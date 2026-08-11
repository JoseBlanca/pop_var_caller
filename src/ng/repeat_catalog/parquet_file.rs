//! The catalog's wire format: a Parquet file with **one row group per contig** and the
//! header in the footer's key-value metadata.
//!
//! Why a column store rather than a compressed TSV is spec §3.5; what it buys here is that
//! the two heaviest readers touch three or four of the seven columns, every column encodes
//! to almost nothing (coordinates ascend within a row group, motifs are drawn from at most
//! 5,460 values), row groups replace a per-contig byte-offset table, and a build that dies
//! mid-write leaves a file with no footer — unreadable rather than short-but-valid.
//!
//! **Arrow types do not leave this module.** Everything above it speaks [`FoundRepeat`] and
//! [`RepeatCatalogHeader`], which is what keeps a format change one module deep.

use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow_array::builder::{
    Float32Builder, Int32Builder, StringBuilder, UInt8Builder, UInt32Builder, UInt64Builder,
};
use arrow_array::{Array, ArrayRef, RecordBatch};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use parquet::arrow::ArrowWriter;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;

use crate::ng::region_typing::segment_criteria::{MinCopies, SsrSegmentCriteria};
use crate::ng::repeat_catalog::criteria::StrRepeatCriteria;
use crate::ng::repeat_catalog::{
    FoundRepeat, RepeatCatalogError, RepeatCatalogHeader, RowsByPeriod, TractSpan,
};
use crate::ng::tandem_repeat::{PeriodRange, ScanParams};
use crate::ng::types::{Bp, ContigId, Motif, Position};

/// The key the header lives under in the file's footer metadata.
pub(crate) const HEADER_KEY: &str = "ng_repeat_catalog";

/// The header encoding's own version, bumped when its grammar changes. Distinct from the
/// tool version: a file may be readable by a later tool that writes the same grammar.
pub(crate) const HEADER_FORMAT_VERSION: u32 = 2;

/// What the writer stamps into the file's `created_by`.
///
/// **Fixed by us on purpose.** The default carries the parquet crate's version, so two
/// builds of the same reference would differ in bytes after a routine dependency bump —
/// and spec §6 asks for a byte-identical file. What wrote the catalog is the header's
/// `tool_version`, which is checked on read.
pub(crate) const CREATED_BY: &str = "ng repeat-catalog";

/// The compression level, fixed for the same reason as [`CREATED_BY`].
const ZSTD_LEVEL: i32 = 3;

/// The seven columns of a catalog row (spec §3.1).
///
/// `contig` is the header's contig index rather than a name: rows address contigs by
/// index, the names are in the header, and inside a row group — one contig — the column is
/// a single repeated value.
pub(crate) fn catalog_schema() -> SchemaRef {
    Arc::new(Schema::new(vec![
        Field::new("contig", DataType::UInt32, false),
        Field::new("detected_start", DataType::UInt64, false),
        Field::new("detected_end", DataType::UInt64, false),
        Field::new("trimmed_start", DataType::UInt64, true),
        Field::new("trimmed_end", DataType::UInt64, true),
        Field::new("period", DataType::UInt8, false),
        Field::new("score", DataType::Int32, false),
        Field::new("motif", DataType::Utf8, false),
        Field::new("purity", DataType::Float32, true),
    ]))
}

/// Writes a catalog: one row group per contig, then the header, then a rename into place.
///
/// **Atomic**: rows go to a sibling `.tmp` file that is renamed only once the footer is
/// written, so a build that dies leaves nothing under the name a reader will look for
/// (spec §6) — and even that leftover `.tmp` has no footer, so it cannot be read as a short
/// catalog.
pub struct RepeatCatalogWriter {
    writer: ArrowWriter<File>,
    schema: SchemaRef,
    final_path: PathBuf,
    temp_path: PathBuf,
    rows_by_period: RowsByPeriod,
}

impl RepeatCatalogWriter {
    /// Open a writer for `path`, creating the temporary file it builds into.
    pub fn create(path: &Path) -> Result<Self, RepeatCatalogError> {
        let temp_path = temp_path_for(path);
        let file = File::create(&temp_path).map_err(|source| RepeatCatalogError::Write {
            path: temp_path.clone(),
            source,
        })?;

        let properties = WriterProperties::builder()
            .set_created_by(CREATED_BY.to_string())
            .set_compression(Compression::ZSTD(
                ZstdLevel::try_new(ZSTD_LEVEL).expect("zstd level 3 is valid"),
            ))
            .set_dictionary_enabled(true)
            .build();

        let schema = catalog_schema();
        let writer = ArrowWriter::try_new(file, Arc::clone(&schema), Some(properties))
            .map_err(|e| unreadable(&temp_path, e))?;

        Ok(Self {
            writer,
            schema,
            final_path: path.to_path_buf(),
            temp_path,
            rows_by_period: RowsByPeriod::default(),
        })
    }

    /// Write one contig's rows as **one row group**, in the order given.
    ///
    /// The caller owns the order (spec §3.1: contig, start, period, end) and the grouping;
    /// this flushes after every call, which is what makes a row group a contig.
    /// An empty slice writes nothing at all rather than an empty row group.
    pub fn write_contig(&mut self, rows: &[FoundRepeat]) -> Result<(), RepeatCatalogError> {
        if rows.is_empty() {
            return Ok(());
        }
        for row in rows {
            self.rows_by_period.count(row.period);
        }
        let batch = rows_to_batch(&self.schema, rows);
        self.writer
            .write(&batch)
            .map_err(|e| unreadable(&self.temp_path, e))?;
        self.writer
            .flush()
            .map_err(|e| unreadable(&self.temp_path, e))?;
        Ok(())
    }

    /// Write the header, close the file and move it into place. Returns what was written,
    /// which is the tally the build command prints (spec §2.6).
    pub fn finish(
        mut self,
        header: &RepeatCatalogHeader,
    ) -> Result<RowsByPeriod, RepeatCatalogError> {
        self.writer
            .append_key_value_metadata(parquet::file::metadata::KeyValue::new(
                HEADER_KEY.to_string(),
                Some(encode_header(header)),
            ));
        self.writer
            .close()
            .map_err(|e| unreadable(&self.temp_path, e))?;

        std::fs::rename(&self.temp_path, &self.final_path).map_err(|source| {
            RepeatCatalogError::Write {
                path: self.final_path.clone(),
                source,
            }
        })?;
        Ok(self.rows_by_period)
    }
}

/// Where a half-written catalog lives until its footer exists.
pub(crate) fn temp_path_for(path: &Path) -> PathBuf {
    let mut name = path.as_os_str().to_os_string();
    name.push(".tmp");
    PathBuf::from(name)
}

fn unreadable(
    path: &Path,
    source: impl std::error::Error + Send + Sync + 'static,
) -> RepeatCatalogError {
    RepeatCatalogError::Unreadable {
        path: path.to_path_buf(),
        source: Box::new(source),
    }
}

/// One contig's rows as a column batch.
fn rows_to_batch(schema: &SchemaRef, rows: &[FoundRepeat]) -> RecordBatch {
    let n = rows.len();
    let mut contig = UInt32Builder::with_capacity(n);
    let mut detected_start = UInt64Builder::with_capacity(n);
    let mut detected_end = UInt64Builder::with_capacity(n);
    let mut trimmed_start = UInt64Builder::with_capacity(n);
    let mut trimmed_end = UInt64Builder::with_capacity(n);
    let mut period = UInt8Builder::with_capacity(n);
    let mut score = Int32Builder::with_capacity(n);
    let mut motif = StringBuilder::with_capacity(n, n * crate::ng::types::MAX_MOTIF_LEN);
    let mut purity = Float32Builder::with_capacity(n);

    for row in rows {
        contig.append_value(row.contig.get());
        detected_start.append_value(row.detected.start.get());
        detected_end.append_value(row.detected.end.get());
        trimmed_start.append_option(row.trimmed.map(|span| span.start.get()));
        trimmed_end.append_option(row.trimmed.map(|span| span.end.get()));
        period.append_value(row.period);
        score.append_value(row.score);
        // A motif is ASCII bases by construction (`Motif::new` takes reference bytes), so
        // this cannot be lossy.
        motif.append_value(String::from_utf8_lossy(row.motif.as_bytes()));
        purity.append_option(row.purity);
    }

    let columns: Vec<ArrayRef> = vec![
        Arc::new(contig.finish()),
        Arc::new(detected_start.finish()),
        Arc::new(detected_end.finish()),
        Arc::new(trimmed_start.finish()),
        Arc::new(trimmed_end.finish()),
        Arc::new(period.finish()),
        Arc::new(score.finish()),
        Arc::new(motif.finish()),
        Arc::new(purity.finish()),
    ];
    RecordBatch::try_new(Arc::clone(schema), columns)
        .expect("the builders were filled in schema order, with one value per row")
}

// ---------------------------------------------------------------------------
// The header, as text in the footer's metadata
// ---------------------------------------------------------------------------

/// Encode the header as tab-separated lines.
///
/// **Text with checked parsing, rather than a derived serialisation**, because two of the
/// values it carries — the period range and the copy-floor table — have checked
/// constructors, and a derive would rebuild them field-by-field and skip the check. Decode
/// goes back through `PeriodRange::new` and `MinCopies::new` (spec §3.4).
pub(crate) fn encode_header(header: &RepeatCatalogHeader) -> String {
    let mut out = String::new();
    out.push_str(&format!("format\t{HEADER_FORMAT_VERSION}\n"));
    out.push_str(&format!("tool_version\t{}\n", header.tool_version));
    out.push_str(&format!("reference_md5\t{}\n", hex(&header.reference_md5)));
    out.push_str(&format!(
        "scan\t{}\t{}\t{}\n",
        header.scan.match_reward, header.scan.mismatch_penalty, header.scan.min_copies
    ));

    let criteria = &header.built_under;
    let class = &criteria.classification;
    let floors: Vec<String> = (1..=crate::ng::types::MAX_MOTIF_LEN as u8)
        .map(|p| class.min_copies.for_period(p).to_string())
        .collect();
    out.push_str(&format!(
        "criteria\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        class.periods.min(),
        class.periods.max(),
        floors.join(","),
        class.min_copies.for_period(u8::MAX), // the beyond-the-table fallback
        class.min_purity,
        class.min_score,
        class.bundle_threshold,
        criteria.min_flank_bp.get(),
        criteria.max_str_len_bp.get(),
    ));

    for (index, contig) in header.contigs.iter().enumerate() {
        out.push_str(&format!(
            "contig\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            contig.name,
            contig.length,
            contig
                .md5
                .map(|d| hex(&d))
                .unwrap_or_else(|| "-".to_string()),
            contig.offset,
            contig.line_bases,
            contig.line_width,
            header.longest_tract_bp.get(index).copied().unwrap_or(0),
        ));
    }
    out
}

/// The inverse of [`encode_header`]. Every failure is a malformed header, which is a
/// malformed file.
pub fn decode_header(text: &str, path: &Path) -> Result<RepeatCatalogHeader, RepeatCatalogError> {
    let malformed = |what: &str| RepeatCatalogError::Unreadable {
        path: path.to_path_buf(),
        source: format!("catalog header: {what}").into(),
    };

    let mut tool_version = None;
    let mut reference_md5 = None;
    let mut scan = None;
    let mut built_under = None;
    let mut contigs = Vec::new();
    let mut longest_tract_bp = Vec::new();

    for line in text.lines().filter(|l| !l.is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "format" => {
                let version: u32 = f
                    .get(1)
                    .and_then(|v| v.parse().ok())
                    .ok_or_else(|| malformed("unreadable format version"))?;
                if version != HEADER_FORMAT_VERSION {
                    return Err(malformed(&format!(
                        "written in header format {version}, this build reads {HEADER_FORMAT_VERSION}"
                    )));
                }
            }
            "tool_version" => tool_version = f.get(1).map(|v| (*v).to_string()),
            "reference_md5" => {
                reference_md5 = f.get(1).and_then(|v| unhex(v));
            }
            "scan" if f.len() >= 4 => {
                scan = Some(ScanParams {
                    match_reward: f[1].parse().map_err(|_| malformed("scan match reward"))?,
                    mismatch_penalty: f[2]
                        .parse()
                        .map_err(|_| malformed("scan mismatch penalty"))?,
                    min_copies: f[3].parse().map_err(|_| malformed("scan min copies"))?,
                });
            }
            "criteria" if f.len() >= 10 => {
                let period_min: u8 = f[1].parse().map_err(|_| malformed("period floor"))?;
                let period_max: u8 = f[2].parse().map_err(|_| malformed("period ceiling"))?;
                let mut by_period = [0u32; crate::ng::types::MAX_MOTIF_LEN];
                let parsed: Vec<&str> = f[3].split(',').collect();
                if parsed.len() != by_period.len() {
                    return Err(malformed("copy-floor table is not six wide"));
                }
                for (slot, text) in by_period.iter_mut().zip(parsed) {
                    *slot = text.parse().map_err(|_| malformed("copy floor"))?;
                }
                let beyond: u32 = f[4].parse().map_err(|_| malformed("copy-floor fallback"))?;

                built_under = Some(StrRepeatCriteria {
                    classification: SsrSegmentCriteria {
                        periods: PeriodRange::new(period_min, period_max)
                            .map_err(|_| malformed("period range"))?,
                        min_copies: MinCopies::new(by_period, beyond),
                        min_purity: f[5].parse().map_err(|_| malformed("purity floor"))?,
                        min_score: f[6].parse().map_err(|_| malformed("score floor"))?,
                        bundle_threshold: f[7].parse().map_err(|_| malformed("bundle radius"))?,
                    },
                    min_flank_bp: Bp(f[8].parse().map_err(|_| malformed("flank floor"))?),
                    max_str_len_bp: Bp(f[9].parse().map_err(|_| malformed("satellite cap"))?),
                });
            }
            "contig" if f.len() >= 8 => {
                contigs.push(crate::ng::reference_info::ContigInfo {
                    name: f[1].to_string(),
                    length: f[2].parse().map_err(|_| malformed("contig length"))?,
                    md5: if f[3] == "-" { None } else { unhex(f[3]) },
                    offset: f[4].parse().map_err(|_| malformed("contig offset"))?,
                    line_bases: f[5].parse().map_err(|_| malformed("contig line bases"))?,
                    line_width: f[6].parse().map_err(|_| malformed("contig line width"))?,
                });
                longest_tract_bp.push(f[7].parse().map_err(|_| malformed("longest tract"))?);
            }
            other => return Err(malformed(&format!("unknown line `{other}`"))),
        }
    }

    Ok(RepeatCatalogHeader {
        contigs,
        longest_tract_bp,
        reference_md5: reference_md5.ok_or_else(|| malformed("no reference digest"))?,
        built_under: built_under.ok_or_else(|| malformed("no criteria"))?,
        scan: scan.ok_or_else(|| malformed("no scan settings"))?,
        tool_version: tool_version.ok_or_else(|| malformed("no tool version"))?,
    })
}

fn hex(digest: &[u8; 16]) -> String {
    digest.iter().map(|b| format!("{b:02x}")).collect()
}

fn unhex(text: &str) -> Option<[u8; 16]> {
    if text.len() != 32 {
        return None;
    }
    let mut out = [0u8; 16];
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = u8::from_str_radix(&text[i * 2..i * 2 + 2], 16).ok()?;
    }
    Some(out)
}

/// Every row of `contig` whose tract overlaps one of `windows`, in file order.
///
/// **The file's own index does the skipping.** The page index — which the writer emits by
/// default, since `EnabledStatistics::Page` is parquet's default — carries each page's
/// minimum and maximum per column, and rows are written in ascending start within a contig,
/// so a filter on `detected_start` / `detected_end` lets the reader drop whole pages without
/// decompressing them. Asking for a megabase of a 90 Mb chromosome reads the pages that
/// megabase touches, not the chromosome.
///
/// `windows` are 1-based inclusive and already widened by the caller (see
/// `RepeatCatalog::rows_of_contig_in`) — this function does no widening of its own, because
/// how far a row can reach is a question about classification, not about the file.
pub fn rows_overlapping(
    path: &Path,
    contig: ContigId,
    windows: &[(u64, u64)],
) -> Result<Vec<FoundRepeat>, Box<dyn std::error::Error + Send + Sync>> {
    use arrow_array::{BooleanArray, UInt32Array, UInt64Array};
    use parquet::arrow::ProjectionMask;
    use parquet::arrow::arrow_reader::{ArrowPredicateFn, ArrowReaderOptions, RowFilter};
    use parquet::file::metadata::PageIndexPolicy;

    let file = File::open(path)?;
    // `Required`, not `Optional`: our own writer always emits the page index, so a file
    // without one is not a catalog this build wrote, and silently falling back to a full
    // scan would turn a wrong file into a slow answer.
    let options = ArrowReaderOptions::new().with_page_index_policy(PageIndexPolicy::Required);
    let builder = ParquetRecordBatchReaderBuilder::try_new_with_options(file, options)?;

    // The contig's own row group, found from the row-group statistics rather than by
    // reading rows.
    let groups: Vec<usize> = (0..builder.metadata().num_row_groups())
        .filter(|group| {
            matches!(
                builder.metadata().row_group(*group).column(0).statistics(),
                Some(stats) if statistics_min_u32(stats) == Some(contig.get())
            )
        })
        .collect();
    if groups.is_empty() {
        return Ok(Vec::new());
    }

    let schema = builder.parquet_schema();
    // Columns 1 and 2 are `detected_start` and `detected_end` (`catalog_schema`).
    let coordinate_columns = ProjectionMask::roots(schema, [0, 1, 2]);
    let windows = windows.to_vec();
    let predicate = ArrowPredicateFn::new(coordinate_columns, move |batch| {
        let contig_column = batch
            .column(0)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .expect("the contig column is UInt32");
        let starts = batch
            .column(1)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .expect("detected_start is UInt64");
        let ends = batch
            .column(2)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .expect("detected_end is UInt64");
        Ok(BooleanArray::from_iter((0..batch.num_rows()).map(|i| {
            let start = starts.value(i);
            let end = ends.value(i);
            Some(
                contig_column.value(i) == contig.get()
                    && windows
                        .iter()
                        .any(|(from, to)| start <= *to && end >= *from),
            )
        })))
    });

    let reader = builder
        .with_row_groups(groups)
        .with_row_filter(RowFilter::new(vec![Box::new(predicate)]))
        .build()?;

    let mut rows = Vec::new();
    for batch in reader {
        let batch = batch?;
        for i in 0..batch.num_rows() {
            rows.push(row_from_batch(&batch, i).ok_or("a row does not match the catalog schema")?);
        }
    }
    Ok(rows)
}

/// The minimum of a column chunk's statistics, when it is an integer column.
fn statistics_min_u32(stats: &parquet::file::statistics::Statistics) -> Option<u32> {
    use parquet::file::statistics::Statistics;
    match stats {
        Statistics::Int32(s) => s.min_opt().map(|v| *v as u32),
        Statistics::Int64(s) => s.min_opt().map(|v| *v as u32),
        _ => None,
    }
}

/// Rebuild one row from the column batch's `i`-th entry. Shared with the reader (B2).
pub fn row_from_batch(batch: &RecordBatch, i: usize) -> Option<FoundRepeat> {
    use arrow_array::{
        Float32Array, Int32Array, StringArray, UInt8Array, UInt32Array, UInt64Array,
    };

    let contig = batch.column(0).as_any().downcast_ref::<UInt32Array>()?;
    let detected_start = batch.column(1).as_any().downcast_ref::<UInt64Array>()?;
    let detected_end = batch.column(2).as_any().downcast_ref::<UInt64Array>()?;
    let trimmed_start = batch.column(3).as_any().downcast_ref::<UInt64Array>()?;
    let trimmed_end = batch.column(4).as_any().downcast_ref::<UInt64Array>()?;
    let period = batch.column(5).as_any().downcast_ref::<UInt8Array>()?;
    let score = batch.column(6).as_any().downcast_ref::<Int32Array>()?;
    let motif = batch.column(7).as_any().downcast_ref::<StringArray>()?;
    let purity = batch.column(8).as_any().downcast_ref::<Float32Array>()?;

    let trimmed = if trimmed_start.is_null(i) || trimmed_end.is_null(i) {
        None
    } else {
        Some(TractSpan {
            start: Position(trimmed_start.value(i)),
            end: Position(trimmed_end.value(i)),
        })
    };

    Some(FoundRepeat {
        contig: ContigId(contig.value(i)),
        detected: TractSpan {
            start: Position(detected_start.value(i)),
            end: Position(detected_end.value(i)),
        },
        trimmed,
        period: period.value(i),
        score: score.value(i),
        motif: Motif::new(motif.value(i).as_bytes()).ok()?,
        purity: (!purity.is_null(i)).then(|| purity.value(i)),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::reference_info::ContigInfo;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    fn header() -> RepeatCatalogHeader {
        RepeatCatalogHeader {
            contigs: vec![
                ContigInfo {
                    name: "chr1".to_string(),
                    length: 1_000,
                    offset: 6,
                    line_bases: 60,
                    line_width: 61,
                    md5: Some([7u8; 16]),
                },
                ContigInfo {
                    name: "chr2".to_string(),
                    length: 500,
                    offset: 1_024,
                    line_bases: 60,
                    line_width: 61,
                    md5: Some([9u8; 16]),
                },
            ],
            reference_md5: [3u8; 16],
            built_under: StrRepeatCriteria::default(),
            scan: ScanParams::default(),
            tool_version: "test-0.1".to_string(),
            longest_tract_bp: vec![0; 2],
        }
    }

    fn row(contig: u32, start: u64, end: u64, period: u8, trimmed: bool) -> FoundRepeat {
        FoundRepeat {
            contig: ContigId(contig),
            detected: TractSpan {
                start: Position(start),
                end: Position(end),
            },
            trimmed: trimmed.then_some(TractSpan {
                start: Position(start),
                end: Position(end),
            }),
            period,
            score: 30,
            motif: Motif::new(b"CAG").expect("valid"),
            purity: trimmed.then_some(0.95),
        }
    }

    fn write_catalog(path: &Path, contigs: &[Vec<FoundRepeat>]) -> RowsByPeriod {
        let mut writer = RepeatCatalogWriter::create(path).expect("writer");
        for rows in contigs {
            writer.write_contig(rows).expect("row group");
        }
        writer.finish(&header()).expect("footer")
    }

    #[test]
    fn the_header_round_trips_through_its_text_encoding() {
        let encoded = encode_header(&header());
        let decoded = decode_header(&encoded, Path::new("test.parquet")).expect("decodes");
        assert_eq!(decoded, header());
    }

    #[test]
    fn a_header_line_the_grammar_does_not_know_is_a_malformed_file() {
        let mut encoded = encode_header(&header());
        encoded.push_str("something_else\t1\n");
        assert!(decode_header(&encoded, Path::new("test.parquet")).is_err());
    }

    #[test]
    fn a_future_header_format_is_refused_rather_than_guessed_at() {
        // Version-agnostic: whatever this build writes, a file one version ahead is refused.
        let encoded = encode_header(&header()).replace(
            &format!("format\t{HEADER_FORMAT_VERSION}"),
            &format!("format\t{}", HEADER_FORMAT_VERSION + 1),
        );
        let err = decode_header(&encoded, Path::new("test.parquet")).expect_err("refused");
        assert!(
            format!("{err}").contains("unreadable or truncated"),
            "unexpected message: {err}"
        );
    }

    #[test]
    fn rows_survive_a_write_and_a_read() {
        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("catalog.parquet");
        let contig0 = vec![row(0, 100, 117, 3, true), row(0, 300, 317, 3, false)];
        let contig1 = vec![row(1, 50, 59, 1, true)];

        let written = write_catalog(&path, &[contig0.clone(), contig1.clone()]);
        assert_eq!(written.total(), 3);
        assert_eq!(written.for_period(3), 2);
        assert_eq!(written.for_period(1), 1);

        let file = File::open(&path).expect("open");
        let reader = ParquetRecordBatchReaderBuilder::try_new(file)
            .expect("parquet")
            .build()
            .expect("reader");
        let mut read_back = Vec::new();
        for batch in reader {
            let batch = batch.expect("batch");
            for i in 0..batch.num_rows() {
                read_back.push(row_from_batch(&batch, i).expect("row"));
            }
        }

        let expected: Vec<FoundRepeat> = contig0.into_iter().chain(contig1).collect();
        assert_eq!(read_back, expected);
    }

    #[test]
    fn each_contig_is_its_own_row_group() {
        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("catalog.parquet");
        write_catalog(
            &path,
            &[
                vec![row(0, 100, 117, 3, true)],
                vec![row(1, 50, 59, 1, true), row(1, 80, 89, 1, true)],
                vec![row(2, 10, 27, 3, true)],
            ],
        );

        let file = File::open(&path).expect("open");
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).expect("parquet");
        let groups = builder.metadata().num_row_groups();
        assert_eq!(groups, 3, "one row group per contig, so a reader can seek");
        assert_eq!(builder.metadata().row_group(1).num_rows(), 2);
    }

    #[test]
    fn the_header_comes_back_from_the_footer() {
        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("catalog.parquet");
        write_catalog(&path, &[vec![row(0, 100, 117, 3, true)]]);

        let file = File::open(&path).expect("open");
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).expect("parquet");
        let kv = builder
            .metadata()
            .file_metadata()
            .key_value_metadata()
            .expect("footer metadata")
            .iter()
            .find(|kv| kv.key == HEADER_KEY)
            .expect("the catalog header")
            .value
            .clone()
            .expect("a value");

        assert_eq!(
            decode_header(&kv, &path).expect("decodes"),
            header(),
            "the header a reader validates against is the one the builder wrote"
        );
    }

    /// Spec §6: the same rows and the same header must give the same bytes, or the
    /// thread-count invariance of §10.5 has nothing to stand on.
    #[test]
    fn two_writes_of_the_same_rows_are_byte_identical() {
        let dir = tempfile::tempdir().expect("tmp");
        let contigs = vec![
            vec![row(0, 100, 117, 3, true), row(0, 300, 317, 3, false)],
            vec![row(1, 50, 59, 1, true)],
        ];

        let first = dir.path().join("first.parquet");
        let second = dir.path().join("second.parquet");
        write_catalog(&first, &contigs);
        write_catalog(&second, &contigs);

        assert_eq!(
            std::fs::read(&first).expect("first"),
            std::fs::read(&second).expect("second"),
            "the writer's codec, level, created_by and row-group boundaries are all fixed"
        );
    }

    /// A build that dies leaves nothing under the name a reader looks for, and what it does
    /// leave has no footer (spec §6).
    #[test]
    fn a_writer_that_never_finishes_leaves_no_readable_catalog() {
        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("catalog.parquet");

        let mut writer = RepeatCatalogWriter::create(&path).expect("writer");
        writer
            .write_contig(&[row(0, 100, 117, 3, true)])
            .expect("row group");
        drop(writer);

        assert!(!path.exists(), "nothing under the final name");
        let leftover = temp_path_for(&path);
        assert!(leftover.exists(), "the half-written file is the .tmp");
        assert!(
            ParquetRecordBatchReaderBuilder::try_new(File::open(&leftover).expect("open")).is_err(),
            "a footerless parquet file cannot be opened, so it cannot read as a short catalog"
        );
    }
}
