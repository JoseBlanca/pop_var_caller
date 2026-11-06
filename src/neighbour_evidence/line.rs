//! **Reading a finished VCF line into the pass, and writing it back re-called**
//! (`doc/devel/ng/spec/neighbour_evidence_pass.md` §2.1, §3.2, §4.1).
//!
//! The pass works on the lines the writer would write, so the same code serves a run with the
//! hidden-duplication filter on, one with it off, and a finished VCF. Reading takes what the walk
//! needs — `POS`, `AF`, each sample's `PL` and `GP` — from a PASS record with rows; any other line
//! passes through untouched. Writing replaces, for every sample the pass re-called, `GT`, `GQ`,
//! `GP` and the repeat counts, recounts `AC` and `AN`, and drops the alternative alleles no
//! sample is called with any more, subsetting every per-allele and per-genotype field as calling
//! does (`vcf_output.md` §7.2). `PL`, `DP`, `QUAL`, `FILTER` and every other field are kept.
//!
//! **The phase** (`neighbour_evidence_phase_output.md` §3): a column the pass phased is written
//! `a|b`, copy 0's allele first, and its block in a `PS` field appended after every other; `PS` is
//! `.` in every other column, and the repeat counts follow the written `GT`'s order.

use super::model::GenotypePairs;
use super::records::{SampleColumns, SiteRecord};
use super::refine::{Outcome, RefinedRecord};
use crate::calling::quality::score_best_genotype;
use crate::vcf::encode::{posteriors_text, written_genotype_quality};
use crate::vcf::line_fields::{
    ALT_COLUMN, COLUMNS_BEFORE_THE_SAMPLES, FORMAT_COLUMN, INFO_COLUMN, LineFields, parse_genotype,
};

/// `INFO` keys with one value per alternative allele (`Number=A` in the header).
const PER_ALTERNATIVE_INFO: [&str; 4] = ["AF", "AC", "MQALT", "MQDIFF"];

/// **What the pass did to the lines it rewrote**, summed over a run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RewriteTally {
    /// Lines that were sites of the walk.
    pub sites: u64,
    /// Sample columns re-called and written with a genotype.
    pub called: u64,
    /// Of those, columns whose reads said nothing — `./.` before the pass (§4.2).
    pub called_without_reads: u64,
    /// Columns of samples without reads re-called and written `./.`: no evidence nearby, `GQ` below
    /// the floor, or probably copying the background (§4.2).
    pub withheld: u64,
    /// Re-called columns whose genotype differs from calling's.
    pub genotypes_changed: u64,
    /// Alternative alleles no sample was called with any more, dropped from `ALT`.
    pub alleles_dropped: u64,
    /// Records left with no alternative allele called, dropped (§4.1, ruled §11 Q2).
    pub records_dropped: u64,
    /// Repeat-count values written `.` because no input column named the allele (§4.1).
    pub repeat_counts_unknown: u64,
}

impl RewriteTally {
    /// Adds `other`'s counts to these. **Destructured exhaustively**, so a count the tally gains
    /// later stops the build here instead of being dropped when lines are rewritten in parallel.
    pub fn add(&mut self, other: RewriteTally) {
        let RewriteTally {
            sites,
            called,
            called_without_reads,
            withheld,
            genotypes_changed,
            alleles_dropped,
            records_dropped,
            repeat_counts_unknown,
        } = other;
        self.sites += sites;
        self.called += called;
        self.called_without_reads += called_without_reads;
        self.withheld += withheld;
        self.genotypes_changed += genotypes_changed;
        self.alleles_dropped += alleles_dropped;
        self.records_dropped += records_dropped;
        self.repeat_counts_unknown += repeat_counts_unknown;
    }
}

/// **The walk's view of a line**, or `None` for a line the pass passes through: not PASS, no
/// alternative allele, no `PL`/`GP`/`AF`, no sample with rows, or genotype rows that are not
/// diploid.
#[must_use]
pub fn site_of(line: &[u8], samples: usize) -> Option<SiteRecord> {
    let fields = LineFields::parse(line)?;
    if fields.column_count() != samples || !fields.is_pass() || fields.fixed[ALT_COLUMN] == "." {
        return None;
    }
    let position: u64 = fields.fixed[1].parse().ok()?;
    let allele_count = fields.fixed[ALT_COLUMN].split(',').count() + 1;
    let genotype_count = allele_count * (allele_count + 1) / 2;
    let phred_at = fields.format_key("PL")?;
    let posterior_at = fields.format_key("GP")?;
    let alt_frequencies: Vec<f64> = info_values(fields.fixed[INFO_COLUMN], "AF")?
        .split(',')
        .map(str::parse)
        .collect::<Result<_, _>>()
        .ok()?;
    if alt_frequencies.len() + 1 != allele_count {
        return None;
    }
    let mut parsed: Vec<Option<(Vec<u32>, Vec<f64>)>> = Vec::with_capacity(samples);
    for column in fields.columns() {
        let values: Vec<&str> = column.split(':').collect();
        let phred = values.get(phred_at).copied().unwrap_or(".");
        let posteriors = values.get(posterior_at).copied().unwrap_or(".");
        if phred == "." || posteriors == "." {
            parsed.push(None);
            continue;
        }
        let phred: Vec<u32> = phred
            .split(',')
            .map(str::parse)
            .collect::<Result<_, _>>()
            .ok()?;
        let posteriors: Vec<f64> = posteriors
            .split(',')
            .map(str::parse)
            .collect::<Result<_, _>>()
            .ok()?;
        if phred.len() != genotype_count || posteriors.len() != genotype_count {
            return None;
        }
        parsed.push(Some((phred, posteriors)));
    }
    if parsed.iter().all(Option::is_none) {
        return None;
    }
    let columns: Vec<Option<SampleColumns<'_>>> = parsed
        .iter()
        .map(|column| {
            column.as_ref().map(|(phred, posteriors)| SampleColumns {
                phred_likelihoods: phred,
                genotype_posteriors: posteriors,
            })
        })
        .collect();
    Some(SiteRecord::new(position, &alt_frequencies, &columns))
}

/// The value of `key` in an `INFO` column, if it is there.
fn info_values<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';').find_map(|entry| {
        entry
            .strip_prefix(key)
            .and_then(|rest| rest.strip_prefix('='))
    })
}

/// **Rewrites a site's line with the pass's result** — `None` when no alternative allele is left
/// called and the record is dropped. `refined` is the line's [`RefinedRecord`]; `pairs` a
/// genotype-order cache the caller keeps.
///
/// # Panics
///
/// If the line is not one [`site_of`] read as a site.
#[must_use]
pub(crate) fn rewrite_line(
    line: &[u8],
    refined: &RefinedRecord,
    pairs: &mut GenotypePairs,
    tally: &mut RewriteTally,
) -> Option<Vec<u8>> {
    let text = std::str::from_utf8(line).expect("a site line is text");
    let fields: Vec<&str> = text.split('\t').collect();
    let alternatives: Vec<&str> = fields[ALT_COLUMN].split(',').collect();
    let allele_count = alternatives.len() + 1;
    let old_pairs: Vec<(usize, usize)> = pairs.of(allele_count).to_vec();
    let genotype_count = old_pairs.len();
    let format: Vec<&str> = fields[FORMAT_COLUMN].split(':').collect();
    let key = |name: &str| format.iter().position(|field| *field == name);
    let genotype_at = key("GT").expect("a site has GT");
    let quality_at = key("GQ");
    let depths_at = key("AD");
    let phred_at = key("PL").expect("a site has PL");
    let posterior_at = key("GP").expect("a site has GP");
    let repeats_at = key("REPCN").or_else(|| key("REPCNB"));
    // `PS` goes after every field the line has, so each field the pass rewrites keeps its place.
    let phase_set_at = key("PS").unwrap_or(format.len());
    let field_count = format.len().max(phase_set_at + 1);
    let mut columns: Vec<Vec<String>> = fields[COLUMNS_BEFORE_THE_SAMPLES..]
        .iter()
        .map(|column| {
            let mut values: Vec<String> = column.split(':').map(str::to_string).collect();
            // A column may drop trailing fields; `PS` must land in its own place. A column with more
            // fields than `FORMAT` is kept as it is.
            if values.len() < field_count {
                values.resize(field_count, ".".to_string());
            }
            values
        })
        .collect();
    tally.sites += 1;

    // Each allele's repeat count, from the first input column that names it (§4.1).
    let mut repeat_counts: Vec<Option<String>> = vec![None; allele_count];
    if let Some(repeats_at) = repeats_at {
        for column in &columns {
            let (Some((first, second)), Some(values)) =
                (parse_genotype(&column[genotype_at]), column.get(repeats_at))
            else {
                continue;
            };
            if values == "." {
                continue;
            }
            for (allele, value) in [first, second].into_iter().zip(values.split(',')) {
                if allele < allele_count && repeat_counts[allele].is_none() && value != "." {
                    repeat_counts[allele] = Some(value.to_string());
                }
            }
        }
    }

    // The genotype every column will carry: the pass's argmax where it called, `./.` where it
    // withheld, calling's own elsewhere.
    let genotypes: Vec<Option<(usize, usize)>> = columns
        .iter()
        .zip(&refined.outcomes)
        .enumerate()
        .map(|(sample, (column, outcome))| match outcome {
            Outcome::Unchanged => parse_genotype(&column[genotype_at]),
            Outcome::Withheld => None,
            Outcome::Called => {
                let posterior =
                    &refined.posteriors[sample * genotype_count..(sample + 1) * genotype_count];
                let (best, _) = score_best_genotype(posterior);
                Some(old_pairs[best.0 as usize])
            }
        })
        .collect();
    let mut carried = vec![false; allele_count];
    for &(first, second) in genotypes.iter().flatten() {
        carried[first] = true;
        carried[second] = true;
    }
    let anyone_called = genotypes.iter().any(Option::is_some);
    if anyone_called && !carried[1..].iter().any(|&carried| carried) {
        tally.records_dropped += 1;
        return None;
    }
    // As calling: with nobody called there is nothing to say an allele is uncarried.
    let kept: Vec<usize> = (0..allele_count)
        .filter(|&allele| allele == 0 || !anyone_called || carried[allele])
        .collect();
    let subsetting = kept.len() < allele_count;
    tally.alleles_dropped += (allele_count - kept.len()) as u64;
    let mut new_index = vec![usize::MAX; allele_count];
    for (index, &allele) in kept.iter().enumerate() {
        new_index[allele] = index;
    }
    // Each kept genotype, in the new table's order, as its index in the old one.
    let kept_genotypes: Vec<usize> = pairs
        .of(kept.len())
        .iter()
        .map(|&(first, second)| {
            let old = (kept[first], kept[second]);
            old_pairs
                .iter()
                .position(|&pair| pair == old)
                .expect("a pair of kept alleles is an old genotype")
        })
        .collect();

    for (sample, column) in columns.iter_mut().enumerate() {
        let outcome = refined.outcomes[sample];
        let recalled = outcome != Outcome::Unchanged;
        // The phase, where the pass wrote one for this very genotype: copy 0's allele first.
        let phase = refined.phases[sample].filter(|phase| {
            outcome == Outcome::Called
                && genotypes[sample]
                    == Some((
                        phase.copy_0_allele.min(phase.copy_1_allele),
                        phase.copy_0_allele.max(phase.copy_1_allele),
                    ))
        });
        let genotype_text = |genotype: Option<(usize, usize)>| match (genotype, phase) {
            (Some(_), Some(phase)) => {
                format!(
                    "{}|{}",
                    new_index[phase.copy_0_allele], new_index[phase.copy_1_allele]
                )
            }
            (Some((first, second)), None) => {
                format!("{}/{}", new_index[first], new_index[second])
            }
            (None, _) => "./.".to_string(),
        };
        if recalled || subsetting {
            column[genotype_at] = genotype_text(genotypes[sample]);
        }
        column[phase_set_at] =
            phase.map_or_else(|| ".".to_string(), |phase| phase.phase_set.to_string());
        let posteriors: Option<Vec<f64>> = if recalled {
            Some(
                refined.posteriors[sample * genotype_count..(sample + 1) * genotype_count].to_vec(),
            )
        } else if subsetting && column[posterior_at] != "." {
            column[posterior_at]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()
                .ok()
        } else {
            None
        };
        if let Some(posteriors) = posteriors {
            let mut kept_posteriors: Vec<f64> =
                kept_genotypes.iter().map(|&old| posteriors[old]).collect();
            let total: f64 = kept_posteriors.iter().sum();
            if total > 0.0 {
                for value in &mut kept_posteriors {
                    *value /= total;
                }
            }
            if let (Some(quality_at), Outcome::Called) = (quality_at, outcome) {
                let (_, quality) = score_best_genotype(&kept_posteriors);
                column[quality_at] = written_genotype_quality(quality).to_string();
            }
            column[posterior_at] = posteriors_text(&kept_posteriors);
        }
        if let (Some(quality_at), Outcome::Withheld) = (quality_at, outcome) {
            column[quality_at] = ".".to_string();
        }
        if subsetting && column[phred_at] != "." {
            let phred: Vec<i64> = column[phred_at]
                .split(',')
                .map(|value| value.parse().expect("PL is whole numbers"))
                .collect();
            let kept_phred: Vec<i64> = kept_genotypes.iter().map(|&old| phred[old]).collect();
            let lowest = kept_phred.iter().copied().min().unwrap_or(0);
            column[phred_at] = join(kept_phred.iter().map(|value| value - lowest));
        }
        if let Some(depths_at) = depths_at.filter(|_| subsetting)
            && column[depths_at] != "."
        {
            let depths: Vec<&str> = column[depths_at].split(',').collect();
            column[depths_at] = kept
                .iter()
                .map(|&allele| depths[allele])
                .collect::<Vec<_>>()
                .join(",");
        }
        if let (Some(repeats_at), true) = (repeats_at, recalled) {
            column[repeats_at] = match genotypes[sample] {
                None => ".".to_string(),
                Some(pair) => phase
                    .map_or([pair.0, pair.1], |phase| {
                        [phase.copy_0_allele, phase.copy_1_allele]
                    })
                    .iter()
                    .map(|&allele| {
                        repeat_counts[allele].clone().unwrap_or_else(|| {
                            tally.repeat_counts_unknown += 1;
                            ".".to_string()
                        })
                    })
                    .collect::<Vec<_>>()
                    .join(","),
            };
        }
        match outcome {
            Outcome::Called => tally.called += 1,
            Outcome::Withheld => tally.withheld += 1,
            Outcome::Unchanged => {}
        }
    }
    // Calling's genotypes, to count the ones the pass moved and the calls made without reads.
    for (sample, column) in fields[COLUMNS_BEFORE_THE_SAMPLES..].iter().enumerate() {
        if refined.outcomes[sample] != Outcome::Called {
            continue;
        }
        let before = parse_genotype(column.split(':').nth(genotype_at).unwrap_or("."));
        if before.is_none() {
            tally.called_without_reads += 1;
        }
        if before != genotypes[sample] {
            tally.genotypes_changed += 1;
        }
    }

    // INFO: AC and AN recounted; per-alternative fields subset.
    let mut copies = vec![0_u64; allele_count];
    let mut called_copies = 0_u64;
    for &(first, second) in genotypes.iter().flatten() {
        copies[first] += 1;
        copies[second] += 1;
        called_copies += 2;
    }
    let info: Vec<String> = fields[INFO_COLUMN]
        .split(';')
        .map(|entry| {
            let (name, value) = entry.split_once('=').unwrap_or((entry, ""));
            match name {
                "AC" => format!(
                    "AC={}",
                    join(kept[1..].iter().map(|&allele| copies[allele]))
                ),
                "AN" => format!("AN={called_copies}"),
                _ if subsetting && PER_ALTERNATIVE_INFO.contains(&name) => {
                    let values: Vec<&str> = value.split(',').collect();
                    format!(
                        "{name}={}",
                        kept[1..]
                            .iter()
                            .map(|&allele| values.get(allele - 1).copied().unwrap_or("."))
                            .collect::<Vec<_>>()
                            .join(",")
                    )
                }
                _ => entry.to_string(),
            }
        })
        .collect();

    let alt = kept[1..]
        .iter()
        .map(|&allele| alternatives[allele - 1])
        .collect::<Vec<_>>()
        .join(",");
    let mut out: Vec<String> = fields[..COLUMNS_BEFORE_THE_SAMPLES]
        .iter()
        .map(|field| (*field).to_string())
        .collect();
    out[ALT_COLUMN] = alt;
    out[INFO_COLUMN] = info.join(";");
    if phase_set_at == format.len() {
        out[FORMAT_COLUMN] = format!("{}:PS", fields[FORMAT_COLUMN]);
    }
    out.extend(columns.iter().map(|column| column.join(":")));
    Some(out.join("\t").into_bytes())
}

fn join<T: std::fmt::Display>(values: impl Iterator<Item = T>) -> String {
    values
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
mod tests;
