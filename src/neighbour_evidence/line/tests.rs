use super::*;

fn line(columns: &[&str], alternatives: &str, info: &str, format: &str) -> Vec<u8> {
    format!(
        "chr1\t100\t.\tA\t{alternatives}\t50\tPASS\t{info}\t{format}\t{}",
        columns.join("\t")
    )
    .into_bytes()
}

fn rewrite(
    line: &[u8],
    outcomes: Vec<Outcome>,
    posteriors: Vec<f64>,
) -> (Option<String>, RewriteTally) {
    let mut tally = RewriteTally::default();
    let rewritten = rewrite_line(
        line,
        &RefinedRecord {
            background_copy: vec![0.0; outcomes.len()],
            phases: vec![None; outcomes.len()],
            outcomes,
            posteriors,
        },
        &mut GenotypePairs::default(),
        &mut tally,
    );
    (
        rewritten.map(|bytes| String::from_utf8(bytes).expect("text")),
        tally,
    )
}

fn sample_columns(text: &str) -> Vec<&str> {
    text.split('\t').skip(9).collect()
}

const BIALLELIC_FORMAT: &str = "GT:GQ:DP:AD:PL:GP";

/// **A PASS record with rows is a site; anything else passes through.**
#[test]
fn only_pass_records_with_rows_are_sites() {
    let columns = [
        "0/0:20:5:5,0:0,15,60:0.99,0.01,0",
        "./.:.:0:0,0:0,0,0:0.8,0.15,0.05",
    ];
    let info = "AF=0.25;AC=0;AN=2";
    let site = site_of(&line(&columns, "C", info, BIALLELIC_FORMAT), 2).expect("a site");
    assert_eq!(site.position, 100);
    assert_eq!(
        site.row(1),
        crate::neighbour_evidence::records::SampleRow::Flat
    );
    let filtered = String::from_utf8(line(&columns, "C", info, BIALLELIC_FORMAT))
        .expect("text")
        .replace("PASS", "hiddenDup");
    assert!(site_of(filtered.as_bytes(), 2).is_none());
    assert!(
        site_of(
            &line(&[".:.:0:.:.:.", ".:.:0:.:.:."], "C", info, BIALLELIC_FORMAT),
            2
        )
        .is_none()
    );
    assert!(site_of(&line(&columns, ".", info, BIALLELIC_FORMAT), 2).is_none());
    assert!(site_of(&line(&columns, "C", info, BIALLELIC_FORMAT), 3).is_none());
}

/// **A sample without reads, re-called, gets the pass's genotype, quality and posterior**, and
/// `AC` and `AN` count it; `PL`, `DP` and `AD` stay as they were.
#[test]
fn a_recalled_sample_gets_genotype_quality_and_posterior() {
    let source = line(
        &[
            "1/1:20:5:0,5:60,15,0:0,0.01,0.99",
            "./.:.:0:0,0:0,0,0:0.8,0.15,0.05",
        ],
        "C",
        "AF=0.5;AC=2;AN=2;DP=5",
        BIALLELIC_FORMAT,
    );
    let (text, tally) = rewrite(
        &source,
        vec![Outcome::Called, Outcome::Called],
        vec![0.0, 0.002, 0.998, 0.01, 0.09, 0.9],
    );
    let text = text.expect("kept");
    assert!(text.contains("\tAF=0.5;AC=4;AN=4;DP=5\t"), "{text}");
    let columns = sample_columns(&text);
    assert_eq!(columns[0], "1/1:27:5:0,5:60,15,0:0,0.002,0.998");
    assert_eq!(columns[1], "1/1:10:0:0,0:0,0,0:0.01,0.09,0.9");
    assert_eq!(tally.called, 2);
    assert_eq!(tally.called_without_reads, 1);
    assert_eq!(tally.genotypes_changed, 1);
}

/// **A withheld sample is `./.` with no quality and the pass's posterior**, and `AN` leaves it
/// out.
#[test]
fn a_withheld_sample_is_a_no_call_with_the_posterior() {
    let source = line(
        &[
            "1/1:20:5:0,5:60,15,0:0,0.01,0.99",
            "./.:.:0:0,0:0,0,0:0.8,0.15,0.05",
        ],
        "C",
        "AF=0.5;AC=2;AN=2",
        BIALLELIC_FORMAT,
    );
    let (text, tally) = rewrite(
        &source,
        vec![Outcome::Unchanged, Outcome::Withheld],
        vec![0.0, 0.0, 0.0, 0.6, 0.3, 0.1],
    );
    let text = text.expect("kept");
    assert!(text.contains("\tAF=0.5;AC=2;AN=2\t"), "{text}");
    assert_eq!(sample_columns(&text)[1], "./.:.:0:0,0:0,0,0:0.6,0.3,0.1");
    assert_eq!(tally.withheld, 1);
}

/// **An alternative no sample is called with any more is dropped as calling drops it** (§4.1):
/// from `ALT` and every per-allele field, every per-genotype field cut to the genotypes over the
/// alleles left — `PL` re-shifted, `GP` renormalised — and `GT` renumbered.
#[test]
fn an_uncalled_alternative_is_dropped_everywhere() {
    let source = line(
        &[
            "0/1:30:6:3,3,0:30,0,30,40,40,80:0,1,0,0,0,0",
            "0/2:10:5:2,0,3:40,30,20,0,30,50:0.1,0,0,0.8,0.1,0",
            "2/2:15:3:0,0,3:90,90,90,30,30,0:0,0,0,0,0.05,0.95",
        ],
        "C,G",
        "AF=0.3,0.2;AC=1,3;AN=6;DP=14;MQALT=60.00,55.00;MQDIFF=0.00,-5.00",
        BIALLELIC_FORMAT,
    );
    let (text, tally) = rewrite(
        &source,
        vec![Outcome::Unchanged, Outcome::Called, Outcome::Called],
        vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
            0.0, 0.0, 0.9, 0.0, 0.1, 0.0, //
            0.0, 0.5, 0.4, 0.0, 0.1, 0.0,
        ],
    );
    let text = text.expect("kept");
    let fields: Vec<&str> = text.split('\t').collect();
    assert_eq!(fields[4], "C");
    assert_eq!(fields[7], "AF=0.3;AC=4;AN=6;DP=14;MQALT=60.00;MQDIFF=0.00");
    let columns = sample_columns(&text);
    assert_eq!(columns[0], "0/1:30:6:3,3:30,0,30:0,1,0");
    assert_eq!(columns[1], "1/1:99:5:2,0:20,10,0:0,0,1");
    assert_eq!(columns[2], "0/1:4:3:0,0:0,0,0:0,0.556,0.444");
    assert_eq!(tally.alleles_dropped, 1);
}

/// **A record left with no alternative called is dropped** (§4.1, ruled §11 Q2) — but one where
/// nobody is called keeps its alleles, as calling's rule does.
#[test]
fn a_record_with_no_alternative_called_is_dropped() {
    let source = line(
        &[
            "0/1:20:5:3,2:20,0,30:0.2,0.8,0",
            "./.:.:0:0,0:0,0,0:0.8,0.15,0.05",
        ],
        "C",
        "AF=0.25;AC=1;AN=2",
        BIALLELIC_FORMAT,
    );
    let (text, tally) = rewrite(
        &source,
        vec![Outcome::Called, Outcome::Called],
        vec![0.6, 0.4, 0.0, 0.9, 0.1, 0.0],
    );
    assert!(text.is_none());
    assert_eq!(tally.records_dropped, 1);
    let (text, _) = rewrite(
        &source,
        vec![Outcome::Withheld, Outcome::Withheld],
        vec![0.6, 0.4, 0.0, 0.9, 0.1, 0.0],
    );
    let text = text.expect("nobody called: nothing dropped");
    assert!(
        text.contains("\tC\t") && text.contains("AC=0;AN=0"),
        "{text}"
    );
}

/// **A re-called repeat's `REPCN` is read off a column naming the same allele**, and is `.` for
/// an allele no input column named.
#[test]
fn repeat_counts_come_from_columns_naming_the_allele() {
    let format = "GT:GQ:DP:AD:PL:GP:REPCN";
    let source = line(
        &[
            "1/1:20:5:0,5:60,15,0:0,0.01,0.99:10,10",
            "./.:.:0:0,0:0,0,0:0.8,0.15,0.05:.",
            "0/0:20:5:5,0:0,15,60:0.99,0.01,0:9,9",
        ],
        "C",
        "AF=0.4;AC=2;AN=4",
        format,
    );
    let (text, _) = rewrite(
        &source,
        vec![Outcome::Unchanged, Outcome::Called, Outcome::Unchanged],
        vec![0.0; 3]
            .into_iter()
            .chain([0.1, 0.8, 0.1])
            .chain([0.0; 3])
            .collect(),
    );
    assert!(sample_columns(&text.expect("kept"))[1].ends_with(":9,10"));

    let source = line(
        &[
            "1/1:20:5:0,5:60,15,0:0,0.01,0.99:10,10",
            "./.:.:0:0,0:0,0,0:0.8,0.15,0.05:.",
        ],
        "C",
        "AF=0.5;AC=2;AN=2",
        format,
    );
    let (text, tally) = rewrite(
        &source,
        vec![Outcome::Unchanged, Outcome::Called],
        vec![0.0, 0.0, 0.0, 0.1, 0.8, 0.1],
    );
    assert!(sample_columns(&text.expect("kept"))[1].ends_with(":.,10"));
    assert_eq!(tally.repeat_counts_unknown, 1);
}
