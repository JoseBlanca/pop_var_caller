use super::*;

/// **The fixed columns split out, the sample columns walked one by one**, and a line without
/// sample columns has none.
#[test]
fn a_line_reads_as_nine_fixed_columns_and_its_samples() {
    let fields = LineFields::parse(b"chr1\t5\t.\tA\tG\t50\tPASS\tAF=0.5\tGT:GQ\t0/1:40\t./.:.")
        .expect("a line");
    assert_eq!(fields.fixed[ALT_COLUMN], "G");
    assert_eq!(fields.fixed[INFO_COLUMN], "AF=0.5");
    assert!(fields.is_pass());
    assert_eq!(fields.format_key("GQ"), Some(1));
    assert_eq!(fields.format_key("PS"), None);
    assert_eq!(fields.columns().collect::<Vec<_>>(), ["0/1:40", "./.:."]);
    assert_eq!(fields.column_count(), 2);

    let sites_only = LineFields::parse(b"chr1\t5\t.\tA\tG\t50\tq10\t.\tGT").expect("a line");
    assert!(!sites_only.is_pass());
    assert_eq!(sites_only.columns().count(), 0);
    assert_eq!(sites_only.column_count(), 0);
    assert!(LineFields::parse(b"chr1\t5\t.\tA").is_none());
    assert!(LineFields::parse(b"chr1\t5\t.\t\xff\tG\t50\tPASS\t.\tGT").is_none());
}

/// **A diploid `GT` reads as its two alleles, ascending**, phased or not and with alleles of
/// several digits; a no-call, a haploid `GT` and a triploid one read as nothing.
#[test]
fn a_genotype_reads_as_two_alleles_or_nothing() {
    assert_eq!(parse_genotype("1/0"), Some((0, 1)));
    assert_eq!(parse_genotype("12|10"), Some((10, 12)));
    assert_eq!(parse_genotype("2/2"), Some((2, 2)));
    for text in ["./.", "1", "0/1/2", ".", "0/."] {
        assert_eq!(parse_genotype(text), None, "{text}");
    }
}
