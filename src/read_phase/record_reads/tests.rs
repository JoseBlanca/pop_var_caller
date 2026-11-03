use super::*;

fn read_allele(sample: u32, read: u64, allele: u32) -> SampleReadAllele {
    SampleReadAllele {
        sample,
        read_allele: ReadAllele {
            read,
            allele,
            error_phred: 30,
        },
    }
}

/// **The merge's alleles reach the record's numbering, and a read allele on an allele selection
/// dropped goes with it** (§2.3, step 1). Merge allele 2 was dropped; 1 became the record's 2.
#[test]
fn a_read_allele_follows_its_allele_through_selection() {
    let record = RecordReadAlleles::from_merge(
        vec![
            read_allele(3, 5, 0),
            read_allele(3, 8, 2),
            read_allele(7, 1, 1),
        ],
        |merge| match merge {
            0 => Some(0),
            1 => Some(2),
            _ => None,
        },
    );
    let flat: Vec<(u32, u64, u32)> = record
        .entries()
        .iter()
        .map(|entry| {
            (
                entry.sample,
                entry.read_allele.read,
                entry.read_allele.allele,
            )
        })
        .collect();
    assert_eq!(flat, vec![(3, 5, 0), (7, 1, 2)]);
    assert_eq!(record.of_sample(7).len(), 1);
    assert!(record.of_sample(4).is_empty());
}

/// **Renumbered with the record's alternatives** (§2.3, step 2): the record dropped its allele 1,
/// so a read allele on it goes and one on 2 becomes 1.
#[test]
fn a_read_allele_follows_its_allele_when_the_record_drops_an_alternative() {
    let mut record = RecordReadAlleles::from_merge(
        vec![
            read_allele(0, 1, 0),
            read_allele(0, 2, 1),
            read_allele(0, 3, 2),
        ],
        Some,
    );
    record.renumber(&[Some(0), None, Some(1)]);
    let flat: Vec<(u64, u32)> = record
        .entries()
        .iter()
        .map(|entry| (entry.read_allele.read, entry.read_allele.allele))
        .collect();
    assert_eq!(flat, vec![(1, 0), (3, 1)]);
}

#[test]
#[should_panic(expected = "read alleles out of order")]
fn samples_out_of_order_are_a_wiring_defect() {
    let _ = RecordReadAlleles::from_merge(vec![read_allele(4, 1, 0), read_allele(2, 1, 0)], Some);
}
