use super::*;

fn read_allele(read: ChainId, allele: u32, error_phred: u8) -> ReadAllele {
    ReadAllele {
        read,
        allele,
        error_phred,
    }
}

/// **The Phred byte**: an error of 1 in 100 is 20, an error of one is 0, a mean too good to matter
/// stops at 60, and a mean that is not a number reads as an error of one.
#[test]
fn a_mean_log_error_becomes_its_rounded_phred() {
    assert_eq!(error_phred_of_mean_ln_error(crate::float::ln(0.01)), 20);
    assert_eq!(error_phred_of_mean_ln_error(crate::float::ln(0.05)), 13);
    assert_eq!(error_phred_of_mean_ln_error(0.0), 0);
    assert_eq!(
        error_phred_of_mean_ln_error(crate::float::ln(1e-9)),
        LARGEST_ERROR_PHRED
    );
    assert_eq!(error_phred_of_mean_ln_error(f64::NAN), 0);
    // An error the walk set to one, the losing mate of an overlapping pair (§2.5).
    assert_eq!(error_phred_of_mean_ln_error(f64::NEG_INFINITY), 60);
    assert_eq!(error_phred_of_mean_ln_error(f64::INFINITY), 0);
}

/// **One entry a read**: a read listed twice with one allele keeps the entry with the lower Phred;
/// a read listed with two alleles is removed; what is left is in read order (§2.2).
#[test]
fn a_read_counts_once_and_a_read_that_showed_two_things_counts_not_at_all() {
    let mut list = vec![
        read_allele(9, 1, 30),
        read_allele(4, 0, 35),
        read_allele(7, 0, 30),
        read_allele(4, 0, 22),
        read_allele(7, 1, 30),
        read_allele(2, 1, 40),
    ];
    settle(&mut list, 0);
    assert_eq!(
        list,
        vec![
            read_allele(2, 1, 40),
            read_allele(4, 0, 22),
            read_allele(9, 1, 30),
        ]
    );
}

/// **Above 32, the reads with the smallest hashed rank are kept, and the same reads whatever was
/// walked first** (§2.4). The second list is the first as a walk would number it after other
/// ground: every id shifted by the same amount, and the region's first id with them. A cap keyed
/// on the raw id would keep different reads.
#[test]
fn the_cap_keeps_the_same_reads_whatever_was_walked_before() {
    let alone: Vec<ReadAllele> = (0..80)
        .map(|rank| read_allele(rank, (rank % 2) as u32, 30))
        .collect();
    let shift = 1_234_567;
    let after_other_ground: Vec<ReadAllele> = alone
        .iter()
        .map(|read_allele| ReadAllele {
            read: read_allele.read + shift,
            ..*read_allele
        })
        .collect();

    let mut kept_alone = alone.clone();
    settle(&mut kept_alone, 0);
    let mut kept_after = after_other_ground;
    settle(&mut kept_after, shift);

    assert_eq!(kept_alone.len(), READ_ALLELES_KEPT);
    assert!(
        kept_alone
            .windows(2)
            .all(|pair| pair[0].read < pair[1].read)
    );
    let ranks_after: Vec<u64> = kept_after.iter().map(|kept| kept.read - shift).collect();
    let ranks_alone: Vec<u64> = kept_alone.iter().map(|kept| kept.read).collect();
    assert_eq!(ranks_alone, ranks_after);

    // The 32 kept are those with the smallest keys.
    let mut keys: Vec<u64> = alone.iter().map(|one| cap_key(one.read, 0)).collect();
    keys.sort_unstable();
    let largest_kept = kept_alone
        .iter()
        .map(|kept| cap_key(kept.read, 0))
        .max()
        .expect("32 kept");
    assert_eq!(largest_kept, keys[READ_ALLELES_KEPT - 1]);

    // And keyed on the raw id, the shifted walk would have kept other reads.
    let mut by_raw_id = alone.iter().map(|one| one.read + shift).collect::<Vec<_>>();
    by_raw_id.sort_unstable_by_key(|read| splitmix64(*read));
    by_raw_id.truncate(READ_ALLELES_KEPT);
    by_raw_id.sort_unstable();
    assert_ne!(
        by_raw_id
            .iter()
            .map(|read| read - shift)
            .collect::<Vec<_>>(),
        ranks_alone
    );
}

/// **Thirty-two or fewer are all kept**, whatever their keys.
#[test]
fn a_list_at_the_cap_is_kept_whole() {
    let mut list: Vec<ReadAllele> = (0..READ_ALLELES_KEPT as u64)
        .map(|read| read_allele(read, 0, 20))
        .collect();
    let before = list.clone();
    settle(&mut list, 0);
    assert_eq!(list, before);
}
