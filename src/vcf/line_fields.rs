//! **A written VCF line, read back** — the fixed columns and the sample columns of one line the
//! writer produced, borrowed, for the stages that rewrite lines after they are encoded: the
//! neighbour-evidence pass and phase from reads alone. One reader, so the two judge a line by the
//! same rules.

/// The `ALT` column's index.
pub(crate) const ALT_COLUMN: usize = 4;
/// The `FILTER` column's index.
pub(crate) const FILTER_COLUMN: usize = 6;
/// The `INFO` column's index.
pub(crate) const INFO_COLUMN: usize = 7;
/// The `FORMAT` column's index.
pub(crate) const FORMAT_COLUMN: usize = 8;
/// The columns before the first sample's: `CHROM` to `FORMAT`.
pub(crate) const COLUMNS_BEFORE_THE_SAMPLES: usize = 9;

/// **One line's columns**, borrowed: the nine fixed ones split out, the sample columns kept as one
/// slice and split only as they are walked.
#[derive(Debug, Clone, Copy)]
pub(crate) struct LineFields<'a> {
    /// `CHROM` to `FORMAT`.
    pub(crate) fixed: [&'a str; COLUMNS_BEFORE_THE_SAMPLES],
    /// Every sample column, tab-separated; empty for a line without them.
    samples: &'a str,
}

impl<'a> LineFields<'a> {
    /// The line's columns, or `None` when it is not text or has fewer than the nine fixed ones.
    #[must_use]
    pub(crate) fn parse(line: &'a [u8]) -> Option<Self> {
        let text = std::str::from_utf8(line).ok()?;
        let mut fields = text.splitn(COLUMNS_BEFORE_THE_SAMPLES + 1, '\t');
        let mut fixed = [""; COLUMNS_BEFORE_THE_SAMPLES];
        for slot in &mut fixed {
            *slot = fields.next()?;
        }
        Some(Self {
            fixed,
            samples: fields.next().unwrap_or(""),
        })
    }

    /// Whether `FILTER` is `PASS`.
    #[must_use]
    pub(crate) fn is_pass(&self) -> bool {
        self.fixed[FILTER_COLUMN] == "PASS"
    }

    /// `FORMAT`'s keys, in order.
    pub(crate) fn format_keys(&self) -> impl Iterator<Item = &'a str> + use<'a> {
        self.fixed[FORMAT_COLUMN].split(':')
    }

    /// Where `key` sits in `FORMAT`, if it is there.
    #[must_use]
    pub(crate) fn format_key(&self, key: &str) -> Option<usize> {
        self.format_keys().position(|name| name == key)
    }

    /// The sample columns, in order.
    pub(crate) fn columns(&self) -> impl Iterator<Item = &'a str> + use<'a> {
        let samples = self.samples;
        (!samples.is_empty())
            .then_some(samples)
            .into_iter()
            .flat_map(|samples| samples.split('\t'))
    }

    /// How many sample columns the line has.
    #[must_use]
    pub(crate) fn column_count(&self) -> usize {
        if self.samples.is_empty() {
            0
        } else {
            self.samples.bytes().filter(|&byte| byte == b'\t').count() + 1
        }
    }
}

/// A written `GT`'s two alleles, ascending, or `None` for a no-call, a haploid `GT` or more than
/// two alleles.
#[must_use]
pub(crate) fn parse_genotype(text: &str) -> Option<(usize, usize)> {
    let mut alleles = text.split(['/', '|']).map(str::parse::<usize>);
    let first = alleles.next()?.ok()?;
    let second = alleles.next()?.ok()?;
    if alleles.next().is_some() {
        return None;
    }
    Some((first.min(second), first.max(second)))
}

#[cfg(test)]
mod tests;
