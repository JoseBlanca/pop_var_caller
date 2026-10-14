# Glossary

*The terms the caller's design documents share, in plain words. A spec that needs a term of its own
defines it where it is used; a term two specs use belongs here. Started 2026-10-08.*

## The genome, as the caller divides it

**Analysed ground.** The stretches of the reference a run was asked to look at: every contig, or the
regions a BED file names (`--regions`). A psp file records the analysed ground it was walked over.

**Typing.** The split of the analysed ground into regions by what the reference sequence is, made
before any read is seen: generic regions, STR segments, STR bundles and satellites
(`typed_regions.md`). It depends on the reference, the repeat catalog and the repeat criteria
(`--min-copies` and its siblings), never on a sample.

**STR segment.** A region the typing calls a short tandem repeat with clean flanks: a motif of one
to six bases repeated enough times, at enough purity, to meet the repeat criteria. In the code it is
`RegionKind::SsrSegment`; some design documents call it a *repeat tract* or a *tract*.

**Generic region.** Any analysed region that is not an STR segment, an STR bundle or a
satellite — ordinary sequence.

**STR bundle.** Two or more repeats within 15 bases of each other, too close for any of them to have
clean flanks (`RegionKind::SsrBundle`). Its region is the hull of its repeats, first start to
furthest end (repeats can nest and overlap). Each repeat is a **member**, carrying its span cut back to whole copies of its motif, and the
motif. A bundle longer than `--max-str-bundle-len` (100 bases by default) is typed as a satellite.
Nothing calls them yet; `str_bundle_loci.md` designs an STR locus per bundle. Earlier documents
call it a *cluster of repeats*.

**Satellite.** A tandem array longer than the longest STR the run types as callable (`--max-str-len`,
100 bases by default), or an STR bundle longer than `--max-str-bundle-len`. Never called.

**Callable genome.** The part of the analysed ground the caller writes calls for: the generic
regions and the STR segments — everything except STR bundles and satellites. The generic
path can walk all of it; the STR path walks the STR segments.

## The two calling paths

**Generic path.** The SNP/indel caller: the pileup walk, which records every position a read covers,
and the calling of the loci built from those records. By default it walks the generic regions and
stops at every STR segment.

**STR path.** The caller of STR segments: a walk that records each segment's reads as whole-segment
spellings, and a read model that allows for stutter — a copy of the motif gained or lost while the
library was made.

**Calling both paths over STR segments.** *Experimental, behind `--tract-calling both`.* The generic
path calls the whole callable genome, STR segments included, and the STR path calls the STR
segments in parallel; at the end the generic call is kept on generic regions, and on each STR
segment the two calls are compared and combined (`tract_dual_path.md`). With the default,
`tract-only`, each STR segment is called by the STR path alone.
