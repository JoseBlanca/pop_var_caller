# Glossary

*The terms the caller's design documents share, in plain words. A spec that needs a term of its own
defines it where it is used; a term two specs use belongs here. Started 2026-10-08.*

## The genome, as the caller divides it

**Analysed ground.** The stretches of the reference a run was asked to look at: every contig, or the
regions a BED file names (`--regions`). A psp file records the analysed ground it was walked over.

**Typing.** The split of the analysed ground into regions by what the reference sequence is, made
before any read is seen: generic regions, STR segments, clusters of repeats and satellites
(`typed_regions.md`). It depends on the reference, the repeat catalog and the repeat criteria
(`--min-copies` and its siblings), never on a sample.

**STR segment.** A region the typing calls a short tandem repeat with clean flanks: a motif of one
to six bases repeated enough times, at enough purity, to meet the repeat criteria. In the code it is
`RegionKind::SsrSegment`; some design documents call it a *repeat tract* or a *tract*.

**Generic region.** Any analysed region that is not an STR segment, a cluster of repeats or a
satellite — ordinary sequence.

**Cluster of repeats.** Several repeats too close together for any of them to have clean flanks
(`RegionKind::SsrBundle`). Nothing calls them today.

**Satellite.** A tandem array longer than the longest STR the run types as callable (`--max-str-len`,
100 bases by default). Never called.

**Callable genome.** The part of the analysed ground the caller writes calls for: the generic
regions and the STR segments — everything except clusters of repeats and satellites. The generic
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
