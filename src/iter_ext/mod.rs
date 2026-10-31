//! Iterator adapters shared across the pipeline.
//!
//! Currently hosts `BufferedPeekable`, a peek+buffer adapter for
//! fallible iterators (`Iterator<Item = Result<T, E>>`). Future
//! cross-stage iterator helpers belong here too.
