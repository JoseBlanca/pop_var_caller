//! **What a `--psp` argument list names** — the one rule three commands share.
//!
//! `call-from-psps`, `estimate-parameters` and (at plan step D1) `regenerate-census` all take the
//! same words: a psp per sample, or a directory holding them, with the flag repeated. **A run
//! types the same thing at every step, so the three had better mean the same thing by it** — and
//! the rule is not obvious enough to be safely written three times: a directory contributes the
//! psps *directly* inside it, sorted by name, and sub-directories are not descended into.
//!
//! **Sorted, because the order is the run's sample order.** It reaches the VCF's sample columns
//! and the parameters file's rows, so two runs naming one directory must open the same cohort the
//! same way whatever order the filesystem answers in.
//!
//! **Flat, because a psp's home is a `generate-psps` `--output-dir`**, which holds one psp a
//! sample and nothing else.

use std::path::PathBuf;

use crate::cli::generate_psps::PSP_FILE_EXTENSION;

/// Why a `--psp` argument could not be turned into a list of files.
///
/// **Not an error type of its own to render**: each command dresses these in its own error enum,
/// because the message a person reads names the command they typed. What is shared is the
/// question, not the wording.
#[derive(Debug)]
pub enum PspArgumentRefusal {
    /// A directory that could not be listed.
    Unlistable {
        /// The directory, as it was given.
        path: PathBuf,
        /// What the filesystem said.
        source: std::io::Error,
    },
    /// A directory holding no psp.
    Empty {
        /// The directory, as it was given.
        path: PathBuf,
    },
}

/// **The psps a command was given, with every directory expanded** — one entry a sample, in the
/// order they were named, and a directory's contents sorted by name.
///
/// # Errors
///
/// [`PspArgumentRefusal::Unlistable`] for a directory that will not list, and
/// [`PspArgumentRefusal::Empty`] for one holding no psp. **A path that is not a directory is
/// taken as given and not checked here**: whether it exists, and whether it is a psp, is the
/// cohort opener's answer, which can say what is wrong with the file rather than only that it is
/// not one.
pub fn psps_named(named: &[PathBuf]) -> Result<Vec<PathBuf>, PspArgumentRefusal> {
    let mut paths = Vec::with_capacity(named.len());
    for argument in named {
        if !argument.is_dir() {
            paths.push(argument.clone());
            continue;
        }
        let mut inside = Vec::new();
        for entry in
            std::fs::read_dir(argument).map_err(|source| PspArgumentRefusal::Unlistable {
                path: argument.clone(),
                source,
            })?
        {
            let entry = entry.map_err(|source| PspArgumentRefusal::Unlistable {
                path: argument.clone(),
                source,
            })?;
            let path = entry.path();
            if path.extension().is_some_and(|it| it == PSP_FILE_EXTENSION) && path.is_file() {
                inside.push(path);
            }
        }
        if inside.is_empty() {
            return Err(PspArgumentRefusal::Empty {
                path: argument.clone(),
            });
        }
        inside.sort();
        paths.extend(inside);
    }
    Ok(paths)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A directory holding `names`, plus one file that is not a psp.
    fn a_directory_of(names: &[&str]) -> tempfile::TempDir {
        let dir = tempfile::tempdir().expect("a temporary directory");
        for name in names {
            std::fs::write(dir.path().join(name), b"").expect("the scratch dir is ours");
        }
        std::fs::write(dir.path().join("notes.txt"), b"").expect("the scratch dir is ours");
        dir
    }

    /// **A directory contributes its psps, sorted, and nothing else.**
    #[test]
    fn a_directory_contributes_its_psps_in_name_order() {
        let dir = a_directory_of(&["zeta.psp", "alpha.psp", "mid.psp"]);

        let named = psps_named(&[dir.path().to_path_buf()]).expect("the directory lists");

        assert_eq!(
            named,
            vec![
                dir.path().join("alpha.psp"),
                dir.path().join("mid.psp"),
                dir.path().join("zeta.psp"),
            ],
        );
    }

    /// **A directory is not descended into**, because a psp's home is one `generate-psps`
    /// output directory and a nested one is another run's.
    #[test]
    fn a_sub_directory_is_not_descended_into() {
        let dir = a_directory_of(&["alpha.psp"]);
        let deeper = dir.path().join("another-run");
        std::fs::create_dir(&deeper).expect("the scratch dir is ours");
        std::fs::write(deeper.join("zeta.psp"), b"").expect("the scratch dir is ours");

        let named = psps_named(&[dir.path().to_path_buf()]).expect("the directory lists");

        assert_eq!(named, vec![dir.path().join("alpha.psp")]);
    }

    /// **A file named directly is taken as given**, whatever it is called and whether or not it
    /// exists — so a person who mistyped a path is told by the cohort opener, which can say what
    /// is wrong with the file.
    #[test]
    fn a_file_named_directly_is_taken_as_given() {
        let named =
            psps_named(&[PathBuf::from("/nowhere/zeta.pspp")]).expect("nothing is checked here");

        assert_eq!(named, vec![PathBuf::from("/nowhere/zeta.pspp")]);
    }

    /// **A directory with no psp in it is refused**, rather than contributing nothing and leaving
    /// the run to say it has no samples.
    #[test]
    fn a_directory_holding_no_psp_is_refused() {
        let dir = a_directory_of(&[]);

        let refusal = psps_named(&[dir.path().to_path_buf()]).expect_err("there is no psp in it");

        assert!(
            matches!(&refusal, PspArgumentRefusal::Empty { path } if path == dir.path()),
            "{refusal:?}",
        );
    }

    /// **Every argument contributes, in the order given** — a directory and a file together, and
    /// the directory's own contents sorted within it.
    #[test]
    fn the_arguments_are_expanded_in_the_order_they_were_given() {
        let dir = a_directory_of(&["zeta.psp", "alpha.psp"]);
        let named = psps_named(&[
            PathBuf::from("/elsewhere/first.psp"),
            dir.path().to_path_buf(),
        ])
        .expect("the directory lists");

        assert_eq!(
            named,
            vec![
                PathBuf::from("/elsewhere/first.psp"),
                dir.path().join("alpha.psp"),
                dir.path().join("zeta.psp"),
            ],
        );
    }
}
