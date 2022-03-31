extern crate core;

use clap::{Arg, Command};
use colored::{ColoredString, Colorize};
use seq_io::fasta::{Reader, Record};
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io;
use std::io::BufRead;

pub type NucGrepResult<T> = Result<T, Box<dyn Error>>;

#[derive(Debug)]
pub struct Config {
    fasta: bool,
    needle: String, // string to be found in file
    file: String,   // "-" is STDIN
    reverse_complement: bool,
    only_reverse_complement: bool,
    allow_non_matching: usize,
    ignore_case: bool,
    headers_only: bool,
}

pub fn search(buff: impl BufRead) -> NucGrepResult<()> {
    Ok(())
}

pub fn search_fasta_stdin(
    needle: &String,
    headers_only: bool,
    search_reverse_complement: bool,
    only_reverse_complement: bool,
) -> NucGrepResult<()> {
    let mut reader = Reader::new(io::stdin());
    let reverse_complement_needle: Option<String> = match search_reverse_complement {
        true => Some(reverse_complement(&needle, None)?),
        false => None,
    };
    while let Some(result) = reader.next() {
        let tmp = result?.clone();
        let fullseq = tmp.seq().iter().map(|&x| x as char).collect::<String>();
        if !only_reverse_complement {
            if fullseq.contains(&*needle) {
                println!(
                    ">{}{}",
                    green(tmp.id().expect("There was a problem with the sequence ID")),
                    if !headers_only {
                        format!("\n{}", fullseq)
                    } else {
                        "".to_string()
                    }
                );
            }
        }
        match reverse_complement_needle {
            //also look for revcomp
            Some(ref p) => {
                if fullseq.contains(&*p) {
                    println!(
                        ">{}{}",
                        green(tmp.id().expect("There was a problem with the sequence ID")),
                        if !headers_only {
                            format!("\n{}", fullseq)
                        } else {
                            "".to_string()
                        }
                    );
                }
            }
            None => {}
        }
    }
    Ok(())
}

pub fn green(s: &str) -> ColoredString {
    format!("{}", s).green().bold()
}

pub fn run(config: Config) -> NucGrepResult<()> {
    //println!("{:#?}", config);
    let mut search_patterns: Vec<String> = Vec::new();
    if config.reverse_complement {
        search_patterns.push(reverse_complement(&config.needle, None)?);
    }
    if !config.only_reverse_complement {
        search_patterns.push(config.needle.clone());
    }
    if config.file == "-" {
        search_fasta_stdin(
            &config.needle,
            config.headers_only,
            config.reverse_complement,
            config.only_reverse_complement,
        )?; //todo refactor
    } else {
        let reader = open(&config.file);
        match reader {
            Err(e) => {
                eprintln!("{}: {}", config.file, e);
                std::process::exit(1)
            }
            Ok(file) => search(file)?,
        }
    }
    Ok(())
}

pub fn parse_args() -> NucGrepResult<Config> {
    let matches = Command::new("nucgrep")
        .version("0.0.1")
        .about(
            "Find sequences in sequences.\n\
   By default, look for PATTERN in each fasta record in FILE. If no file is specified, use STDIN.\n\n\t -- under construction --",
        )
        .arg(
            Arg::new("needle")
                .short('p')
                .long("--pattern")
                .help(
                    "PATTERN to look for in FILE, e.g. a nucleotide sequence like 'AaTGATAcGGCGg'",
                )
                .required(true)
                .takes_value(true)
                .value_name("PATTERN"),
        )
        .arg(
            Arg::new("path")
                .allow_invalid_utf8(true)
                .default_value("-")
                .value_name("FILE"),
        )
        .arg(
            Arg::new("reverse_complement")
                .short('r')
                .long("--reverse-complement")
                .takes_value(false)
                .help("Also show matches for the reverse complement of PATTERN"),
        )
        .arg(
            Arg::new("allow_non_matching")
                .long("allow-non-matching")
                .short('N')
                .takes_value(true)
                .required(false)
                .value_name("N")
                .default_value("0")
                .help("Maximum number of allowed non matching characters"),
        )
        .arg(
            Arg::new("only_reverse_complement")
                .long("reverse-complement-only")
                .short('R')
                .help("Show only matches for the reverse complement of PATTERN")
                .takes_value(false),
        )
        .arg(
            Arg::new("ignore_case")
                .long("--ignore-case")
                .short('i')
                .help("Ignore case\nE.g. find match 'aTgA' in FILE for PATTERN 'ATGA' or 'atga' and vice versa"),
        )
        .arg(Arg::new("headers_only").long("--headers-only").short('H').help("Only show headers for records that match").takes_value(false))
        .get_matches();
    Ok(Config {
        fasta: true, //matches.is_present("fasta"),
        needle: matches.value_of("needle").map(String::from).unwrap(),
        file: matches.value_of_lossy("path").map(String::from).unwrap(),
        reverse_complement: matches.is_present("reverse_complement")
            || matches.is_present("only_reverse_complement"),
        only_reverse_complement: matches.is_present("only_reverse_complement"),
        allow_non_matching: matches
            .value_of("allow_non_matching")
            .unwrap()
            .parse::<usize>()
            .map_err(|e| {
                eprintln!(
                    "{}\n{}",
                    e, "\nHint: --allow-non-matching takes a number as argument"
                );
                std::process::exit(2)
            })
            .unwrap(),
        ignore_case: matches.is_present("ignore_case"),
        headers_only: matches.is_present("headers_only"),
    })
}

pub fn open(path: &String) -> NucGrepResult<Box<dyn BufRead>> {
    if path == "-" {
        Ok(Box::new(std::io::BufReader::new(std::io::stdin())))
    } else {
        Ok(Box::new(std::io::BufReader::new(File::open(path)?)))
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct NucleotideComplementError;

impl Display for NucleotideComplementError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "invalid nucleotide")
    }
}

impl Error for NucleotideComplementError {}

pub fn reverse_complement(
    s: &str,
    is_dna: Option<bool>,
) -> Result<String, NucleotideComplementError> {
    let is_dna_or_guess = match is_dna {
        None => {
            if s.find('u').is_some() || s.find('U').is_some() {
                if s.find('t').is_some() || s.find('T').is_some() {
                    return Err(NucleotideComplementError); // u and t in one sequence
                } else {
                    false
                }
            } else {
                true
            }
        }
        Some(b) => b,
    };
    let complement = s
        .chars()
        .map(|c| match c {
            'A' => {
                if is_dna_or_guess {
                    'T'
                } else {
                    'U'
                }
            }
            'a' => {
                if is_dna_or_guess {
                    't'
                } else {
                    'u'
                }
            }
            'G' => 'C',
            'g' => 'c',
            'T' => 'A',
            't' => 'a',
            'U' => 'A',
            'u' => 'a',
            'C' => 'G',
            'c' => 'g',
            'N' => 'N',
            'n' => 'n',
            'W' => 'W',
            'w' => 'w',
            'S' => 'S',
            's' => 's',
            'M' => 'K',
            'm' => 'k',
            'K' => 'M',
            'k' => 'm',
            'R' => 'Y',
            'r' => 'y',
            'Y' => 'R',
            'y' => 'r',
            'B' => 'V',
            'b' => 'v',
            'D' => 'H',
            'd' => 'h',
            'H' => 'D',
            'h' => 'd',
            'V' => 'B',
            'v' => 'b',
            '-' => '-',
            _ => '?',
        })
        .collect::<String>();
    if complement.contains('?') {
        Err(NucleotideComplementError)
    } else {
        let rev: String = complement.chars().rev().collect();

        Ok(rev)
    }
}

mod tests {
    use crate::{reverse_complement, NucleotideComplementError};

    #[test]
    pub fn test_revcomp() {
        let input = "ATG";
        let expected = "CAT";
        assert_eq!(reverse_complement(&input, Some(true)).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_lower() {
        let input = "atgcn";
        let expected = "ngcat";
        assert_eq!(reverse_complement(&input, Some(true)).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_lower_rna() {
        let input = "augcn";
        let expected = "ngcau";
        assert_eq!(reverse_complement(&input, Some(false)).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_upper() {
        let input = "GCTA";
        let expected = "TAGC";
        assert_eq!(reverse_complement(&input, Some(true)).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_upper_rna() {
        let input = "GCUA";
        let expected = "UAGC";
        assert_eq!(reverse_complement(&input, None).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_mixed_rna() {
        let input = "GcuA";
        let expected = "UagC";
        assert_eq!(reverse_complement(&input, None).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_guess() {
        let input = "GGGGaaaaaaaatttatatat";
        let expected = "atatataaattttttttCCCC";
        assert_eq!(reverse_complement(&input, None).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_guess2() {
        let input = "GGGGaaaaaaaauuuauauau";
        let expected = "auauauaaauuuuuuuuCCCC";
        assert_eq!(reverse_complement(&input, None).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_rna() {
        let input = "AUG";
        let expected = "CAU";
        assert_eq!(reverse_complement(&input, Some(false)).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_lower_iupac() {
        let input = "vnhdbatgcywsmkrd";
        let expected = "hymkswrgcatvhdnb";
        assert_eq!(reverse_complement(&input, None).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_upper_iupac() {
        let input = "ATGCNHVBYWSMKRD";
        let expected = "HYMKSWRVBDNGCAT";
        assert_eq!(reverse_complement(&input, None).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_mixed_case_iupac() {
        let input = "ATGCNhVBYWsMKRdD";
        let expected = "HhYMKsWRVBdNGCAT";
        assert_eq!(reverse_complement(&input, None).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_fail() {
        let input = "blarghblergh!";
        let expected = Err(NucleotideComplementError);
        assert_eq!(reverse_complement(&input, None), expected);
    }
    #[test]
    pub fn test_revcomp_fail_rna_dna() {
        let input = "atgcu";
        let expected = Err(NucleotideComplementError);
        assert_eq!(reverse_complement(&input, None), expected);
    }
    #[test]
    pub fn test_revcomp_default_to_dna() {
        let input = "aAaAAAA";
        let expected = "TTTTtTt";
        assert_eq!(reverse_complement(&input, None).unwrap(), expected);
    }
    #[test]
    pub fn test_revcomp_force_to_rna() {
        let input = "aAaAAAA";
        let expected = "UUUUuUu";
        assert_eq!(reverse_complement(&input, Some(false)).unwrap(), expected);
    }
}
