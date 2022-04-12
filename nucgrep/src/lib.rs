pub mod search_fasta;

extern crate core;

use clap::{Arg, Command};
use colored::{ColoredString, Colorize};
use regex::{Regex, RegexBuilder};
use seq_io::fasta::{Reader, Record, RefRecord};
use std::collections::HashSet;
use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::fs::{read, File};
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
    linewrap: Option<usize>,
}

#[derive(Debug)]
struct NucGrepMatch {
    start: usize,
    end: usize,
}
pub fn search_fasta<T: std::io::Read>(
    config: &Config,
    needle: &Regex,
    mut reader: Reader<T>,
) -> NucGrepResult<()> {
    while let Some(result) = reader.next() {
        let tmp = result?.clone();
        let fullseq = tmp
            .seq()
            .iter()
            .map(|&x| x as char)
            .filter(|c| !c.is_whitespace())
            .collect::<String>();

        //refactor
        if config.allow_non_matching == 0 {
            let res = search_nonfuzzy(&fullseq, &needle, &config);
            if res.is_some() {
                if config.headers_only {
                    println!(">{}", tmp.id()?)
                } else if config.linewrap.is_some() {
                    println!(
                        ">{}\n{}",
                        tmp.id()?,
                        textwrap::fill(&*res.unwrap(), config.linewrap.unwrap())
                    );
                } else {
                    println!("{}{}", format!(">{}\n", tmp.id()?), res.unwrap())
                }
            }
        }
        // fuzzy matching  // under construction
        else {
            let res = crate::search_fasta::align_new(&tmp, &config);

            if res.is_ok() {
                if let Some(result) = res? {
                    if config.headers_only {
                        if !result.is_empty() {
                            println!(">{}", tmp.id()?)
                        }
                    } else if config.linewrap.is_some() && !result.is_empty() {
                        println!(
                            ">{}\n{}",
                            tmp.id()?,
                            textwrap::fill(&*result, config.linewrap.unwrap())
                        );
                    } else {
                        if !result.is_empty() {
                            println!("{}{}", format!(">{}\n", tmp.id()?), result)
                        }
                    }
                }
            }
        }
    }
    Ok(())
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

    let pattern = match config.reverse_complement {
        true => {
            format!(
                r"{}{}",
                reverse_complement(&config.needle, None)?,
                if !config.only_reverse_complement {
                    format!("|{}", config.needle)
                } else {
                    "".to_string()
                }
            )
        }
        false => config.needle.clone(),
    };
    let needle_regex = match config.ignore_case {
        true => RegexBuilder::new(&*pattern)
            .case_insensitive(true)
            .build()
            .expect("Invalid Regex"),
        false => RegexBuilder::new(&*pattern)
            .case_insensitive(false)
            .build()
            .expect("Invalid Regex"),
    };

    if config.file == "-" {
        let reader = Reader::new(std::io::stdin());
        search_fasta(&config, &needle_regex, reader)?
    } else {
        let reader = Reader::from_path(&config.file);
        match reader {
            Err(e) => {
                eprintln!("{}: {}", config.file, e);
                std::process::exit(1)
            }
            Ok(reader) => search_fasta(&config, &needle_regex, reader)?,
        }
    }
    Ok(())
}

#[derive(Debug)]
struct InvalidValueNonMatchingAllowedLargerEqualPattern;

impl Display for InvalidValueNonMatchingAllowedLargerEqualPattern {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "number of allowed mismatches is longer than or as long as pattern"
        )
    }
}

impl Error for InvalidValueNonMatchingAllowedLargerEqualPattern {}

pub fn parse_args() -> NucGrepResult<Config> {
    let matches = Command::new("nucgrep")
        .version("0.0.1")
        .about(
            "Find sequences in sequences.\nBy default, look for PATTERN in each fasta record in FILE. If no file is specified, use STDIN.\n\n\t -- under construction --",
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
                .help("Maximum number of allowed non-matching characters (Caution: experimental & slow!)"),
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
        .arg(Arg::new("wrap_lines").long("--wrap-lines").short('w').help("Wrap output every C characters").takes_value(true).value_name("C").required(false))
        .get_matches();

    let number_non_matching_allowed = matches
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
        .unwrap();
    let pattern_length = matches.value_of("needle").map(String::from).unwrap().len();
    if pattern_length <= number_non_matching_allowed {
        return Err(Box::new(InvalidValueNonMatchingAllowedLargerEqualPattern));
    }
    let wrap = match matches.value_of("wrap_lines") {
        None => None,
        Some(num) => Some(
            num.parse::<usize>()
                .map_err(|e| {
                    eprintln!(
                        "{}\n{}",
                        e, "\nHint: --wrap-lines takes a number as argument"
                    );
                    std::process::exit(2)
                })
                .unwrap(),
        ),
    };
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
        linewrap: wrap,
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

pub fn highlight_match(haystack: &str, needle: &str) -> NucGrepResult<String> {
    //let result = String::new();
    let needle_regex = RegexBuilder::new(needle).case_insensitive(true).build()?;
    let mut display_seq = haystack.clone();

    let mut nmatches: Vec<NucGrepMatch> = Vec::new();

    let mut found = HashSet::new();

    for i in needle_regex.find_iter(&haystack) {
        nmatches.push(NucGrepMatch {
            start: i.start(),
            end: i.end(),
        });
        found.insert(String::from(i.as_str()));
    }
    let mut offset: usize = 0;
    let mut result = String::from("");
    for m in &nmatches {
        let l = m.end - m.start;
        result.push_str(&*format!(
            "{}{}",
            &display_seq[offset..m.start],
            if &haystack[m.start..m.end].to_uppercase() == &needle.to_uppercase() {
                haystack[m.start..m.end].color("blue").bold()
            } else {
                haystack[m.start..m.end].color("yellow").bold()
            }
        ));
        offset = m.start + l;
    }
    result.push_str(&haystack[offset..]);

    Ok(result)
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

pub fn search_fuzzy(record: &RefRecord, config: &Config) -> NucGrepResult<Option<String>> {
    let windows_size = config.needle.len() - config.allow_non_matching; //todo: ending? //indels

    let fullseq = record
        .seq()
        .iter()
        .map(|&x| x as char)
        .filter(|c| !c.is_whitespace())
        .collect::<String>(); //rm whitespace

    let revcomp_needle = Some(reverse_complement(&*config.needle, None)?);
    let searchpattern_forward = config.needle.chars().map(|c| c as u8).collect::<Vec<u8>>();
    let searchpattern_revcomp = revcomp_needle
        .unwrap()
        .chars()
        .map(|c| c as u8)
        .collect::<Vec<u8>>();
    //-----
    let mut was_found = false;
    let mut foundpatterns: Vec<String> = Vec::new();
    let mut result = String::new();

    //only forward
    if !config.reverse_complement {
        for w in fullseq //mismatches
            .chars()
            .map(|c| c as u8)
            .collect::<Vec<u8>>()
            .windows(windows_size)
        {
            //todo
            if generic_levenshtein::distance(w, &searchpattern_forward[..])
                <= config.allow_non_matching
            {
                was_found = true;
                foundpatterns.push(w.iter().map(|&c| c as char).collect::<String>());
            }
        }
        if was_found {
            let assembled = foundpatterns.join("|");
            result.push_str(&*format!("{}", highlight_match(&fullseq, &assembled)?));
        } //todo
    }
    //only reverse complement
    else if config.only_reverse_complement {
        for w in fullseq
            .chars()
            .map(|c| c as u8)
            .collect::<Vec<u8>>()
            .windows(windows_size)
        {
            //todo
            if generic_levenshtein::distance(w, &searchpattern_revcomp[..])
                <= config.allow_non_matching
            {
                was_found = true;
                foundpatterns.push(w.iter().map(|&c| c as char).collect::<String>());
            }
        }
        if was_found {
            let assembled = foundpatterns.join("|");
            result.push_str(&*format!("{}", highlight_match(&fullseq, &assembled)?));
        } //todo
    }
    //both
    else {
        if config.reverse_complement {
            foundpatterns = Vec::new();
            result = String::new();

            // search for both
            for w in fullseq
                .chars()
                .map(|c| c as u8)
                .collect::<Vec<u8>>()
                .windows(windows_size)
            {
                //todo
                if generic_levenshtein::distance(w, &searchpattern_forward[..])
                    <= config.allow_non_matching
                {
                    was_found = true;
                    foundpatterns.push(w.iter().map(|&c| c as char).collect::<String>());
                }
                if generic_levenshtein::distance(w, &searchpattern_revcomp[..])
                    <= config.allow_non_matching
                {
                    was_found = true;
                    foundpatterns.push(w.iter().map(|&c| c as char).collect::<String>());
                }
            }
            if was_found {
                let assembled = foundpatterns.join("|");
                result.push_str(&*format!("{}", highlight_match(&fullseq, &assembled)?));
            } //todo
        }
    }

    Ok(Some(result))
}

pub fn search_nonfuzzy(haystack: &String, needle: &Regex, config: &Config) -> Option<String> {
    let mut result = String::new();
    let mut display_seq = haystack.clone();
    let mut nmatches: Vec<NucGrepMatch> = Vec::new();
    let mut found = HashSet::new();

    for i in needle.find_iter(&*haystack) {
        nmatches.push(NucGrepMatch {
            start: i.start(),
            end: i.end(),
        });
        found.insert(String::from(i.as_str()));
    }
    if found.is_empty() {
        return None;
    }

    let mut offset: usize = 0;
    for m in &nmatches {
        let l = m.end - m.start;
        result.push_str(&*format!(
            "{}{}",
            &display_seq[offset..m.start],
            if &haystack[m.start..m.end].to_uppercase() == &config.needle.to_uppercase() {
                haystack[m.start..m.end].color("green").bold()
            } else {
                haystack[m.start..m.end].color("purple").bold()
            }
        ));
        offset = m.start + l;
    }
    result.push_str(&haystack[offset..]);

    Some(result)
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

    #[test]
    pub fn test_levenshtein_one_mismatch() {
        let maxdist: usize = 3;
        let expected = true;
        let mut result = false;
        let input = "atgcta";
        let needle = "tcta";
        let n: Vec<char> = needle.chars().collect();
        let tmp = input.chars().collect::<Vec<char>>();
        for t in tmp.windows(needle.len() - maxdist) {
            let dist = generic_levenshtein::distance(t, &n);
            if dist <= maxdist {
                println!("{:?} matches {:?}", &t, &n);
                result = true;
            }
        }
        assert_eq!(result, expected);
    }

    #[test] //todo
    pub fn test_levenshtein_one_deletion_in_haystack_short() {
        let maxdist: usize = 1;
        let expected = true;
        let mut result = false;
        let input = "atgta";
        let needle = "gcta";
        let n: Vec<char> = needle.chars().collect();
        let tmp = input.chars().collect::<Vec<char>>();
        for t in tmp.windows(needle.len() - maxdist) {
            let dist = generic_levenshtein::distance(t, &n);
            if dist <= maxdist {
                println!("{:?} matches {:?}", &t, &n);
                result = true;
            }
        }
        assert_eq!(result, expected);
    }

    #[test] //todo
    pub fn test_levenshtein_one_deletion_in_haystack_short_dist_3() {
        let maxdist: usize = 3;
        let expected = true;
        let mut result = false;
        let input = "atgtaaaaaaatangagtttagactagnatag";
        let needle = "atgtaaaaaaaaaat";
        let n: Vec<char> = needle.chars().collect();
        let tmp = input.chars().collect::<Vec<char>>();
        for t in tmp.windows(needle.len() - maxdist) {
            let dist = generic_levenshtein::distance(t, &n);
            if dist <= maxdist {
                println!("{:?} matches {:?}", &t, &n);
                result = true;
            }
        }
        assert_eq!(result, expected);
    }

    #[test]
    pub fn test_levenshtein_one_insertion_in_haystack() {
        let maxdist: usize = 3;
        let expected = true;
        let mut result = false;
        let input = "atgctaatgatgatg";
        let needle = "tcta";
        let n: Vec<char> = needle.chars().collect();
        let tmp = input.chars().collect::<Vec<char>>();
        for t in tmp.windows(needle.len() - maxdist) {
            let dist = generic_levenshtein::distance(t, &n);
            println!("{:?} cmp {:?} -> {:?}", &t, &n, dist);
            if dist <= maxdist {
                println!("{:?} matches {:?}", &t, &n);
                result = true;
            }
        }
        assert_eq!(result, expected);
    }
    //todo more tests...
    #[test]
    pub fn test_levenshtein_insertion_in_haystack_1() {
        let maxdist: usize = 1;
        let expected = true;
        let mut result = false;
        let input = "actgt";
        let needle = "acgt";
        let n: Vec<char> = needle.chars().collect();
        let tmp = input.chars().collect::<Vec<char>>();
        for t in tmp.windows(needle.len() - maxdist) {
            //deletion in haystack
            let dist = generic_levenshtein::distance(t, &n);
            println!("{:?} cmp {:?} -> {:?}", &t, &n, dist);
            if dist <= maxdist {
                println!("{:?} matches {:?}", &t, &n);
                result = true;
            }
        }
        for t in tmp.windows(needle.len() + maxdist) {
            // insertion in haystack
            let dist = generic_levenshtein::distance(t, &n);
            println!("{:?} cmp {:?} -> {:?}", &t, &n, dist);
            if dist <= maxdist {
                println!("{:?} matches {:?}", &t, &n);
                result = true;
            }
        }
        for t in tmp.windows(needle.len()) {
            // only mismatches
            let dist = generic_levenshtein::distance(t, &n);
            println!("{:?} cmp {:?} -> {:?}", &t, &n, dist);
            if dist <= maxdist {
                println!("{:?} matches {:?}", &t, &n);
                result = true;
            }
        }
        assert_eq!(result, expected);
    }

    #[test]
    pub fn test_levenshtein_insertion_in_haystack_2() {
        let maxdist: usize = 2;
        let expected = true;
        let mut result = false;
        let input = "acnntgt";
        let needle = "acgt";
        let n: Vec<char> = needle.chars().collect();
        let tmp = input.chars().collect::<Vec<char>>();
        for t in tmp.windows(needle.len() - maxdist) {
            //deletion in haystack
            let dist = generic_levenshtein::distance(t, &n);
            println!("{:?} cmp {:?} -> {:?}", &t, &n, dist);
            if dist <= maxdist {
                println!("{:?} matches {:?}", &t, &n);
                result = true;
            }
        }
        for t in tmp.windows(needle.len() + maxdist) {
            // insertion in haystack
            let dist = generic_levenshtein::distance(t, &n);
            println!("{:?} cmp {:?} -> {:?}", &t, &n, dist);
            if dist <= maxdist {
                println!("{:?} matches {:?}", &t, &n);
                result = true;
            }
        }
        for t in tmp.windows(needle.len()) {
            // only mismatches
            let dist = generic_levenshtein::distance(t, &n);
            println!("{:?} cmp {:?} -> {:?}", &t, &n, dist);
            if dist <= maxdist {
                println!("{:?} matches {:?}", &t, &n);
                result = true;
            }
        }
        assert_eq!(result, expected);
    }

    /*
        #[test]
        pub fn test_levenshtein_insertion_in_haystack_3() {
            let maxdist: usize = 3;
            let expected = true;
            let mut result = false;
            let input =  "acnnggtacgt";
            let needle = "acgtacgt";
            let n: Vec<char> = needle.chars().collect();
            let tmp = input.chars().collect::<Vec<char>>();
            for t in tmp.windows(needle.len()-maxdist) { //deletion in haystack
                let dist = generic_levenshtein::distance(t, &n);
                println!("{:?} cmp {:?} -> {:?}", &t, &n, dist);
                if dist <= maxdist {
                    println!("{:?} matches {:?}", &t, &n);
                    result = true;
                }
            }
            for t in tmp.windows(needle.len()+maxdist) { // insertion in haystack
                let dist = generic_levenshtein::distance(t, &n);
                println!("{:?} cmp {:?} -> {:?}", &t, &n, dist);
                if dist <= maxdist {
                    println!("{:?} matches {:?}", &t, &n);
                    result = true;
                }
            }
            for t in tmp.windows(needle.len()) { // only mismatches
                let dist = generic_levenshtein::distance(t, &n);
                println!("{:?} cmp {:?} -> {:?}", &t, &n, dist);
                if dist <= maxdist {
                    println!("{:?} matches {:?}", &t, &n);
                    result = true;
                }
            }
            assert_eq!(result, expected);
        }


    */
}
//todo pretty print
