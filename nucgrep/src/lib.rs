use clap::{Arg, Command};
use std::error::Error;
use std::fs::File;
use std::io::BufRead;

pub type NucGrepResult<T> = Result<T, Box<dyn Error>>;

#[derive(Debug)]
pub struct Config {
    fasta: bool,
    needle: String, // string to be found in file
    file: String,   // "-" is STDIN
    reverse_complement: bool,
}

pub fn run(config: Config) -> NucGrepResult<()> {
    println!("{:#?}", config);
    let mut reader = open(config.file)?;
    let mut buff = String::new();
    while let Ok(bytes_read) = reader.read_line(&mut buff) {
        println!("{} {}", bytes_read, buff);
        buff.clear();
        if bytes_read == 0 {
            break;
        }
    }
    Ok(())
}

pub fn parse_args() -> NucGrepResult<Config> {
    let matches = Command::new("nucgrep")
        .version("0.1.0")
        .about(
            "Find sequences in sequences.\n\
    By default, look for PATTERN in each line of FILE. If no file is specified, use STDIN.\n",
        )
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .takes_value(false)
                .help("Use fasta format records instead of lines"),
        )
        .arg(
            Arg::new("needle")
                .short('p')
                .long("--pattern")
                .help("PATTERN to look for in FILE, e.g. a nucleotide sequence like 'AaTGATAcGGCGg'")
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
                .help("Maximum number of allowed non matching characters"),
        )
        .arg(
            Arg::new("only_reverse_complement")
                .long("only-reverse-complement")
                .short('R')
                .help("Show only matches for the reverse complement of PATTERN")
                .takes_value(false),
        )
        .get_matches();
    Ok(Config {
        fasta: matches.is_present("fasta"),
        needle: matches.value_of("needle").map(String::from).unwrap(),
        file: matches.value_of_lossy("path").map(String::from).unwrap(),
        reverse_complement: matches.is_present("reverse_complement"),
    })
}

pub fn open(path: String) -> NucGrepResult<Box<dyn BufRead>> {
    if path == "-" {
        Ok(Box::new(std::io::BufReader::new(std::io::stdin())))
    } else {
        Ok(Box::new(std::io::BufReader::new(File::open(path)?)))
    }
}
