use assert_cmd::Command;
use predicates::prelude::*;
use rand::{distributions::Alphanumeric, Rng};
use std::fs;

type TestResult = Result<(), Box<dyn std::error::Error>>;

const PRG: &str = "nucgrep";
const EMPTY: &str = "tests/inputs/empty.txt";

const FASTA_INVALID: &str = "tests/inputs/invalid.fa";

const FASTA_ONE_RECORD_PER_LINE: &str = "tests/inputs/fasta_one_record_per_line.fa";
const FASTA_WITH_LINEBREAKS: &str = "tests/inputs/fasta_with_line_breaks.fa";
const FASTA_ONE_RECORD_PER_LINE_MC: &str = "tests/inputs/fasta_one_record_per_line_mc.fa"; //mixed case
const FASTA_REV_ORIG: &str = "tests/inputs/rev_orig.fa";

const FASTA_ONE_RECORD_PER_LINE_OUT: &str = "tests/expected/fasta_one_record_per_line.fa.out";
const FASTA_WITH_LINEBREAKS_OUT: &str = "tests/expected/fasta_with_line_breaks.fa.out";
const FASTA_ONE_RECORD_PER_LINE_MC_OUT: &str = "tests/expected/fasta_one_record_per_line_mc.fa.out"; //mixed case

const FASTA_REV_ORIG_OUT: &str = "tests/expected/rev_orig.fa.out";
const FASTA_REV_OUT: &str = "tests/expected/rev.fa.out";
const FASTA_ORIG_OUT: &str = "tests/expected/orig.fa.out";

const PATTERN_FOUND_ALL_UPPER: &str = "TACGGCCT";
const PATTERN_FOUND_MC: &str = "TCGTgtgCAG";
const PATTERN_FOUND_MC_ACROSS_LINEBREAK: &str = "TACGGCCT";

// --------------------------------------------------

fn run(args: &[&str], expected_file: &str) -> TestResult {
    let expected = fs::read_to_string(expected_file)?;
    Command::cargo_bin(PRG)?
        .args(args)
        .assert()
        .success()
        .stdout(expected);
    Ok(())
}

fn run_stdin(args: &[&str], input_file: &str, expected_file: &str) -> TestResult {
    let input = fs::read_to_string(input_file)?;
    let expected = fs::read_to_string(expected_file)?;
    Command::cargo_bin(PRG)?
        .write_stdin(input)
        .args(args)
        .assert()
        .success()
        .stdout(expected);
    Ok(())
}
fn gen_bad_file() -> String {
    loop {
        let filename = rand::thread_rng()
            .sample_iter(&Alphanumeric)
            .take(7)
            .map(char::from)
            .collect();

        if fs::metadata(&filename).is_err() {
            return filename;
        }
    }
}
// -- read file

/* #[test]
fn empty() -> TestResult {
    run(&[EMPTY, "--pattern", "ATG"], "tests/expected/empty.txt.out")
}*/
//todo

/*#[test]
fn simple() -> TestResult {
    run(&[EMPTY, "--pattern", "ATG"], "tests/expected/simple.fa.out")
} */
//todo

/* #[test]
fn simple_reverse_complement() -> TestResult {
    run(
        &[EMPTY, "--pattern", "CAT"],
        "tests/expected/simple.fa.reverse_complement.out",
    )
} */
//todo

/*
#[test]
fn dies_fasta_invalid_file() -> TestResult {
    let expected = format!("Error FASTA parse error: expected '>'");
    Command::cargo_bin(PRG)?
        .arg(FASTA_INVALID)
        .arg("--pattern")
        .arg("ATG")
        .arg("--fasta")
        .assert()
        .failure()
        .stderr(predicate::str::is_match(expected)?);
    Ok(())
}*/
//todo

#[test]
fn dies_bad_file() -> TestResult {
    let bad = gen_bad_file();
    let expected = format!("{}: .* [(]os error 2[)]", bad);
    Command::cargo_bin(PRG)?
        .arg(bad)
        .arg("--pattern")
        .arg("ABC")
        .assert()
        .failure()
        .stderr(predicate::str::is_match(expected)?);
    Ok(())
}

// --- stdin ---
#[test]
fn dies_fasta_invalid_stdin() -> TestResult {
    let input = fs::read_to_string(FASTA_INVALID)?;
    let expected = "Error FASTA parse error: expected '>'";
    Command::cargo_bin(PRG)?
        .write_stdin(input)
        .arg("--pattern")
        .arg("A")
        .assert()
        .failure()
        .stderr(predicate::str::is_match(expected)?);
    Ok(())
}

#[test]
fn dies_invalid_num_mismatch_args_stdin() -> TestResult {
    let input = fs::read_to_string(FASTA_ONE_RECORD_PER_LINE)?;
    let expected = format!("invalid digit found in string");
    Command::cargo_bin(PRG)?
        .write_stdin(input)
        .arg("--pattern")
        .arg("ATG")
        .arg("--allow-non-matching")
        .arg("A")
        .assert()
        .failure()
        .stderr(predicate::str::is_match(expected)?);
    Ok(())
}
#[test]
fn dies_allowed_mismatch_greater_pattern_length_args_stdin() -> TestResult {
    let input = fs::read_to_string(FASTA_ONE_RECORD_PER_LINE)?;
    let expected =
        format!("Error number of allowed mismatches is longer than or as long as pattern");
    Command::cargo_bin(PRG)?
        .write_stdin(input)
        .arg("--pattern")
        .arg("ATG")
        .arg("--allow-non-matching")
        .arg("4")
        .assert()
        .failure()
        .stderr(predicate::str::is_match(expected)?);
    Ok(())
}
#[test]
fn dies_allowed_mismatch_equal_pattern_length_args_stdin() -> TestResult {
    let input = fs::read_to_string(FASTA_ONE_RECORD_PER_LINE)?;
    let expected =
        format!("Error number of allowed mismatches is longer than or as long as pattern");
    Command::cargo_bin(PRG)?
        .write_stdin(input)
        .arg("--pattern")
        .arg("ATG")
        .arg("--allow-non-matching")
        .arg("3")
        .assert()
        .failure()
        .stderr(predicate::str::is_match(expected)?);
    Ok(())
}

#[test]
fn empty_stdin() -> TestResult {
    run_stdin(
        &["--pattern", "ATG"],
        "tests/inputs/empty.txt",
        "tests/expected/empty.txt.out",
    )
}

#[test]
fn simple_stdin() -> TestResult {
    run_stdin(
        &["--pattern", "ATG"],
        "tests/inputs/simple.fa",
        "tests/expected/simple.fa.out",
    )
}
#[test]
fn simple_reverse_complement_stdin() -> TestResult {
    run_stdin(
        &["--pattern", "CAT", "--reverse-complement"],
        "tests/inputs/simple.fa",
        "tests/expected/simple.fa_reverse_complement.out",
    )
}
#[test]
fn simple_headers_only_stdin() -> TestResult {
    run_stdin(
        &["--pattern", "ATG", "--headers-only"],
        "tests/inputs/simple.fa",
        "tests/expected/simple.fa.headers_only.out",
    )
}
#[test]
fn simple_headers_only_reverse_complement_stdin() -> TestResult {
    run_stdin(
        &["--pattern", "CAT", "--headers-only", "--reverse-complement"],
        "tests/inputs/simple.fa",
        "tests/expected/simple.fa.headers_only_reverse_complement.out",
    )
}
#[test]
fn simple_reverse_complement_only_stdin() -> TestResult {
    run_stdin(
        &["--pattern", "GTAGTAGTAG", "--reverse-complement-only"],
        FASTA_REV_ORIG,
        FASTA_REV_OUT,
    )
}
#[test]
fn simple_reverse_complement_and_orig_stdin() -> TestResult {
    run_stdin(
        &["--pattern", "GTAGTAGTAG", "-r"],
        FASTA_REV_ORIG,
        FASTA_REV_ORIG_OUT,
    )
}
#[test]
fn simple_orig_stdin() -> TestResult {
    run_stdin(&["--pattern", "GTAGTAGTAG"], FASTA_REV_ORIG, FASTA_ORIG_OUT)
}

/* #[test]
fn simple_reverse_complement() -> TestResult {
    run(
        &[EMPTY, "--pattern", "CAT"],
        "tests/expected/simple.fa.reverse_complement.out",
    )
} */
//todo

// tests line break
// test ignore case
// test 1 mismatch
// test 1 insertion
// test 1 deletion
// test multi mismatches/indels
// test fail more mismatches than pattern length
