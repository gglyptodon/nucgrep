use crate::{Config, NucGrepResult};
use bio::alignment::AlignmentOperation;
use seal::pair::{
    Alignment, AlignmentSet, InMemoryAlignmentMatrix, NeedlemanWunsch, SmithWaterman, Step,
};
//
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;
use bio::scores::blosum62;
//

use colored::{ColoredString, Colorize};

use seq_io::fasta::{Record, RefRecord};

pub fn align_new(record: &RefRecord, config: &Config) -> NucGrepResult<Option<String>> {
    let x: Vec<u8> = record
        .seq()
        .iter()
        .map(|&c| {
            if config.ignore_case {
                c.to_ascii_uppercase() as u8
            } else {
                c as u8
            }
        })
        .filter(|c| !c.to_ascii_uppercase().is_ascii_whitespace())
        .collect();

    let mut ys: Vec<Vec<u8>> = Vec::new();
    let y: Vec<u8> = config
        .needle
        .clone()
        .chars()
        .into_iter()
        .map(|c| {
            if config.ignore_case {
                c.to_ascii_uppercase() as u8
            } else {
                c as u8
            }
        })
        .collect();
    if !config.only_reverse_complement {
        ys.push(y);
    }
    if config.reverse_complement {
        let y_rev_comp: Vec<u8> = crate::reverse_complement(&config.needle, None)?
            .chars()
            .into_iter()
            .map(|c| {
                if config.ignore_case {
                    c.to_ascii_uppercase() as u8
                } else {
                    c as u8
                }
            })
            .collect();
        ys.push(y_rev_comp);
    }

    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    let scoring = Scoring::new(-5, -1, &score) // Gap open, gap extend and match score function
        .yclip(0) // Clipping penalty for x set to 'negative infinity', hence global in x
        .xclip(0); // Clipping penalty for y set to 0, hence local in y

    let mut results: Vec<String> = Vec::new();
    let mut result = String::new();
    let mut aligner = Aligner::with_scoring(scoring);
    for y in ys {
        let alignment = aligner.custom(&*x, &*y); // The custom aligner invocation
        let num_matches = alignment
            .operations
            .iter()
            .filter(|&o| o == &AlignmentOperation::Match)
            .count();
        if num_matches < config.needle.len() - config.allow_non_matching {
            results.push(String::new());
            //return Ok(Some(String::new()));
        }

        let mut x_str = &*x[..].iter().map(|u| *u as char).collect::<String>();
        result.push_str(&x_str[0..alignment.xstart]);
        let mut found = x_str[alignment.xstart..=alignment.xend].to_string();
        found = found.color("blue").to_string();
        result.push_str(&*found);
        if !alignment.xend + 1 >= x.len() {
            result.push_str(
                &*x[alignment.xend + 1..]
                    .iter()
                    .map(|u| *u as char)
                    .collect::<String>(),
            ); //todo
        }
        results.push(result.clone());
    }
    //todo merge
    let tmp = results.join("\n---todo---\n");
    Ok(Some(tmp))
}

pub fn align_sw(record: &RefRecord, config: &Config) -> NucGrepResult<Option<String>> {
    let mut result = String::new();
    let str_x: String = record
        .seq()
        .iter()
        .map(|&c| c as char)
        .filter(|c| !c.is_whitespace())
        .collect();
    let str_y = config.needle.clone();

    let strategy = SmithWaterman::new(2, -1, -1, -1);

    let sequence_x: Vec<char> = str_x.chars().collect();
    let sequence_y: Vec<char> = str_y.chars().collect();
    let set: AlignmentSet<InMemoryAlignmentMatrix> =
        AlignmentSet::new(sequence_x.len(), sequence_y.len(), strategy, |x, y| {
            sequence_x[x] == sequence_y[y]
        })
        .unwrap();
    if set.local_score() >= ((config.needle.len()) * 2 - 3 * config.allow_non_matching) as isize {
        //todo
        let mut ident = 0;
        let mut col = String::new();
        let mut last = 0;
        // let pos:Vec<Step> = set.local_alignment().steps().collect();
        // println!("pos: {:?}, last:{:?}", pos, pos.last().unwrap());
        for s in set.local_alignment().steps() {
            match s {
                Step::Align { x, y } => {
                    if sequence_x[x] == sequence_y[y] {
                        ident += 1;
                        col.push_str(&String::from(sequence_x[x]).color("blue").to_string());
                        last = x;
                    } else {
                        col.push_str(&String::from(sequence_x[x]).color("red").to_string());
                        last = x;
                    }
                }
                Step::Delete { x } => {
                    col.push_str(&String::from(sequence_x[x]).to_string());
                    last = x;
                }
                Step::Insert { y } => {}
            }
        }
        if ident >= config.needle.len() - config.allow_non_matching {
            result.push_str(&*col); //&*str_x.as_str().color("red").to_string())
            result.push_str(&*sequence_x[last + 1..].iter().collect::<String>())
            /* println!(
                "{}",
                record
                    .id()
                    .map_err(|e| eprintln!(
                        "something went wrong while getting the fasta header: {}",
                        e
                    ))
                    .unwrap()
            );*/
        }
    }

    return Ok(Some(result)); //todo ends are cut off...
}
