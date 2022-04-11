use crate::{Config, NucGrepResult};
use seal::pair::{
    Alignment, AlignmentSet, InMemoryAlignmentMatrix, NeedlemanWunsch, SmithWaterman, Step,
};
use colored::{ColoredString, Colorize};

use seq_io::fasta::{Record, RefRecord};
pub fn align_sw(record: &RefRecord, config: &Config) -> NucGrepResult<Option<String>> {
    let mut result = String::new();
    let str_x: String = record.seq().iter().map(|&c| c as char).filter(|c|!c.is_whitespace()).collect();
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

        for s in set.local_alignment().steps() {
            match s {
                Step::Align { x, y } => {
                    if sequence_x[x] == sequence_y[y] {
                        ident += 1;
                        col.push_str(&String::from(sequence_x[x]).color("green").to_string())
                    }
                },
                Step::Delete {x} => {col.push_str(&String::from(sequence_x[x]).to_string())},
                Step::Insert { y } => {}
            }
        }
        if ident >= config.needle.len() - config.allow_non_matching {
            result.push_str(&*col)//&*str_x.as_str().color("red").to_string())
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

    return Ok(Some(result));//todo ends are cut off...
}