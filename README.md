# nucgrep
like grep, but for nucleotide sequences

```
nucgrep 0.0.1
Find sequences in sequences.
By default, look for PATTERN in each fasta record in FILE. If no file is specified, use STDIN.

         -- under construction --

USAGE:
    nucgrep [OPTIONS] --pattern <PATTERN> [FILE]

ARGS:
    <FILE>    [default: -]

OPTIONS:
    -h, --help                       Print help information
    -H, --headers-only               Only show headers for records that match
    -i, --ignore-case                Ignore case
                                     E.g. find match 'aTgA' in FILE for PATTERN 'ATGA' or 'atga' and
                                     vice versa
    -N, --allow-non-matching <N>     Maximum number of allowed non matching characters [default: 0]
    -p, --pattern <PATTERN>          PATTERN to look for in FILE, e.g. a nucleotide sequence like
                                     'AaTGATAcGGCGg'
    -r, --reverse-complement         Also show matches for the reverse complement of PATTERN
    -R, --only-reverse-complement    Show only matches for the reverse complement of PATTERN
    -V, --version                    Print version information
```