use nucgrep::run;

fn main() {
    if let Err(e) = nucgrep::parse_args().and_then(run) {
        eprintln!("Error {}", e);
        std::process::exit(1);
    }
}
