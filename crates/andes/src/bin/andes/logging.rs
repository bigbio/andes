//! The single stderr log format of andes.
//!
//! Every line is `[   12.3s] LEVEL stage: message`: the time since the run started, the
//! level (`ERROR`, `WARN`, `INFO`, `DEBUG`), the stage the line belongs to and the message.
//! The stage is the `target:` of the log call, or the last segment of the calling module's
//! path when no target is given. A message with several lines gets the prefix on each line.
//!
//! Library crates log through the `log` facade and never print to stderr themselves. A call
//! below the active level costs one relaxed atomic load and formats nothing. Each line is
//! formatted once and written with one `write_all` on the locked stderr, so lines from
//! different threads never interleave.

use std::fmt::Write as _;
use std::io::Write as _;
use std::sync::OnceLock;
use std::time::Instant;

use log::{Level, LevelFilter, Log, Metadata, Record};

static START: OnceLock<Instant> = OnceLock::new();
static LOGGER: StderrLogger = StderrLogger;

struct StderrLogger;

impl Log for StderrLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= log::max_level()
    }

    fn log(&self, record: &Record) {
        if !self.enabled(record.metadata()) {
            return;
        }
        let secs = START.get().map_or(0.0, |t| t.elapsed().as_secs_f64());
        let text = render(
            secs,
            record.level(),
            stage_of(record.target()),
            record.args(),
        );
        let _ = std::io::stderr().lock().write_all(text.as_bytes());
    }

    fn flush(&self) {}
}

/// The stage name of a log target: an explicit target as given, a module path by its last
/// segment (`search::refinement` → `refinement`).
fn stage_of(target: &str) -> &str {
    target.rsplit("::").next().unwrap_or(target)
}

fn level_word(level: Level) -> &'static str {
    match level {
        Level::Error => "ERROR",
        Level::Warn => "WARN",
        Level::Info => "INFO",
        Level::Debug => "DEBUG",
        Level::Trace => "TRACE",
    }
}

/// Formats one record as complete lines, the prefix repeated on each line of a
/// multi-line message.
fn render(secs: f64, level: Level, stage: &str, args: &std::fmt::Arguments) -> String {
    let message = match args.as_str() {
        Some(s) => std::borrow::Cow::Borrowed(s),
        None => std::borrow::Cow::Owned(args.to_string()),
    };
    let message = message.trim_end_matches('\n');
    let mut out = String::with_capacity(message.len() + 32);
    for line in message.split('\n') {
        let _ = writeln!(
            out,
            "[{secs:>7.1}s] {:<5} {stage}: {line}",
            level_word(level)
        );
    }
    out
}

/// The level selected by `--quiet` (warnings and errors) and `--verbose` (adds debug
/// lines); `--quiet` wins when both are given.
pub(crate) fn level_from_flags(quiet: bool, verbose: bool) -> LevelFilter {
    if quiet {
        LevelFilter::Warn
    } else if verbose {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    }
}

/// Installs the stderr logger and starts the run clock. Calling it again only changes
/// the level.
pub(crate) fn init(level: LevelFilter) {
    START.get_or_init(Instant::now);
    let _ = log::set_logger(&LOGGER);
    log::set_max_level(level);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_line_carries_time_level_and_stage() {
        let text = render(
            12.34,
            Level::Warn,
            "calibration",
            &format_args!("too few PSMs"),
        );
        assert_eq!(text, "[   12.3s] WARN  calibration: too few PSMs\n");
    }

    #[test]
    fn every_line_of_a_multi_line_message_gets_the_prefix() {
        let text = render(0.0, Level::Info, "summary", &format_args!("a\nb\n"));
        assert_eq!(
            text,
            "[    0.0s] INFO  summary: a\n[    0.0s] INFO  summary: b\n"
        );
    }

    #[test]
    fn a_module_path_target_is_named_by_its_last_segment() {
        assert_eq!(stage_of("search::refinement"), "refinement");
        assert_eq!(stage_of("calibration"), "calibration");
    }

    #[test]
    fn quiet_wins_over_verbose() {
        assert_eq!(level_from_flags(true, true), LevelFilter::Warn);
        assert_eq!(level_from_flags(false, true), LevelFilter::Debug);
        assert_eq!(level_from_flags(false, false), LevelFilter::Info);
    }
}
