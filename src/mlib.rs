#![allow(dead_code, unused_imports)]

mod logger {
    use chrono::Local;
    use ctor::ctor;
    use log::{Level, LevelFilter, Log, Metadata, Record};
    use std::io::Write;

    struct SimpleLogger;

    impl Log for SimpleLogger {
        fn enabled(&self, metadata: &Metadata) -> bool {
            metadata.level() <= Level::Info // 你可以换成 Debug
        }

        fn log(&self, record: &Record) {
            if self.enabled(record.metadata()) {
                let tid = format!("{:?}", std::thread::current().id());
                let tid_str = tid
                    .strip_prefix("ThreadId(")
                    .and_then(|s| s.strip_suffix(")"))
                    .unwrap_or(&tid);
                let timestamp = Local::now().format("%Y-%m-%d %H:%M:%S%.3f");

                let mut stderr = std::io::stderr().lock();
                let _ = writeln!(
                    stderr,
                    "[{tid} {level} {timestamp}] {msg}",
                    tid = tid_str,
                    level = record.level(),
                    timestamp = timestamp,
                    msg = record.args()
                );
            }
        }

        fn flush(&self) {}
    }

    static LOGGER: SimpleLogger = SimpleLogger;

    #[ctor]
    fn init_logger() {
        log::set_logger(&LOGGER)
            .map(|()| log::set_max_level(LevelFilter::Info))
            .expect("Failed to init logger");
    }
}

mod resource {
    use ctor::ctor;
    use lazy_static::lazy_static;
    use libc::{RUSAGE_SELF, getrusage, rusage};
    use std::{env, fmt::Write, fs, mem::MaybeUninit, time::Instant};

    lazy_static! {
        static ref START_TIME: Instant = Instant::now();
    }

    #[ctor]
    fn init_timer() {
        lazy_static::initialize(&START_TIME);
    }

    fn usage() -> rusage {
        unsafe {
            let mut usage = MaybeUninit::uninit();
            if getrusage(RUSAGE_SELF, usage.as_mut_ptr()) != 0 {
                panic!("getrusage failed");
            }
            usage.assume_init()
        }
    }

    fn cputime() -> i64 {
        let u = usage();
        u.ru_utime.tv_sec + u.ru_stime.tv_sec
    }

    fn realtime() -> u64 {
        START_TIME.elapsed().as_secs()
    }

    fn peakrss() -> i64 {
        usage().ru_maxrss
    }

    pub fn resource_str() -> String {
        let mut s = String::with_capacity(1024);
        let version_file = concat!(env!("OUT_DIR"), "/VERSION");
        if let Ok(v) = fs::read_to_string(version_file) {
            writeln!(&mut s, "Version: {v}").unwrap();
        }
        s.push_str("CMD:");
        for arg in env::args() {
            write!(&mut s, " {arg}").unwrap();
        }
        writeln!(
            &mut s,
            "\nReal time: {} sec; CPU: {} sec; Peak RSS: {:.3} GB",
            realtime(),
            cputime(),
            peakrss() as f64 / 1024.0 / 1024.0
        )
        .unwrap();
        s
    }

    //calculate interval time
    pub struct IntervalTimer {
        last: Instant,
        label: String,
        enabled: bool,
    }

    impl IntervalTimer {
        pub fn new() -> Self {
            Self::with_label("Interval")
        }

        pub fn with_label(label: impl Into<String>) -> Self {
            Self {
                last: Instant::now(),
                label: label.into(),
                enabled: true,
            }
        }

        pub fn with_label_disabled(label: impl Into<String>) -> Self {
            Self {
                last: Instant::now(),
                label: label.into(),
                enabled: false,
            }
        }

        pub fn tick(&mut self) -> f32 {
            let now = Instant::now();
            let delta = now.duration_since(self.last);
            self.last = now;
            (delta.as_secs_f64() * 1000.0).floor() as f32 / 1000.0
        }

        pub fn set_enabled(&mut self, enable: bool) {
            self.enabled = enable;
        }
    }

    impl Drop for IntervalTimer {
        fn drop(&mut self) {
            if self.enabled {
                let elapsed = self.tick();
                eprintln!("[{}] elapsed: {:.3} ms", self.label, elapsed);
            }
        }
    }
}

mod fileio {
    use std::{
        fs::{self, File},
        io::{BufReader, BufWriter},
        path::{Path, PathBuf},
    };

    pub struct FileIO {
        path: PathBuf,
    }

    impl FileIO {
        pub fn new<P: Into<PathBuf>>(path: P) -> Self {
            Self { path: path.into() }
        }

        pub fn writer(&self) -> BufWriter<File> {
            BufWriter::new(File::create(&self.path).unwrap_or_else(|e| {
                panic!("Failed to create file '{}': {}", self.path.display(), e)
            }))
        }

        pub fn reader(&self) -> BufReader<File> {
            BufReader::new(
                File::open(&self.path).unwrap_or_else(|e| {
                    panic!("Failed to read file '{}': {}", self.path.display(), e)
                }),
            )
        }

        pub fn set_done(&self) {
            let done_path = self.done_path();
            File::create(&done_path).unwrap_or_else(|e| {
                panic!(
                    "Failed to create done file '{}': {}",
                    done_path.display(),
                    e
                )
            });
        }

        pub fn check_done(&self) -> bool {
            self.done_path().exists() && self.path.exists()
        }

        pub fn remove(&self) {
            if self.path.exists() {
                fs::remove_file(&self.path).unwrap_or_else(|e| {
                    panic!("Failed to remove file '{}': {}", self.path.display(), e)
                });
            }

            let done_path = self.done_path();
            if done_path.exists() {
                fs::remove_file(&done_path).unwrap_or_else(|e| {
                    panic!(
                        "Failed to remove done file '{}': {}",
                        done_path.display(),
                        e
                    )
                });
            }
        }

        fn done_path(&self) -> PathBuf {
            let file_name = self
                .path
                .file_name()
                .expect("Invalid file path: missing file name");
            let parent = self.path.parent().unwrap_or_else(|| Path::new("."));

            parent.join(format!(".{}.done", file_name.to_string_lossy()))
        }
    }
}

mod updater {
    use reqwest::blocking::Client;

    fn latestver(owner: &str, repo: &str) -> Option<String> {
        let url = format!(
            "https://api.github.com/repos/{}/{}/releases/latest",
            owner, repo
        );

        let client = Client::builder()
            .timeout(std::time::Duration::from_secs(1))
            .build()
            .ok()?; 
        let resp = client
            .get(&url)
            .header("User-Agent", "rust-update-checker")
            .send()
            .ok()?
            .text()
            .ok()?;
        
        let key = "\"tag_name\"";
        let pos = resp.find(key)?; 
        let rest = &resp[pos + key.len()..];

        let colon_pos = rest.find(':')?;
        let rest = &rest[colon_pos + 1..];

        let first_quote = rest.find('"')?;
        let rest = &rest[first_quote + 1..];

        let second_quote = rest.find('"')?;
        let tag = &rest[..second_quote];

        Some(tag.trim_start_matches('v').to_string())
    }

    pub fn update_checker(owner: &str, repo: &str) {
        let ver = env!("CARGO_PKG_VERSION");
        if let Some(latest) = latestver(owner, repo) && latest != ver {
            eprintln!(
                "\x1b[35m Please update to the latest version: {}, current version: {} \x1b[0m",
                latest, ver
            );
        }
    }

    pub fn update_checker_default() {
        let owner = env!("GIT_OWNER");
        let repo = env!("GIT_REPO");
        if owner != "unknown" && repo != "unknown" {
            update_checker(owner, repo);
        }
    }
}

pub use fileio::FileIO;
pub use resource::{IntervalTimer, resource_str};
pub use updater::update_checker;
