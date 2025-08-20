#![allow(dead_code)]
use byte_unit::Byte;
use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};
use libc_stdhandle::stdout;
use path_absolutize::*;
use std::{
    ffi::CString,
    io::{Error, ErrorKind, Result},
    os::unix::ffi::OsStrExt,
    path::Path,
};

const VERSION: &str = include_str!(concat!(env!("OUT_DIR"), "/VERSION"));

#[derive(Clone, Debug)]
pub struct Option {
    pub cmd: String,

    pub genome: String,
    pub bam: String,
    pub mode: String,
    pub thread: usize,     //-t
    pub seed: u64,         //-s
    pub k: usize,          //-k
    pub mapq: u8,          //-q
    pub map_len: usize,    //-l
    pub map_fra: f32,      //-l
    pub batch_size: usize, //-b
    pub read_type: String, //-r

    pub depth: usize,            //-d
    pub base: usize,             //-D
    pub median_len: u32,         //-l
    pub min_len: u32,            //-m
    pub max_len: u32,            //-M
    pub noise: f64,              //-n
    pub temperature: f64,        //-e
    pub disable_chimeric: bool,  //-c
    pub disable_unaligned: bool, //-u
    pub global_error_rate: bool, //-g
}

impl Option {
    pub fn new() -> Option {
        Option::default()
    }

    pub fn from_args() -> Option {
        let opt = Option::default();
        let args = Command::new("CycSim")
			.version(VERSION)
			.about("A context-based long-read simulator")
			.arg_required_else_help(true)
			.subcommand(
				Command::new("train")
				.about("Train to extract comprehensive read features")
				.arg(
					Arg::new("bam")
						.value_name("read.bam")
						.value_parser(|x: &str| to_abspath_string(x, true))
						.required(true)
						.help("input read alignments in BAM/SAM format."),
				)
				.arg(Arg::new("genome")
					.value_name("genome.fa")
					.value_parser(|x: &str| to_abspath_string(x, true))
					.required(true)
					.help("reference genome in fasta [GZIP] format."),
				)
				.arg(
					Arg::new("read_type")
						.short('r')
						.value_name("TYPE")
						// .default_value(opt.read_type.to_string())
						.required(true)
						.value_parser(["hifi", "nanopore"])
						.help("specify read type, 'nanopore' for all electrical signal sequencing."),
				)
				.arg(
					Arg::new("thread")
						.short('t')
						.value_name("INT")
						.default_value(opt.thread.to_string())
						.value_parser(value_parser!(usize))
						.help("number of threads."),
				)
				.arg(
					Arg::new("k")
						.short('k')
						.value_name("INT")
						.default_value(opt.k.to_string())
						.value_parser(value_parser!(usize))
						.help("size of k-mers (5 < k ≤ 16)."),
				)
				.arg(
					Arg::new("mapq")
						.short('q')
						.long("mapq")
						.value_name("INT")
						.default_value(opt.mapq.to_string())
						.value_parser(value_parser!(u8))
						.help("minimum mapping quality, ignore alignments with mapping quality below this value."),
				).arg(
					Arg::new("map_len")
						.short('m')
						.long("map_len")
						.value_name("INT.FLOAT")
						.default_value((opt.map_len as f32 + opt.map_fra).to_string())
						.value_parser(value_parser!(f32))
						.help("minimum mapping length, ignore alignments shorter than max(INT, FLOAT × read_length),\ndoes not affect read length feature extraction."),
				)
				.arg(
					Arg::new("batch_size")
						.short('b')
						.long("batch_size")
						.value_name("INT[M|G]")
						.default_value(opt.batch_size.to_string())
						.value_parser(|x: &str| Byte::parse_str(x, true).map(|x| x.as_u64() as usize))
						.help("load INT bases into RAM at once for characterization."),
				)
				.arg(
					Arg::new("out")
						.short('o')
						.long("out")
						.value_name("FILE")
						.default_value("stdout")
						.value_parser(|x: &str| {
							if x != "stdout" {
								freopen_stdout(x)
							} else {
								Ok(())
							}
						})
						.help("output file."),
				)
			)
			.subcommand(
				Command::new("sim")
				.about("Simulate long reads based on a trained model")
				.arg(
					Arg::new("mode")
						.value_name("train.mode")
						.value_parser(|x: &str| to_abspath_string(x, true))
						.required(true)
						.help("input trained model file."),
				)
				.arg(Arg::new("genome")
					.value_name("genome.fa")
					.value_parser(|x: &str| to_abspath_string(x, true))
					.required(true)
					.help("reference genome in fasta [GZIP] format."),
				)
				.arg(
					Arg::new("thread")
						.short('t')
						.long("thread")
						.value_name("INT")
						.default_value(opt.thread.to_string())
						.value_parser(value_parser!(usize))
						.help("number of threads."),
				)
				.arg(
					Arg::new("seed")
						.short('s')
						.long("seed")
						.value_name("INT")
						.default_value(opt.seed.to_string())
						.value_parser(value_parser!(u64))
						.help("seed for the pseudo-random number generator (set to 0 to use a random seed).\nNote: even with the same seed, the order of simulated reads may vary due to internal multi-threading."),
				)
				.arg(
					Arg::new("depth")
						.short('d')
						.long("depth")
						.value_name("INT")
						.default_value(opt.depth.to_string())
						.value_parser(value_parser!(usize))
						.help("approximate total coverage depth to simulate."),
				)
				.arg(
					Arg::new("base")
						.short('D')
						.long("base")
						.value_name("INT")
						// .default_value(opt.base.to_string())
						.value_parser(|x: &str| Byte::parse_str(x, true).map(|x| x.as_u64() as usize))
						.help("approximate total number of bases to simulate. If not specified, initialized from `-d`."),
				).arg(
					Arg::new("median_len")
						.short('l')
						.long("median_len")
						.value_name("INT")
						// .default_value(opt.median_len.to_string())
						.value_parser(|x: &str| Byte::parse_str(x, true).map(|x| x.as_u64() as u32))
						.help("median aligned length of simulated reads. By default, inferred from training mode."),
				)
				.arg(
					Arg::new("min_len")
						.short('m')
						.long("min_len")
						.value_name("INT")
						.default_value(opt.min_len.to_string())
						.value_parser(|x: &str| Byte::parse_str(x, true).map(|x| x.as_u64() as u32))
						.help("minimum read aligned length."),
				)
				.arg(
					Arg::new("max_len")
						.short('M')
						.long("max_len")
						.value_name("INT")
						.default_value(opt.max_len.to_string())
						.value_parser(|x: &str| Byte::parse_str(x, true).map(|x| x.as_u64() as u32))
						.help("maximum read aligned length."),
				).arg(
					Arg::new("noise")
						.short('n')
						.long("noise")
						.value_name("FLOAT")
						.default_value(opt.noise.to_string())
						.value_parser(value_parser!(f64))
						.help("the proportion of random error bases introduced for a known error type."),
				)
				.arg(
					Arg::new("temperature")
						.short('e')
						.long("temperature")
						.value_name("FLOAT")
						// .default_value(opt.temperature.to_string())
						.value_parser(value_parser!(f64))
						.help("adjusts the sampling temperature (>0.), lower means more errors, higher means fewer."),
				)
				.arg(
					Arg::new("disable_chimeric")
						.short('c')
						.hide(true)
						.long("disable_chimeric")
						.help("disable simulation of chimeric reads.")
						.action(ArgAction::SetTrue),
				).arg(
					Arg::new("disable_unaligned")
						.short('u')
						.long("disable_unaligned")
						.help("disable simulation of unaligned regions at both ends of the read.")
						.action(ArgAction::SetTrue),
				)
				.arg(
					Arg::new("global_error_rate")
						.short('g')
						.long("global_error_rate")
						.help("use a global, context-independent error rate instead of modeling context-dependent error distributions.")
						.action(ArgAction::SetTrue),
				)
				.arg(
					Arg::new("out")
						.short('o')
						.long("out")
						.value_name("FILE")
						.default_value("stdout")
						.value_parser(|x: &str| {
							if x != "stdout" {
								freopen_stdout(x)
							} else {
								Ok(())
							}
						})
						.help("output file."),
				)
			)
			.subcommand(
				Command::new("rhq")
				.about("Set MAPQ to 61 for primary alignments of reads mapping equally to homozygous regions")
				.arg(
					Arg::new("bam")
						.value_name("read.bam")
						.value_parser(|x: &str| to_abspath_string(x, true))
						.required(true)
						.help("input read alignments in BAM/SAM format."),
				)
				.arg(Arg::new("genome")
					.value_name("genome.fa")
					.value_parser(|x: &str| to_abspath_string(x, true))
					.required(true)
					.help("reference genome in fasta [GZIP] format."),
				)
				.arg(
					Arg::new("thread")
						.short('t')
						.value_name("INT")
						.default_value(opt.thread.to_string())
						.value_parser(value_parser!(usize))
						.help("number of threads."),
				).arg(
					Arg::new("map_fra")
						.short('m')
						.long("map_fra")
						.value_name("FLOAT")
						.default_value(opt.map_fra.to_string())
						.value_parser(value_parser!(f32))
						.help("minimum mapping fraction, ignore alignments shorter than FLOAT × read_length."),
				)
				.arg(
					Arg::new("out")
						.short('o')
						.long("out")
						.value_name("FILE")
						.default_value("stdout")
						.value_parser(|x: &str| {
							if x != "stdout" {
								freopen_stdout(x)
							} else {
								Ok(())
							}
						})
						.help("output file."),
				)
			)
			.get_matches();
        opt.update(args)
    }

    fn update(self, args: ArgMatches) -> Option {
        match args.subcommand() {
            Some(("train", args)) => {
                let mut args = args.to_owned();
                let map_len_ = args.remove_one::<f32>("map_len").unwrap();

                Option {
                    //safely unwrap, becasue the default values have been set
                    cmd: "train".to_string(),
                    genome: args
                        .remove_one::<String>("genome")
                        .expect("Missing sequences file!"),
                    bam: args
                        .remove_one::<String>("bam")
                        .expect("Missing sequences file!"),
                    read_type: args.remove_one::<String>("read_type").unwrap(),
                    thread: args.remove_one::<usize>("thread").unwrap(),
                    k: args
                        .remove_one::<usize>("k")
                        .filter(|&k| k <= 16 && k > 5)
                        .expect("k must be between 6 and 16 (inclusive)!"),
                    mapq: args.remove_one::<u8>("mapq").unwrap(),
                    map_len: map_len_.trunc() as usize,
                    map_fra: map_len_.fract(),
                    batch_size: args.remove_one::<usize>("batch_size").unwrap(),
                    ..Default::default()
                }
            }
            Some(("sim", args)) => {
                let mut args = args.to_owned();
                let mut opt = Option {
                    //safely unwrap, becasue the default values have been set
                    cmd: "sim".to_string(),
                    genome: args
                        .remove_one::<String>("genome")
                        .expect("Missing sequences file!"),
                    mode: args
                        .remove_one::<String>("mode")
                        .expect("Missing sequences file!"),
                    thread: args.remove_one::<usize>("thread").unwrap(),
                    seed: args.remove_one::<u64>("seed").unwrap(),
                    depth: args.remove_one::<usize>("depth").unwrap(),
                    noise: args.remove_one::<f64>("noise").unwrap().max(1e-6),
                    disable_chimeric: args.get_flag("disable_chimeric"),
                    disable_unaligned: args.get_flag("disable_unaligned"),
                    global_error_rate: args.get_flag("global_error_rate"),
                    ..Default::default()
                };

                if let Some(v) = args.remove_one::<usize>("base") {
                    opt.base = v;
                }

                if let Some(v) = args.remove_one::<u32>("median_len") {
                    opt.median_len = v;
                }
                if let Some(v) = args.remove_one::<u32>("min_len") {
                    opt.min_len = v;
                }
                if let Some(v) = args.remove_one::<u32>("max_len") {
                    opt.max_len = v;
                }

                if let Some(v) = args.remove_one::<f64>("temperature") {
                    opt.temperature = v.max(0.);
                }
                opt
            }
            Some(("rhq", args)) => {
                let mut args = args.to_owned();
                Option {
                    cmd: "rhq".to_string(),
                    genome: args
                        .remove_one::<String>("genome")
                        .expect("Missing sequences file!"),
                    bam: args
                        .remove_one::<String>("bam")
                        .expect("Missing sequences file!"),
                    thread: args.remove_one::<usize>("thread").unwrap(),
                    map_fra: args.remove_one::<f32>("map_fra").unwrap(),
                    ..Default::default()
                }
            }
            _ => unreachable!("Subcommand required"),
        }
    }
}

impl Default for Option {
    fn default() -> Self {
        Option {
            cmd: String::new(),
            genome: String::new(),
            bam: String::new(),
            mode: String::new(),
            read_type: String::new(),
            seed: 0, //42,
            thread: 3,
            k: 11,
            mapq: 0,
            map_len: 2000,
            map_fra: 0.8,
            batch_size: 1000000000,
            depth: 30,
            base: 0,
            median_len: 0,
            min_len: 1000,
            max_len: 500000,
            noise: 1e-6,
            temperature: 0.,
            disable_chimeric: false,
            disable_unaligned: false,
            global_error_rate: false,
        }
    }
}

fn to_abspath_string<P: AsRef<Path>>(path: P, check_exist: bool) -> Result<String> {
    let abs_path = path.as_ref().absolutize()?.to_path_buf();

    if abs_path.exists() || !check_exist {
        Ok(abs_path.to_string_lossy().to_string())
    } else {
        Err(Error::new(
            ErrorKind::NotFound,
            format!("{abs_path:?} does not exist!"),
        ))
    }
}

fn freopen_stdout(path: &str) -> Result<()> {
    let path = Path::new(path).absolutize()?.to_path_buf();
    if path.exists() {
        return Err(Error::new(
            ErrorKind::AlreadyExists,
            format!("{path:?} already exists!"),
        ));
    }

    let mode = CString::new("w")?;
    let c_path = CString::new(path.as_os_str().as_bytes())?;

    if unsafe { libc::freopen(c_path.as_ptr(), mode.as_ptr(), stdout()) }.is_null() {
        return Err(Error::other(format!("Failed to freopen: {path:?}")));
    }
    Ok(())
}
