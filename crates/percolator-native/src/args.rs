use std::path::PathBuf;

pub struct Args {
    pub seed: u64,
    pub target: PathBuf,
    pub decoy: PathBuf,
    pub weights: Option<PathBuf>,
    pub tdc: bool,
    pub threads: usize,
    pub static_model: bool,
    pub init_weights: Option<PathBuf>,
    pub overrule: bool,
    pub pin: PathBuf,
}

const VALUE_FLAGS: [&str; 6] = [
    "--seed",
    "--results-psms",
    "--decoy-results-psms",
    "--weights",
    "--num-threads",
    "--init-weights",
];

const SWITCHES: [&str; 5] = [
    "--only-psms",
    "-Y",
    "--post-processing-tdc",
    "--static",
    "--override",
];

fn bounded(flag: &str, value: &str, min: u64, max: u64) -> Result<u64, String> {
    match value.parse::<u64>() {
        Ok(v) if (min..=max).contains(&v) => Ok(v),
        _ => Err(format!(
            "{flag} must be an integer from {min} to {max}, got '{value}'"
        )),
    }
}

pub fn parse(args: &[String]) -> Result<Args, String> {
    let mut seed = 1u64;
    let mut target = None;
    let mut decoy = None;
    let mut weights = None;
    let mut tdc = false;
    let mut threads = 3usize;
    let mut static_model = false;
    let mut init_weights = None;
    let mut overrule = false;
    let mut pin = None;
    let mut i = 0usize;
    while i < args.len() {
        let arg = args[i].as_str();
        i += 1;
        let (flag, inline) = match arg.split_once('=') {
            Some((f, v)) if f.starts_with("--") => (f, Some(v.to_string())),
            _ => (arg, None),
        };
        if VALUE_FLAGS.contains(&flag) {
            let value = match inline {
                Some(v) => v,
                None => {
                    let v = args
                        .get(i)
                        .ok_or_else(|| format!("argument {flag} needs a value"))?;
                    i += 1;
                    v.clone()
                }
            };
            match flag {
                "--seed" => seed = bounded(flag, &value, 1, 20000)?,
                "--num-threads" => threads = bounded(flag, &value, 1, 128)? as usize,
                "--results-psms" => target = Some(PathBuf::from(value)),
                "--decoy-results-psms" => decoy = Some(PathBuf::from(value)),
                "--weights" => weights = Some(PathBuf::from(value)),
                _ => init_weights = Some(PathBuf::from(value)),
            }
        } else if SWITCHES.contains(&arg) {
            match arg {
                "-Y" | "--post-processing-tdc" => tdc = true,
                "--static" => static_model = true,
                "--override" => overrule = true,
                _ => {}
            }
        } else if arg.starts_with('-') && arg.len() > 1 {
            return Err(format!("unknown argument: {arg}"));
        } else if i == args.len() {
            pin = Some(PathBuf::from(arg));
        } else {
            return Err(format!(
                "unexpected argument: {arg} (the PIN file comes last)"
            ));
        }
    }
    if static_model && init_weights.is_none() {
        return Err(
            "Error: the --static option requires the --init-weights option to be specified."
                .to_string(),
        );
    }
    if !static_model && init_weights.is_some() {
        return Err("argument --init-weights is only supported with --static".to_string());
    }
    Ok(Args {
        seed,
        target: target.ok_or("argument --results-psms is required")?,
        decoy: decoy.ok_or("argument --decoy-results-psms is required")?,
        weights,
        tdc,
        threads,
        static_model,
        init_weights,
        overrule,
        pin: pin.ok_or("no PIN file given: the PIN path comes last")?,
    })
}
