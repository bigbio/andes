//! Value parsers for the search command-line arguments: numeric ranges, tolerances with units,
//! charge and isotope-error ranges, fractions and crosslinker lists. `cli.rs` declares the
//! arguments, keeps the enumerated-value parser `parse_enum`, and re-exports these parsers by
//! name; config files and `--print-config` reuse them so every path reads the same text.
use model::Tolerance;

pub(crate) fn dash_sep_index(s: &str) -> Option<usize> {
    let b = s.as_bytes();
    (1..b.len()).find(|&i| b[i] == b'-' && b[i - 1].is_ascii_digit())
}

pub(crate) fn parse_int_range<T>(s: &str, label: &str) -> Result<(T, T), String>
where
    T: std::str::FromStr + PartialOrd + std::fmt::Display + Copy,
{
    let trimmed = s.trim();
    let (lo_s, hi_s) = if let Some((a, b)) = trimmed.split_once("..") {
        (a.trim(), b.trim())
    } else if let Some(idx) = dash_sep_index(trimmed) {
        (trimmed[..idx].trim(), trimmed[idx + 1..].trim())
    } else {
        return Err(format!(
            "invalid {label} `{s}`: expected MIN..MAX (or MIN-MAX)"
        ));
    };
    let lo: T = lo_s
        .parse()
        .map_err(|_| format!("invalid {label} `{s}`: bad MIN `{lo_s}`"))?;
    let hi: T = hi_s
        .parse()
        .map_err(|_| format!("invalid {label} `{s}`: bad MAX `{hi_s}`"))?;
    if lo > hi {
        return Err(format!("invalid {label} `{s}`: MIN {lo} > MAX {hi}"));
    }
    Ok((lo, hi))
}

pub(crate) const MAX_SUPPORTED_CHARGE: u8 = 50;

pub(crate) fn parse_charge_range(s: &str) -> Result<(u8, u8), String> {
    let (lo, hi) = parse_int_range::<u8>(s, "charge")?;
    if lo < 1 {
        return Err(format!(
            "invalid charge `{s}`: minimum charge must be >= 1 (got {lo})"
        ));
    }
    if hi > MAX_SUPPORTED_CHARGE {
        return Err(format!(
            "invalid charge `{s}`: maximum charge {hi} exceeds supported bound {MAX_SUPPORTED_CHARGE}"
        ));
    }
    Ok((lo, hi))
}

pub(crate) fn parse_isotope_error_range(s: &str) -> Result<(i8, i8), String> {
    parse_int_range::<i8>(s, "isotope-error")
}

pub(crate) fn parse_precursor_tol(s: &str) -> Result<Tolerance, String> {
    let t = s.trim();
    let lower = t.to_ascii_lowercase();
    let (num_str, is_ppm) = if let Some(n) = lower.strip_suffix("ppm") {
        (n, true)
    } else if let Some(n) = lower.strip_suffix("da") {
        (n, false)
    } else {
        return Err(format!(
            "invalid --precursor-tol `{s}`: expected VALUE+unit, e.g. 20ppm or 0.02da"
        ));
    };
    let v: f64 = num_str
        .trim()
        .parse()
        .map_err(|_| format!("invalid --precursor-tol `{s}`: bad number `{num_str}`"))?;
    if !v.is_finite() || v <= 0.0 {
        return Err(format!(
            "invalid --precursor-tol `{s}`: value must be finite and > 0 (got {v})"
        ));
    }
    Ok(if is_ppm {
        Tolerance::Ppm(v)
    } else {
        Tolerance::Da(v)
    })
}

pub(crate) fn parse_open_search_delta(s: &str) -> Result<(f64, f64), String> {
    let t = s.trim();
    let (min_s, max_s) = t.split_once("..").ok_or_else(|| {
        format!("invalid --open-search-delta `{s}`: expected MIN..MAX in Da, e.g. -100..300")
    })?;
    let parse = |v: &str| -> Result<f64, String> {
        v.trim()
            .parse::<f64>()
            .map_err(|_| format!("invalid --open-search-delta `{s}`: bad number `{}`", v.trim()))
    };
    let (lo, hi) = (parse(min_s)?, parse(max_s)?);
    if !lo.is_finite() || !hi.is_finite() || lo >= hi {
        return Err(format!(
            "invalid --open-search-delta `{s}`: need finite MIN < MAX (got {lo}..{hi})"
        ));
    }
    if lo > 0.0 || hi < 0.0 {
        return Err(format!(
            "invalid --open-search-delta `{s}`: the window must contain 0 (got {lo}..{hi})"
        ));
    }
    if hi - lo > 2000.0 {
        return Err(format!(
            "invalid --open-search-delta `{s}`: width {} Da exceeds the 2000 Da limit",
            hi - lo
        ));
    }
    Ok((lo, hi))
}

pub(crate) fn parse_unit_fraction(s: &str) -> Result<f64, String> {
    let v: f64 = s
        .trim()
        .parse()
        .map_err(|_| format!("invalid value `{s}`: expected a number in [0, 1]"))?;
    if !v.is_finite() || !(0.0..=1.0).contains(&v) {
        return Err(format!(
            "invalid value `{s}`: must be finite and within [0, 1] (got {v})"
        ));
    }
    Ok(v)
}

pub(crate) fn parse_non_negative(s: &str) -> Result<f64, String> {
    let v: f64 = s
        .trim()
        .parse()
        .map_err(|_| format!("invalid value `{s}`: expected a number >= 0"))?;
    if !v.is_finite() || v < 0.0 {
        return Err(format!(
            "invalid value `{s}`: must be finite and >= 0 (got {v})"
        ));
    }
    Ok(v)
}

pub(crate) fn parse_non_negative_f32(s: &str) -> Result<f32, String> {
    let v = parse_non_negative(s)?;
    if v > f32::MAX as f64 {
        return Err(format!("invalid value `{s}`: out of range"));
    }
    Ok(v as f32)
}

pub(crate) fn parse_positive_tol(s: &str) -> Result<f64, String> {
    let v: f64 = s
        .trim()
        .parse()
        .map_err(|_| format!("invalid tolerance `{s}`: expected a positive number"))?;
    if !v.is_finite() || v <= 0.0 {
        return Err(format!(
            "invalid tolerance `{s}`: must be finite and > 0 (got {v})"
        ));
    }
    Ok(v)
}

/// --dia-interference: `auto`, or `F:DRT` with a shared-fragment fraction in (0, 1] and an RT
/// window >= 0 s.
pub(crate) fn parse_dia_interference(s: &str) -> Result<crate::dia::Interference, String> {
    if s.trim().eq_ignore_ascii_case("auto") {
        return Ok(crate::dia::Interference::Auto);
    }
    let err = || {
        format!(
            "{s}: expected auto or F:DRT, a fragment fraction in (0, 1] and seconds, e.g. 0.4:4"
        )
    };
    let (f, drt) = s.split_once(':').ok_or_else(err)?;
    let (f, drt): (f64, f64) = (
        f.trim().parse().map_err(|_| err())?,
        drt.trim().parse().map_err(|_| err())?,
    );
    if !(f > 0.0 && f <= 1.0 && drt >= 0.0 && drt.is_finite()) {
        return Err(err());
    }
    Ok(crate::dia::Interference::Fixed(f, drt))
}

/// The crosslinker items of an `--xl` value: split at `;`, trimmed, empty items dropped.
pub(crate) fn split_xl(s: &str) -> Vec<&str> {
    s.split(';').map(str::trim).filter(|x| !x.is_empty()).collect()
}

pub(crate) fn parse_xl_one(s: &str) -> Result<xl::crosslinker::Crosslinker, String> {
    xl::spec::XlDecl::parse(s)
        .and_then(|d| d.resolve(model::unimod::bundled()))
        .map_err(|e| e.to_string())
}

pub(crate) fn parse_xl(s: &str) -> Result<crate::cli::XlSet, String> {
    let linkers = split_xl(s)
        .into_iter()
        .map(parse_xl_one)
        .collect::<Result<Vec<_>, _>>()?;
    crate::cli::XlSet::new("--xl", linkers)
}
