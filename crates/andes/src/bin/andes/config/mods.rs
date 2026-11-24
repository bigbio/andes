//! Reads the modification and crosslinker declarations of a `--config` file.
//!
//! A modification or crosslinker is written as one `Name (sites)` string or as a map with
//! the fields of the declaration; both forms resolve to the same values as the command-line
//! flags.

pub fn mod_list(
    key: &str,
    items: &[serde_yaml::Value],
    fixed: bool,
) -> Result<crate::cli::ModList, String> {
    use model::modification::spec::ModDecl;
    let mut mods = Vec::new();
    for (i, item) in items.iter().enumerate() {
        let at = |e: String| format!("{key}[{i}]: {e}");
        let decl = match item {
            serde_yaml::Value::String(s) => ModDecl::parse(s).map_err(|e| at(e.to_string()))?,
            serde_yaml::Value::Mapping(map) => mod_decl(map).map_err(at)?,
            _ => {
                return Err(at(
                    "expected `Name (sites)` or a map with `name` or `mass` and `sites`".into(),
                ))
            }
        };
        mods.extend(
            decl.resolve(fixed, model::unimod::bundled())
                .map_err(|e| at(e.to_string()))?,
        );
    }
    Ok(crate::cli::ModList(mods))
}

fn mod_decl(map: &serde_yaml::Mapping) -> Result<model::modification::spec::ModDecl, String> {
    use serde_yaml::Value;
    let text = |key: &str, v: &Value| -> Result<String, String> {
        match v {
            Value::String(s) => Ok(s.clone()),
            _ => Err(format!("`{key}` must be a string")),
        }
    };
    let numbers = |key: &str, v: &Value| -> Result<Vec<f64>, String> {
        let one = |v: &Value| v.as_f64().ok_or_else(|| format!("`{key}` takes numbers"));
        match v {
            Value::Sequence(seq) => seq.iter().map(one).collect(),
            v => Ok(vec![one(v)?]),
        }
    };
    let mut d = model::modification::spec::ModDecl::default();
    let mut sites = None;
    for (k, v) in map {
        let key = k.as_str().ok_or("modification keys are strings")?;
        match key {
            "name" => d.name = Some(text(key, v)?),
            "mass" => {
                d.mass = Some(
                    v.as_f64()
                        .ok_or("`mass` is a number (the mass shift in Da)")?,
                )
            }
            "sites" => sites = Some(text(key, v)?),
            "neutral_loss" => d.neutral_loss = numbers(key, v)?,
            "loss_class" => d.loss_class = Some(text(key, v)?),
            "restrict_to" => d.restrict_to = Some(text(key, v)?),
            "cleavage" => d.cleavage = Some(text(key, v)?),
            "diagnostic_ions" => d.diagnostic_ions = numbers(key, v)?,
            other => {
                return Err(format!(
                    "unknown key `{other}` (expected name, mass, sites, neutral_loss, \
                     loss_class, restrict_to, cleavage, diagnostic_ions)"
                ))
            }
        }
    }
    match (sites, d.name.clone(), d.mass) {
        (Some(s), _, _) => d.sites = s,
        (None, Some(name), None) => {
            let parsed = model::modification::spec::ModDecl::parse(&name).map_err(|e| e.to_string())?;
            d.name = parsed.name;
            d.sites = parsed.sites;
        }
        (None, _, _) => return Err("`sites` is missing".into()),
    }
    Ok(d)
}

pub fn xl_linker(key: &str, v: &serde_yaml::Value) -> Result<xl::crosslinker::Crosslinker, String> {
    use serde_yaml::Value;
    use xl::spec::XlDecl;
    let at = |e: String| format!("{key}: {e}");
    let db = model::unimod::bundled();
    match v {
        Value::String(s) => XlDecl::parse(s)
            .and_then(|d| d.resolve(db))
            .map_err(|e| at(e.to_string())),
        Value::Mapping(map) if map.contains_key("intact") => {
            let r: model::crosslinkers::Reagent =
                serde_yaml::from_value(v.clone()).map_err(|e| at(e.to_string()))?;
            xl::spec::custom(&r, db).map_err(|e| at(e.to_string()))
        }
        Value::Mapping(map) => {
            let mut name = None;
            let mut sites = None;
            for (k, val) in map {
                let text = || {
                    val.as_str().map(str::to_string).ok_or_else(|| {
                        at(format!("`{}` must be a string", k.as_str().unwrap_or("?")))
                    })
                };
                match k.as_str() {
                    Some("name") => name = Some(text()?),
                    Some("sites") => sites = Some(text()?),
                    other => {
                        return Err(at(format!(
                            "unknown key `{}` (a named reagent takes name and sites; a custom \
                             reagent also takes intact, stubs, doublets, monolinks, aliases, \
                             xlmod, source)",
                            other.unwrap_or("?")
                        )))
                    }
                }
            }
            let name = name.ok_or_else(|| at("`name` is missing".into()))?;
            let mut d = XlDecl::parse(&name).map_err(|e| at(e.to_string()))?;
            if sites.is_some() {
                d.sites = sites;
            }
            d.resolve(db).map_err(|e| at(e.to_string()))
        }
        _ => Err(at("expected `Name (sites)` or a map with `name`".into())),
    }
}

/// The crosslinkers of the `search.xl` value: one declaration, or a list of them.
pub fn xl_set(v: &serde_yaml::Value) -> Result<crate::cli::XlSet, String> {
    let linkers = match v {
        serde_yaml::Value::Sequence(items) => items
            .iter()
            .map(|i| xl_linker("search.xl", i))
            .collect::<Result<Vec<_>, _>>()?,
        _ => vec![xl_linker("search.xl", v)?],
    };
    crate::cli::XlSet::new("search.xl", linkers)
}
