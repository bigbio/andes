//! Instrument detection from the mzML header: the instrument model term, or the mass
//! analyzer of the configuration that the MS2 scans reference.

use std::collections::HashMap;
use std::io::BufRead;

use quick_xml::{events::Event, Reader};

use model::InstrumentType;

use super::attr_str;
use super::cv::{
    CV_ANALYZER_FTICR, CV_ANALYZER_ION_TRAP, CV_ANALYZER_LINEAR_ION_TRAP, CV_ANALYZER_ORBITRAP,
    CV_ANALYZER_QUAD_ION_TRAP, CV_ANALYZER_RADIAL_LIT, CV_ANALYZER_TOF, CV_MODEL_ORBITRAP_ASTRAL,
    CV_MODEL_ORBITRAP_FUSION, CV_MODEL_Q_EXACTIVE, CV_MODEL_Q_EXACTIVE_HF,
    CV_MODEL_Q_EXACTIVE_HF_X, CV_MODEL_Q_EXACTIVE_PLUS, CV_MS_LEVEL,
};

pub fn detect_instrument_type<R: BufRead>(reader: R) -> Option<InstrumentType> {
    let mut xml = Reader::from_reader(reader);
    xml.trim_text(true);

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum S {
        Outside,
        ReferenceableParamGroupList,
        ReferenceableParamGroup,
        InstrumentConfigurationList,
        InstrumentConfiguration,
        ComponentListAnalyzer,
        Run,
        Spectrum,
        Scan,
    }

    let mut state = S::Outside;
    let mut buf: Vec<u8> = Vec::with_capacity(4096);

    let mut ic_map: HashMap<String, InstrumentType> = HashMap::new();
    let mut current_ic_id: Option<String> = None;
    let mut current_ic_type: Option<InstrumentType> = None;
    let mut current_ic_analyzer: Option<InstrumentType> = None;

    let mut rpg_map: HashMap<String, InstrumentType> = HashMap::new();
    let mut current_rpg_id: Option<String> = None;
    let mut current_rpg_type: Option<InstrumentType> = None;

    let mut default_ic_ref: Option<String> = None;

    let mut ms2_counts: HashMap<InstrumentType, usize> = HashMap::new();
    let mut current_spec_is_ms2: Option<bool> = None;
    let mut current_spec_ic_ref: Option<String> = None;
    let mut ms2_seen: usize = 0;

    const MAX_PEEK: usize = 64;

    loop {
        buf.clear();
        let event = match xml.read_event_into(&mut buf) {
            Ok(e) => e,
            Err(_) => break,
        };
        match event {
            Event::Eof => break,

            Event::Start(ref e) => {
                let tag = e.local_name().as_ref().to_owned();
                match tag.as_slice() {
                    b"referenceableParamGroupList" if state == S::Outside => {
                        state = S::ReferenceableParamGroupList;
                    }
                    b"referenceableParamGroup" if state == S::ReferenceableParamGroupList => {
                        current_rpg_id = attr_str(e, b"id");
                        current_rpg_type = None;
                        state = S::ReferenceableParamGroup;
                    }
                    b"instrumentConfigurationList" if state == S::Outside => {
                        state = S::InstrumentConfigurationList;
                    }
                    b"instrumentConfiguration" if state == S::InstrumentConfigurationList => {
                        current_ic_id = attr_str(e, b"id");
                        current_ic_type = None;
                        current_ic_analyzer = None;
                        state = S::InstrumentConfiguration;
                    }
                    b"analyzer" if state == S::InstrumentConfiguration => {
                        state = S::ComponentListAnalyzer;
                    }
                    b"run" if state == S::Outside => {
                        default_ic_ref = attr_str(e, b"defaultInstrumentConfigurationRef");
                        state = S::Run;
                    }
                    b"spectrum" if state == S::Run => {
                        current_spec_is_ms2 = None;
                        current_spec_ic_ref = None;
                        state = S::Spectrum;
                    }
                    b"scan" if state == S::Spectrum => {
                        if let Some(r) = attr_str(e, b"instrumentConfigurationRef") {
                            current_spec_ic_ref = Some(r);
                        }
                        state = S::Scan;
                    }
                    _ => {}
                }
            }

            Event::Empty(ref e) => {
                let tag = e.local_name().as_ref().to_owned();
                if tag == b"scan" && state == S::Spectrum {
                    if let Some(r) = attr_str(e, b"instrumentConfigurationRef") {
                        current_spec_ic_ref = Some(r);
                    }
                }
                if tag == b"referenceableParamGroupRef" && state == S::InstrumentConfiguration {
                    if let Some(r) = attr_str(e, b"ref") {
                        if let Some(&t) = rpg_map.get(&r) {
                            current_ic_type = Some(t);
                        }
                    }
                }
                if tag == b"cvParam" {
                    let acc = attr_str(e, b"accession").unwrap_or_default();
                    match state {
                        S::ComponentListAnalyzer => {
                            let typ = match acc.as_str() {
                                CV_ANALYZER_ORBITRAP => Some(InstrumentType::QExactive),
                                CV_ANALYZER_FTICR => Some(InstrumentType::HighRes),
                                CV_ANALYZER_TOF => Some(InstrumentType::TOF),
                                CV_ANALYZER_ION_TRAP
                                | CV_ANALYZER_QUAD_ION_TRAP
                                | CV_ANALYZER_RADIAL_LIT
                                | CV_ANALYZER_LINEAR_ION_TRAP => Some(InstrumentType::LowRes),
                                _ => None,
                            };
                            if let Some(t) = typ {
                                if current_ic_type.is_none() {
                                    current_ic_type = Some(t);
                                }
                                current_ic_analyzer.get_or_insert(t);
                            }
                        }
                        S::InstrumentConfiguration => {
                            let cv_name = attr_str(e, b"name").unwrap_or_default();
                            if let Some(t) = instrument_model_type(acc.as_str(), &cv_name) {
                                current_ic_type = Some(t);
                            }
                        }
                        S::ReferenceableParamGroup => {
                            let cv_name = attr_str(e, b"name").unwrap_or_default();
                            if current_rpg_type.is_none() {
                                if let Some(t) = instrument_model_type(acc.as_str(), &cv_name) {
                                    current_rpg_type = Some(t);
                                }
                            }
                        }
                        S::Spectrum => {
                            if acc == CV_MS_LEVEL {
                                let val = attr_str(e, b"value").unwrap_or_default();
                                current_spec_is_ms2 = Some(val == "2");
                            }
                        }
                        _ => {}
                    }
                }
            }

            Event::End(ref e) => {
                let tag = e.local_name().as_ref().to_owned();
                match tag.as_slice() {
                    b"referenceableParamGroup" if state == S::ReferenceableParamGroup => {
                        if let (Some(id), Some(t)) =
                            (current_rpg_id.take(), current_rpg_type.take())
                        {
                            rpg_map.insert(id, t);
                        }
                        state = S::ReferenceableParamGroupList;
                    }
                    b"referenceableParamGroupList" if state == S::ReferenceableParamGroupList => {
                        state = S::Outside;
                    }
                    b"analyzer" if state == S::ComponentListAnalyzer => {
                        state = S::InstrumentConfiguration;
                    }
                    b"instrumentConfiguration" if state == S::InstrumentConfiguration => {
                        let ic_type = match current_ic_analyzer.take() {
                            Some(InstrumentType::LowRes) => Some(InstrumentType::LowRes),
                            _ => current_ic_type.take(),
                        };
                        if let (Some(id), Some(t)) = (current_ic_id.take(), ic_type) {
                            ic_map.insert(id, t);
                        }
                        state = S::InstrumentConfigurationList;
                    }
                    b"instrumentConfigurationList" if state == S::InstrumentConfigurationList => {
                        state = S::Outside;
                    }
                    b"scan" if state == S::Scan => {
                        state = S::Spectrum;
                    }
                    b"spectrum" if state == S::Spectrum => {
                        let is_ms2 = current_spec_is_ms2.unwrap_or(false);
                        if is_ms2 {
                            let ic_ref = current_spec_ic_ref
                                .clone()
                                .or_else(|| default_ic_ref.clone());
                            if let Some(r) = ic_ref {
                                if let Some(&t) = ic_map.get(&r) {
                                    *ms2_counts.entry(t).or_insert(0) += 1;
                                }
                            }
                            ms2_seen += 1;
                            if ms2_seen >= MAX_PEEK {
                                break;
                            }
                        }
                        current_spec_is_ms2 = None;
                        current_spec_ic_ref = None;
                        state = S::Run;
                    }
                    b"run" if state == S::Run => {
                        state = S::Outside;
                    }
                    _ => {}
                }
            }

            _ => {}
        }
    }

    if !ms2_counts.is_empty() {
        return ms2_counts
            .iter()
            .max_by_key(|(&k, &n)| (n, std::cmp::Reverse(k as u8)))
            .map(|(&t, _)| t);
    }

    if let Some(r) = default_ic_ref.as_ref() {
        if let Some(&t) = ic_map.get(r) {
            return Some(t);
        }
    }

    if ic_map.len() == 1 {
        return ic_map.into_values().next();
    }

    None
}

fn instrument_model_type(acc: &str, name: &str) -> Option<InstrumentType> {
    match acc {
        CV_MODEL_ORBITRAP_ASTRAL => Some(InstrumentType::OrbitrapAstral),
        CV_MODEL_Q_EXACTIVE
        | CV_MODEL_Q_EXACTIVE_HF
        | CV_MODEL_Q_EXACTIVE_HF_X
        | CV_MODEL_Q_EXACTIVE_PLUS
        | CV_MODEL_ORBITRAP_FUSION => Some(InstrumentType::QExactive),
        _ => {
            if name.to_ascii_lowercase().contains("astral") {
                Some(InstrumentType::OrbitrapAstral)
            } else {
                None
            }
        }
    }
}
