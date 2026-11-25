//! Tests of labelled-pair quantification of a run.

use super::fit::{decoy_probes, fit_at_identified_profile, fit_folded, learn_width, probe_of};
use super::hits::{heavy_proline_delta, labelled, satellite_on_another_channel, ChannelHit, Group};
use super::proteins::{protein_ratio_tsv, protein_ratios};
use super::settings::{max_label_sites, search_label_mods};
use quant::label::{ChannelEnvelope, LabelSpec};
use quant::Ms1RunIndex;

use super::*;

#[test]
fn arg_pro_satellite_is_folded_into_its_channel() {
    use model::scan::Ms1Scan;
    let spec = LabelSpec::parse("silac2").unwrap();
    let arg10 = spec.channels[1].label_at(Site::Residue(b'R')).unwrap();
    let d_pro = heavy_proline_delta(arg10).unwrap();
    assert!((d_pro - 6.013809).abs() < 1e-4, "{d_pro}");
    let light = Formula::peptide(b"PEPTIDER".iter().map(|&r| (r, 0.0)));
    let sites = spec.sites();
    let counts: Vec<u32> = sites
        .iter()
        .map(|s| u32::from(*s == Site::Residue(b'R')))
        .collect();
    let envs: Vec<ChannelEnvelope> = (0..2)
        .map(|c| spec.channel_envelope(c, &sites, &counts, &light))
        .collect();
    let z = 2.0;
    let ml = 478.73;
    let mh = ml + 10.008269 / z;
    let mp = mh + d_pro / z;
    let species = [
        (ml, envs[0], 1e6),
        (mh, envs[1], 0.8e6),
        (mp, envs[1], 0.2e6),
    ];
    let scans: Vec<Ms1Scan> = (0..200)
        .map(|i| {
            let rt = i as f64;
            let g = (-(rt - 100.0f64).powi(2) / 32.0).exp();
            let mut peaks: Vec<(f64, f32)> =
                (0..40).map(|k| (300.0 + k as f64 * 23.1, 50.0)).collect();
            for (mz, env, amp) in species {
                for off in -2..=5 {
                    let v = amp * env.at(off) * g;
                    if v > 50.0 {
                        peaks.push((mz + off as f64 * model::mass::ISOTOPE / z, v as f32));
                    }
                }
            }
            peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
            Ms1Scan {
                rt,
                peaks,
                id: format!("scan={i}"),
                ..Default::default()
            }
        })
        .collect();
    let ms1 = Ms1RunIndex::new(scans);
    let mut g = Group {
        key: "PEPTIDER".into(),
        sequence: "PEPTIDER".into(),
        charge: 2,
        proteins: Vec::new(),
        channel_mz: vec![ml, mh],
        envelopes: envs.clone(),
        deltas: vec![0.0, 10.008269],
        faims_cv: None,
        satellites: vec![(1, mp)],
        identity_decoy: false,
    };
    assert!(!satellite_on_another_channel(&g, 10.0));
    let p = probe_of(&g, 100.0, 4.0, 4.0, &[0.0, 0.0]);
    let (fit, converted) = fit_folded(&ms1, &g, &p, &FitParams::default()).unwrap();
    let (r, _) = log2_ratio(&fit.channels[1], &fit.channels[0]).unwrap();
    assert!(r.abs() < 0.02, "{r}");
    assert!((converted[1].unwrap() - 0.2).abs() < 0.01, "{converted:?}");
    g.satellites.clear();
    let p = probe_of(&g, 100.0, 4.0, 4.0, &[0.0, 0.0]);
    let (fit, _) = fit_folded(&ms1, &g, &p, &FitParams::default()).unwrap();
    let (r, _) = log2_ratio(&fit.channels[1], &fit.channels[0]).unwrap();
    assert!((r - 0.8f64.log2()).abs() < 0.02, "{r}");
}

#[test]
fn missing_channel_scores_like_its_decoys_next_to_an_interferer() {
    use model::scan::Ms1Scan;
    let env = ChannelEnvelope::new(&[0.55, 0.3, 0.11, 0.04], 0, 1.0, 0, 1.0);
    let z = 2.0;
    let ml = 512.27;
    let g = Group {
        key: "PEPTIDEK".into(),
        sequence: "PEPTIDEK".into(),
        charge: 2,
        proteins: Vec::new(),
        channel_mz: vec![ml, ml + 8.014199 / z],
        envelopes: vec![env, env],
        deltas: vec![0.0, 8.014199],
        faims_cv: None,
        satellites: Vec::new(),
        identity_decoy: false,
    };
    let ids = vec![Some(ChannelHit { hi: 0, rt: 100.0 }), None];
    let probe = probe_of(&g, 100.0, 4.0, 4.0, &[0.0, 0.0]);
    let params = FitParams::default();
    let run = |mz: f64| {
        let scans: Vec<Ms1Scan> = (0..600)
            .map(|k| {
                let rt = k as f64 / 3.0;
                let mut peaks: Vec<(f64, f32)> =
                    (0..60).map(|j| (350.0 + 11.3 * j as f64, 50.0)).collect();
                for (m, apex, amp) in [(ml, 100.0, 1e6), (mz, 106.0, 5e6)] {
                    let e = (-(rt - apex) * (rt - apex) / 32.0f64).exp();
                    for off in 0..4 {
                        let v = amp * env.at(off) * e;
                        if v > 50.0 {
                            peaks.push((m + off as f64 * model::mass::ISOTOPE / z, v as f32));
                        }
                    }
                }
                peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
                Ms1Scan {
                    rt,
                    peaks,
                    id: format!("scan={k}"),
                    ..Default::default()
                }
            })
            .collect();
        Ms1RunIndex::new(scans)
    };
    let decoy_mz = decoy_probes(&g, &probe, &[1], 10.0, 7, 1)[0].channel_mz[1];
    assert!((decoy_mz - g.channel_mz[1]).abs() > 0.01);

    let at_heavy = run(g.channel_mz[1]);
    let (free, _) = fit_folded(&at_heavy, &g, &probe, &params).unwrap();
    assert!(free.t0 > 102.0, "free apex {}", free.t0);
    let ((target, _), _) = fit_at_identified_profile(&at_heavy, &g, &probe, &ids, &params).unwrap();
    assert!((target.t0 - 100.0).abs() < 0.5, "apex {}", target.t0);

    let at_decoy = run(decoy_mz);
    let (_, pinned) = fit_at_identified_profile(&at_decoy, &g, &probe, &ids, &params).unwrap();
    let dp = decoy_probes(&g, &pinned, &[1], 10.0, 7, 1).remove(0);
    assert_eq!(dp.channel_mz[1], decoy_mz);
    let (decoy, _) = fit_folded(&at_decoy, &g, &dp, &params).unwrap();

    let (st, sd) = (
        channel_score(&target.channels[1]),
        channel_score(&decoy.channels[1]),
    );
    assert!(
        (st - sd).abs() <= 0.02 * st.max(sd) + 1e-9,
        "target {st} vs decoy {sd}"
    );
    assert!(
        (target.channels[1].apex_offset - decoy.channels[1].apex_offset).abs() < 0.2,
        "{} vs {}",
        target.channels[1].apex_offset,
        decoy.channels[1].apex_offset
    );
    let (r, _) = log2_ratio(&target.channels[0], &decoy.channels[0]).unwrap();
    assert!(r.abs() < 0.02, "{r}");
}

#[test]
fn learned_width_is_stable() {
    use model::scan::Ms1Scan;
    let env = ChannelEnvelope::new(&[0.55, 0.3, 0.11, 0.04], 0, 1.0, 0, 1.0);
    let n_groups = 60;
    let (sigma, z) = (4.0, 2.0);
    let groups: Vec<Group> = (0..n_groups)
        .map(|i| {
            let ml = 400.0 + 7.3 * i as f64;
            Group {
                key: format!("PEPTIDE{i}K"),
                sequence: format!("PEPTIDE{i}K"),
                charge: 2,
                proteins: Vec::new(),
                channel_mz: vec![ml, ml + 8.014199 / z],
                envelopes: vec![env, env],
                deltas: vec![0.0, 8.014199],
                faims_cv: None,
                satellites: Vec::new(),
                identity_decoy: false,
            }
        })
        .collect();
    let apex = |i: usize| 40.0 + 9.0 * i as f64;
    let scans: Vec<Ms1Scan> = (0..1800)
        .map(|k| {
            let rt = k as f64 / 3.0;
            let mut peaks: Vec<(f64, f32)> =
                (0..60).map(|j| (350.0 + 11.3 * j as f64, 50.0)).collect();
            for (i, g) in groups.iter().enumerate() {
                let e = (-(rt - apex(i)).powi(2) / (2.0 * sigma * sigma)).exp();
                for (c, amp) in [(0, 1e6), (1, 5e5)] {
                    for off in 0..4 {
                        let v = amp * env.at(off) * e;
                        if v > 50.0 {
                            let m = g.channel_mz[c] + off as f64 * model::mass::ISOTOPE / z;
                            peaks.push((m, v as f32));
                        }
                    }
                }
            }
            peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
            Ms1Scan {
                rt,
                peaks,
                id: format!("scan={k}"),
                ..Default::default()
            }
        })
        .collect();
    let ms1 = Ms1RunIndex::new(scans);
    let anchors: Vec<(usize, f64)> = (0..n_groups).map(|i| (i, apex(i))).collect();
    let all = learn_width(&ms1, &groups, &anchors, 10.0);
    let fewer = learn_width(&ms1, &groups, &anchors[1..], 10.0);
    assert!((all - sigma).abs() < 0.04 * sigma, "learned {all}");
    assert!((all - fewer).abs() < 0.02 * all, "{all} vs {fewer}");
}

#[test]
fn satellite_overlap_check_covers_the_modelled_isotopes() {
    let iso = model::mass::ISOTOPE / 2.0;
    let group = |channel_mz: Vec<f64>, satellites: Vec<(usize, f64)>| Group {
        key: "PEPTIDER".into(),
        sequence: "PEPTIDER".into(),
        charge: 2,
        proteins: Vec::new(),
        envelopes: vec![ChannelEnvelope::new(&[1.0], 0, 1.0, 0, 1.0); channel_mz.len()],
        deltas: vec![0.0; channel_mz.len()],
        channel_mz,
        faims_cv: None,
        satellites,
        identity_decoy: false,
    };
    assert!(satellite_on_another_channel(
        &group(vec![500.0, 510.0], vec![(0, 510.0 - iso)]),
        10.0
    ));
    assert!(satellite_on_another_channel(
        &group(vec![500.0, 510.0], vec![(0, 510.0 + 5.0 * iso)]),
        10.0
    ));
    assert!(satellite_on_another_channel(
        &group(vec![500.0, 520.0], vec![(0, 540.0), (1, 540.0 + iso)]),
        10.0
    ));
    assert!(!satellite_on_another_channel(
        &group(vec![500.0, 520.0], vec![(0, 540.0), (1, 560.0)]),
        10.0
    ));
}

#[test]
fn four_channels_with_arg_pro_satellites_are_fitted() {
    use model::scan::Ms1Scan;
    let spec = LabelSpec::parse("[][Arg6][R+8.014199:6C2N][Arg10]").unwrap();
    let light = Formula::peptide(b"PEPTIDER".iter().map(|&r| (r, 0.0)));
    let sites = spec.sites();
    let counts: Vec<u32> = sites
        .iter()
        .map(|s| u32::from(*s == Site::Residue(b'R')))
        .collect();
    let envs: Vec<ChannelEnvelope> = (0..4)
        .map(|c| spec.channel_envelope(c, &sites, &counts, &light))
        .collect();
    let z = 2.0;
    let ml = 478.73;
    let deltas = [0.0, 6.020129, 8.014199, 10.008269];
    let mz: Vec<f64> = deltas.iter().map(|d| ml + d / z).collect();
    let sats: Vec<(usize, f64)> = (1..4)
        .map(|c| {
            let arg = spec.channels[c].label_at(Site::Residue(b'R')).unwrap();
            (c, mz[c] + heavy_proline_delta(arg).unwrap() / z)
        })
        .collect();
    let amps = [1e6, 0.5e6, 0.25e6, 2e6];
    let mut species: Vec<(f64, ChannelEnvelope, f64)> = vec![(mz[0], envs[0], amps[0])];
    for &(c, mp) in &sats {
        species.push((mz[c], envs[c], 0.8 * amps[c]));
        species.push((mp, envs[c], 0.2 * amps[c]));
    }
    let scans: Vec<Ms1Scan> = (0..200)
        .map(|i| {
            let rt = i as f64;
            let g = (-(rt - 100.0f64).powi(2) / 32.0).exp();
            let mut peaks: Vec<(f64, f32)> =
                (0..40).map(|k| (300.0 + k as f64 * 23.1, 50.0)).collect();
            for &(m, env, amp) in &species {
                for off in -2..=5 {
                    let v = amp * env.at(off) * g;
                    if v > 50.0 {
                        peaks.push((m + off as f64 * model::mass::ISOTOPE / z, v as f32));
                    }
                }
            }
            peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
            Ms1Scan {
                rt,
                peaks,
                id: format!("scan={i}"),
                ..Default::default()
            }
        })
        .collect();
    let ms1 = Ms1RunIndex::new(scans);
    let g = Group {
        key: "PEPTIDER".into(),
        sequence: "PEPTIDER".into(),
        charge: 2,
        proteins: Vec::new(),
        channel_mz: mz.clone(),
        envelopes: envs.clone(),
        deltas: deltas.to_vec(),
        faims_cv: None,
        satellites: sats,
        identity_decoy: false,
    };
    assert!(satellite_on_another_channel(&g, 10.0));
    let p = probe_of(&g, 100.0, 4.0, 4.0, &[0.0; 4]);
    assert_eq!(p.channel_mz.len(), 7);
    let (fit, converted) = fit_folded(&ms1, &g, &p, &FitParams::default()).unwrap();
    for c in 1..4 {
        let (r, _) = log2_ratio(&fit.channels[c], &fit.channels[0]).unwrap();
        let want = (amps[c] / amps[0]).log2();
        assert!((r - want).abs() < 0.05, "channel {c}: {r} vs {want}");
        assert!((converted[c].unwrap() - 0.2).abs() < 0.03, "{converted:?}");
    }
}

#[test]
fn terminal_mod_leaves_the_residue_label_site() {
    use model::amino_acid::AminoAcid;
    use model::modification::{ModLocation, Modification, ResidueSpec};
    use model::peptide::Peptide;
    let spec = LabelSpec::parse("silac2").unwrap();
    let m = |name: &str, delta: f64, r: u8, loc: ModLocation| -> &'static Modification {
        Box::leak(Box::new(Modification {
            name: name.into(),
            mass_delta: delta,
            residue: ResidueSpec::Specific(r),
            location: loc,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: Vec::new(),
            parts: Vec::new(),
        }))
    };
    let lys8 = spec.channels[1].label_at(Site::Residue(b'K')).unwrap();
    let arg10 = spec.channels[1].label_at(Site::Residue(b'R')).unwrap();
    let acetyl = m("Acetyl", 42.010565, b'K', ModLocation::ProtNTerm);
    let heavy_k = m(
        &format!("{}+Acetyl", lys8.name),
        lys8.delta + 42.010565,
        b'K',
        ModLocation::ProtNTerm,
    );
    let heavy_r = m(&arg10.name, arg10.delta, b'R', ModLocation::Anywhere);
    let cand = |first: &'static Modification, last: Option<&'static Modification>| {
        let mut res: Vec<AminoAcid> = b"KPEPTIDER"
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap().clone())
            .collect();
        res[0] = res[0].clone().with_mod(first);
        if let Some(l) = last {
            res[8] = res[8].clone().with_mod(l);
        }
        Candidate {
            peptide: Peptide::new(res, b'_', b'A'),
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: true,
            is_protein_c_term: false,
        }
    };
    let light = labelled(&cand(acetyl, None), &spec).unwrap();
    let heavy = labelled(&cand(heavy_k, Some(heavy_r)), &spec).unwrap();
    assert_eq!(light.channel, 0);
    assert_eq!(heavy.channel, 1);
    let sorted = |v: &[Site]| {
        let mut v = v.to_vec();
        v.sort_by_key(|s| format!("{s:?}"));
        v
    };
    assert_eq!(sorted(&light.instances), sorted(&heavy.instances));
    assert_eq!(light.instances.len(), 2);
    assert_eq!(light.key, heavy.key);
    assert!((light.light_mass - heavy.light_mass).abs() < 1e-6);
}

#[test]
fn protein_ratio_is_the_median_of_peptide_medians() {
    let row = |pep: &str, z: u8, v: f64, status: &str| output::LabelRow {
        run: "r".into(),
        peptide: pep.into(),
        sequence: pep.into(),
        charge: z,
        proteins: vec!["P1".into()],
        channel: 2,
        channel_name: "heavy".into(),
        label: None,
        status: status.into(),
        spec_id: None,
        psm_q_value: None,
        rt_psm: None,
        labelled_peptidoform: None,
        mz: 0.0,
        rt_apex: 0.0,
        rt_shift: 0.0,
        width: 0.0,
        quantity: 0.0,
        quantity_se: 0.0,
        quantity_upper: 0.0,
        llr: 0.0,
        envelope_cos: 0.0,
        coherence: 0.0,
        m1_ratio: 0.0,
        mass_err_ppm: Some(0.0),
        channel_q_value: None,
        group_llr: 0.0,
        transfer: false,
        transfer_q_value: None,
        log2_ratio: Some(v),
        log2_ratio_se: None,
        log2_ratio_normalized: Some(v - 0.5),
        converted_fraction: None,
        trace_quantity: 0.0,
    };
    let rows = vec![
        row("A", 2, 0.5, "identified"),
        row("A", 3, 1.5, "identified"),
        row("B", 2, 3.0, "match_channel"),
        row("C", 2, -1.0, "identified"),
    ];
    let r = protein_ratios(&rows, |_| None, 2);
    assert_eq!(r.len(), 1);
    assert_eq!(r[0].ratio_count, 3);
    assert_eq!(r[0].inferred_count, 1);
    assert!((r[0].log2_ratio - 1.0).abs() < 1e-12);
    assert!((r[0].log2_ratio_normalized.unwrap() - 0.5).abs() < 1e-12);
    assert!(protein_ratios(&rows[3..], |_| None, 2).is_empty());
    let tsv = protein_ratio_tsv(&r);
    assert!(tsv
        .lines()
        .nth(1)
        .unwrap()
        .starts_with("P1\tP1\tr\t2\theavy\t1.000000"));
}

#[test]
fn max_sites_follow_missed_cleavages() {
    use model::enzyme::Enzyme;
    let db = model::ProteinDb {
        proteins: vec![model::Protein {
            accession: "P1".into(),
            description: String::new(),
            sequence: b"GGLRLRLLRAGLRGLRLGGGGGGGGGGG".to_vec(),
            annotation: None,
        }],
    };
    let s = LabelSpec::parse("silac2").unwrap();
    assert_eq!(max_label_sites(&s, &[Enzyme::Trypsin], 2, 2, &db, 50), 3);
    assert_eq!(max_label_sites(&s, &[Enzyme::Trypsin], 1, 2, &db, 50), 3);
    let d = LabelSpec::parse("dimethyl3").unwrap();
    assert_eq!(max_label_sites(&d, &[Enzyme::Trypsin], 2, 1, &db, 50), 3);
    let l = LabelSpec::parse("[][Leu3]").unwrap();
    assert_eq!(max_label_sites(&l, &[Enzyme::Trypsin], 2, 2, &db, 50), 7);
    assert_eq!(max_label_sites(&l, &[Enzyme::Trypsin], 2, 2, &db, 6), 4);
    assert_eq!(max_label_sites(&s, &[Enzyme::LysC], 2, 2, &db, 50), 3 + 5);
    assert_eq!(max_label_sites(&s, &[Enzyme::Trypsin], 0, 2, &db, 50), 5);
}

#[test]
fn search_mods_name_channels() {
    let d = LabelSpec::parse("dimethyl3").unwrap();
    let (mods, channels) = search_label_mods(&d);
    assert_eq!(mods.len(), 6);
    assert_eq!(channels.len(), 3);
    assert_eq!(channels[0], ("Dimethyl".to_string(), 1));
    let s = LabelSpec::parse("silac2").unwrap();
    let (mods, channels) = search_label_mods(&s);
    assert_eq!(mods.len(), 2);
    assert!(channels.iter().all(|c| c.1 == 2));
}
