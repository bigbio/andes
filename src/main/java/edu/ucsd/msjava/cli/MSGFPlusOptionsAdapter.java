package edu.ucsd.msjava.cli;

import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.params.ParamManager.ParamNameEnum;
import edu.ucsd.msjava.params.Parameter;

/**
 * Phase 1 adapter: populates a {@link ParamManager} from a parsed
 * {@link MSGFPlusOptions} by round-tripping each set field through the
 * canonical string form that the existing
 * {@link Parameter#parse(String)} hierarchy expects.
 *
 * This deliberately reuses the legacy parsing logic so Phase 1 is
 * behavior-preserving. Phase 3 deletes the {@code params.Parameter}
 * hierarchy and replaces this adapter with direct construction of the
 * downstream {@code SearchParams}.
 *
 * Returns {@code null} on success, or a human-readable error string
 * matching the format used by {@link ParamManager#parseParams(String[])}.
 */
public final class MSGFPlusOptionsAdapter {

    private MSGFPlusOptionsAdapter() {}

    /**
     * Populate {@code paramManager} (already initialized via
     * {@link ParamManager#addMSGFPlusParams()}) with values from
     * {@code opts}. Caller is responsible for calling
     * {@link ParamManager#isValid()} afterwards if final validation
     * is desired (this method also runs it as the last step).
     */
    public static String adapt(MSGFPlusOptions opts, ParamManager paramManager) {
        String err;

        // Files / paths
        if ((err = setIfPresent(paramManager, ParamNameEnum.CONFIGURATION_FILE,
                opts.configFile == null ? null : opts.configFile.getPath())) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.SPECTRUM_FILE,
                opts.spectrumFile == null ? null : opts.spectrumFile.getPath())) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.DB_FILE,
                opts.databaseFile == null ? null : opts.databaseFile.getPath())) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.SEARCH_OUTPUT_FILE,
                opts.outputFile == null ? null : opts.outputFile.getPath())) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.MOD_FILE,
                opts.modificationFile == null ? null : opts.modificationFile.getPath())) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.DD_DIRECTORY,
                opts.dbIndexDir == null ? null : opts.dbIndexDir.getPath())) != null) return err;

        // Plain strings / domain strings parsed by ToleranceParameter / RangeParameter / EnumParameter
        if ((err = setIfPresent(paramManager, ParamNameEnum.DECOY_PREFIX, opts.decoyPrefix)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.PRECURSOR_MASS_TOLERANCE, opts.precursorTolerance)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.ISOTOPE_ERROR, opts.isotopeErrorRange)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.OUTPUT_FORMAT, opts.outputFormat)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.PRECURSOR_CAL, opts.precursorCalMode)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.MS_LEVEL, opts.msLevel)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.SPEC_INDEX, opts.specIndexRange)) != null) return err;

        // Integer-valued flags (enum + numeric)
        if ((err = setIfPresent(paramManager, ParamNameEnum.PRECURSOR_MASS_TOLERANCE_UNITS, opts.precursorToleranceUnits)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.NUM_THREADS, opts.numThreads)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.NUM_TASKS, opts.numTasks)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.MIN_SPECTRA_PER_THREAD, opts.minSpectraPerThread)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.VERBOSE, opts.verbose)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.TDA_STRATEGY, opts.tdaStrategy)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.FRAG_METHOD, opts.fragMethodId)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.INSTRUMENT_TYPE, opts.instrumentTypeId)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.ENZYME_ID, opts.enzymeId)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.PROTOCOL_ID, opts.protocolId)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.ENZYME_SPECIFICITY, opts.numTolerableTermini)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.MIN_PEPTIDE_LENGTH, opts.minPeptideLength)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.MAX_PEPTIDE_LENGTH, opts.maxPeptideLength)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.MIN_CHARGE, opts.minCharge)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.MAX_CHARGE, opts.maxCharge)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.NUM_MATCHES_SPEC, opts.numMatchesPerSpec)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.ADD_FEATURES, opts.addFeatures)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.MAX_MISSED_CLEAVAGES, opts.maxMissedCleavages)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.MAX_NUM_MODS, opts.maxNumMods)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.ALLOW_DENSE_CENTROIDED_PEAKS, opts.allowDenseCentroidedPeaks)) != null) return err;

        // Hidden integer flags
        if ((err = setIfPresent(paramManager, ParamNameEnum.EDGE_SCORE, opts.edgeScore)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.MIN_NUM_PEAKS, opts.minNumPeaks)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.NUM_ISOFORMS, opts.numIsoforms)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.IGNORE_MET_CLEAVAGE, opts.ignoreMetCleavage)) != null) return err;
        if ((err = setIfPresent(paramManager, ParamNameEnum.MIN_DE_NOVO_SCORE, opts.minDeNovoScore)) != null) return err;

        // Doubles
        if ((err = setIfPresent(paramManager, ParamNameEnum.CHARGE_CARRIER_MASSES, opts.chargeCarrierMass)) != null) return err;

        return paramManager.isValid();
    }

    private static String setIfPresent(ParamManager paramManager, ParamNameEnum name, String value) {
        if (value == null) return null;
        Parameter p = paramManager.getParameter(name.getKey());
        if (p == null) return "Internal error: parameter not registered: -" + name.getKey();
        String err = p.parse(value);
        if (err != null) {
            return "Invalid value for parameter -" + name.getKey() + ": " + value + "\n        (" + err + ")";
        }
        p.setValueAssigned();
        return null;
    }

    private static String setIfPresent(ParamManager paramManager, ParamNameEnum name, Integer value) {
        if (value == null) return null;
        return setIfPresent(paramManager, name, value.toString());
    }

    private static String setIfPresent(ParamManager paramManager, ParamNameEnum name, Double value) {
        if (value == null) return null;
        return setIfPresent(paramManager, name, value.toString());
    }
}
