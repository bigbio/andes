package edu.ucsd.msjava.msutil;

import java.util.Comparator;

/**
 * Representation of a peak in a spectrum object.
 *
 * @author Sangtae Kim
 */
public class Peak implements Comparable<Peak> {

    private int charge = 1;
    private float mz;
    private float intensity;

    private int index = -1;
    private int rank = 151;

    public Peak(float mz, float intensity, int charge) {
        this.mz = mz;
        this.intensity = intensity;
        this.charge = charge;
    }

    public int getIndex() {
        return index;
    }

    public float getMz() {
        return mz;
    }

    /** Returns (m/z - H) * charge: the de-charged monoisotopic mass. */
    public float getMass() {
        Float monoMass = (mz - (float)Composition.ChargeCarrierMass()) * (float)charge;
        if (monoMass > 0)
            return monoMass;
        else
            return 0;
    }


    public float getIntensity() {
        return intensity;
    }

    public int getCharge() {
        return this.charge;
    }

    public Peak getShiftedPeak(float mz) {
        Peak newPeak = new Peak(mz, this.intensity, this.charge);
        newPeak.rank = this.rank;
        newPeak.index = this.index;
        return newPeak;
    }

    public void setRank(int rank) {
        this.rank = rank;
    }

    public int getRank() {
        return rank;
    }

    /**
     * Given the parent mass return the mass of the uncharged complement peak.
     * This assumes that the parent mass has no charge (H).
     *
     * @param parentMass the deprotonated and decharged parent mass
     * @return the deprotonated and decharged complement mass
     */
    public float getComplementMass(float parentMass) {
        return parentMass - getMass();
    }


    public void setIntensity(float intensity) {
        this.intensity = intensity;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public void setMz(float mz) {
        this.mz = mz;
    }

    public void setCharge(int charge) {
        this.charge = charge;
    }

    public float toUnitTolerance(float ppmTolerance) {
        return getMass() * ppmTolerance / Constants.MILLION;
    }

    /**
     * Compares this peak to another peak by mass. If the masses are equal,
     * compare by intensity.
     */
    public int compareTo(Peak p) {
        if (mz > p.mz) return 1;
        if (p.mz > mz) return -1;

        if (intensity > p.intensity) return 1;
        if (p.intensity > intensity) return -1;

        return 0;
    }


    @Override
    public int hashCode() {
        return (int) (mz + intensity + charge);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Peak)
            return equals((Peak) obj);
        return false;
    }

    public boolean equals(Peak p) {
        // this might not be a good idea for floats
        return mz == p.mz && intensity == p.intensity && charge == p.charge;
    }


    public static float getAbsoluteMassDiff(Peak p1, Peak p2) {
        return Math.abs(p1.mz - p2.mz);
    }

    @Override
    public String toString() {
        return mz + " " + intensity;
    }

    public Peak clone() {
        Peak p = new Peak(mz, intensity, charge);
        p.index = index;
        p.rank = rank;
        return p;
    }


    public static class IntensityComparator implements Comparator<Peak> {

        public int compare(Peak p1, Peak p2) {
            if (p1.intensity > p2.intensity) return 1;
            if (p2.intensity > p1.intensity) return -1;

            if (p1.mz > p2.mz) return 1;
            if (p2.mz > p1.mz) return -1;

            return 0;
        }

        public boolean equals(Peak p1, Peak p2) {
            // float exact equality intentional: these are cached values, not computed
            return p1.mz == p2.mz && p1.intensity == p2.intensity;
        }
    }

    public static class MassComparator implements Comparator<Peak> {

        public int compare(Peak p1, Peak p2) {
            return p1.compareTo(p2);
        }

        public boolean equals(Peak p1, Peak p2) {
            return p1.equals(p2);
        }

    }

    public Peak duplicate(float offset) {
        float mzOffset = offset / this.charge;
        return new Peak(mz + mzOffset, this.intensity, this.charge);
    }

}





