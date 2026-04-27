package edu.ucsd.msjava.msutil;


/** Root class for anything that has a mass. */
public abstract class Matter implements Comparable<Matter> {

    public abstract float getMass();

    public double getAccurateMass() {
        return getMass();
    }

    public abstract int getNominalMass();

    public int compareTo(Matter other) {
        if (this.getMass() > other.getMass()) return 1;
        if (other.getMass() > this.getMass()) return -1;
        return 0;
    }

    public String toString() {
        return String.format("[%.2f]", getMass());
    }

    public abstract boolean equals(Object obj);
}
