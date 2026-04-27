package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.Matter;

public record ProfilePeak<T extends Matter>(T node, float probability) implements Comparable<ProfilePeak<T>> {

    public T getNode() { return node; }
    public float getProbability() { return probability; }

    @Override
    public int compareTo(ProfilePeak<T> p) {
        return node.compareTo(p.node);
    }
}
