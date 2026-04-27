package edu.ucsd.msjava.msutil;

import java.util.Comparator;

/** Generic ordered pair. */
public class Pair<A, B> {

    private A first;
    private B second;

    public Pair(A first, B second) {
        super();
        this.first = first;
        this.second = second;
    }

    public int hashCode() {
        int hashFirst = first != null ? first.hashCode() : 0;
        int hashSecond = second != null ? second.hashCode() : 0;

        return (hashFirst + hashSecond) * hashSecond + hashFirst;
    }

    public boolean equals(Object other) {
        if (other instanceof Pair<?, ?>) {
            Pair<?, ?> otherPair = (Pair<?, ?>) other;
            return
                    ((this.first == otherPair.first ||
                            (this.first != null && otherPair.first != null &&
                                    this.first.equals(otherPair.first))) &&
                            (this.second == otherPair.second ||
                                    (this.second != null && otherPair.second != null &&
                                            this.second.equals(otherPair.second))));
        }

        return false;
    }

    public String toString() {
        return "(" + first + ", " + second + ")";
    }

    public A getFirst() {
        return first;
    }

    public void setFirst(A first) {
        this.first = first;
    }

    public B getSecond() {
        return second;
    }

    public void setSecond(B second) {
        this.second = second;
    }

    public static class PairComparator<A extends Comparable<? super A>, B extends Comparable<? super B>> implements Comparator<Pair<A, B>> {
        boolean useSecondForComprison;

        public PairComparator() {
            this(false);
        }

        public PairComparator(boolean useSecondForComprison) {
            this.useSecondForComprison = useSecondForComprison;
        }

        public int compare(Pair<A, B> p1, Pair<A, B> p2) {
            if (!useSecondForComprison)
                return p1.getFirst().compareTo(p2.getFirst());
            else
                return p1.getSecond().compareTo(p2.getSecond());
        }
    }

    public static class PairReverseComparator<A extends Comparable<? super A>, B extends Comparable<? super B>> implements Comparator<Pair<A, B>> {
        boolean useSecondForComprison;

        public PairReverseComparator() {
            this(false);
        }

        public PairReverseComparator(boolean useSecondForComprison) {
            this.useSecondForComprison = useSecondForComprison;
        }

        public int compare(Pair<A, B> p1, Pair<A, B> p2) {
            if (!useSecondForComprison)
                return p2.getFirst().compareTo(p1.getFirst());
            else
                return p2.getSecond().compareTo(p1.getSecond());
        }
    }
}
