package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.Mass;
import edu.ucsd.msjava.msutil.Matter;

import java.util.ArrayList;

public class MassListComparator<T extends Matter> {
    ArrayList<T> massList1;
    ArrayList<T> massList2;

    // massList1 and massList2 must be sorted
    public MassListComparator(ArrayList<T> massList1, ArrayList<T> massList2) {
        this.massList1 = massList1;
        this.massList2 = massList2;
    }

    public MatchedPair[] getMatchedList(Tolerance tolerance) {
        int i1 = 0, i2 = 0;
        ArrayList<MatchedPair> matches = new ArrayList<MatchedPair>();

        float m1, m2;
        while (i1 < massList1.size() && i2 < massList2.size()) {
            m1 = massList1.get(i1).getMass();
            m2 = massList2.get(i2).getMass();
            float tol = tolerance.getToleranceAsDa(m1);
            if (m2 <= m1 - tol) {
                i2++;
                continue;
            }
            // m2 > m1-tolerance
            if (m2 < m1 + tol) {
                matches.add(new MatchedPair<T>(massList1.get(i1), massList1.get(i2)));
                if (i1 == massList1.size() - 1)
                    i2++;
                else if (i2 == massList2.size() - 1)
                    i1++;
                else {
                    if (massList1.get(i1 + 1).getMass() < massList2.get(i2 + 1).getMass())
                        i1++;
                    else
                        i2++;
                }
            } else    // m2 >= m1+tolerance
            {
                i1++;
            }
        }
        return matches.toArray(new MatchedPair[0]);
    }


    public static class MatchedPair<T extends Matter> {
        T m1, m2;

        public MatchedPair(T m1, T m2) {
            this.m1 = m1;
            this.m2 = m2;
        }

        public T getMass1() {
            return m1;
        }

        public T getMass2() {
            return m2;
        }
    }

}
