package edu.ucsd.msjava.msgf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

public class AAFrequencyCounter {
    Histogram<String> frequencyTable;
    int nMer;
    int sizeNMer;

    public AAFrequencyCounter() {
        frequencyTable = new Histogram<String>();
        sizeNMer = 0;
    }

    public void setNMer(int nMer) {
        this.nMer = nMer;
    }

    public void readFromFreqFile(String fileName) {
        BufferedReader in = null;
        try {
            in = new BufferedReader(new FileReader(fileName));
            String s;

            s = in.readLine();
            String[] token = s.split("\t");
            assert (token[0].equalsIgnoreCase("n"));
            this.nMer = Integer.parseInt(token[1]);

            s = in.readLine();
            token = s.split("\t");
            assert (token[0].equalsIgnoreCase("size"));
            this.sizeNMer = Integer.parseInt(token[1]);

            while ((s = in.readLine()) != null) {
                token = s.split("\t");
                assert (token.length == 2);
                frequencyTable.put(token[0], Integer.parseInt(token[1]));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void readFromFasta(String fileName) {
        BufferedReader in = null;
        try {
            in = new BufferedReader(new FileReader(fileName));
            String s;
            while ((s = in.readLine()) != null) {
                if (s.startsWith(">"))
                    continue;
                StringBuffer buf = new StringBuffer();
                for (int i = 0; i < s.length(); i++) {
                    if (i >= nMer) {
                        frequencyTable.add(buf.toString());
                        sizeNMer++;
                        buf.deleteCharAt(0);
                    }
                    buf.append(s.charAt(i));
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static float getRandomFrequency(String str) {
        float uniFreq = 0.05f;
        int numLI = 0;
        for (int i = 0; i < str.length(); i++)
            if (str.charAt(i) == 'L' || str.charAt(i) == 'I')
                numLI++;
        return (float) (Math.pow(2, numLI) * Math.pow(uniFreq, str.length()));
    }

    public float getFrequency(String str) {
        ArrayList<String> strSet = new ArrayList<String>();
        strSet.add(str);
        for (int i = 0; i < str.length(); i++) {
            char c = str.charAt(i);
            if (c == 'L') {
                int size = strSet.size();
                for (int j = 0; j < size; j++) {
                    String s = strSet.get(j);
                    strSet.add(s.substring(0, i) + "I" + s.substring(i + 1));
                }
            } else if (c == 'I') {
                int size = strSet.size();
                for (int j = 0; j < size; j++) {
                    String s = strSet.get(j);
                    strSet.add(s.substring(0, i) + "L" + s.substring(i + 1));
                }
            }
        }
        int occ = 0;
        for (String s : strSet)
            occ += getOccurrence(s);
        return occ / (float) sizeNMer;
    }

    public int getOccurrence(String str) {
        Integer occ = frequencyTable.get(str);
        if (occ == null)
            return 0;
        else
            return occ;
    }

}
