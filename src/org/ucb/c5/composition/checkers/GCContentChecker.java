package org.ucb.c5.composition.checkers;

/**
 * Checks a sequence for its GC content, ensuring it remains within
 * desired bounds (for instance, 30% to 70%).
 *
 * @author J. Christopher Anderson
 */
public class GCContentChecker {

    private double lowerBound;
    private double upperBound;

    public void initiate(double lowerBound, double upperBound) {
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
    }

    /**
     *
     * @param inSeq DNA sequence
     * @return false if GC content is outside the desired bounds
     */
    public boolean run(String inSeq) {
        int countGC = 0;

        for (char base : inSeq.toCharArray()) {
            if (base == 'G' || base == 'C') {
                countGC++;
            }
        }

        double gcContent = ((double) countGC / inSeq.length()) * 100;

        if (gcContent < lowerBound || gcContent > upperBound) {
            return false;
        }

        return true;
    }
}
