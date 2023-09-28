package org.ucb.c5.composition.checkers;

import org.ucb.c5.sequtils.HairpinCounter;

/**
 * Checks a sequence for secondary structure in windows using
 * HairpinCounter
 * 
 * It checks windows of 20 bp staggered by 10 bp for secondary structure
 * 
 * @author J. Christopher Anderson
 */
public class HairpinChecker {
    private HairpinCounter hpc;

    public void initiate() throws Exception {
        hpc = new HairpinCounter();
        hpc.initiate();
    }

    /**
     * 
     * @param inSeq DNA sequence
     * @return false if fails test
     * @throws Exception 
     */
    public boolean run(String inSeq) throws Exception {
        
        for(int start = 0; start < inSeq.length() - 20; start = start + 10) {
            int end = Math.min(inSeq.length()-1,start + 20);
            double d = hpc.run(inSeq.substring(start, end));
            if (d > 20) {
                return false;
            }
        }
        
        return true;
    }
}
