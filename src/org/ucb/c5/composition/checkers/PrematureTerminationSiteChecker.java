package org.ucb.c5.composition.checkers;

import org.ucb.c5.sequtils.RevComp;

import java.util.ArrayList;
import java.util.List;

/**
 * Checks a sequence for potential premature Termination Sites
 * I didn't end up building out this code because I struggled to find a good database with prematureTerminationSites, but left it in
 * incase somebody else wanted to do this eventually
 *
 * @author Aditya Vunnum
 */
public class PrematureTerminationSiteChecker {
    private RevComp revcomp;
    private List<String> terminators;

    public void initiate() {
        revcomp = new RevComp();
        revcomp.initiate();

        // Populate terminator sequences (just as an example, actual sequences might be different)
        terminators = new ArrayList<>();
        terminators.add("GCUGCUGCU");
        terminators.add("UGUGUGUGU");
        // ... add more terminator sequences here based on the patterns you described.
    }

    /**
     * Checks a DNA sequence for terminator Strings
     *
     * @param dnaseq
     * @return true if passes; false if contains a terminator sequence
     */
    public boolean run(String dnaseq) {
        String rc = revcomp.run(dnaseq);
        String combined = dnaseq + "x" + rc;
        combined = combined.toUpperCase();

        for (String term : terminators) {
            if (combined.contains(term)) {
//                System.out.println(term);
                return false;
            }
        }

        return true;
    }

    public static void main(String[] args) {
        PrematureTerminationSiteChecker checker = new PrematureTerminationSiteChecker();
        checker.initiate();
        boolean result = checker.run("GCUGCUGCU");  //returns false due to terminator sequence
        System.out.println(result);
    }
}

