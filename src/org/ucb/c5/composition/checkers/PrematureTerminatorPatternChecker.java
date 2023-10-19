package org.ucb.c5.composition.checkers;
import java.util.regex.*;

/**
 * Checks to see if there are any premature termination sites
 * Generalized Regex pattern from a much more complicated analysis of premature termination sites.
 * This code certainly doesn't catch all, but will catch basic premature termination sites
 * Papers referenced:
 *
 *  Yves d'Aubenton Carafa, Edward Brody, Claude Thermes (1990).
 *  "Prediction of rho-independent Escherichia coli transcription terminators: A statistical analysis of their RNA stem-loop structures".
 *  Journal of Molecular Biology, Volume 216, Issue 4, Pages 835-858.
 *  DOI: 10.1016/S0022-2836(99)80005-9
 *
 *  Jason M. Peters, Abbey D. Vangeloff, and Robert Landick (2011).
 *  "Bacterial Transcription Terminators: The RNA 3′-End Chronicles".
 *  Journal of Molecular Biology, Volume 412, Issue 5, Pages 793-813.
 *  Published online 2011 Mar 23.
 *  DOI: 10.1016/j.jmb.2011.03.036
 *  PMCID: PMC3622210
 *  NIHMSID: NIHMS284190
 *  PMID: 21439297
 *
 * @author Aditya Vunnum
 */
public class PrematureTerminatorPatternChecker {
    // Regular expressions for the various terminator parts
    private static final String HAIRPIN_STEM = "(GC){5,8}"; // 5 to 8 GC pairs
    private static final String LOOP = "[ACGT]{3,4}"; // 3 or 4 nucleotides
    private static final String T_TRACT = "T{2,}"; // At least 2 Ts

    private static final int WINDOW_SIZE = 55; // size of the window
    private static final int SLIDE_AMOUNT = 37; // amount to slide the window each iteration

    public boolean run(String dnaSequence) {
        int sequenceLength = dnaSequence.length();
        int startIndex = 0;

        if (sequenceLength < WINDOW_SIZE) {
            return !checkTerminatorInWindow(dnaSequence);
        }

        while (startIndex + WINDOW_SIZE <= sequenceLength) {
            String window = dnaSequence.substring(startIndex, startIndex + WINDOW_SIZE);
            if (checkTerminatorInWindow(window)) {
                return false; // Found a premature terminator, return false
            }
            startIndex += SLIDE_AMOUNT;
        }
        return true; // No premature terminator found, return true
    }

    private boolean checkTerminatorInWindow(String window) {
        String patternString = HAIRPIN_STEM + LOOP + T_TRACT;
        Pattern pattern = Pattern.compile(patternString);
        Matcher matcher = pattern.matcher(window);

        return matcher.find();
    }

    public static void main(String[] args) {
        String dna = "GCGCGCGCGCGCGCGCACGTTTTT"; // replace with your actual DNA sequence
        PrematureTerminatorPatternChecker checker = new PrematureTerminatorPatternChecker();
        System.out.println(checker.run(dna));
    }
}
