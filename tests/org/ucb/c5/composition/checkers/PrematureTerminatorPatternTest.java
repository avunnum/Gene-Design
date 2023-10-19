package org.ucb.c5.composition.checkers;

import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.*;

/**
 * Testing PrematureTerminatorPatternChecker based off of random and known sequences to confirm proper function.
 *
 * @author Aditya Vunnum
 */
public class PrematureTerminatorPatternTest {
    @Test
    public void testPrematureTerminationSequence() throws Exception {
        String dna = "GCGCGCGCGCGCGCGCACGTTTTT"; // replace with your actual DNA sequence
        PrematureTerminatorPatternChecker checker = new PrematureTerminatorPatternChecker();
        assertFalse(checker.run(dna));
    }

    @Test
    public void testRandomSequence() throws Exception {
        String dna = "aactagactacgaagaattcacgcagctacggaagccgctgaacagccgg";
        PrematureTerminatorPatternChecker checker = new PrematureTerminatorPatternChecker();
        assertTrue(checker.run(dna.toUpperCase()));

        String longDNA = "aaaaactgcgaaaggtcggatcagtgcccccaattgatgaacgagtgctccgcgtaacta\n" +
                "ctccaactggagcatgctactggagtgtggggcccacagggtacaaatgcgttcgcaatg\n" +
                "gcccaagtctggggagtcttctcatagtattgatattcctaagtgcacatcaccggttct\n" +
                "cctcgtgaaatcgcatttctaccgacccgagcagaagcgcatgtcggatcaaccaaatgt\n" +
                "gccacagggtcgttttaagcactagtgcgaaagtgtagcctgcgtggccatacattccac\n" +
                "gttaacattccaggacataaagggtctctccgagcaggagctcgtataaatcgaaagacg\n" +
                "gtaaggtcaattgccattactccaggcagcagctttgtcgtgaaccgcagtacctagctc\n" +
                "tccaatgcgtaaggttcgtagcgtaagaacccgagcagatcgcgggagcgggaagataga\n" +
                "gtgcggtgtgccacttagttgctcaatcctatgcattgcgattccgattgaagaagctac\n" +
                "ttcatattcccagaatatgtaactgaggagagaccaggggctccgccacacgaatatgaa\n" +
                "cggatgttacgtcaagccttgtctaactcgttagggtcttgtaatgtgcggccgaaaacc\n" +
                "tgtgacctaagtgtttgaagccattgagcggtttgcgccgcctaaggcacgaccagtatg\n" +
                "cgtatggctgactctgtaaaaaaccatatacccttctacggtggggtagttggattcttt\n" +
                "agtataccgttatgaataga";
        assertTrue(checker.run(longDNA.toUpperCase()));
    }

    @Test
    public void testKnownPrematureTerminationSites() throws Exception {
        String terminationSequence = "GCACAGGAATTGTTCTGATTAAAAAGGCGCGCGCGCGCGCGCCACTTTTCAGTTTGCTGAC";
        PrematureTerminatorPatternChecker checker = new PrematureTerminatorPatternChecker();
        assertFalse(checker.run(terminationSequence));

        String anotherTerminationSequence = "GCACTCTCGCGCGCGCGCACTTTCCCGCTTCGGTTTTTTTTTATGCGCGCATCTATCGCGCG";
        assertFalse(checker.run(anotherTerminationSequence));
    }
}


