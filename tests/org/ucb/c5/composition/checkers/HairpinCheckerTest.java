package org.ucb.c5.composition.checkers;

import org.junit.Test;

import static org.junit.Assert.*;

public class HairpinCheckerTest {
    @Test
    public void testTetraloops() throws Exception {
        HairpinChecker hairpinChecker = new HairpinChecker();
        hairpinChecker.initiate();
        // UUCG is very stable
        // UNCG - stable. UUCG (part of this group) is very stable
        // UNAC - stable
        // GNRA - stable
        // CUYG - somewhat stable
        // ANYA - somewhat stable
        boolean b = hairpinChecker.run("AAAAAAAAAAAAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("AAAAAAATTCGAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("AATTCGAAAATTCGAA");
        assertFalse(b);
        b = hairpinChecker.run("AAAAAAATCCGAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("ATCCGAAATCCGAA");
        assertFalse(b);
        b = hairpinChecker.run("AAAAAAATACGAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("ATACGAAATACGAA");
        assertFalse(b);
        b = hairpinChecker.run("AAAAAAATGCGAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("AAAAAAATTACAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("AAAAAAATAACAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("AAAAAAATGACAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("AAAAAAATCACAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("CCCCCCCGTAACCCCCCCCCCCC");
        assertTrue(b);
        b = hairpinChecker.run("AAAAAAAAAAACTTGAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("AAAAAAAAAAACTCGAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("AAAAAAAAAAAATTAAAAAAAAA");
        assertTrue(b);
    }

    @Test
    public void testPseudoknots() throws Exception {
        HairpinChecker hairpinChecker = new HairpinChecker();
        hairpinChecker.initiate();
        boolean b = hairpinChecker.run("AAAAAAGGGGGGTTTTTTGGGGGG");
        assertTrue(b);
        b = hairpinChecker.run("AAAAAAGGGGGGTTTTTTCCCCCC");
        assertFalse(b);

        // A real life pseudoknot
        b = hairpinChecker.run("TTTAAACTGATACGGGGTATCAGTCAGGCTCGGCTGGTACCCCTTGCAAAGCGAGCC");
        assertFalse(b);
        b = hairpinChecker.run("TTTAAACTGATACAAAATATCAGTCAGGCTCGGCTGGTACCCCTTGCAAAGCGAGCC");
        assertTrue(b);
    }
    
    @Test
    public void testSimple() throws Exception {
        HairpinChecker hairpinChecker = new HairpinChecker();
        hairpinChecker.initiate();
        boolean b = hairpinChecker.run("CCCAAAAAAAGGGAAAAAAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("AAACGCGAAAAAACGCGAAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("AAACACGAAAAAACGTGAAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("AAACCCCCAAAAAGGGGGAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("AAACCCCCAAAAAGGGGGAAA");
        assertFalse(b);

        b = hairpinChecker.run("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        assertTrue(b);
        b = hairpinChecker.run("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
        assertTrue(b);
        b = hairpinChecker.run("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
        assertTrue(b);
        b = hairpinChecker.run("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        assertTrue(b);
    }
}
