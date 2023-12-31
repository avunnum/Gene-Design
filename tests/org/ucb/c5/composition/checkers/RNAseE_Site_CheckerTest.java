package org.ucb.c5.composition.checkers;

import org.junit.BeforeClass;
import org.junit.Test;

public class RNAseE_Site_CheckerTest {

    private static RNAseE_Site_Checker checker;

    @BeforeClass
    public static void setUpClass() throws Exception {

        checker = new RNAseE_Site_Checker();
        checker.initiate();
    }

    /**
     * Inputs known sequences and checks for RNAseE Sites. False if match is present.
     * Written By Naveen Kumaran
     * @throws Exception
     */
    
    //JCA:  I silenced this test.  I think this Promoter seq does have a site

//    @Test
//    public void CaseSensitivityTest() throws Exception {
//        //Tests whether checker only looks at Upper Case Letters.
//        //Promoter region obtained from CompositionToDNA.
//        String Promoter = "ttatgacaacttgacggctacatcattcactttttcttcacaaccggcacggaactcgctcgggctggccccggtgcattttttaaatacccgcgagaaatagagttgatcgtcaaaaccaacattgcgaccgacggtggcgataggcatccgggtggtgctcaaaagcagcttcgcctggctgatacgttggtcctcgcgccagcttaagacgctaatccctaactgctggcggaaaagatgtgacagacgcgacggcgacaagcaaacatgctgtgcgacgctggcgatatcaaaattgctgtctgccaggtgatcgctgatgtactgacaagcctcgcgtacccgattatccatcggtggatggagcgactcgttaatcgcttccatgcgccgcagtaacaattgctcaagcagatttatcgccagcagctccgaatagcgcccttccccttgcccggcgttaatgatttgcccaaacaggtcgctgaaatgcggctggtgcgcttcatccgggcgaaagaaccccgtattggcaaatattgacggccagttaagccattcatgccagtaggcgcgcggacgaaagtaaacccactggtgataccattcgcgagcctccggatgacgaccgtagtgatgaatctctcctggcgggaacagcaaaatatcacccggtcggcaaacaaattctcgtccctgatttttcaccaccccctgaccgcgaatggtgagattgagaatataacctttcattcccagcggtcggtcgataaaaaaatcgagataaccgttggcctcaatcggcgttaaacccgccaccagatgggcattaaacgagtatcccggcagcaggggatcattttgcgcttcagccatacttttcatactcccgccattcagagaagaaaccaattgtccatattgcatcagacattgccgtcactgcgtcttttactggctcttctcgctaaccaaaccggtaaccccgcttattaaaagcattctgtaacaaagcgggaccaaagccatgacaaaaacgcgtaacaaaagtgtctataatcacggcagaaaagtccacattgattatttgcacggcgtcacactttgctatgccatagcatttttatccataagattagcggatcctacctgacgctttttatcgcaactctctactgtttctccatacccgtttttttgggctagc";
//        boolean result = checker.run(Promoter);
//        assert(result == true);
//    }

    //JCA:  I silenced this test.  I think this Promoter seq does have a site
    
//    @Test
//    public void PositiveSequenceTest() throws Exception {
//        //Inputs a known sequence that DOESN'T contain any RNAseE binding sites.
//        //Should return positive.
//        String DNA = "ttatgacaacttgacggctacatcattcactttttcttcacaaccggcacggaactcgctcgggctggccccggtgcattttttaaatacccgcgagaaatagagttgatcgtcaaaaccaacattgcgaccgacggtggcgataggcatccgggtggtgctcaaaagcagcttcgcctggctgatacgttggtcctcgcgccagcttaagacgctaatccctaactgctggcggaaaagatgtgacagacgcgacggcgacaagcaaacatgctgtgcgacgctggcgatatcaaaattgctgtctgccaggtgatcgctgatgtactgacaagcctcgcgtacccgattatccatcggtggatggagcgactcgttaatcgcttccatgcgccgcagtaacaattgctcaagcagatttatcgccagcagctccgaatagcgcccttccccttgcccggcgttaatgatttgcccaaacaggtcgctgaaatgcggctggtgcgcttcatccgggcgaaagaaccccgtattggcaaatattgacggccagttaagccattcatgccagtaggcgcgcggacgaaagtaaacccactggtgataccattcgcgagcctccggatgacgaccgtagtgatgaatctctcctggcgggaacagcaaaatatcacccggtcggcaaacaaattctcgtccctgatttttcaccaccccctgaccgcgaatggtgagattgagaatataacctttcattcccagcggtcggtcgataaaaaaatcgagataaccgttggcctcaatcggcgttaaacccgccaccagatgggcattaaacgagtatcccggcagcaggggatcattttgcgcttcagccatacttttcatactcccgccattcagagaagaaaccaattgtccatattgcatcagacattgccgtcactgcgtcttttactggctcttctcgctaaccaaaccggtaaccccgcttattaaaagcattctgtaacaaagcgggaccaaagccatgacaaaaacgcgtaacaaaagtgtctataatcacggcagaaaagtccacattgattatttgcacggcgtcacactttgctatgccatagcatttttatccataagattagcggatcctacctgacgctttttatcgcaactctctactgtttctccatacccgtttttttgggctagcacgaatccgactctcacattatcaggggtataaaaATGGCGCTGTTGTCAAGTTCACTTAGCAGCCAAATCCCAACCGGTTCCCACCCCCTCACGCACACACAGTGTATACCGCATTTCTCGACCACTATAAACGCCGGTATCAGCGCCGGTAAGCCTCGATCTTTTTACTTGAGGTGGGGCAAAGGTTCCAACAAGATTATTGCATGTGTCGGGGAGGGGACTACTTCTCTGCCTTATCAAAGCGCTGAGAAGACCGATAGTTTGTCAGCTCCCACACTCGTGAAAAGGGAGTTTCCGCCAGGGTTCTGGAAAGATCATGTTATCGATTCACTGACGAGTAGCCACAAGGTCAGCGCGGCAGAGGAAAAACGGATGGAAACATTAATAAGCGAAATCAAGAATATCTTTCGGTCAATGGGCTATGGGGAAACCAATCCGTCGGCGTACGATACAGCTTGGGTCGCTCGCATACCAGCTGTCGATGGATCAGAACACCCCGAGTTCCCAGAGACATTAGAGTGGATTCTCCAAAATCAGCTCAAAGATGGTAGCTGGGGGGAAGGCTTCTACTTCCTTGCCTACGACAGGATTCTGGCAACGTTGGCGTGCATAATAACGTTAACTCTTTGGCGGACCGGAGAAACCCAAATAAGGAAAGGCATTGAGTTCTTTAAAACTCAGGCCGGAAAGATTGAAGATGAAGCCGATTCTCACAGGCCGTCTGGGTTCGAGATCGTCTTCCCGGCGATGTTAAAAGAGGCTAAAGTTCTGGGGTTAGACTTACCCTACGAATTGCCATTTATAAAGCAGATCATCGAAAAACGGGAGGCTAAGTTGGAACGATTACCTACTAATATATTGTACGCCCTGCCAACGACGCTACTCTACTCACTTGAGGGATTGCAGGAGATCGTTGATTGGGAAAAGATAATTAAACTACAATCCAAGGACGGAAGTTTCTTAACGAGCCCCGCTTCGACGGCCGCGGTTTTCATGCGGACTGGGAATAAGAAGTGCCTTGAGTTCCTTAACTTCGTCCTTAAGAAATTCGGAAACCACGTTCCATGTCACTACCCCCTTGACTTGTTCGAACGTTTATGGGCGGTAGACACGGTTGAGCGGTTGGGCATCGACCACCATTTTAAGGAGGAAATTAAGGATGCACTCGATTACGTATATTCGCACTGGGATGAACGCGGGATTGGGTGGGCCAGAGAGAACCCTATTCCAGATATTGACGACACCGCAATGGGCTTAAGAATATTACGGCTGCACGGCTATAACGTGAGCTCAGACGTGTTGAAGACTTTTCGTGACGAGAATGGGGAGTTTTTCTGCTTCCTCGGGCAAACACAAAGAGGCGTCACAGATATGTTAAACGTTAACCGCTGTTCGCATGTGGCGTTCCCGGGTGAAACAATCATGCAGGAAGCTAAATTATGCACGGAGCGATATCTTCGAAACGCACTAGAGGATGTAGGTGCGTTCGACAAATGGGCTCTGAAAAAGAACATCCGGGGCGAAGTAGAGTATGCATTAAAGTACCCGTGGCACCGCTCCATGCCTAGGTTAGAAGCTCGTAGTTATATAGAGCACTATGGGCCTAACGACGTGTGGCTGGGGAAAACTATGTACATGATGCCTTATATTAGTAACCTGAAATATCTTGAGTTGGCTAAGCTGGATTTCAATCATGTCCAGTCGCTACACCAAAAGGAGCTTAGAGACTTACGTCGATGGTGGAAATCGTCGGGTCTGTCAGAGCTGAAGTTTACAAGAGAACGCGTGACTGAGATTTACTTCTCGGCGGCCTCTTTCATATTTGAACCGGAATTTGCTACTTGTAGGGATGTCTACACAAAAATCAGTATATTTACGGTGATTCTGGACGATCTCTATGATGCCCATGGGACTCTGGACAACCTTGAGCTTTTCAGTGAGGGTGTCAAACGGTGGGATTTATCGTTAGTAGATCGTATGCCGCAGGATATGAAAATTTGCTTCACAGTTCTGTACAATACGGTTAACGAGATCGCTGTGGAAGGCAGAAAACGCCAAGGGCGCGATGTACTCGGATACATTAGAAACGTCCTCGAAATCTTACTCGCGGCACATACCAAGGAGGCAGAATGGTCGGCCGCTCGTTACGTTCCTTCCTTTGATGAATATATTGAGAACGCGTCCGTTTCGATCTCTCTCGGAACGTTGGTTCTCATATCTGTTCTCTTCACTGGGGAAATTCTCACCGATGATGTATTATCTAAGATTGGCCGTGGATCGAGGTTTCTACAGCTTATGGGACTGACAGGGCGCCTCGTTAATGACACTAAGACCTATGAAGCAGAAAGAGGTCAGGGCGAGGTAGCGTCGGCAGTCCAATGCTATATGAAGGAACATCCCGAGATATCCGAGGAGGAAGCGTTAAAACATGTGTATACAGTTATGGAAAATGCCCTAGACGAGTTGAACAGAGAATTTGTCAACAACCGTGATGTCCCTGACTCGTGCAGGCGCTTGGTTTTTGAAACAGCTAGGATAATGCAGTTATTTTATATGGAGGGGGACGGGCTGACGTTAAGTCACGAAATGGAAATAAAAGAGCACGTAAAGAATTGTCTATTTCAGCCGGTGGCTgttgaccgagcactgtgattttttgaggtaacaagATGTATCCCTTCATACGCACGGCCAGAATGACGGTCTGCGCGAAGAAACACGTGCACCTTACAAGGGACGCGGCAGAGCAACTGCTCGCCGATATCGACAGACGCCTGGATCAACTCTTACCAGTGGAGGGTGAACGTGACGTTGTCGGGGCAGCCATGCGGGAAGGGGCTCTAGCTCCGGGCAAAAGAATAAGGCCCATGCTATTACTTCTCACTGCGCGTGATCTGGGGTGTGCAGTCTCTCATGATGGCCTATTGGACCTTGCCTGCGCGGTTGAGATGGTTCACGCCGCGAGTTTAATCCTTGATGACATGCCCTGTATGGATGACGCTAAACTCAGGCGCGGGCGGCCTACCATTCATAGCCACTATGGTGAGCACGTGGCGATCTTAGCTGCCGTCGCACTACTTTCGAAGGCGTTCGGTGTCATAGCTGATGCGGACGGACTTACGCCATTAGCTAAGAACCGGGCAGTATCGGAATTGAGTAATGCGATCGGTATGCAGGGGCTTGTCCAAGGCCAGTTCAAGGACCTGAGCGAGGGAGACAAACCTAGATCCGCTGAAGCGATCCTCATGACCAACCATTTCAAAACATCGACGCTCTTCTGTGCCTCTATGCAGATGGCTTCAATAGTGGCAAATGCGTCTTCTGAAGCGCGCGACTGCTTGCACCGATTCTCCTTGGACTTGGGTCAGGCCTTCCAACTTTTAGATGACCTCACAGACGGTATGACCGACACTGGGAAGGACAGTAACCAGGATGCAGGGAAATCAACACTTGTGAACCTATTAGGTCCCAGGGCGGTTGAAGAACGGCTTCGGCAACACCTCCAGCTCGCCTCTGAACACCTCTCCGCGGCTTGCCAACATGGCCATGCTACCCAACATTTTATCCAAGCCTGGTTCGACAAGAAACTTGCCGCAGTTAGCTGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTT";
//        boolean result = checker.run(DNA);
//        assert(result == true);
//    }

    @Test
    public void AnotherPositiveSequenceTest() throws Exception {
       //Inputs blank DNA sequence which should return RNAseEsite free.
        // Checks the used algorith to see if size of string affects results.
        String DNA = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        boolean result = checker.run(DNA);
        assert(result == true);
    }

    @Test
    public void NegativeSequenceTest() throws Exception {
        //Inputs a known sequence that does contain RNAse Binding sites.
        // Should return as false.
        String DNA = "ttatgacaacttgacggctacatcattcactttttcttcacaaccggcacggaactcgctcgggctggccccggtgcattttttaaatacccgcgagaaatagagttgatcgtcaaaaccaacattgcgaccgacggtggcgataggcatccgggtggtgctcaaaagcagcttcgcctggctgatacgttggtcctcgcgccagcttaagacgctaatccctaactgctggcggaaaagatgtgacagacgcgacggcgacaagcaaacatgctgtgcgacgctggcgatatcaaaattgctgtctgccaggtgatcgctgatgtactgacaagcctcgcgtacccgattatccatcggtggatggagcgactcgttaatcgcttccatgcgccgcagtaacaattgctcaagcagatttatcgccagcagctccgaatagcgcccttccccttgcccggcgttaatgatttgcccaaacaggtcgctgaaatgcggctggtgcgcttcatccgggcgaaagaaccccgtattggcaaatattgacggccagttaagccattcatgccagtaggcgcgcggacgaaagtaaacccactggtgataccattcgcgagcctccggatgacgaccgtagtgatgaatctctcctggcgggaacagcaaaatatcacccggtcggcaaacaaattctcgtccctgatttttcaccaccccctgaccgcgaatggtgagattgagaatataacctttcattcccagcggtcggtcgataaaaaaatcgagataaccgttggcctcaatcggcgttaaacccgccaccagatgggcattaaacgagtatcccggcagcaggggatcattttgcgcttcagccatacttttcatactcccgccattcagagaagaaaccaattgtccatattgcatcagacattgccgtcactgcgtcttttactggctcttctcgctaaccaaaccggtaaccccgcttattaaaagcattctgtaacaaagcgggaccaaagccatgacaaaaacgcgtaacaaaagtgtctataatcacggcagaaaagtccacattgattatttgcacggcgtcacactttgctatgccatagcatttttatccataagattagcggatcctacctgacgctttttatcgcaactctctactgtttctccatacccgtttttttgggctagcacgaatccgactctcacattatcaggggtataaaaATGGCGCTGTTGTCAAGTTCACTTAGCAGCCAAATCCCAACCGGTTCCCACCCCCTCACGCACACACAGTGTATACCGCATTTCTCGACCACTATAAACGCCGGTATCAGCGCCGGTAAGCCTCGATCTTTTTACTTGAGGTGGGGCAAAGGTTCCAACAAGATTATTGCATGTGTCGGGGAGGGGACTACTTCTCTGCCTTATCAAAGCGCTGAGAAGACCGATAGTTTGTCAGCTCCCACACTCGTGAAAAGGGAGTTTCCGCCAGGGTTCTGGAAAGATCATGTTATCGATTCACTGACGAGTAGCCACAAGGTCAGCGCGGCAGAGGAAAAACGGATGGAAACATTAATAAGCGAAATCAAGAATATCTTTCGGTCAATGGGCTATGGGGAAACCAATCCGTCGGCGTACGATACAGCTTGGGTCGCTCGCATACCAGCTGTCGATGGATCAGAACACCCCGAGTTCCCAGAGACATTAGAGTGGATTCTCCAAAATCAGCTCAAAGATGGTAGCTGGGGGGAAGGCTTCTACTTCCTTGCCTACGACAGGATTCTGGCAACGTTGGCGTGCATAATAACGTTAACTCTTTGGCGGACCGGAGAAACCCAAATAAGGAAAGGCATTGAGTTCTTTAAAACTCAGGCCGGAAAGATTGAAGATGAAGCCGATTCTCACAGGCCGTCTGGGTTCGAGATCGTCTTCCCGGCGATGTTAAAAGAGGCTAAAGTTCTGGGGTTTTTTTTACCCTACGAATTGCCATTTATAAAGCAGATCATCGAAAAACGGGAGGCTAAGTTGGAACGATTACCTACTAATATATTGTACGCCCTGCCAACGACGCTACTCTACTCACTTGAGGGATTGCAGGAGATCGTTGATTGGGAAAAGATAATTAAACTACAATCCAAGGACGGAAGTTTCTTAACGAGCCCCGCTTCGACGGCCGCGGTTTTCATGCGGACTGGGAATAAGAAGTGCCTTGAGTTCCTTAACTTCGTCCTTAAGAAATTCGGAAACCACGTTCCATGTCACTACCCCCTTGACTTGTTCGAACGTTTATGGGCGGTAGACACGGTTGAGCGGTTGGGCATCGACCACCATTTTAAGGAGGAAATTAAGGATGCACTCGATTACGTATATTCGCACTGGGATGAACGCGGGATTGGGTGGGCCAGAGAGAACCCTATTCCAGATATTGACGACACCGCAATGGGCTTAAGAATATTACGGCTGCACGGCTATAACGTGAGCTCAGACGTGTTGAAGACTTTTCGTGACGAGAATGGGGAGTTTTTCTGCTTCCTCGGGCAAACACAAAGAGGCGTCACAGATATGTTAAACGTTAACCGCTGTTCGCATGTGGCGTTCCCGGGTGAAACAATCATGCAGGAAGCTAAATTATGCACGGAGCGATATCTTCGAAACGCACTAGAGGATGTAGGTGCGTTCGACAAATGGGCTCTGAAAAAGAACATCCGGGGCGAAGTAGAGTATGCATTAAAGTACCCGTGGCACCGCTCCATGCCTAGGTTAGAAGCTCGTAGTTATATAGAGCACTATGGGCCTAACGACGTGTGGCTGGGGAAAACTATGTACATGATGCCTTATATTAGTAACCTGAAATATCTTGAGTTGGCTAAGCTGGATTTCAATCATGTCCAGTCGCTACACCAAAAGGAGCTTAGAGACTTACGTCGATGGTGGAAATCGTCGGGTCTGTCAGAGCTGAAGTTTACAAGAGAACGCGTGACTGAGATTTACTTCTCGGCGGCCTCTTTCATATTTGAACCGGAATTTGCTACTTGTAGGGATGTCTACACAAAAATCAGTATATTTACGGTGATTCTGGACGATCTCTATGATGCCCATGGGACTCTGGACAACCTTGAGCTTTTCAGTGAGGGTGTCAAACGGTGGGATTTATCGTTAGTAGATCGTATGCCGCAGGATATGAAAATTTGCTTCACAGTTCTGTACAATACGGTTAACGAGATCGCTGTGGAAGGCAGAAAACGCCAAGGGCGCGATGTACTCGGATACATTAGAAACGTCCTCGAAATCTTACTCGCGGCACATACCAAGGAGGCAGAATGGTCGGCCGCTCGTTACGTTCCTTCCTTTGATGAATATATTGAGAACGCGTCCGTTTCGATCTCTCTCGGAACGTTGGTTCTCATATCTGTTCTCTTCACTGGGGAAATTCTCACCGATGATGTATTATCTAAGATTGGCCGTGGATCGAGGTTTCTACAGCTTATGGGACTGACAGGGCGCCTCGTTAATGACACTAAGACCTATGAAGCAGAAAGAGGTCAGGGCGAGGTAGCGTCGGCAGTCCAATGCTATATGAAGGAACATCCCGAGATATCCGAGGAGGAAGCGTTAAAACATGTGTATACAGTTATGGAAAATGCCCTAGACGAGTTGAACAGAGAATTTGTCAACAACCGTGATGTCCCTGACTCGTGCAGGCGCTTGGTTTTTGAAACAGCTAGGATAATGCAGTTATTTTATATGGAGGGGGACGGGCTGACGTTAAGTCACGAAATGGAAATAAAAGAGCACGTAAAGAATTGTCTATTTCAGCCGGTGGCTgttgaccgagcactgtgattttttgaggtaacaagATGTATCCCTTCATACGCACGGCCAGAATGACGGTCTGCGCGAAGAAACACGTGCACCTTACAAGGGACGCGGCAGAGCAACTGCTCGCCGATATCGACAGACGCCTGGATCAACTCTTACCAGTGGAGGGTGAACGTGACGTTGTCGGGGCAGCCATGCGGGAAGGGGCTCTAGCTCCGGGCAAAAGAATAAGGCCCATGCTATTACTTCTCACTGCGCGTGATCTGGGGTGTGCAGTCTCTCATGATGGCCTATTGGACCTTGCCTGCGCGGTTGAGATGGTTCACGCCGCGAGTTTAATCCTTGATGACATGCCCTGTATGGATGACGCTAAACTCAGGCGCGGGCGGCCTACCATTCATAGCCACTATGGTGAGCACGTGGCGATCTTAGCTGCCGTCGCACTACTTTCGAAGGCGTTCGGTGTCATAGCTGATGCGGACGGACTTACGCCATTAGCTAAGAACCGGGCAGTATCGGAATTGAGTAATGCGATCGGTATGCAGGGGCTTGTCCAAGGCCAGTTCAAGGACCTGAGCGAGGGAGACAAACCTAGATCCGCTGAAGCGATCCTCATGACCAACCATTTCAAAACATCGACGCTCTTCTGTGCCTCTATGCAGATGGCTTCAATAGTGGCAAATGCGTCTTCTGAAGCGCGCGACTGCTTGCACCGATTCTCCTTGGACTTGGGTCAGGCCTTCCAACTTTTAGATGACCTCACAGACGGTATGACCGACACTGGGAAGGACAGTAACCAGGATGCAGGGAAATCAACACTTGTGAACCTATTAGGTCCCAGGGCGGTTGAAGAACGGCTTCGGCAACACCTCCAGCTCGCCTCTGAACACCTCTCCGCGGCTTGCCAACATGGCCATGCTACCCAACATTTTATCCAAGCCTGGTTCGACAAGAAACTTGCCGCAGTTAGCTGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTT";
        boolean result = checker.run(DNA);
        assert(result == false);
    }

    @Test
    public void AnotherNegativeSequenceTest() throws Exception {
        //another negative test just using a smaller snipet of DNA.
        //Determines whether size of string plays an issue.
        String DNA = "ACGATTTTACAAGTTTTG";
        boolean result = checker.run(DNA);
        assert(result == false);
    }
}
