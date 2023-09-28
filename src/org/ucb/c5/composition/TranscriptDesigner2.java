package org.ucb.c5.composition;

import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.composition.model.Transcript;
import org.ucb.c5.sequtils.HairpinCounter;

import java.util.*;

/**
 * This reverse translates a protein sequence to a DNA and chooses an RBS. It
 * uses the highest CAI codon for each amino acid in the specified Host.
 *
 * @author J. Christopher Anderson
 */
public class TranscriptDesigner2 {

    private Map<Character, String> aminoAcidToCodon;
    private RBSChooser rbsChooser;

    private HairpinCounter hairpinCounter;

    public void initiate() throws Exception {
        //Initialize the RBSChooser
        rbsChooser = new RBSChooser();  //Delete the 2 to use old algorithm
        rbsChooser.initiate();

        //Initialize hairpin chooser
        hairpinCounter = new HairpinCounter();
        hairpinCounter.initiate();
        
        //Construct a map between each amino acid and the highest-CAI codon for E coli
        aminoAcidToCodon = new HashMap<>();

        aminoAcidToCodon.put('A', "GCG");
        aminoAcidToCodon.put('C', "TGC");
        aminoAcidToCodon.put('D', "GAT");
        aminoAcidToCodon.put('E', "GAA");
        aminoAcidToCodon.put('F', "TTC");
        aminoAcidToCodon.put('G', "GGT");
        aminoAcidToCodon.put('H', "CAC");
        aminoAcidToCodon.put('I', "ATC");
        aminoAcidToCodon.put('K', "AAA");
        aminoAcidToCodon.put('L', "CTG");
        aminoAcidToCodon.put('M', "ATG");
        aminoAcidToCodon.put('N', "AAC");
        aminoAcidToCodon.put('P', "CCG");
        aminoAcidToCodon.put('Q', "CAG");
        aminoAcidToCodon.put('R', "CGT");
        aminoAcidToCodon.put('S', "TCT");
        aminoAcidToCodon.put('T', "ACC");
        aminoAcidToCodon.put('V', "GTT");
        aminoAcidToCodon.put('W', "TGG");
        aminoAcidToCodon.put('Y', "TAC");
    }

    public Transcript run(String peptide, Set<RBSOption> ignores) throws Exception {
        // Choose codons for each amino acid
        String[] codons = new String[peptide.length()];
        for(int i=0; i<peptide.length(); i++) {
            char aa = peptide.charAt(i);
            String codon = aminoAcidToCodon.get(aa);
            codons[i] = codon;
        }

        // Break the sequence into 3aa windows (9 nt each) and optimize them individually
        List<String> optimizedCodons = new ArrayList<>();
        int i = 0;

        while (i < codons.length - 2) { // Ensure there's room for a window
            List<String> windowOptions = new ArrayList<>();

            // Enumerate the window combined with its upstream and downstream sequences
            for (int j = 0; j < 3; j++) {
                String currentCombination = "";
                if (i > 0 && optimizedCodons.size() > 0) {
                    currentCombination += optimizedCodons.get(optimizedCodons.size() - 1); // Upstream
                }
                currentCombination += codons[i + j]; // Current window
                if (i + j + 1 < codons.length) {
                    currentCombination += codons[i + j + 1]; // Downstream
                }
                windowOptions.add(currentCombination);
            }

            // Rank each combination based on secondary structure
            String bestCombination = rankBySecondaryStructure(windowOptions);

            // Retain the middle 3 codons
            String middleThree = bestCombination.substring(3, 6); // Assuming each codon is a single char; adjust if needed
            optimizedCodons.add(middleThree);

            i += 3; // Move to the next window
        }

        // Convert optimizedCodons to an array for consistency with the original structure
        String[] finalCodons = optimizedCodons.toArray(new String[0]);

        // Choose an RBS
        StringBuilder cds = new StringBuilder();
        for(String codon : finalCodons) {
            cds.append(codon);
        }
        RBSOption selectedRBS = rbsChooser.run(cds.toString(), ignores);

        // Construct the Transcript and return it
        Transcript out = new Transcript(selectedRBS, peptide, finalCodons);
        return out;
    }

    private String rankBySecondaryStructure(List<String> windowOptions) throws Exception {
        double bestScore = Double.MAX_VALUE; // Assuming a lower deltaG score is better
        String bestSequence = null;

        for (String seq : windowOptions) {
            double currentScore = hairpinCounter.run(seq);
            if (currentScore < bestScore) {
                bestScore = currentScore;
                bestSequence = seq;
            }
        }

        return bestSequence;
    }
}
