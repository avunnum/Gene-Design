package org.ucb.c5.composition;

import org.ucb.c5.composition.checkers.*;
import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.composition.model.Transcript;
import org.ucb.c5.sequtils.HairpinCounter;
import org.ucb.c5.utils.FileUtils;

import java.util.*;

/**
 * This reverse translates a protein sequence to a DNA and chooses an RBS.
 * It uses guided random monte carlo to generate random sequences within sliding windows iteratively
 * Checks for forbidden sequences, secondary structure, repetitive sequences, high GC content, and potential termination sites
 *
 * @author J. Christopher Anderson
 * @author Aditya Vunnum
 */
public class TranscriptDesigner2 {

    private RBSChooser2 rbsChooser2;
    private HairpinChecker hairpinChecker;
    private GCContentChecker gcContentChecker;
    private ForbiddenSequenceChecker forbiddenSequenceChecker;
    private RepeatSequenceChecker repeatSequenceChecker;
    private PrematureTerminatorPatternChecker prematureTerminatorPatternChecker;
    public Map<String, List<String>> codon_usage;
    private int seed = 56372;
    private Random rand;

    public void initiate() throws Exception {
        // Reading in codon usage file
        String codon_usage_file = FileUtils.readResourceFile("composition/data/codon_usage.txt");
        String[] lines = codon_usage_file.split("\\r|\\r?\\n");

        // Maps each amino acid to all the codons that code for it
        codon_usage = new HashMap<>();

        // maps amino acids to all codons, should add up to 100.
        for (int i = 0; i < lines.length; i++) {
            String line = lines[i];
            if(line.isEmpty()) {
                continue;
            }
            String[] tabs = line.split("\t");
            String codon = tabs[0];
            String aminoAcid = tabs[1];
            int frequency  = (int) Math.round(Double.parseDouble(tabs[2]) * 100);
            codon_usage.putIfAbsent(aminoAcid, new ArrayList<>());
            for (int j = 0; j < frequency; j++) {
                codon_usage.get(aminoAcid).add(codon);
            }
        }
        // I and R both were at 99, so added to account for edge case of choosing 99
        codon_usage.get("I").add("ATA");
        codon_usage.get("R").add("AGG");

        rbsChooser2 = new RBSChooser2();  //Delete the 2 to use old algorithm
        rbsChooser2.initiate();

        hairpinChecker = new HairpinChecker();
        hairpinChecker.initiate();

        gcContentChecker = new GCContentChecker();
        gcContentChecker.initiate(30.0, 70.0);

        forbiddenSequenceChecker = new ForbiddenSequenceChecker();
        forbiddenSequenceChecker.initiate();

        repeatSequenceChecker = new RepeatSequenceChecker();
        prematureTerminatorPatternChecker = new PrematureTerminatorPatternChecker();

    }

    public Transcript run(String peptide, Set<RBSOption> ignores) throws Exception {
        boolean entireSequenceChecker = true;
        int numTries = 0;
        rand = new Random(seed);
        String testedCDS = null;
        while (entireSequenceChecker) {
            StringBuilder optimizedCDS = new StringBuilder();

            for (int i = 0; i < peptide.length(); i += 3) {
                String currentWindow = getCurrentWindow(peptide, i);
                String preamble = getPreamble(optimizedCDS.toString());
                String downstream = getDownstream(peptide, i + 3);

                String optimizedWindow = optimizeWindow(currentWindow, preamble, downstream);
                optimizedCDS.append(optimizedWindow);
            }
            numTries++;

            // checking the entire cds rather than windows to look for problems at the junctions between windows
            if ((forbiddenSequenceChecker.run(optimizedCDS.toString()) && hairpinChecker.run(optimizedCDS.toString())
                    && prematureTerminatorPatternChecker.run(optimizedCDS.toString())) || numTries == 100) {
                entireSequenceChecker = false;
                testedCDS = optimizedCDS.toString();
            }
        }

        RBSOption selectedRBS = rbsChooser2.run(testedCDS, ignores);

        //changing cds to an array of codons
        String[] codons = stringToCodon(testedCDS);

        //Construct the Transcript and return it
        Transcript out = new Transcript(selectedRBS, peptide, codons);
        return out;
    }

    private String getCurrentWindow(String peptide, int start) {
        if (start + 3 > peptide.length()) {
            return peptide.substring(start);
        }
        return peptide.substring(start, start + 3);
    }

    private String getPreamble(String optimizedCDS) {
        int length = optimizedCDS.length();
        return (length >= 9) ? optimizedCDS.substring(length - 9, length) : "";
    }

    private String getDownstream(String peptide, int downstreamStart) {
        // returns the 6 amino acids downstream
        if (downstreamStart + 5 >= peptide.length()) return "";
        String downstreamPeptide = peptide.substring(downstreamStart, downstreamStart + 6);
        return downstreamPeptide;
    }

    private String optimizeWindow(String window, String preamble, String downstream) throws Exception {
        List<String> guidedRandomSequences = generateGuidedRandomSequences(window + downstream);
        for (String seq : guidedRandomSequences) {
            String testSeq = preamble + seq;
            if (forbiddenSequenceChecker.run(testSeq) && hairpinChecker.run(testSeq) && gcContentChecker.run(testSeq) && repeatSequenceChecker.run(testSeq)) {
                // return the first 3 codons, which is the current window
                if (seq.length() < 9) {
                    return seq;
                }
                return seq.substring(0, 9);
            }
        }
        // just return last sequence otherwise
        if (guidedRandomSequences.get(guidedRandomSequences.size() - 1).length() < 9) {
            return guidedRandomSequences.get(guidedRandomSequences.size() - 1);
        } else {
            return guidedRandomSequences.get(guidedRandomSequences.size() - 1).substring(0, 9);
        }
    }

    private List<String> generateGuidedRandomSequences(String windowAndDownstream) {
        List<String> sequences = new ArrayList<>();

        for (int j = 0; j < 30; j++) {
            StringBuilder currentSequence = new StringBuilder();

            for (int i = 0; i < windowAndDownstream.length(); i++) {
                String aminoAcid = windowAndDownstream.substring(i,i+1);
                int index = rand.nextInt(100);
                List<String> aminoAcidCodons = codon_usage.get(aminoAcid);
                String chosenCodon = aminoAcidCodons.get(index);
                currentSequence.append(chosenCodon);
            }

            sequences.add(currentSequence.toString());
        }
        return sequences;
    }

    private String[] stringToCodon(String cds) {
        int codonCount = cds.length() / 3;
        String[] codons = new String[codonCount];
        for (int i = 0; i < codonCount; i++) {
            codons[i] = cds.substring(i * 3, (i + 1) * 3);
        }
        return codons;
    }
}
