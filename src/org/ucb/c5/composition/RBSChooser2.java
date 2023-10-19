package org.ucb.c5.composition;

import java.util.*;

import javassist.Translator;
import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.sequtils.CalcEditDistance;
import org.ucb.c5.sequtils.HairpinCounter;
import org.ucb.c5.sequtils.Translate;
import org.ucb.c5.utils.FileUtils;

/**
 * Second generation RBSChooser algorithm
 *
 * Employs a list of genes and their associated ribosome binding sites for
 * highly-expressed proteins in E. coli.
 *
 * @author J. Christopher Anderson
 * @author Aditya Vunnum
 */
public class RBSChooser2 {
    private List<RBSOption> rbss;
    private Translate translator;
    private HairpinCounter hairpinCounter;

    /**
     * Takes into two data tables, constructs a list of all potential rbs options
     * Each rbs option tells us the name of the gene, rbs sequence, cds sequence, and first 6 aas of cds
     */
    public void initiate() throws Exception {
        // Function used to get the amino acid residues from the cds sequence
        translator = new Translate();
        translator.initiate();

        // Instantiate and initiate HairpinCounter
        hairpinCounter = new HairpinCounter();
        hairpinCounter.initiate();

        //Read the data file
        String coli_data = FileUtils.readResourceFile("composition/data/coli_genes.txt");
        String rbs_options = FileUtils.readResourceFile("composition/data/rbs_options.txt");

        //Split it into lines (fancy Regex is to support Mac and Windows and Excel variations) - from lecture
        String[] lines = coli_data.split("\\r|\\r?\\n");

        // initializing map that maps gene name to the cds region
        // We will use this map to link coli_genes table to the rbs_options table
        Map<String, String> coli_genes = new HashMap<>();

        //Scan through each line of coli_genes.txt
        for (int i = 0; i < lines.length; i++) {
            String line = lines[i];
            if(line.isEmpty()) {continue;}
            //Split the line between tabs
            String[] tabs = line.split("\t");

            // tabs[1] contains the name of the gene
            // tabs[6] contains the cds of the gene
            coli_genes.put(tabs[1], tabs[6]);
        }

        //Splitting rbs_options into lines
        String[] rbs_lines = rbs_options.split("\\r|\\r?\\n");

        // initializing the list of all the rbs options from the rbs_options table
        rbss = new ArrayList<>();

        //Scan through each line of rbs_options.txt
        for (int i = 0; i < rbs_lines.length; i++) {
            String rbs_line = rbs_lines[i];
            if (rbs_line.isEmpty()) {continue;}
            //Split the rbs_line between white spaces
            String[] rbs_tabs = rbs_line.split("\\s+");;

            String geneName = rbs_tabs[0];
            String rbs = rbs_tabs[1];
            // extracting the gene cds from the map created from coli_genes table
            String cds = coli_genes.get(geneName);

            // Translate the DNA sequence to amino acid sequence.
            String first6aas = translator.run(cds.substring(0, 18));

            // I was unsure how to get a description, so just used name again to satisfy RBSOption constructor below
            String description = geneName;

            RBSOption singleRbsOpt = new RBSOption(geneName, description, rbs, cds, first6aas);
            rbss.add(singleRbsOpt);

        }
    }

    /**
     * Provided an ORF of sequence 'cds', this computes the best ribosome
     * binding site to use from a list of options.
     *
     * It also permits a list of options to exclude.
     *
     * @param cds The DNA sequence, ie ATGCATGAT...
     * @param ignores The list of RBS's to exclude
     * @return the best RBSOption based on the hairpin distance of rbs / beginning of cds and similarity between first 6 AAs
     * @throws Exception
     */
    public RBSOption run(String cds, Set<RBSOption> ignores) throws Exception {
        // ALL FUNCTIONS USED ARE EXPLAINED BELOW

        // Given a list of rbs options to ignore, returns only valid rbs options
        List<RBSOption> valid_rbss = filterRBSOptions(ignores);

        // Calculates hairpin scores and returns the 10 rbs options with the lowest hairpin scores (based on rbs + first 30 residues of cds)
        Set<RBSOption> lowestHairpinScores = getTop10LowestHairpinScores(valid_rbss, cds);

        // Calculates edit distance between first 6 AAs of each rbs and given cds, then chooses closest option
        RBSOption closestRBSOption = getClosestRBSOption(lowestHairpinScores, cds);

        return closestRBSOption;
    }

    /**
     * Provided list of ribosome binding sites to ignore, filters out these sites from
     * our current list of RBS from the rbs_options.txt file
     *
     * @param ignores The list of RBS's to exclude
     * @return filtered list of ribosome binding site options
     * @throws Exception
     */
    public List<RBSOption> filterRBSOptions(Set<RBSOption> ignores) throws Exception {
        List<RBSOption> filteredOptions = new ArrayList<>();
        Iterator<RBSOption> var3 = this.rbss.iterator();

        // iterates through our list of ribosome binding sites and removes sites that are also in ignores
        while (var3.hasNext()) {
            RBSOption rbsopt = var3.next();
            if (!ignores.contains(rbsopt)) {
                filteredOptions.add(rbsopt);
            }
        }

        if (filteredOptions.isEmpty()) {
            throw new Exception("All RBSOption objects are in the ignores set!");
        }

        return filteredOptions;
    }

    /**
     * Provided list of valid ribosome binding sites, counts the hairpins and determines
     * the 10 ribosome binding sites with the lowest scores
     *
     * @param rbsOptions the filtered list of ribosome binding sites
     * @return set containing the ribosome binding sites with the lowest hairpin scores
     * @throws Exception
     */
    public Set<RBSOption> getTop10LowestHairpinScores(List<RBSOption> rbsOptions, String cds) throws Exception {
        // Use a list to store both the RBSOption and its score (utility class ScoredRBSOption below this method)
        List<ScoredRBSOption> scoredOptions = new ArrayList<>();

        // for each rbs option, calculate hairpin score and add to list
        for (RBSOption rbsOption : rbsOptions) {
            String hairpin_region = rbsOption.getRbs() + cds.substring(0,30);
            double score = hairpinCounter.run(hairpin_region);
            scoredOptions.add(new ScoredRBSOption(rbsOption, score));
        }

        // Sort the list based on the scores
        Collections.sort(scoredOptions);

        // Return a set of the top 10 RBSOption with the lowest scores
        Set<RBSOption> top10 = new HashSet<>();
        for (int i = 0; i < 10 && i < scoredOptions.size(); i++) {
            top10.add(scoredOptions.get(i).getRbsOption());
        }

        return top10;
    }

    /**
     * Utility class to associate an RBSOption with its score
     *
     * @implements Comparable
     */
    private static class ScoredRBSOption implements Comparable<ScoredRBSOption> {
        private final RBSOption rbsOption;
        private final double score;

        // create a new class linking hairpin score to each rbs option
        public ScoredRBSOption(RBSOption rbsOption, double score) {
            this.rbsOption = rbsOption;
            this.score = score;
        }

        public RBSOption getRbsOption() {
            return rbsOption;
        }

        // need to override java's normal compareTo method to compare the score values
        @Override
        public int compareTo(ScoredRBSOption o) {
            return Double.compare(this.score, o.score);
        }
    }

    /**
     * Given a set of ribosome binding sites, calculates the edit distance between the first 6 AAs of the rbs's gene's ORF
     * and the first 6 AAs of the given open reading frame.
     *
     * @param topRBSOptions the filtered list of ribosome binding sites with low hairpin counter scores
     * @return RBSOption object containing the best ribosome binding site for the given sequence
     * @throws Exception
     */
    public RBSOption getClosestRBSOption(Set<RBSOption> topRBSOptions, String cds) throws Exception {
        String cdsFirst6AAs = translator.run(cds.substring(0,18));  // Extract first 6 AAs of CDS

        // Initiate edit distance calculator
        CalcEditDistance calculator = new CalcEditDistance();
        calculator.initiate();

        // initialize a variable for closest rbs to update iteratively
        RBSOption closestRBS = null;
        // Similarly, initiate a distance variable starting at a high value to compare each rbs sequence with. Updates iteratively
        int closestDistance = Integer.MAX_VALUE;  // Set to maximum possible integer value

        for (RBSOption rbsOption : topRBSOptions) {
            String rbsFirst6AAs = rbsOption.getFirst6aas();
            int distance = calculator.run(cdsFirst6AAs, rbsFirst6AAs);

            // if the current ribosome binding site is more similar than the one previously saved, update
            if (distance < closestDistance) {
                closestDistance = distance;
                closestRBS = rbsOption;
            }
        }

        // This will return the RBSOption with the closest 6 AAs to the CDS!
        return closestRBS;
    }


    public static void main(String[] args) throws Exception {
        //Create an example
        String cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATTATCACCATCGTTATATGACACAGCTTGGGTGGCTATGGTGCCGGAGCGTAGTTCTTCTCAA";

        //Initiate the chooser
        RBSChooser2 chooser = new RBSChooser2();
        chooser.initiate();

        //Make the first choice with an empty Set of ignores
        Set<RBSOption> ignores = new HashSet<>();
        RBSOption selected1 = chooser.run(cds, ignores);

        //Add the first selection to the list of things to ignore
        ignores.add(selected1);

        //Choose again with an ignore added
        RBSOption selected2 = chooser.run(cds, ignores);

        //Print out the two options, which should be different
        System.out.println("CDS starts with:");
        System.out.println(cds.substring(0, 18));
        System.out.println();
        System.out.println("Selected1:\n");
        System.out.println(selected1.toString());
        System.out.println();
        System.out.println("Selected2:\n");
        System.out.println(selected2.toString());
    }
}

