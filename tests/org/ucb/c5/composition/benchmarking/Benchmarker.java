package org.ucb.c5.composition.benchmarking;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;
import org.ucb.c5.composition.CompositionToDNA;
import org.ucb.c5.composition.TranscriptDesigner;
import org.ucb.c5.composition.TranscriptDesigner2;
import org.ucb.c5.composition.checkers.ForbiddenSequenceChecker;
import org.ucb.c5.composition.model.Composition;
import org.ucb.c5.composition.model.Construct;
import org.ucb.c5.composition.model.Host;
import org.ucb.c5.composition.model.Transcript;
import org.ucb.c5.sequtils.HairpinCounter;
import org.ucb.c5.sequtils.Translate;
import org.ucb.c5.utils.FileUtils;

/**
 *
 * @author J. Christopher Anderson
 */
public class Benchmarker {

    private static Logger logger;
    private static FileHandler msgHandler;

    static {
        try {
            File logfil = new File("BenchmarkerLog.txt");
            FileUtils.writeFile("", logfil.getAbsolutePath());
            msgHandler = new FileHandler(logfil.getAbsolutePath(), true);
            logger = Logger.getLogger("C5");
            logger.setLevel(Level.ALL);
            logger.addHandler(msgHandler);
            msgHandler.setFormatter(new Formatter() {
                @Override
                public String format(LogRecord record) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(record.getMessage());
                    sb.append("\n");
                    return sb.toString();
                }
            });
        } catch (Exception ex) {
            System.err.println("Unable to initiate logger");
            System.exit(0);
        }
    }

    public static void main(String[] args) throws Exception {
        Translate translate = new Translate();
        translate.initiate();
        ForbiddenSequenceChecker fsc = new ForbiddenSequenceChecker();
        fsc.initiate();
        HairpinCounter haircount = new HairpinCounter();
        haircount.initiate();

        double score = 0;

        try {
            //Does TranscriptDesigner compile, instantiate, and initiate?
            TranscriptDesigner2 td2 = new TranscriptDesigner2();
            score++;
            logger.info("*  TranscriptDesigner2 instantiates   +1");
            td2.initiate();
            score++;
            logger.info("*  TranscriptDesigner2 initiates   +1");

            //Does TranscriptDesigner run a single simple-but-long sequence and return non-empty objects?
            String firstMyosinResult = null;
            try {
                Transcript myosin_Result = td2.run("MPKPVANQEDEDPTPYLFVSLEQRRIDQSKPYDSKKSCWIPDEKEGYLLGEIKATKGDIVSVGLQGGEVRDIKSEKVEKVNPPKFEKIEDMADMTVLNTPCVLHNLRQRYYAKLIYTYSGLFCVAINPYKRYPVYTNRCAKMYRGKRRNEVPPHIFAISDGAYVDMLTNHVNQSMLITGESGAGKTENTKKVIAYFATVGASKKTDEAAKSKGSLEDQVVQTNPVLEAFGNAKTVRNDNSSRFGKFIRIHFGPTGKLAGADIETYLLEKARVISQQSLERSYHIFYQIMSGSVPGVKEMVFLGQHIGDYPGICQGKTRIPGVNDGEEFELTDQAFDILGFTKQEKEDVYRITAAVMHMGGMKFKQRGREEQAEQDGEEEGGRVSKLFGCDTAELYKNLLKPRIKVGNEFVTQGRNVQQVTNSIGALCKGVFDRLFKWLVKKCNETLDTQQKRQHFIGVLDIAGFEIFEYNGFEQLCINFTNEKLQQFFNHHMFVLEQEEYKREGIDWAFIDFGMDLLACIDLIEKPMGILSILEEESMFPKATDQTFSEKLTNTHLGKSAPFQKPKPPKPGQQAAHFAIAHYAGCVSYNITGWLEKNKDPLNDTVVDQFKKSQNKLLIEIFADHAGQSGGGEQAKGGRGKKGGGFATVSSAYKEQLNSLMTTLRSTQPHFVRCIIPNEMKQPGVVDAHLVMHQLTCNGVLEGIRICRKGFPNRMMYPDFKMRYKIMCPKLLQGVEKDKKATEIIIKFIDLPEDQYRLGNTKVFFRAGVLGQMEEFRDERLGKIMSWMQAWARGYLSRKGFKKLQEQRVALKVVQRNLRKYLQLRTWPWYKLWQKVKPLLNVSRIEDEIARLEEKAKKAEELHAAEVKVRKELEALNAKLLAEKTALLDSLSGEKGALQDYQERNAKLTAQKNDLENQLRDIQERLTQEEDARNQLFQQKKKADQEISGLKKDIEDLELNVQKAEQDKATKDHQIRNLNDEIAHQDELINKLNKEKKMQGETNQKTGEELQAAEDKINHLNKVKAKLEQTLDELEDSLEREKKVRGDVEKSKRKVEGDLKLTQEAVADLERNKKELEQTIQRKDKELSSITAKLEDEQVVVLKHQRQIKELQARIEELEEEVEAERQARAKAEKQRADLARELEELGERLEEAGGATSAQIELNKKREAELSKLRRDLEEANIQHESTLANLRKKHNDAVAEMAEQVDQLNKLKAKAEHDRQTCHNELNQTRTACDQLGRDKAAQEKIAKQLQHTLNEVQSKLDETNRTLNDFDASKKKLSIENSDLLRQLEEAESQVSQLSKIKISLTTQLEDTKRLADEESRERATLLGKFRNLEHDLDNLREQVEEEAEGKADLQRQLSKANAEAQVWRSKYESDGVARSEELEEAKRKLQARLAEAEETIESLNQKCIGLEKTKQRLSTEVEDLQLEVDRANAIANAAEKKQKAFDKIIGEWKLKVDDLAAELDASQKECRNYSTELFRLKGAYEEGQEQLEAVRRENKNLADEVKDLLDQIGEGGRNIHEIEKARKRLEAEKDELQAALEEAEAALEQEENKVLRAQLELSQVRQEIDRRIQEKEEEFENTRKNHQRALDSMQASLEAEAKGKAEALRMKKKLEADINELEIALDHANKANAEAQKNIKRYQQQLKDIQTALEEEQRARDDAREQLGISERRANALQNELEESRTLLEQADRGRRQAEQELADAHEQLNEVSAQNASISAAKRKLESELQTLHSDLDELLNEAKNSEEKAKKAMVDAARLADELRAEQDHAQTQEKLRKALEQQIKELQVRLDEAEANALKGGKKAIQKLEQRVRELENELDGEQRRHADAQKNLRKSERRVKELSFQSEEDRKNHERMQDLVDKLQQKIKTYKRQIEEAEEIAALNLAKFRKAQQELEEAEERADLAEQAISKFRAKGRAGSVGRGASPAPRATSVRPQFDGLAFPPRFDLAPENEF", new HashSet<>());
                int length = myosin_Result.getCodons().length;  //Will cause Exception if empty
                score++;
                logger.info("*  TranscriptDesigner2 runs   +1");
                firstMyosinResult = myosin_Result.toSeq();
            } catch (Exception err) {
                logger.warning("!!! Unable to run myosin");
            }

            //Check for determinacy
            try {
                Transcript myosin_Result = td2.run("MPKPVANQEDEDPTPYLFVSLEQRRIDQSKPYDSKKSCWIPDEKEGYLLGEIKATKGDIVSVGLQGGEVRDIKSEKVEKVNPPKFEKIEDMADMTVLNTPCVLHNLRQRYYAKLIYTYSGLFCVAINPYKRYPVYTNRCAKMYRGKRRNEVPPHIFAISDGAYVDMLTNHVNQSMLITGESGAGKTENTKKVIAYFATVGASKKTDEAAKSKGSLEDQVVQTNPVLEAFGNAKTVRNDNSSRFGKFIRIHFGPTGKLAGADIETYLLEKARVISQQSLERSYHIFYQIMSGSVPGVKEMVFLGQHIGDYPGICQGKTRIPGVNDGEEFELTDQAFDILGFTKQEKEDVYRITAAVMHMGGMKFKQRGREEQAEQDGEEEGGRVSKLFGCDTAELYKNLLKPRIKVGNEFVTQGRNVQQVTNSIGALCKGVFDRLFKWLVKKCNETLDTQQKRQHFIGVLDIAGFEIFEYNGFEQLCINFTNEKLQQFFNHHMFVLEQEEYKREGIDWAFIDFGMDLLACIDLIEKPMGILSILEEESMFPKATDQTFSEKLTNTHLGKSAPFQKPKPPKPGQQAAHFAIAHYAGCVSYNITGWLEKNKDPLNDTVVDQFKKSQNKLLIEIFADHAGQSGGGEQAKGGRGKKGGGFATVSSAYKEQLNSLMTTLRSTQPHFVRCIIPNEMKQPGVVDAHLVMHQLTCNGVLEGIRICRKGFPNRMMYPDFKMRYKIMCPKLLQGVEKDKKATEIIIKFIDLPEDQYRLGNTKVFFRAGVLGQMEEFRDERLGKIMSWMQAWARGYLSRKGFKKLQEQRVALKVVQRNLRKYLQLRTWPWYKLWQKVKPLLNVSRIEDEIARLEEKAKKAEELHAAEVKVRKELEALNAKLLAEKTALLDSLSGEKGALQDYQERNAKLTAQKNDLENQLRDIQERLTQEEDARNQLFQQKKKADQEISGLKKDIEDLELNVQKAEQDKATKDHQIRNLNDEIAHQDELINKLNKEKKMQGETNQKTGEELQAAEDKINHLNKVKAKLEQTLDELEDSLEREKKVRGDVEKSKRKVEGDLKLTQEAVADLERNKKELEQTIQRKDKELSSITAKLEDEQVVVLKHQRQIKELQARIEELEEEVEAERQARAKAEKQRADLARELEELGERLEEAGGATSAQIELNKKREAELSKLRRDLEEANIQHESTLANLRKKHNDAVAEMAEQVDQLNKLKAKAEHDRQTCHNELNQTRTACDQLGRDKAAQEKIAKQLQHTLNEVQSKLDETNRTLNDFDASKKKLSIENSDLLRQLEEAESQVSQLSKIKISLTTQLEDTKRLADEESRERATLLGKFRNLEHDLDNLREQVEEEAEGKADLQRQLSKANAEAQVWRSKYESDGVARSEELEEAKRKLQARLAEAEETIESLNQKCIGLEKTKQRLSTEVEDLQLEVDRANAIANAAEKKQKAFDKIIGEWKLKVDDLAAELDASQKECRNYSTELFRLKGAYEEGQEQLEAVRRENKNLADEVKDLLDQIGEGGRNIHEIEKARKRLEAEKDELQAALEEAEAALEQEENKVLRAQLELSQVRQEIDRRIQEKEEEFENTRKNHQRALDSMQASLEAEAKGKAEALRMKKKLEADINELEIALDHANKANAEAQKNIKRYQQQLKDIQTALEEEQRARDDAREQLGISERRANALQNELEESRTLLEQADRGRRQAEQELADAHEQLNEVSAQNASISAAKRKLESELQTLHSDLDELLNEAKNSEEKAKKAMVDAARLADELRAEQDHAQTQEKLRKALEQQIKELQVRLDEAEANALKGGKKAIQKLEQRVRELENELDGEQRRHADAQKNLRKSERRVKELSFQSEEDRKNHERMQDLVDKLQQKIKTYKRQIEEAEEIAALNLAKFRKAQQELEEAEERADLAEQAISKFRAKGRAGSVGRGASPAPRATSVRPQFDGLAFPPRFDLAPENEF", new HashSet<>());
                String secondMyosinResult = myosin_Result.toSeq();
                if (!firstMyosinResult.equals(secondMyosinResult)) {
                    score = score - 5;
                    logger.warning("!!! Failed determinacy test -5");
                }
            } catch (Exception err) {
                score = score - 5;
                logger.warning("!!! Unable to run determinacy test -5");
            }

            //Does TranscriptDesigner run 20 ~100 aa MG1655-derived proteins?
            String[] prot20 = new String[20];   //Protein sequences   (MSKE...)
            String[] cds20 = new String[20];    //CDS sequences (ATG...)
            Transcript[] tran20 = new Transcript[20];   //Transcripts (RBS + CDS)
            long[] timings = new long[20];      //Run time in millis

            prot20[0] = "MSKEHTTEHLRAELKSLSDTLEEVLSSSGEKSKEELSKIRSKAEQALKQSRYRLGETGDAIAKQTRVAAARADEYVRENPWTGVGIGAAIGVVLGVLLSRR";
            prot20[1] = "MATLLQLHFAFNGPFGDAMAEQLKPLAESINQEPGFLWKVWTESEKNHEAGGIYLFTDEKSALAYLEKHTARLKNLGVEEVVAKVFDVNEPLSQINQAKLA";
            prot20[2] = "MSNQFGDTRIDDDLTLLSETLEEVLRSSGDPADQKYVELKARAEKALDDVKKRVSQASDSYYYRAKQAVYRADDYVHEKPWQGIGVGAAVGLVLGLLLARR";
            prot20[3] = "MAKGQSLQDPFLNALRRERVPVSIYLVNGIKLQGQIESFDQFVILLKNTVSQMVYKHAISTVVPSRPVSHHSNNAGGGTSSNYHHGSSAQNTSAQQDSEETE";
            prot20[4] = "MQPNDITFFQRFQDDILAGRKTITIRDESESHFKTGDVLRVGRFEDDGYFCTIEVTATSTVTLDTLTEKHAEQENMTLTELKKVIADIYPGQTQFYVIEFKCL";
            prot20[5] = "MLTVIAEIRTRPGQHHRQAVLDQFAKIVPTVLKEEGCHGYAPMVDCAAGVSFQSMAPDSIVMIEQWESIAHLEAHLQTPHMKAYSEAVKGDVLEMNIRILQPGI";
            prot20[6] = "MMIRERIEEKLRAAFQPVFLEVVDESYRHNVPAGSESHFKVVLVSDRFTGERFLNRHRMIYSTLAEELSTTVHALALHTYTIKEWEGLQDTVFASPPCRGAGSIA";
            prot20[7] = "MAEWSGEYISPYAEHGKKSEQVKKITVSIPLKVLKILTDERTRRQVNNLRHATNSELLCEAFLHAFTGQPLPDDADLRKERSDEIPEAAKEIMREMGINPETWEY";
            prot20[8] = "MNDSEFHRLADQLWLTIEERLDDWDGDSDIDCEINGGVLTITFENGSKIIINRQEPLHQVWLATKQGGYHFDLKGDEWICDRSGETFWDLLEQAATQQAGETVSFR";
            prot20[9] = "MAKNRSRRLRKKMHIDEFQELGFSVAWRFPEGTSEEQIDKTVDDFINEVIEPNKLAFDGSGYLAWEGLICMQEIGKCTEEHQAIVRKWLEERKLDEVRTSELFDVWWD";
            prot20[10] = "MFGKGGLGNLMKQAQQMQEKMQKMQEEIAQLEVTGESGAGLVKVTINGAHNCRRVEIDPSLLEDDKEMLEDLVAAAFNDAARRIEETQKEKMASVSSGMQLPPGFKMPF";
            prot20[11] = "METTKPSFQDVLEFVRLFRRKNKLQREIQDVEKKIRDNQKRVLLLDNLSDYIKPGMSVEAIQGIIASMKGDYEDRVDDYIIKNAELSKERRDISKKLKAMGEMKNGEAK";
            prot20[12] = "MSFFISDAVAATGAPAQGSPMSLILMLVVFGLIFYFMILRPQQKRTKEHKKLMDSIAKGDEVLTNGGLVGRVTKVAENGYIAIALNDTTEVVIKRDFVAAVLPKGTMKAL";
            prot20[13] = "MAESFTTTNRYFDNKHYPRGFSRHGDFTIKEAQLLERHGYAFNELDLGKREPVTEEEKLFVAVCRGEREPVTEAERVWSKYMTRIKRPKRFHTLSGGKPQVEGAEDYTDSDD";
            prot20[14] = "MKKIDAIIKPFKLDDVREALAEVGITGMTVTEVKGFGRQKGHTELYRGAEYMVDFLPKVKIEIVVPDDIVDTCVDTIIRTAQTGKIGDGKIFVFDVARVIRIRTGEEDDAAI";
            prot20[15] = "MMKKSILAFLLLLTSSAAALAAPQVITVSRFEVGKDKWAFNREEVMLTCRPGNALYVINPSTLVQYPLNDIAQKEVASGKTNAQPISVIQIDDPNNPGEKMSLAPFIERAEKLC";
            prot20[16] = "MRIFVYGSLRHKQGNSHWMTNAQLLGDFSIDNYQLYSLGHYPGAVPGNGTVHGEVYRIDNATLAELDALRTRGGEYARQLIQTPYGSAWMYVYQRPVDGLKLIESGDWLDRDK";
            prot20[17] = "MTIVRIDAEARWSDVVIHNNTLYYTGVPENLDADAFEQTANTLAQIDAVLEKQGSNKSSILDATIFLADKNDFAAMNKAWDAWVVAGHAPVRCTVQAGLMNPKYKVEIKIVAAV";
            prot20[18] = "MQKIVIVANGAPYGSSSSESLFNSLRLAIALREQESNLDLRLFLMSDAVTAGLRGQKPGEGYNIQQMLEILTAQNVPVKLCKTCTDGRGISTLPLIDGVEIGTLVELAQWTLSADKVLTF";
            prot20[19] = "MVTLYGIKNCDTIKKKKARRWLEANNIDYRFHDYRVDGLDSELLNDFINELGWEALLNTRGTTWRKLDETTRNKITDAASAAALMTEMPAIIKRPLLCVPGKPMLLGFSDSSYQQFFHEV";

            boolean runs20 = true;

            for (int i = 0; i < 20; i++) {
                String cds = prot20[i];
                long startTime = System.currentTimeMillis();
                try {
                    Transcript t = td2.run(cds, new HashSet<>());
                    StringBuilder sb = new StringBuilder();
                    for (String codon : t.getCodons()) {
                        sb.append(codon);
                    }
                    cds20[i] = sb.toString();
                    tran20[i] = t;
                } catch (Exception err) {
                    logger.warning("!!! Exception when running cds #:  " + i + "\n" + prot20[i]);
                    runs20 = false;
                }
                long stopTime = System.currentTimeMillis();
                timings[i] = stopTime - startTime;
            }

            if (runs20) {
                score++;
                logger.info("*  Ran without Exception on 20 E. coli proteins  +1");
            }

            //Examining the results of the 20 protein sequences
            boolean nonenull = true;
            boolean alltranslate = true;
            boolean allstartandstop = true;

            for (int i = 0; i < 20; i++) {
                if (cds20[i] == null) {
                    nonenull = false;
                    logger.warning("!!! Protein returned null:  " + i + "\n" + prot20[i]);
                    continue;
                }

                if (!translate.run(cds20[i]).equals(prot20[i])) {
                    alltranslate = false;
                    logger.warning("!!! Protein doesn't match translation: " + i + "\n" + prot20[i]);
                }

                if (!cds20[i].startsWith("ATG")) {
                    allstartandstop = false;
                    logger.warning("!!! CDS doesn't begin with start codon for: " + i + "\n" + prot20[i]);
                }

                if (tran20[i].toSeq().endsWith("TAA") || tran20[i].toSeq().endsWith("TGA") || tran20[i].toSeq().endsWith("TAG")) {
                } else {
                    allstartandstop = false;
                    logger.warning("!!! CDS doesn't end with stop codon for:  " + i + "\n" + prot20[i]);
                }
            }

            //Are all 20 returned Transcripts non-null?
            if (nonenull) {
                score += 2;
                logger.info("*  All 20 E. coli proteins return non-null  +2");
            }

            //Do the 20 returned sequences translate as original protein?
            if (alltranslate) {
                score += 5;
                logger.info("*  All 20 designed CDS translate correctly  +5");
            }

            //Do the 20 returned sequences begin with a start and end with a stop?
            if (allstartandstop) {
                score += 5;
                logger.info("*  All 20 designed CDS have start and stop  +5");
            }

            //Do the 20 returned sequences pass ForbiddenSequenceChecker?
            {
                int numpass = 0;
                int numfail = 0;
                for (int i = 0; i < 20; i++) {
                    if (cds20[i] == null) {
                        continue;
                    }
                    boolean result = fsc.run(cds20[i]);
                    if (result) {
                        numpass++;
                    } else {
                        numfail++;
                        logger.warning("!!! CDS has forbidden sequences for: " + i + "\n" + prot20[i]);
                    }
                }

                if (numfail == 0 && numpass > 10) {
                    score += 5;
                    logger.info("*  Algorithm generates no forbidden sequences   +5");
                }
                if (numpass == 20) {
                    score += 5;
                    logger.info("*  Algorithm generates no forbidden sequences for all 20   +5");
                }
            }

            //Do the 20 returned sequences have fewer hairpins than original algorithm?
            {
                //Use the old algorithm to compare
                TranscriptDesigner td = new TranscriptDesigner();
                td.initiate();

                int numpass = 0;
                for (int i = 0; i < 20; i++) {
                    if (cds20[i] == null) {
                        continue;
                    }

                    //Calculate hairpin count for the new algorithm's result
                    double hcNew = haircount.run(tran20[i].toSeq());

                    //Calculate hairpin count for the old algorithm
                    Transcript oldT = td.run(prot20[i], new HashSet<>());
                    double hcOld = haircount.run(oldT.toSeq());

                    //Decrement the old score to require improvement
                    hcOld = hcOld - 1;

                    //Compare to the counts
                    if (hcNew < hcOld) {
                        numpass++;
                    } else {
                        logger.warning("!!! CDS has secondary structure for:  " + i + "\n" + prot20[i]);
                    }
                }

                //increments score by <=10
                score += .5 * numpass;
                logger.info("*  Algorithm generates fewer hairpins for " + numpass + " sequences   +" + 0.5 * numpass);
            }

            //Does it use different codons?
            Map<String, Integer> codonsToCount = new HashMap<>();
            for (int i = 0; i < 20; i++) {
                String[] codons = tran20[i].getCodons();
                for (int x = 0; x < codons.length; x++) {
                    String codon = codons[x];
                    Integer count = codonsToCount.get(codon);
                    if (count == null) {
                        count = 0;
                    }
                    count++;
                    codonsToCount.put(codon, count);
                }
            }
            logger.info("*  This many codons were used: " + codonsToCount.size());
            if (codonsToCount.size() < 25) {
                score = score - 5;
                logger.warning("!!! Algorithm is not using alternate codons -5");
            }

            //Is it better than random at avoiding rare AGG codons?
            int argCount = 0;
            Integer aggCount = codonsToCount.get("AGG");
            if (aggCount == null) {
                aggCount = 0;
            }
            argCount += aggCount;

            Integer agaCount = codonsToCount.get("AGA");
            if (agaCount != null) {
                argCount += agaCount;
            }

            Integer cgaCount = codonsToCount.get("CGA");
            if (cgaCount != null) {
                argCount += cgaCount;
            }

            Integer cgtCount = codonsToCount.get("CGT");
            if (cgtCount != null) {
                argCount += cgtCount;
            }

            Integer cgcCount = codonsToCount.get("CGC");
            if (cgcCount != null) {
                argCount += cgcCount;
            }

            Integer cggCount = codonsToCount.get("CGG");
            if (cggCount != null) {
                argCount += cggCount;
            }

            if (aggCount >= argCount / 6) {
                score -= 5;
                logger.warning("!!! Algorithm is not reducing AGG codons -5");
            }

            //Is the maximum run time acceptable?
            float maxtime = 10000;  //10 sec in millis
            boolean timeok = true;
            for (int i = 0; i < 20; i++) {
                if (timings[i] > maxtime) {
                    timeok = false;
                    break;
                }
            }
            if (timeok) {
                score += 2;
                logger.info("*  Algorithm has an acceptable run time   +2");
            }

            //Does CompositionToDNA run Compositions with 20 proteins?
            CompositionToDNA c2d = new CompositionToDNA();
            c2d.initiate();
            List<String> twentyProts = new ArrayList<>();
            for (int i = 0; i < 20; i++) {
                twentyProts.add(prot20[i]);
            }
            String promoter = "ttatgacaacttgacggctacatcattcactttttcttcacaa";
            String terminator = "TGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCC";
            String protein = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVE";
            Composition comp = new Composition(Host.Ecoli, promoter, twentyProts, terminator);
            String design = null;
            try {
                Construct dna = c2d.run(comp);
                design = dna.toSeq();
                logger.info("*  Algorithm handled a 20 protein composition");
            } catch (Exception err) {
                score -= 2;
                logger.warning("!!! Algorithm could not handle a 20 protein composition  -2");
            }

            //Does the returned Construct pass checkers including junction tests?
            if (design != null) {
                boolean result = fsc.run(design);
                if (result) {
                    score += 2;
                    logger.info("*  20 protein construct has no forbidden sequences   +2");
                } else {
                    logger.warning("!!! 20 protein construct has forbidden sequences");
                }
            }

            logger.info("*  Benchmarker completed.  Score: " + score);
        } catch (Exception err) {
            logger.severe("!!!  Exception encountered.  Score: " + score);
            System.exit(0);
        }
    }
}
