package picard.analysis.library;

import com.google.gson.Gson;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.reflect.TypeToken;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Tuple;
import org.apache.commons.math3.distribution.BetaDistribution;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Alpha;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Type;
import java.nio.file.Files;
import java.util.*;

/**
 * <p>Checks the library identities between a collection of read groups. For every pair of read groups, a same-
 * library / different-library decision is generated. A corresponding log-likelihood ratio (same-lib / diff-lib)
 * quantifies the confidence in the decision. False decisions are written to a MetricsFile.</p>
 * <p>
 * This implementation relies on same-library and different-library similarity value distributions pre-trained against
 * trusted data. Currently, the user must supply appropriate shape parameters for these distributions. We intend to
 * preload sensible default values (while retaining option for user-override).</p>
 * <p>
 * The tool currently relies on feature inputs generated by {@link FindReadGroupDistributions}, stored in Json format.
 * Ultimately, we will wish to merge these two programs. Caching of values will be retained, as it is useful in some
 * circumstances.</p>
 * <p>
 * This is an alpha implementation. Further development is expected. In particular, more tests are required.</p>
 *
 * @author Bradley Taylor
 */
@CommandLineProgramProperties(
        usage = " <p>Checks the library identities between a collection of read groups. For every pair of read groups, a same-" +
                " library / different-library decision is generated. A corresponding log-likelihood ratio (same-lib / diff-lib)" +
                " quantifies the confidence in the decision. False decisions are written to a MetricsFile.\n" +
                " </p>" +
                " <p>This implementation relies on same-library and different-library similarity value distributions pre-trained against" +
                " trusted data. Currently, the user must supply appropriate shape parameters for these distributions. We intend to" +
                " preload sensible default values (while retaining option for user-override).\n" +
                " </p>" +
                " <p>The tool currently relies on feature inputs generated by {@link FindReadGroupDistributions}, stored in Json format." +
                " Ultimately, we will wish to merge these two programs. Caching of values will be retained, as it is useful in some" +
                " circumstances." +
                " </p>" +
                " <p>This is an alpha implementation. Further development is expected. In particular, more tests are required.</p>",
        usageShort = "Checks the library identities between a collection of read groups.",
        programGroup =Alpha.class
)
public class CheckLibraryIdentity extends CommandLineProgram {

    @Option(shortName = "RGK",
            doc = "Read group key file, output by {@link FindReadGroupDistributions}. A json file used to translate TABLE_INPUT. Contains a mapping of read-group-index : { SAMReadGroupRecord }"
    )
    public File READ_GROUP_KEY;

    @Option(shortName = "TI",
            doc = "JSON file output by {@link FindReadGroupDistributions}. Contains observed patterns of insert distribution, encoded as a mapping of [read-group-indices] : duplicate-set-count")
    public File TABLE_INPUT;

    // TODO add an enum where you can specify data type (RNA, Exome, PCR-plus Genome, PCR-free Genome). Include a datastructure, keyed off of this enum, that presets the appropriate values for beta distributions. Make it mutex with all the beta dist + threshold stuff below.

    @Option(shortName = "AD",
            doc = "Pretrained alpha parameter for the different-library beta distribution.")
    public Double ALPHA_DIFF;

    @Option(shortName = "BD",
            doc = "Pretrained beta parameter for the different-library beta distribution.")
    public Double BETA_DIFF;

    @Option(shortName = "AS",
            doc = "Pretrained alpha parameter for the same-library beta distribution.")
    public Double AlPHA_SAME;

    @Option(shortName = "BS",
            doc = "Pretrained beta parameter for the same-library beta distribution.")
    public Double BETA_SAME;

    @Option(shortName = "LT",
            doc = "Pretrained log-likelihood ratio threshold.")
    public Double LLR_THRESHOLD = 0.0;

    @Option(shortName = "OM",
            doc = "Include in the detailed metrics output read group pairs whose inferred library identity matches existing annotations. By default, only mismatches are output.")
    public Boolean OUTPUT_MATCHES;

    @Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Name of the metrics output file.")
    public File OUTPUT;

    /** A key associating a full SAMReadGroupRecord with a numerical index */
    private Map<Integer, SAMReadGroupRecord> readGroups = new HashMap<>();

    /**
     * A map of read group petterns to duplicate set counts. Key is a tuple of read group indices, from readGroups.
     * (1,1,1,2) : 314 indicates 314 duplicate sets in the original data set had 4 reads, 3 from RG1 and 1 from RG2
     * */
    private Map<ArrayList<Integer>, Integer> featureData = new HashMap<>();

    private static final Gson gson = new Gson();


    // Stock main method
    public static void main(final String[] args) {
        new CheckLibraryIdentity().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        readInputs();
        if (readGroups.size() < 2) {
            throw new PicardException("Read group key must contain at least two read groups.");
        }

        final BetaDistribution diffLibraryDistribution = new BetaDistribution(ALPHA_DIFF, BETA_DIFF);
        final BetaDistribution sameLibraryDistribution = new BetaDistribution(AlPHA_SAME, BETA_SAME);

        final Map<Tuple<Integer, Integer>, Double> similarities = calculateReadGroupSimilarities();
        if (similarities.size() == 0) {
            throw new PicardException("Obtained empty pairwise similarities. Input must have at least one pair of read groups. Both must be represented in table data.");
        }

        final MetricsFile<ReadGroupSimilarityMetrics,?> metricsFile = getMetricsFile();
        updateReadGroupSimilarityMetrics(metricsFile, similarities, diffLibraryDistribution, sameLibraryDistribution);
        metricsFile.write(OUTPUT);

        return 0;
    }

    /**
     * Calculates a similarity value for every pair of read groups.
     *
     * A rough, inefficient implementation of the Exact Pairs Affinity custom similarity measure.
     *      Method requires testing and optimization.
     *
     * Similarity values are capped slightly below 1, so that they can be handled by beta distributions
     *
     * @return Similarity map
     */
    private Map<Tuple<Integer, Integer>, Double> calculateReadGroupSimilarities() {
        final Map<Tuple<Integer, Integer>, Double> similarities = new HashMap<>(); // todo - set initial size. It's predictable.

        // Collect the raw counts
        for (Map.Entry<ArrayList<Integer>, Integer> entry : featureData.entrySet()) {
            List<Integer> groups = entry.getKey();
            Double count = entry.getValue().doubleValue();

            if (groups.size() == 2) {
                Integer group1 = groups.get(0);
                Integer group2 = groups.get(1);
                Tuple<Integer, Integer> groupKey = new Tuple<>(group1, group2);
                similarities.merge(groupKey, count, Double::sum);
            }
        }

        // Normalize different-read-group counts according to same-read-group counts. NB- we never update the same-read-group counts
        for (int group1 = 0; group1 < readGroups.size()-1; group1++) {
            for (int group2 = (group1+1); group2 < readGroups.size(); group2++) {
                Tuple<Integer, Integer> abPattern = new Tuple<>(group1, group2);
                Double aa = similarities.get(new Tuple<>(group1, group1));
                Double ab = similarities.get(abPattern);
                Double bb = similarities.get(new Tuple<>(group2, group2));

                Double ExactPairsAffinity = ab / (2 * Math.sqrt(aa * bb));
                similarities.put(abPattern, ExactPairsAffinity);

                System.out.println(group1);
                System.out.println(aa);
                System.out.println(group2);
                System.out.println(bb);
            }
        }

        // Just throw out the same-group entries. We don't want them around to confuse things
        for (int group = 0; group < readGroups.size(); group++) {
            similarities.remove(new Tuple<>(group, group));
        }
        return similarities;
    }

    /**
     * This method does the work of evaluating observed similarities against the same/different library model,
     *      obtaining log-likelihood ratios, and generating same/different decisions.
     *
     * Results are written to a metrics file.
     *      NB: Currently outputting a detailed metric listing every false decision.
     *          Will develop a summary metric that aggregates false-decisions by read group
     *
     * @param metricsFile  MetricsFile for ReadGroupSimilarityMetrics
     * @param similarities A map associating every pair of read groups with a similarity value
     * @param diffLibraryDistribution Beta distribution pre-trained across known different-library pairs
     * @param sameLibraryDistribution Beta distribution pre-trained across known same-library pairs
     */
    private void updateReadGroupSimilarityMetrics(MetricsFile metricsFile, Map<Tuple<Integer, Integer>, Double> similarities, BetaDistribution diffLibraryDistribution, BetaDistribution sameLibraryDistribution) {
        for (Map.Entry<Tuple<Integer, Integer>, Double> entry : similarities.entrySet()) {
            Tuple<Integer, Integer> groups = entry.getKey();
            SAMReadGroupRecord group1 = readGroups.get(groups.a);
            SAMReadGroupRecord group2 = readGroups.get(groups.b);

            Double similarity = entry.getValue();

            Double logLikelihoodDiffLibrary = Math.log(diffLibraryDistribution.density(similarity)); // NB: natural log
            Double logLikelihoodSameLibrary = Math.log(sameLibraryDistribution.density(similarity));
            Double logLikelihoodRatio = logLikelihoodSameLibrary - logLikelihoodDiffLibrary;
            Boolean sameLibrary = logLikelihoodRatio > LLR_THRESHOLD;

            Boolean expectedSame = group1.getLibrary().equals(group2.getLibrary());
            Boolean observedLibraryIdentityMatchesExpected = sameLibrary.equals(expectedSame);

            if (OUTPUT_MATCHES || !observedLibraryIdentityMatchesExpected) {
                ReadGroupSimilarityMetrics metrics = new ReadGroupSimilarityMetrics();
                metrics.SIMILARITY = similarity;
                metrics.READ_GROUP_ID_1 = group1.getPlatformUnit(); // full-flowcell-barcode.lane.molecular-barcode-sequence.
                metrics.READ_GROUP_ID_2 = group2.getPlatformUnit();
                metrics.SAMPLE_1 = group1.getSample();
                metrics.SAMPLE_2 = group2.getSample();
                metrics.EXPECTED_LIBRARY_1 = group1.getLibrary();
                metrics.EXPECTED_LIBRARY_2 = group2.getSample();
                metrics.LLR_SAME_LIBRARY = logLikelihoodRatio;
                metrics.OBSERVED_SAME_LIBRARY = sameLibrary;
                metrics.OBSERVED_MATCHES_EXPECTED = observedLibraryIdentityMatchesExpected;
                metricsFile.addMetric(metrics);
            }
        }
    }

    // Read json inputs created by FindReadGroupDistributions.
    private void readInputs() {
        IOUtil.assertFileIsReadable(READ_GROUP_KEY);
        IOUtil.assertFileIsReadable(TABLE_INPUT);
        try {
            String tableJson = Files.readAllLines(TABLE_INPUT.toPath()).get(0);
            parseTableJson(tableJson);

            String keyJson = Files.readAllLines(READ_GROUP_KEY.toPath()).get(0);
            parseReadGroupKeyJson(keyJson);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private void parseReadGroupKeyJson(String keyJson) {
        JsonParser parser = new JsonParser();
        JsonObject job = parser.parse(keyJson).getAsJsonObject();
        for (Map.Entry<String, JsonElement> entry : job.entrySet()) {
            Integer id = gson.fromJson(entry.getKey(), Integer.class);

            String recordId = entry.getValue().getAsJsonObject().get("mReadGroupId").getAsString();
            SAMReadGroupRecord record = new SAMReadGroupRecord(recordId);

            // NB: Can probably get away with just LB and SM here, but I am including everything from our actual production readgroups for completeness. Can reduce when we're sure what isn't needed.
            JsonObject recordAttributes = entry.getValue().getAsJsonObject().get("mAttributes").getAsJsonObject();
            record.setLibrary(recordAttributes.get("LB").getAsString());
            record.setPlatform(recordAttributes.get("PL").getAsString());
            record.setPlatformUnit(recordAttributes.get("PU").getAsString());
//                record.setRunDate(recordAttributes.get("DT").getAsString()); //TODO- Needs conversion to Date. This field is in the read group record, but it's not really used by the rest of the program. Add if convenient.
            record.setSample(recordAttributes.get("SM").getAsString());
            record.setSequencingCenter(recordAttributes.get("CN").getAsString());

            readGroups.put(id, record);
        }
    }

    private void parseTableJson(String tableJson) {
        Type listType = new TypeToken<List<Integer>>() {}.getType();
        JsonParser parser = new JsonParser();
        JsonObject job = parser.parse(tableJson).getAsJsonObject();
        for (Map.Entry<String,JsonElement> entry : job.entrySet()) {
            ArrayList<Integer> groups = gson.fromJson(entry.getKey(), listType);
            Integer count = (entry.getValue().getAsInt());
            featureData.put(groups, count);
        }
    }
}
