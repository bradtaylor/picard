package picard.analysis.library;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Created by bradt on 11/20/16. Detailed metrics information when inferring whether read groups belong to the same
 * library. See {@link CheckLibraryIdentity}.
 */
public class ReadGroupSimilarityMetrics extends MetricBase {
    public String READ_GROUP_ID_1;
    public String READ_GROUP_ID_2;
    public String SAMPLE_1;
    public String SAMPLE_2;
    public String EXPECTED_LIBRARY_1;
    public String EXPECTED_LIBRARY_2;
    public Double SIMILARITY;
    public Double LLR_SAME_LIBRARY;
    public Boolean OBSERVED_SAME_LIBRARY;
    public Boolean OBSERVED_MATCHES_EXPECTED;
}
