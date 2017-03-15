package picard.filter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;
import picard.PicardException;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

/**
 * Filter fragments or pairs with low mapping quality.
 *
 * If a record is paired, it is considered a match if either
 *   its quality OR the quality of its mate is < the given threshold.
 *
 * The mate's mapping quality is derived from the 'MQ' SAM tag. This tag is required for paired records.
 */
public class PairsMappingQualityFilter implements SamRecordFilter {

    private int minimumMappingQuality = Integer.MIN_VALUE;

    public PairsMappingQualityFilter(final int minimumMappingQuality) {
        this.minimumMappingQuality = minimumMappingQuality;
    }

    @Override
    public boolean filterOut(final SAMRecord record) {
        if (record.getReadPairedFlag()) {
            if (record.getMateUnmappedFlag()) { // NB: Unmapped mates have mapping quality 0
                return true;
            }
            Object mqTag = record.getAttribute("MQ");
            if (mqTag == null) {
                throw new PicardException("Found a paired record with no 'MQ' SAM tag: " + record.getReadName());
            } else {
                int mateMappingQuality = ((Number) record.getAttribute("MQ")).intValue();
                return (record.getMappingQuality() < this.minimumMappingQuality || mateMappingQuality < minimumMappingQuality);
            }
        } else {
            return record.getMappingQuality() < this.minimumMappingQuality;
        }
    }

    /**
     * Method not implemented.
     *
     * This filter is not intended to operate on arbitrary pairs of records.
     * Pairwise filtering is handled by the single-SAMRecord method via the 'MQ' SAM tag.
     */
    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
        throw new NotImplementedException();
    }
}
