package picard.analysis.library;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by bradt on 1/25/17.
 */
public class FindReadGroupDistributionsTest {

    protected final static int DEFAULT_BASE_QUALITY = 10;

    @Test
    public void testReadGroupIdManager() {
        // Create a list of read groups with IDs "0" to "5"
        List<SAMReadGroupRecord> readGroups = new ArrayList<>();
        for (int i = 0; i < 6; i++) {
            SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(Integer.toString(i));
            readGroups.add(readGroupRecord);
        }
        // Pack read groups into a SAM header
        SAMFileHeader header = new SAMFileHeader();
        header.setReadGroups(readGroups);

        // Construct a ReadGroupIdManager, which will populate its internal maps with read groups from the header
        FindReadGroupDistributions.ReadGroupIdManager readGroupIdManager = new FindReadGroupDistributions.ReadGroupIdManager(header);

        Assert.assertFalse(readGroupIdManager.getReadGroupRecordMap().isEmpty());
        Assert.assertFalse(readGroupIdManager.getReadGroupRecordMap().keySet().contains(-1));

        // Check that read group IDs returned by the Manager match what is expected
        for (short i = 0; i < 6; i++) {
            Assert.assertEquals(readGroupIdManager.getGroupFromId(i).getId(), Short.toString(i));
            Assert.assertEquals(readGroupIdManager.getIdOfGroup(readGroups.get(i)), i);
        }
    }

    @Test
    public void testTwoMappedPairs() {
        final FindReadGroupDistributionsTester tester = new FindReadGroupDistributionsTester();
        tester.addMappedPair(1, 1, 100, false, false, DEFAULT_BASE_QUALITY);
        tester.addMappedPair(1, 1, 100, true, true, DEFAULT_BASE_QUALITY); // duplicate!!!
        tester.runTest();
    }
}
