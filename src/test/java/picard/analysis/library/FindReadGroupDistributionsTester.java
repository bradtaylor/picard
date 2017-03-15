package picard.analysis.library;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import org.testng.Assert;
import picard.sam.testers.SamFileTester;

import java.io.File;
import java.util.List;

/**
 * Created by bradt on 1/28/17.
 */
public class FindReadGroupDistributionsTester extends SamFileTester {
    private String outputBasename;
    private File tableOutput;
    private File readGroupKey;

    @Override
    public String getCommandLineProgramName() {
        return FindReadGroupDistributions.class.getSimpleName();
    }

    // NB: InsertComparator.scoringStrategy is same as SAMRecordSetBuilder.DEFAULT_DUPLICATE_SCORING_STRATEGY
    public FindReadGroupDistributionsTester() {
        super(50, true);

        // Presumably can do some add-args here based on parameters we pass in, if we need args other than default
    }

    public String[] makePicardCommandLineArgs(final List<String> args) {
        final String[] picardCommandLineArgs = new String[args.size()];
        int i = 0;
        for (final String arg : args) {
            picardCommandLineArgs[i++] = arg;
        }
        return picardCommandLineArgs;
    }

    // NB: This is duplicated directly from SamFileTester, where it has private access.
    // I do not love this; it was done out of expediency when this was part of a separate repo.
    // TODO - modify SamFileTester
    private File createInputFile() {
        // Create the input file
        final File input = new File(getOutputDir(), "input.sam");
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(samRecordSetBuilder.getHeader(), true, input);
        samRecordSetBuilder.getRecords().forEach(writer::addAlignment);
        writer.close();
        return input;
    }


    /**
     * Main test execution method
     * NB: I am not using the typical SamFileTester run-test here for two reasons:
     *
     * 1) Output is a string prefix for naming two other outputs, not a File as in SamFileTester
     * 2) Because I wish to capture an instance of the CLP class, such that I can access the protected summary stats it maintains.
     */
    @Override
    public void runTest() {

        final File input = createInputFile();

        setOutputs();

        addArg("INPUT=" + input.getAbsoluteFile());
        addArg("OUTPUT=" + outputBasename);

        FindReadGroupDistributions finder = new FindReadGroupDistributions();
        Assert.assertEquals(finder.instanceMain(makePicardCommandLineArgs(getArgs())), 0);

    }

    private void setOutputs() {
        outputBasename = (new File(getOutputDir(), "output")).getAbsolutePath();
        tableOutput = new File(outputBasename + FindReadGroupDistributions.tableOutputSuffix);
        readGroupKey = new File(outputBasename + FindReadGroupDistributions.readGroupKeySuffix);
    }

    @Override
    protected void test() {}
}
