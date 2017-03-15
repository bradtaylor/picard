/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.analysis.library;


import com.google.gson.Gson;
import htsjdk.samtools.*;
import htsjdk.samtools.filter.*;
import htsjdk.samtools.util.*;
import picard.PicardException;
import picard.analysis.replicates.IndependentReplicateMetric;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Alpha;
import picard.filter.PairsMappingQualityFilter;
import picard.filter.ProperPairsFilter;
import picard.sam.markduplicates.util.AbstractMarkDuplicatesCommandLineProgram;
import picard.sam.markduplicates.util.AbstractOpticalDuplicateFinderCommandLineProgram;
import picard.sam.markduplicates.util.LibraryIdGenerator;
import picard.sam.markduplicates.util.ReadEnds;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.*;


/**
 * <p>A CLP that records the pattern of duplicate dispersal across read groups. Consumes a set of BAM inputs, assembles
 * duplicate sets, and checks the read group of each set member. Sets are represented as a pattern of read groups,
 * and multiple sets displaying the same pattern are collapsed to form a count map. Patterns are represented by a list
 * of numerical index keys. A read group key file associates each index with full {@link SAMReadGroupRecord} information.
 * This file is output alongside the pattern information, both in Json format.
 * <p>
 * The intended purpose of collecting this information is to provide features for the {@link CheckLibraryIdentity}
 * CLP. Caching feature information is convenient during development, and is needed if model fitting etc. will be applied
 * in an alternate environment. However, the ultimate intention is to merge these two programs.</p>
 * <p>
 * Multiple filters are applied to remove problematic reads/pairs before duplicate set sorting. In addition, the method tracks
 * and removes optical duplicates from the final read group pattern outputs. Both of these features may be disabled.</p>
 * <p>
 * NOTE: This class is very much in alpha stage, and still under heavy development.</p>
 *
 * @author Bradley Taylor
 *
 */

@CommandLineProgramProperties(
        usage = "<p>A CLP that records the pattern of duplicate dispersal across read groups. Consumes a set of BAM inputs, assembles" +
                " duplicate sets, and checks the read group of each set member. Sets are represented as a pattern of read groups," +
                " and multiple sets displaying the same pattern are collapsed to form a count map. Patterns are represented by a list" +
                " of numerical index keys. A read group key file associates each index with full {@link SAMReadGroupRecord} information." +
                " This file is output alongside the pattern information, both in Json format." +
                " <p>" +
                " The intended purpose of collecting this information is to provide features for the {@link CheckLibraryIdentity}" +
                " CLP. Caching feature information is convenient during development, and is needed if model fitting etc. will be applied" +
                " in an alternate environment. However, the ultimate intention is to merge these two programs.</p>" +
                " <p>" +
                " Multiple filters are applied to remove problematic reads/pairs before duplicate set sorting. In addition, the method tracks" +
                " and removes optical duplicates from the final read group pattern outputs. Both of these features may be disabled.</p>" +
                " <p>" +
                " NOTE: This class is very much in alpha stage, and still under heavy development.</p>",
        usageShort = "TBD",
        programGroup = Alpha.class
)
public class FindReadGroupDistributions extends AbstractOpticalDuplicateFinderCommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "One or more input SAM or BAM files to analyze. Must be coordinate sorted.", minElements = 1) // TODO- is that true about sorting?
    public List<File> INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Basename for the table and read-group-key output json files.")
    public String OUTPUT;

    @Option(shortName = "MQ", doc = "Minimal value for the mapping quality of the reads to be examined.", optional = true)
    public Integer MINIMUM_MQ = 30;

    @Option(doc = "Number of duplicate sets to examine before stopping. Duplicate set formation (inc. sorting) will " +
            "still occur; this only halts the subsequent walk-through.", optional = true)
    public Integer STOP_AFTER = 0;

    @Option(doc = "Whether to pre-filter reads before capturing read group patterns. Defaults to true, as doing so should improve accuracy.", optional = true)
    public boolean FILTER_READS = true;

    private static final Log log = Log.getInstance(FindReadGroupDistributions.class);
    private static final Gson gson = new Gson();

    private LibraryIdGenerator libraryIdGenerator;

    protected static String readGroupKeySuffix = ".read_groups.json";
    protected static String tableOutputSuffix = ".table.json";

    protected int totalReads = 0;
    protected int duplicateSetsExamined = 0;
    protected double setsToReadsRatio = 0.0;
    protected int readsExamined = 0;
    protected int opticalDuplicateReadsIgnored = 0;

    @Override
    protected String[] customCommandLineValidation() {
        if (READ_NAME_REGEX == null) {
            log.warn("Skipped optical duplicate cluster discovery; results may be biased by noise introduced by on-sequencer processes following library construction.");
        }
        return super.customCommandLineValidation();
    }

    // Read ends is abstract, but we don't need any additional functionality
    private class ReadEndsForFindReadGroupDistributions extends ReadEnds {}

    /**
     * A BiMap for short ID : SAMReadGroupRecord group, implemented as two co-initialized HashMaps.
     * NB: Guava has a true BiMap, but I did not want to add the dependency.
     */
    protected static class ReadGroupIdManager {
        private final Map<Short, SAMReadGroupRecord> readGroupRecordMap = new HashMap<>();
        private final Map<SAMReadGroupRecord, Short> readGroupIndexMap = new HashMap<>();

        ReadGroupIdManager(SAMFileHeader header) {
            final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
            short index = -1; // read groups should be zero indexed for building distance matrices :)
            for (SAMReadGroupRecord readGroup : readGroups) {
                this.readGroupRecordMap.put(++index, readGroup);
                this.readGroupIndexMap.put(readGroup, index);
            }
        }

        short getIdOfGroup(SAMReadGroupRecord readGroup) {
            return this.readGroupIndexMap.get(readGroup);
        }

        SAMReadGroupRecord getGroupFromId(short id) {
            return this.readGroupRecordMap.get(id);
        }

        Map<Short, SAMReadGroupRecord> getReadGroupRecordMap() {
            return this.readGroupRecordMap;
        }
    }

    @Override
    protected int doWork() {

        // Get a merging iterator for all of our inputs, and it's merged header, in order to create a duplicate set iterator
        final MergingSamRecordIterator mergingFilteringIterator = getMergingFilteringIterator();
        final SAMFileHeader mergedHeader = mergingFilteringIterator.getMergedHeader();
        log.info("Finished opening outputs and constructing underlying iterators");
        final DuplicateSetIterator duplicateSets = new DuplicateSetIterator(mergingFilteringIterator,
                                                                            mergedHeader,
                                                                            false,
//                                                                            new SAMRecordDuplicateComparator()); TODO: replace InsertComparator with an appropriately modified form of SAMRecordDuplicateComparator
                                                                            new InsertComparator());

        System.out.println(mergedHeader.getReadGroups().toString());

        // Map each read group from the merged header to a numerical index, to serve as a lookup table for the feature data.
        final ReadGroupIdManager groupInfo = new ReadGroupIdManager(mergedHeader);
        writeReadGroupInformation(groupInfo.getReadGroupRecordMap());

        System.out.println(groupInfo.getReadGroupRecordMap().toString());
        System.out.println(groupInfo.getReadGroupRecordMap().entrySet().toString());

        // Feature data is a map associating (A) an observed pattern of read group IDs, to (B) the count of how often that pattern was observed.
        // Keys will a list of shorts, corresponding to the lookup-table keys we built above. So 'RG1,RG1,RG2,RG4' would be '0,0,1,3'
        // We want the pattern of read groups, so duplicates are allowed. The group-index list is sorted.
        // NB: There are a number of ways to reduce the representation of the RGID pattern. This is a first pass.
        Map<List<Short>, Integer> featureData = new HashMap<>();

        // Associates each library name (string) with a libraryId (short). We don't need all the metrics accessors, but the class is required for the OpticalDuplicateTracker
        libraryIdGenerator = new LibraryIdGenerator(mergedHeader);

        // Set up logging information
        final ProgressLogger progress = new ProgressLogger(log, 10000000, "examined", "duplicate sets");

        log.info("Starting iteration on duplicate sets");

        Set<String> opticalDuplicateMatesSeen = new HashSet<>();

        while (duplicateSets.hasNext()) {
            final DuplicateSet set = duplicateSets.next();
            log.debug("NEW SET");

            // Do some logging
            duplicateSetsExamined++;
            final int setSize = set.size();
            totalReads += setSize;

            // Get the set's representative record (i.e. the non-duplicate). Mostly for logging, but also as a keeper to pass to the OpticalDuplicateFinder
            // Best-record selection uses a more stringent comparison based on duplication score, but orientations are collapsed.
            // NB: Even in our optical duplicate marking we don't particularly care about which record this is.
            //     A duplicate set can contain multiple subsets of optical duplicates, and for each optical duplicate subset, the finder arbitrarily picks a non-optical-duplicate
            //     The representative is passed in to ensure that the read we wish to call not-a-duplicate is not then called an OPTICAL duplicate. We are not concerned about this, because...
            // NB: Optical duplicates are required to be of the same read group
            final SAMRecord setRep = set.getRepresentative();
            progress.record(setRep);

            // Set up for optical duplicate tracking
            final boolean doOpticalDuplicateTracking = (this.READ_NAME_REGEX != null) &&
                    isPairedAndBothMapped(setRep); // &&
            final Set<String> duplicateReadEndsSeen = new HashSet<>();
            final List<ReadEnds> duplicateReadEnds = new ArrayList<>();
            final List<ReadEnds> allReadEnds = new ArrayList<>();


            for (final SAMRecord read : set.getRecords()) {

                final ReadEndsForFindReadGroupDistributions readEnd = getProperlyConfigeredReadEnds(read, groupInfo, libraryIdGenerator);

                allReadEnds.add(readEnd); // NB: The full list of records, whether we consider them for optical duplicates or not

                String whichMate = read.getFirstOfPairFlag() ? "first" : "second";
                log.debug("read name: " + read.getReadName());
                log.debug("read is " + whichMate + " of pair");

                // To track optical duplicates, store a set of locations for each read. We expect each end of a mapped pair to be in a separate duplicate set.
                //  We care about orientation relative to the first end of the pair for optical duplicate tracking, which is more stringent than PCR duplicate tracking.
                if (doOpticalDuplicateTracking &&
                        isPairedAndBothMapped(read) &&
                        !duplicateReadEndsSeen.contains(read.getReadName())) {
                    duplicateReadEnds.add(readEnd);
                    duplicateReadEndsSeen.add(read.getReadName());
                }
            }
            log.debug("set size: " + set.size());
            log.debug("duplicateReadEnds size: " + duplicateReadEnds.size());
            log.debug("duplicateReadEndsSeen size: " + duplicateReadEndsSeen.size());
            log.debug("allReadEnds size: " + allReadEnds.size());

            // Track the optical duplicates. NB: This sets the transient isOpticalDuplicate flag on each ReadEnds.
            if (this.READ_NAME_REGEX != null && 1 < duplicateReadEnds.size()) {
                AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(duplicateReadEnds, duplicateReadEnds.get(0), opticalDuplicateFinder, libraryIdGenerator);
            }

            // Collect the observed pattern of read groups and the count of reads, ignoring those marked as optical duplicates
            List<Short> readGroupIndexes = new ArrayList<>();
            int libraryRecordsExamined = 0;
            for (ReadEnds end : allReadEnds) {
                if (end.isOpticalDuplicate) {
                    opticalDuplicateReadsIgnored++;

                } else {
                    short readGroupIndex = end.getReadGroup();
                    if (readGroupIndex >= 0) { // Ignore reads with -1 for read group Id. This most likely occurs because phys info could not be parsed from read name.
                        readGroupIndexes.add(end.getReadGroup());
                        libraryRecordsExamined++;
                    }
                }
            }
            readsExamined += libraryRecordsExamined;

            // Add or update the read count associated with this particular pattern of read groups
            Collections.sort(readGroupIndexes);
            log.debug("read group pattern is: " + readGroupIndexes);
            if (!readGroupIndexes.isEmpty()) {
                featureData.merge(readGroupIndexes, libraryRecordsExamined, Integer::sum);
            }

            if (STOP_AFTER > 0 && progress.getCount() > STOP_AFTER) break;
        }
        setsToReadsRatio = (duplicateSetsExamined / (double)totalReads);
        log.info("Reads examined: " + totalReads);
        log.info("Duplicate sets examined: " + duplicateSetsExamined);
        log.info("Sets / Read ratio: " + setsToReadsRatio);
        log.info("Library records examined: " + readsExamined);
        log.info("Optical duplicates ignored: " + opticalDuplicateReadsIgnored);
        log.info("Found " + (libraryIdGenerator.getNumberOfOpticalDuplicateClusters()) + " optical duplicate clusters.");
        writeTableOutput(featureData);
        duplicateSets.close();

        return 0;
    }

    /*
     * Build a ReadEnds, setting group and library info, as well as the orientation info needed for optical duplicate tracking.
     */
    private ReadEndsForFindReadGroupDistributions getProperlyConfigeredReadEnds(SAMRecord read, ReadGroupIdManager groupInfo, LibraryIdGenerator libraryIdGenerator) {
        final ReadEndsForFindReadGroupDistributions readEnd = new ReadEndsForFindReadGroupDistributions();

        if (read.getFirstOfPairFlag()) {
            readEnd.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(read.getReadNegativeStrandFlag(), read.getMateNegativeStrandFlag());
        } else {
            readEnd.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(read.getMateNegativeStrandFlag(), read.getReadNegativeStrandFlag());
        }
        if (READ_NAME_REGEX==null || opticalDuplicateFinder.addLocationInformation(read.getReadName(), readEnd)) {
            SAMReadGroupRecord readGroup = read.getReadGroup();
            if (null != read.getReadGroup()) {
                final short libIndex = libraryIdGenerator.getLibraryId(read);
                readEnd.setLibraryId(libIndex);
                readEnd.setReadGroup(groupInfo.getIdOfGroup(readGroup)); // NB: Store an ID from our map instead of the actual group
            }
        }
        return readEnd;
    }

    private static boolean isPairedAndBothMapped(final SAMRecord record) {
        return record.getReadPairedFlag() &&
                !record.getReadUnmappedFlag() &&
                !record.getMateUnmappedFlag();
    }

    /**
     * For each input, create a SamReader and obtain a FilteringSamIterator.
     * Pass the all the FilteringSamIterators into a MergingSamRecordIterator
     *
     * Simultaneously, get the merged header.
     * @return a SamHeaderAndIterator containing the merged header and merging iterator
     */
    private MergingSamRecordIterator getMergingFilteringIterator() {
        IOUtil.assertFilesAreReadable(INPUT);
        final List<SAMFileHeader> headers = new ArrayList<>(INPUT.size());
        Map<SamReader, CloseableIterator<SAMRecord>> iterators = new HashMap<>();
        for (final File inFile : INPUT) {
            final SamReader reader = SamReaderFactory.makeDefault().open(inFile);
            final SAMRecordIterator samRecordIterator = reader.iterator();
            final FilteringSamIterator filteredSamRecordIterator = new FilteringSamIterator(samRecordIterator, new AggregateFilter(getSamRecordFilters()));
            iterators.put(reader, filteredSamRecordIterator);

            final SAMFileHeader header = reader.getFileHeader();
            headers.add(header);
        }
        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(headers.get(0).getSortOrder(), headers, false);
        final MergingSamRecordIterator iterator = new MergingSamRecordIterator(headerMerger, iterators, false);
        log.info("Constructed filtering merging SAM/BAM iterator header");
        return iterator;
    }

    /*
     * return the SAM record filters we wish to apply when iterating over reads.
     *
     * return an empty list if filtering is disabled at the command line
     *
     * Our desired filters treat pairs the same way,
     *   i.e. if either end is filtered out, it's mate must be filtered out as well.
     *   This despite our data likely being in coordinate order.
     *   We are interested in the whole molecule. We don't want to aberrantly create fragments
     *
     * Filters:
     *   Exclude any reads not aligned in a proper pair
     *   No secondary or supplementary alignments (not independent evidence of the insert molecule!)
     *   Mapping Quality, which is used as a means of excluding so-called "reference duplicates". Preferable over an interval filter because MQ is simple and reference-agnostic
     *     For justification, see: https://github.com/broadinstitute/palantir/blob/71b5b9345fd6105b1783a0a4514aef5adcbc35c1/Analysis/duplicate_set_sizes/Readme.md
     */
    private List<SamRecordFilter> getSamRecordFilters() {
        if (!FILTER_READS) {
            return new ArrayList<>();
        }

        return CollectionUtil.makeList(
                new FailsVendorReadQualityFilter(), // NB: Both reads in a pair will posses a 'failed vendor quality' flag if either is low-quality.
                new SecondaryOrSupplementaryFilter(),
                new ProperPairsFilter(),
                new PairsMappingQualityFilter(MINIMUM_MQ)
        );
    }

    /*
     * Write full read group header entries in JSON format to an OUTPUT-derived file, as a key to the table data
     */
    private void writeReadGroupInformation(final Map<Short, SAMReadGroupRecord> readGroups) {
        log.info("Writing read group information");
        String output = gson.toJson(readGroups);
        writeStringOutput(output, new File(OUTPUT + readGroupKeySuffix));
    }

    /*
     * Convert {read group indexes -> count} map into JSON representation and write to an OUTPUT-derived file
     */
    private void writeTableOutput(final Map<List<Short>, Integer> featureData) {
            log.info("Writing table data");
            String output = gson.toJson(featureData);
            writeStringOutput(output, new File(OUTPUT + tableOutputSuffix));
    }

    /*
     * Write generic string data to a given file
     */
    private void writeStringOutput(String output, File file) {
        IOUtil.assertFileIsWritable(file);
        try {
            final FileOutputStream fileOutputStream = new FileOutputStream(file);
            final OutputStreamWriter outputStreamWriter = new OutputStreamWriter(fileOutputStream);
            outputStreamWriter.append(output);
            outputStreamWriter.close();
            fileOutputStream.close();
        } catch (Exception e) {
            throw new PicardException("Error writing table output: ", e);
        }
    }

    public static void main(final String[] args) {
        new FindReadGroupDistributions().instanceMainWithExit(args);
    }
}
