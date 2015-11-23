/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package org.embl.gbcs.je.jeduplicates;

import htsjdk.samtools.DuplicateScoringStrategy;
import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.SortingLongCollection;
import htsjdk.samtools.util.StringUtil;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.embl.cg.utilitytools.utils.FileUtil;
import org.embl.gbcs.je.Jexception;

import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.SamOrBam;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.AbstractMarkDuplicatesCommandLineProgram;
import picard.sam.markduplicates.util.LibraryIdGenerator;
import picard.sam.markduplicates.util.ReadEnds;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicatesMap;

/**
 * A better duplication marking algorithm that handles all cases including clipped
 * and gapped alignments ; and molecular barcodes
 *
 * @author Charles Girardot (based on original code from Tim Fennell)
 */
@CommandLineProgramProperties(
		usage = "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules taking into"
				+ " account molecular barcodes (Unique Molecular Identifiers or UMIs) found in read header. " 
				+ "All records are then either written to the output file with the duplicate records flagged or trashed.\n"
				+"Example :\n"
				+"\t je markdupes INPUT=file_with_dupes.bam OUPUT=result.bam MM=1", 
		usageShort = "je markdupes INPUT=file_with_dupes.bam OUPUT=result.bam MM=1", 
		programGroup = SamOrBam.class
		)
public class MarkDuplicatesWithMolecularCode extends AbstractMarkDuplicatesCommandLineProgram {
	
	//a constant used to store all reads (from the same positional duplicate group) that have too many Ns in their molecular code 
	private static final String DUPES_GROUP_NAME_UNDEF = "UNDEF";


	private final Log log = Log.getInstance(MarkDuplicatesWithMolecularCode.class);

	
	@Option(shortName = "MM", optional = false,
			printOrder = 2010, //let the first two options from AbstractMarkDuplicatesCommandLineProgram go first ie use order > 2000 and < 3000 
			doc="Number of MisMatches (inclusive) to still consider two Unique Molecular Identifiers (UMIs) identical i.e. "
					+ "this option buffers for sequencing errors."
					+ "Indeed, in case of a sequencing error, 2 duplicate reads would not be considered duplicates anymore."
					+ "Note that N are not considered mismatches during comparison ie ATTNGG and NTTANG are seen as"
					+ " the same barcode and these two reads would be flagged duplicates."
					+ "This option takes a single value even when several barcodes are present (see SLOTS). "
					+ "Note that when declaring several barcodes (see SLOTS) AND providing a predefined set "
					+ "of barcodes (see BC option), the MM value is applicable in each lookup. When a predefined set "
					+ "of barcodes is NOT given, the different barcodes (SLOTS) are concatenated first and the MM value"
					+ " is therefore considered *overall* as the concatenated code is seen as a unique code.\n"
					+ "MM=null is like MM=0\n"
					+ "Use the minimum Hamming distance of the original barcode set (if applicable)."
					
			)
	public Integer MISMATCHES = null;
	
	@Option(shortName = "MAX_N", optional = true,
			printOrder = 2020,
			doc="Maximum number of Ns a molecular code can contain (inclusive). Above this value, reads are placed in a UNDEF group."
					+ "More precisely, these 'too degenarate' codes will not :\n"
					+ "\t * be compared to the list of predefined codes [predefined code list situation ie BC option given] nor \n"
					+ "\t * be considered as a potential independent code [no predefined code list situation ie BC option not given]\n"
					+"Default value is the MISMATCHES number.\n "
					+ "Note that when declaring several barcodes (see SLOTS) AND providing a predefined set "
					+ "of barcodes (see BC option), the MAX_N value is applicable to each barcode. When a predefined set "
					+ "of barcodes is NOT given, the different barcodes (SLOTS) are concatenated first and the MAX_N value"
					+ " is therefore considered *overall*."
					
			)
	public Integer MAX_NUMBER_OF_N = null;
	
	@Option(shortName = "SLOTS", optional = true,
			printOrder = 2030,
			doc="Where to find the UMIs (and only the UMIs) in the read name once read name has been tokenized using the SPLIT character (e.g. ':'). \n"
					+ "By default, the UMI is considered to be found at the end of the read header i.e. after the last ':'. "
					+ "Use this option to indicate other or additional UMI positions (e.g. multiple UMIs present in read header.\n"
					+ "IMPORTANT : counting starts at 1 and negative numbers can be used to start counting from the end.\n"
					+ "For example, consider the following read name that lists 3 different barcodes in the end : \n"
					+ "\t HISEQ:44:C6KC0ANXX:8:2112:20670:79594:CGATGTTT:GATCCTAG:AAGGTACG \n"
					+ "\t to indicate that the three barcodes are molecular codes, use \n"
					+ "\t\tSLOTS=-1 SLOTS=-2 SLOTS=-3\n"
					+ "\t if only the 2 last ones should be considered (the third one being a sample encoding barcode), use \n"
					+ "\t\tSLOTS=-1 SLOTS=-2\n"
					)
	public List<Integer> SLOTS = null;
	
	//if user deactivate UMIs using explicit SLOTS=null, we ll turn this internal flag to false
	private boolean USE_UMIS = true;
	
	@Option(shortName = "BC", optional = true,
			printOrder = 2040,
			doc="Pre-defined list of UMIs that can be expected. "
					+ "Format: one column text file, one barcode per line. "
					+ "All UMIs MUST have the same length. "
			)
	public File BARCODE_FILE = null;
	
	
	@Option(shortName="T", optional = true,
			printOrder = 5050,
			doc="Should barcode information be removed from read names in the output BAM ? This is usefull to save storage space.\n" 
			)
	public boolean TRIM_HEADERS = false;

	@Option(shortName = "TSLOTS", optional = true,
			printOrder = 5060,
			doc="Where to find *all* barcode(s) (i.e. sample encoding and UMIs) in the read name once has been tokenized"
					+ " using the SPLIT character (e.g. ':'). \n"
					+ "This option is only considered when TRIM_HEADERS=true. "
					+ "When TSLOTS is ommited while TRIM_HEADERS=true, the values of SLOTS apply.\n"
					+ "IMPORTANT : counting starts at 1 and negative numbers can be used to start counting from the end.\n"
					+ "See SLOT help for examples." 
			)
	public List<Integer> TSLOTS = null;
	
	@Option(shortName = "SPLIT", optional = true,
			printOrder = 5070,
			doc="Character to use to split up the read header line, default is ':'." 
					)
	public String SPLIT_CHAR = ":";
	
	
	@Option(shortName = "MAX_FILE_HANDLES",
			printOrder = 20000, //throw this option to the end
			doc = "Maximum number of file handles to keep open when spilling read ends to disk. " +
					"Set this number a little lower than the per-process maximum number of file that may be open. " +
			"This number can be found by executing the 'ulimit -n' command on a Unix system.")
	public int MAX_FILE_HANDLES_FOR_READ_ENDS_MAP = 8000;

	@Option(printOrder = 21000, //throw this option to the end
			doc = "This number, plus the maximum RAM available to the JVM, determine the memory footprint used by " +
			"some of the sorting collections.  If you are running out of memory, try reducing this number.")
	public double SORTING_COLLECTION_SIZE_RATIO = 0.25;

	
	
	/*
	 * the MolecularBarcodeFinder to extract molecular code from read names 
	 */
	private MolecularBarcodeFinder molecularBarcodeFinder = null;
	
	/*
	 * Set of "random" code to be expected, will be filled only when BC=file is given 
	 */
	private Set<String> predefinedBarcodes = null;
	// all barcode have same length
	private Integer barcodeLength = null;
	
	private SortingCollection<ReadEndsForMarkDuplicatesWithMolecularCode> pairSort;
	private SortingCollection<ReadEndsForMarkDuplicatesWithMolecularCode> fragSort;
	private SortingLongCollection duplicateIndexes;
	private int numDuplicateIndices = 0;

	private LibraryIdGenerator libraryIdGenerator = null; // this is initialized in buildSortedReadEndLists

	
	
	public MarkDuplicatesWithMolecularCode() {
		super();
		DUPLICATE_SCORING_STRATEGY = ScoringStrategy.SUM_OF_BASE_QUALITIES;
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		new MarkDuplicatesWithMolecularCode().instanceMainWithExit(args);
	}

	/**
	 * Main work method.  Reads the BAM file once and collects sorted information about
	 * the 5' ends of both ends of each read (or just one end in the case of pairs).
	 * Then makes a pass through those determining duplicates before re-reading the
	 * input file and writing it out with duplication flags set correctly.
	 */
	protected int doWork() {
		IOUtil.assertInputsAreValid(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		IOUtil.assertFileIsWritable(METRICS_FILE);

		/*
		 * Molecular Barcode Support 
		 * Initialize the MolecularBarcodeFinder
		 */
		if(SLOTS == null || SLOTS.size() == 0){
			SLOTS = Arrays.asList(new Integer[]{-1});
		}
		else if(SLOTS.size() == 1 && SLOTS.get(0) == 0){
			//user wants to run without UMIs
			SLOTS = new ArrayList<Integer>(); 
			USE_UMIS = false;
		}
		
		if(TSLOTS == null || TSLOTS.size() == 0){
			TSLOTS = SLOTS;
			if(!USE_UMIS){
				//force TRIM_HEADERS to false to avoid problems
				TRIM_HEADERS = false;
			}
		}
		if(USE_UMIS){
			molecularBarcodeFinder = new MolecularBarcodeFinder(SPLIT_CHAR, SLOTS.toArray(new Integer[SLOTS.size()]));
		}else{
			molecularBarcodeFinder = new MolecularBarcodeFinder(SPLIT_CHAR, null);
		}
		
		if(MAX_NUMBER_OF_N == null){
			MAX_NUMBER_OF_N = MISMATCHES;
		}
		
		/*
		 * Do we have a pre-defined list of "random" barcode ?
		 */
		if(BARCODE_FILE!=null && USE_UMIS){
			IOUtil.assertFileIsReadable(BARCODE_FILE);
			
			//load up barcodes 
			try {
				predefinedBarcodes = new TreeSet<String>( FileUtil.getTokenListFromSingleColumnFile(BARCODE_FILE.getAbsolutePath()) );
				barcodeLength = predefinedBarcodes.iterator().next().length();
				for (String bc : predefinedBarcodes) {
					if(bc.length() != barcodeLength){
						throw new IllegalArgumentException("Barcodes read from "+BARCODE_FILE.getName()+"do not have all the same length.");
					}
				}
			} catch (IOException e) {
				 throw new IllegalArgumentException("Cannot read barcode list from "+BARCODE_FILE.getName()+". Check format (one column text file, one barcode per line).");
			}
			
			//sounds good make sure SLOTS all give expected length UMIs or die
			try{
				validateUMISlots(barcodeLength);
			}catch(Jexception jexc){
				//something is wrong with extracted UMIs
				log.error(jexc, jexc.getMessage());
				//exit
				return 1;
			}
			
		}		
		
		reportMemoryStats("Start of doWork");
		log.info("Reading input file and constructing read end information.");
		try{
			buildSortedReadEndLists();
		}catch(Jexception jexc){
			//something is wrong with extracted UMIs
			log.error(jexc, jexc.getMessage());
			//exit
			return 1;
		}
		reportMemoryStats("After buildSortedReadEndLists");
		generateDuplicateIndexes();
		reportMemoryStats("After generateDuplicateIndexes");
		log.info("Marking " + this.numDuplicateIndices + " records as duplicates.");

		if (this.READ_NAME_REGEX == null) {
			log.warn("Skipped optical duplicate cluster discovery; library size estimation may be inaccurate!");
		} else {
			log.info("Found " + (this.libraryIdGenerator.getNumberOfOpticalDuplicateClusters()) + " optical duplicate clusters.");
		}

		final SamHeaderAndIterator headerAndIterator = openInputs();
		final SAMFileHeader header = headerAndIterator.header;

		final SAMFileHeader outputHeader = header.clone();
		outputHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
		for (final String comment : COMMENT) outputHeader.addComment(comment);

		// Key: previous PG ID on a SAM Record (or null).  Value: New PG ID to replace it.
		final Map<String, String> chainedPgIds = getChainedPgIds(outputHeader);

		final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(outputHeader,
				true,
				OUTPUT);

		// Now copy over the file while marking all the necessary indexes as duplicates
		long recordInFileIndex = 0;
		long nextDuplicateIndex = (this.duplicateIndexes.hasNext() ? this.duplicateIndexes.next() : -1);

		final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Written");
		final CloseableIterator<SAMRecord> iterator = headerAndIterator.iterator;
		/*
		 * help for read name trimming (option T) if needed
		 */
		Set<Integer> readNameSlotsToIgnore = new TreeSet<Integer>();
		
		while (iterator.hasNext()) {
			final SAMRecord rec = iterator.next();
			if (!rec.isSecondaryOrSupplementary()) {
				final String library = libraryIdGenerator.getLibraryName(header, rec);
				DuplicationMetrics metrics = libraryIdGenerator.getMetricsByLibrary(library);
				if (metrics == null) {
					metrics = new DuplicationMetrics();
					metrics.LIBRARY = library;
					libraryIdGenerator.addMetricsByLibrary(library, metrics);
				}

				// First bring the simple metrics up to date
				if (rec.getReadUnmappedFlag()) {
					++metrics.UNMAPPED_READS;
				} else if (!rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) {
					++metrics.UNPAIRED_READS_EXAMINED;
				} else {
					++metrics.READ_PAIRS_EXAMINED; // will need to be divided by 2 at the end
				}


				if (recordInFileIndex == nextDuplicateIndex) {
					rec.setDuplicateReadFlag(true);

					// Update the duplication metrics
					if (!rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) {
						++metrics.UNPAIRED_READ_DUPLICATES;
					} else {
						++metrics.READ_PAIR_DUPLICATES;// will need to be divided by 2 at the end
					}

					// Now try and figure out the next duplicate index
					if (this.duplicateIndexes.hasNext()) {
						nextDuplicateIndex = this.duplicateIndexes.next();
					} else {
						// Only happens once we've marked all the duplicates
						nextDuplicateIndex = -1;
					}
				} else {
					rec.setDuplicateReadFlag(false);
				}
			}
			recordInFileIndex++;

			/*
			 * With this option the user wants to eliminate the barcode slots from the read name.
			 * This is handy to clip back all (molecular) barcodes from the name 
			 */
			if(TRIM_HEADERS){
				String readName = rec.getReadName();
				String[] tokens = readName.split(SPLIT_CHAR);
				if(readNameSlotsToIgnore.size() == 0){
					//set
					for (int slot : TSLOTS) {
						if(slot>=0){
							readNameSlotsToIgnore.add(slot);
						}else{
							//from end
							readNameSlotsToIgnore.add( (tokens.length + slot));
						}
					}
				}
				StringBuffer nname = new StringBuffer();
				boolean addTok = false;
				for (int i = 0; i < tokens.length; i++) {
					if(!readNameSlotsToIgnore.contains(i)){
						if(addTok){
							nname.append(":");
						}else{
							addTok=true;
						}
						nname.append(tokens[i]);
					}
				}
				rec.setReadName(nname.toString());
			} 
			
			
			if (!this.REMOVE_DUPLICATES || !rec.getDuplicateReadFlag()) {
				if (PROGRAM_RECORD_ID != null) {
					rec.setAttribute(SAMTag.PG.name(), chainedPgIds.get(rec.getStringAttribute(SAMTag.PG.name())));
				}
				out.addAlignment(rec);
				progress.record(rec);
			}
		}

		// remember to close the inputs
		iterator.close();

		this.duplicateIndexes.cleanup();

		reportMemoryStats("Before output close");
		out.close();
		reportMemoryStats("After output close");

		// Write out the metrics
		finalizeAndWriteMetrics(libraryIdGenerator);

		return 0;
	}

	
	
	/**
	 * When using a UMI file, checks that all UMI slots contain UMI of expected length or throws a Jexception
	 */
	private void validateUMISlots(int UMILength) {
		
		final SamHeaderAndIterator headerAndIterator = openInputs();
		final CloseableIterator<SAMRecord> iterator = headerAndIterator.iterator;
		final SAMRecord rec = iterator.next();
		String readname = rec.getReadName();
		iterator.close();
		
		//now check
		if( !molecularBarcodeFinder.validateUMILength(readname, UMILength)){
			throw new Jexception( "One or more of the UMIs indicated by the SLOTS option (SLOTS=["+org.embl.cg.utilitytools.utils.StringUtil.mergeIntList(SLOTS, ",")+"]"
					+") do not have the expected length of "+UMILength+". Please check your command line options." );
		}	
		
	}

	
	
	/**
	 * package-visible for testing
	 */
	long numOpticalDuplicates() { return ((long) this.libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap().getSumOfValues()); } // cast as long due to returning a double

	/** Print out some quick JVM memory stats. */
	private void reportMemoryStats(final String stage) {
		System.gc();
		final Runtime runtime = Runtime.getRuntime();
		log.info(stage + " freeMemory: " + runtime.freeMemory() + "; totalMemory: " + runtime.totalMemory() +
				"; maxMemory: " + runtime.maxMemory());
	}

	/**
	 * Goes through all the records in a file and generates a set of ReadEndsForMarkDuplicatesWithMolecularCode objects that
	 * hold the necessary information (reference sequence, 5' read coordinate) to do
	 * duplication, caching to disk as necessary to sort them.
	 */
	private void buildSortedReadEndLists() {
		final int maxInMemory = (int) ((Runtime.getRuntime().maxMemory() * SORTING_COLLECTION_SIZE_RATIO) / ReadEndsForMarkDuplicatesWithMolecularCode.SIZE_OF);
		log.info("Will retain up to " + maxInMemory + " data points before spilling to disk.");

		this.pairSort = SortingCollection.newInstance(ReadEndsForMarkDuplicatesWithMolecularCode.class,
				new ReadEndsForMarkDuplicatesWithMolecularCodeCodec(),
				new ReadEndsMDComparator(),
				maxInMemory,
				TMP_DIR);
 
		this.fragSort = SortingCollection.newInstance(ReadEndsForMarkDuplicatesWithMolecularCode.class,
				new ReadEndsForMarkDuplicatesWithMolecularCodeCodec(),
				new ReadEndsMDComparator(),
				maxInMemory,
				TMP_DIR);

		final SamHeaderAndIterator headerAndIterator = openInputs();
		final SAMFileHeader header = headerAndIterator.header;
		final ReadEndsForMarkDuplicatesMap tmp = new DiskBasedReadEndsForMarkDuplicatesWithMolecularCodeMap(MAX_FILE_HANDLES_FOR_READ_ENDS_MAP);
		long index = 0;
		final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");
		final CloseableIterator<SAMRecord> iterator = headerAndIterator.iterator;

		if (null == this.libraryIdGenerator) {
			this.libraryIdGenerator = new LibraryIdGenerator(header);
		}

		while (iterator.hasNext()) {
			final SAMRecord rec = iterator.next();

			// This doesn't have anything to do with building sorted ReadEnd lists, but it can be done in the same pass
			// over the input
			if (PROGRAM_RECORD_ID != null) {
				// Gather all PG IDs seen in merged input files in first pass.  These are gathered for two reasons:
				// - to know how many different PG records to create to represent this program invocation.
				// - to know what PG IDs are already used to avoid collisions when creating new ones.
				// Note that if there are one or more records that do not have a PG tag, then a null value
				// will be stored in this set.
				pgIdsSeen.add(rec.getStringAttribute(SAMTag.PG.name()));
			}

			if (rec.getReadUnmappedFlag()) {
				if (rec.getReferenceIndex() == -1) {
					// When we hit the unmapped reads with no coordinate, no reason to continue.
					break;
				}
				// If this read is unmapped but sorted with the mapped reads, just skip it.
			} else if (!rec.isSecondaryOrSupplementary()) {
				final ReadEndsForMarkDuplicatesWithMolecularCode fragmentEnd = buildReadEnds(header, index, rec);
				this.fragSort.add(fragmentEnd);

				if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
					final String key = rec.getAttribute(ReservedTagConstants.READ_GROUP_ID) + ":" + rec.getReadName();
					ReadEndsForMarkDuplicates pairedEnds = tmp.remove(rec.getReferenceIndex(), key);

					// See if we've already seen the first end or not
					if (pairedEnds == null) {
						pairedEnds = buildReadEnds(header, index, rec);
						tmp.put(pairedEnds.read2ReferenceIndex, key, pairedEnds);
					} else {
						final int sequence = fragmentEnd.read1ReferenceIndex;
						final int coordinate = fragmentEnd.read1Coordinate;

						// Set orientationForOpticalDuplicates, which always goes by the first then the second end for the strands.  NB: must do this
						// before updating the orientation later.
						if (rec.getFirstOfPairFlag()) {
							pairedEnds.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(rec.getReadNegativeStrandFlag(), pairedEnds.orientation == ReadEnds.R);
						} else {
							pairedEnds.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(pairedEnds.orientation == ReadEnds.R, rec.getReadNegativeStrandFlag());
						}

						// If the second read is actually later, just add the second read data, else flip the reads
						if (sequence > pairedEnds.read1ReferenceIndex ||
								(sequence == pairedEnds.read1ReferenceIndex && coordinate >= pairedEnds.read1Coordinate)) {
							pairedEnds.read2ReferenceIndex = sequence;
							pairedEnds.read2Coordinate = coordinate;
							pairedEnds.read2IndexInFile = index;
							pairedEnds.orientation = ReadEnds.getOrientationByte(pairedEnds.orientation == ReadEnds.R,
									rec.getReadNegativeStrandFlag());
						} else {
							pairedEnds.read2ReferenceIndex = pairedEnds.read1ReferenceIndex;
							pairedEnds.read2Coordinate = pairedEnds.read1Coordinate;
							pairedEnds.read2IndexInFile = pairedEnds.read1IndexInFile;
							pairedEnds.read1ReferenceIndex = sequence;
							pairedEnds.read1Coordinate = coordinate;
							pairedEnds.read1IndexInFile = index;
							pairedEnds.orientation = ReadEnds.getOrientationByte(rec.getReadNegativeStrandFlag(),
									pairedEnds.orientation == ReadEnds.R);
						}

						pairedEnds.score += DuplicateScoringStrategy.computeDuplicateScore(rec, this.DUPLICATE_SCORING_STRATEGY);
						this.pairSort.add((ReadEndsForMarkDuplicatesWithMolecularCode)pairedEnds);
					}
				}
			}

			// Print out some stats every 1m reads
			++index;
			if (progress.record(rec)) {
				log.info("Tracking " + tmp.size() + " as yet unmatched pairs. " + tmp.sizeInRam() + " records in RAM.");
			}
		}

		log.info("Read " + index + " records. " + tmp.size() + " pairs never matched.");
		iterator.close();

		// Tell these collections to free up memory if possible.
		this.pairSort.doneAdding();
		this.fragSort.doneAdding();
	}

	/** Builds a read ends object that represents a single read. */
	private ReadEndsForMarkDuplicatesWithMolecularCode buildReadEnds(final SAMFileHeader header, final long index, final SAMRecord rec) {
		final ReadEndsForMarkDuplicatesWithMolecularCode ends = new ReadEndsForMarkDuplicatesWithMolecularCode();
		ends.read1ReferenceIndex = rec.getReferenceIndex();
		ends.read1Coordinate = rec.getReadNegativeStrandFlag() ? rec.getUnclippedEnd() : rec.getUnclippedStart();
		ends.orientation = rec.getReadNegativeStrandFlag() ? ReadEnds.R : ReadEnds.F;
		ends.read1IndexInFile = index;
		ends.score = DuplicateScoringStrategy.computeDuplicateScore(rec, this.DUPLICATE_SCORING_STRATEGY);

		// Doing this lets the ends object know that it's part of a pair
		if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
			ends.read2ReferenceIndex = rec.getMateReferenceIndex();
		}

		// Fill in the library ID
		ends.libraryId = libraryIdGenerator.getLibraryId(rec);

		// Fill in the location information for optical duplicates
		if (this.opticalDuplicateFinder.addLocationInformation(rec.getReadName(), ends)) {
			// calculate the RG number (nth in list)
			ends.readGroup = 0;
			final String rg = (String) rec.getAttribute("RG");
			final List<SAMReadGroupRecord> readGroups = header.getReadGroups();

			if (rg != null && readGroups != null) {
				for (final SAMReadGroupRecord readGroup : readGroups) {
					if (readGroup.getReadGroupId().equals(rg)) break;
					else ends.readGroup++;
				}
			}
		}
		
		/*
		 * Molecular barcode support. 
		 * Fill in the barcode from header 
		 */
		this.molecularBarcodeFinder.addMolecularBarcodes(rec.getReadName(), ends);

		return ends;
	}

	/**
	 * Goes through the accumulated ReadEndsForMarkDuplicatesWithMolecularCode objects and determines which of them are
	 * to be marked as duplicates.
	 *
	 * @return an array with an ordered list of indexes into the source file
	 */
	private void generateDuplicateIndexes() {
		// Keep this number from getting too large even if there is a huge heap.
		final int maxInMemory = (int) Math.min((Runtime.getRuntime().maxMemory() * 0.25) / SortingLongCollection.SIZEOF,
				(double) (Integer.MAX_VALUE - 5));
		log.info("Will retain up to " + maxInMemory + " duplicate indices before spilling to disk.");
		this.duplicateIndexes = new SortingLongCollection(maxInMemory, TMP_DIR.toArray(new File[TMP_DIR.size()]));

		ReadEndsForMarkDuplicatesWithMolecularCode firstOfNextChunk = null;
		final List<ReadEndsForMarkDuplicatesWithMolecularCode> nextChunk = new ArrayList<ReadEndsForMarkDuplicatesWithMolecularCode>(200);

		// First just do the pairs
		log.info("Traversing read pair information and detecting duplicates.");
		for (final ReadEndsForMarkDuplicatesWithMolecularCode next : this.pairSort) {
			if (firstOfNextChunk == null) {
				firstOfNextChunk = next;
				nextChunk.add(firstOfNextChunk);
			} else if (areComparableForDuplicates(firstOfNextChunk, next, true)) {
				nextChunk.add(next);
			} else {
				if (nextChunk.size() > 1) {
					markDuplicatePairs(nextChunk);
				}

				nextChunk.clear();
				nextChunk.add(next);
				firstOfNextChunk = next;
			}
		}
		if (nextChunk.size() > 1) markDuplicatePairs(nextChunk);
		this.pairSort.cleanup();
		this.pairSort = null;

		// Now deal with the fragments
		log.info("Traversing fragment information and detecting duplicates.");
		boolean containsPairs = false;
		boolean containsFrags = false;

		for (final ReadEndsForMarkDuplicatesWithMolecularCode next : this.fragSort) {
			if (firstOfNextChunk != null && areComparableForDuplicates(firstOfNextChunk, next, false)) {
				nextChunk.add(next);
				containsPairs = containsPairs || next.isPaired();
				containsFrags = containsFrags || !next.isPaired();
			} else {
				if (nextChunk.size() > 1 && containsFrags) {
					markDuplicateFragments(nextChunk, containsPairs);
				}

				nextChunk.clear();
				nextChunk.add(next);
				firstOfNextChunk = next;
				containsPairs = next.isPaired();
				containsFrags = !next.isPaired();
			}
		}
		markDuplicateFragments(nextChunk, containsPairs);
		this.fragSort.cleanup();
		this.fragSort = null;

		log.info("Sorting list of duplicate records.");
		this.duplicateIndexes.doneAddingStartIteration();
	}

	private boolean areComparableForDuplicates(final ReadEndsForMarkDuplicatesWithMolecularCode lhs, final ReadEndsForMarkDuplicatesWithMolecularCode rhs, final boolean compareRead2) {
		boolean retval = lhs.libraryId == rhs.libraryId &&
				lhs.read1ReferenceIndex == rhs.read1ReferenceIndex &&
				lhs.read1Coordinate == rhs.read1Coordinate &&
				lhs.orientation == rhs.orientation;

		if (retval && compareRead2) {
			retval = lhs.read2ReferenceIndex == rhs.read2ReferenceIndex &&
					lhs.read2Coordinate == rhs.read2Coordinate;
		}

		return retval;
	}

	private void addIndexAsDuplicate(final long bamIndex) {
		this.duplicateIndexes.add(bamIndex);
		++this.numDuplicateIndices;
	}

	/**
	 * Takes a list of ReadEndsForMarkDuplicatesWithMolecularCode objects and removes from it all objects that should
	 * not be marked as duplicates.  This assumes that the list contains objects representing pairs.
	 *
	 * @param duplicates list of reads seen as duplicates before molecular code inspection
	 */
	private void markDuplicatePairs(final List<ReadEndsForMarkDuplicatesWithMolecularCode> duplicates) {
		/*
		 * Addition for molecular code support (Charles Girardot)
		 * we split the list based on the molecular codes
		 */
		
		
		Map<String, List<ReadEndsForMarkDuplicatesWithMolecularCode>> groupsOfDupes = null;
		if(USE_UMIS){
			groupsOfDupes = splitDuplicatesByMolecularCodeGroup(duplicates, MISMATCHES, MAX_NUMBER_OF_N);
		}else{
			//we put all of the reads into the DUPES_GROUP_NAME_UNDEF: since this is the only group, one read will be kept as normal (see text below)
			groupsOfDupes = new HashMap<String, List<ReadEndsForMarkDuplicatesWithMolecularCode>>();
			groupsOfDupes.put(DUPES_GROUP_NAME_UNDEF, duplicates);
		}
		/*
		 *  At this point the map can contain an UNDEF category. We do NOT want to let this UNDEF category become a independent group of dupes
		 *  Indeed, imagine the following situation. After positional inspection, we have many (100s of 1000s) small groups with say 3 REAL dupes. 
		 *  But imagine that we systematically have 1 of these 3 dupes with a crapy molecular code (ie ends up UNDEF) and the other 2 reads with the 
		 *  same barcodes eg 'AAAA' => we therefore have at this point 2 groups, one for AAAA and one UNDEF
		 *  If we continue like this we will end up defining 2 individual reads instead of 1 for these 100s of 1000s groups. This can artificially 
		 *  inflate the number f non-dupes reads A LOT. 
		 *  For this reason I prefer to be conservative and considered ALL these degenerate reads as duplicates. The only exception would be if the 
		 *  map only contain one group of the category UNDEF  
		 */
	
		boolean undef_is_only_group = false;
		if(groupsOfDupes.size() == 1 && groupsOfDupes.containsKey(DUPES_GROUP_NAME_UNDEF)){
			undef_is_only_group = true;
		}
		for(Entry<String, List<ReadEndsForMarkDuplicatesWithMolecularCode>> e: groupsOfDupes.entrySet()){
			List<ReadEndsForMarkDuplicatesWithMolecularCode> list = e.getValue();
			boolean isUNDEFGroup = e.getKey().equals(DUPES_GROUP_NAME_UNDEF);
			
			short maxScore = 0;
			ReadEndsForMarkDuplicatesWithMolecularCode best = null;

			/** All read ends should have orientation FF, FR, RF, or RR **/
			/** we only look up the best one if this is not the UNDEF group OR UDEF is the only returned group**/
			if(!isUNDEFGroup || undef_is_only_group){
				for (final ReadEndsForMarkDuplicatesWithMolecularCode end : list) {
					if (end.score > maxScore || best == null) {
						maxScore = end.score;
						best = end;
					}
				}
			}
			
			//flag all reads duplicates but best
			for (final ReadEndsForMarkDuplicatesWithMolecularCode end : list) {
				if (best == null || end != best) {
					addIndexAsDuplicate(end.read1IndexInFile);
					addIndexAsDuplicate(end.read2IndexInFile);
				}
			}

			if (this.READ_NAME_REGEX != null) {
				AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(list, opticalDuplicateFinder, libraryIdGenerator);
			}
		}
	}

	/**
	 * Takes a list of ReadEndsForMarkDuplicatesWithMolecularCode objects and removes from it all objects that should
	 * not be marked as duplicates.  This will set the duplicate index for only list items are fragments.
	 *
	 * @param list
	 * @param containsPairs true if the list also contains objects containing pairs, false otherwise.
	 */
	protected void markDuplicateFragments(final List<ReadEndsForMarkDuplicatesWithMolecularCode> duplicates, final boolean containsPairs) {
		if (containsPairs) {
			for (final ReadEndsForMarkDuplicatesWithMolecularCode end : duplicates) {
				if (!end.isPaired()) addIndexAsDuplicate(end.read1IndexInFile);
			}
		} else {
			
			/*
			 * Addition for molecular code support (Charles Girardot)
			 * we split the list based on the molecular codes
			 */
			
			Map<String, List<ReadEndsForMarkDuplicatesWithMolecularCode>> groupsOfDupes = splitDuplicatesByMolecularCodeGroup(duplicates, MISMATCHES, MAX_NUMBER_OF_N);
			/*
			 *  At this point the map can contain an UNDEF category. We do NOT want to let this UNDEF category become a independent group of dupes
			 *  Indeed, imagine the following situation. After positional inspection, we have many (100s of 1000s) small groups with say 3 REAL dupes. 
			 *  But imagine that we systematically have 1 of these 3 dupes with a crapy molecular code (ie ends up UNDEF) and the other 2 reads with the 
			 *  same barcodes eg 'AAAA' => we therefore have at this point 2 groups, one for AAAA and one UNDEF
			 *  If we continue like this we will end up defining 2 individual reads instead of 1 for these 100s of 1000s groups. This can artificially 
			 *  inflate the number f non-dupes reads A LOT. 
			 *  For this reason I prefer to be conservative and considered ALL these degenerate reads as duplicates. The only exception would be if the 
			 *  map only contain one group of the category UNDEF  
			 */
		
			boolean undef_is_only_group = false;
			if(groupsOfDupes.size() == 1 && groupsOfDupes.containsKey(DUPES_GROUP_NAME_UNDEF)){
				undef_is_only_group = true;
			}
			
			
			for(Entry<String, List<ReadEndsForMarkDuplicatesWithMolecularCode>> e: groupsOfDupes.entrySet()){
				List<ReadEndsForMarkDuplicatesWithMolecularCode> list = e.getValue();
				boolean isUNDEFGroup = e.getKey().equals(DUPES_GROUP_NAME_UNDEF);
				
				short maxScore = 0;
				ReadEndsForMarkDuplicatesWithMolecularCode best = null;

				/** All read ends should have orientation FF, FR, RF, or RR **/
				/** we only look up the best one if this is not the UNDEF group OR UDEF is the only returned group**/
				if(!isUNDEFGroup || undef_is_only_group){
					for (final ReadEndsForMarkDuplicatesWithMolecularCode end : list) {
						if (end.score > maxScore || best == null) {
							maxScore = end.score;
							best = end;
						}
					}
				}
				
				//flag all reads duplicates but best
				for (final ReadEndsForMarkDuplicatesWithMolecularCode end : list) {
					if (best == null || end != best) {
						addIndexAsDuplicate(end.read1IndexInFile);
					}
				}
			
			}
		}
	}

	/**
	 * Split a list of duplicates computed based on their position into sub-lists based on molecular barcodes proximity
	 * @param duplicates
	 * @return the lists of duplicates with their key molecular code. Note the special {@link MarkDuplicatesWithMolecularCode#DUPES_GROUP_NAME_UNDEF} 
	 * key used for reads with too degenerate molecualr code (too many Ns)  
	 * @author girardot@embl.de
	 */
	private <T extends MolecularCodedReadEnds> Map<String, List<T>> splitDuplicatesByMolecularCodeGroup(
			List<T> duplicates, int mismatches, int maxNNumber) {

		Map<String, List<T>> uniqcode2dupes = null;
		if(predefinedBarcodes!=null && predefinedBarcodes.size()>0){
			uniqcode2dupes = splitDuplicatesByMolecularCodeGroupWithPredefinedCodeList(duplicates, mismatches, maxNNumber, barcodeLength, predefinedBarcodes);
		}else{
			uniqcode2dupes = splitDuplicatesByMolecularCodeGroupWithoutPredefinedCodeList(duplicates, mismatches, maxNNumber);
		}
		
		
		return uniqcode2dupes;
	}	

	
	/**
	 * Inspects a list of {@link MolecularCodedReadEnds} and distribute them into subgroups based on their molecular codes.
	 * Each created group can then be regarded as a real list of duplicates. 
	 * This method tries its best to discover the real initial barcodes. 
	 * If you know the exact list of expected molecualr codes, you should instead use 
	 * {@link MarkDuplicatesWithMolecularCode#splitDuplicatesByMolecularCodeGroupWithPredefinedCodeList(List, int, int, int, Set)}
	 *   
	 * 
	 * @param duplicates the list of {@link MolecularCodedReadEnds} to split into subgroups based on molecular codes
	 * @param mismatches the number of mismatches (inclusive) to still consider 2 sequences identical. 
	 * @param maxNNumber the max number of N a code can have. Passed this number of N in code, the {@link MolecularCodedReadEnds} 
	 * are collected into a unique 'UNDEF' group 
	 * @return the lists of duplicates with their key molecular code
	 * @author girardot@embl.de 
	 */
	protected <T extends MolecularCodedReadEnds> Map<String, List<T>> splitDuplicatesByMolecularCodeGroupWithoutPredefinedCodeList(
				List<T> duplicates, int mismatches, int maxNNumber) {
			
		//build a map with all codes
		HashMap<String, List<T>> code2dupes = new HashMap<String, List<T>>();
		HashMap<String, Integer> code2numberOfNs = new HashMap<String, Integer>(); // number of N in the code sequence
 		for (T end : duplicates) {
 			String code = end.getMolecularCode();
			if(!code2dupes.containsKey(code)){
				code2dupes.put(code, new ArrayList<T>());
				code2numberOfNs.put(code, StringUtils.countMatches(code, "N"));
			}
			code2dupes.get(code).add(end);
		}
 		//counts
 		Map<String, Integer> counts = new HashMap<String, Integer>();
 		for (Entry<String, List<T>> e : code2dupes.entrySet()) {
 			counts.put(e.getKey(), e.getValue().size());
		}
 		
 		List<String> keys = new ArrayList<String>(code2dupes.keySet());
 		/*
 		 * sort the codes : this gives a list with sequence ordered in (1)decreasing order of N they contain and, for equivalent number of N, (2) the sequence occurring the most
 		 * ie the first seq is the most occurring sequence with the lowest number of N 
 		 */
 		Collections.sort(keys, new CompareSequenceByNNumAndOccurence(counts));
		
 		/*
		 * follow the sorted list and add as a new entry in uniqcode2dupes if the code does not match any of the already added code
		 * use a linked hashmap to keep traversing the map in an ordered way 
		 */
 		Map<String, List<T>> uniqcode2dupes = Collections.synchronizedMap( new LinkedHashMap<String, List<T>>() );
 		Map<String, List<String>> uniqcodeAliases = new HashMap<String, List<String>>();
 		for (String code : keys) {
			//compute 
 			if(uniqcode2dupes.isEmpty()) {
 				/*
 				 * This is the first code we inspect ie it contains the least observed number of N (due to sorting)
 				 * A special case is when this first code is already too degenerate ie contains more Ns than allowed
 				 * In this situation, no need to continue processing => all reads will end up in UNDEF category and we can therefore directly 
 				 * return 
 				 */
 				if(code2numberOfNs.get(code) > maxNNumber){
 					return Collections.singletonMap(DUPES_GROUP_NAME_UNDEF, duplicates);
 				}else{
 					uniqcode2dupes.put(code, code2dupes.get(code));
 					List<String> aliases = new ArrayList<String>();
 					aliases.add(code);
 					uniqcodeAliases.put(code, aliases);
 					continue;
 				}
 			}
 			
 			if(code2numberOfNs.get(code) > maxNNumber){
 				//goes in the UNDEF group
 				if(!uniqcode2dupes.containsKey(DUPES_GROUP_NAME_UNDEF))
 					uniqcode2dupes.put(DUPES_GROUP_NAME_UNDEF, code2dupes.get(code));
 				else
 					uniqcode2dupes.get(DUPES_GROUP_NAME_UNDEF).addAll(code2dupes.get(code));
 				
 				continue;
			}
 			//check with each of the code already added in uniqcode2dupes, stop checking as soon as a 'merge case is found' 
 			String mergeWith = null;
 			LOOP : for (String uniqCode : uniqcode2dupes.keySet()) {
 				//check each aliases
 				for(String k : uniqcodeAliases.get(uniqCode)){
 					if(countMismatches(code, k) <= mismatches){
 						//this is a merge situation
 						mergeWith = uniqCode;
 						break LOOP;
 					}
 				}
			}
 			if(mergeWith != null){
 				uniqcode2dupes.get(mergeWith).addAll(code2dupes.get(code));
 				/*
 				 * only use sequence without Ns as aliases (because N are not considered as a mismatch).
 				 * We will indeed end up with very bad degenerate sequences that will match any subsequences if we dont do this  
 				 */
 				
 				if(code2numberOfNs.get(code) == 0)
 					uniqcodeAliases.get(mergeWith).add(code);
 			}else{
 				uniqcode2dupes.put(code, code2dupes.get(code));
 				List<String> aliases = new ArrayList<String>();
 				aliases.add(code);
 				uniqcodeAliases.put(code, aliases);
 			}
		}
 		
 		
		return uniqcode2dupes;
	}
	
	/**
	 * Inspects a list of {@link MolecularCodedReadEnds} and distribute them into subgroups based on the proximity of their 
	 * molecular codes to a predefined list of expected codes.
	 * If you don't know the exact list of expected molecular codes, you should instead use 
	 * {@link MarkDuplicatesWithMolecularCode#splitDuplicatesByMolecularCodeGroupWithoutPredefinedCodeList(List, int, int)}
	 *   
	 * 
	 * @param duplicates the list of {@link MolecularCodedReadEnds} to split into subgroups based on molecular codes
	 * @param mismatches the number of mismatches (inclusive) to still consider 2 sequences identical. 
	 * @param maxNNumber the max number of N a code can have. Passed this number of N in code, the {@link MolecularCodedReadEnds} 
	 * are collected into a unique 'UNDEF' group 
	 * @param barcodeLength the barcode length. This is then used to call {@link MolecularCodedReadEnds#getMolecularCodes(int)} and
	 *  each returned code will be used to look up the predefined list independently. Finally, all identified predefined codes are sorted 
	 *  and concatenated into a unique code. This unique code is then considered to form a duplicate group.
	 *  @param predefinedBarcodes the list of predefined codes used to tag the reads 'molecularly'   
	 * @return the lists of duplicates with their key molecular code
	 * @author girardot@embl.de 
	 */
	protected <T extends MolecularCodedReadEnds> Map<String, List<T>> splitDuplicatesByMolecularCodeGroupWithPredefinedCodeList(
			List<T> duplicates, int mismatches, int maxNNumber, int barcodeLength, Set<String> predefinedBarcodes) {
		
		Map<String, List<T>> uniqcode2dupes = new HashMap<String, List<T>>();
		
		//for each read, find the matching predefined code for each molecular code (as many as SLOTS size)
		for (T end : duplicates) {
			List<String> codes = end.getMolecularCodes(barcodeLength);
			List<String> identifiedcodes = new ArrayList<String>();
			boolean goesToUndefGroup = false;
			for (String code : codes) {
				String best = null;
				//is the code too degenerate to be considered ?
				if( StringUtils.countMatches(code, "N") > maxNNumber ){
					goesToUndefGroup = true;
				}else{
					for(String k : predefinedBarcodes){
						if(countMismatches(code, k) <= mismatches){
							best = k;
							break; // we break as soon as we found a match in MISMATCHES range
						}
					}
					if(best == null){
						//this is far apart from expected, all reads in this situation are placed in the same same group to be seen as dups
						goesToUndefGroup = true;
					}else{
						identifiedcodes.add( best );
					}
				}
			}
			
			/*
			 * turn list into a unique code, sort first (?) then concat. unless goes to UNDEF group
			 * The sorting is not obvious : if I sort then A:B is same as B:A (from read header) which supposes that 
			 * order is not important. Since UMIs are placed in read headers in a defined order (i.e. UMI from read 1 then from read 2)
			 *  and there is no particular reason to consider that UMI in read 1 is the same as in read 2, 
			 *  we should preserve the order and consider A:B is not B:A => so don't sort    
			 */
			String uniqCode = DUPES_GROUP_NAME_UNDEF;
			if(!goesToUndefGroup){
				// Collections.sort(identifiedcodes); //don t sort ie see above
				uniqCode = org.embl.cg.utilitytools.utils.StringUtil.mergeList(identifiedcodes, ":");
			}
			
			//register
			if(!uniqcode2dupes.containsKey(uniqCode)){
				uniqcode2dupes.put(uniqCode, new ArrayList<T>());
			}
			uniqcode2dupes.get(uniqCode).add(end);
		}
		
		
		return uniqcode2dupes;
		
	}

	/**
	 * 
	 * Compare barcode sequence to bases from read
	 * @param barcode1 the barcode 1 
	 * @param barcode2 the barcode 2 
	 * @return how many bases are different, 'N' always considered as a similar base 
	 * 
	 */
	protected int countMismatches(final String barcode1, final String barcode2 ) {
		 byte[] barcodeBytes1 = StringUtil.stringToBytes(barcode1);
		 byte[] barcodeBytes2 = StringUtil.stringToBytes(barcode2);
		int numMismatches = 0;
		
		for (int i = 0; i < barcodeBytes1.length; ++i) {
			if (!SequenceUtil.isNoCall(barcodeBytes1[i]) && !SequenceUtil.isNoCall(barcodeBytes2[i]) ) { //returns true if the value of base represents a no call
				//if bases in barcode1 and barcode2 are different => increase mismatch count
				if (!SequenceUtil.basesEqual(barcodeBytes1[i], barcodeBytes2[i])){
					++numMismatches;
				}
			}
		}

		return numMismatches;
	}

	
	
	
	/** Comparator for ReadEndsForMarkDuplicatesWithMolecularCode that orders by read1 position then pair orientation then read2 position. */
	static class ReadEndsMDComparator implements Comparator<ReadEndsForMarkDuplicatesWithMolecularCode> {
		public int compare(final ReadEndsForMarkDuplicatesWithMolecularCode lhs, final ReadEndsForMarkDuplicatesWithMolecularCode rhs) {
			int retval = lhs.libraryId - rhs.libraryId;
			if (retval == 0) retval = lhs.read1ReferenceIndex - rhs.read1ReferenceIndex;
			if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
			if (retval == 0) retval = lhs.orientation - rhs.orientation;
			if (retval == 0) retval = lhs.read2ReferenceIndex - rhs.read2ReferenceIndex;
			if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
			if (retval == 0) retval = (int) (lhs.read1IndexInFile - rhs.read1IndexInFile);
			if (retval == 0) retval = (int) (lhs.read2IndexInFile - rhs.read2IndexInFile);
			return retval;
		}
	}
}

