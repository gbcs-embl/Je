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
package org.embl.gbcs.je.retag;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.embl.cg.utilitytools.utils.parser.csv.CSVLine;
import org.embl.cg.utilitytools.utils.parser.csv.CSVParser;
import org.embl.cg.utilitytools.utils.parser.csv.InvalidCSVSetUpException;
import org.embl.gbcs.je.JeUtils;
import org.embl.gbcs.je.Jexception;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

/**
 * 
 * @author Charles Girardot 
 */
@CommandLineProgramProperties(
		usage = "Extracts barcode and UMI sequence(s) embedded in read names and tag reads with proper BAM tag.", 
		usageShort = "je tag INPUT=file.bam OUTPUT=result.bam ", 
		programGroup = SamOrBam.class
		)
public class TagFromReadName extends CommandLineProgram {
	/*
	 * Relevant SAM tag to use :
	 * - BC for sample barcode, raw or corrected, with QT to store its quality string
	 * - RX : Sequence bases of the (possibly corrected) unique molecular identifier
	 * - QX : Quality score of the unique molecular identifier in the RX tag
	 * - OX : Original unique molecular barcode bases
	 * - BZ : Phred quality of the unique molecular barcode bases in the OX tag
	 * - MI : Molecular identifier; a string that uniquely identifies the molecule from which the record was derived
	 */
	
	public static final String BARCODE_UMI_SEPARATOR = "-";


	private final Log log = Log.getInstance(TagFromReadName.class);


	// input fastq files
	@Option(shortName="I", 
			optional = false,
			printOrder=10,
			doc="Input SAM/BAM file"
			)
	public File INPUT;
	
	@Option(shortName="O", 
			optional = false,
			printOrder=20,
			doc="Output SAM/BAM file"
			)
	public File OUTPUT;
	
	@Option(optional = true,
			printOrder = 30,
			doc="Where to find the BARCODE(s) in the read name once read name has been tokenized using the SPLIT character (e.g. ':'). \n"
					+ "This option can be specified multiple time when multiple BARCODES are present in read header.\n"
					+ "Counting starts at 1 and negative numbers can be used to start counting from the end (last token is '-1').\n"
					+ "BARCODE(s) extracted from read name are used to assemble a 'BC:Z:GATCCTAG' tag (BC is default, see BC_TAG). \n"
					+ "Following SAM specifications, in the case of multiple barcodes, all the barcodes are concatenated using \n"
					+ "a hyphen ('-') between the barcodes. The order of concatenation follows the order of BC_SLOT in the command line.\n"
					+ "For example, consider the following read name that lists 3 different barcodes in the end : \n"
					+ "\t HISEQ:44:C6KC0ANXX:8:2112:20670:79594:CGATGTTT:GATCCTAG:AAGGTACG \n"
					+ "\t to indicate that the three slots contain barcodes, use \n"
					+ "\t\t BC_SLOT=-1 BC_SLOT=-2 BC_SLOT=-3 ; which will result as BC:Z:AAGGTACG-GATCCTAG-CGATGTTT\n"
					+ "\t if only the 2 last ones should be considered, use \n"
					+ "\t\t BC_SLOT=-1 BC_SLOT=-2 ; which will result as BC:Z:AAGGTACG-GATCCTAG"
					+ "\t Note that BC_SLOT order matters as : \n"
					+ "\t\t BC_SLOT=-2 BC_SLOT=-1 ; would result as BC:Z:GATCCTAG-AAGGTACG"
					)
	public List<Integer> BC_SLOT = null;
	
	@Option(optional = true,
			printOrder = 40,
			doc="Where to find the UMI(s) in the read name once read name has been tokenized using the SPLIT character (e.g. ':'). \n"
					+ "This option can be specified multiple time when multiple UMIs are present in read header.\n"
					+ "Counting starts at 1 and negative numbers can be used to start counting from the end (last token is '-1').\n"
					+ "UMI(s) extracted from read name are used to assemble both a RX and OX tag e.g. 'RX:Z:GATCCTAG' tag. \n"
					+ "Following SAM specifications, in the case of multiple UMIs, all the UMIs are concatenated using \n"
					+ "a hyphen ('-'). The order of concatenation follows the order of UMI_SLOT in the command line.\n"
					+ "For example, consider the following read name that lists 3 different sequence codes in the end : \n"
					+ "\t HISEQ:44:C6KC0ANXX:8:2112:20670:79594:CGATGTTT:GATCCTAG:AAGGTACG \n"
					+ "\t to indicate that the 2 last slots contain UMIs, use \n"
					+ "\t\t UMI_SLOT=-1 UMI_SLOT=-2 ; which will result as RX:Z:AAGGTACG-GATCCTAG"
					+ "\t Note that UMI_SLOT order matters as : \n"
					+ "\t\t UMI_SLOT=-2 UMI_SLOT=-1 ; would result as RX:Z:GATCCTAG-AAGGTACG"
					)
	public List<Integer> UMI_SLOT = null;
	
	
	@Option(shortName = "CC", optional = true,
			printOrder = 50,
			doc="Group of barcodes (clusters) for on-the-fly correction of the BARCODEs. This file contains mapping between the barcodes sequences\n"
					+ "as originally found in the FASTQ reads (and in the read name of the input BAM) and the corrected sample barcode\n"
					+ "sequences (e.g. as identified by clustering with mismatches using tools like starcode or vsearch ) \n"
					+ "that should be written in the BC tag in place of the original sequence.\n"
					+ "If provided, the BARCODE sequence(s) extracted from read name are converted to the 'real' barcodes according\n"
					+ "to the mapping described in this file. If the sequence is not found in the supplied mapping file, the read is either\n"
					+ "trashed or kept (according to option KEEP_UNASSIGNED_BARCODES), in which case the value defined by UNASSIGNED_BARCODE_VALUE\n"
					+ "is used.\n"
					+ "Format: two column text file, one cluster per line with the real barcode in the first line and the comma separated\n"
					+ "list of codes in the second column i.e. :\n"
					+ "\t\t ACTGTAC \tACTCTAC,TCTGTAC,ACTGTAG  \n"
					+ "All the codes MUST have the same length"
			)
	public File CLUSTERED_CODES_FILE = null;
	
	@Option(shortName="KUP", optional = true,
			printOrder = 60,
			doc="Should read be keep when no mapping was defined for the orginal barcode sequence in provided CLUSTERED_CODES_FILE.\n"
					+ "If false, the read is not written in output file." 
			)
	public boolean KEEP_UNASSIGNED_BARCODES = true;
	
	@Option(shortName = "UBV", optional = true,
			printOrder = 70,
			doc="Value to use for the BARCODE tag when CLUSTERED_CODES_FILE was provided and no mapping was defined for the orginal barcode sequence." 
					)
	public String UNASSIGNED_BARCODE_VALUE = "NA";
	
	@Option(shortName="T", optional = true,
			printOrder = 80,
			doc="Should barcode/UMIs information be removed from read names in the output BAM ? " 
			)
	public boolean TRIM_HEADERS = false;

	@Option(shortName = "TSLOTS", optional = true,
			printOrder = 90,
			doc="Where to find *all* barcode(s) and UMIs in the read name once has been tokenized using the SPLIT character (e.g. ':'). \n"
					+ "This option is only considered when TRIM_HEADERS=true and should only be used when UMI_SLOT and BC_SLOT do not\n"
					+ "describe all the slots that should be trimmed. When TSLOTS is ommited while TRIM_HEADERS=true, the values\n"
					+ "of UMI_SLOT and BC_SLOT apply.\n"
					+ "IMPORTANT : counting starts at 1 and negative numbers can be used to start counting from the end.\n"
					+ "See UMI_SLOT help for examples." 
			)
	public List<Integer> TSLOTS = null;
	
	
	@Option(shortName = "SPLIT", optional = true,
			printOrder = 130,
			doc="Character to use to split up the read header line, default is ':'." 
					)
	public String SPLIT_CHAR = ":";
	
	@Option(optional = true,
			printOrder = 140,
			doc="SAM Tag to use to store barcode(s) sequences extracted from barcode slots (BC by default)."
					+ " Do not change unless you have good reasons to."
					)
	public String BC_TAG = JeUtils.SAMTAG_BC;
	
	@Option(optional = true,
			printOrder = 141,
			doc="SAM Tag to use to store barcode(s) quality score extracted from barcode slots (QT by default)."
					+ " Do not change unless you have good reasons to."
					)
	public String QT_TAG = JeUtils.SAMTAG_QT;
	
	
	@Option(optional = true,
			printOrder = 144,
			doc="Should the RX (and QX when relevant) SAM Tag(s) be used to store UMI(s)"
					+ " sequence (and quality) extracted from UMI slots. "
					+ "Set to FALSE if you don't want these tags to be set. "
					)
	public boolean WITH_RX = true;
	
	@Option(optional = true,
			printOrder = 144,
			doc="Should the OX (and BZ when relevant) SAM Tag(s) be used to store UMI(s)"
					+ " sequence (and quality) extracted from UMI slots. "
					+ "Set to FALSE if you don't want these tags to be set. "
					)
	public boolean WITH_OX = true;
	
	
	@Option(optional = true,
			printOrder = 146,
			doc="SAM Tag to use to store original UMI(s) sequence extracted from UMI slots (instead of RX / OX)"
	)
	public String UMI_SEQ_TAG = null;
	
	@Option(optional = true,
			printOrder = 148,
			doc="SAM Tag to use to store original UMI(s) sequence extracted from UMI slots (instead of QX / BZ)"
					)
	public String UMI_QUAL_TAG = null;
	
	
	
	@Option(shortName="ARG", optional = true,
			printOrder = 150,
			doc="Should a read group be created for each barcode. This option is only considered when providing a CLUSTERED_CODES_FILE." 
			)
	public boolean ADD_RG = false;

	
	@Option(doc = "Read Group platform (e.g. illumina, solid) ; only considered when RG=true")
    public String RGPL;
	
	@Option(doc = "Read Group program group; only considered when RG=true", optional = true)
    public String RGPG;
	
	
	@Option(shortName = StandardOptionDefinitions.PROGRAM_RECORD_ID_SHORT_NAME,
            doc = "The program record ID for the @PG record(s) created by this program. Set to null to disable " +
                    "PG record creation.  This string may have a suffix appended to avoid collision with other " +
                    "program record IDs.",
            optional = true)
    public String PROGRAM_RECORD_ID = "TagFromReadName";

	@Option(shortName = "PG_VERSION",
            doc = "Value of VN tag of PG record to be created. If not specified, the version will be detected automatically.",
            optional = true)
    public String PROGRAM_GROUP_VERSION;

	@Option(shortName = "PG_COMMAND",
            doc = "Value of CL tag of PG record to be created. If not supplied the command line will be detected automatically.",
            optional = true)
    public String PROGRAM_GROUP_COMMAND_LINE;

	@Option(shortName = "PG_NAME",
            doc = "Value of PN tag of PG record to be created.")
    public String PROGRAM_GROUP_NAME = getClass().getSimpleName();

	@Option(shortName = "CO",
            doc = "Comment(s) to include in the output file's header.",
            optional = true)
    public List<String> COMMENT = new ArrayList<String>();

	
	
	
	/** The program groups that have been seen during the course of examining the input records. */
    protected final Set<String> pgIdsSeen = new HashSet<String>();
	
	
	public TagFromReadName() {
		super();
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		new TagFromReadName().instanceMainWithExit(args);
	}

	
	
	/**
	 * Main work method.  Reads the BAM file once and collects sorted information about
	 * the 5' ends of both ends of each read (or just one end in the case of pairs).
	 * Then makes a pass through those determining duplicates before re-reading the
	 * input file and writing it out with duplication flags set correctly.
	 */
	protected int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		
		/*
		 * Check that at least BC_SLOT or UMI_SLOT is not empty
		 */
		if( (BC_SLOT ==null || BC_SLOT.isEmpty()) && (UMI_SLOT ==null || UMI_SLOT.isEmpty()) )
			throw new IllegalArgumentException("At least one of BC_SLOT or UMI_SLOT must be provided.");
		
		
		boolean extracting_barcodes = BC_SLOT !=null && !BC_SLOT.isEmpty();
		boolean extracting_UMIs = UMI_SLOT !=null && !UMI_SLOT.isEmpty();
		
		if(TSLOTS == null || TSLOTS.size() == 0){
			TSLOTS = new ArrayList<Integer>();
			if(extracting_barcodes)
				TSLOTS.addAll(BC_SLOT);
			if(extracting_UMIs)
				TSLOTS.addAll(UMI_SLOT);
			
		}
		
		/*
		 * Do we have a pre-defined list of barcodes ?
		 */
		Map<String, String> sequence2barcode = null;
		Set<String> correctedcodes = null;
		if(CLUSTERED_CODES_FILE!=null){
			IOUtil.assertFileIsReadable(CLUSTERED_CODES_FILE);
			
			//load up barcodes 
			CSVParser p = null;
			sequence2barcode = new HashMap<String, String>();
			correctedcodes = new TreeSet<String>();
			try {
				p = new CSVParser(-1,-1,false);
				Iterator<CSVLine> lines = p.iterator(CLUSTERED_CODES_FILE.getAbsolutePath(), "\t");
				int n = 0;
				int barcodeLength = -1;
				while (lines.hasNext()) {
					CSVLine l = lines.next();
					n++;
					String code = l.getValue(0);
					correctedcodes.add(code);
					sequence2barcode.put(code, code); //make sure code is in there
					if(barcodeLength < 0 ){
						barcodeLength = code.length();
					}else if(code.length() != barcodeLength){
						throw new IllegalArgumentException("Barcode '"+code+"' read from "+CLUSTERED_CODES_FILE.getName()+
								" in col 1 at line "+n+" do not have the same length as previously read codes ("+barcodeLength+" bases)");
					}
					for(String s :l.getValue(1).split(",")){
						if(s.length() != barcodeLength){
							throw new IllegalArgumentException("Barcode '"+s+"' read from "+CLUSTERED_CODES_FILE.getName()+
									" in col 2 at line "+n+" do not have the same length as previously read codes ("+barcodeLength+" bases)");
						}
						sequence2barcode.put(s, code);
					}
				}
			} catch (InvalidCSVSetUpException e) {
				throw new IllegalArgumentException("Cannot read barcode mappings from "+CLUSTERED_CODES_FILE.getName()+": file does not seem to have the expected format.");

			} catch (IOException e) {
				 throw new IllegalArgumentException("Cannot read barcode mappings from "+CLUSTERED_CODES_FILE.getName()+". Check format.");
			} finally{
				if(p!=null)
					p.terminate();
			}
		}		
		boolean correcting_barcodes = sequence2barcode != null;
		
		boolean adding_readgroup = correcting_barcodes && ADD_RG;
		
		reportMemoryStats("Start of doWork");
		SamReaderFactory readerFactory = SamReaderFactory.makeDefault();
        SamReader reader = readerFactory.enable(SamReaderFactory.Option.EAGERLY_DECODE).open(SamInputResource.of(INPUT)) ;
        final SAMFileHeader header = reader.getFileHeader();

		
		final SAMFileHeader outputHeader = header.clone();
		outputHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
		
		
		/*
		 * create all needed read groups
		 */
		Map<String, SAMReadGroupRecord> readgroups = null;
		if(adding_readgroup){
			readgroups = new HashMap<String, SAMReadGroupRecord>();
			int i=0;
			for(String code : correctedcodes){
				i++;
				String RGID = ""+i;
				final SAMReadGroupRecord rg = new SAMReadGroupRecord(RGID);
				rg.setPlatform(RGPL);
				if (RGPG != null) rg.setProgramGroup(RGPG);
				rg.setKeySequence(code);
				readgroups.put(RGID,  rg);
			}
			
			outputHeader.setReadGroups(new ArrayList<SAMReadGroupRecord>(readgroups.values()));
		}

		//add comment if any
		for (final String comment : COMMENT) outputHeader.addComment(comment);

		// Key: previous PG ID on a SAM Record (or null).  Value: New PG ID to replace it.
		final Map<String, String> chainedPgIds = getChainedPgIds(outputHeader);

				
		//open output
		final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(outputHeader,
				true,
				OUTPUT);

		final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Written");
		final CloseableIterator<SAMRecord> iterator = reader.iterator();
		
		
		/*
		 * help for read name trimming (option T) if needed
		 */
		Set<Integer> readNameSlotsToIgnore = new TreeSet<Integer>();
		Boolean barcodeWithQual = null; //not init
		Boolean umiWithQual = null; //not init
		while (iterator.hasNext()) {
			final SAMRecord rec = iterator.next();
			
			String readName = rec.getReadName();
			String[] tokens = readName.split(SPLIT_CHAR);
			
			/*
			 * With this option the user wants to eliminate the barcode slots from the read name.
			 * This is handy to clip back all (molecular) barcodes from the name 
			 */
			if(TRIM_HEADERS){
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
			
			/*
			 * Extract barcodes if required
			 */
			if(extracting_barcodes){
				if(barcodeWithQual == null) {
					//do we have quality numbers with the barcodes ? Lets supposed not and check if we got numbers in the result string 
					barcodeWithQual = extractCodes(tokens, BC_SLOT).getReadString().matches(".*[0-9]+.*");
				}
				//get barcodes
				FastqRecord _barcode = barcodeWithQual ? extractCodesWithQual(tokens, BC_SLOT) : extractCodes(tokens, BC_SLOT);
				String barcode = _barcode.getReadString();
				//optionally correct 
				if(correcting_barcodes){
					barcode = sequence2barcode.get(barcode);
					if(barcode == null){
						if(!KEEP_UNASSIGNED_BARCODES)
							continue;
						
						barcode = UNASSIGNED_BARCODE_VALUE;
					}
					else if(adding_readgroup){
						rec.setAttribute(SAMTag.RG.name(), readgroups.get(barcode));
					}
				}
				//add TAG 
				rec.setAttribute(this.BC_TAG, barcode);
				//and quality
				if(barcodeWithQual)
					rec.setAttribute(this.QT_TAG, _barcode.getBaseQualityString());
					
			}
			
			if(extracting_UMIs){
				if(umiWithQual == null) {
					//do we have quality numbers with the umis ? Lets supposed not and check if we got numbers in the result string 
					umiWithQual = extractCodes(tokens, UMI_SLOT).getReadString().matches(".*[0-9]+.*");
				}
				//get umis
				FastqRecord _umi = umiWithQual ? extractCodesWithQual(tokens, UMI_SLOT) : extractCodes(tokens, UMI_SLOT);
				String umi = _umi.getReadString();
				
				//add TAG 
				if(WITH_RX) {
					rec.setAttribute("RX", umi);
					//and quality
					if(umiWithQual)
						rec.setAttribute("QX", _umi.getBaseQualityString());
						
				}
				if(WITH_OX) {
					rec.setAttribute("OX", umi);
					//and quality
					if(umiWithQual)
						rec.setAttribute("BZ", _umi.getBaseQualityString());
						
				}
				if(UMI_SEQ_TAG!=null){
					rec.setAttribute(UMI_SEQ_TAG, umi);
					//and quality
					if(umiWithQual)
						rec.setAttribute(UMI_QUAL_TAG, _umi.getBaseQualityString());
						
				}
			}
			
			
			if (PROGRAM_RECORD_ID != null) {
				rec.setAttribute(SAMTag.PG.name(), chainedPgIds.get(rec.getStringAttribute(SAMTag.PG.name())));
			}
			out.addAlignment(rec);
			progress.record(rec);
		}

		// remember to close the inputs
		iterator.close();

		reportMemoryStats("Before output close");
		out.close();
		reportMemoryStats("After output close");

		return 0;
	}

	
	
	static Pattern codeWithQuality = Pattern.compile("^([a-Z]+)([0-9]*)$");
	static Matcher codeWithQualityM = codeWithQuality.matcher("");
	
	/**
	 * extract all the barcodes from indicated slots
	 * @param tokens
	 * @param slots
	 * @return
	 */
	private FastqRecord extractCodesWithQual(String[] tokens, List<Integer> slots) {
		StringBuffer molCode = new StringBuffer();
		StringBuffer molQual = new StringBuffer();
		boolean empty = true;
		for (int slot : slots) {
			String umi = null;
			if(slot>=0){
				umi = tokens[slot];
			}else{
				//from end
				umi = tokens[tokens.length + slot];
			}
			//extract sequence and qual scores
			codeWithQualityM.reset(umi);
			if(!codeWithQualityM.matches()){
				throw new Jexception("String extracted from read name slot "+slot+" does not comply with expected 'barcode with quality' format '[a-Z]+[0-9]*' :"+umi);
			}
			String umiSeq = codeWithQualityM.group(1);
			String umiQual = codeWithQualityM.group(2);
			if(umiQual.length() != (2 * umiSeq.length()) ){
				throw new Jexception("Lenght of code quality string must be twice the size of the code but is not the case in "+umi+" (slot "+slot+")");
			}
			if(!empty) {
				molCode.append(BARCODE_UMI_SEPARATOR);
				molQual.append(" ");
			}
			molCode.append(umiSeq);
			molQual.append(  toBytesThenPhred(umiQual) );
			empty = false;
		}
		return new FastqRecord("", molCode.toString(), "", molQual.toString() );
	}
	
	public String toBytesThenPhred(String s) {
		byte [] arr = new byte [s.length()/2];
		int i =0;
		for(String t : s.split("(?<=\\G.{2})")) { // the regex splits string into pairs of char
			arr[i] = Byte.parseByte(t);
			i++;
		}
		return SAMUtils.phredToFastq(arr);
	}

	/**
	 * extract all the barcodes from indicated slots
	 * @param tokens
	 * @param slots
	 * @return
	 */
	private FastqRecord extractCodes(String[] tokens, List<Integer> slots) {
		StringBuffer molCode = new StringBuffer();
		boolean empty = true;
		for (int slot : slots) {
			String umi = null;
			if(slot>=0){
				umi = tokens[slot];
			}else{
				//from end
				umi = tokens[tokens.length + slot];
			}
			
			if(!empty)
				molCode.append(BARCODE_UMI_SEPARATOR);
			
			molCode.append(umi);
			empty = false;
		}
		return new FastqRecord("", molCode.toString(), "","" );
	}

	/**
     * We have to re-chain the program groups based on this algorithm.  This returns the map from existing program group ID
     * to new program group ID.
     */
    protected Map<String, String> getChainedPgIds(final SAMFileHeader outputHeader) {
        final Map<String, String> chainedPgIds;
        // Generate new PG record(s)
        if (PROGRAM_RECORD_ID != null) {
            final SAMFileHeader.PgIdGenerator pgIdGenerator = new SAMFileHeader.PgIdGenerator(outputHeader);
            if (PROGRAM_GROUP_VERSION == null) {
                PROGRAM_GROUP_VERSION = this.getVersion();
            }
            if (PROGRAM_GROUP_COMMAND_LINE == null) {
                PROGRAM_GROUP_COMMAND_LINE = this.getCommandLine();
            }
            chainedPgIds = new HashMap<String, String>();
            for (final String existingId : this.pgIdsSeen) {
                final String newPgId = pgIdGenerator.getNonCollidingId(PROGRAM_RECORD_ID);
                chainedPgIds.put(existingId, newPgId);
                final SAMProgramRecord programRecord = new SAMProgramRecord(newPgId);
                programRecord.setProgramVersion(PROGRAM_GROUP_VERSION);
                programRecord.setCommandLine(PROGRAM_GROUP_COMMAND_LINE);
                programRecord.setProgramName(PROGRAM_GROUP_NAME);
                programRecord.setPreviousProgramGroupId(existingId);
                outputHeader.addProgramRecord(programRecord);
            }
        } else {
            chainedPgIds = null;
        }
        return chainedPgIds;
    }
	
	/** Print out some quick JVM memory stats. */
	private void reportMemoryStats(final String stage) {
		System.gc();
		final Runtime runtime = Runtime.getRuntime();
		log.info(stage + " freeMemory: " + runtime.freeMemory() + "; totalMemory: " + runtime.totalMemory() +
				"; maxMemory: " + runtime.maxMemory());
	}
}

