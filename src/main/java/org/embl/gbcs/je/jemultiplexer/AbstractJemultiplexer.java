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
package org.embl.gbcs.je.jemultiplexer;




import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.FastqQualityFormat;
import htsjdk.samtools.util.QualityEncodingDetector;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SolexaQualityConverter;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.embl.cg.utilitytools.utils.ExceptionUtil;
import org.embl.cg.utilitytools.utils.FileUtil;
import org.embl.cg.utilitytools.utils.StringUtil;
import org.embl.gbcs.embase.api.model.NGSLibrary;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.sam.FastqToSam;

public abstract class AbstractJemultiplexer extends CommandLineProgram {
	private static Logger log = LoggerFactory.getLogger(AbstractJemultiplexer.class);


	/*
	 * constants
	 */
	
	//sample name for unassigned read 
	protected static final String UNASSIGNED = "UNASSIGNED";

	
	/**
	 * Defaults
	 */
	protected static final BarcodePosition DEFAULT_BARCODE_READ_POS = BarcodePosition.BOTH;
	protected static final BarcodePosition DEFAULT_BARCODE_FOR_SAMPLE_MATCHING = BarcodePosition.BOTH;
	protected static final boolean DEFAULT_REDUNDANT_BARCODES = true;
	protected static final boolean DEFAULT_STRICT = false;
	protected static final Integer DEFAULT_MAX_MISMATCHES = 1;
	protected static final Integer DEFAULT_MIN_MISMATCH_DELTA = 1;
	protected static final Integer DEFAULT_MIN_BASE_QUALITY = 10;
	protected static final Integer DEFAULT_XTRIMLEN = 0;
	protected static final Integer DEFAULT_ZTRIMLEN = 0;
	protected static final boolean DEFAULT_CLIP_BARCODE = true;
	protected static final boolean DEFAULT_ADD_BARCODE_TO_HEADER = true;
	protected static final FastqQualityFormat DEFAULT_QUALITY_FORMAT = FastqQualityFormat.Standard; 
	protected static final boolean DEFAULT_GZIP_OUTPUTS = true;
	protected static final boolean DEFAULT_KEEP_UNASSIGNED_READ = true;
	protected static final boolean DEFAULT_WRITER_FACTORY_USE_ASYNC_IO = true;
	protected static final boolean DEFAULT_ENSURE_IDENTICAL_HEADER_NAMES = true;
	protected static final String DEFAULT_UNASSIGNED_FILE_NAME_2 = "unassigned_2.txt";
	protected static final String DEFAULT_UNASSIGNED_FILE_NAME_1 = "unassigned_1.txt";
	protected static final String DEFAULT_METRICS_FILE_NAME = "jemultiplexer_out_stats.txt";
	
	
	
	/*
	 * We define here all common options
	 */
	
	@Option(shortName="F1", optional = false,
			printOrder=10,
			doc="Input fastq file (optionally gzipped) for single end data, or first read in paired end data.\n"
			)
	public File FASTQ_FILE1;

	@Option(shortName="F2", optional = true,
			printOrder=20,
			doc="Input fastq file (optionally gzipped) for the second read of paired end data.\n"
			)
	public File FASTQ_FILE2 = null;

	/*
	 * Barcode definition options
	 */
	@Option(shortName="BF", optional = false,  mutex={"USE_EMBASE"},
			printOrder=30,
			doc="Barcode file describing sequence list and sample names. "
					+"Tab-delimited file with 2 columns, with the sample in col1 and the corresponding barcode in col2.\n" +
					"Simple barcode file format : 2 tab-delimited colums\n"+
					"\t"+"If multiple barcode map to the same sample, either line can be duplicated e.g.\n" +
					"\t"+"\t"+"sample1\tATAT\n" +
					"\t"+"\t"+"sample1\tGAGG\n" +
					"\t"+"\t"+"sample2\tCCAA\n" +
					"\t"+"\t"+"sample2\tTGTG\n" +
					"\t"+"Or barcodes can be combined using the OR operator '|' i.e. the file above can be re-written like\n " +
					"\t"+"\t"+"sample1\tATAT|GAGG\n" +
					"\t"+"\t"+"sample2\tCCAA|TGTG\n" +
					"\t"+"Finally, for the special situation of paired-end data in which barcodes differ at both ends " +
					"(ie BPOS=BOTH BRED=false BM=BOTH , see BRED option description), barcodes for read_1 and read_2 can be" +
					" distinguished using a ':' separator i.e. \n" +
					"\t"+"\t"+"sample1\tATAT:GAGG\n" +
					"\t"+"\t"+"sample2\tCCAA:TGTG\n" +
					"\t"+"This above syntax means that sample 1 is encoded with ATAT barcode at read_1 AND GAGG barcode at read_2. " +
					"Note that you can still combine barcodes using | e.g. \n"+
					"\t"+"\t"+"sample1\tATAT|GAGG:CCAA|TGTG\n" +
					"\t"+"would mean that sample 1 is mapped by the combination of barcode: ATAT OR GAGG at read_1 AND CCAA OR TGTG at read_2.\n"+
					"Extended barcode file format : 3 (single-end) or 4 (paired-end) tab-delimited colums\n"+
					"\t"+"same as the simple barcode file format but the extra columns contains the file name(s) to use to name output files." +
					" A unique extra column is expected for single-end while 2 extra columns are expected for paired-end. In case, lines are duplicated (multiple barcodes" +
					"mapping the same sample), the same file name should be indicated in the third (and fourth) column(s). \n"+
					"\t"+"\t"+"sample1\tATAT\tspl1_1.txt.gz\tspl1_2.txt.gz\n" +
					"\t"+"\t"+"sample1\tGAGG\tspl1_1.txt.gz\tspl1_2.txt.gz\n" +
					"\t"+"\t"+"sample2\tCCAA\tspl2_1.txt.gz\tspl2_2.txt.gz\n" +
					"\t"+"Or\n" +
					"\t"+"\t"+"sample1 \t ATAT|GAGG:CCAA|TGTG \t spl1_1.txt.gz \t spl1_2.txt.gz\n" +
					"Ns in barcode sequence are allowed and are used to flag positions that should be ignored in sample matching \n"+
					"\t i.e. they will be clipped off the read sequence (like in iCLIP protocol)."
			)
	public File BARCODE_FILE = null;
	

	/*
	 * This option cannot be used when using INDEX_FILE.
	 */
	@Option(shortName="EM", optional = false, mutex={"BARCODE_FILE"},
			printOrder=400,
			doc="Enables emBASE mode i.e fetch information from emBASE and place demultiplexed files directly in emBASE repository structure.\n" +
					"This option is mutually exclusive with BARCODE_FILE.\n" +
					"Note : this option forces O=null GZ=true UN=true UF1=null UF2=null STATS_ONLY=false (all other user options supported).\n" 
			)
	public boolean USE_EMBASE = false;
	
	


	@Option(shortName="S", optional = true,
			printOrder=50,
			doc="For paired-end data and when two distinct barcodes/indices are used to encode samples,\n" +
					" this option tells if both barcodes should resolve to the same sample.\n"
					+" When true and if only one of the two reads has a barcode match, the read pair is 'unassigned'.\n"
					+" When false and if only one of the two reads has a barcode match, the read pair is assigned to the\n"
					+" corresponding sample\n"
					+ "When false and reads resolve to different samples, the read pair with the lowest mismatch is chosen\n "
					+ "(read pair is 'unassigned' in case of mismatch equality).\n"
			)
	public boolean STRICT = false;

	protected Integer BCLEN_1 = 0;
	protected Integer BCLEN_2 = 0;

	
	/*
	 * barcode matching options 
	 */

	@Option(shortName = "MM", optional = true,
			printOrder=60,
			doc="Maximum mismatches for a barcode to be considered a match. In situations where both barcodes are used for"
					+ " sample matching i.e. BPOS=BOTH BM=BOTH (or 2 INDEX_FILE given), two distinct\n" +
					" values can be given here using the syntax MM=X:Z where X and Z are 2 integers to use for read_1 and read_2 respectively.\n"
					+ "MM=null is like MM=0\n"
			)
	public String MAX_MISMATCHES = DEFAULT_MAX_MISMATCHES.toString();
	public Integer MAX_MISMATCHES_1 = DEFAULT_MAX_MISMATCHES;
	public Integer MAX_MISMATCHES_2 = DEFAULT_MAX_MISMATCHES;

	@Option(shortName = "MMD", optional = true,
			printOrder=70,
			doc="Minimum difference between the number of mismatches against the best and the second best barcode. When MMD is not respected, "+
					"the read remains unassigned.\n" +
					"When two distinct barcodes are used for sample matching (dual encoding), two distinct" +
					" values can be given using the syntax MMD=X:Z where X and Z are 2 integers to use for first (e.g. from read_1 or index_1)\n"+
					"MMD=null is like MMD=0\n"
			)
	public String MIN_MISMATCH_DELTA = DEFAULT_MIN_MISMATCH_DELTA.toString();
	public Integer MIN_MISMATCH_DELTA_1 = DEFAULT_MIN_MISMATCH_DELTA;
	public Integer MIN_MISMATCH_DELTA_2 = DEFAULT_MIN_MISMATCH_DELTA;


	@Option(shortName="Q", optional = true,
			printOrder=80,
			doc="Minimum base quality during barcode matching: bases which quality is less than this cutoff are always considered as a mismatch." +
					"When two distinct barcodes are used for sample matching (dual encoding), two distinct" +
					" values can be given using the syntax Q=X:Z where X and Z are 2 integers to use for first (e.g. from read_1 or index_1)"
					+ " and second barcode (e.g. from read_2 or index_2) respectively.\n"
					+ "Q=null is like Q=0.\n" 
			)
	public String MIN_BASE_QUALITY = DEFAULT_MIN_BASE_QUALITY.toString();
	public Integer MIN_BASE_QUALITY_1 = DEFAULT_MIN_BASE_QUALITY;
	public Integer MIN_BASE_QUALITY_2 = DEFAULT_MIN_BASE_QUALITY;


	/*
	 * read trimming options 
	 */


	@Option(shortName = "XT", optional = true,
			printOrder=90,
			doc = "Optional extra number of base to be trimmed right after the barcode (only used if CLIP_BARCODE=true). \n" +
					"When running paired-end, two distinct values can be given using the syntax XT=X:Z where X and Z are 2 integers" +
					" to use for read_1 and read_2 respectively. Note that even when BPOS=READ_1 or BPOS=READ_2, a X:Y synthax can be given to " +
					"trim the read w/o barcode as to end up with reads of the same length (note that this can also be operated using ZT). " +
					"If a unique value is given, e.g. XT=1, while running paired-end the following rule applies :\n " +
					"\t(1) BPOS=READ_1 or BPOS=READ_2, no trim is applied at the read w/o barcode \n"
					+ "\t(2) BPOS=BOTH, the value is used for both reads.\n"+
					"Note that XT=null is like XT=0.\n"
			)
	public String XTRIMLEN= DEFAULT_XTRIMLEN.toString() ; //must be further inspected to set XTRIMLEN_1 & XTRIMLEN_2
	public Integer XTRIMLEN_1=null ; // we do not set defaults as having null value here is part of the reasoning !!
	public Integer XTRIMLEN_2=null ; // we do not set defaults as having null value here is part of the reasoning !!

	@Option(shortName = "ZT", optional = true,
			printOrder=100,
			doc = "Optional extra number of bases to be trimmed from the read end i.e. 3' end.\n" +
					"When running paired-end, two distinct values can be given here using the syntax ZT=X:Z where X and Z are 2 integers" +
					" to use for read_1 and read_2 respectively. Note that even when BPOS=READ_1 or BPOS=READ_2, a X:Y synthax can be given to " +
					"trim the read w/o barcode as to end up with reads of the same length (note that this can also be operated using XT). " +
					"Note that if a single value is passed, the value always applies to both reads in paired-end mode without further consideration.\n"
					+ "ZT=null is like ZT=0.\n"
			)
	public String ZTRIMLEN=DEFAULT_ZTRIMLEN.toString() ;
	public Integer ZTRIMLEN_1=DEFAULT_ZTRIMLEN ;
	public Integer ZTRIMLEN_2=DEFAULT_ZTRIMLEN ;


	@Option(shortName="C", optional = true,
			printOrder=110,
			doc="Clip barcode sequence from read sequence, as well as XTRIMLEN (and ZTRIMLEN) bases if applicable, before writing to output file.\n"
					+" If false, reads are written without modification to output file. \n" +
					"Apply to both barcodes when BPOS=BOTH.\n"
			)
	public boolean CLIP_BARCODE = DEFAULT_CLIP_BARCODE;


	
	@Option(shortName = "RCHAR", optional = true,
			printOrder=120,
			doc="Replace spaces in read name/header using provided character. This is particularly handy when you need to retain" 
					+"\t ADDed barcode in read name/header during mapping (everything after space in read name is usually clipped in BAM files)." 
					+"\t"+"For example, with RCHAR=':' :\n"
					+"\t"+"\t" +"'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:'\n"
					+"\t" +"becomes\n"
					+"\t"+"\t" +"'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965:2:N:0:BARCODE'\n"
			)
	public String READ_NAME_REPLACE_CHAR = null;
	
	/*
	 * More general options below
	 */

	@Option(shortName="V", optional = true,
			printOrder=130,
			doc="A value describing how the quality values are encoded in the fastq.  Either 'Solexa' for pre-pipeline 1.3 " +
					"style scores (solexa scaling + 66), 'Illumina' for pipeline 1.3 and above (phred scaling + 64) or 'Standard' for phred scaled " +
					"scores with a character shift of 33.  If this value is not specified (or 'null' is given), the quality format will be detected.\n"
			)
	public FastqQualityFormat QUALITY_FORMAT = DEFAULT_QUALITY_FORMAT;

	
	
	@Option(shortName = "O", optional = true,
			printOrder=140,
			doc="Output directory. By default, output files are written in running directory.\n")
	public File OUTPUT_DIR = null;

	@Option(shortName="UN", optional = true,
			printOrder=150,
			doc="Should un-assigned reads be saved in files or simply ignored. File names are automatically created"
					+ " or can be given using UF1 & UF2 options.\n"
				)
	public boolean KEEP_UNASSIGNED_READ = DEFAULT_KEEP_UNASSIGNED_READ;


	@Option(shortName="UF1", optional=true, 
			printOrder=160,
			doc="Name of the file in which to write unassigned reads from FILE1." 
					+"Either a name (in which case the file will be created in the output dir) or full path.\n"
			)
	public String UNASSIGNED_FILE_NAME_1 = DEFAULT_UNASSIGNED_FILE_NAME_1;
	public File unassignedFile1 = null; //final File object to use with UNASSIGNED_FILE1, set in cmdline validation if necessary

	@Option(shortName="UF2", optional=true,
			printOrder=170,
			doc="Name of the file in which to write unassigned reads from FILE2."
					+"Either a name (in which case the file will be created in the output dir) or full path.\n"
			)
	public String UNASSIGNED_FILE_NAME_2 = DEFAULT_UNASSIGNED_FILE_NAME_2;
	public File unassignedFile2 = null;  //final File object to use with UNASSIGNED_FILE2, set in cmdline validation if necessary

	
	@Option(shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME, optional = true,
			printOrder=180,
			doc="File name where to write demultiplexing statistics. "
					+"Either a name (in which case the file will be created in the output dir) or an absolute path.\n"
			)
	public String METRICS_FILE_NAME = DEFAULT_METRICS_FILE_NAME;
	File metricsFile = null; //final File object to use with METRICS_FILE_NAME, set in cmdline validation 

	/**
	 * This is a tab-delimited 7 columns file:
	 * - Read_1_Header and Read_2_Header : fastq headers of read 1 (fwd read ie from fastq file 1) and 2 (rev read ie from fastq file 2)
	 * - seq1 and seq2 : seq as extracted from the reads 1 and 2 
	 * - bc1 and bc2 : barcode after remapping ; null if no remapping could be made
	 * - sample : final decision ie what sample was used for this read pair
	 */
	@Option(shortName = "DIAG", optional = true,
			printOrder=190,
			doc="Name for a barcode match reporting file (not generated by default)."
					+"Either a name (in which case the file will be created in the output dir) or full path."
					+" This file will contain a line per read pair with the barcode best matching the read subsequence "
					+"or 'null' when no match is found according to matching parameters ; and the final selected sample."
					+" This file is useful for debugging or further processing in case both ends are barcoded.\n"
					+ "N.B: this file will have a size of about one of the fastq input files.")
	String BARCODE_DIAG_FILE = null;
	File bcDiagFile = null; //final File object to use with BARCODE_DIAG_FILE, set in cmdline validation if necessary


	@Option(optional = true,
			printOrder=200,
			doc="Allows to overwrite existing files (system rights still apply).\n"
			)
	public boolean FORCE = false;
	
	
	@Option(shortName="GZ", optional = true,
			printOrder=210,
			doc="Compress output files using gzip.\n"
			)
	public boolean GZIP_OUTPUTS = DEFAULT_GZIP_OUTPUTS;


	@Option(optional=true,
			printOrder=220,
			doc="Change the default extension of created fastq files, eg 'fastqsanger'. By default uses the "
					+ "file extension from input fastq file. If result file names are given in the barcode file, "
					+ "this option is only used to adapt the unassigned file names. When using compression, a .gz is "
					+ "always appended to file names and should not be specified in FASTQ_FILE_EXTENSION i.e. \n"
					+ "use FASTQ_FILE_EXTENSION=fastq and NOT FASTQ_FILE_EXTENSION=fastq.gz\n"
			)
	public String FASTQ_FILE_EXTENSION = null;
	

	@Option(shortName="ASYNC", optional = true,
			printOrder=230,
			doc="Use one thread per Fastq Writer.\n")
	public boolean WRITER_FACTORY_USE_ASYNC_IO = DEFAULT_WRITER_FACTORY_USE_ASYNC_IO;

	@Option(optional = true,
			printOrder=240,
			doc="Only produces metric and diagnostic reports i.e. no output fastq file produced.\n")
	public boolean STATS_ONLY = false;



	protected static final SolexaQualityConverter solexaQualityConverter = SolexaQualityConverter.getSingleton();

	/**
	 * object in charge of validating the barcode file and holds useful maps ; to be set at cmd line validation
	 */
	protected BarcodeValidator barcodeValidator = null; 
	
	
	
	/*
	 * test mode : tell the code to stop as soon as the doWork() is called i.e. allows to test the parsing in order to check 
	 * the status after parsing.
	 */
	protected Boolean TEST_MODE_STOP_AFTER_PARSING = false;
	
	/**
	 * A boolean to indicate that sample barcodes come in separate index file(s) ie like in Illumina indexing
	 */
	protected Boolean USE_SAMPLE_INDEX_FILES = null;
	
	/**
	 * Are we running SE or PE ? Option set at cmd line parsing time
	 */
	protected Boolean RUNNING_PAIRED_END = false;

	
	/*
	 * barcodes as byte arrays for both 
	 */
	protected byte[][] barcodeBytes1 = null;
	protected byte[][] barcodeBytes2 = null;

	/*
	 * factory to create fastq writers
	 */
	protected JemultiplexerFastqWriterFactory fastqFactory = null;
	
	
	public AbstractJemultiplexer() {
		super();
		/*
		 * try load properties
		 */
		try {
			ApplicationConfiguration.init();
		} catch (IOException e) {
			log.error("Failed to read properties from property file.");
			throw new RuntimeException(e);
		}
	}

	/*
	 * Abstract methods that give access to command line arguments that we had to move (the cmd line arg I mean)
	 *  in implementing classes. Indeed, the argument have the same name but are used differently depending on the 
	 *  implementing class. Thus we needed to move them in implementing class to provide a documentation matching the 
	 *  situation. Note declaring these options here and overwriting them in implementing classes does not work  
	 *  But we still need access to these values here, this is obtained through abstract method 
	 */
	
	public abstract BarcodePosition getBarcodeReadPosOptionValue();
	public abstract boolean getRedundantBarcodesOptionValue();
	public abstract BarcodePosition getBarcodeForSampleMatchingOptionValue();
	public abstract String getBCLenOptionValue();
	public abstract boolean getAddBarcodeToHeaderOptionValue();
	
	
	/*
	 * VALIDATION methods
	 */
	
	
	protected List<String> validateFASTQFiles() {

		List<String> messages = new ArrayList<String>();

		if(!FASTQ_FILE1.exists())
			messages.add("FASTQ_FILE1 does not exist :"+FASTQ_FILE1.getAbsolutePath());

		if(FASTQ_FILE1.exists() && !FASTQ_FILE1.canRead())
			messages.add("Unsufficient rights to read FASTQ_FILE1 :"+FASTQ_FILE1.getAbsolutePath());

		if(FASTQ_FILE2 != null){

			if(!FASTQ_FILE2.exists())
				messages.add("FASTQ_FILE2 does not exist :"+FASTQ_FILE2.getAbsolutePath());

			if(FASTQ_FILE2.exists() &&!FASTQ_FILE2.canRead())
				messages.add("Unsufficient rights to read FASTQ_FILE2 :"+FASTQ_FILE2.getAbsolutePath());

			if(FASTQ_FILE2.getAbsolutePath().equals(FASTQ_FILE1.getAbsolutePath()))
				messages.add("FASTQ_FILE1 and FASTQ_FILE2 are the same file !");
		}
		return messages;
	}
	
	

	/**
	 * Perform validation of the barcode positioning options 
	 * 
	 * @return the list of error messages or empty list
	 */
	protected abstract List<String> validateBarcodePositions();

	
	/**
	 * Perform barcode file validation
	 * 
	 * @return the list of error messages or empty list
	 */
	protected abstract List<String> validateBarcodes();
	
	/**
	 * Tells if options that optionally accept 2 values separated by a ':' can actually 
	 * accept two values *in the current setup*
	 * 
	 * @return null if OK or a error message
	 */
	protected abstract String validateOptionAcceptsTwoValues(String optionName, String optionValue, boolean failIfPEwithSingleBarcode);
	

	/** 
	 * Called at cmd line validation when the USE_EMBASE switch is used 
	 * 
	 * @return the list of error messages or empty list
	 */
	protected abstract List<String> prepareForEmBASE() ;
	
	
	/**
	 * 
	 * To be call by implementating classes
	 * Validate that all BARCODEs are the same length and unique
	 * 
	 *
	 * @return null if command line is valid.  If command line is invalid, returns an array of error message
	 *         to be written to the appropriate place.
	 */
	protected String[] customCommonCommandLineValidation(List<String> messages) {
		
		
		/*
		 * Check barcode option are consitent
		 * 
		 */
		List<String> _mess = validateBarcodePositions();
		if(_mess.size() > 0){
			//fail directly
			messages.addAll(_mess);
			return messages.toArray(new String[messages.size()]);
		}
			
		
		/*
		 * Intercept USE_EMBASE to initialize BF and other things
		 */
		if(USE_EMBASE){
			//use embase mode not supported with index files yet
//			if(USE_SAMPLE_INDEX_FILES){
//				messages.add("[with sample index files] USE_EMBASE mode not currently supported with index files");
//				return messages.toArray(new String[messages.size()]);
//			}
			messages.addAll( prepareForEmBASE() );
		}
		
		
		
		/*
		 * parse options that can have 'X:Y' syntax
		 * We first set these options blindly, later check is performed for those options that need it ie:
		 * - BCLEN* => make sure to have non 0 values at ends with barcodes and 0 values at ends w/o barcodes
		 * - XTRIMLEN => must be further inspected i.e. make sure we don t have a '1' for ends w/o 
		 * barcodes UNLESS explicitly set by user in the form X:Y 
		 */
		int [] arr = null;
		String BCLEN = getBCLenOptionValue(); // ask implementing class its BCLEN arg
		arr = setOption("BCLEN", BCLEN,  messages, true);
		if(arr.length == 2){ 
			//BCLEN was given in cmd line as it is null by default
			BCLEN_1 = arr[0];
			log.info("setting BCLEN_1 to "+BCLEN_1+" (BCLEN="+BCLEN);
			BCLEN_2 = arr[1];
			log.info("setting BCLEN_2 to "+BCLEN_2+" (BCLEN was "+BCLEN);
		}else{
			//BCLEN was NOT given in cmd line ... let it init with barcode file...
		}
				
		arr = setOption("MAX_MISMATCHES", MAX_MISMATCHES,  messages, true);
		if(arr.length == 2){
			MAX_MISMATCHES_1 = arr[0];
			MAX_MISMATCHES_2 = arr[1];
		}
		
		arr = setOption("MIN_MISMATCH_DELTA", MIN_MISMATCH_DELTA, messages, true);
		if(arr.length == 2){
			MIN_MISMATCH_DELTA_1 = arr[0];
			MIN_MISMATCH_DELTA_2 = arr[1];
		}
		arr = setOption("MIN_BASE_QUALITY", MIN_BASE_QUALITY,  messages, true);
		if(arr.length == 2){
			MIN_BASE_QUALITY_1 = arr[0];
			MIN_BASE_QUALITY_2 = arr[1];
		}
		
		arr = setOption("XTRIMLEN", XTRIMLEN,  messages, false);
		if(arr.length == 2){
			XTRIMLEN_1 = arr[0];
			XTRIMLEN_2 = arr[1];
		}
		
		arr = setOption("ZTRIMLEN", ZTRIMLEN, messages, false);
		if(arr.length == 2){
			ZTRIMLEN_1 = arr[0];
			ZTRIMLEN_2 = arr[1];
		}
		
		//check XTRIMLEN (not needed when using INDEX_FILE as we already checked it before)
		if(!USE_SAMPLE_INDEX_FILES && RUNNING_PAIRED_END && getBarcodeReadPosOptionValue() != BarcodePosition.BOTH){
			//set XTRIMLEN to 0 for the end w/o barcode UNLESS set by user in command line
			if(!XTRIMLEN.contains(":")){
				//user did not set 2 values => make it 0 
				if(getBarcodeReadPosOptionValue() == BarcodePosition.READ_1)
					XTRIMLEN_2 = 0;
				else
					XTRIMLEN_1 = 0;
			}
		}
		
		
		messages.addAll( validateBarcodes() );

		/*
		 * check output dir
		 */
		if(!USE_EMBASE){
			
			if(OUTPUT_DIR == null){
				OUTPUT_DIR = new File(System.getProperty("user.dir"));
			}
			
			if(!OUTPUT_DIR.exists()){
				log.info("Attempting to create output directory : "+OUTPUT_DIR.getAbsolutePath());
				FileUtil.checkWritableDir(OUTPUT_DIR, true);
			}
		}

		/*
		 * Check quality format
		 */
		if (QUALITY_FORMAT == null) {
			FastqReader reader = new FastqReader(FASTQ_FILE1);
			QUALITY_FORMAT = QualityEncodingDetector.detect(100000, reader);
			log.info(String.format("Auto-detected quality format as: %s.", QUALITY_FORMAT));
		}

		

		if (messages.size() == 0) {
			
			/*
			 * All is ok, finish understanding the command line options 
			 */

			if(METRICS_FILE_NAME == null) METRICS_FILE_NAME=DEFAULT_METRICS_FILE_NAME; //this is user mistake, reset the default name
			metricsFile = setFile(METRICS_FILE_NAME);

			if(BARCODE_DIAG_FILE != null)
				bcDiagFile = setFile(this.BARCODE_DIAG_FILE);

			if(KEEP_UNASSIGNED_READ){
				if(UNASSIGNED_FILE_NAME_1 == null) UNASSIGNED_FILE_NAME_1=DEFAULT_UNASSIGNED_FILE_NAME_1; //this is user mistake, reset the default name
				if(UNASSIGNED_FILE_NAME_2 == null && RUNNING_PAIRED_END) UNASSIGNED_FILE_NAME_2=DEFAULT_UNASSIGNED_FILE_NAME_2; //this is user mistake, reset the default name

				if(FASTQ_FILE_EXTENSION != null){
					//is the unassignedFile1 the default name ? If so , we replace the extension ; else user knows what he is doing ...
					if(UNASSIGNED_FILE_NAME_1.equals(DEFAULT_UNASSIGNED_FILE_NAME_1)){
						UNASSIGNED_FILE_NAME_1 = DEFAULT_UNASSIGNED_FILE_NAME_1.replace(".txt", "."+FASTQ_FILE_EXTENSION);
						if(RUNNING_PAIRED_END){
							//blindly force name for rev 
							UNASSIGNED_FILE_NAME_2 = DEFAULT_UNASSIGNED_FILE_NAME_2.replace(".txt", "."+FASTQ_FILE_EXTENSION);
						}
					}
				} 
				
				
				if(GZIP_OUTPUTS && ! UNASSIGNED_FILE_NAME_1.endsWith(".gz"))
					UNASSIGNED_FILE_NAME_1 += ".gz";

				unassignedFile1 = setFile(UNASSIGNED_FILE_NAME_1);
				if(RUNNING_PAIRED_END){
					if(GZIP_OUTPUTS && ! UNASSIGNED_FILE_NAME_2.endsWith(".gz"))
						UNASSIGNED_FILE_NAME_2 += ".gz";
					unassignedFile2 = setFile(UNASSIGNED_FILE_NAME_2);
				}
			}

			
			return null;
		}
		
		return messages.toArray(new String[messages.size()]);
	}


	/**
	 * builds a file to use for barcode file
	 * 
	 * @param fastqFile the fastq file (mate 1 in case of PE) 
	 * @param username the user how is about to demultiplex
	 * @param d the dir to use or null to use system defaults
	 * @return
	 */
	protected File getBarcodeFileTmpLocation(File fastqFile, String username, File d) {
		String tmpdir = d == null ? System.getProperty("java.io.tmpdir") : d.getAbsolutePath();
		if(!StringUtil.isValid(tmpdir) || !new File(tmpdir).exists())
			tmpdir = ".";
		
		File f= new File(
				tmpdir,
				FileUtil.removeExtension( fastqFile.getName() ) +"_"+username+"_barcodes.txt"
				);
		log.debug("barcode file to use is "+f.getAbsolutePath());
		return f;
	}
	

	/**
	 * Deals with integer options that can take X:Y values ; also check that given values are >= 0 
	 * @param optionName name of the option (for error messages)
	 * @param optionValue as given by user ; if blank returned array is empty
	 * @param messages the error message list to use in case of error
	 * @param failIfPEwithSingleBarcode
	 * @return an int array with the two parsed values for read_1 and read_2. 
	 * Array is empty is nothing should be made ie in such case one should check if messages size increased
	 */
	protected int[] setOption(String optionName, String optionValue, List<String> messages, boolean failIfPEwithSingleBarcode) {
		
		int[] arr = new int[2];
		
		if(StringUtils.isBlank(optionValue))
			return arr;
		
		String [] values = new String[2]; 
		if(optionValue.contains(":")){
			//validate 
			String _mess = validateOptionAcceptsTwoValues(optionName, optionValue, failIfPEwithSingleBarcode);
			if(_mess!=null){
				messages.add(_mess);
				return arr ;
			}
			//split values
			values = optionValue.split(":");
			if(values.length!=2){
				messages.add("Invalid value used for "+optionName+" : "+optionValue +". The X:Z synthax can contain only a single ':' ");
				return arr;
			}
		}
		else{
			values[0] = optionValue;
			values[1] = optionValue;
		}
		
		//now convert to int
		int optionRead1 = -1; 
		int optionRead2 = -1;
		try{
			optionRead1 = Integer.parseInt(values[0]);
		}catch (NumberFormatException e) {
			messages.add("Invalid value used for "+optionName+" : "+optionValue +" => "+values[0]+" is not an integer");
		}
		
		try{
			optionRead2 = Integer.parseInt(values[1]);
		}catch (NumberFormatException e) {
			messages.add("Invalid value used for "+optionName+" : "+optionValue +" => "+values[1]+" is not an integer");
		}
		
		log.debug(optionName+" => [1]="+optionRead1+" ; [2]="+optionRead2);
		
		if(optionRead1<0 || optionRead2<0)
			messages.add("Invalid value used for "+optionName+" : "+optionValue +" => value(s) must be >= 0.");
		
		if(RUNNING_PAIRED_END && getBarcodeReadPosOptionValue()==BarcodePosition.BOTH && getRedundantBarcodesOptionValue() && optionRead1 != optionRead2 ){
			log.warn("POTENTIAL OPTION MISTAKE : Found 2 different values for option "+optionName+" (given cmd line value was '"+optionValue
					+"') while REDUNDANT_BARCODES=true: "+optionRead1+" for read_1 and "+optionRead2+" for read_2. Code continues for now.");
		}
		
		return new int []{optionRead1, optionRead2};
	}


	protected File setFile(String fnameOrPath) {
		if(fnameOrPath.contains("/")){
			//it is a path, use file path as provided
			return new File(fnameOrPath);
		}
		else{
			//it is just a name
			return new File(OUTPUT_DIR, fnameOrPath);
		}

	}

	
	
	/**
	 * 
	 * Find out what file should be used for a given library. This method does NOT fetches existing file
	 * but builds a file path to use by convention i.e. when you want to automatically demultiplex fastq file
	 *  
	 * 
	 * @param lib the lib for which we need to know the file to use to write demultiplexed reads.
	 * Note that lib should have its final sampleDir correctly set 
	 * @param readPairNumber 1 or 2 (for single end, it is always 1) 
	 * @param gzip whether output should be gzipped
	 * @return the File to use
	 */
	protected File getDemultiplexedFile(NGSLibrary lib, int readPairNumber, boolean gzip) {
		
		//for now simply build up a name with lib name
		File f =  new File(lib.getSampleDir(), lib.getName()+"_"+lib.getSampleDir().getParentFile().getName()+"_"+readPairNumber+".txt"+(gzip?".gz":"")); 
		log.debug("Demultiplexed file for "+lib.getName()+" (read pair "+readPairNumber+") should be "+f.getAbsolutePath());
		return f;
	}

	/*
	 * 
	 * COMMON DEMULTIPLEXING METHODS
	 * 
	 * 
	 */
	
	
	/**
	 * internal variable storing the barcode to sample-name map, 
	 * initialized after successful validation in doWork() 
	 */
	protected Map<String, String> bc2sample = null;
	
	/**
	 * internal variable storing the sample-name to list of barcodes map, 
	 * initialized after successful validation in doWork() 
	 */
	protected Map<String, Set<String>> sample2bcSet = null;
	
	/**
	 * internal variable storing the sample-name to fastq writer map, 
	 * initialized after successful validation in doWork() 
	 */
	protected HashMap<String, FastqWriter> writers_fwd = null;
	protected HashMap<String, FastqWriter> writers_rev = null;
	
	
	/**
	 * internal variable storing map storing how many read pair are assigned to each sample 
	 * initialized after successful validation in doWork() 
	 */
	protected HashMap<String, Integer> sampleCountMap = null;
	
	@Override
	protected int doWork() {
		
		/*
		 * We perform here all common initialization work 
		 */
		
		if(TEST_MODE_STOP_AFTER_PARSING)
			return 0;
		
		//should we write a barcode match result file ?
		boolean writingBarcodeDiagFile = (this.bcDiagFile != null);

		
		/*
		 * Initialize all barcode maps 
		 */
		barcodeValidator.parse();
	
		bc2sample = barcodeValidator.barcode2sample;
		sample2bcSet = barcodeValidator.sample2barcodes;
		//fill in barcodeBytes array
		List<String> barcodeList1 = new ArrayList<String>(barcodeValidator.barcodeSetRead1);
		barcodeBytes1 = new byte[barcodeList1.size()][];
		for (int i = 0; i < barcodeList1.size(); i++) {
			barcodeBytes1[i] = htsjdk.samtools.util.StringUtil.stringToBytes(barcodeList1.get(i));
		}
		List<String> barcodeList2 = new ArrayList<String>(barcodeValidator.barcodeSetRead2);
		barcodeBytes2 = new byte[barcodeList2.size()][];
		for (int i = 0; i < barcodeList2.size(); i++) {
			barcodeBytes2[i] = htsjdk.samtools.util.StringUtil.stringToBytes(barcodeList2.get(i));
		}

		/*
		 * prepare FASTQ file writers
		 */
		fastqFactory = new JemultiplexerFastqWriterFactory();
		fastqFactory.setUseAsyncIo(WRITER_FACTORY_USE_ASYNC_IO);

		//init a writer and counter for each sample
		writers_fwd = new HashMap<String, FastqWriter>();
		writers_rev = new HashMap<String, FastqWriter>();
		
		//map storing how many read pair are assigned to each sample
		sampleCountMap = new HashMap<String, Integer>(); 

		//do we use long barcode format?
		Map<String, String> sample2file1 = null;
		Map<String, String> sample2file2 = null;
		if(barcodeValidator.isLongFormat && !STATS_ONLY){
			//yes load file name to use
			sample2file1 = barcodeValidator.sample2file1;
			if(RUNNING_PAIRED_END)
				sample2file2 = barcodeValidator.sample2file2;
		}
		for (Entry<String, Set<String>> e : sample2bcSet.entrySet()) {
			String sample = e.getKey();
			sampleCountMap.put(sample, 0);
			List<String> barcodes = new ArrayList<String>( e.getValue() );
			//default names
			String fname_fwd = sample + "_"+StringUtil.mergeList(barcodes, "-")+"_1." + (FASTQ_FILE_EXTENSION == null ? "txt" : FASTQ_FILE_EXTENSION);
			String fname_rev = sample + "_"+StringUtil.mergeList(barcodes, "-")+"_2." + (FASTQ_FILE_EXTENSION == null ? "txt" : FASTQ_FILE_EXTENSION);
			//overwrite with user provided file names if available
			fname_fwd = sample2file1 == null || !sample2file1.containsKey(sample) ? fname_fwd : sample2file1.get(sample);
			fname_rev = sample2file2 == null || !sample2file2.containsKey(sample) ? fname_rev : sample2file2.get(sample);
			if(GZIP_OUTPUTS){
				//make sure files ends with .gz
				if(!fname_fwd.endsWith(".gz"))
					fname_fwd += ".gz";
				if(fname_rev!=null && !fname_rev.endsWith(".gz"))
					fname_rev += ".gz";
			}
			if(!STATS_ONLY){
				log.info("Registering writers for "+sample+" :");
				log.info("\tFile 1 Name : "+fname_fwd);
				File fwdFile = null;
				if(looksLikeAPath(fname_fwd)){
					fwdFile = new File(fname_fwd);
				}else{
					log.info("Using OUTPUT_DIR : "+OUTPUT_DIR.getAbsolutePath());
					fwdFile = new File(OUTPUT_DIR, fname_fwd);
				}
				//check if file exists
				if(fwdFile.exists()){
					if(!FORCE){
						log.error(
								"\nOutput file already exists : " + fwdFile.getAbsolutePath() +
								"\nNote that other output files not indicated here might also exist."+
								"\nPlease use FORCE=true to allow overwrite."
								);
						System.exit(1); // abnormal exit
					}else{
						//delete existing file
						if(!fwdFile.delete()){
							log.error(
									"\nFailed to remove existing output file (FORCE=true) : " + fwdFile.getAbsolutePath() +
									"\nPlease make sure you have write access to *all* existing output files (adapt file rights or delete them first)."
									);
							System.exit(1); // abnormal exit
						}else{
							log.info("\tExisting output file successfully deleted : "+fwdFile.getAbsolutePath());
						}
					}
				}
				//register writer 
				writers_fwd.put(sample, fastqFactory.newWriter(fwdFile, GZIP_OUTPUTS, CREATE_MD5_FILE));
				
				if(RUNNING_PAIRED_END){
					log.info("\tFile 2 Name : "+fname_rev);
					File revFile = null;
					if(looksLikeAPath(fname_rev)){
						revFile = new File(fname_rev);
					}else{
						revFile = new File(OUTPUT_DIR, fname_rev);
					}
					//check if file exists
					if(revFile.exists() ){
						if(!FORCE){
							log.error(
									"\nOutput file already exists : " + revFile.getAbsolutePath() +
									"\nNote that other output files not indicated here might also exist."+
									"\nPlease use FORCE=true to allow overwrite."
									);
							System.exit(1); // abnormal exit
						}else{
							//delete existing file
							if(!revFile.delete()){
								log.error(
										"\nFailed to remove existing output file (FORCE=true) : " + revFile.getAbsolutePath() +
										"\nPlease make sure you have write access to *all* existing output files (adapt file rights or delete them first)."
										);
								System.exit(1); // abnormal exit
							}else{
								log.info("\tExisting output file successfully deleted : "+revFile.getAbsolutePath());
							}
						}
					}
					//register writer 
					writers_rev.put(sample, fastqFactory.newWriter(revFile, GZIP_OUTPUTS, CREATE_MD5_FILE));
					
				}
			}
		}
		//add writers for the unassigned reads if required
		if(KEEP_UNASSIGNED_READ && !STATS_ONLY){
			log.info("Registering writers unassigned reads :");
			log.info("\tFile 1 Name : "+unassignedFile1);
			
			writers_fwd.put(UNASSIGNED, fastqFactory.newWriter(unassignedFile1, GZIP_OUTPUTS, CREATE_MD5_FILE));
			if(RUNNING_PAIRED_END){
				writers_rev.put(UNASSIGNED, fastqFactory.newWriter(unassignedFile2, GZIP_OUTPUTS, CREATE_MD5_FILE));
				log.info("\tFile 2 Name : "+unassignedFile2);
			}
		}

		/*
		 * Should we init the diag file ?
		 */
		PrintWriter diagPW = null;
		if( writingBarcodeDiagFile ){
			try{
				diagPW = new PrintWriter(bcDiagFile);
			} catch (Exception e) {
				log.error(ExceptionUtil.getStackTrace(e));
				throw new RuntimeException(e);
			}
			/*
			 * init headers with SE  :
			 * - Read_Header (is read header as found in FASTQ) ; 
			 * - "seq" : read sub-sequence to match against barcodes ; 
			 * - "MM_Best", "MM_Second": Mismatch number with best and second best match
			 * - "bc" : chosen barcode sequence ; this is ofcourse function of global options
			 * - "sample" : sample ultimately selected ie the read should be found in the corresponding file
			 * For PE case only :
			 * - same as above but duplicated with ending _1 and _2 for read 1 and 2, respectively ; 
			 * - "smpl_1" and "smpl_2" : sample match for bc1 and bc2
			 */
			String[] headers = new String []{"Read_Header", "seq", "MM_Best", "MM_Second", "bc", "sample"};
			if(RUNNING_PAIRED_END)
				headers = new String []{"Read_1_Header", "Read_2_Header", "seq1", "seq2" , "MM_Best_1", "MM_Second_1", "bc1", "smpl_1", "MM_Best_2", "MM_Second_2", "bc2", "smpl_2", "sample"};

			diagPW.println(StringUtil.mergeArray(headers,"\t"));
		}

		
		/*
		 * call the case specific implementation
		 */
		return doDemultiplexingWork(diagPW);
	}
	
	protected abstract int doDemultiplexingWork(PrintWriter diagPW); 
	

	
	/**
	 * Print the final demultiplexing report
	 * The method prints reports in {@link AbstractJemultiplexer#metricsFile}
	 * and individual sample counts are from {@link AbstractJemultiplexer#sampleCountMap}
	 * map that is filled up by the demultiplexing implementation
	 * @param demultiplexer the actual implementation calling the method
	 * @param demultiplexerName the implementation name as known by the user i.e. 'demultiplex' or 'demultiplex-illu' 
	 * @param total_read_count 
	 * @param unassigned number of read (pair) unassigned
	 * @param assigned  number of read (pair) assigned
	 */
	/**
	 * @param demultiplexer
	 * @param demultiplexerName
	 * @param total_read_count
	 * @param unassigned
	 * @param assigned
	 */
	protected void printMetricFile(AbstractJemultiplexer demultiplexer, String demultiplexerName, int total_read_count, int unassigned,
			int assigned) {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(metricsFile);
			
				
			pw.println("## Created on " + new Date() + " by je "+demultiplexerName +" ["+demultiplexer.getClass()+"] version "+getVersion());
			pw.println("## Command line : "+demultiplexer.getCommandLine());           
           
			pw.println("Processed Reads (pairs)"+"\t"+ total_read_count);
			pw.println("Assigned Reads (pairs)"+"\t"+ assigned );
			pw.println("Unassigned Reads (pairs)"+"\t"+ unassigned);
			pw.println("# Individual sample read (pair) counts :");
			for (Entry<String, Integer> e : sampleCountMap.entrySet()) {
				pw.println(e.getKey()+"\t"+ e.getValue() );
			}

		} catch (Exception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			throw new RuntimeException(e);
		} finally {
			if(pw!=null)
				pw.close();
		}
	}
	
	
	
	
	/**
	 * 
	 * Finds best matching barcode according to global settings i.e. 
	 * MINIMUM_BASE_QUALITY
	 * MIN_MISMATCH_DELTA
	 * MAX_MISMATCHES
	 * 
	 * @param barcodePosition tell what barcode we are currently matching i.e. either READ_1 or READ_2 
	 * @param readSubsequence the read subsequence to match against barcodes
	 * @param qualitySubsequence the quality subsequence corresponding to s
	 * @param bc2sample
	 * @param sample2bcSet
	 * @return a {@link BarcodeMatch} describing the match or null if no match was found
	 */
	protected BarcodeMatch findBarcode(
			final BarcodePosition barcodePosition,
			final String readSubsequence,
			final String qualitySubsequence,
			final Map<String, String> bc2sample,
			final Map<String, Set<String>> sample2bcSet) {

		if(readSubsequence == null)
			return null;
		
		if(barcodePosition == BarcodePosition.BOTH)
			throw new RuntimeException("findBarcode method can t be called with barcodePosition = BOTH");

		//fecth set of allowed barcodes at this positin
		Set<String> possibleBarcodes = barcodePosition == BarcodePosition.READ_1 ? barcodeValidator.barcodeSetRead1 : barcodeValidator.barcodeSetRead2;
		
		boolean barcodesHaveNs = barcodeValidator.areNsInBarcodes(barcodePosition);
		/*
		 * we first check if the subsequence matches a known barcode, in which case, we directly take this as a result, 
		 * BUT this cannot work if we have Ns in barcodes (like in iCLIP), to speed up code we directly use a predefined 
		 */
		if(!barcodesHaveNs && possibleBarcodes.contains(readSubsequence)){
			BarcodeMatch bm = new BarcodeMatch();
			bm.matched = true;
			bm.barcode = readSubsequence;
			bm.mismatches = 0;
			bm.mismatchesToSecondBest = readSubsequence.length();
			return bm;
		}else{
			BarcodeMatch bm =  findBestBarcode(barcodePosition, readSubsequence, qualitySubsequence); 
			if(bm.matched)
				return bm;

		}
		return null;
	}



	/**
	 * Find the best barcode match for the given read sequence
	 * @param barcodePosition the {@link BarcodePosition} 
	 * @param readSubsequence portion of read containing barcode
	 * @param qualityScoreString the string representing the quality scores for the read subsequence
	 * @return a BarcodeMatch holding the matched status and extra info when matched is true
	 */
	protected BarcodeMatch findBestBarcode(final BarcodePosition barcodePosition, final String readSubsequence, final String qualityScoreString) {

		//turn to bytes
		byte [] subseqBytes = htsjdk.samtools.util.StringUtil.stringToBytes(readSubsequence);
		byte[] qualBytes = convertQualityStringToByteArray(qualityScoreString);

		int numMismatchesInBestBarcode = readSubsequence.length() + 1; //init with max mismatch num + 1
		int numMismatchesInSecondBestBarcode = readSubsequence.length() + 1;
		//find best and second best barcode
		byte[] bestBarcode = null;
		
		//look up adequate variables
		byte [][] barcodeBytes = barcodePosition == BarcodePosition.READ_1 ? barcodeBytes1 : barcodeBytes2;
		int minQuality = barcodePosition == BarcodePosition.READ_1 ? MIN_BASE_QUALITY_1 : MIN_BASE_QUALITY_2;
		int maxMismatches = barcodePosition == BarcodePosition.READ_1 ? MAX_MISMATCHES_1 : MAX_MISMATCHES_2;
		int minMMDelta = barcodePosition == BarcodePosition.READ_1 ? MIN_MISMATCH_DELTA_1 : MIN_MISMATCH_DELTA_2;
		
		for (int i = 0; i < barcodeBytes.length; i++) {
			final int numMismatches = countMismatches(barcodeBytes[i], subseqBytes, qualBytes, minQuality);
			if (numMismatches < numMismatchesInBestBarcode) {
				if (bestBarcode != null) {
					numMismatchesInSecondBestBarcode = numMismatchesInBestBarcode;
				}
				numMismatchesInBestBarcode = numMismatches;
				bestBarcode = barcodeBytes[i];
			} else if (numMismatches < numMismatchesInSecondBestBarcode) {
				numMismatchesInSecondBestBarcode = numMismatches;
			}
		}


		final boolean matched = (bestBarcode != null &&
				numMismatchesInBestBarcode <= maxMismatches &&
				numMismatchesInSecondBestBarcode - numMismatchesInBestBarcode >= minMMDelta)
				;


		final BarcodeMatch match = new BarcodeMatch();
		match.matched = matched;

		if (matched) {
			match.mismatches = numMismatchesInBestBarcode;
			match.mismatchesToSecondBest = numMismatchesInSecondBestBarcode;
			match.barcode = htsjdk.samtools.util.StringUtil.bytesToString( bestBarcode );
		}


		return match;
	}



	/**
	 * This considers a byte offset of 33 i.e. Phred+33 (Illumina compatible)
	 * Method inspired from {@link FastqToSam}
	 * @param qualityScores the quality scire string from Illumina
	 * @return
	 */
	protected byte[] convertQualityStringToByteArray(final String qualityScores) {
		byte [] qualBytes = htsjdk.samtools.util.StringUtil.stringToBytes(qualityScores);
		convertQuality(qualBytes, QUALITY_FORMAT);
		return qualBytes;
	}

	/** 
	 * Based on the type of quality scores coming in, converts them to a numeric byte[] in phred scale. 
	 */
	void convertQuality(final byte[] quals, final FastqQualityFormat version) {
		switch (version)  {
		case Standard:
			SAMUtils.fastqToPhred(quals);
			break ;
		case Solexa:
			solexaQualityConverter.convertSolexaQualityCharsToPhredBinary(quals);
			break ;
		case Illumina:
			solexaQualityConverter.convertSolexa_1_3_QualityCharsToPhredBinary(quals);
			break ;
		}
	}



	/**
	 * 
	 * Compare barcode sequence to bases from read
	 * @param barcodeBytes the barcode as array of bytes
	 * @param readSubsequence the subsequence as array of bytes
	 * @param qualities the base quality score in phred scale (+33)
	 * @param minimumBaseQuality an int in [0,40]
	 * @return how many bases did not match , 'N' or base with quality less than MINIMUM_BASE_QUALITY always considered
	 * as a mismatch 
	 * 
	 */
	protected  int countMismatches(final byte[] barcodeBytes, final byte[] readSubsequence, final byte[] qualities, final int minimumBaseQuality) {
		int numMismatches = 0;
		// Read sequence and barcode length may not be equal, so we just use the shorter of the two ; I don t see why this should happen but it is free to do this !
		final int basesToCheck = Math.min(barcodeBytes.length, readSubsequence.length);
		for (int i = 0; i < basesToCheck; ++i) {
			/*
			 * we first need to check this is a valid position in barcode (barcode can contain N to specify this position should be ignored)
			 * If the bc contains a N, we skip
			 */
			if (!SequenceUtil.isValidBase(barcodeBytes[i])) { // Returns true if the byte is in [acgtACGT].
				continue;
			}
			
			if (!SequenceUtil.isNoCall(readSubsequence[i])) { //returns true if the value of base represents a no call
				//if bases in barcode and readsubsequence are different => increase mismatch count ; else check quality
				if (!SequenceUtil.basesEqual(barcodeBytes[i], readSubsequence[i])){
					++numMismatches;
				}
				else if (qualities != null){
					final int uQual = qualities[i] & 0xff; //from FastqToSam.java ; this is kind of bit masking 
					if(uQual < minimumBaseQuality)
						++numMismatches;
				}
			}else{
				//the read subsequence contains a 'N' at this position => this is a mismatch 
				++numMismatches;
			}
		}

		return numMismatches;
	}

	
	/**
	 * @param r the read to write
	 * @param bcLen the barcode length
	 * @param xtrim how many extra base should be clipped after the barcode ; 0 or more
	 * @param ztrim how many extra base should be clipped from the bread's end ; 0 or more
	 * @param writer the writer
	 * @param barcode the barcode sequence extracted from this read. Note that the caller is responsible to give the 
	 *  sequence as extracted from the read or the barcode sequence after matching the extracted sequence against list of expected barcode
	 * @param addBarcodeToHeader whether the barcode should be added to the read header. If true, the string ':barcode' is added 
	 * with the starting ':' added only if current read header does not end with ':'. For example :
	 * '@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:'
	 * becomes
	 * '@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:BARCODE'
	 * @param clipHeaderAfterSpace if true the '2:N:0' part is removed before adding the barcodes
	 * @return the trimmed {@link FastqRecord} ; mostly for testing purposes
	 */
	protected FastqRecord write(final FastqRecord r, int bcLen, int xtrim, int ztrim, final FastqWriter writer, String barcode, Boolean addBarcodeToHeader, Boolean clipHeaderAfterSpace) {
		
		log.debug("Writing read [bcLen="+bcLen+", xtrim="+xtrim+", ztrim="+ztrim+", barcode="+(barcode==null?"NULL":barcode)
				+", addBarcodeToHeader="+addBarcodeToHeader.toString()+", clipHeaderAfterSpace="+clipHeaderAfterSpace.toString()+"]");
		
		int l = bcLen + xtrim ;
		
		String header = null;
		
		if(clipHeaderAfterSpace){
			header = r.getReadHeader().split("\\s+")[0] + " ";
		}else{
			header = r.getReadHeader();
			
			if(addBarcodeToHeader && barcode!=null && !header.endsWith(":"))
				header+= ":";
				
		}
		log.debug(READ_NAME_REPLACE_CHAR);
		if(addBarcodeToHeader){
			if(barcode!=null){
				header+= barcode;
			}  
				
			if(READ_NAME_REPLACE_CHAR != null){
				header = header.replaceAll(" ", READ_NAME_REPLACE_CHAR);
			}
		}
		log.debug(header);
		log.debug("will clip "+l+"bases of "+r.getReadString());

		String subseq_trimmed = r.getReadString();
		String qualstring_trimmed = r.getBaseQualityString();
		if(CLIP_BARCODE && ztrim > 0){
			subseq_trimmed = r.getReadString().substring(l, subseq_trimmed.length() - ztrim) ;
			qualstring_trimmed = r.getBaseQualityString().substring(l, qualstring_trimmed.length() - ztrim) ;
		}
		else if(CLIP_BARCODE){
			subseq_trimmed = r.getReadString().substring(l) ;
			qualstring_trimmed = r.getBaseQualityString().substring(l) ;
		}
		else if(ztrim > 0){
			subseq_trimmed = r.getReadString().substring(0, subseq_trimmed.length() - ztrim) ;
			qualstring_trimmed = r.getBaseQualityString().substring(0, qualstring_trimmed.length() - ztrim) ;
		}


		FastqRecord trimmed = new FastqRecord(
				header,//seqHeaderPrefix
				subseq_trimmed, //seqLine,
				r.getBaseQualityHeader(),//qualHeaderPrefix,
				qualstring_trimmed // qualLine
				);

		if(writer != null) //mostly here for testing ie to be able to call method w/o writer
			writer.write(trimmed);

		return trimmed;
	}
	
	
	/**
	 * @param r
	 * @param ztrim
	 * @param header
	 * @param l
	 * @return
	 */
	protected FastqRecord prepareTrimmedRead(FastqRecord r, Integer ztrim,
			String header, int l) {
		
		String subseq_trimmed = r.getReadString();
		String qualstring_trimmed = r.getBaseQualityString();
		if(CLIP_BARCODE && ztrim > 0){
			subseq_trimmed = r.getReadString().substring(l, subseq_trimmed.length() - ztrim) ;
			qualstring_trimmed = r.getBaseQualityString().substring(l, qualstring_trimmed.length() - ztrim) ;
		}
		else if(CLIP_BARCODE){
			subseq_trimmed = r.getReadString().substring(l) ;
			qualstring_trimmed = r.getBaseQualityString().substring(l) ;
		}
		else if(ztrim > 0){
			subseq_trimmed = r.getReadString().substring(0, subseq_trimmed.length() - ztrim) ;
			qualstring_trimmed = r.getBaseQualityString().substring(0, qualstring_trimmed.length() - ztrim) ;
		}
		

		FastqRecord trimmed = new FastqRecord(
				header,//seqHeaderPrefix
				subseq_trimmed, //seqLine,
				r.getBaseQualityHeader(),//qualHeaderPrefix,
				qualstring_trimmed // qualLine
				);
		return trimmed;
	}

	

	
	
	/**
	 * Checks if a given file name looks like a path
	 * @param fname
	 * @return
	 */
	protected boolean looksLikeAPath(String fname) {
		
		if(fname.contains("/")){
			File f = new File(fname);
			return f.getParentFile().exists() && f.getParentFile().isDirectory();
		}
		return false;
	}
	
}
