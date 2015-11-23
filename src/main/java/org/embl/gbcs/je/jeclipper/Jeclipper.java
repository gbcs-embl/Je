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
package org.embl.gbcs.je.jeclipper;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.lang3.StringUtils;
import org.embl.cg.utilitytools.utils.ExceptionUtil;
import org.embl.cg.utilitytools.utils.FileUtil;
import org.embl.gbcs.je.jemultiplexer.BarcodePosition;
import org.embl.gbcs.je.jemultiplexer.JemultiplexerFastqWriterFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Illumina;




/**
 * 
 * Jeclipper is a utility to clip molecular barcodes off FASTQ sequence.
 * Jeclipper is NOT a demultiplexer ie it writes a single result file (2 in case of paired end). 
 * Clipped sequences can be added to read headers 
 * 
 * @author girardot
 *
 */
@CommandLineProgramProperties(
		usage = "\tInput fastq file(s) can be in gzip compressed format (end in .gz). \n"
				+"\tBy default, output file(s) are gzipped and have names following the pattern :\n"
				+"\t\t'<inputfile_name>_clipped.<inputfile_extension>[.gz]'\n" 
				+"\tSee help for a detailled description of all options.\n"+
				"Example: \n"
				+"\tje clip F1=file.fastq.gz LEN=8 O=/path/to/resultdir/\n",
				usageShort = "je clip F1=file.fastq.gz LEN=8 O=/path/to/resultdir/", 
				programGroup = Illumina.class 
		)
public class Jeclipper extends CommandLineProgram {
	private static Logger log = LoggerFactory.getLogger(Jeclipper.class);

	
	
	/**
	 * Defaults
	 */
	protected static final BarcodePosition DEFAULT_BARCODE_READ_POS = BarcodePosition.BOTH;
	protected static final Integer DEFAULT_XTRIMLEN = 0;
	protected static final Integer DEFAULT_ZTRIMLEN = 0;
	protected static final boolean DEFAULT_ADD_BARCODE_TO_HEADER = true;
	protected static final String DEFAULT_READ_NAME_REPLACE_CHAR = ":";
	
	/**
	 * when DEFAULT_ADD_BARCODE_TO_HEADER is false, clipped barcodes are written to a specific file which 
	 * default name is DEFAULT_BARCODE_RESULTFILENAME
	 */
	protected static final String DEFAULT_BARCODE_RESULTFILENAME = "clipped_barcodes.txt"; 
	protected static final boolean DEFAULT_GZIP_OUTPUTS = true;
	
	
	
	/**
	 * test mode : tell the code to stop as soon as the doWork() is called i.e. allows to test the parsing in order to check 
	 * the status after parsing.
	 */
	protected boolean TEST_MODE_STOP_AFTER_PARSING = false;
	


	@Option(shortName="F1", optional = false, 
			printOrder=10, 
			doc="Input fastq file (optionally gzipped) for single end data, or first read in paired end data."
			)
	public File FASTQ_FILE1;

	@Option(shortName="F2", optional = true, 
			printOrder=20, 
			doc="Input fastq file (optionally gzipped) for the second read of paired end data."
			)
	public File FASTQ_FILE2 = null;

	/**
	 * Are we running SE or PE ? Option set at cmd line parsing time
	 */
	protected boolean RUNNING_PAIRED_END = false;

	
	@Option(shortName = "LEN", optional = false, 
			printOrder=30, 
			doc = "Length of the barcode sequences. When BARCODE_READ_POS == BOTH, two distinct" +
					" lengths can be provided using the syntax LEN=X:Z where X and Z are 2 integers representing" +
					" the barcode length for read_1 and read_2 respectively.\n"
			)
	public String BCLEN=null ; //BCLEN parsed to set bclen1 and bclen2
	protected Integer BCLEN_1 = 0;
	protected Integer BCLEN_2 = 0;
	
	@Option(shortName="BPOS", optional = true, 
			printOrder=40, 
			doc="Reads containing the sequence (i.e. UMIs) to clip:\n"
					+ "* READ_1 (beginning of read from FASTQ_FILE_1),\n"
					+ "* READ_2 (beginning of read from FASTQ_FILE_2),\n"
					+ "* BOTH (beginning of both reads).\n"
					+ "Automatically set to READ_1 in single end mode and BOTH in paired end mode. "
					+ "Actually not relevant for single end data\n"
			)
	public BarcodePosition BARCODE_READ_POS = DEFAULT_BARCODE_READ_POS ;


	@Option(shortName="ADD", optional = false, 
			printOrder=50, 
			doc="Should clipped UMIs be added to the read header (at the end); "
					+ "apply to both barcodes when BPOS=BOTH.\n"
					+"\t"+"If ADD=true, the string ':barcode' is added at the end of the read header with " 
					+"a ':' added only if current read header does not end with ':'.\n"
					+"\t"+"If both reads of the pair contains a UMI (i.e. BARCODE_READ_POS == BOTH), " 
					+"the UMI from the second read is also added to the read header. \n"
					+ "\t Else, the header of the read without UMI receives the UMI from the other read.\n"
					+"\t"+"For example :\n"
					+"\t"+"\t" +"'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:'\n"
					+"\t" +"becomes\n"
					+"\t"+"\t" +"'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:BARCODE'\n"
			)
	public boolean ADD_BARCODE_TO_HEADER  = DEFAULT_ADD_BARCODE_TO_HEADER;

	@Option(shortName = "RCHAR", optional = true,
			printOrder=70,
			doc="Replace spaces in read name/header using provided character.\n"
					+"This is needed when you need to retain ADDed barcode in read name/header " 
					+"during mapping as everything after space in read name is usually clipped in BAM files.\n" 
					+"For example, with RCHAR=':' :\n"
					+"\t"+"\t" +"'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 1:N:0:'\n"
					+"\t" +"becomes\n"
					+"\t"+"\t" +"'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965:1:N:0:BARCODE'\n"
			)
	public String READ_NAME_REPLACE_CHAR = DEFAULT_READ_NAME_REPLACE_CHAR;
	
	@Option(shortName="SAME_HEADERS", optional = true,
			printOrder=60,
			doc= "Makes sure headers of both reads of a pair are identical.\n"
					+ "Read name (or headers) will follow the pattern (for both reads of a pair) :\n"
					+ "\t'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 CLIPPED_SEQ_FROMREAD1:CLIPPED_SEQ_FROMREAD2 \n"
					+ "This option only makes sense in paired end mode and ADD=true.Some (if not all) mappers will indeed complain when "
					+ "read headers of a read pair are not identical.\n"
					+ "When SAME_HEADERS=FALSE and the RCHAR is used, read headers look like this :\n"
					+ "\t\tHISEQ:44:C6KC0ANXX:5:1101:1491:1994:1:N:0:TGGAGTAG\n"
					+ "\t\tHISEQ:44:C6KC0ANXX:5:1101:1491:1994:3:N:0:CGTTGTAT\n"
					+ "SAME_HEADERS=true will instead generates the following identical header for both reads :\n"
					+ "\tHISEQ:44:C6KC0ANXX:5:1101:1491:1994:TGGAGTAG:CGTTGTAT\n"
					+ "Note that we also clipped the useless '1:N:0' amd '3:N:0' as they also result in different headers\n"
					+ "Important : this option will force RCHAR=: UNLESS you specify RCHAR=null ; in which case a space will be preserved ie : \n"
					+ "\tHISEQ:44:C6KC0ANXX:5:1101:1491:1994 TAGAACAC:TGGAGTAG:CGTTGTAT"
					
			)
	public boolean ENSURE_IDENTICAL_HEADER_NAMES  = true;

	
	@Option(shortName = "XT", optional = true,
			printOrder=80,
			doc = "Optional extra number of base(s) to be trimmed right after the barcode. These extra bases are not added to read headers.\n" +
					"When running paired-end, two distinct values can be given using the syntax XT=X:Z where X and Z are 2 integers" +
					" to use for read_1 and read_2 respectively. Note that even when BPOS=READ_1 or BPOS=READ_2, a X:Y synthax can be given to " +
					"trim the read w/o barcode to end up with reads of identical length (note that this can also be operated using ZT). " +
					"If a unique value is given, e.g. XT=1, while running paired-end the following rule applies :\n " +
					"\t(1) BPOS=READ_1 or BPOS=READ_2, no trim is applied at the read w/o barcode \n"
					+ "\t(2) BPOS=BOTH, the value is used for both reads.\n"+
					"Note that XT=null is like XT=0.\n"
			)
	public String XTRIMLEN= DEFAULT_XTRIMLEN.toString() ; //must be further inspected to set XTRIMLEN_1 & XTRIMLEN_2
	public Integer XTRIMLEN_1=null ; // we do not set defaults as having null value here is part of the reasoning !!
	public Integer XTRIMLEN_2=null ; // we do not set defaults as having null value here is part of the reasoning !!

	@Option(shortName = "ZT", optional = true,
			printOrder=90,
			doc = "Optional extra number of bases to be trimmed from the read end i.e. 3' end. These extra bases are not added to read headers.\n" +
					"When running paired-end, two distinct values can be given here using the syntax ZT=X:Z where X and Z are 2 integers" +
					" to use for read_1 and read_2 respectively. Note that even when BPOS=READ_1 or BPOS=READ_2, a X:Y synthax can be given to " +
					"trim the read w/o barcode as to end up with reads of the same length (note that this can also be operated using XT). " +
					"Note that if a single value is passed, the value always applies to both reads in paired-end mode without further consideration.\n"
			)
	public String ZTRIMLEN=DEFAULT_ZTRIMLEN.toString() ;
	public Integer ZTRIMLEN_1=DEFAULT_ZTRIMLEN ;
	public Integer ZTRIMLEN_2=DEFAULT_ZTRIMLEN ;

	
	/*
	 * More general options below
	 */
	@Option(shortName = "OF1", optional = true,
			printOrder=100,
			doc="Optional result file name for Read_1 fastq file (i.e. F1) . Default name is built after input file name with the following pattern : <inputfile_name>_clipped.<inputfile_extension>[.gz]"
			+"\nCan either be a name (in which case the file will be created in the output dir) or a full path.\n"
			)
	public String RESULT_FILENAME_1 = null;
	File RESULT_FASTQ_FILE1 = null;  

	@Option(shortName = "OF2", optional = true,
			printOrder=110,
			doc="Optional result file name for Read_2 fastq file (i.e. F2) . Default name is built after input file name with the following pattern : <inputfile_name>_clipped.<inputfile_extension>[.gz]"
			+"\nCan either be a name (in which case the file will be created in the output dir) or a full path.\n"
			)
	public String RESULT_FILENAME_2 = null;
	File RESULT_FASTQ_FILE2 = null;  

	@Option(shortName = "BF", optional = true,
			printOrder=200,
			doc="Optional file name where to write clipped barcodes, default name is "+DEFAULT_BARCODE_RESULTFILENAME+".\n"
					+"Can either be a name (in which case the file will be created in the output dir) or a full path.\n"
					+"This file is automatically created if ADD=FALSE i.e. even if this option is not provided by user (and always created if this option is given).\n"
					+ "File format is tab delimited with : \n"
					+ "\t read header (col 1), barcode from read_1 (col 2), barcode quality from read_1 (col 2) and barcode + quality from read_2 (col 4 and 5 respectively) when relevant.\n"

			)
	public String BARCODE_RESULT_FILENAME = null;
	File barcodeFile = null; //final File object to use with METRICS_FILE_NAME, set in cmdline validation 

	
	@Option(shortName = "O", optional = true,
			printOrder=130, 
			doc="Where to write output file(s). By default, these are written in running directory.\n")
	public File OUTPUT_DIR = null;

	@Option(optional = true,
			printOrder=140,
			doc="Allows overwriting existing files.\n" 
			)
	public boolean FORCE = false;
	
	
	@Option(shortName="GZ", optional = true,
			printOrder=150,
			doc="Compress output s_l_t_barcode.txt files using gzip and append a .gz extension to the filenames.\n")
	public boolean GZIP_OUTPUTS = DEFAULT_GZIP_OUTPUTS;


	public Jeclipper(){
		super();
	}


	/**
	 * 
	 * @return null if command line is valid.  If command line is invalid, returns an array of error message
	 *         to be written to the appropriate place.
	 */
	@Override
	protected String[] customCommandLineValidation() {
		final ArrayList<String> messages = new ArrayList<String>();


		/*
		 * Check mandatory input files
		 */
		if(!FASTQ_FILE1.exists()){
			messages.add("FASTQ_FILE1 does not exist :"+FASTQ_FILE1.getAbsolutePath());
		}else if(!FASTQ_FILE1.canRead()){
			messages.add("Unsufficient rights to read FASTQ_FILE1 :"+FASTQ_FILE1.getAbsolutePath());
		}
		
		if(FASTQ_FILE2 != null){
			RUNNING_PAIRED_END = true;

			if(!FASTQ_FILE2.exists()){
				messages.add("FASTQ_FILE2 does not exist :"+FASTQ_FILE2.getAbsolutePath());
			}else if(!FASTQ_FILE2.canRead()){
				messages.add("Unsufficient rights to read FASTQ_FILE2 :"+FASTQ_FILE2.getAbsolutePath());
			}
			
			if(FASTQ_FILE2.getAbsolutePath().equals(FASTQ_FILE1.getAbsolutePath())){
				messages.add("FASTQ_FILE1 and FASTQ_FILE2 are the same file !");
			}
		}
		
		if(BARCODE_READ_POS == BarcodePosition.NONE){
			messages.add("Jeclipper cannot run with BPOS=NONE !");
		}
		
		if(!RUNNING_PAIRED_END){
			//force to the only option
			BARCODE_READ_POS=BarcodePosition.READ_1;
		}

		
		/*
		 * parse options that can have 'X:Y' syntax
		 * We first set these options blindly, later check is performed for those options that need it ie:
		 * - BCLEN* => make sure to have non 0 values at ends with barcodes and 0 values at ends w/o barcodes
		 * - XTRIMLEN => must be further inspected i.e. make sure we don t have a '1' for ends w/o 
		 * barcodes UNLESS explicitly set by user in the form X:Y 
		 */
		int [] arr = null;
		arr = setOption("BCLEN", BCLEN,  messages, true);
		if(arr.length == 2){ 
			//BCLEN was given in cmd line as it is null by default
			if(BARCODE_READ_POS != BarcodePosition.READ_2){
				BCLEN_1 = arr[0];
				log.debug("setting BCLEN_1 to "+BCLEN_1+" (BCLEN="+BCLEN);
			}
			if(BARCODE_READ_POS != BarcodePosition.READ_1){
				BCLEN_2 = arr[1];
				log.debug("setting BCLEN_2 to "+BCLEN_2+" (BCLEN was "+BCLEN);
			}
		}else{
			//BCLEN was NOT given in cmd line which is not possible since it is mandatory
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
		
		//check XTRIMLEN 
		if(RUNNING_PAIRED_END && BARCODE_READ_POS != BarcodePosition.BOTH){
			//set XTRIMLEN to 0 for the end w/o barcode UNLESS set by user in command line
			if(!XTRIMLEN.contains(":")){
				//user did not set 2 values => make it 0 
				if(BARCODE_READ_POS == BarcodePosition.READ_1)
					XTRIMLEN_2 = 0;
				else
					XTRIMLEN_1 = 0;
			}
		}
		
		if(READ_NAME_REPLACE_CHAR!=null){
			if(!org.embl.cg.utilitytools.utils.StringUtil.isValid(READ_NAME_REPLACE_CHAR)){
				messages.add("The option READ_NAME_REPLACE_CHAR is empty, for special characters like ';', you might need to escape them with a backslash.");
			}
			log.debug("RCHAR is '"+READ_NAME_REPLACE_CHAR+"'");
		}
		
		
		//only continue if all seems ok 
		if(messages.isEmpty()){

			if(OUTPUT_DIR == null){
				OUTPUT_DIR = new File(System.getProperty("user.dir"));
			}
			/*
			 * check output dir
			 */
			if(!OUTPUT_DIR.exists()){
				log.info("Attempting to create output directory : "+OUTPUT_DIR.getAbsolutePath());
				try{
					FileUtil.checkWritableDir(OUTPUT_DIR, true);
				}catch(Exception e){
					//creation failed 
					messages.add("Creation of result directory failed! Please check you have sufficient privileges to create this dir.");
				}
			}


			/*
			 * Does user want the clipped barcodes in the read headers or in a separate file ?  
			 * We create one automatically if ADD_BARCODE_TO_HEADER is false
			 */
			if(!ADD_BARCODE_TO_HEADER || BARCODE_RESULT_FILENAME!=null){
				String fname = DEFAULT_BARCODE_RESULTFILENAME;
				if(BARCODE_RESULT_FILENAME != null) 
					fname = BARCODE_RESULT_FILENAME;
				if(GZIP_OUTPUTS){
					fname += fname.endsWith("gz") ? "":".gz";
				}

				barcodeFile = setFile(fname);
			}

			if(!StringUtils.isBlank(RESULT_FILENAME_1)){
				RESULT_FASTQ_FILE1 = setFile(this.RESULT_FILENAME_1);
			}else{
				String defaultname = FileUtil.removeExtension(FASTQ_FILE1.getName()) + "_clipped."+ FileUtil.getExtensionZipOrGZIPSafe(FASTQ_FILE1.getName()) ;
				if(GZIP_OUTPUTS){
					defaultname += defaultname.endsWith("gz") ? "":".gz";
				}
				RESULT_FASTQ_FILE1 = setFile( defaultname );
			}

			//check the file does not exist yet
			if(RESULT_FASTQ_FILE1.exists()){
				if(!FORCE){
					String mess = "\nOutput file already exists : " + RESULT_FASTQ_FILE1.getAbsolutePath() +
							"\nPlease use FORCE=true to overwrite existing file(s).";
					messages.add(mess);
				}else if(!RESULT_FASTQ_FILE1.delete()){
					String mess = "Failed to remove existing output file (FORCE=true) : " + RESULT_FASTQ_FILE1.getAbsolutePath() +
							"\nPlease check you have enought privileges to perform this operation.";
					messages.add(mess);
				}
			}

			if(RUNNING_PAIRED_END){
				if(RESULT_FILENAME_2 != null){
					RESULT_FASTQ_FILE2 = setFile(this.RESULT_FILENAME_2);
				}else{
					String defaultname = FileUtil.removeExtension(FASTQ_FILE2.getName()) + "_clipped."+ FileUtil.getExtensionZipOrGZIPSafe(FASTQ_FILE2.getName()) ;
					if(GZIP_OUTPUTS){
						defaultname += defaultname.endsWith("gz") ? "":".gz";
					}
					RESULT_FASTQ_FILE2 = setFile( defaultname );
				}

				//check the file does not exist yet
				if(RESULT_FASTQ_FILE2.exists()){
					if(!FORCE){
						String mess = "\nOutput file already exists : " + RESULT_FASTQ_FILE2.getAbsolutePath() +
								"\nPlease use FORCE=true to overwrite existing file(s).";
						messages.add(mess);
					}else if(!RESULT_FASTQ_FILE2.delete()){
						String mess = "Failed to remove existing output file (FORCE=true) : " + RESULT_FASTQ_FILE2.getAbsolutePath() +
								"\nPlease check you have enought privileges to perform this operation.";
						messages.add(mess);
					}
				}
			} 
		}
		
		if(messages.isEmpty()){
			return null;
		}
		return messages.toArray(new String[messages.size()]);
	}



	
	@Override
	protected int doWork() {
		
		log.info("[options valid, start processing]");
		
		if(TEST_MODE_STOP_AFTER_PARSING)
			return 0;
		
	
		/*
		 * prepare FASTQ file writers
		 */
		JemultiplexerFastqWriterFactory fastqFactory = new JemultiplexerFastqWriterFactory();
		fastqFactory.setUseAsyncIo(false);

		/*
		 * init writers for output files
		 * at this point, the file path has been already set and the overwrite status checked (if needed)
		 * we can blindly init writers 
		 */
		FastqWriter fw1 = fastqFactory.newWriter(RESULT_FASTQ_FILE1, GZIP_OUTPUTS, CREATE_MD5_FILE);
		log.debug("[OK] writer on result file "+RESULT_FASTQ_FILE1.getAbsolutePath());
		FastqWriter fw2 = null;
		if(RUNNING_PAIRED_END){
			fw2 = fastqFactory.newWriter(RESULT_FASTQ_FILE2, GZIP_OUTPUTS, CREATE_MD5_FILE);
			log.debug("[OK] writer on result file "+RESULT_FASTQ_FILE2.getAbsolutePath());
		}	
		
		
		
		/*
		 * Should we init a barcode file ?
		 */
		PrintWriter barcodeW = null;
		if( ! ADD_BARCODE_TO_HEADER ){
			try{
				barcodeW = GZIP_OUTPUTS ? new PrintWriter( new GZIPOutputStream(new FileOutputStream(barcodeFile)) ) : new PrintWriter(barcodeFile);
			} catch (Exception e) {
				log.error(ExceptionUtil.getStackTrace(e));
				throw new RuntimeException(e);
			}
			
			String[] headers = new String []{"READ_HEADER", "BC_SEQ_READ1", "BC_QUAL_READ1"};
			if(RUNNING_PAIRED_END)
				headers = new String []{"READ_HEADER", "BC_SEQ_READ1", "BC_QUAL_READ1", "BC_SEQ_READ2", "BC_QUAL_READ2"};

			barcodeW.println(org.embl.cg.utilitytools.utils.StringUtil.mergeArray(headers,"\t"));
		}

		FastqReader fqr1 = new FastqReader(FASTQ_FILE1);
		Iterator<FastqRecord> it1 = fqr1.iterator();
		FastqReader fqr2 = null;
		Iterator<FastqRecord> it2 = null;

		
		if(RUNNING_PAIRED_END){
			fqr2 = new FastqReader(FASTQ_FILE2);
			it2 = fqr2.iterator();
		}

		
		/*
		 * Read files and write
		 */
		int cnt = 0; //read pair counter
		while(it1.hasNext()){
			FastqRecord r1 = it1.next();
			FastqRecord r2 = null; 
			if(RUNNING_PAIRED_END){
				r2 = it2.next();
			}
			cnt++;
			if(cnt % (1000*1000) == 0){
				log.info(cnt+" reads processed");

			}

			/*
			 * Operate barcode EXTRACTION 
			 */
			String s1 = null;
			String q1 = null;
			String s2 = null;
			String q2 = null;


			//do we have a barcode at read 1 : yes if SE or PE with BC position != read2 
			if(!RUNNING_PAIRED_END || BARCODE_READ_POS != BarcodePosition.READ_2){
				//read 1 has a barcode, extract it
				s1 = r1.getReadString().substring(0, BCLEN_1).toUpperCase();
				q1 = r1.getBaseQualityString().substring(0, BCLEN_1);
				log.debug("Got bc sequence from read1 : "+s1);
			}
			//do we have a barcode at read 2 ? => yes is PE with BC position != read1 
			if(RUNNING_PAIRED_END && BARCODE_READ_POS != BarcodePosition.READ_1){
				s2 = r2.getReadString().substring(0, BCLEN_2);
				q2 = r2.getBaseQualityString().substring(0, BCLEN_2);
				log.debug("Got bc sequence from read2 : "+s2);
			}

			
			/*
			 * write to barcode file if requested
			 */
			if(barcodeW!=null){
				String [] values = null;
				if(RUNNING_PAIRED_END){
					values = new String[]{
							r1.getReadHeader().split("\\s+")[1], 
							s1, 
							q1,
							s2,
							q2
					};
				}else{
					values = new String[]{
							r1.getReadHeader().split("\\s+")[1], 
							s1, 
							q1
					};
				}
				barcodeW.println(org.embl.cg.utilitytools.utils.StringUtil.mergeArray(values,"\t"));
			}
			
			/*
			 * write reads to fastq files
			 */
			String codeForRead1 = s1;
			String codeForRead2 = s2;
			
			// if ENSURE_IDENTICAL_HEADER_NAMES is true we need to make sure the codes are identical for each read
			if(RUNNING_PAIRED_END && ADD_BARCODE_TO_HEADER && ENSURE_IDENTICAL_HEADER_NAMES){
				if(BARCODE_READ_POS == BarcodePosition.BOTH){
					// merge both code and init both BCLEN_1 and BCLEN_2 to it
					String merged = s1 + ":" + s2;
					codeForRead1 = merged;
					codeForRead2 = merged;
				}else{
					// make sure both BCLEN_1 and BCLEN_2 have the value
					if(BARCODE_READ_POS == BarcodePosition.READ_1){
						codeForRead2 = s1;
					}else{
						codeForRead1 = s2;
					}
				}
				 
			}
			
			
			write(r1, BCLEN_1, XTRIMLEN_1, ZTRIMLEN_1, fw1, codeForRead1, ADD_BARCODE_TO_HEADER);
			if(RUNNING_PAIRED_END){
				write(r2, BCLEN_2 , XTRIMLEN_2 , ZTRIMLEN_2, fw2, codeForRead2, ADD_BARCODE_TO_HEADER );
			}

			
		}
			
		//close all writers
		fw1.close();
		fqr1.close();
		if(fw2!=null){
			fw2.close();
			fqr2.close();
		}
		if(barcodeW!=null) barcodeW.close(); 
		
		//all fine
		return 0;

	}


	/**
	 * @param r the read to write
	 * @param bcLen the barcode length
	 * @param xtrim how many extra base should be clipped after the barcode ; 0 or more
	 * @param writer the writer
	 * @param barcode the barcode sequence extracted from this read. Note that the caller is responsible to give the 
	 *  sequence as extracted from the read or the barcode sequence after matching the extracted sequence against list of expected barcode
	 * @param addBarcodeToHeader whether the barcode should be added to the read header. If true, the string ':barcode' is added 
	 * with the starting ':' added only if current read header does not end with ':'. For example :
	 * '@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:'
	 * becomes
	 * '@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:BARCODE'
	 * @return the trimmed {@link FastqRecord} ; mostly for testing purposes
	 */
	protected FastqRecord write(final FastqRecord r, int bcLen, int xtrim, int ztrim, final FastqWriter writer, String barcode, Boolean addBarcodeToHeader) {
		
		log.debug("Writing read [bcLen="+bcLen+", xtrim="+xtrim+", ztrim="+ztrim+", barcode="+(barcode==null?"NULL":barcode)+", addBarcodeToHeader="+addBarcodeToHeader.toString()+"]");
	
		int l = bcLen + xtrim ;
		
		String header = r.getReadHeader();
		
		if(RUNNING_PAIRED_END && ENSURE_IDENTICAL_HEADER_NAMES && addBarcodeToHeader){
			//clip off header after space
			header = r.getReadHeader().split("\\s+")[0]+" "+(barcode!=null?barcode:"");
			if(READ_NAME_REPLACE_CHAR != null){
				header = header.replaceAll(" ", READ_NAME_REPLACE_CHAR);
			}
		}
		else if(addBarcodeToHeader){
			header+= (header.endsWith(":")? "":":")+(barcode!=null?barcode:"");
			if(READ_NAME_REPLACE_CHAR != null){
				header = header.replaceAll(" ", READ_NAME_REPLACE_CHAR);
			}
		}
		log.debug("will clip "+l+" bases of "+r.getReadString());

		String subseq_trimmed = r.getReadString();
		String qualstring_trimmed = r.getBaseQualityString();
		if(ztrim > 0){
			subseq_trimmed = r.getReadString().substring(l, subseq_trimmed.length() - ztrim) ;
			qualstring_trimmed = r.getBaseQualityString().substring(l, qualstring_trimmed.length() - ztrim) ;
		}
		else{
			subseq_trimmed = r.getReadString().substring(l) ;
			qualstring_trimmed = r.getBaseQualityString().substring(l) ;
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
	 * Deals with integer options that can take X:Y values ; also check that given values are >= 0 
	 * @param optionName name of the option (for error messages)
	 * @param optionValue as given by user ; if blank returned array is empty
	 * @param messages the error message list to use in case of error
	 * @param failIfPEwithSingleBarcode 
	 * @return an int array with the two parsed values for read_1 and read_2. 
	 * Array is empty is nothing should be made ie in such case one should check if messages size increased
	 */
	private int[] setOption(String optionName, String optionValue, ArrayList<String> messages, boolean failIfPEwithSingleBarcode) {
		
		int[] arr = new int[2];
		
		if(htsjdk.samtools.util.StringUtil.isBlank(optionValue))
			return arr;
		
		String [] values = new String[2]; 
		if(optionValue.contains(":")){
			 if(!RUNNING_PAIRED_END || (failIfPEwithSingleBarcode && BARCODE_READ_POS!=BarcodePosition.BOTH)){
				messages.add("Invalid value used for "+optionName+" : "+optionValue +". The X:Z synthax can only be used for paired-end data with BPOS=BOTH");
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
		
		return new int []{optionRead1, optionRead2};
	}


	private File setFile(String fnameOrPath) {
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
	 * @param args
	 */
	public static void main(final String[] argv) {
		new Jeclipper().instanceMainWithExit(argv); //eventually execute the doWork() method
	}
	
	

}
