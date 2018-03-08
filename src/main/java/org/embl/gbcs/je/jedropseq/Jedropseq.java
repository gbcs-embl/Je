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
package org.embl.gbcs.je.jedropseq;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.lang3.StringUtils;
import org.embl.cg.utilitytools.utils.FileUtil;
import org.embl.gbcs.je.JemultiplexerFastqWriterFactory;
import org.embl.gbcs.je.ReadLayoutConsumer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Illumina;




/**
 * 
 * Jedropseq is a utility to clip cell barcode and molecular barcode off reads gained from a DROP-seq protocol.
 * Jedropseq only accepts paired end data as cell barcodes and UMI sequences are clipped from the read1 and added to read name of read2. 
 * Jedropseq is NOT a demultiplexer ie it writes a single result file (2 in case of paired end). 
 * 
 * 
 * @author girardot
 *
 */
@CommandLineProgramProperties(
		usage ="Reformat Drop-seq files into a single fastq file. DROP-seq produces 2 fastq files : one (F1) contains the cell barcode (usually 12 bp, LEN)"+
				" followed by a UMI (usually 8bp, ULEN) while the second (F2) contains the RNA sequence. The output file is similar to F2 but holds the "+
				"parsed barcode and UMI in read names.\n"+ 
				"Input fastq file(s) can be in gzip compressed format (end in .gz). See help for a detailled description of all options.\n" +
				"\n"+
				"Example: \n"+
				"\t"+"je dropseq F1=file.fastq.gz F2=file.fastq.gz LEN=12 ULEN=8 O=/path/to/resultdir/out.fastq.gz\n",
				usageShort = "je dropseq F1=file_1.fastq.gz F2=file_2.fastq.gz LEN=12 ULEN=8 O=/path/to/resultdir/file.fastq.gz", 
				programGroup = Illumina.class 
		)
public class Jedropseq extends CommandLineProgram {
	private static Logger log = LoggerFactory.getLogger(Jedropseq.class);

	
	
	/**
	 * Defaults
	 */	
	protected static final String DEFAULT_READ_NAME_REPLACE_CHAR = ":";
	protected static final boolean DEFAULT_GZIP_OUTPUTS = true;
	
	
	
	/**
	 * test mode : tell the code to stop as soon as the doWork() is called i.e. allows to test the parsing in order to check 
	 * the status after parsing.
	 */
	public boolean TEST_MODE_STOP_AFTER_PARSING = false;
	


	@Option(shortName="F1", optional = false, 
			printOrder=10, 
			doc="Input fastq file (optionally gzipped) for first read. This read contains the cell barcode followed by the UMI"
			)
	public File FASTQ_FILE1;

	@Option(shortName="F2", optional = false, 
			printOrder=20, 
			doc="Input fastq file (optionally gzipped) for the second read."
			)
	public File FASTQ_FILE2 = null;

	
	@Option(shortName = "LEN", optional = false, 
			printOrder=30, 
			doc = "Length of the cell barcode sequence in read 1.\n"
			)
	public Integer BCLEN=null ; 
	 
	
	@Option(shortName = "ULEN", optional = false, 
			printOrder=40, 
			doc = "Length of the UMI sequence in read 1 found right after the cell barcode.\n"
			)
	public Integer UCLEN=null ; 
	
	
	@Option(shortName="WQ", optional = true,
			printOrder=35,
			doc="Should quality string of barcode and UMI also be injected in read names.\n"+
			"If true, the quality string is translated into 2 digits number and a e.g. UMI will look like\n"+
					"\t"+" '...:ATGCAT333423212322:...' instead of '...:ATGCAT:...'\n"+
			"This option is particularly useful with the retag module that knows how to extract quality numbers into BAM tags."
				)
	public boolean WITH_QUALITY_IN_READNAME = false;

	
	@Option(shortName = "N", optional = true, 
			printOrder=50, 
			doc = "Maximum number of N's in the cell barcode sequence. If the cell barcode has this number or more N in the sequence, the read is ignored.\n"
			)
	public Integer MAX_N = 6 ; 
	
	
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
	
		
	/*
	 * More general options below
	 */
	@Option(shortName = "O", optional = false,
			printOrder=100,
			doc="Result file name with headers modified.\n"
			+"\nCan either be a name (in which case the file will be created in the output dir) or a full path.\n"
			)
	public String RESULT_FILENAME_1 = null;
	File RESULT_FASTQ_FILE1 = null;  

	
	@Option(optional = true,
			printOrder=140,
			doc="Allows overwriting existing files.\n" 
			)
	public boolean FORCE = false;
	
	
	@Option(shortName="GZ", optional = true,
			printOrder=150,
			doc="Compress output using gzip and append a .gz extension to the result filename if necessary.\n")
	public boolean GZIP_OUTPUTS = DEFAULT_GZIP_OUTPUTS;


	public Jedropseq(){
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
		
		if(!FASTQ_FILE2.exists()){
			messages.add("FASTQ_FILE2 does not exist :"+FASTQ_FILE2.getAbsolutePath());
		}else if(!FASTQ_FILE2.canRead()){
			messages.add("Unsufficient rights to read FASTQ_FILE2 :"+FASTQ_FILE2.getAbsolutePath());
		}
		
		if(FASTQ_FILE2.getAbsolutePath().equals(FASTQ_FILE1.getAbsolutePath())){
			messages.add("FASTQ_FILE1 and FASTQ_FILE2 are the same file !");
		}
		
	
		if(READ_NAME_REPLACE_CHAR!=null){
			if(!org.embl.cg.utilitytools.utils.StringUtil.isValid(READ_NAME_REPLACE_CHAR)){
				messages.add("The option READ_NAME_REPLACE_CHAR is empty, for special characters like ';', you might need to escape them with a backslash.");
			}
			log.debug("RCHAR is '"+READ_NAME_REPLACE_CHAR+"'");
		}
		
		
		//only continue if all seems ok 
		if(messages.isEmpty()){



			if(!StringUtils.isBlank(RESULT_FILENAME_1)){
				RESULT_FASTQ_FILE1 = setFile(this.RESULT_FILENAME_1);
			}else{
				String defaultname = FileUtil.removeExtension(FASTQ_FILE2.getName()) + "_tagged."+ FileUtil.getExtensionZipOrGZIPSafe(FASTQ_FILE2.getName()) ;
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
		FastqWriter writer = fastqFactory.newWriter(RESULT_FASTQ_FILE1, GZIP_OUTPUTS, CREATE_MD5_FILE);
		log.debug("[OK] writer on result file "+RESULT_FASTQ_FILE1.getAbsolutePath());
		

		FastqReader fqr1 = new FastqReader(FASTQ_FILE1);
		Iterator<FastqRecord> it1 = fqr1.iterator();
		FastqReader fqr2 = new FastqReader(FASTQ_FILE2);
		Iterator<FastqRecord> it2 = fqr2.iterator();		
		
		/*
		 * Read files and write
		 */
		int cnt = 0; //read pair counter
		int cntTooManyNs = 0;
		int written = 0;
		while(it1.hasNext()){
			FastqRecord r1 = it1.next();
			FastqRecord r2 = it2.next();

			cnt++;
			if(cnt % (1000*1000) == 0){
				log.info(cnt+" read pairs processed");

			}

			/*
			 * Operate barcode EXTRACTION 
			 */
			String cellBarcodeSeq = r1.getReadString().substring(0, BCLEN).toUpperCase();
			String cellBarcodeQual = r1.getBaseQualityString().substring(0, BCLEN);
			String umiSeq = r1.getReadString().substring(BCLEN, BCLEN+UCLEN).toUpperCase();
			String umiQual = r1.getBaseQualityString().substring(BCLEN, BCLEN+UCLEN);
			
			
			if( StringUtils.countMatches(cellBarcodeSeq, "N") >= MAX_N ){
				cntTooManyNs++;
				continue;
			}
			
			if(WITH_QUALITY_IN_READNAME) {
				//add the converted quality
				cellBarcodeSeq += ReadLayoutConsumer.qualityToNumberString( SAMUtils.fastqToPhred(cellBarcodeQual) );
				umiSeq += ReadLayoutConsumer.qualityToNumberString( SAMUtils.fastqToPhred(umiQual) );
			}
			
			//write read2 with augmented header
			String header = r2.getReadName().split("\\s")[0];
			header = header + this.READ_NAME_REPLACE_CHAR + cellBarcodeSeq + this.READ_NAME_REPLACE_CHAR + umiSeq;
					
			FastqRecord modifed = new FastqRecord(
				header,//seqHeaderPrefix
				r2.getReadString(), //seqLine,
				r2.getBaseQualityHeader(),//qualHeaderPrefix,
				r2.getBaseQualityString() // qualLine
				);

			if(writer != null){ //mostly here for testing ie to be able to call method w/o writer
				writer.write(modifed);
				written++;
			}
		}

		//close all writers
		writer.close();
		fqr1.close();		
		fqr2.close();
		
		//all fine
		log.info("Read pairs processed : "+"\t"+cnt);
		log.info("Reads written : "+"\t"+written);
		log.info("Reads ignored (too many Ns in cell barcode) : "+"\t"+cntTooManyNs);
		return 0;

	}

	private File setFile(String fnameOrPath) {
		if(fnameOrPath.contains("/")){
			//it is a path, use file path as provided
			return new File(fnameOrPath);
		}
		else{
			//it is just a name
			return new File(fnameOrPath);
		}

	}

	
	/**
	 * @param args
	 */
	public static void main(final String[] argv) {
		new Jedropseq().instanceMainWithExit(argv); //eventually execute the doWork() method
	}
	
	

}
