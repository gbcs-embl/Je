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

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.embl.cg.utilitytools.utils.ExceptionUtil;
import org.embl.cg.utilitytools.utils.FileUtil;
import org.embl.cg.utilitytools.utils.StringUtil;
import org.embl.cg.utilitytools.utils.parser.csv.CSVLine;
import org.embl.cg.utilitytools.utils.parser.csv.InvalidHeaderException;
import org.embl.gbcs.embase.api.EmBaseDatabaseFacade;
import org.embl.gbcs.embase.api.Queries;
import org.embl.gbcs.embase.api.exception.EmBASEConnectionException;
import org.embl.gbcs.embase.api.model.Item;
import org.embl.gbcs.embase.api.model.NGSLibrary;
import org.embl.gbcs.je.ApplicationConfiguration;
import org.embl.gbcs.je.BarcodeMatch;
import org.embl.gbcs.je.Je;
import org.embl.gbcs.je.JeUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Illumina;




@CommandLineProgramProperties(
		usage = "\tFastq files demultiplexer for in-line barcoded Illumina Fastq files.\n"
				+ "\tInput fastq files can be in gzip compressed (end in .gz). \n" 
				+ "\tBy default output files are gzipped and have names following the pattern \n"
				+ "\t'<samplename>_<barcode[-barcode2...-barcodeN]>_[1|2].txt[.gz]' unless you gave file\n"
				+ "\t names to use within the barcode description file.\n" 
				+"Example:\n"
				+"\t je demultiplex F1=fastq_1.txt.gz BF=barcodes.bs O=/path/to/jemultiplexer-results/ ",
				usageShort = "Demultiplex fastq files with reads which first bases hold a barcode", 
				programGroup = Illumina.class 
		)
public class Jemultiplexer extends AbstractJemultiplexer {
	private static Logger log = LoggerFactory.getLogger(Jemultiplexer.class);

	
	
	
	@Option(shortName="BPOS", optional = true,
			printOrder=41,
			doc="For paired-end data, where to expect the barcode(s) :\n"
					+ "\t * READ_1 (beginning of read from FASTQ_FILE_1), \n"
					+" \t * READ_2 (beginning of read from FASTQ_FILE_2),\n"
					+ "\t * BOTH (beginning of both reads).\n"+
					"Automatically set to READ_1 in single end mode.\n"
			)
	public BarcodePosition BARCODE_READ_POS = DEFAULT_BARCODE_READ_POS ;

	
	@Option(shortName = "LEN", optional = true,
			printOrder=43,
			doc = "Length of the barcode sequences, optional. Taken from barcode file when not given.\n" +
					"In situations where BARCODE_READ_POS == BOTH AND REDUNDANT_BARCODES=false, two distinct" +
					" length can be provided using the syntax LEN=X:Z where X and Z are 2 integers representing" +
					" the barcode length for read_1 and read_2 respectively.\n"
			)
	public String BCLEN=null ; //BCLEN parsed to set bclen1 and bclen2
	
	
	@Option(shortName="BM", optional = true,
			printOrder=45,
			doc="Indicates which barcode(s) should be used for sample lookup\n" +
					"Automatically set to READ_1 in single end mode. \n" +
					"For paired-end data and when BARCODE_READ_POS == BOTH, which barcode should be used to resolve sample :\n"
					+"\t"+"- use BM=READ_1 (beginning of read from FASTQ_FILE_1) if only this read should be used for sample matching,\n"
					+"\t"+"- use BM=READ_2 (beginning of read from FASTQ_FILE_2) if only this read should be used for sample matching,\n"
					+"\t"+"- use BM=BOTH (beginning of both reads) if both should be used \n"
					+ "When BM=BOTH, the behaviour is different based on the value of REDUNDANT_BARCODES : \n"
					+"\t"+"\t"+"If REDUNDANT_BARCODES=true, the two barcodes are considered to map to the same sample and 'Je demultiplex'" +
					" uses the two barcodes according to the STRICT value.\n" +
					"\t"+"\t"+"If REDUNDANT_BARCODES=false, the barcode file should map a couple of barcode to each sample (e.g. sample1 => AGAGTG:TTGATA)" +
					" and 'Je demultiplex' needs both barcodes to find the relevant sample. Note that this is the only situation in which all barcode matching options" +
					" (MM, MMD, Q) accept different values for both barcodes in the form X:Z where X and Z are 2 integers.\n"
					
			)
	public BarcodePosition BARCODE_FOR_SAMPLE_MATCHING = BarcodePosition.BOTH;

	
	
	@Option(shortName="BRED", optional = true,
			printOrder=47,
			doc="This option only applies for paired-end data with BARCODE_READ_POS set to 'BOTH'\n" +
					"Indicates if both read's barcodes encode redundant information or if barcodes are supposed to "
					+ "be identical at both ends (or to resolve to the same sample when a pool of barcodes is used per sample).\n "+
					"\tWhen REDUNDANT_BARCODES=false, the 2 barcodes potentially encode\n" +
					" different information. For example, only one of the barcodes encodes the sample identity while \n" +
					"the second barcode might be a random barcode (UMI) to tell apart PCR artefacts from real duplicates.\n" +
					"Another example is when both barcodes should be used in a combined fashion to resolve the sample.\n" +
					"In the first example, you should use BPOS=BOTH BRED=false BM=READ_1.\n" +
					"In the second example, you should have BPOS=BOTH BRED=false BM=BOTH. \n"
					+ "Note that with BPOS=BOTH BRED=true BM=BOTH, the behavior would be different as 'demultiplex' would then check the STRICT option to perform sample resolution.\n" +
					"Importantly, when BARCODE_READ_POS (BPOS) == BOTH AND REDUNDANT_BARCODES=false, BLEN, barcode matching options" +
					" (MM, MMD, Q) and read trimming/clipping options (XT, ZT) accept different values for both barcodes in " +
					"the form X:Z where X and Z are 2 integers.\n"
					
			)
	public boolean REDUNDANT_BARCODES = DEFAULT_REDUNDANT_BARCODES;

	
	
	@Option(shortName="ADD", optional = true,
			printOrder=112,
			doc="Add barcode at the end of the read header. Apply to both barcodes when BPOS=BOTH.\n"
					+"\t"+"If true, the string ':barcode' is added at the end of the read header with " 
					+"a ':' added only if current read header does not end with ':'.\n"
					+"\t"+"If both reads of the pair have a barcode (i.e. BARCODE_READ_POS == BOTH), then" 
					+"the second read also has its own matched barcode written. Else, the read without a barcode "
					+"receives the barcode from the barcoded read.\n"
					+"\t"+"For example :\n"
					+"\t"+"\t" +"'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:'\n"
					+"\t" +"becomes\n"
					+"\t"+"\t" +"'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:BARCODE'\n"
					+"\n\tWhen barcodes containing random positions, i.e. 'N', (for example like " +
					"\tin the iCLIP protocol) or are UMIs, the added sequence is the sequence clipped from the read and NOT the matched barcode.\n"
			)
	public boolean ADD_BARCODE_TO_HEADER  = DEFAULT_ADD_BARCODE_TO_HEADER;

	
	@Option(shortName="SAME_HEADERS", optional = true,
			printOrder=115,
			doc= "Makes sure that headers of both reads of a pair are identical, using the following read header pattern (for both reads of a pair) :\n"
					+ "\t\t'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 SAMPLEBARCODE_READ1:SAMPLEBARCODE_READ2(if applicable)':CLIPPED_SEQ_FROMREAD1:CLIPPED_SEQ_FROMREAD2"
					+ "This option only makes sense in paired end mode and ADD=true. "
					+ "Some (if not all) mappers will indeed complain when"
					+ " the read headers are not identical. When molecular barcodes are present in reads (either as"
					+ " additional barcodes or as degenerate barcodes ie with 'N') and the RCHAR is used, you will end"
					+ " with (problematic) read headers like this :\n"
					+ "\t\tHISEQ:44:C6KC0ANXX:5:1101:1491:1994:1:N:0:TAGAACAC:TGGAGTAG\n"
					+ "\t\tHISEQ:44:C6KC0ANXX:5:1101:1491:1994:3:N:0:TAGAACAC:CGTTGTAT\n"
					+ "SAME_HEADERS=true will instead generates the following identical header for both reads :\n"
					+ "\t\tHISEQ:44:C6KC0ANXX:5:1101:1491:1994:TAGAACAC:TGGAGTAG:CGTTGTAT\n"
					+ "\tNote that we also clipped the useless '1:N:0' and '3:N:0' has they will also result in generating different headers.\n"
					+ "\t\t Important : this option will force RCHAR=: UNLESS you specify RCHAR=null ; in which case a space will be preserved ie : \n"
					+ "\t\tHISEQ:44:C6KC0ANXX:5:1101:1491:1994 TAGAACAC:TGGAGTAG:CGTTGTAT\n"
					
			)
	public boolean ENSURE_IDENTICAL_HEADER_NAMES  = DEFAULT_ENSURE_IDENTICAL_HEADER_NAMES;

	
	
	
	public Jemultiplexer(){
		super();
		USE_SAMPLE_INDEX_FILES = false;
	}

	
	@Override
	protected String[] customCommandLineValidation() {
		
		/*
		 * run custom validation that does NOT depend on other validation 
		 *  
		 */
		ArrayList<String> messages = new ArrayList<String>();
		
		
		/*
		 * validate FASTQ files
		 */
		messages.addAll( validateFASTQFiles());
		
		if(FASTQ_FILE2 != null)
			RUNNING_PAIRED_END = true;

		
		if(ENSURE_IDENTICAL_HEADER_NAMES && RUNNING_PAIRED_END ){
			// check RCHAR : if this is not explicitely given in cmd line, then turn it to :
			String _cmdline = Arrays.toString(this.getCommandLineParser().getArgv());
			log.debug(Arrays.toString(this.getCommandLineParser().getArgv()));
			if(!_cmdline.contains("READ_NAME_REPLACE_CHAR") && !_cmdline.contains("RCHAR")){
				READ_NAME_REPLACE_CHAR = ":";
			}else{
				if(READ_NAME_REPLACE_CHAR == null || READ_NAME_REPLACE_CHAR.equalsIgnoreCase("NULL"))
					READ_NAME_REPLACE_CHAR = " ";
			}
		}
		
		if(READ_NAME_REPLACE_CHAR!=null){
			if(!StringUtil.isValid(READ_NAME_REPLACE_CHAR)){ // this is not failing for a space, this is wanted
				messages.add("The option READ_NAME_REPLACE_CHAR is empty, for special characters like ';', you might need to escape them with a backslash.");
			}
			log.debug("RCHAR is "+READ_NAME_REPLACE_CHAR);
		}
		
		return customCommonCommandLineValidation(messages);
		
	}
	
	
	@Override
	protected int doDemultiplexingWork(PrintWriter diagPW) {
		//should we write a barcode match result file ?
		boolean writingBarcodeDiagFile = (this.bcDiagFile != null);

		
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
		int unassigned = 0; //count how many pairs end up unassigned 
		int assigned = 0; //count how many pairs end up with a sample assignment
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
			BarcodeMatch bc1 = null;
			String s2 = null;
			String q2 = null;
			BarcodeMatch bc2 = null;
			

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
			 * Operate barcode MATCHING
			 */
			//do we need to match barcode 1? 
			//yes if (1) single end or (2) there is a barcode in READ1 position AND it is used for sample matching
			if(!RUNNING_PAIRED_END || (BARCODE_READ_POS != BarcodePosition.READ_2 && BARCODE_FOR_SAMPLE_MATCHING != BarcodePosition.READ_2 )){
				bc1 = findBarcode(BarcodePosition.READ_1, s1, q1, bc2sample, sample2bcSet);
				log.debug("bc1 match : "+ (bc1 == null?"NULL":bc1.toString()));
			}
			
			//do we need to match barcode 2? 
			// yes if (1) running paired-end AND there is a barcode in READ2 position AND it is used for sample matching
			// 
			if(	(RUNNING_PAIRED_END && BARCODE_READ_POS != BarcodePosition.READ_1 && BARCODE_FOR_SAMPLE_MATCHING != BarcodePosition.READ_1)){
				bc2 = findBarcode(BarcodePosition.READ_2, s2, q2, bc2sample, sample2bcSet);
				log.debug("bc2 match : "+ (bc2 == null?"NULL":bc2.toString()));
			}
			
						
			/*
			 * Special case
			 *  we have 2 barcode but only one is used for sample lookup, the second one being most likely a random tag 
			 * for PCR duplicate assessment => we still need to fill in bc1 or bc2
			 * 
			 * Indeed, if index files are used then bc are set correctly and we can t be in a weird situation (also
			 * remember that using index files is exclusive ie no barcode for sample matching can be in the reads
			 *  anymore, only molecular indices)
			 * 
			 */
			if(RUNNING_PAIRED_END && BARCODE_READ_POS == BarcodePosition.BOTH && BARCODE_FOR_SAMPLE_MATCHING != BarcodePosition.BOTH){
				if(BARCODE_FOR_SAMPLE_MATCHING == BarcodePosition.READ_1){
					//set bc2
					log.debug("Setting bc2 (not used in spl matching) to its seq : "+s2);
					bc2 = new BarcodeMatch();
					bc2.matched = true;
					bc2.barcode = s2;
					bc2.mismatches = 0;
					bc2.mismatchesToSecondBest = s2.length();
				}
				else{
					log.debug("Setting bc1 (not used in spl matching) to its seq : "+s1);
					bc1 = new BarcodeMatch();
					bc1.matched = true;
					bc1.barcode = s1;
					bc1.mismatches = 0;
					bc1.mismatchesToSecondBest = s1.length();
				}
			}


			
			
			
			/*
			 * Operate the barcode to sample matching
			 * At this point writers are set to null. Writers will be set only if barcodes
			 * can be mapped to a sample. 
			 */
			FastqWriter fwdW = null;
			FastqWriter revW = null;
			String spl1 = null; //the sample read_1 resolves to
			String spl2 = null; //the sample read_2 resolves to
			String chosen_sample = null; //final sample to attribute the read to
	
			
			/*
			 * Sample Lookup : set the fwdW, revW writers and return an array with samples variables (in this order) : spl1, spl2 and chosen_sample
			 * 
			 * but only if sample lookup is successful, else variables remain null 
			 * 
			 * The method delegates to methods specific for when using INDEX_FILE or not  
			 */
			 
			String [] arr = sampleLookupWithoutIndexFiles( bc1, bc2, bc2sample);
			spl1 = arr[0];
			spl2 = arr[1];
			chosen_sample = arr[2];
			
			if(chosen_sample!=null){
				fwdW = writers_fwd.get(chosen_sample);
				if(RUNNING_PAIRED_END)
					revW = writers_rev.get(chosen_sample);
			}
			
			log.debug("spl1="+spl1);
			log.debug("spl2="+spl2);
			log.debug("chosen_sample="+chosen_sample);
			if(bc1!=null)
				log.debug("bc1="+bc1.toString());
			else
				log.debug("bc1=NULL");
			
			if(bc2!=null)
				log.debug("bc2="+bc2.toString());
			else
				log.debug("bc2=NULL");
			
			//if chosen_sample == null, the read pair is unassigned
			boolean is_unassigned = false;
			if(chosen_sample == null){
				is_unassigned = true;
				unassigned++;
				if(KEEP_UNASSIGNED_READ){
					fwdW = writers_fwd.get(UNASSIGNED);
					if(RUNNING_PAIRED_END) 
						revW = writers_rev.get(UNASSIGNED);
				}
			}
			else{
				assigned++;
				//increase count for sample
				int c = sampleCountMap.get(chosen_sample);
				c++;
				sampleCountMap.put(chosen_sample, c);
			}

			//write diag no matter if read pair is assigned to a sample or not
			if(writingBarcodeDiagFile){
				String [] vals = null;
				if(RUNNING_PAIRED_END){
					//	headers = new String []{"Read_1_Header", "Read_2_Header", "seq1", "seq2" , "MM_Best_1", "MM_Second_1", "bc1", "smpl_1", "MM_Best_2", "MM_Second_2", "bc2", "smpl_2", "sample"};
					String mm1 = bc1 == null ? "" : bc1.mismatches+"";
					String mms1 = bc1 == null ? "" : bc1.mismatchesToSecondBest+"";
					String _bc1 =  bc1 == null ? "" : bc1.barcode;
					String mm2 = bc2 == null ? "" : bc2.mismatches+"";
					String mms2 = bc2 == null ? "" : bc2.mismatchesToSecondBest+"";
					String _bc2 =  bc2 == null ? "" : bc2.barcode;
					vals = new String []{
							r1.getReadHeader(), r2.getReadHeader(), s1, s2, 
							mm1 , mms1 , _bc1, spl1 , 
							mm2 , mms2 , _bc2, spl2 ,
							chosen_sample};
				}else{
					//headers =  new String []{"Read_Header", "seq", "MM_Best", "MM_Second", "bc", "sample"};
					String mm1 = bc1 == null ? "" : bc1.mismatches+"";
					String mms1 = bc1 == null ? "" : bc1.mismatchesToSecondBest+"";
					String _bc1 =  bc1 == null ? "" : bc1.barcode;

					vals = new String []{r1.getReadHeader(), s1, mm1 , mms1 , _bc1, chosen_sample};
				}

				diagPW.println(StringUtil.mergeArray(vals,"\t"));
			}

			if(!KEEP_UNASSIGNED_READ && chosen_sample == null){
				//then unassigned reads should not be saved, skip read writing
				continue;
			}

			/*
			 * write reads to out files
			 */
			if(!STATS_ONLY){
				/*
				 * special situation where barcodes have N's : give extracted sequence from read, not the matched barcodes
				 */
				if( !is_unassigned && barcodeValidator.barcodesHaveNs ){
					if(RUNNING_PAIRED_END){
						String bcForRead1 = (bc1 == null ? null : s1);
						String bcForRead2 = (bc2 == null ? null : s2);
						if(BARCODE_READ_POS == BarcodePosition.READ_2 ){
							bcForRead1 = (bc2 == null ? null : s2);
						}
						else if(BARCODE_READ_POS==BarcodePosition.READ_1){
							bcForRead2= (bc1 == null ? null : s1);
						}
						else if(BARCODE_READ_POS==BarcodePosition.BOTH && ENSURE_IDENTICAL_HEADER_NAMES){
							//we build the same bc for both
							if(bcForRead1!=null && bcForRead2!=null){
								bcForRead1 = bcForRead1+":"+bcForRead2;
							}else if(bcForRead2!=null ){
								bcForRead1 = bcForRead2;
							}
							bcForRead2= bcForRead1; // make bcForRead2 the same 
						}
						write(r1, BCLEN_1, XTRIMLEN_1, ZTRIMLEN_1, fwdW, bcForRead1, ADD_BARCODE_TO_HEADER, ENSURE_IDENTICAL_HEADER_NAMES);
						write(r2, BCLEN_2 , XTRIMLEN_2 , ZTRIMLEN_2, revW, bcForRead2, ADD_BARCODE_TO_HEADER, ENSURE_IDENTICAL_HEADER_NAMES );
					}else{
						//single end
						write(r1, BCLEN_1, XTRIMLEN_1, ZTRIMLEN_1, fwdW, (bc1 == null ? null : s1), ADD_BARCODE_TO_HEADER, false);
					}
				}else if(!is_unassigned){
					//normal situation
					if(RUNNING_PAIRED_END){
						String bcForRead1 = (bc1 == null ? null : bc1.barcode);
						String bcForRead2 = (bc2 == null ? null : bc2.barcode);
						if(BARCODE_READ_POS == BarcodePosition.READ_2 ){
							bcForRead1 = (bc2 == null ? null : bc2.barcode);
						}
						else if(BARCODE_READ_POS==BarcodePosition.READ_1){
							bcForRead2= (bc1 == null ? null : bc1.barcode);
						}
						else if(BARCODE_READ_POS==BarcodePosition.BOTH && ENSURE_IDENTICAL_HEADER_NAMES){
							//we build the same bc for both
							if(bcForRead1!=null && bcForRead2!=null){
								bcForRead1 = bcForRead1+":"+bcForRead2;
							}else if(bcForRead2!=null ){
								bcForRead1 = bcForRead2;
							}
							bcForRead2= bcForRead1;// make bcForRead2 the same 
						}
						write(r1, BCLEN_1, XTRIMLEN_1, ZTRIMLEN_1, fwdW, bcForRead1, ADD_BARCODE_TO_HEADER, ENSURE_IDENTICAL_HEADER_NAMES);
						write(r2, BCLEN_2 , XTRIMLEN_2 , ZTRIMLEN_2, revW, bcForRead2, ADD_BARCODE_TO_HEADER,ENSURE_IDENTICAL_HEADER_NAMES );
					}else{
						//single end
						write(r1, BCLEN_1, XTRIMLEN_1, ZTRIMLEN_1, fwdW, (bc1 == null ? null : bc1.barcode), ADD_BARCODE_TO_HEADER,false);
					}
				}else{
					/*
					 *  is_unassigned == TRUE situation
					 *  we want to print unmodified reads
					 */
					fwdW.write(r1);
					if(RUNNING_PAIRED_END){
						revW.write(r2);	
					}
				}

			}

		}

		
		//close readers
		try {
			fqr1.close();
			if(RUNNING_PAIRED_END && fqr2!=null)
				fqr2.close();
			
		} catch (Exception e) {
			// ignore
		}
		
		
		//close all writers
		if(!STATS_ONLY){
			for (FastqWriter w : writers_fwd.values()) {
				if(w!=null) w.close(); 
			}
			if(RUNNING_PAIRED_END){
				for (FastqWriter w : writers_rev.values()) {
					if(w!=null) w.close();
				}
			}
		}
		if(writingBarcodeDiagFile){
			if(diagPW!=null) diagPW.close();
		}

		//print counts
		printMetricFile(this, Je.COMMAND_MULTIPLEX, cnt, unassigned, assigned);

		//all fine

		return 0;

	}


	

	@Override
	public BarcodePosition getBarcodeReadPosOptionValue() {
		return BARCODE_READ_POS;
	}

	@Override
	public boolean getRedundantBarcodesOptionValue() {
		return REDUNDANT_BARCODES;
	}

	@Override
	public BarcodePosition getBarcodeForSampleMatchingOptionValue(){
		return BARCODE_FOR_SAMPLE_MATCHING;
	}
	
	@Override
	public  String getBCLenOptionValue(){
		return BCLEN;
	}
	
	@Override
	public  boolean getAddBarcodeToHeaderOptionValue(){
		return ADD_BARCODE_TO_HEADER;
	}




	/* (non-Javadoc)
	 * @see org.embl.gbcs.jemultiplexer.AbstractJemultiplexer#validateBarcodes()
	 */
	@Override
	protected List<String> validateBarcodes() {
		
		List<String> messages = new ArrayList<String>();
		
		/*
		 * Validate the barcode file format
		 */
		barcodeValidator = null; 
		if(BARCODE_FILE==null){
			if(USE_EMBASE)
				messages.add("Running EMBASE Mode: No Barcode file was created; this is most likely due to previous errors.");
			else
				messages.add("No Barcode file was given in command line!"); //this is normally NOT possible
		} 
		else if(!BARCODE_FILE.exists()){
			messages.add("BARCODE FILE does not exist :"+BARCODE_FILE.getAbsolutePath());
		}
		else if(!BARCODE_FILE.canRead()){
			messages.add("Unsufficient rights to read BARCODE FILE :"+BARCODE_FILE.getAbsolutePath());
		}
		else{
			barcodeValidator = new BarcodeValidator(BARCODE_FILE, this);
			//validate file format
			List<String> _messages = barcodeValidator.validate();

			if(_messages.size() != 0){
				messages.addAll(_messages);

				//end of the story
				return messages;
			}

			/*
			 * first inspect a delicate situation: when we have 2 barcodes but only one is meant for barcode matching
			 * and the BC file does not have a bc1:bc2 syntax (which is expected in case the second barcode is 
			 * a random one) 
			 * in this situation, we must have BCLEN defined from cmd line for BOTH barcodes 
			 * 
			 * this situation is only possible when no index file is given 
			 * 
			 */
			if(	 RUNNING_PAIRED_END 
					&& BARCODE_READ_POS==BarcodePosition.BOTH 
					&& !REDUNDANT_BARCODES && BARCODE_FOR_SAMPLE_MATCHING!= BarcodePosition.BOTH 
					&& !barcodeValidator.hasBarcodesForBothEnds ){

				if(BCLEN_1 == 0 && BCLEN_2 == 0){
					//we dont fail ONLY if CLIP barcode is false and ADD is false as well
					if(CLIP_BARCODE || ADD_BARCODE_TO_HEADER){
						messages.add("BCLEN for BOTH barcode must be specified in the command line (e.g. BCLEN=6:8) when" +
								" BARCODE_READ_POS==BOTH REDUNDANT_BARCODES=false BARCODE_FOR_SAMPLE_MATCHING="+BARCODE_FOR_SAMPLE_MATCHING.toString()+
								" and the barcode file does not distinguish between the barcodes of both ends (using bc1:bc2 syntax). " +
								"The situation is too ambigous.");
					}
				}

			}


			/*
			 * set barcode length if needed ; if still equals to 0 => BCLEN option was not
			 * given in cmd line ; get value from barcodeValidator if applicable 
			 */
			if(BCLEN_1 == 0){  
				if(!RUNNING_PAIRED_END || BARCODE_READ_POS!=BarcodePosition.READ_2)
					BCLEN_1 = barcodeValidator.bclen1;
			}
			else if (BCLEN_1 != 0 && RUNNING_PAIRED_END && BARCODE_READ_POS==BarcodePosition.READ_2){
				//it should be 0 if no barcode is described at this end
				BCLEN_1 = 0;
			}


			/*
			 * set barcode length if needed ; if still equals to 0 => BCLEN option was not given 
			 * in cmd line ; get value from barcodeValidator if applicable
			 */
			if(BCLEN_2 == 0 && RUNNING_PAIRED_END && (BARCODE_READ_POS==BarcodePosition.READ_2 || BARCODE_READ_POS == BarcodePosition.BOTH)){ 
				BCLEN_2 = barcodeValidator.bclen2;
			}
			else if (BCLEN_2 != 0 && RUNNING_PAIRED_END && BARCODE_READ_POS==BarcodePosition.READ_1){
				//it should be 0 if no barcode is described at this end
				BCLEN_2 = 0;
			}

			//more over, it does not make sense to have different BCLEN if barcodes are "redundant"
			if(RUNNING_PAIRED_END 
					&& BARCODE_READ_POS==BarcodePosition.BOTH && BARCODE_FOR_SAMPLE_MATCHING==BarcodePosition.BOTH 
					&& REDUNDANT_BARCODES && BCLEN_1 != BCLEN_2 ){
				messages.add("It does not make sense to have different barcode length (BCLEN1="+BCLEN_1+", BCLEN2="+BCLEN_2+") "
						+" when REDUNDANT_BARCODES=true ");
			}

		}

		return messages;
	}
	

	
	@Override
	protected String validateOptionAcceptsTwoValues(String optionName, String optionValue, boolean failIfPEwithSingleBarcode){
		String mess = null;
		if(!RUNNING_PAIRED_END || (failIfPEwithSingleBarcode && BARCODE_READ_POS!=BarcodePosition.BOTH)){
			mess = "Invalid value used for "+optionName+" : "+optionValue +". The X:Z synthax can only be used for paired-end data with BPOS=BOTH";
			
		}
		return mess;
	}
	
	@Override
	protected List<String> validateBarcodePositions(){
		if(RUNNING_PAIRED_END){
			//if we have a unique barcode then BARCODE_FOR_SAMPLE_MATCHING value should match value of BARCODE_READ_POS
			if(BARCODE_READ_POS!=BarcodePosition.BOTH)
				BARCODE_FOR_SAMPLE_MATCHING = BARCODE_READ_POS;
			
		}
		else{
			BARCODE_READ_POS=BarcodePosition.READ_1;
			BARCODE_FOR_SAMPLE_MATCHING=BarcodePosition.READ_1;
		}
		
		return new ArrayList<String>();
	}
	
	
	/**
	 * 
	 * Sample matching in the normal case of not using index files ie barcode embedded in read(s)
	 * 
	 * @param bc1
	 * @param bc2
	 * @param bc2sample
	 * @return a String array containing, in this order, the new values for 'spl1', 'spl2' and 'chosen_sample'
	 */
	private String [] sampleLookupWithoutIndexFiles(BarcodeMatch bc1, BarcodeMatch bc2, Map<String, String> bc2sample){

		
		String spl1 = null;
		String spl2 = null;
		String chosen_sample = null;
		
		//if running SE and we could match up a barcode
		if(!RUNNING_PAIRED_END && bc1!=null &&bc1.matched){
			chosen_sample = bc2sample.get(bc1.barcode);
			
			//stop method
			return new String[]{chosen_sample, null, chosen_sample};
		}


		/*
		 * if only one barcode should be used for sample look up (note: no matter if barcodes are present at both ends)
		 */				
		if(BARCODE_FOR_SAMPLE_MATCHING != BarcodePosition.BOTH ){
			BarcodeMatch bm = BARCODE_FOR_SAMPLE_MATCHING==BarcodePosition.READ_1 ? bc1:bc2;
			String spl = (bm==null || !bm.matched) ? null : bc2sample.get(bm.barcode);

			if(spl!=null){
				chosen_sample = spl;
			}
			
			return new String[]{
					BARCODE_FOR_SAMPLE_MATCHING==BarcodePosition.READ_1 ? chosen_sample:null, 
							BARCODE_FOR_SAMPLE_MATCHING==BarcodePosition.READ_2 ? chosen_sample:null, 
									chosen_sample};
		}
		/*
		 * if both barcodes are redundant
		 */
		else if(BARCODE_FOR_SAMPLE_MATCHING == BarcodePosition.BOTH && REDUNDANT_BARCODES == true){
			spl1 = (bc1==null || !bc1.matched) ? null : bc2sample.get(bc1.barcode);
			spl2 = (bc2==null || !bc2.matched) ? null : bc2sample.get(bc2.barcode);

			if(spl1 != null &&  spl2 != null){
				if(spl1.equals(spl2)){
					chosen_sample=spl1;
				}else if(!STRICT){
					//which is the best match ?
					if(bc1.mismatches > bc2.mismatches){
						//bc2 is best
						chosen_sample=spl2;
					}
					else if(bc1.mismatches < bc2.mismatches){
						//bc1 is best
						chosen_sample=spl1;
					}
					else{
						//look int delta : we don t do this anymore in v1.2
						// this was a bug as it always resolve to the same sample when barcodes have no mismatches (see above)  
//						if(bc1.mismatchesToSecondBest!=bc2.mismatchesToSecondBest){
//							chosen_sample = bc1.mismatchesToSecondBest>bc2.mismatchesToSecondBest? spl1:spl2;
//						}
					}
				}
			}
			else if(!STRICT && spl1 != null){
				chosen_sample=spl1;
			}
			else if(!STRICT && spl2 != null){
				chosen_sample=spl2;					
			}

			
			return new String[]{spl1, spl2, chosen_sample};

		}
		/*
		 * if both barcode are needed to lookup sample
		 */
		else if(BARCODE_FOR_SAMPLE_MATCHING == BarcodePosition.BOTH && REDUNDANT_BARCODES == false){
			//we need a match for both 
			if(bc1!=null && bc1.matched && bc2!=null && bc2.matched){
				spl1="NA";
				spl2="NA";
				//then the lookup is made with concatenated barcodes 
				chosen_sample = bc2sample.get(bc1.barcode+bc2.barcode);
				
			}
			
			return new String[]{spl1, spl2, chosen_sample};
		}
		else{
			//unexpected situation reflecting a bug
			throw new RuntimeException("Don t know how to rematch barcode(s) to sample : unexpected code dead end reached ! Please report.");
		}


	}
	

	
	@Override
	protected List<String> prepareForEmBASE() {
		//override whatever user could have set for some options
		//
		
		//error messages to be accumulated
		List<String> messages = new ArrayList<String>();
		
		
		OUTPUT_DIR = null;
		//we keep unassigned reads
		KEEP_UNASSIGNED_READ = true;
		//and call files accordign to input files, simply adding '_unassigned-reads' before extension
		UNASSIGNED_FILE_NAME_1= FileUtil.removeExtension( FASTQ_FILE1.getName() ) + "_unassigned-reads"+".txt.gz" ;
		UNASSIGNED_FILE_NAME_2= RUNNING_PAIRED_END ? FileUtil.removeExtension( FASTQ_FILE2.getName() ) + "_unassigned-reads"+".txt.gz" : null;
		GZIP_OUTPUTS=true;
		STATS_ONLY=false;
		CREATE_MD5_FILE = true;
		
		/*
		 * Get information from emBASE 
		 */		
		
		//Determine the caller
		String username = System.getProperty("user.name");
		//check user is not trying to trick me with a -Duser.name=another user!
		String whoami = JeUtils.whoami();
		if(!username.equalsIgnoreCase(whoami)){
			messages.add("User "+username+" does not match 'whoami' ("+whoami+") ! Are you trying to hack the system ?");
			//we cannot continue
			return messages;
		}
		
		EmBaseDatabaseFacade facade = null;
		try {
			facade = EmBaseDatabaseFacade.createEmBaseDatabaseFacade(
					username,
					ApplicationConfiguration.DB_URL,
					ApplicationConfiguration.DB_USER,
					ApplicationConfiguration.DB_PWD
					);
		} catch (EmBASEConnectionException e) {
			log.error(ExceptionUtil.getStackTrace(e));
			messages.add("Error when trying to connect emBASE DB : "+e.getMessage());
			//we cannot continue
			return messages;
		}
		
		//what s the lane number ? first get teh NGS assay
		Item ngsAssay = facade.getNGSAssay(FASTQ_FILE1);
		List<CSVLine> rbaDetails = null;
		if(ngsAssay == null){
			messages.add("No NGS Assay found in emBASE for lane file : "+FASTQ_FILE1.getAbsolutePath()+"\nAre you sure this file has been already imported in emBASE ? If so, please contact emBASE admins.");
		}else{
			//then fetch details
			ArrayList<Item> _list = new ArrayList<Item>();
			_list.add(ngsAssay);
			/*
			 * rbaDetails contains values for all org.embl.gbcs.embase.api.Queries.HEADER_*
			 */
			rbaDetails = facade.fetchRBADetailsForNGSAssays(_list).values().iterator().next();
			if(rbaDetails == null || rbaDetails.size()==0){
				messages.add("No RawBioAssay found in emBASE for NGS Assay "+ngsAssay.getName()+" (id="+ngsAssay.getId()+") i.e. lane file : "+FASTQ_FILE1.getAbsolutePath()+"\nAre you sure this file has been successfully imported in emBASE ? If so, please contact emBASE admins.");
			}
		}
		
		if(messages.size()!=0){
			//we stop here
			return messages;
		}
		
		//get lane number : any RBA can be used for this
		Integer laneNumber = null;
		try {
			laneNumber = rbaDetails.get(0).getValueAsInt(Queries.HEADER_Lane);
		} catch (InvalidHeaderException e1) {
			//not possible if the code has been run at least once!
			throw new RuntimeException(e1);
		}
		
		// Check that expected directories exist
		File runDir = FASTQ_FILE1.getParentFile();
		//get lane dir
		File laneDir = new File (runDir, EmBaseDatabaseFacade.getLaneDirName(laneNumber));
		//then check if expected dir exists for stats and unassigned reads
		File runStatDir = new File (laneDir, EmBaseDatabaseFacade.STAT_DIR_NAME);
		
		if(!runStatDir.exists() || !runStatDir.canWrite()){
			if(!TEST_MODE_STOP_AFTER_PARSING)
				messages.add("The expected '"+EmBaseDatabaseFacade.STAT_DIR_NAME+"' dir in the lane dir '"+laneDir.getAbsolutePath()+"' is either not existing or not writable");
			else{
				METRICS_FILE_NAME = new File(EmBaseDatabaseFacade.TMP_DIR, FileUtil.removeExtension( FASTQ_FILE1.getName() ) +"_jemultiplexing_metrics.txt").getAbsolutePath();
				log.info("TEST MODE: METRICS_FILE_NAME set to "+METRICS_FILE_NAME);
			}
		}else{
			//set the full path so the remaining of the code leaves it like this
			METRICS_FILE_NAME = new File(runStatDir, FileUtil.removeExtension( FASTQ_FILE1.getName() ) +"_jemultiplexing_metrics.txt").getAbsolutePath();
			log.info("METRICS_FILE_NAME set to "+METRICS_FILE_NAME);
		}
		
		File unassignedReadsDir = new File (laneDir, EmBaseDatabaseFacade.UNASSIGNED_DIR_NAME);
		if(!unassignedReadsDir.exists() || !unassignedReadsDir.canWrite()){
			if(!TEST_MODE_STOP_AFTER_PARSING)
				messages.add("The expected '"+EmBaseDatabaseFacade.UNASSIGNED_DIR_NAME+"' dir in the run dir '"+laneDir.getAbsolutePath()+"' is either not existing or not writable");
			else{
				OUTPUT_DIR = EmBaseDatabaseFacade.TMP_DIR;
				log.info("TEST MODE: OUTPUT_DIR set to "+OUTPUT_DIR.getAbsolutePath());
			}
		}else{
			 //set output dir to where the input files are located => this will be the place for unassigned reads 
			OUTPUT_DIR = unassignedReadsDir;
			log.info("OUTPUT_DIR set to "+OUTPUT_DIR.getAbsolutePath());
		}

		
		/*
		 * Check caller can actually demultiplex this dataset
		 * The user should have read access on the NGS assay
		 */
		try {
			if(!facade.canDemultiplex(FASTQ_FILE1, username)){
				messages.add("According to emBASE DB, the user "+username+" is not allowed to read NGS Assay corresponding to file "+FASTQ_FILE1.getAbsolutePath());
			}else{
				log.debug("User allowed to demultiplex");
			}
		}catch (Exception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			messages.add("Error when trying to check read rights from emBASE DB : "+e.getMessage());
		}
		
		
		/*
		 * Get all libs 
		 */
		List<NGSLibrary> libraries = null;
		try {
			libraries = facade.fetchLibraries(FASTQ_FILE1);
			facade.findSampleStorageDir(laneDir, libraries);
			
			//and check dirs have been set for all libs
			for (NGSLibrary lib : libraries) {
				if(lib.getSampleDir() == null || !lib.getSampleDir().exists()){
					String mess = "No library dir for sample/library "+lib.getName()+" is available.\nThe emBASE storage architecture is not ready. " +
							"Are you sure this lane has been successfully uploaded to emBASE ?\nIf so, please contact emBASE admins.";
					log.debug(mess);
					if(TEST_MODE_STOP_AFTER_PARSING){
						log.debug("TEST MODE : Setting fake dir for "+lib.getName());
						lib.setSampleDir ( EmBaseDatabaseFacade.getFakeSampleDir(lib.getName(), lib.getId()) );//does not exist, just fake for test
					}
					else{
						messages.add(mess);
					}
				}
				else {
					//update lib.getSampleDir with fastq sub-dir
					File fastqDir = new File(lib.getSampleDir(), "fastq");
					if(!fastqDir.exists() || !fastqDir.canWrite()){
						String mess = "The 'fastq' dir under the library dir  "+lib.getSampleDir().getAbsolutePath()
								+" is either not existing OR not writable.\n";
						if(!fastqDir.exists()){
							mess+="The storage architecture is not ready. Are you sure this lane has been successfully uploaded to emBASE ?\nIf so, please contact emBASE admins.";
						}
						else{
							mess+="The fastq storage file has been already locked to prevent further modifications.\n" +
									"Please contact emBASE admins if you think this should not be the case.";
						}
						log.debug(mess);
						messages.add(mess);
					}else{
						lib.setSampleDir( fastqDir );
					}
				}
			}
		} catch (Exception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			messages.add("Error when trying to fetch barcodes from emBASE DB : "+e.getMessage());
		}
		
		if(messages.size()!=0){
			//we stop here
			return messages;
		}
		
		/*
		 * All seems perfect, create barcode file 
		 * Check demultiplexed files are not already found
		 */
		PrintWriter pw = null;
		try{
			BARCODE_FILE = getBarcodeFileTmpLocation(FASTQ_FILE1, username, runStatDir);
			if(BARCODE_FILE.exists()){
				log.warn("The barcode file already exists and will be overwritten : "+BARCODE_FILE.getAbsolutePath());
			}
			pw = new PrintWriter(BARCODE_FILE);
			for (NGSLibrary lib : libraries) {
				//File fastq1 = getDemultiplexedFile(lib, FileUtil.removeExtension(FASTQ_FILE1.getName()) , GZIP_OUTPUTS);
				//File fastq2 = RUNNING_PAIRED_END ? getDemultiplexedFile(lib, FileUtil.removeExtension(FASTQ_FILE2.getName()), GZIP_OUTPUTS) : null;
				
				File fastq1 = getDemultiplexedFile(lib, 1 , GZIP_OUTPUTS);
				File fastq2 = getDemultiplexedFile(lib, 2 , GZIP_OUTPUTS);
				
				
				/*
				 * check if file already exists , the test on the length is meant to avoid considering a 'fake' file as a real file
				 * we indeed create such files to prepare the architecture
				 */
				if(fastq1.exists() && fastq1.length() > 10*1000L && !FORCE){ 
					if(fastq1.canWrite()){
						messages.add("Demultiplexed file "+fastq1.getAbsolutePath()+" already exists. You need to run with FORCE=true to overwrite existing files.");
					}else{
						messages.add("Demultiplexed file "+fastq1.getAbsolutePath()+" already exists BUT you do NOT have write permission on it.\n" +
								"You first need to solve the write permission issue, then use FORCE=true to overwrite existing files.");
					}
				}
				else if(fastq1.exists() && FORCE && fastq1.canWrite()){
					log.warn("Existing demultiplexed file will be overritten (option FORCE=true) : "+fastq1.getAbsolutePath());
				}
				
				if(RUNNING_PAIRED_END && fastq2.exists() & fastq2.length() > 10*1000L && !FORCE){
					if(fastq2.canWrite()){
						messages.add("Demultiplexed file "+fastq2.getAbsolutePath()+" already exists. You need to run with FORCE=true to overwrite existing files.");
					}else{
						messages.add("Demultiplexed file "+fastq2.getAbsolutePath()+" already exists BUT you do NOT have write permission on it.\n" +
								"You first need to solve the write permission issue, then use FORCE=true to overwrite existing files.");
					}
				}else if(RUNNING_PAIRED_END && fastq2.exists() && FORCE && fastq2.canWrite()){
					log.warn("Existing demultiplexed file will be overritten (option FORCE=true) : "+fastq2.getAbsolutePath());
				}
				
				pw.print(lib.getName()+"\t"+lib.getBarcode()+"\t"+fastq1.getAbsolutePath());
				if(RUNNING_PAIRED_END)
					pw.print("\t"+fastq2.getAbsolutePath());
				pw.println();
			}
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}finally{
			if(pw!=null)
				pw.close();
		}
		
		//then ok just let the code go on
		
		return messages;
	}


	/**
	 * @param args
	 */
	public static void main(final String[] argv) {
		new Jemultiplexer().instanceMainWithExit(argv); //eventually execute the doWork() method
	}
	
}
