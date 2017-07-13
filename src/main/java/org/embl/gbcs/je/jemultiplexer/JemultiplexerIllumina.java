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
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.embl.cg.utilitytools.utils.StringUtil;
import org.embl.gbcs.je.BarcodeMatch;
import org.embl.gbcs.je.Je;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Illumina;


@CommandLineProgramProperties(
		usage = "\tFastq files demultiplexer for Illumina Fastq files using Illumina Index files.\n"
				+ "\tFastq files (reads and index) can be in gzip compressed (end in .gz). \n" 
				+ "\tBy default output files are gzipped and have names following the pattern \n"
				+ "\t'<samplename>_<barcode[-barcode2...-barcodeN]>_[1|2].txt[.gz]' unless you gave file\n"
				+ "\t names to use within the barcode description file.\n"
				+ "Example : \n"
				+ "\t je demultiplex-illu F1=fastq_1.txt.gz I1=index_1.txt.gz BF=barcodes.bs O=~/Desktop/test-jemultiplexer2/ ",
				usageShort = "Demultiplex fastq files with indices provided in separate index files",
				programGroup = Illumina.class 
		)
public class JemultiplexerIllumina extends AbstractJemultiplexer {
	private static Logger log = LoggerFactory.getLogger(JemultiplexerIllumina.class);

	@Option(shortName="I1", optional = false,
			printOrder=22,
			doc="Fastq file for index 1 (barcode) reads, optionally gzipped.\n" 
			)
	public File INDEX_FILE1;

	@Option(shortName="I2", optional = true,
			printOrder=24,
			doc="Fastq file for index 2 (barcode) reads, optionally gzipped.\n" +
					"A INDEX_FILE1 MUST be provided when INDEX_FILE2 is given. This situation corresponds to Illumina dual indexing.\n"
					
			)
	public File INDEX_FILE2 = null;

	
	@Option(shortName="BPOS", optional = true,
			printOrder=41,
			doc="Indicates the location of additional barcodes present in the read(s). Setting this option implies setting the LEN option.\n"+
					"\tImportantly, these additional barcodes must not encode sample identity information but used for \n"+
					"\te.g. molecular barcoding (UMIs) or for any purpose other than sample identity encoding."
			)
	public BarcodePosition BARCODE_READ_POS = DEFAULT_BARCODE_READ_POS ;

	@Option(shortName = "LEN", optional = true,
			printOrder=43,
			doc = "Length of the additional barcodes present in the read(s) as indicated by the BPOS option. " +
					"Two distinct length can be provided using the syntax LEN=X:Z where X and Z are 2 integers representing" +
					" the barcode length for read_1 and read_2 respectively.\n"+
					"Only relevant when BPOS != NONE."
			)
	public String BCLEN=null ; //BCLEN parsed to set bclen1 and bclen2
	
	
	
	@Option(shortName="BRED", optional = true,
			printOrder=45,
			doc="This option only applies for paired-end data with *both* INDEX_FILE1 and INDEX_FILE2 provided.\n" +
					"Indicates if both index barcodes encode redundant information i.e. if both barcodes are supposed" +
					" to be identical (or resolve to the same sample when a pool of barcodes is used per sample).\n "+
					"\tWhen BRED=true, the STRICT option guides the sample lookup behavior" +
					"\tWhen BRED=false, barcodes are combined prior to sample lookup.\n"
			)
	public boolean REDUNDANT_BARCODES = DEFAULT_REDUNDANT_BARCODES;

	
	
	@Option(shortName="ADD", optional = true,
			printOrder=112,
			doc="Add matched barcode at the end of the read header. Applies to both index when INDEX_FILE2 is also provided.\n" +
					"\t"+"First the sample encoding barcodes from I1 (and I2 when relevant) are added to the read headers like \n"+
					"\t"+"\t"+"'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:I1_BARCODE:I2_BARCODE'\n"+
					"\t"+"Then, if BPOS!=NONE, the additional barcodes (UMIs) clipped from the read(s) are added to their own header, like \n"+
					"\t"+"\t"+"'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:I1_BARCODE:I2_BARCODE:CLIPPED_SEQ_FROMREAD'\n"
			)
	public boolean ADD_BARCODE_TO_HEADER  = DEFAULT_ADD_BARCODE_TO_HEADER;


	@Option(shortName="SAME_HEADERS", optional = true,
			printOrder=115,
			doc= "Makes sure that headers of both reads of a pair are identical, using the following read header pattern (for both reads of a pair) :\n"
					+ "\t\t'@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 I1_BARCODE:I2_BARCODE(if applicable)':CLIPPED_SEQ_FROMREAD1:CLIPPED_SEQ_FROMREAD2 \n"
					+ "This option only makes sense in paired end mode and ADD=true. "
					+ "Some (if not all) mappers will indeed complain when"
					+ " the read headers are not identical. When molecular barcodes are"
					+ " present in reads and the RCHAR is used, you will end"
					+ " with (problematic) read headers like this :\n"
					+ "\t\tHISEQ:44:C6KC0ANXX:5:1101:1491:1994:1:N:0:TAGAACAC:TGGAGTAG\n"
					+ "\t\tHISEQ:44:C6KC0ANXX:5:1101:1491:1994:3:N:0:TAGAACAC:CGTTGTAT\n"
					+ "SAME_HEADERS=true will instead genetates the following identical header for both reads :\n"
					+ "\t\tHISEQ:44:C6KC0ANXX:5:1101:1491:1994:TAGAACAC:TGGAGTAG:CGTTGTAT\n"
					+ "Note that we also clipped the useless '1:N:0' and '3:N:0' has they will also result in generating different headers\n"
					+ "\t Important : this option will force RCHAR=: UNLESS you specify RCHAR=null ; in which case a space will be preserved ie : \n"
					+ "\t\tHISEQ:44:C6KC0ANXX:5:1101:1491:1994 TAGAACAC:TGGAGTAG:CGTTGTAT\n"
			)
	public boolean ENSURE_IDENTICAL_HEADER_NAMES  = DEFAULT_ENSURE_IDENTICAL_HEADER_NAMES;



	public JemultiplexerIllumina() {
		super();
		USE_SAMPLE_INDEX_FILES = true;
	}
	
	
	@Override	
	protected String[] customCommandLineValidation() {
		/*
		 * run custom validation that does NOT depend on other validation 
		 *  
		 */
		List<String> messages = new ArrayList<String>();
		

		/*
		 * validate FASTQ files
		 */
		messages.addAll( validateFASTQFiles());
		
		if(FASTQ_FILE2 != null)
			RUNNING_PAIRED_END = true;

		
		/*
		 * Validate INDEX files when applicable 
		 */
		messages.addAll( validateIndexFiles());

		
		if(ENSURE_IDENTICAL_HEADER_NAMES && RUNNING_PAIRED_END ){
			// check RCHAR : if this is not explicitely given in cmd line, then turn it to :
			String _cmdline = Arrays.toString(this.getCommandLineParser().getArgv());
			if(!_cmdline.contains("READ_NAME_REPLACE_CHAR") && ! _cmdline.contains("RCHAR")){
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

		
		/*
		 * we still need to init write for unassigned barcodes
		 * will be init below only if applicable
		 */
		FastqWriter writer_index1 = null;
		FastqWriter writer_index2 = null;
		if(KEEP_UNASSIGNED_READ && !STATS_ONLY){
			if(INDEX_FILE1!=null){
				String unIndex1Fname = "unassigned_"+INDEX_FILE1.getName() ;
				//make sure a gz is there if needed
				if(GZIP_OUTPUTS && !unIndex1Fname.toUpperCase().endsWith(".GZ")){
					unIndex1Fname += ".gz";
				}

				//should not be the case at this point but...
				writer_index1 = fastqFactory.newWriter(new File(unassignedFile1.getParentFile(), unIndex1Fname), GZIP_OUTPUTS, CREATE_MD5_FILE);
			}
			if(INDEX_FILE2!=null){
				String unIndex2Fname = "unassigned_"+INDEX_FILE2.getName() ;
				//make sure a gz is there if needed
				if(GZIP_OUTPUTS && !unIndex2Fname.toUpperCase().endsWith(".GZ")){
					unIndex2Fname += ".gz";
				}

				//should not be the case at this point but...
				writer_index2 = fastqFactory.newWriter(new File(unassignedFile2.getParentFile(), unIndex2Fname), GZIP_OUTPUTS, CREATE_MD5_FILE);
			}

		}
		
		
		/*
		 * init readers
		 */
		FastqReader fqr1 = new FastqReader(FASTQ_FILE1);
		Iterator<FastqRecord> it1 = fqr1.iterator();
		FastqReader fqr2 = null;
		Iterator<FastqRecord> it2 = null;

		FastqReader idxr1 = null;
		Iterator<FastqRecord> idxit1 = null;

		FastqReader idxr2 = null;
		Iterator<FastqRecord> idxit2 = null;

		
		if(RUNNING_PAIRED_END){
			fqr2 = new FastqReader(FASTQ_FILE2);
			it2 = fqr2.iterator();
		}
		
		idxr1 = new FastqReader(INDEX_FILE1);
		idxit1 = idxr1.iterator();

		if(INDEX_FILE2!=null){
			idxr2 = new FastqReader(INDEX_FILE2);
			idxit2 = idxr2.iterator();
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
			
			FastqRecord idx_r1 = idxit1.next();
			FastqRecord idx_r2 = null;
			
			s1 = idx_r1.getReadString();
			q1 = idx_r1.getBaseQualityString();
			log.debug("BC1="+s1);
			if(INDEX_FILE2!=null){
				idx_r2 = idxit2.next();
				s2 = idx_r2.getReadString();
				q2 = idx_r2.getBaseQualityString();
				log.debug("BC2="+s2);
			}
			
			
			/*
			 * Operate barcode MATCHING
			 */
			bc1 = findBarcode(BarcodePosition.READ_1, s1, q1, bc2sample, sample2bcSet);
			log.debug("bc1 match : "+ (bc1 == null?"NULL":bc1.toString()));
			
			//do we need to match barcode 2? 
			if(	INDEX_FILE2!=null){ 
				bc2 = findBarcode(BarcodePosition.READ_2, s2, q2, bc2sample, sample2bcSet);
				log.debug("bc2 match : "+ (bc2 == null?"NULL":bc2.toString()));
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
			 
			String [] arr = sampleLookupWithIndexFiles( bc1, bc2, bc2sample);
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
				 * up to now, index files can only barcode encoding sample information ie no Ns can be found in there
				 * But there might be additional barcodes in the reads to extract
				 */

				if(!is_unassigned && RUNNING_PAIRED_END){
					String bcForRead1 = (bc1 == null ? null : bc1.barcode);
					String bcForRead2 = (bc2 == null ? null : bc2.barcode);
					if(INDEX_FILE2 == null ){
						//then it is the same barcode for read2 as read1
						bcForRead1 = (bcForRead1 == null ? "" : bcForRead1);
						bcForRead2= bcForRead1;
					}else {
						//when 2 INDEX_FILES and BRED=false => both barcodes should be added to each read like I1:I2
						String concatBCs =  (bcForRead1 == null ? null : bcForRead1)+":"+(bcForRead2 == null ? null : bcForRead2);
						bcForRead1 = concatBCs;
						bcForRead2 = concatBCs;
					}

					if(ENSURE_IDENTICAL_HEADER_NAMES){
						writePEWithIdenticalBarcodesForIndexFiles(r1, r2, bcForRead1, bcForRead2, BCLEN_1, BCLEN_2, XTRIMLEN_1, XTRIMLEN_2, ZTRIMLEN_1, ZTRIMLEN_2, fwdW, revW, ADD_BARCODE_TO_HEADER);
					}
					else{
						//previous code
						writeWithIndexFiles(r1, BCLEN_1, XTRIMLEN_1, ZTRIMLEN_1, fwdW, bcForRead1, ADD_BARCODE_TO_HEADER);
						writeWithIndexFiles(r2, BCLEN_2 , XTRIMLEN_2 , ZTRIMLEN_2, revW, bcForRead2, ADD_BARCODE_TO_HEADER);

					}

				}else if(!is_unassigned){
					//single end
					// write(r1, BCLEN_1, XTRIMLEN_1, ZTRIMLEN_1, fwdW, (bc1 == null ? null : bc1.barcode), ADD_BARCODE_TO_HEADER, false); //bug reported by Mark Heron
					writeWithIndexFiles(r1, BCLEN_1, XTRIMLEN_1, ZTRIMLEN_1, fwdW, (bc1 == null ? null : bc1.barcode), ADD_BARCODE_TO_HEADER);
				}else{
					/*
					 *  is_unassigned == TRUE situation
					 *  we want to print unmodified reads
					 */
					fwdW.write(r1);
					if(RUNNING_PAIRED_END){
						revW.write(r2);	
					}
					//index 
					writer_index1.write(idx_r1);
					if(INDEX_FILE2!=null){
						writer_index2.write(idx_r2);
					}
				}
			}

		}

		
		//close readers
		try {
			fqr1.close();
			if(RUNNING_PAIRED_END && fqr2!=null)
				fqr2.close();
			if(idxr1!=null)
				idxr1.close();
			if(RUNNING_PAIRED_END && idxr2!=null)
				idxr2.close();
			
		} catch (Exception e) {
			// TODO: handle exception
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
		printMetricFile(this, Je.COMMAND_MULTIPLEX_ILLUMINA, cnt, unassigned, assigned);

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
		return BarcodePosition.NONE;
	}
	
	@Override
	public  String getBCLenOptionValue(){
		return BCLEN;
	}
	
	@Override
	public  boolean getAddBarcodeToHeaderOptionValue(){
		return ADD_BARCODE_TO_HEADER;
	}
	
	@Override
	protected List<String> validateBarcodes() {

		List<String> messages = new ArrayList<String>();

		/*
		 * Validate the barcode file format
		 */
		barcodeValidator = null; 
		if(BARCODE_FILE==null){
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
			/*
			 * when using index files, the lenght are barcodes found in the BC file has not been checked against what is found in the index files !
			 */
			_messages.addAll(
					barcodeValidator.checkBarcodesLengthWithIndexFiles(INDEX_FILE1, INDEX_FILE2)
					);


			if(_messages.size() != 0){
				messages.addAll(_messages);
				//end of the story
				return messages;
			}

			/*
			 * set barcode length if needed ; if still equals to 0 => BCLEN option was not
			 * given in cmd line ; get value from barcodeValidator if applicable 
			 */
			if(BCLEN_1 == 0){  
				if (BARCODE_READ_POS == BarcodePosition.READ_1 || BARCODE_READ_POS == BarcodePosition.BOTH){
					messages.add("barcode lenght for READ 1 cannot be 0 with BPOS="+BARCODE_READ_POS);
				}
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
				if(BARCODE_READ_POS == BarcodePosition.READ_2 || BARCODE_READ_POS == BarcodePosition.BOTH){
					messages.add("barcode lenght for READ 2 cannot be 0 with BPOS="+BARCODE_READ_POS);
				}
			}
			else if (BCLEN_2 != 0 && RUNNING_PAIRED_END && BARCODE_READ_POS==BarcodePosition.READ_1){
				//it should be 0 if no barcode is described at this end
				BCLEN_2 = 0;
			}

		}


		return messages;
	}

	
	
	protected List<String> validateIndexFiles() {
		
		List<String> messages = new ArrayList<String>();
			
		// INDEX_FILE1 is mandatory
		if(!INDEX_FILE1.exists())
			messages.add("INDEX_FILE1 does not exist :"+INDEX_FILE1.getAbsolutePath());

		if(INDEX_FILE1.exists() &&!INDEX_FILE1.canRead())
			messages.add("Unsufficient rights to read INDEX_FILE1 :"+INDEX_FILE1.getAbsolutePath());

		if(INDEX_FILE1.getAbsolutePath().equals(FASTQ_FILE1.getAbsolutePath()))
			messages.add("FASTQ_FILE1 and INDEX_FILE1 are the same file !");

		if(RUNNING_PAIRED_END && INDEX_FILE1.getAbsolutePath().equals(FASTQ_FILE2.getAbsolutePath()))
			messages.add("FASTQ_FILE2 and INDEX_FILE1 are the same file !");

		
		
		if(INDEX_FILE2 != null){
			
			if(INDEX_FILE1 == null)
				messages.add("INDEX_FILE2 provided while INDEX_FILE1 was not given : can't give INDEX_FILE2 wihtout INDEX_FILE1");
			
			if(!INDEX_FILE2.exists())
				messages.add("INDEX_FILE2 does not exist :"+INDEX_FILE2.getAbsolutePath());

			if(INDEX_FILE2.exists() &&!INDEX_FILE2.canRead())
				messages.add("Unsufficient rights to read INDEX_FILE2 :"+INDEX_FILE2.getAbsolutePath());
			
			if(INDEX_FILE2.getAbsolutePath().equals(INDEX_FILE1.getAbsolutePath()))
				messages.add("INDEX_FILE1 and INDEX_FILE2 are the same file !");
			
			if(INDEX_FILE2.getAbsolutePath().equals(FASTQ_FILE1.getAbsolutePath()))
				messages.add("FASTQ_FILE1 and INDEX_FILE2 are the same file !");
			
			if(INDEX_FILE2.getAbsolutePath().equals(FASTQ_FILE2.getAbsolutePath()))
				messages.add("FASTQ_FILE2 and INDEX_FILE2 are the same file !");
			
			
		}
		
		return messages;
	}
	
	
	@Override
	protected String validateOptionAcceptsTwoValues(String optionName, String optionValue, boolean failIfPEwithSingleBarcode){
		String mess = null;
		if(failIfPEwithSingleBarcode && !optionName.equalsIgnoreCase("BCLEN")){
			mess = "Invalid value used for "+optionName+" : "+optionValue +". The X:Z synthax can only be used with 2 sample index files";
		}
		else if(!RUNNING_PAIRED_END || (failIfPEwithSingleBarcode && BARCODE_READ_POS!=BarcodePosition.BOTH)){
			mess = "Invalid value used for "+optionName+" : "+optionValue +". The X:Z synthax can only be used for paired-end data with BPOS=BOTH";
		}
		return mess;
	}
	
	
	@Override
	protected List<String> validateBarcodePositions(){
		List<String> messages = new ArrayList<String>();

		/*
		 * if BPOS not given then turn it to NONE
		 */
		String _cmdline = Arrays.toString(this.getCommandLineParser().getArgv());
		log.debug(_cmdline);
		if(!_cmdline.contains("BPOS") && ! _cmdline.contains("BARCODE_READ_POS")){
			log.debug("BPOS absent from command line : turning it to NONE");
			BARCODE_READ_POS = BarcodePosition.NONE;
		}
		
		/*
		 * if 2 index files were given and BRED not given, turn it to false
		 */
		if(INDEX_FILE2 !=null && !_cmdline.contains("BRED") && ! _cmdline.contains("REDUNDANT_BARCODES")){
			log.debug("[with sample index files] BRED not given in cmd line : turning it to false");
			REDUNDANT_BARCODES = false;
		}
		
		//BCLEN must be set to null unless BPOS != NONE as BCLEN is used for internal barcodes when using index files
		if(BARCODE_READ_POS == BarcodePosition.NONE){
			BCLEN = null;
		}
		
		
		if(!RUNNING_PAIRED_END){
			if(BARCODE_READ_POS==BarcodePosition.BOTH || BARCODE_READ_POS==BarcodePosition.READ_2){
				messages.add("Can't have barcode(s) in position '"+BARCODE_READ_POS+"' in single end sequencing mode" );
			}
		}
		
		/*
		 * if barcodes are also in the reads, BLEN must be given
		 */
		if(BARCODE_READ_POS != BarcodePosition.NONE && StringUtils.isBlank(BCLEN)){
			messages.add("BCLEN must be given to describe the length of barcode(s) in read(s) : '"+BARCODE_READ_POS );
		}
			
		//change defaults of XT 
		if(!_cmdline.contains("XT") && !_cmdline.contains("XTRIMLEN")){
			//user did not specify this option, it must then be defaulted to 0 
			this.XTRIMLEN = "0";
		}
		
		
		return messages;
	}

	
	/**
	 * 
	 * sample matching in case of use of index files (I1 and I2)
	 * 
	 * @param bc1
	 * @param bc2
	 * @param bc2sample
	 * 
	 * @return a String array containing, in this order, the new values for 'spl1', 'spl2' and 'chosen_sample'
	 */
	private String [] sampleLookupWithIndexFiles(BarcodeMatch bc1, BarcodeMatch bc2, 
			Map<String, String> bc2sample){
		
		String spl1 = null;
		String spl2 = null;
		String chosen_sample = null;
		
		if(bc1!=null)
			log.debug("[sampleLookupWithIndexFiles] bc1="+bc1.toString());
		else
			log.debug("[sampleLookupWithIndexFiles] bc1=NULL");
		
		if(bc2!=null)
			log.debug("[sampleLookupWithIndexFiles] bc2="+bc2.toString());
		else
			log.debug("[sampleLookupWithIndexFiles] bc2=NULL");
		
		//only one barcode then it MUST be bc1 as INDEX_FILE2 == null
		if(INDEX_FILE2 == null){
			chosen_sample = (bc1==null || !bc1.matched) ? null : bc2sample.get(bc1.barcode); 
			return new String[]{chosen_sample, null, chosen_sample};
		}
		
		/*
		 * We then have 2 barcodes.
		 * Case 1 : both barcodes are redundant
		 */
		if(REDUNDANT_BARCODES == true){
			log.debug("[sampleLookupWithIndexFiles] 2 index file with redundant barcodes");
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
						//look int delta
						if(bc1.mismatchesToSecondBest!=bc2.mismatchesToSecondBest){
							chosen_sample = bc1.mismatchesToSecondBest>bc2.mismatchesToSecondBest? spl1:spl2;
						}
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
		 * Case 2: both barcodes are needed to lookup sample
		 */
		else {
			log.debug("[sampleLookupWithIndexFiles] 2 index file with non-redundant barcodes");
			//we need a match for both 
			if(bc1!=null && bc1.matched && bc2!=null && bc2.matched){
				spl1="NA";
				spl2="NA";
				//then the lookup is made with concatenated barcodes 
				String _bc = bc1.barcode+bc2.barcode;
				chosen_sample = bc2sample.get(_bc);
				
			}
			
			return new String[]{spl1, spl2, chosen_sample};
		}
		
		
	}
	
	/**
	 * Writes the FastqRecord in the case of using index files
	 * 
	 * @param r the read to write, as read from input FASTQ file
	 * @param bcLen the barcode length. In this case, this is the length of an extra barcode found at the read start
	 * @param xtrim how many extra base should be clipped after the barcode ; 0 or more
	 * @param writer the writer
	 * @param barcode the barcode sequence from this read. Note that the caller is responsible to give the 
	 *  sequence as extracted from the read or the barcode sequence after matching the extracted sequence against list of expected barcode
	 * @param addBarcodeToHeader whether the barcode should be added to the read header. If true, the string ':barcode' is added 
	 * with the starting ':' added only if current read header does not end with ':'. For example :
	 * '@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:'
	 * becomes
	 * '@D3FCO8P1:178:C1WLBACXX:7:1101:1836:1965 2:N:0:BARCODE'
	 * @return the trimmed {@link FastqRecord} ; mostly for testing purposes
	 */
	private FastqRecord writeWithIndexFiles(FastqRecord r, int bcLen,
			int xtrim, int ztrim, FastqWriter writer, String barcode,
			boolean addBarcodeToHeader) {
		
		
		int l = bcLen + xtrim ; //safe for bcLen = 0
		log.debug(l+"");
		//create the new header
		String header = r.getReadHeader();
		if(addBarcodeToHeader){
			String extractedBCSeq = bcLen > 0 ? ":"+r.getReadString().substring(0, bcLen).toUpperCase() : "";
			
			header+= (header.endsWith(":")? "":":")+(barcode!=null?barcode:"")+extractedBCSeq;
			if(READ_NAME_REPLACE_CHAR != null){
				header = header.replaceAll(" ", READ_NAME_REPLACE_CHAR);
			}
		}

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
	 * Writes a pair of FastqRecord in the case of using index files making sure that the generated read headers are identical for both reads
	 * 
	 * @param r1 the fwd read to write, as read from input FASTQ file
	 * @param r2 the rev read to write, as read from input FASTQ file
	 * @param bcForRead1 the sample barcode sequence to add to the FWD read header. Note that the caller is responsible to give the 
	 *  sequence as extracted from the index files. This is I1 or I1:I2. 
	 * @param bcForRead2 the sample barcode sequence to add to the REV read header. Note that the caller is responsible to give the 
	 *  sequence as extracted from the index files. This is I1 or I1:I2. 
	 * @param bcLen1 the length of a potential extra barcode in read 1.
	 * @param bcLen2 the length of a potential extra barcode in read 2.
	 * @param xtrim1 how many extra base should be clipped after the barcode from read 1 ; 0 or more
	 * @param xtrim2 how many extra base should be clipped after the barcode from read 2 ; 0 or more
	 * @param ztrim1 how many extra base should be clipped from read 1 end; 0 or more
	 * @param ztrim2 how many extra base should be clipped from read 2 end; 0 or more
	 * @param fwdW the FWD (read1) writer
	 * @param revW the REV (read2) writer
	 * @param addBarcodeToHeader whether the barcode should be added to the read header. 
	 * @return the trimmed reads {@link FastqRecord} as an array with the fwd read in first position; mostly for testing purposes
	 */
	private FastqRecord[] writePEWithIdenticalBarcodesForIndexFiles(
			FastqRecord r1, FastqRecord r2,
			String bcForRead1, String bcForRead2, 
			Integer bcLen1, Integer bcLen2, 
			Integer xtrim1, Integer xtrim2,
			Integer ztrim1, Integer ztrim2, 
			FastqWriter fwdW, FastqWriter revW, 
			boolean addBarcodeToHeader) {
		
		
		/*
		 * create the new header 
		 */
		
		//clip off the '1:N:0'
		String commonHeader = r1.getReadHeader().split("\\s+")[0]; //take any of the two, once '1:N:0' is clipped they are similar
		if(READ_NAME_REPLACE_CHAR == null){
			//then we add back a space
			commonHeader += " ";
		}
		
		//build the extra common header : I1:[I2]:[BC1]:[BC2]
		if(addBarcodeToHeader){
			String commonBC = "";
			
			//make sure bcForRead1 and bcForRead2 are identical 
			if(!bcForRead1.equals(bcForRead2)){
				commonBC = bcForRead1+":"+bcForRead2;
			}else{
				commonBC = bcForRead1;
			}
			commonHeader+= (commonHeader.endsWith(":")? "":":")+(commonBC!=null?commonBC:"");
			//add up the extra barcodes
			String extractedBCSeq1 = bcLen1 > 0 ? ":"+r1.getReadString().substring(0, bcLen1).toUpperCase() : "";
			String extractedBCSeq2 = bcLen2 > 0 ? ":"+r2.getReadString().substring(0, bcLen2).toUpperCase() : "";
			//add it
			commonHeader+=extractedBCSeq1+extractedBCSeq2;
		}

		
		/*
		 * prepare reads for writing
		 */
		int l1 = bcLen1 + xtrim1 ; //safe for bcLen1 = 0
		int l2 = bcLen2 + xtrim2 ; //safe for bcLen2 = 0
		log.debug("l1="+l1+" ; l2="+l2);
		
		
		FastqRecord trimmed1 = prepareTrimmedRead(r1, ztrim1, commonHeader, l1);
		FastqRecord trimmed2 = prepareTrimmedRead(r2, ztrim2, commonHeader, l2);

		if(fwdW != null) //mostly here for testing ie to be able to call method w/o writer
			fwdW.write(trimmed1);
		
		if(revW != null) //mostly here for testing ie to be able to call method w/o writer
			revW.write(trimmed2);
		

		return new FastqRecord[]{trimmed1, trimmed2};
	}

	
	
	@Override
	protected List<String> prepareForEmBASE() {
		//use embase mode not supported with index files yet
		List<String> messages = new ArrayList<String>();
		//TODO : why is this commented ?
		//messages.add("USE_EMBASE mode not currently supported in Illumina demultiplexing using Index files");

		return messages;
	}


	/**
	 * @param args
	 */
	public static void main(final String[] argv) {
		new JemultiplexerIllumina().instanceMainWithExit(argv); //eventually execute the doWork() method
	}
	
}
