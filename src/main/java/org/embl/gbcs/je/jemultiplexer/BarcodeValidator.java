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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.embl.cg.utilitytools.utils.ExceptionUtil;
import org.embl.cg.utilitytools.utils.FileUtil;
import org.embl.cg.utilitytools.utils.parser.csv.CSVLine;
import org.embl.cg.utilitytools.utils.parser.csv.CSVParser;
import org.embl.cg.utilitytools.utils.parser.csv.InvalidCSVSetUpException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * 
 * Validator for barcode file
 * 
 * @author girardot
 *
 */
public class BarcodeValidator {

	private static Logger log = LoggerFactory.getLogger(BarcodeValidator.class);

	/**
	 * The barcode file to validate
	 */
	protected File file;
	
	/**
	 * The demultiplexer implementation calling this validator 
	 */
	protected AbstractJemultiplexer demultiplexer;
	
	/**
	 * Set after validation call.
	 * If true, the file also defines output file names
	 */
	public Boolean isLongFormat = null; 
	
	/**
	 * indicates if a barcode for each read is described (PE mode)
	 * set according to file parsing 
	 * only set in PE situation
	 */
	public Boolean hasBarcodesForBothEnds = null; 
	
	/**
	 * length for barcode 1
	 */
	public Integer bclen1 = null;
	
	/**
	 * length for barcode 2, if applicable
	 */
	public Integer bclen2 = null;
	
	/**
	 * A logical indicating if the barcodes have Ns in their sequence (indicating degenerate position)
	 * this happens in very rare occasions like iCLIP protocol
	 */
	public Boolean barcodesHaveNs = false;
	
	/**
	 * barcode to sample mapping, for each read end. Only relevant ones will be populated 
	 */
	public Map<String, String> barcode2sample = new HashMap<String, String>();
	
	/**
	 * sample name to barcode list
	 */
	public Map<String, Set<String>> sample2barcodes = new HashMap<String, Set<String>>();
	
	/**
	 * sample name to result file name for read1 
	 * populated only if isLongFormat=true 
	 */
	public Map<String, String> sample2file1 = new HashMap<String, String>();
	
	/**
	 * sample name to result file name for read2, in PE mode 
	 * populated only if isLongFormat=true 
	 */
	public Map<String, String> sample2file2 = new HashMap<String, String>();
	
	/**
	 * barcode sets that can be found in read 1
	 */
	protected Set<String> barcodeSetRead1 = null;
	
	/**
	 * barcode sets that can be found in read 2, if applicable
	 */
	protected Set<String> barcodeSetRead2 = null;
	
	/**
	 * are Ns found in barcodes found in read 1 
	 */
	protected Boolean areNsInBarcodeSetRead1 = null;
	
	/**
	 * are Ns found in barcodes found in read 2, if applicable
	 */
	protected Boolean areNsInBarcodeSetRead2 = null;
	
	/**
	 * @param file barcode fiel to validate
	 * @param demultiplexer will be used to access command line arguments
	 */
	public BarcodeValidator(File file, AbstractJemultiplexer demultiplexer){
		this.file = file;
		this.demultiplexer = demultiplexer;
		//when using idx files, the BCLEN, if given, relates to extra barcodes found in reads ie not to barcodes described in the bc file
		if(!demultiplexer.USE_SAMPLE_INDEX_FILES){
			this.bclen1 = demultiplexer.BCLEN_1;
			this.bclen2 = demultiplexer.BCLEN_2;
		}
	}
	
	
	/**
	 * @param barcodePosition one of {@link BarcodePosition#READ_1} or {@link BarcodePosition#READ_2}
	 * @return null is barcodeSetRead1 or barcodeSetRead2 have not been initialized
	 */
	public Boolean areNsInBarcodes(BarcodePosition barcodePosition) {
		if(barcodePosition == BarcodePosition.READ_1){
			if(barcodeSetRead1 == null)
				return null;
			
			if(areNsInBarcodeSetRead1 == null){
				areNsInBarcodeSetRead1 = false;
				for(String bc : barcodeSetRead1){
					if(bc.contains("N")){
						areNsInBarcodeSetRead1 = true;
						break;
					}
				}
			}
			return areNsInBarcodeSetRead1;
		}
		else if(barcodePosition == BarcodePosition.READ_2){
			if(barcodeSetRead2 == null)
				return null;
			
			if(areNsInBarcodeSetRead2 == null){
				areNsInBarcodeSetRead2 = false;
				for(String bc : barcodeSetRead2){
					if(bc.contains("N")){
						areNsInBarcodeSetRead2 = true;
						break;
					}
				}
			}
			return areNsInBarcodeSetRead2;
		}
		return null;
	}
	
	
	/**
	 * Parse the barcode file after successful validation of the file.
	 * Should be called ONLY if validate() did not returned errors messages
	 * (behavior not guaranteed otherwise) 
	 */
	public void parse(){
		
		/*
		 * Load basic barcode => sample name mapping in which barcode can be a complex format
		 */
		HashMap<String, String> bc2sample = new HashMap<String, String>();
		try {
			HashMap<String, String> _bc2sample = FileUtil.getMapFromFile(file.getAbsolutePath(), 0, 1, 0, "\\s+", false);
			//filter, make sure no space are in sample and barcode are uppercase
			for (Entry<String, String> e : _bc2sample.entrySet()) {
				bc2sample.put(e.getKey().toUpperCase(), e.getValue().replaceAll("\\s+", "_"));
			}
		} catch (IOException e) {
			log.error("IO error while reading barcode file "+file.getAbsolutePath());
			log.error(ExceptionUtil.getStackTrace(e));
			throw new RuntimeException(e);
		}catch (Exception e) {
			//should not happen as file was checked at parsing time
			log.error(ExceptionUtil.getStackTrace(e));
			throw new RuntimeException(e);
		} 
		
		/*
		 * Compute final maps
		 */
		
		if(
				demultiplexer.RUNNING_PAIRED_END 
				&& (demultiplexer.getBarcodeReadPosOptionValue() == BarcodePosition.BOTH || (demultiplexer.USE_SAMPLE_INDEX_FILES && ((JemultiplexerIllumina)demultiplexer).INDEX_FILE2 != null) ) 
				&& demultiplexer.getRedundantBarcodesOptionValue() == false
				){
			/*
			 * barcodes might define a read separator':' if both barcodes should be used for look up
			 */
			if(demultiplexer.getBarcodeForSampleMatchingOptionValue() == BarcodePosition.BOTH || (demultiplexer.USE_SAMPLE_INDEX_FILES && ((JemultiplexerIllumina)demultiplexer).INDEX_FILE2 != null)){
				
				initBarcodeSeparateMaps(bc2sample);
			}else{
				//only one bc is for sample lookup while the other one is eg a random sequence
				initBarcodeMaps(bc2sample);
			}
		}
		else{
			initBarcodeMaps(bc2sample);
		}
		
	}
	
	


	/**
	 * This method is called in the special situation where both ends barcode should be used to resolve the sample.
	 * i.e. when RUNNING_PAIRED_END=true BARCODE_READ_POS=BOTH REDUNDANT_BARCODES=false
	 * In this situation, both barcode are concatenated (not joining character) BC_READ1 + BC_READ2 
	 * e.g. ATAGCA:GCGGCA => ATAGCAGCGGCA and used this unique longer barcode will be used to populate maps 
	 * @param bc2sample the "barcode => sample" name mapping in which barcode is not parsed (i.e. as read from file)
	 *  and defines barcode(s) for both ends i.e. in a format like bc1[|bc]*:bc2[|bc]*
	 */
	private void initBarcodeSeparateMaps(Map<String, String> bc2sample) {
		
		barcodeSetRead1 = new TreeSet<String>();
		barcodeSetRead2 = new TreeSet<String>();
		
		Map<String, String> m = new HashMap<String, String>();
		for (Entry<String, String> e : bc2sample.entrySet()) {
			String sample = e.getValue();
			String [] bcArr = e.getKey().split(":");
			//should not happen but...
			if(bcArr.length!=2)
				throw new RuntimeException("Splitting "+e.getKey()+" by ':' returned "+bcArr.length+" tokens instead of 2. At this step, barcode read from file are supposed to comply to the bc1[|bc]*:bc2[|bc]* format.");
		
			String [] bc1arr = bcArr[0].split("\\|");
			String [] bc2arr = bcArr[1].split("\\|");
			
			for (String bc1 : bc1arr) {
				barcodeSetRead1.add(bc1);
				for (String bc2 : bc2arr) {
					barcodeSetRead2.add(bc2);
					m.put(bc1+bc2, sample);
				}
			}
		}
		
		//first the bc key might be of the format bc1|bc2 => split it and explode the map
		barcode2sample = new HashMap<String, String>();
		for (Entry<String, String> e : m.entrySet()) {
			for(String bc : e.getKey().split("\\|")){
				barcode2sample.put(bc, e.getValue());
				log.debug("PUT in barcode2sample : '"+bc+"' => '"+e.getValue()+"'");
			}
		}

		//build sample2barcodes maps
		sample2barcodes = new HashMap<String, Set<String>>();
		for (Entry<String, String> e : barcode2sample.entrySet()) {
			if(!sample2barcodes.containsKey(e.getValue()))
				sample2barcodes.put(e.getValue(), new TreeSet<String>());
			sample2barcodes.get(e.getValue()).add(e.getKey());
		}

		initSampleFileMaps();
	}


	/**
	 * This method is called in all simple situation ie a single barcode (list) is used to lookup sample
	 *  @param bc2sample the "barcode => sample" name mapping in which barcode is not parsed (i.e. as read from file)
	 *  but complies to a simple situation format like bc1[|bc]*
	 */
	private void initBarcodeMaps(Map<String, String> bc2sample) {
		//note both mapps are init with SAME content to ease further processing
		barcodeSetRead1 = new TreeSet<String>();
		barcodeSetRead2 = new TreeSet<String>();
		
		//first the bc key might be of the format bc1|bc2 => split it and explode the map
		barcode2sample = new HashMap<String, String>();
		for (Entry<String, String> e : bc2sample.entrySet()) {
			for(String bc : e.getKey().split("\\|")){
				//there is still the option that one bc is fully random and described as NNN:realbarcode or realbarcode:NNN
				if(bc.contains(":")){
					int realBCIdx = demultiplexer.getBarcodeForSampleMatchingOptionValue() == BarcodePosition.READ_1 ? 0:1;
					bc = bc.split(":")[realBCIdx]; //this will drop the NNNN part
				}
				//populate both sets
				barcodeSetRead1.add(bc);
				barcodeSetRead2.add(bc);
				
				barcode2sample.put(bc, e.getValue());
			}
		}
		
		//build sample2barcodes maps
		sample2barcodes = new HashMap<String, Set<String>>();
		for (Entry<String, String> e : barcode2sample.entrySet()) {
			if(!sample2barcodes.containsKey(e.getValue()))
				sample2barcodes.put(e.getValue(), new TreeSet<String>());
			sample2barcodes.get(e.getValue()).add(e.getKey());
		}

		initSampleFileMaps();

	}

	private void initSampleFileMaps(){
		//load sample file mapping
		if(isLongFormat){
			try {
				HashMap<String, String> _sample2file = FileUtil.getMapFromFile(file.getAbsolutePath(), 0, 0, 2, "\\s+", false);
				sample2file1 = new HashMap<String, String>();
				//filter, make sure no space are in sample and file names
				for (Entry<String, String> e : _sample2file.entrySet()) {
					sample2file1.put(e.getKey().replaceAll("\\s+", "_"), e.getValue().replaceAll("\\s+", "_"));
				}
				//second file ?
				if(demultiplexer.RUNNING_PAIRED_END){
					_sample2file = FileUtil.getMapFromFile(file.getAbsolutePath(), 0, 0, 3, "\\s+", false);
					sample2file2 = new HashMap<String, String>();
					//filter, make sure no space are in sample and file names
					for (Entry<String, String> e : _sample2file.entrySet()) {
						sample2file2.put(e.getKey().replaceAll("\\s+", "_"), e.getValue().replaceAll("\\s+", "_"));
					}
				}

			} catch (Exception e) {
				log.error(ExceptionUtil.getStackTrace(e));
				throw new RuntimeException(e);
			}
		}
	}
	

	/**
	 * @return a list of error messages, which is empty if no pbs were found
	 */
	public List<String> validate(){
		List<String> messages = new ArrayList<String>();
		if(demultiplexer.RUNNING_PAIRED_END){
			validatePE(messages);
		}else{
			validateSE(messages);
		}
		
		return messages;
	}



	private void validateSE(List<String> messages) {
		isLongFormat = null;
		/*
		 * read barcodes and check they are unique and of adequate length
		 */
		Set<String> bcSet = new TreeSet<String>();
		CSVParser p = null;
		try {
			p = new CSVParser(-1,-1,false);
			Iterator<CSVLine> lines = p.iterator(file.getAbsolutePath(), "\\s+");
			while (lines.hasNext()) {
				CSVLine line = lines.next();
				log.debug("barcode file has "+line.getColNum() +" columns");
				if(isLongFormat == null){
					if(line.getColNum()>3 || line.getColNum()<2)
						messages.add("Invalid barcode file format : barcode file for single end must have 2 or 3 columns, not "+line.getColNum());
					
					isLongFormat = line.getColNum() == 3;
				}
				String barcodeLine = line.getValue(1);
				if(!checkBarcodeValueSyntaxSE(barcodeLine, messages)){
					continue;
				}
				
				//ok 
				String[] bcs = barcodeLine.split("\\|");
				for(String barcode : bcs){
					//check bc is unique
					if(bcSet.contains(barcode))
						messages.add("Barcode found multiple times :"+barcode);

					int l = barcode.length();
					//should we init BCLEN_1?
					if(bclen1==null || bclen1 == 0){
						log.info("Setting barcode length with first barcode length read from barcode file : "+l);
						bclen1 = l;
					}

					if(l!=bclen1)
						messages.add("Barcode "+barcode+" does not have the expected size of "+bclen1);
					
					bcSet.add(barcode);
				}
			}

			//do we use long format?
			if(isLongFormat){
				//yes try read to make sure no parsing error occurs
				try {
					//get mapping sample => file name for read 1
					FileUtil.getMapFromFile(file.getAbsolutePath(), 0, 0, 2, "\\s+", false);
				} catch (Exception e) {
					log.error(ExceptionUtil.getStackTrace(e));
					messages.add("Barcode file with long format ("+file.getAbsolutePath()+") is malformed : "+e.getMessage());
				}
			}

		} catch (InvalidCSVSetUpException e) {
			log.error(ExceptionUtil.getStackTrace(e));
			messages.add("Barcode file is malformed: "+e.getMessage());
		} catch (IOException e) {
			log.error(ExceptionUtil.getStackTrace(e));
			messages.add("IO Error while barcode file "+file.getAbsolutePath()+" : "+e.getMessage());
		} finally{
			if(p!=null)
				p.terminate();
		}
	}



	private void validatePE(List<String> messages) {
		isLongFormat = null;
		/*
		 * read barcodes and check they are unique and of adequate length
		 */
		CSVParser p = null;
		try {
			p = new CSVParser(-1,-1,false);
			Iterator<CSVLine> lines = p.iterator(file.getAbsolutePath(), "\\s+");
			hasBarcodesForBothEnds = null; //according to file parsing
			
			Set<String> bcSet1 = new TreeSet<String>();
			Set<String> bcSet2 = new TreeSet<String>();
			//init with cmd line attr for the end that shoudl be used for smpl lookup
			Integer _bclen = demultiplexer.getBarcodeForSampleMatchingOptionValue() == BarcodePosition.READ_2 ? bclen2 : bclen1;
			
			while (lines.hasNext()) {
				CSVLine line = lines.next();
				
				if(isLongFormat == null){
					if(line.getColNum()!=2 && line.getColNum()!=4)
						messages.add("Invalid barcode file format : barcode file for paired end must have 2 or 4 columns, not "+line.getColNum());
					
					isLongFormat = line.getColNum() == 4;
				}
				String barcodeLine = line.getValue(1);
				log.debug(barcodeLine);
				if(!checkBarcodeValueSyntaxPE(barcodeLine, messages)){ //this checks compliance of having a ':' with BARCODE_READ_POS option value
					log.warn(barcodeLine+" does not look like a valid barcode expression");
					continue;
				}
				String [] barcodeStringArr = barcodeLine.split(":");
				//check format coherence across the file 
				if(hasBarcodesForBothEnds == null){
					hasBarcodesForBothEnds = barcodeStringArr.length == 2;
				}
				if((barcodeStringArr.length == 2) != hasBarcodesForBothEnds){
					//this check is still needed as checkBarcodeValueSyntaxPE does not check if lines are coherent altogether (one with : then one without...)
					if(hasBarcodesForBothEnds)
						messages.add("Barcode value :"+barcodeLine+" does not define barcode for both read ends as in previous lines (no ':' found) ; in line "+line.merge("\t"));
					else
						messages.add("Barcode value :"+barcodeLine+" defines barcode for both read ends while it was not the case in previous lines (':' found) ; in line "+line.merge("\t"));
					
					continue;
				}
				//ok , check barcodes
				if(hasBarcodesForBothEnds){
					for (int i = 0; i < barcodeStringArr.length; i++) {
						String[] bcs = barcodeStringArr[i].split("\\|");
						Set<String> bcSet = i==0 ? bcSet1 : bcSet2;
						for(String barcode : bcs){
							//do not check bc is unique unless one of them is not for sample lookup OR barcodes are redundant
							if(!demultiplexer.USE_SAMPLE_INDEX_FILES && 
									(demultiplexer.getBarcodeForSampleMatchingOptionValue()!=BarcodePosition.BOTH || demultiplexer.getRedundantBarcodesOptionValue() )
									){
								if(bcSet.contains(barcode)){
									//this is a problem only when this is the barcode for sample matching 
									if( 
											(i == 0 && demultiplexer.getBarcodeForSampleMatchingOptionValue()==BarcodePosition.READ_2) ||  
											(i == 1 && demultiplexer.getBarcodeForSampleMatchingOptionValue()==BarcodePosition.READ_1)  
									){
										log.debug("IGNORING (bc not used for sample matching): Barcode found multiple times : '"+barcode+"' for read_"+(i+1)+" end.");
									}else{
										messages.add("Barcode found multiple times : '"+barcode+"' for read_"+(i+1)+" end.");
									}
								}
							}
							int l = barcode.length();
							//should we init bclen_*?
							Integer blen = i==0 ? bclen1 : bclen2;
							if(blen==null || blen == 0){
								log.info("Setting barcode length for read_"+(i+1)+" end with first barcode length read from barcode file : "+l);
								blen = l;
								if(i==0)
									bclen1 = l;
								else
									bclen2 = l;
							}

							 
							if(l!=blen)
								messages.add("Barcode "+barcode+" does not have the expected size of "+blen+" for read_"+(i+1)+" end.");

							bcSet.add(barcode);
						}
					}
					
					if(demultiplexer.getRedundantBarcodesOptionValue() == true){
						//highly suspect !
						if(!demultiplexer.USE_SAMPLE_INDEX_FILES)
							log.warn("POSSIBLE ERROR : using different barcodes at both ends while REDUNDANT_BARCODES is set to TRUE !! ");
						else
							log.warn("POSSIBLE ERROR : using different barcodes in index files while REDUNDANT_BARCODES is set to TRUE !! ");
					}
					
				}else{
					//work with bcSet1 and bclen1 but adapt later on
					String[] bcs = barcodeStringArr[0].split("\\|");
					log.info("Init barcode len with "+_bclen);
					for(String barcode : bcs){
						log.info(barcode);
						//check bc is unique
						if(bcSet1.contains(barcode))
							messages.add("Barcode found multiple times :"+barcode);

						int l = barcode.length();
						//should we init bclen_*?
						if(_bclen==null || _bclen == 0){
							log.info("Setting barcode length with first barcode length read from barcode file : "+l);
							_bclen = l;
						}

						
						if(l!=_bclen){
							messages.add("Barcode "+barcode+" does not have the expected size of "+_bclen);
						}

						bcSet1.add(barcode);
					}
					bclen1=_bclen;
				}
				
			}
			// all data is now in the bclen1 and bcSet1 => should we reverse ie if it is READ_2 
			boolean reverseToread2 = false;
			/*
			 * note that the null testing here is meant for cases where no line of the barcode file was valid => the hasBarcodesForBothEnds
			 * could therefore not be set
			 */
			if(hasBarcodesForBothEnds!=null && !hasBarcodesForBothEnds){
				if(demultiplexer.getBarcodeReadPosOptionValue() == BarcodePosition.BOTH){
					if(demultiplexer.getBarcodeForSampleMatchingOptionValue() == BarcodePosition.BOTH ){
						//copy values to read 2 as well
						bcSet2 = bcSet1;
						bclen2 = bclen1;
					}else if (demultiplexer.getBarcodeForSampleMatchingOptionValue() == BarcodePosition.READ_2){
						reverseToread2 = true;
					}
				}else if (demultiplexer.getBarcodeReadPosOptionValue() ==  BarcodePosition.READ_2){
					reverseToread2 = true;
				}
			}
			//should we put value into read slots ?
			if(reverseToread2){
				bcSet2 = new TreeSet<String>(bcSet1);
				bclen2 = bclen1.intValue();
				//cancel read 1 values unless there is a rdm bc
				bclen1= (demultiplexer.getBarcodeReadPosOptionValue() == BarcodePosition.BOTH && demultiplexer.getBarcodeForSampleMatchingOptionValue() == BarcodePosition.READ_2)? demultiplexer.BCLEN_1: 0;
				bcSet1 = new TreeSet<String>();
			}

			//do we use long format?
			if(isLongFormat){
				//yes try read to make sure no parsing error occurs
				try {
					//get mapping sample => file name for read 1
					FileUtil.getMapFromFile(file.getAbsolutePath(), 0, 0, 2, "\\s+", false);
					//get mapping sample => file name for read 2
					FileUtil.getMapFromFile(file.getAbsolutePath(), 0, 0, 3, "\\s+", false);
				} catch (Exception e) {
					log.error(ExceptionUtil.getStackTrace(e));
					messages.add("Barcode file with long format ("+file.getAbsolutePath()+") is malformed : "+e.getMessage());
				}
			}

		} catch (InvalidCSVSetUpException e) {
			log.error(ExceptionUtil.getStackTrace(e));
			messages.add("Barcode file is malformed: "+e.getMessage());
		} catch (IOException e) {
			log.error(ExceptionUtil.getStackTrace(e));
			messages.add("IO Error while barcode file "+file.getAbsolutePath()+" : "+e.getMessage());
		} finally{
			if(p!=null)
				p.terminate();
		}
	}
	

	/**
	 * Check that the value given for barcode complies to the syntax bc[|bc2]*[:bc[|bc2]*]* ie
	 * bc or bc|bc2|bc3 or bc:bc2 or bc|bc2|bc3:bc4|bc5 where ':' delimits barcode for read_1
	 *  (left part) and read_2 (right part). And '|' delimits different possible barcodes for 
	 *  the same read end (multiple barcodes resolving to the same sample) 
	 * 
	 * @param barcode the  barcode character sequence to validate
	 * @param messages the list of error messages, any error message should be added in this list
	 * @return true if looks ok or false if format does not comply
	 */
	private boolean checkBarcodeValueSyntaxSE(String barcode, List<String> messages){
		if(barcode.contains(":")){
			messages.add("Invalid value used for barcode : "+barcode +". The X:Z synthax can only be used for paired-end data");
			return false;
		}
		
		Pattern p  = Pattern.compile("^[NATGCnatgc\\|]+$");
		Matcher m = p.matcher(barcode);
		if(!m.matches()){
			messages.add("Invalid value used for barcode value : "+barcode +". Barcodes can only contains ATGCN letters and the '|' separator (case insensitive) ");
			return false;
		}
		if(!barcodesHaveNs && barcode.toUpperCase().contains("N"))
			barcodesHaveNs = true;
		
		return true;
	}
	
	/**
	 * Check that the value given for barcode complies to the syntax bc[|bc2]*[:bc[|bc2]*]* ie
	 * bc or bc|bc2|bc3 or bc:bc2 or bc|bc2|bc3:bc4|bc5 where ':' delimits barcode for read_1
	 *  (left part) and read_2 (right part). And '|' delimits different possible barcodes for 
	 *  the same read end (multiple barcodes resolving to the same sample) 
	 * 
	 * @param barcode the  barcode character sequence to validate
	 * @param messages the list of error messages, any error message should be added in this list
	 * @return true if looks ok or false if format does not comply
	 */
	private boolean checkBarcodeValueSyntaxPE(String barcode, List<String> messages){
		String [] values = new String[2]; 
		if(barcode.contains(":")){
			if(
					(demultiplexer.getBarcodeReadPosOptionValue()!=BarcodePosition.BOTH && !demultiplexer.USE_SAMPLE_INDEX_FILES)
					|| ( demultiplexer.USE_SAMPLE_INDEX_FILES && ((JemultiplexerIllumina)demultiplexer).INDEX_FILE2 == null )
					){
				if(demultiplexer.USE_SAMPLE_INDEX_FILES)
					messages.add("Invalid value used for barcode : "+barcode +". The X:Z synthax can only be used for paired-end data with two Illumina INDEX FILEs");
				else
					messages.add("Invalid value used for barcode : "+barcode +". The X:Z synthax can only be used for paired-end data with BPOS=BOTH");
				return false;
			}else if (demultiplexer.getBarcodeForSampleMatchingOptionValue()!=BarcodePosition.BOTH && !demultiplexer.USE_SAMPLE_INDEX_FILES){
				//barcodes for both ends should only be given if both ends are used for sample lookup ; but it is not necessarily a mistake
				log.warn("Possible mistake with barcode expression "+barcode+" : It does not make great sense to provide barcodes for both ends while BARCODE_FOR_SAMPLE_MATCHING="+demultiplexer.getBarcodeForSampleMatchingOptionValue());
			}
			//split values
			values = barcode.split(":");
			if(values.length!=2){
				messages.add("Invalid value used for barcode value : "+barcode +". The X:Z synthax can contain only a single ':' ");
				return false;
			}
		}else{
			if(
					demultiplexer.getBarcodeForSampleMatchingOptionValue()==BarcodePosition.BOTH && demultiplexer.getRedundantBarcodesOptionValue() == false 
					&& 
					( (demultiplexer.getBarcodeReadPosOptionValue()==BarcodePosition.BOTH && !demultiplexer.USE_SAMPLE_INDEX_FILES) || (demultiplexer.USE_SAMPLE_INDEX_FILES && ((JemultiplexerIllumina)demultiplexer).INDEX_FILE2 != null ) )
					){
				//problem : use must give separate barcodes in such situation
				String _mess = "Invalid barcode expression "+barcode+" : Expecting barcodes for BOTH ends (syntax bc1:bc2) when using BARCODE_READ_POS=BOTH BARCODE_FOR_SAMPLE_MATCHING=BOTH REDUNDANT_BARCODES=false.";
				messages.add(_mess);
				log.debug(_mess);
				return false;
			}
			values = new String[]{barcode}; 
		}
		
		Pattern p  = Pattern.compile("^[NATGCnatgc\\|]+$");
		for (String s : values) {
			Matcher m = p.matcher(s);
			if(!m.matches()){
				log.debug(s+" is NOT a valid barcode expression");
				messages.add("Invalid value used for barcode value : "+barcode +". Barcodes can only contains NATGC letters and the '|' separator (case insensitive) ");
				return false;
			}else{
				log.debug(s+" is valid barcode expression");
				if(!barcodesHaveNs && s.toUpperCase().contains("N"))
					barcodesHaveNs = true;
				
			}
		}
		
		return true;
	}


	/**
	 * 
	 * Checks that the bclen1 and bclen2 extracted from BC file match the size of sequences in index file(s)
	 * 
	 * @param indexFile1
	 * @param indexFile2
	 * @return
	 */
	public List<String> checkBarcodesLengthWithIndexFiles(
			File indexFile1, File indexFile2) {
		
		List<String> messages = new ArrayList<String>();
		
		int sizeInFile = getRealReadLength(indexFile1);
		if(sizeInFile != bclen1){
			messages.add("Index reads in file "+indexFile1.getAbsolutePath()+" are "+sizeInFile+" bases long while barcodes described in barcode file are "+bclen1+" bases long. These two must be of same size !");
		}
		
		if(indexFile2!=null){
			sizeInFile = getRealReadLength(indexFile2);
			int _bclen2 = (bclen2 == null || bclen2 == 0) ? bclen1: bclen2; //if BC are redundant in the index files, then bclen2 might not have been initialized
			if(sizeInFile != _bclen2){
				messages.add("Index reads in file "+indexFile2.getAbsolutePath()+" are "+sizeInFile+" bases long while barcodes described in barcode file are "+_bclen2+" bases long. These two must be of same size !");
			}
		}

		return messages;
	}


	protected int getRealReadLength(File fastqFile) {
    	
		final FastqReader r = new FastqReader(fastqFile);
		
		try {
			FastqRecord rec = r.iterator().next();
			return rec.length();
		} finally {
			if(r!=null)
				r.close();
				
		}
    	
	}

	
}
