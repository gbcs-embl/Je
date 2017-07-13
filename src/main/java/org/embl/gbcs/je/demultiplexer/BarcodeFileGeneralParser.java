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
package org.embl.gbcs.je.demultiplexer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.embl.cg.utilitytools.utils.CollectionUtils;
import org.embl.cg.utilitytools.utils.FileUtil;
import org.embl.cg.utilitytools.utils.StringUtil;
import org.embl.cg.utilitytools.utils.parser.csv.CSVLine;
import org.embl.cg.utilitytools.utils.parser.csv.InvalidHeaderException;
import org.embl.gbcs.je.Jexception;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * In this format, the file structure is governed with headers:<br />
 * <ol>
 * <li>the mandatory'SAMPLE' column lists the sample names</li>
 * <li>the mandatory 'BARCODEn' columns list the matching BARCODE from the BARCODEn slot (where n is a number). It is mandatory
 *  to have as many 'BARCODEn' columns as described BARCODE slots. 
 *  Here again, barcodes can be combined using the OR operator '|' 
 *  </li>
 *  <li>the optional 'OUTn' columns (where n is a number) list the output file names for this sample and matching output number. 
 *  When used, one must define as many OUTn columns as defined OUTPUT_LAYOUT on the command line </li>
 * 
 * @author girardot
 *
 */
public class BarcodeFileGeneralParser  {

	protected static final String HEADER_OUT = "OUT";

	protected static final String HEADER_BARCODE = "BARCODE";

	protected static final String HEADER_SAMPLE = "SAMPLE";

	private static Logger log = LoggerFactory.getLogger(BarcodeFileGeneralParser.class);

	/**
	 * regex of the header line
	 */
	public static String headerLineRegex = "^" + HEADER_SAMPLE + "(\\s"+HEADER_BARCODE+"\\d+)+(\\s"+HEADER_OUT+"\\d+)*$";
	
	/**
	 * The barcode file to validate
	 */
	protected File file;
	
	
	/**
	 * Maps the sample names to their barcode combination i.e. each sample is mapped to one or more barcodes 
	 * (hence the use of a Set) from each BARCODE slot (hence the use of a List which size is the number of BARCODE
	 *  columns in the barcode file)
	 */
	protected Map<String, List<Set<String>>> sample2BarcodeSets = new HashMap<String, List<Set<String>>>();
	
	/**
	 * length for the different barcode i.e. each BARCODE column contains barcodes of identical length which length
	 * can differ between 2 BARCODE columns.
	 * List order follows the BARCODE column idx ; which is the 'n' from the BARCODEn headers or the order in file
	 */
	protected List<Integer> barcodeLengths ;
	
	
	/**
	 * A logical indicating if the barcodes have Ns in their sequence (indicating degenerate position)
	 * A boolean is needed for each BARCODE column
	 */
	protected List<Boolean> barcodesHaveNs ;
	
	/**
	 * Indicates if some barcodes are defined as a mix i.e. ACTAC|ATGTC|GTCAG
	 */
	protected Boolean foundMixedBarcodes = false;
	
	/**
	 * Sample name to result file list which order follows the 'n' from OUTn headers (or the order in file) 
	 * populated only if OUT headers are present 
	 */
	protected Map<String, List<String>> sample2outfiles ;
	
	
	
	
	private Pattern headerPattern = Pattern.compile("^("+HEADER_SAMPLE+"|"+HEADER_BARCODE+"|"+HEADER_OUT+")(\\d*)$");
	private Matcher headerMatcher = headerPattern.matcher("");	
	private boolean hasOutHeaders = false;
	private Integer maxBarcodeIndex = null;
	private Integer maxOutIndex = null;

	
	public BarcodeFileGeneralParser(File file) {
		if(file == null)
			throw new Jexception("Barcode file is null");
		
		if(!file.canRead())
			throw new Jexception("Barcode file "+file.getAbsolutePath()+" cannot be read (rights issue ?).");
		
		this.file = file;
	}

	
	
	public void parse() throws InvalidHeaderException {
		try {
			for( CSVLine l : FileUtil.loadLinesFromCSVHeadedFile(file, "\t", true)){
				
				/*
				 * Is this the header line ?
				 */
				if(l.isHeaderLine()){
					if(!isValidHeaderLine(l.merge("\t"))){
						throw new InvalidHeaderException("Headers in the barcode file do not comply to specifications (please check manual): "+l.merge("\t"));
					}
					//we initialise few things
					boolean foundSampleH= false;
					boolean foundBarcodeH= false;
					Set<Integer> barcodeIndices = new TreeSet<Integer>();
					Set<Integer> outIndices = new TreeSet<Integer>();
					
					for(String h : l.getHeaders()){
						headerMatcher.reset(h);
						if(!headerMatcher.matches()){
							throw new InvalidHeaderException("Found unexpected header in barcode file (please check manual): "+h);
						}
						String t = headerMatcher.group(1);						
						if(t.equalsIgnoreCase(HEADER_SAMPLE)){
							if(foundSampleH)
								throw new InvalidHeaderException("Found twice the header "+h+" in barcode file. Please correct");
							foundSampleH = true;
						}
						else if(t.equalsIgnoreCase(HEADER_BARCODE)){
							try{
								int n = Integer.parseInt(headerMatcher.group(2));
								foundBarcodeH = true;
								if(barcodeIndices.contains(n))
									throw new InvalidHeaderException("Found twice the header "+h+" in barcode file. Please correct.");
								barcodeIndices.add(n);
							}catch(NumberFormatException nfe){
								throw new InvalidHeaderException("Found "+HEADER_BARCODE+" without mandatory index in barcode file: "+ h);
							}
						}
						else if(t.equalsIgnoreCase(HEADER_OUT)){
							try{
								int n = Integer.parseInt(headerMatcher.group(2));
								hasOutHeaders = true;
								if(outIndices.contains(n))
									throw new InvalidHeaderException("Found twice the header "+h+" in barcode file. Please correct.");

								outIndices.add(n);
							}catch(NumberFormatException nfe){
								throw new InvalidHeaderException("Found "+HEADER_OUT+" without mandatory index in barcode file: "+ h);
							}
							
						}
						else{
							//is it even possible to get here ...
							throw new InvalidHeaderException("Invalid header found in barcode file : "+t);
						}
					}
					
					/*
					 *  headers look good, few more checks:
					 *  # did we get all mandatory headers
					 *  # check that the  BARCODE and OUT headers are number like 1,2.3 ...n
					 */
					if(!foundSampleH)
						throw new InvalidHeaderException("Did not find mandatory header : "+HEADER_SAMPLE);
					if(!foundBarcodeH)
						throw new InvalidHeaderException("Did not find mandatory header : "+HEADER_BARCODE);
					
					this.maxBarcodeIndex = checkHeaderNumbering(HEADER_BARCODE, barcodeIndices);
					// init 
					Integer[] integers = new Integer[this.maxBarcodeIndex];
					Arrays.fill(integers, 0);
					this.barcodeLengths = Arrays.asList(integers); 
					Boolean[] bools = new Boolean[this.maxBarcodeIndex];
					Arrays.fill(bools, false);
					this.barcodesHaveNs = Arrays.asList(bools); 
					if(hasOutHeaders){
						this.maxOutIndex = checkHeaderNumbering(HEADER_OUT, outIndices);
					}
					continue;
				}
				
				/*
				 * process lines
				 */
				String smpl = l.getValue(HEADER_SAMPLE);
				for (int i = 1 ; i <= maxBarcodeIndex ; i++) {
					checkBarcodeLengthAndRegsiter(smpl, l.getValue(HEADER_BARCODE+i), i);
				}
				
				if(this.hasOutHeaders){
					if(sample2outfiles == null){
						sample2outfiles = new HashMap<String, List<String>>();
					}
					if(!sample2outfiles.containsKey(smpl)){
						sample2outfiles.put(smpl, new ArrayList<String>(this.maxOutIndex));
					}
					
					for (int i = 1 ; i <= maxOutIndex ; i++) {
						String fpath = l.getValue(HEADER_OUT+i); //this can be just a file name
						if(StringUtils.isBlank(fpath))
							throw new Jexception("no value for sample "+smpl+" in "+HEADER_OUT+i);
						sample2outfiles.get(smpl).set((i-1), fpath);
					}
				}
				
			}
		} catch (IOException e) {
			new RuntimeException(e);
		} 
		
	}



	/**
	 * @param l
	 * @return
	 */
	public static boolean isValidHeaderLine(String line) {
		return Pattern.matches(BarcodeFileGeneralParser.headerLineRegex, line);
	}



	/**
	 * @param sample the sample name for this line
	 * @param barcodes can be a single or multiple barcodes separated with |
	 * @param i
	 */
	private void checkBarcodeLengthAndRegsiter(String sample, String barcodes, int i) {
		Set<String> allPossiblesBarcodes = new TreeSet<String>();
		if(!foundMixedBarcodes)
			foundMixedBarcodes = barcodes.contains("|");
		
		for(String b: barcodes.split("\\|")){
			//init if necessary
			if(barcodeLengths.get(i - 1) == 0){
				barcodeLengths.set((i - 1), b.length());
			}
			
			//check length is same as before
			if(b.length() != barcodeLengths.get(i - 1)){
				throw new Jexception("Error in barcode file : found barcode of length "+b.length()+" ("+b+") in column "
						+HEADER_BARCODE+i+" while previous barcodes were "+barcodeLengths.get(i - 1)+" bases long. All barcodes must have same length in a given column.");
			}
			
			//check barcode only contains allowed letters 
			if(!Pattern.matches("^[ATGCUNatugcn]+$", b)){
				throw new Jexception("Barcode "+b+" contains inavlid letters (in "+HEADER_BARCODE+i+" column)");
			}
			
			//save N status
			if(false == barcodesHaveNs.get(i-1) ){
				barcodesHaveNs.set((i-1), containsN(b));
			}
			
			//register
			allPossiblesBarcodes.add(b);
		}
		
		//save
		if(!sample2BarcodeSets.containsKey(sample)){
			//create
			 List<Set<String>> _l = new ArrayList<Set<String>>();
			 for (int j = 0; j < this.maxBarcodeIndex; j++) {
				_l.add(new TreeSet<String>());
			}
			sample2BarcodeSets.put(sample, _l);
		}
		sample2BarcodeSets.get(sample).set((i-1), allPossiblesBarcodes);
	}



	/**
	 * @param b
	 * @return
	 */
	protected static boolean containsN(String b) {
		return Pattern.matches("[Nn]", b);
	}



	/**
	 * check that the indices form a series n(i+1) = n(i) + 1 
	 * @param headertype
	 * @param ids
	 * @return return the max index
	 * @throws InvalidHeaderException if the indices do not form a series n(i+1) = n(i) + 1 
	 */
	protected static int checkHeaderNumbering(String headertype, Set<Integer> ids) throws InvalidHeaderException {
		Integer min = null;
		Integer max = null;
		for (Integer i : ids) {
			if(min == null){
				min = i; 
				max = i;
				continue;
			}
			if(i>max) max = i;
			if(i<min) min = i;
		}
		// starts at 1 ?
		if(min != 1)
			throw new InvalidHeaderException("Header for "+headertype+" must start with number 1 (not "+min+")");
		
		// start at expected number accordign to id size ?
		int expectedEnd = min + ids.size() - 1;
		if(max != expectedEnd){
			ArrayList<Integer> _ids = new ArrayList<Integer>(ids);
			Collections.sort(_ids);
			throw new InvalidHeaderException("Indices for the "+ids.size()+" "+headertype
					+" headers do not form a continous sequence of integers from 1 to "+expectedEnd+" : "+StringUtil.mergeIntList(_ids,  ","));
		}
		
		return max;
	}




	/**
	 * @return the file
	 */
	public File getFile() {
		return file;
	}


	/**
	 * @return the sample2BarcodeSets
	 */
	public Map<String, List<Set<String>>> getSample2BarcodeSets() {
		return sample2BarcodeSets;
	}


	/**
	 * @return the barcodeLengths
	 */
	public List<Integer> getBarcodeLengths() {
		return barcodeLengths;
	}


	/**
	 * @return the barcodesHaveNs
	 */
	public List<Boolean> getBarcodesHaveNs() {
		return barcodesHaveNs;
	}


	/**
	 * @return the sample2outfiles or null if not OUT headers were found
	 */
	public Map<String, List<String>> getSample2outfiles() {
		return sample2outfiles;
	}


	/**
	 * @param file the file to set
	 */
	public void setFile(File file) {
		this.file = file;
	}


	/**
	 * @param sample2BarcodeSets the sample2BarcodeSets to set
	 */
	public void setSample2BarcodeSets(
			Map<String, List<Set<String>>> sample2BarcodeSets) {
		this.sample2BarcodeSets = sample2BarcodeSets;
	}


	/**
	 * @param barcodeLengths the barcodeLengths to set
	 */
	public void setBarcodeLengths(List<Integer> barcodeLengths) {
		this.barcodeLengths = barcodeLengths;
	}


	/**
	 * @param barcodesHaveNs the barcodesHaveNs to set
	 */
	public void setBarcodesHaveNs(List<Boolean> barcodesHaveNs) {
		this.barcodesHaveNs = barcodesHaveNs;
	}


	/**
	 * @param sample2outfiles the sample2outfiles to set
	 */
	public void setSample2outfiles(Map<String, List<String>> sample2outfiles) {
		this.sample2outfiles = sample2outfiles;
	}



	public boolean useBarcodeMix() {
		return foundMixedBarcodes;
	}



	public int getBarcodeColumnNumber() {
		return this.barcodeLengths.size();
	}


}
