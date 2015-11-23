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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.embl.gbcs.je.Jexception;

/**
 *
 * This class is responsible to extract the molecular barcode(s) from the read headers.
 * @author girardot
 *
 */
public class MolecularBarcodeFinder {

	
	protected static final Pattern UMI_VALIDATION_PATTERN = Pattern.compile("[ATGCNatgcn]+");
	protected static Matcher UMI_VALIDATION_MATCHER = null;
	
	//init
	static{
		UMI_VALIDATION_MATCHER = UMI_VALIDATION_PATTERN.matcher("");
	}
	
	/**
	 * length of the different UMIs, this is initialized using the first read and later used to validate
	 * all reads. An error is launched if different length are observed
	 * Let to null ir null is used to launch initialisation 
	 */
	protected static Integer[] UMI_LENGTH = null;
	
	
	private Integer [] slots = new Integer[0];
	private String splitter  = ":";
	private boolean HAS_UMIS = true;
	
	public MolecularBarcodeFinder( String splitter, Integer[] slots){
		this.splitter = splitter;
		if(this.slots==null){
			HAS_UMIS = false;
		}else if(this.slots.length==1 && this.slots[0] == 0){
			//slot=0 means no slots
			HAS_UMIS = false;
		}else{
			this.slots = slots;
		}
	}
	
	
	public void addMolecularBarcodes(String readName, ReadEndsForMarkDuplicatesWithMolecularCode ends){
		
		if(!HAS_UMIS){
			UMI_LENGTH = new Integer[]{0};
			return; // do nothing
		}
		
		String[] tokens = readName.split(splitter);
		StringBuffer molCode = new StringBuffer();
		boolean init_umi_len = (UMI_LENGTH == null);
		if(init_umi_len){
			UMI_LENGTH= new Integer[slots.length];
		}
		int umi_len_cnt = 0;
		for (int slot : slots) {
			String umi = null;
			if(slot>=0){
				umi = tokens[slot];
			}else{
				//from end
				umi = tokens[tokens.length + slot];
			}
			//init if required 
			if(init_umi_len){
				UMI_LENGTH[umi_len_cnt]  = umi.length();
			}else if(UMI_LENGTH[umi_len_cnt] != umi.length()){
				//length not valide
				throw new Jexception("Invalid UMI length in slot "+slot+" : found a UMI with length "+umi.length()+" while previous UMIs from this slot had length of "+UMI_LENGTH[umi_len_cnt]);
			}
			//validate UMI 			
			UMI_VALIDATION_MATCHER.reset(umi);
			if(!UMI_VALIDATION_MATCHER.matches()){
				throw new Jexception("Invalid format for UMI in slot "+slot+" : found UMI '"+umi+"' while expecting UMI matching pattern "+UMI_VALIDATION_PATTERN.toString());
			}
			molCode.append(umi);
			umi_len_cnt++;
		}
		ends.setMolecularCode( molCode.toString() );
	}
	
	/**
	 * Helper method to make sure that all UMI slots contain a UMI of expected length
	 * @param readName
	 * @param expectedLength
	 * @return false if one of the slots does not contain sequence of expected length
	 */
	public boolean validateUMILength(String readName, int expectedLength){
		
		if(!HAS_UMIS){
			return true;
		}
		
		String[] tokens = readName.split(splitter);

		boolean isvalid = true;
		for (int slot : slots) {
			String umi = null;
			if(slot>=0){
				umi = tokens[slot]; 
			}else{
				//from end
				umi = tokens[tokens.length + slot];
			}
			if(umi.length() != expectedLength){
				isvalid = false;
				break;
			}
		}
		return isvalid;
	}
	
	
}
