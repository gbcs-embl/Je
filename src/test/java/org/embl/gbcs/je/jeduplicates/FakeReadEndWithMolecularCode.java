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

import htsjdk.samtools.util.StringUtil;

import java.util.ArrayList;
import java.util.List;

import org.embl.gbcs.je.jeduplicates.MolecularCodedReadEnds;

public class FakeReadEndWithMolecularCode implements MolecularCodedReadEnds {

	protected byte[] molecularCodeArr  = null;
	protected Integer uniqcode = null; // a code with potenial base errors (including N)
	public String originalcode = null; // the code this barcode shoudl match to 

	public FakeReadEndWithMolecularCode(int uniqcode, String originalcode, String code){
		this.uniqcode = uniqcode;
		this.originalcode = originalcode;
		molecularCodeArr = StringUtil.stringToBytes(code);
	}


	public void setMolecularCode(String code){
		molecularCodeArr = StringUtil.stringToBytes(code);
	}
	
	public void setMolecularCode(byte[] molecularCodeArr){
		this.molecularCodeArr = molecularCodeArr;
	}

	/**
	 * @return the molecularCode
	 */
	public String getMolecularCode() {
		return StringUtil.bytesToString(molecularCodeArr);
	}
	
	/**
	 * @param codelen length of an individual code. Each code is supposed to have the length
	 * @return the molecularCodes
	 */
	public List<String> getMolecularCodes(int codelen) {
		int codenumber = molecularCodeArr.length / codelen;
		List<String> codes = new ArrayList<String>(codenumber);
		for(int i = 0 ; i < codenumber; i++){
			codes.add(
					StringUtil.bytesToString(molecularCodeArr, i*codelen, codelen)
					);
		}
		return codes;
	}
	

	/**
	 * @return the molecularCodeArr
	 */
	public byte[] getMolecularCodeArr() {
		return molecularCodeArr;
	}


}
