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

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import htsjdk.samtools.util.StringUtil;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;

public class ReadEndsForMarkDuplicatesWithMolecularCode2 extends
		ReadEndsForMarkDuplicates implements MolecularCodedReadEnds{

	/*
    What do we need to store you ask?  Well, we need to store:
       - byte: orientation
       - short: libraryId, readGroup, tile, x, y, score
       - int: read1ReferenceIndex, read1Coordinate, read2ReferenceIndex, read2Coordinate, code length
       - long: read1IndexInFile, read2IndexInFile
       - byte array for code : fix it to 20, should be a max here
     */
	public static final int SIZE_OF = (1 * 1) + (5 * 2) + (4 * 5) + (8 * 2) + 1
            + 8 + // last 8 == reference overhead
            + 20 + //arbitrary byte array for code
            13; // This is determined experimentally with JProfiler
	
	protected byte[] molecularCodeArr  = null;
	
	public void setMolecularCode(String code){
		molecularCodeArr = StringUtil.stringToBytes(code.toUpperCase());
	}
	
	/**
	 * This is only meant to be used by the Codec. We don t want people to use this ie 
	 * we want to make sure they use the string version so that we ensure codes are set upper case 
	 * @param molecularCodeArr
	 */
	protected void setMolecularCode(byte[] molecularCodeArr){
		this.molecularCodeArr = molecularCodeArr;
	}

	/**
	 * @return the molecularCode
	 */
	public String getMolecularCode() {
		return StringUtil.bytesToString(molecularCodeArr);
	}
	
	/**
	 * @param codelen length of an individual code. Each code is supposed to have the same length
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

	public String toString(){
		String s = "RG="+readGroup+" ; "+tile+":"+x+":"+y+" ; orientation="+orientation+" ; paired="+(isPaired()?"1":"0")+" \n "
				+"   R1-ref="+read1ReferenceIndex + ", R1-coor="+read1Coordinate +", idx="+read1IndexInFile+" \n "
				+"   R2-ref="+read2ReferenceIndex + ", R2-coor="+read2Coordinate +", idx="+read2IndexInFile+" \n ";
				
		return s;
	}
	
}
