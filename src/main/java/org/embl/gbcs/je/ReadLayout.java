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
package org.embl.gbcs.je;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.embl.cg.utilitytools.utils.StringUtil;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Describe the read layout i.e. indicates where barcode, UMI, fixed bases and sample sequences are located. 
 * Barcodes, UMIs and read sequence are specified using a "<BLOCKCODE:length>" while fixed bases are depicted using ATCGN directly.
 * The following BLOCKCODE are supported :    
 * <ul>
 * <li>BARCODE, expects a fixed length e.g. <BARCODEn:6> ; where 'n' is an optional number to uniquely identify this barcode slot </li>
 * <li>UMI, expects a fixed length e.g. <UMIn:8> ; where 'n' is an optional number to uniquely identify this umi slot</li>
 * <li>SAMPLE, expects a fixed length or "x" e.g. <SAMPLEn:30> or <SAMPLE:x> ; where 'n' is an optional number to uniquely identify this sample sequence slot
 * When no length ('x') is specified, all the sequence till the end is considered. 
 * When a length is provided, the exact sequence length is considered (any extra bases are discarded). 
 * Negative values are also supported to indicate all but the last x bases ; this is only accepted when the SAMPLE block is the last one
 * </li>
 * </ul>
 * 
 * <br/>
 * <h2>Note on the optional 'n' in blocks.</h2> 
 * When multiple Barcode slots are described (potentially across different read layout), Je needs to understand how to use these barcodes. 
 * 
 * Let s talk examples : <br />
 * 
 * In standard PE protocol, the barcode is found at the beginning of both reads and is identical between read1 and read2 (REDUNDANT_BARCODES=True).
 * In such a situation, you can simply give READ_LAYOUT=1:<BARCODE:6><SAMPLE:x> and READ_LAYOUT=2:<BARCODE:6><SAMPLE:x> which is 
 * equivalent to READ_LAYOUT=1:<BARCODE1:6><SAMPLE:x> and READ_LAYOUT=2:<BARCODE1:6><SAMPLE:x> or simply READ_LAYOUT=<BARCODE:6><SAMPLE:x>. 
 * <br/>
 * Note that the "1:", "2:", ... "i:" number preceding the read layout. This is used to match the read layout to the correct FASTQ file (passed like FASTQ=1:fastq_1.txt.gz, FASTQ=2:fastq_2.txt.gz, ...)
 * An alternative is to provide the read layout with the FASTQ i.e.  FASTQ=1:<BARCODE1:6><SAMPLE:x>:fastq_1.txt.gz 
 * <br/>
 * When this is not the case and the barcodes are different, all barcodes are then used to look up sample name and the barcode file must
 *  have BARCODEn headers i.e. matching the number of defined BARCODE ; for example if 2 barcode blocks were defined like
 * <br/>
 * READ_LAYOUT=1:<BARCODE1:6><SAMPLE:x> and READ_LAYOUT=2:<BARCODE2:6><SAMPLE:x>
 * <br/>
 * The barcode file format is expected to be like :
 * <br />
 * SAMPLE   BARCODE1   BARCODE2
 * <br />
 * sample1  ATGCGC     TTCGAA
 * sample2  GCGCTA     AACTGA
 * ...
 * 
 * <br/>
 * Important:
 * <ul>
 * <li> The sample sequence must always come after UMI and/or Barcode blocks.</li> 
 * <li>The Layout accepts multiple blocks for Barcode and UMI but not for the sample sequence </li>
 * </ul>
 * <br/>
 * 
 * Example of read layout :
 * <ul>
 *   <li> 
 *   	NN<BARCODE:6>N<UMI:8>A<SAMPLE:x> would instruct Je to ignore the first two bases, use the first 6 bases to match the sample barcode, 
 *      ignore a base (N), use the next 8 bases as a UMI, ignore the next A and use the remaining of the sequence as the sample sequence  
 *   </li> 
 * </ul>
 * 
 * Note that in a paired-end setup, two read layouts should be assembled when the read structures are different ; otherwise a unique layout is sufficient     
 * 
 * @author girardot
 *
 */
public class ReadLayout {
	private static Logger log = LoggerFactory.getLogger(ReadLayout.class);

	
	protected static final String LABEL_FOR_SAMPLE = "Sample";
	protected static final String LABEL_FOR_BARCODE = "Barcode";
	protected static final String LABEL_FOR_UMI = "UMI";

	protected static final String CHAR_FOR_SAMPLE = "X";
	protected static final String CHAR_FOR_BARCODE = "B";
	protected static final String CHAR_FOR_UMI = "I";

	
	public static final String BLOCKCODE_BARCODE = "BARCODE";
	public static final String BLOCKCODE_UMI = "UMI";
	public static final String BLOCKCODE_SAMPLE = "SAMPLE";
	
	protected String layout;
	protected boolean hasUMIBlock;
	protected boolean hasBarcodeBlock;
	protected boolean hasSampleBlock;
	
	protected List<Integer> umiLength;
	protected List<Integer> bcLength;
	//protected List<Integer> sampleSequenceLength; 
	protected Integer sampleSequenceLength; //a unique sample sequence block is currently supported
	
	protected List<Integer> umiStart;
	protected List<Integer> bcStart;
	//protected List<Integer> sampleSequenceStart;
	protected Integer sampleSequenceStart; //a unique sample sequence block is currently supported
	
	/*
	 * BLOCK Id i.e. the number associated to the blocks i.e. is 2 from <SAMPLE2:x> or <UMI2:6> ;
	 *  defaults to '1' when no id is found e.g <SAMPLE:x>
	 *  The maps map the block ID to its position in the layout
	 */
	protected Integer sampleBlockId; 
	protected Map<Integer, Integer> umiBlockId2BlockPosition;
	protected Map<Integer, Integer> bcBlockId2BlockPosition;
	
	public ReadLayout(String layout){
		this.layout = layout;
		parseLayout();
	}
	
	/**
	 * Process the layout and initialize useful variables
	 */
	public void parseLayout(){
		Pattern umiBlockP = Pattern.compile("<"+BLOCKCODE_UMI+"(\\d?):(\\d+)>");
		Pattern bcBlockP = Pattern.compile("<"+BLOCKCODE_BARCODE+"(\\d?):(\\d+)>");
		Pattern readBlockP = Pattern.compile("<"+BLOCKCODE_SAMPLE+"(\\d?):(-?[0-9x]+)>");
		
		String flatlayout = new String(layout); //copy the layout
		log.debug(flatlayout);
		/*
		 * Look for the different blocks, the idea is to replace each block in the layout with a flatlayout ie :
		 * <br />
		 * NN<BARCODE:6>N<UMI:8>A<SAMPLE:x>
		 * <br /> then becomes
		 * NNBBBBBBNIIIIIIIIAXXXXX
		 * <br />
		 * 
		 */
				
		
		// look for UMI blocks
		int i = 0;
		umiLength = new ArrayList<Integer>();
		umiStart = new ArrayList<Integer>();
		umiBlockId2BlockPosition = new HashMap<Integer, Integer>();
		while(flatlayout.contains(BLOCKCODE_UMI)){
			hasUMIBlock = true;
			int l = getBlockLength(umiBlockP, LABEL_FOR_UMI, flatlayout);
			int blockNumber = getBlockNumber(umiBlockP, LABEL_FOR_UMI, flatlayout);
			log.debug("got block len of "+l);
			flatlayout = replaceBlock(flatlayout, umiBlockP, CHAR_FOR_UMI, l);
			umiLength.add( l );
			umiBlockId2BlockPosition.put(blockNumber, i);
			i++;
			log.debug(flatlayout);
		}
		
		// look for BARCODE blocks
		i = 0;
		bcLength = new ArrayList<Integer>();
		bcStart = new ArrayList<Integer>();
		bcBlockId2BlockPosition = new HashMap<Integer, Integer>();
		while(flatlayout.contains(BLOCKCODE_BARCODE)){
			hasBarcodeBlock = true;
			int l = getBlockLength(bcBlockP, LABEL_FOR_BARCODE, flatlayout);
			int blockNumber = getBlockNumber(bcBlockP, LABEL_FOR_BARCODE, flatlayout);
			flatlayout = replaceBlock(flatlayout, bcBlockP, CHAR_FOR_BARCODE, l);
			bcLength.add( l );
			bcBlockId2BlockPosition.put(blockNumber, i);
			i++;
			log.debug(flatlayout);
		}
		
		// look for sample sequence blocks
		if(flatlayout.contains(BLOCKCODE_SAMPLE)){
			hasSampleBlock = true;
			Integer l = getBlockLength(readBlockP, LABEL_FOR_SAMPLE, flatlayout);
			sampleBlockId = getBlockNumber(readBlockP, LABEL_FOR_SAMPLE, flatlayout);
			flatlayout = replaceBlock(flatlayout, readBlockP, CHAR_FOR_SAMPLE, (l==null || l <=0 ? 1:l));
			sampleSequenceLength = l;
			if(l==null || l <=0 ){
				//the block had to be the last one ; check this 
				if(flatlayout.contains("<")){
					throw new ReadLayoutMalformedException(
							"The use of "
									+ (l==null?" the 'x'" : " a negative ("+l+")")
									+" value is only accepted if the "+BLOCKCODE_SAMPLE+" block is the last one; which is not the case. \n"
									+ "Currently flatten to "+flatlayout+")", layout);
				}
				//ok then 
			}
		}
		
		/*
		 * At this step all block replacement should have occurred and the flat layout should only be made of ACTGUNXIB letters
		 */
		log.debug(flatlayout);
		if(!Pattern.matches("[ACTGUN"+CHAR_FOR_BARCODE+CHAR_FOR_SAMPLE+CHAR_FOR_UMI+"]+", flatlayout)){
			throw new ReadLayoutMalformedException(
					"Read layout contains unexpected characters! \n"
							+ "End of layout flattening is : "+flatlayout+")", layout);
		}
		
		if(!hasBarcodeBlock && !hasSampleBlock && !hasUMIBlock){
			throw new ReadLayoutMalformedException("Not a single block found in layout!", layout);
		}
		
		
		//extract start positions and block names now
		if(hasUMIBlock){
			for (int j = 0; j < umiLength.size(); j++) {
				int len = umiLength.get(j);
				String repeat = String.format("%0" + len + "d", 0).replace("0",CHAR_FOR_UMI);
				int sta = flatlayout.indexOf(repeat, (j == 0 ? 0 : umiStart.get(j-1)+umiLength.get(j-1)) );
				log.debug(LABEL_FOR_UMI+" of length "+len+" starts at "+sta);
				umiStart.add( sta );
			}
		}
		
		if(hasBarcodeBlock){
			for (int j = 0; j < bcLength.size(); j++) {
				int len = bcLength.get(j);
				String repeat = String.format("%0" + len + "d", 0).replace("0",CHAR_FOR_BARCODE);
				int sta = flatlayout.indexOf(repeat, (j == 0 ? 0 : bcStart.get(j-1)+bcLength.get(j-1)) );
				log.debug(LABEL_FOR_BARCODE+" of length "+len+" starts at "+sta);
				bcStart.add( sta );
			}
		}
		
		if(hasSampleBlock){
			sampleSequenceStart = flatlayout.indexOf(CHAR_FOR_SAMPLE);
			log.debug(LABEL_FOR_SAMPLE+" of length "+sampleSequenceLength+" starts at "+sampleSequenceStart);
		}
		
		//final check : the sample block must be the last one for layout with a sample block
		if(hasSampleBlock){
			for(int s: umiStart){

				if(s > sampleSequenceStart ){
					throw new ReadLayoutMalformedException("The "+BLOCKCODE_SAMPLE+" block is not the last 3' block of your layout while it must : at least one "+BLOCKCODE_UMI+" block is found after", layout);
				}
			}
			for(int s: bcStart){
				if(s > sampleSequenceStart ){
					throw new ReadLayoutMalformedException("The "+BLOCKCODE_SAMPLE+" block is not the last 3' block of your layout while it must : at least one "+BLOCKCODE_BARCODE+" block is found after", layout);
				}
			}
		}
	}
	
	private String replaceBlock(String flatlayout, Pattern pat,
			String charToRepeat, Integer len) {
		String repeat = String.format("%0" + len + "d", 0).replace("0",charToRepeat);
		Matcher m = pat.matcher(flatlayout);
		flatlayout = m.replaceFirst(repeat);
		return flatlayout;
	}

	/**
	 * @param pat the pattern 
	 * @param blockName for error reporting only, use a user meaningful name here 
	 * @return null if block has unlimited length
	 * @throws ReadLayoutMalformedException is pattern can't be matched or is invalid
	 */
	private Integer getBlockLength(Pattern pat, String blockName, String flatlayout) {
		Matcher m = pat.matcher(flatlayout);
		if(!m.find()){
			throw new ReadLayoutMalformedException("Malformed read layout : block for "+blockName+" could not be found ! ", layout);
		}
		
		Integer l = null;  
		String token = m.group(2);
		try {
			l = Integer.parseInt(token);
		} catch (NumberFormatException e) {
			// if it is 'x' it is fine
			if( !token.equalsIgnoreCase("x")){
				//then it is a format error
				throw new ReadLayoutMalformedException("Malformed read block : length in block for "+blockName+" should be specified with a valid number or 'x'! ", layout);
			}
		}
		return l;
		
	}
	
	/**
	 * @param pat the pattern 
	 * @param blockName for error reporting only, use a user meaningful name here 
	 * @return the number found next to the bllock name i.e. 2 from <UMI2:6>; or 1 if block has no number
	 * @throws ReadLayoutMalformedException is pattern can't be matched or is invalid
	 */
	private Integer getBlockNumber(Pattern pat, String blockName, String flatlayout) {
		Matcher m = pat.matcher(flatlayout);
		if(!m.find()){
			throw new ReadLayoutMalformedException("Malformed read layout : block for "+blockName+" could not be found ! ", layout);
		}
		
		Integer l = null;  
		String token = m.group(1);
		if(token == null || token.isEmpty())
			return 1;
		
		try {
			l = Integer.parseInt(token);
		} catch (NumberFormatException e) {
			//then it is a format error
			throw new ReadLayoutMalformedException("Malformed read block : block number for "+blockName+" should be specified with a valid number or absent ", layout);
		}
		return l;
		
	}

	/**
	 * @return true if this layout has a UMI block
	 */
	public boolean containsUMI(){
		return hasUMIBlock;
	}
	
	/**
	 * @return true if this layout has a Barcode block
	 */
	public boolean containsBarcode(){
		return hasBarcodeBlock;
	}
	
	/**
	 * @return true if this layout has a Sample block
	 */
	public boolean containsSampleSequence(){
		return hasSampleBlock;
	}
	
	
	/**
	 * @return the number of UMI block number found in this layout
	 */
	public int umiBlockNumber(){
		if(!hasUMIBlock)
			return 0;
		
		return umiStart.size();
	}
	
	
	/**
	 * @return the number of BARCODE block number found in this layout
	 */
	public int barcodeBlockNumber(){
		if(!hasBarcodeBlock)
			return 0;
		
		return bcStart.size();
	}
	
	/**
	 * Extract the subsequence(s) corresponding to the UMI blocks in the read layout
	 * 
	 * @param read the whole read or quality string
	 * @return the sequences corresponding to the UMI or null if this layout has no UMI block
	 */
	public String[] extractUMIs(String read){
		if(!hasUMIBlock)
			return null;
		
		String [] umis = new String[umiStart.size()];
		for (int i = 0; i < umiStart.size(); i++) {
			int start = umiStart.get(i);
			int len = umiLength.get(i);
			umis[i] = read.substring(start, start+len);
		}
		
		return umis;
	}
	
	/**
	 * Extract the subsequence(s) corresponding to the UMI blocks in the read layout
	 * and merge them in a unique String (following the 5' to 3' order on the layout)
	 * 
	 * @param read the whole read or quality string
	 * @return the sequences corresponding to the UMI or null if this layout has no UMI block
	 */
	public String extractUMI(String read){
		if(!hasUMIBlock)
			return null;
		
		return StringUtil.mergeArray( extractUMIs(read) , "" );
	}
	
	/**
	 * Extract the subsequence(s) corresponding to a particular UMI block identified by its ID i.e. '3' from <UMI3:6> 
	 * in the read layout . 
	 * 
	 * @param read the whole read or quality string
	 * @param umiBlockId the UMI block id i.e. '3' from a layout like <UMI3:6> ; this is NOT the 5' to 3' 
	 * order number in which UMI blocks were described in the layout. 
	 * For example, the <UMI3:6> slot is the "first" UMI block in the layout "<BARCODE2:6><UMI3:6><UMI2:6><SAMPLE:x>"
	 * but is identified with the umiBlock ID 3 (i.e. UMI3). 
	 * Remember that block Id applies to all provided layout i.e. if you described twice UMI3 in different layouts 
	 * given to Je ; Je understands that these two blocks should contain the same sequence in a read set     
	 * @return the sequences corresponding to the UMI named with umiBlock ID or null if this layout has no UMI block or no UMI identified with this number
	 *   
	 */
	public String extractUMI(String read, int umiBlockId){
		if(!hasUMIBlock ){
			log.error("NO UMI block in this layout : "+this.layout);
			return null;
		}
		if(!umiBlockId2BlockPosition.containsKey(umiBlockId)){
			log.error("NO UMI block with ID "+umiBlockId +" in this layout : "+this.layout);
			return null;
		}
		
		int i = umiBlockId2BlockPosition.get(umiBlockId);
		int start = umiStart.get(i);
		int len = umiLength.get(i);
		return read.substring(start, start+len);
	}
	
	/**
	 * Extract the subsequence(s) corresponding to a particular BARCODE block identified by its Id i.e. '3' from <BARCODE3:6> 
	 * in the read layout . 
	 * 
	 * @param read the whole read or quality string
	 * @param bcBlockId the BARCODE block id i.e. '3' from a layout like <BARCODE3:6> ; this is NOT the 5' to 3' 
	 * order number in which BARCODE blocks were described in the layout. 
	 * For example, the <BARCODE3:6> slot is the "first" BARCODE block in the layout "<UMI2:6><BARCODE3:6><BARCODE2:6><SAMPLE:x>"
	 * but is identified with the ID 3 (i.e. BARCODE3). 
	 * Remember that block Id applies to all provided layout i.e. if you described twice BARCODE3 in different layouts 
	 * given to Je ; Je understands that these two blocks should contain the same sequence in a read set     
	 * @return the sequences corresponding to the BARCODE with ID or null if this layout has no BARCODE block or no BARCODE identified with this number
	 *   
	 */
	public String extractBarcode(String read, int bcBlockId){
		if(!hasBarcodeBlock ){
			log.error("NO BARCODE block in this layout : "+this.layout);
			return null;
		}
		if(!bcBlockId2BlockPosition.containsKey(bcBlockId)){
			log.error("NO BARCODE block with ID "+bcBlockId +" in this layout : "+this.layout);
			return null;
		}
		
		int i = bcBlockId2BlockPosition.get(bcBlockId);
		int start = bcStart.get(i);
		int len = bcLength.get(i);
		return read.substring(start, start+len);
	}
	
	
	/**
	 * 
	 * Extract the subsequence(s) corresponding to the BARCODE blocks in the read layout
	 * 
	 * @param read the whole read or quality string
	 * @return the barcode subsequences or null if this layout has no BARCODE block
	 */
	public String[] extractBarcodes(String read){
		if(!hasBarcodeBlock)
			return null;
		String [] bcs = new String[bcStart.size()];
		for (int i = 0; i < bcStart.size(); i++) {
			int start = bcStart.get(i);
			int len = bcLength.get(i);
			bcs[i] = read.substring(start, start+len);
		}
		
		return bcs;
	}
	
	
	/**
	 * Extract the subsequence(s) corresponding to the BARCODE blocks in the read layout
	 * and merge them in a unique String (following the 5' to 3' order on the layout)
	 * 
	 * @param read the whole read or quality string
	 * @return the sequences corresponding to the BARCODE or null if this layout has no BARCODE block
	 */
	public String extractBarcode(String read){
		if(!hasBarcodeBlock)
			return null;
		return StringUtil.mergeArray( extractBarcodes(read) , "" );
	}
	
	
	/**
	 * @param read the whole read or quality string
	 * @return the read subsequence or null if this layout has no SAMPLE block
	 */
	public String extractSampleSequence(String read){
		if(!hasSampleBlock)
			return null;
		
		if(sampleSequenceLength == null){
			//take all 
			return read.substring(sampleSequenceStart);
		}
		else if(sampleSequenceLength>0){
			return read.substring(sampleSequenceStart, sampleSequenceStart+sampleSequenceLength);
		}
		else{
			return read.substring(sampleSequenceStart, read.length()+sampleSequenceLength);
		}
	}
	
}
