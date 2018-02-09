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
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.embl.cg.utilitytools.utils.StringUtil;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Describe the read layout i.e. indicates where barcode, UMI, fixed bases and sample sequences are located. 
 * Barcodes, UMIs and read sequence are specified using a "<BLOCKCODEn:length>" while fixed bases are depicted using ATCGN directly.
 * The BLOCKCODE 'BARCODE', 'UMI' and 'SAMPLE' are supported with :  
 * <ol>
 * <li> 'n' is a number to uniquely identify this barcode slot ; please starts numbering at 1. Note that Je expects 
 * continuous series of indices for each block (1,2,3 ...n). 'n' can be omitted only if a unique block of a particular
 *  BLOCKCODE is found across all read layout</li>
 * <li> Obvioulsy, when a length is provided, the exact sequence length is considered (any extra bases are discarded) </li> 
 * <li> When no length ('x') is specified, all the sequence till the end is considered ; it only possible to use the 'x' 
 * shortcut in the last block of a layout</li>
 * <li> When a negative value is given in place of length (e.g. '<BLOCKCODEn:-2>'), all but the last x (2 in the 
 * '<BLOCKCODEn:-2>' example) bases ; a negative length value is only accepted in the last block of a layout</li>
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

	protected static final String CHAR_FOR_SAMPLE = "S";
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
	protected List<Integer> sampleSequenceLength; 
	
	
	protected List<Integer> umiStart;
	protected List<Integer> bcStart;
	protected List<Integer> sampleSequenceStart;
	
	
	/*
	 * BLOCK Id i.e. the number associated to the blocks i.e. is 2 from <SAMPLE2:x> or <UMI2:6> ;
	 *  defaults to '1' when no id is found e.g <SAMPLE:x>
	 *  The maps map the block ID to its position in the layout
	 */ 
	protected Map<Integer, Integer> sampleBlockId2BlockPosition;
	protected Map<Integer, Integer> umiBlockId2BlockPosition;
	protected Map<Integer, Integer> bcBlockId2BlockPosition;
	protected List<Integer> sampleBlockIdOrdered;
	protected List<Integer> umiBlockIdOrdered;
	protected List<Integer> bcBlockIdOrdered;
	
	protected Pattern blockP = Pattern.compile("<("+BLOCKCODE_SAMPLE+"|"+BLOCKCODE_BARCODE+"|"+BLOCKCODE_UMI+")(\\d?):(-?[0-9x]+)>");
	protected String layoutValidationRegex = "^([a-zA-Z]*<("+BLOCKCODE_SAMPLE+"|"+BLOCKCODE_BARCODE+"|"+BLOCKCODE_UMI+")(\\d?):(-?[0-9x]+)>[a-zA-Z]*)+$";
	
	
	public ReadLayout(String layout){
		this.layout = layout;
		parseLayout();
	}
	
	/**
	 * change
	 * Process the layout and initialize useful variables
	 */
	public void parseLayout(){
		
		if(!Pattern.matches(layoutValidationRegex, layout)){
			throw new LayoutMalformedException("Layout does not match validation regex '"+layoutValidationRegex+"'", layout);
		}
		
		String flatlayout = new String(layout); //copy the layout
		log.debug(flatlayout);
		
		/*
		 * Look for the different blocks, the idea is to replace each block in the layout with a flatlayout ie :
		 * <br />
		 * NN<BARCODE:6>N<UMI:8>A<SAMPLE:x>
		 * <br /> then becomes
		 * NNBBBBBBNIIIIIIIIASSSSS
		 * <br />
		 * 
		 */		
		
		// init all 
		int umiCnt = 0;
		umiLength = new ArrayList<Integer>();
		umiStart = new ArrayList<Integer>();
		umiBlockId2BlockPosition = new HashMap<Integer, Integer>();
		umiBlockIdOrdered = new ArrayList<Integer>();
		int smpCnt = 0;
		sampleSequenceLength = new ArrayList<Integer>();
		sampleSequenceStart = new ArrayList<Integer>();
		sampleBlockId2BlockPosition = new HashMap<Integer, Integer>();
		sampleBlockIdOrdered = new ArrayList<Integer>();
		int bcCnt = 0;
		bcLength = new ArrayList<Integer>();
		bcStart = new ArrayList<Integer>();
		bcBlockId2BlockPosition = new HashMap<Integer, Integer>();
		bcBlockIdOrdered = new ArrayList<Integer>();
		
		Matcher layoutMatcher = blockP.matcher(layout);
		boolean previousLengthWasIncomplete = false;
		while(layoutMatcher.find()){
			//we get a new block mathc
			if(previousLengthWasIncomplete){
				throw new LayoutMalformedException(
						"Undefined (i.e. 'x') or negative block length are only allowed in the last block.", layout);
			}
			String blockType = layoutMatcher.group(1);
			int blockNumber = getBlockNumber(blockType, layoutMatcher.group(2));
			Integer l = getBlockLength(blockType, layoutMatcher.group(3));
			if(l == null || l < 0){
				previousLengthWasIncomplete = true;
			}
			String wholeMatch = layoutMatcher.group(0);
			log.debug(
					wholeMatch+ " => \n" + 
							"  Type = "+blockType+"\n" +
							"  Idx = '"+blockNumber+"'\n" +
							"  Len = "+l+"\n" 
							);
			if(blockType.equals(BLOCKCODE_BARCODE)){
				hasBarcodeBlock = true;
				bcLength.add( l );
				bcBlockId2BlockPosition.put(blockNumber, bcCnt);
				bcBlockIdOrdered.add(blockNumber);
				bcCnt++;
				flatlayout = replaceBlock(flatlayout, wholeMatch, CHAR_FOR_BARCODE, (l==null || l <=0 ? 1:l));
			}
			else if(blockType.equals(BLOCKCODE_UMI)){
				hasUMIBlock = true;
				umiLength.add( l );
				umiBlockId2BlockPosition.put(blockNumber, umiCnt);
				umiBlockIdOrdered.add(blockNumber);
				umiCnt++;
				flatlayout = replaceBlock(flatlayout, wholeMatch, CHAR_FOR_UMI, (l==null || l <=0 ? 1:l));
			}
			else if(blockType.equals(BLOCKCODE_SAMPLE)){
				hasSampleBlock = true;
				sampleSequenceLength.add( l );
				sampleBlockId2BlockPosition.put(blockNumber, smpCnt);
				sampleBlockIdOrdered.add(blockNumber);
				smpCnt++;
				flatlayout = replaceBlock(flatlayout, wholeMatch, CHAR_FOR_SAMPLE, (l==null || l <=0 ? 1:l));
			}
			else{
				throw new Jexception("Unknown block type in read layout : "+blockType);
			}
			log.debug(flatlayout);
		}
	
		
		/*
		 * At this step all block replacement should have occurred and the flat layout should only be made of ACTGUNXIB letters
		 */
		log.debug(flatlayout);
		if(!Pattern.matches("[ACTGUN"+CHAR_FOR_BARCODE+CHAR_FOR_SAMPLE+CHAR_FOR_UMI+"]+", flatlayout)){
			throw new LayoutMalformedException(
					"Read layout contains unexpected characters! \n"
							+ "End of layout flattening is : "+flatlayout+")", layout);
		}
		
		if(!hasBarcodeBlock && !hasSampleBlock && !hasUMIBlock){
			throw new LayoutMalformedException("Not a single block found in layout!", layout);
		}
		
		
		//extract start positions and block names now
		if(hasUMIBlock){
			for (int j = 0; j < umiLength.size(); j++) {
				Integer len = umiLength.get(j);
				int sta = nextStartInLayout(flatlayout, CHAR_FOR_UMI, len, (j == 0 ? 0 : umiStart.get(j-1)+umiLength.get(j-1)) );
				log.debug(LABEL_FOR_UMI+" of length "+len+" starts at "+sta);
				umiStart.add( sta );
			}
		}
		
		if(hasBarcodeBlock){
			for (int j = 0; j < bcLength.size(); j++) {
				Integer len = bcLength.get(j);
				int sta = nextStartInLayout(flatlayout, CHAR_FOR_BARCODE, len, (j == 0 ? 0 : bcStart.get(j-1)+bcLength.get(j-1)) );
				log.debug(LABEL_FOR_BARCODE+" of length "+len+" starts at "+sta);
				bcStart.add( sta );
			}
		}
		
		if(hasSampleBlock){
			for (int j = 0; j < sampleSequenceLength.size(); j++) {
				Integer len = sampleSequenceLength.get(j);
				int sta = nextStartInLayout(flatlayout, CHAR_FOR_SAMPLE, len, (j == 0 ? 0 : sampleSequenceStart.get(j-1)+sampleSequenceLength.get(j-1)) );
				log.debug(LABEL_FOR_SAMPLE+" of length "+len+" starts at "+sta);
				sampleSequenceStart.add( sta );
			}
		}

	}

	private int nextStartInLayout(String flatlayout, String blockChar,
			Integer len, int from) {
		int l = (len!=null && len > 0 ? len : 1);
		String repeat = String.format("%0" + l + "d", 0).replace("0",blockChar);
		return flatlayout.indexOf(repeat, from );
	}

	/**
	 * @param flatlayout
	 * @param l
	 */
	protected void validateBlockLengthWithLayout(String flatlayout, Integer l) {
		if(l==null || l <=0 ){
			//the block had to be the last one ; check this 
			if(flatlayout.contains("<")){
				throw new LayoutMalformedException(
						"The use of "
								+ (l==null?" the 'x'" : " a negative ("+l+")")
								+" value is only accepted if the block is the last one; which is not the case. \n"
								+ "Currently flatten to "+flatlayout+")", layout);
			}
			//ok then 
		}
	}
	

	
	private String replaceBlock(String flatlayout, String toReplace,
			String charToRepeat, Integer len) {
		String repeat = String.format("%0" + len + "d", 0).replace("0",charToRepeat);
		return flatlayout.replace(toReplace, repeat);
	}

	/**
	 * @param blockName for error reporting only, use a user meaningful name here 
	 * @param token the string extracted from layout for the length 
	 * @return null if block has unlimited length
	 * @throws LayoutMalformedException is pattern can't be matched or is invalid
	 */
	private Integer getBlockLength(String blockName, String token) {	
		Integer l = null;  
		try {
			l = Integer.parseInt(token);
		} catch (NumberFormatException e) {
			// if it is 'x' it is fine
			if( !token.equalsIgnoreCase("x")){
				//then it is a format error
				throw new LayoutMalformedException("Malformed read block : length in block for "+blockName+" should be specified with a valid number or 'x'! ", layout);
			}
		}
		return l;
		
	}
	
	/**
	 * @param token the string extracted from layout for the block index 
	 * @param blockName for error reporting only, use a user meaningful name here 
	 * @return the number found next to the block name i.e. 2 from <UMI2:6>; or 1 if block has no number
	 * @throws LayoutMalformedException is pattern can't be matched or is invalid
	 */
	private Integer getBlockNumber(String blockName, String token) {
		Integer l = null;  	
		if(token == null || token.isEmpty())
			return 1;
		
		try {
			l = Integer.parseInt(token);
		} catch (NumberFormatException e) {
			//then it is a format error
			throw new LayoutMalformedException("Malformed read block : block number for "+blockName+" should be specified with a valid number or absent ", layout);
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
	 * @return the unique set of SAMPLE ids or an empty set 
	 */
	public Set<Integer> getSampleBlockUniqueIds(){
		if(!containsSampleSequence())
			return new TreeSet<Integer>();
		return this.sampleBlockId2BlockPosition.keySet();
	}
	
	
	/**
	 * @return the SAMPLE ids in the order they are found (from 5 to 3')
	 */
	public List<Integer> getOrderedSampleBlockUniqueIds(){
		if(!containsBarcode())
			return new ArrayList<Integer>();
		
		return this.sampleBlockIdOrdered;
	}
	
	/**
	 * @return the unique set of BARCODE ids or an empty set 
	 */
	public Set<Integer> getBarcodeBlockUniqueIds(){
		if(!containsBarcode())
			return new TreeSet<Integer>();
		return this.bcBlockId2BlockPosition.keySet();
	}
	
	/**
	 * @return the BARCODE ids in the order they are found (from 5 to 3')
	 */
	public List<Integer> getOrderedBarcodeBlockUniqueIds(){
		if(!containsBarcode())
			return new ArrayList<Integer>();
		
		return this.bcBlockIdOrdered;
	}
	
	
	/**
	 * @return the unique set of UMI ids or an empty set 
	 */
	public Set<Integer> getUMIBlockUniqueIds(){
		if(!containsUMI())
			return new TreeSet<Integer>();
		return this.umiBlockId2BlockPosition.keySet();
	}
	
	/**
	 * @return the UMI ids in the order they are found (from 5 to 3')
	 */
	public List<Integer> getOrderedUMIBlockUniqueIds(){
		if(!containsUMI())
			return new ArrayList<Integer>();
		
		return this.umiBlockIdOrdered;
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
	 * @return the number of SAMPLE block number found in this layout
	 */
	public int sampleBlockNumber(){
		if(!hasSampleBlock)
			return 0;
		
		return sampleSequenceStart.size();
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
			Integer len = umiLength.get(i);
			umis[i] = extract(read, start, len);
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
		Integer len = umiLength.get(i);
		return extract(read, start, len);
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
		Integer len = bcLength.get(i);
		return extract(read, start, len);
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
			Integer len = bcLength.get(i);
			bcs[i] = extract(read, start, len);
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
	 * Extract the subsequence(s) corresponding to a particular SAMPLE block identified by its Id i.e. '3' from <SAMPLE3:6> 
	 * in the read layout . 
	 * 
	 * @param read the whole read or quality string
	 * @param bcBlockId the SAMPLE block id i.e. '3' from a layout like <SAMPLE3:6> ; this is NOT the 5' to 3' 
	 * order number in which SAMPLE blocks were described in the layout. 
	 * For example, the <SAMPLE3:6> slot is the "first" BARCODE block in the layout "<UMI2:6><SAMPLE3:6><BARCODE2:6><SAMPLE:x>"
	 * but is identified with the ID 3 (i.e. SAMPLE3). 
	 * Remember that block Id applies to all provided layout i.e. if you described twice SAMPLE3 in different layouts 
	 * given to Je ; Je understands that these two blocks should contain the same sequence in a read set     
	 * @return the sequences corresponding to the SAMPLE with ID or null if this layout has no SAMPLE block or no SAMPLE identified with this number
	 *   
	 */
	public String extractSample(String read, int bcBlockId){
		if(!hasSampleBlock ){
			log.error("NO SAMPLE block in this layout : "+this.layout);
			return null;
		}
		if(!sampleBlockId2BlockPosition.containsKey(bcBlockId)){
			log.error("NO SAMPLE block with ID "+bcBlockId +" in this layout : "+this.layout);
			return null;
		}
		
		int i = sampleBlockId2BlockPosition.get(bcBlockId);
		int start = sampleSequenceStart.get(i);
		Integer len = sampleSequenceLength.get(i);
		
		return extract(read, start, len);
	}
	
	/**
	 * 
	 * Extract the subsequence(s) corresponding to the SAMPLE blocks in the read layout
	 * 
	 * @param read the whole read or quality string
	 * @return the sample subsequences or null if this layout has no SAMPLE block
	 */
	public String[] extractSamples(String read){
		if(!hasSampleBlock)
			return null;
		String [] _sampleSeqs = new String[sampleSequenceStart.size()];
		for (int i = 0; i < sampleSequenceStart.size(); i++) {
			int start = sampleSequenceStart.get(i);
			Integer len = sampleSequenceLength.get(i);
			_sampleSeqs[i] = extract(read, start, len);
		}
		
		return _sampleSeqs;
	}
	
	
	/**
	 * Extract the subsequence(s) corresponding to the SAMPLE blocks in the read layout
	 * and merge them in a unique String (following the 5' to 3' order on the layout)
	 * 
	 * @param read the whole read or quality string
	 * @return the sequences corresponding to the SAMPLE or null if this layout has no SAMPLE block
	 */
	public String extractSample(String read){
		if(!hasSampleBlock)
			return null;
		return StringUtil.mergeArray( extractSamples(read) , "" );
	}
	
	
	private String extract(String read, int start, Integer len) {
		if(len == null){
			//take all 
			return read.substring(start);
		}
		else if(len>0){
			return read.substring(start, start+len);
		}
		else{
			return read.substring(start, read.length()+len);
		}
	}
	
}
