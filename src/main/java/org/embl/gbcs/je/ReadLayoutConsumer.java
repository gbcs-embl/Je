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

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;

import java.util.ArrayList;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ReadLayoutConsumer {
	private static Logger log = LoggerFactory.getLogger(ReadLayoutConsumer.class);
	

	private static final byte BYTECODE_SAMPLE = 0;
	private static final byte BYTECODE_BARCODE = 1;
	private static final byte BYTECODE_UMI = 2;
	
	
	ArrayList<Byte> slotCodes = new ArrayList<Byte>();
	ArrayList<Integer> slotIdx = new ArrayList<Integer>();
	ArrayList<Set<Integer>> layoutIndicesToUseForSlots = new ArrayList<Set<Integer>>();
	ReadLayout [] readLayouts;
	String outPutLayout;
	
	/**
	 * @param outPutLayout in short format
	 * @param readLayouts all ordered layout (order is as the reads are read from files)
	 */
	public ReadLayoutConsumer(String outPutLayout, ReadLayout [] readLayouts){
		this.outPutLayout = outPutLayout;
		this.readLayouts = readLayouts;

		Pattern sub = Pattern.compile("([BUS])(\\d+)");
		Matcher subMatcher = sub.matcher("");

		Pattern p = Pattern.compile("([BUS]\\d+)");
		Matcher m = p.matcher(outPutLayout);
		log.debug("ReadLayoutConsumer created with "+outPutLayout);
		while(m.find()){
			String _g = m.group();
			log.debug(_g);
			subMatcher.reset(_g);
			subMatcher.matches();
			String slotType = subMatcher.group(1);
			Integer slotId = Integer.parseInt(subMatcher.group(2));

			//save the slot code
			byte bCode ;
			if(slotType.equalsIgnoreCase("B"))
				bCode = BYTECODE_BARCODE;
			else if(slotType.equalsIgnoreCase("U"))
				bCode = BYTECODE_UMI;
			else 
				bCode = BYTECODE_SAMPLE;

			log.debug("slot type "+slotType+" => "+bCode);
			log.debug("slot id => "+slotId);

			slotCodes.add(bCode);

			//save the slot id
			slotIdx.add(slotId);

			//identify the ReadLayout to use to extract this slot
			Set<Integer> allPossibleLayoutIdxToUse = new TreeSet<Integer>();
			for (int i = 0; i < readLayouts.length; i++) {
				ReadLayout rl = readLayouts[i];

				boolean canBeUsed = false;
				switch (bCode) {
				case BYTECODE_BARCODE:
					canBeUsed = rl.containsBarcode() && rl.getBarcodeBlockUniqueIds().contains(slotId); 
					break;
				case BYTECODE_UMI:
					canBeUsed = rl.containsUMI() && rl.getUMIBlockUniqueIds().contains(slotId);
					break;
				case BYTECODE_SAMPLE:
					canBeUsed = rl.containsSampleSequence() && rl.getSampleBlockUniqueIds().contains(slotId);
					break;
				default:
					throw new LayoutMalformedException("unknown block code in output layout :"+bCode, outPutLayout);
				}

				if(canBeUsed){
					allPossibleLayoutIdxToUse.add(i);
					log.debug("ReadLayout idx "+i+" can be used for lookup");

				}
			}
			if(allPossibleLayoutIdxToUse.isEmpty())
				throw new LayoutMalformedException("no read layout identified for block "+slotType+" with index "+slotIdx, this.outPutLayout);
			layoutIndicesToUseForSlots.add(allPossibleLayoutIdxToUse);
		}
	}


	/**
	 * Assemble a read name by concatenating the output layout to the original read name.  
	 * Concatenation is made by inserting a readNameDelimitor between each added slot 
	 *  
	 * 
	 * @param reads the reads in order matching that of the {@link ReadLayout} array used at construction
	 * @param useReadSequenceForBarcodes dictates what to write in the read header layouts of the {@link FastqWriterLayout} for each BARCODE.
	 * When false, the matched barcode is used. When true, the exact read sequence extracted from the barcode slot is written
	 * @param m a {@link SampleMatch} holding all the barcode matches
	 * 
	 * @return
	 */
	public String assembleNewReadName(FastqRecord [] reads, boolean[] useReadSequenceForBarcodes, SampleMatch m){

		String newname = reads[0].getReadName().split("\\s")[0];
		if(newname.endsWith(FastqWriterLayout.readNameDelimitor))
			newname = newname.substring(0, newname.length()-1);
		
		log.debug("assembling read name with pattern "+this.outPutLayout);
		for (int i = 0; i < slotCodes.size(); i++) {
			byte slotTypeCode = slotCodes.get(i);
			int slotIdx = this.slotIdx.get(i); 
			log.debug("assembling read name adding slot code "+slotTypeCode+" with id "+slotIdx);
			/*
			 * when a slot can be obtained from different reads (e.g. redundant barcode), keep the one with best overall quality
			 */
			String subseq = null;
			int bestQual = 0;
			if(!useReadSequenceForBarcodes[slotIdx-1] && slotTypeCode == BYTECODE_BARCODE){
				// we init the subseq with the matched barcode directly
				subseq = m.getBarcodeMatches().get(slotIdx).barcode;
			}else{
				for(int rlIdx : layoutIndicesToUseForSlots.get(i)){

					ReadLayout rl = readLayouts[rlIdx];
					FastqRecord readForLayout = reads[rlIdx];

					String _subseq = null;
					String _subqual = null;
					switch (slotTypeCode) {
					case BYTECODE_BARCODE:
						_subseq  = rl.extractBarcode(readForLayout.getReadString(), slotIdx);
						_subqual = rl.extractBarcode(readForLayout.getBaseQualityString(), slotIdx);
						break;
					case BYTECODE_UMI:
						_subseq = rl.extractUMI(readForLayout.getReadString(), slotIdx);
						_subqual = rl.extractUMI(readForLayout.getBaseQualityString(), slotIdx);
						break;
					default:
						_subseq = rl.extractSample(readForLayout.getReadString(), slotIdx);
						_subqual = rl.extractSample(readForLayout.getBaseQualityString(), slotIdx);
						break;
					}
					int _qualsum = overallQuality( SAMUtils.fastqToPhred(_subqual) );
					if(subseq == null || _qualsum > bestQual){
						subseq = _subseq;
						bestQual = _qualsum;
					}
				}
			}
			//concatenate to the growing name
			newname += FastqWriterLayout.readNameDelimitor + subseq;
			log.debug("header is now : "+newname);
		}

		return newname;
	}

	/**
	 * Assemble a FastqRecord using original read name and quality header. 
	 * The read sequence (and associated quality string) are assembled according to the 
	 * output layout given to this consumer 
	 * 
	 * @param reads the reads in order matching that of the {@link ReadLayout} array used at construction
	 * @return
	 */
	public FastqRecord assembleNewRead(FastqRecord [] reads){

		String newseq = "";
		String newqual = "";

		log.debug("#####  assembling read according to layout "+this.outPutLayout);
		for (int i = 0; i < slotCodes.size(); i++) {
			byte slotTypeCode = slotCodes.get(i);
			int slotIdx = this.slotIdx.get(i); 

			log.debug("gettign info for slot code "+slotTypeCode+" with idx "+slotIdx);
			/*
			 * when a slot can be obtained from different reads (e.g. redundant barcode), keep the one with best overall quality
			 */
			String subseq = null;
			String subqual = null;
			int bestQual = 0;
			
			for(int rlIdx : layoutIndicesToUseForSlots.get(i)){
				ReadLayout rl = readLayouts[rlIdx];
				FastqRecord readForLayout = reads[rlIdx];
				String _subseq, _subqual = null;
				switch (slotTypeCode) {
				case BYTECODE_BARCODE:
					_subseq  = rl.extractBarcode(readForLayout.getReadString(), slotIdx);
					_subqual = rl.extractBarcode(readForLayout.getBaseQualityString(), slotIdx);
					break;
				case BYTECODE_UMI:
					_subseq = rl.extractUMI(readForLayout.getReadString(), slotIdx);
					_subqual = rl.extractUMI(readForLayout.getBaseQualityString(), slotIdx);
					break;
				default:
					_subseq = rl.extractSample(readForLayout.getReadString(), slotIdx);
					_subqual = rl.extractSample(readForLayout.getBaseQualityString(), slotIdx);
					break;
				}
				int _qualsum = overallQuality( SAMUtils.fastqToPhred(_subqual) );
				
				if(subseq == null || _qualsum > bestQual){
					subseq = _subseq;
					subqual = _subqual;
					bestQual = _qualsum;
				}
			}
			
			//concatenate to the growing newseq/newqual
			newseq+=subseq;
			newqual+=subqual;
			
		}
		//log.debug(newseq);
		//log.debug(newqual);
		//we borrow read name and qual header from first read, which must always be here
		return new FastqRecord(reads[0].getReadName().split("\\s")[0], newseq, reads[0].getBaseQualityHeader(), newqual);
	}

	private int overallQuality(byte [] qualScores) {
		int s = 0;
		for (byte b : qualScores) {
			s+= b;
		}
		return s;
	}

}


