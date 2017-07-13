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

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import picard.util.MathUtil;
import htsjdk.samtools.fastq.FastqRecord;

/**
 * 
 * Describes how and what to write in FASTQ output files when demultiplexing input FASTQ files
 * More precisely, it describes the output file layout(s) using the slots defined in read layouts e.g. '<BARCODE1><UMI1><UMI2>' or '<SAMPLE1>'  
 * Each FastqWriterLayout needs two of these descriptors : 
 *   1. One describing how to write the read sequence e.g. '<SAMPLE1>' (only writes the sample sequence) or '<BARCODE1><SAMPLE1>' to also 
 *      keep the barcode in the output read sequence (and qualities)
 *   2. A second one describing how to write the read name (header) e.g. '<BARCODE1><UMI1><UMI2>' to add the barcode and two extracted UMIs 
 *      in the final read name, in addition to the original read name (ie header up to the space). Here each written slot is separated with ':' by default
 *      
 * 
 * Note that a short layout format can also be used like 'B1', 'U2', 'S1'' instead of '<BARCODE1>' , '<UMI2>' and '<SAMPLE1>' ; respectively. 
 * For example, 'B1U1U2' is the same as '<BARCODE1><UMI1><UMI2>'. 
 * 
 *  Technically speaking, the short layout format is the only one used.  
 * 
 * @author girardot
 *
 */
public class FastqWriterLayout {

	private static Logger log = LoggerFactory.getLogger(FastqWriterLayout.class);
	
	private static final String LONG_LAYOUT_REGEX = "^(<?(BARCODE|UMI|SAMPLE)\\d+>?)+$";
	private static final String SHORT_LAYOUT_REGEX = "^([BUS]\\d+)+$";
	
	
	
	/**
	 * char to use to delineate slots in read name ; if needed
	 */
	protected static String readNameDelimitor = ":";
	
	
	
	/**
	 * Layout for writing the read sequence ; in short format 
	 */
	protected String readSequenceLayout;
	
	/**
	 * Consumer associated with the readSequenceLayout
	 */
	protected ReadLayoutConsumer sequenceConsumer;
	
	/**
	 * Layout for writing the read name ; in short format 
	 */
	protected String readNameLayout;
	
	/**
	 * Consumer associated with the readNameLayout
	 */
	protected ReadLayoutConsumer readNameConsumer ;
	
	
	
	/**
	 * All the {@link ReadLayout} defined in the demultiplexing ; used to find how to extract needed information  
	 */
	protected ReadLayout [] readLayouts;
	
	/**
	 * @param readSequenceLayout
	 * @param readNameLayout can be null when the read name should be reused unmodified
	 * @param readLayouts
	 */
	public FastqWriterLayout(final String readSequenceLayout, final String readNameLayout, final ReadLayout [] readLayouts) {
		
		this.readNameLayout = (StringUtils.isBlank(readNameLayout) ? null : convertToShortLayout(readNameLayout));
		this.readSequenceLayout = convertToShortLayout(readSequenceLayout);
		this.readLayouts = readLayouts;
		init(); //build all maps for easy lookup
	}
	
	/**
	 * @param readSequenceLayout
	 * @param readNameLayout can be null when the read name should be reused unmodified
	 * @param readLayout
	 */
	public FastqWriterLayout(final String readSequenceLayout, final String readNameLayout, final ReadLayout readLayout) {
		this(readSequenceLayout, readNameLayout, new ReadLayout[]{ readLayout });
	}

	
	/**
	 * @param layout
	 * @return
	 */
	private String convertToShortLayout(final String layout) {
		log.debug("given layout : "+layout);
		if(StringUtils.isBlank(layout))
			return layout;
		String shortLayout = layout;
		if(!Pattern.matches(SHORT_LAYOUT_REGEX, layout)){
			if(!Pattern.matches(LONG_LAYOUT_REGEX, layout))
				throw new LayoutMalformedException("FASTQ Output Layout does not match expected short ("+SHORT_LAYOUT_REGEX+") nor long ("+LONG_LAYOUT_REGEX+") formats", layout);
			
			//convert to short
			shortLayout = shortLayout.replaceAll("<", "");
			shortLayout = shortLayout.replaceAll(">", "");
			shortLayout = shortLayout.replaceAll("ARCODE", "");
			shortLayout = shortLayout.replaceAll("MI", "");
			shortLayout = shortLayout.replaceAll("AMPLE", "");
		}
		log.debug("short layout : "+shortLayout);
		return shortLayout;
		
	}


	/**
	 * Assemble the {@link FastqRecord} that should be written in the output file according to the layout(s)
	 * @param reads the {@link FastqRecord} from the input fastq files in the order matching the {@link ReadLayout} given at construction 
	 * @param useReadSequenceForBarcodes dictates what to write in the read header layouts of the {@link FastqWriterLayout} for each BARCODE.
	 * When false, the matched barcode is used. When true, the exact read sequence extracted from the barcode slot is written
	 * @param m a {@link SampleMatch} holding all the barcode matches
	 * @return
	 */
	public FastqRecord assembleRecord( FastqRecord[] reads, boolean[] useReadSequenceForBarcodes, SampleMatch m ){
		
		FastqRecord rec = sequenceConsumer.assembleNewRead(reads);
		String name = rec.getReadName(); 
		if(readNameConsumer != null) 
			name = readNameConsumer.assembleNewReadName(reads, useReadSequenceForBarcodes, m);
		
		FastqRecord ass = new FastqRecord(name, rec.getReadString(), rec.getBaseQualityHeader(), rec.getBaseQualityString());
		log.debug("Assembled read for output using layout [NameLayout="+this.readNameLayout+" ; SequenceLayout="+this.readSequenceLayout+"] => \n"+ass.toFastQString());
		return ass;
	}
	
	/**
	 * Convenient wrapper for single end configuration
	 * @param read the {@link FastqRecord} from the input fastq file
	 * @param useReadSequenceForBarcodes dictates what to write in the read header layouts of the {@link FastqWriterLayout} for each BARCODE.
	 * When false, the matched barcode is used. When true, the exact read sequence extracted from the barcode slot is written
	 * @param m a {@link SampleMatch} holdign all the barcode matches  
	 * @return
	 */
	public FastqRecord assembleRecord( FastqRecord read, boolean[] useReadSequenceForBarcodes, SampleMatch m ){
		return assembleRecord(new FastqRecord[]{read}, useReadSequenceForBarcodes, m);
	}
	
	
	
	/**
	 * Checks that this layout can be produced from a list of {@link ReadLayout}
	 * @param readLayouts
	 * @return
	 */
	protected void init(){
		
		
		/*
		 * Process (short format) layout for easy output assembly at FASTQ writing time
		 */
		
		// do for read seq
		if(!Pattern.matches(SHORT_LAYOUT_REGEX, this.readSequenceLayout)){
			throw new LayoutMalformedException("FASTQ Output Layout for read sequence does not match expected short format (regex is :"+SHORT_LAYOUT_REGEX+")", this.readSequenceLayout);
		}
		sequenceConsumer = new ReadLayoutConsumer(this.readSequenceLayout, this.readLayouts);

		// do for read name if not null
		if( this.readNameLayout != null){
			if(!Pattern.matches(SHORT_LAYOUT_REGEX, this.readNameLayout)){
				throw new LayoutMalformedException("FASTQ Output Layout for read name does not match expected short format (regex is :"+SHORT_LAYOUT_REGEX+")", this.readNameLayout);
			}
			readNameConsumer = new ReadLayoutConsumer(this.readNameLayout, this.readLayouts);
		}
	}
	
	
	/**
	 * @param reads
	 * @param readLayouts
	 * @return
	 */
	public static Map<Integer, List<FastqRecord>> extractBarcodeSlots(FastqRecord[] reads, ReadLayout[] readLayouts) {
		Map<Integer, List<FastqRecord>> m = new HashMap<Integer, List<FastqRecord>>();
		
		/*
		 * for each read layout
		 */
		for (int i = 0; i < readLayouts.length; i++) {
			//find those with BARCODE slot(s)
			if(readLayouts[i].containsBarcode()){
				//extract the BARCODE subsequence
				String [] sequences = readLayouts[i].extractBarcodes(reads[i].getReadString());
				//and the corresponding quality strings
				String [] qualities = readLayouts[i].extractBarcodes(reads[i].getBaseQualityString());
				//and the BARCODE idx corresponding to them
				List<Integer> barcodeBlockIds = readLayouts[i].getOrderedBarcodeBlockUniqueIds();
				//save each subsequence/quality as a FastqRecord in the list corresponding to the BARCODE idx
				for (int j = 0; j < barcodeBlockIds.size(); j++) {
					int blockId = barcodeBlockIds.get(j);
					if(!m.containsKey(blockId))
						m.put(blockId, new ArrayList<FastqRecord>());
					
					FastqRecord fr = new FastqRecord(null, sequences[j], null, qualities[j]);
					log.debug("Extracted Barcode : "+fr.toFastQString());
					m.get(blockId).add(
							fr
							);
				}
			}
		}
		return m;
	}


	/**
	 * @return the readNameDelimitor
	 */
	public String getReadNameDelimitor() {
		return readNameDelimitor;
	}


	/**
	 * @param readNameDelimitor the readNameDelimitor to set
	 */
	public void setReadNameDelimitor(String readNameDelimitor) {
		this.readNameDelimitor = readNameDelimitor;
	}

	/**
	 * @return the readSequenceLayout in short format
	 */
	public String getReadSequenceLayout() {
		return readSequenceLayout;
	}

	/**
	 * @return the readNameLayout in short format
	 */
	public String getReadNameLayout() {
		return readNameLayout;
	}


	
}
