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
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.FastqQualityFormat;

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
 * Note that in case of barcode, one might want to write the barcode or the read sequence corresponding to the looked up sample barcode.  
 * 
 * The possible keys are :<br/>
 * <ul>
 * <li>SAMPLEn  : refers to the SAMPLE slot with idx 'n' defined in the {@link ReadLayout} objects</li>
 * <li>UMIn     : refers to the UMI slot with idx 'n' defined in the {@link ReadLayout} objects</li>
 * <li>BARCODEn : refers to the sample barcode resolved from the read sequence found in the of the BARCODE slot with idx 'n' defined in the {@link ReadLayout} objects</li>
 * <li>READBARn : refers to the read sequence found in the BARCODE slot with idx 'n' defined in the {@link ReadLayout} objects ; this is only valida in read name layout 
 * i.e. read sequence always contain original sequence</li>
 * </ul> 
 * 
 * Note that a short layout format can also be used like 'B1', 'U2', 'S1' or 'R1' instead of '<BARCODE1>' , '<UMI2>' , '<SAMPLE1>' and <READBAR>; respectively. 
 * For example, 'B1U1U2' is the same as '<BARCODE1><UMI1><UMI2>'. 
 * 
 * Technically speaking, the short layout format is the only one used.  
 * 
 * @author girardot
 *
 */
public class FastqWriterLayout {

	private static Logger log = LoggerFactory.getLogger(FastqWriterLayout.class);
	
	public static final String DEFAULT_READNAME_DELIMITOR = ":";
	private static final String LONG_LAYOUT_REGEX = "^(<?(BARCODE|UMI|SAMPLE|READBAR)\\d+>?)+$";
	private static final String SHORT_LAYOUT_REGEX = "^([BUSR]\\d+)+$";
	
	
	
	/**
	 * char to use to delineate slots in read name ; if needed
	 */
	protected String readNameDelimitor = DEFAULT_READNAME_DELIMITOR;
	
	/**
	 * Should the quality string be injected into read name together with READBAR and UMI slots ? 
	 */
	protected boolean withQualityInReadName = false;
	

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
	 * the {@link FastqQualityFormat} of the input fastq files
	 */
	protected FastqQualityFormat fastqQualityFormat = null;
	

	
	
	/**
	 * @param readSequenceLayout the string representation of the output layout to use for the read sequence e.g. "S1"
	 * @param readNameLayout the string representation of the output layout to use for the read name e.g. "B1U1S1" 
	 * @param readLayouts the ordered {@link ReadLayout} objects defining how input fastq files are formatted
	 * @param withQualityInReadName indicates if the Barcode/UMI quality should be injected in the read name together with their sequence
	 * @param readNameDelimitor the character to use to split up the read name (':' is the default)
	 * @param convertBarcodeToReadbar  if true all BARCODE slots are converted to READBAR in the readNameLayout (BARCODE == READBAR in readSequenceLayout)
	 * @param fastqQualityFormat the {@link FastqQualityFormat} of the input fastq files
	 */
	public FastqWriterLayout(
			final String readSequenceLayout, 
			final String readNameLayout, 
			final ReadLayout [] readLayouts, 
			final boolean withQualityInReadName, 
			final String readNameDelimitor, 
			final boolean convertBarcodeToReadbar, 
			final FastqQualityFormat fastqQualityFormat) {
		
		this.readNameLayout = (StringUtils.isBlank(readNameLayout) ? null : convertToShortLayout(readNameLayout));
		this.readSequenceLayout = convertToShortLayout(readSequenceLayout);
		this.readLayouts = readLayouts;
		this.withQualityInReadName = withQualityInReadName;
		this.readNameDelimitor = readNameDelimitor;
		if(convertBarcodeToReadbar && readNameLayout!=null) {
			this.readNameLayout = this.readNameLayout.replaceAll("B", "R");
		}
		this.fastqQualityFormat = fastqQualityFormat;
		init(); //build all maps for easy lookup
	}
	

	/**
	 * @param readSequenceLayout the string representation of the output layout to use for the read sequence e.g. "S1"
	 * @param readNameLayout the string representation of the output layout to use for the read name e.g. "B1U1S1" 
	 * @param readLayouts the ordered {@link ReadLayout} objects defining how input fastq files are formatted
	 * @param withQualityInReadName indicates if the Barcode/UMI quality should be injected in the read name together with their sequence
	 * @param readNameDelimitor the character to use to split up the read name (':' is the default)
	 * @param fastqQualityFormat the {@link FastqQualityFormat} of the input fastq files
	 */
	public FastqWriterLayout(final String readSequenceLayout, final String readNameLayout, final ReadLayout [] readLayouts, final boolean withQualityInReadName, final String readNameDelimitor, final FastqQualityFormat fastqQualityFormat) {
		this(readSequenceLayout, readNameLayout, readLayouts, withQualityInReadName, readNameDelimitor, false, fastqQualityFormat);
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
			shortLayout = shortLayout.replaceAll("EADBAR", "");
		}
		log.debug("short layout : "+shortLayout);
		return shortLayout;
		
	}

	
	/**
	 * Assemble the {@link FastqRecord} that should be written in the output file according to the layout(s) ; this method should be used when a barcode has been matched
	 * 
	 * @param reads the {@link FastqRecord} from the input fastq files in the order matching the {@link ReadLayout} given at construction
	 * @param sampleMatch a {@link SampleMatch} holding all the barcode matches
	 * @return
	 */
	public FastqRecord assembleRecord( FastqRecord[] reads, SampleMatch sampleMatch ){
		
		FastqRecord rec = sequenceConsumer.assembleNewRead(reads);
		String name = rec.getReadName(); 
		if(readNameConsumer != null) 
			name = readNameConsumer.assembleNewReadName(reads, sampleMatch);
		
		FastqRecord ass = new FastqRecord(name, rec.getReadString(), rec.getBaseQualityHeader(), rec.getBaseQualityString());
		log.debug("Assembled read for output using layout [NameLayout="+this.readNameLayout+" ; SequenceLayout="+this.readSequenceLayout+"] => \n"+ass.toFastQString());
		return ass;
	}
	
	/**
	 * Convenient wrapper for single end configuration
	 * @param read the {@link FastqRecord} from the input fastq file 
	 * @param m a {@link SampleMatch} holding all the barcode matches  
	 * @return
	 */
	public FastqRecord assembleRecord( FastqRecord read, SampleMatch m ){
		return assembleRecord(new FastqRecord[]{read}, m);
	}
	
	
	
	/**
	 * 
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
			readNameConsumer = new ReadLayoutConsumer(this.readNameLayout, this.readLayouts, this.withQualityInReadName , this.readNameDelimitor, this.fastqQualityFormat);
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
	 * @return the withQualityInReadName
	 */
	public boolean isWithQualityInReadName() {
		return withQualityInReadName;
	}

	/**
	 * @param withQualityInReadName the withQualityInReadName to set
	 */
	public void setWithQualityInReadName(boolean withQualityInReadName) {
		this.withQualityInReadName = withQualityInReadName;
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
