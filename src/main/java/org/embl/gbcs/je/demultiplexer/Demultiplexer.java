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

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.FastqQualityFormat;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SolexaQualityConverter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import org.embl.cg.utilitytools.utils.CollectionUtils;
import org.embl.gbcs.je.BarcodeMatch;
import org.embl.gbcs.je.FastqWriterLayout;
import org.embl.gbcs.je.JeUtils;
import org.embl.gbcs.je.JemultiplexerFastqWriterFactory;
import org.embl.gbcs.je.Jexception;
import org.embl.gbcs.je.ReadLayout;
import org.embl.gbcs.je.SampleMatch;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Demultiplexer {

	private static Logger log = LoggerFactory.getLogger(Demultiplexer.class);

	
	//sample name for unassigned read 
	protected static final String UNASSIGNED = "__UNASSIGNED__";

	
	/*
	 * ordered list of the FASTQ in files
	 */
	File [] fastqInFiles; 
	
	/*
	 * ordered list of the ReadLayout ie matching fastqInFiles order
	 */
	ReadLayout [] readLayouts; 
	
	/*
	 * FASTQ Quality format of the input FASTQ
	 */
	FastqQualityFormat fastqQualityFormat;
	
	
	/*
	 * sample unique name to ordered list (ie matching the BARCODEn slots of the ReadLayouts) of barcode sets (only when multiple redundant barcodes were given)
	 */
	Map<String, List<Set<String>>> sample2BarcodeSets;
	
	/*
	 * for each sample gives the ordered list of output files
	 */
	Map<String, List<File>> sample2outputFileList;
	
	/*
	 * what to write in output
	 */
	FastqWriterLayout [] outputFastqLayout;
	
	/*
	 * if the FASTQ output file should be gzipped. Note that we do not rely on the file name 
	 * for this to be compatible with situations where the file name is e.g. just a number 
	 */
	boolean gzipOutput;
	
	/**
	 * if md5sum shoujdl be created for each created FASTQ output
	 */
	boolean createMD5Sum;
	
	/*
	 * output files for unassigned reads, ordered as input FASTQ (same number expected)
	 */
	List<File> unassignedFiles;
	
	/*
	 * output metrics file
	 */
	File metricFile; 
	
	/*
	 * factory to create fastq writers
	 */
	protected JemultiplexerFastqWriterFactory fastqFactory = null;
	
	/**
	 * The command line caller
	 */
	protected Jedemultiplex jedemultiplexer;
	
	
	/**
	 * internal variable storing map storing how many read pair are assigned to each sample 
	 * initialized after successful validation in doWork() 
	 */
	protected HashMap<String, Integer> sampleCountMap = null;
	
	
	/**
	 * 
	 * Create a Demultiplexer. Note that option logics should have occurred before as only minimal checks happens here.
	 * 
	 * In particular, it is expected that :
	 *    - no sample output file exist e.g. 'overwrite' option should have been dealt with before and existing files already removed
	 *    -  
	 * @param jedemultiplexer the command line caller
	 * @param fastqInFiles ordered list of the FASTQ in files
	 * @param readLayouts ordered list of the ReadLayout matching FASTQ file order in fastqInFiles
	 * @param sample2BarcodeSets ordered list (ie matching the BARCODEn slots of the ReadLayouts) of
	 *  barcode sets (only when multiple redundant barcodes were given) for each sample
	 * @param sample2outputFileList ordered list of output files for each sample 
	 * @param outputFastqLayout array of output FASTQ file layout (same number as the number of list of output files per sample)
	 * @param unassignedFiles output files for unassigned reads, ordered as input FASTQ (same number expected) ; or null to ignore unassigned reads
	 * @param metricFile output metrics file
	 * @param fastqQualityFormat FASTQ Quality format of the input FASTQ
	 * @param gzipOutput
	 * @param createMD5Sum
	 */
	public Demultiplexer(
			Jedemultiplex jedemultiplexer,
			File [] fastqInFiles, 
			ReadLayout [] readLayouts,
			Map<String, List<Set<String>>> sample2BarcodeSets, 
			Map<String, List<File>> sample2outputFileList, 
			FastqWriterLayout [] outputFastqLayout, 
			List<File> unassignedFiles,
			File metricFile, 
			FastqQualityFormat fastqQualityFormat, 
			boolean gzipOutput, 
			boolean createMD5Sum
			) {
		
		this.jedemultiplexer = jedemultiplexer;
		
		if(fastqInFiles == null || fastqInFiles.length == 0)
			throw new Jexception("no input FASTQ files provided");
		
		if(readLayouts == null || readLayouts.length == 0)
			throw new Jexception("no read layouts provided");
		
		if(readLayouts.length != fastqInFiles.length)
			throw new Jexception("The number of read layout must be the same as FASTQ input files but was "+readLayouts.length + "and "+ fastqInFiles.length+" (respectively)");
		
		if(metricFile == null)
			throw new Jexception("no metric file provided (null)");
		
		if(fastqQualityFormat == null)
			throw new Jexception("no fastq quality format provided (null)");
		
		if(sample2BarcodeSets == null || sample2BarcodeSets.size() == 0)
			throw new Jexception("no sample2BarcodeSets provided (null or empty)");
		
		if(sample2outputFileList == null || sample2outputFileList.size() == 0)
			throw new Jexception("no sample2outputFileList provided (null or empty)");
		
		if(sample2BarcodeSets.size() != sample2outputFileList.size())
			throw new Jexception("The number samples described in sample2BarcodeSets and sample2outputFileList maps is not the same : "+sample2BarcodeSets.size() + "and "+ sample2outputFileList.size()+" (respectively)");
		
		//check all entries in sample2BarcodeSets have the same number of barcode set ; which must match the overall number of barcode slots across the read layouts
		int barcodeBlockUniqueIdNumber = JeUtils.barcodeSlotCount(readLayouts);
		for(Entry<String, List<Set<String>>> e : sample2BarcodeSets.entrySet()){
			if(e.getValue().size() != barcodeBlockUniqueIdNumber)
				throw new Jexception("The number of barcode set described for sample "+e.getKey()+" ("+e.getValue().size() +") does not match the overall number of BARCODE slot across the read layouts ("+barcodeBlockUniqueIdNumber+")");
		}
		
		//check we have the same number of output file for each sample ; and that file do not already exists
		int expectedOutputNumber = sample2outputFileList.values().iterator().next().size();
		for(Entry<String, List<File>> e : sample2outputFileList.entrySet()){
			if(e.getValue().size() != expectedOutputNumber)
				throw new Jexception("The number of output files described for sample "+e.getKey()+" ("+e.getValue().size() +") does not match the number of previously observed file output ("+expectedOutputNumber+")");
			for (File _f : e.getValue()) {
				if(_f.exists()){
					throw new Jexception("Output file already exists ("+_f.getAbsolutePath()+") please make sure to delete existing files before launching me!");
				}
			}
		}
		
		/*
		 * default checks occurred, we should be pretty ok now
		 */
		
		this.fastqInFiles = fastqInFiles;			
		this.readLayouts = readLayouts;
		this.sample2BarcodeSets = sample2BarcodeSets;
		this.sample2outputFileList = sample2outputFileList;
		this.outputFastqLayout = outputFastqLayout;
		this.unassignedFiles = unassignedFiles;
		this.metricFile = metricFile;
		this.fastqQualityFormat = fastqQualityFormat;
		
		//map storing how many read pair are assigned to each sample
		sampleCountMap = new HashMap<String, Integer>(); 

	}
	
	
	





	/**
	 * @param max_mismatches one value for each defined BARCODE slot (in the read layouts)  
	 * @param min_mismatch_deltas one value for each defined BARCODE slot (in the read layouts)
	 * @param min_base_qualities one value for each defined BARCODE slot (in the read layouts)
	 * @param useReadSequenceForBarcodes dictates what to write in the read header layouts of the {@link FastqWriterLayout}.
	 * When false, the matched barcode is used. When true, the exact read sequence extracted from the barcode slot is written  
	 * @param strict how to handle barcode with redundant slots
	 * @param asyncWrite whether we should use async FASTQ writers
	 * @param diagnosticFile if not null Je writes info on sample matching process 
	 */
	public void run(int[] max_mismatches, int[] min_mismatch_deltas, int[] min_base_qualities, boolean[] useReadSequenceForBarcodes, boolean strict, boolean asyncWrite, File diagnosticFile){
		
		/*
		 * Initialize all barcode maps 
		 */
		
		// a map to get the list of all possible barcodes for a given slots
		Map<Integer, Set<String>> barcodeSetBySlotId = new LinkedHashMap<Integer, Set<String>>();
		//list of barcode length ordered by slot id
		List<Integer> orderedBarcodeLengths = new ArrayList<Integer>();
		for(Entry<String, List<Set<String>>> e : sample2BarcodeSets.entrySet()){
			for (int i = 0; i < e.getValue().size(); i++) {
				//set of redundant barcodes for this sample and for the BARCODE slot i
				Set<String> _bcs = e.getValue().get(i); 
				int bcIdx = i+1;
				if(!barcodeSetBySlotId.containsKey(bcIdx)){
					barcodeSetBySlotId.put(bcIdx, new TreeSet<String>());
					orderedBarcodeLengths.add(_bcs.iterator().next().length());
				}
				barcodeSetBySlotId.get(bcIdx).addAll(_bcs);
			}
		}
		
		
		// the same map as above but with byte[][] arrays
		Map<Integer, byte[][]> barcodeBytesBySlotId = new HashMap<Integer, byte[][]>();
		for (Entry<Integer, Set<String>> e : barcodeSetBySlotId.entrySet()) {
			int bcSlotId = e.getKey();
			List<String> tmp = new ArrayList<String>(e.getValue());
			byte[][] arr = new byte[tmp.size()][];
			for (int i = 0; i < tmp.size(); i++) {
				arr[i] = htsjdk.samtools.util.StringUtil.stringToBytes(tmp.get(i));
			}
			barcodeBytesBySlotId.put(bcSlotId, arr);
		}
		
		/*
		 * build  a barcode hashcode to sample name map for easy look up after barcode identification
		 * for this, we produce every single possible combination of barcode from all slots, hash it and save this
		 * 
		 */
		Map<Integer, String> barcodehash2sample = new HashMap<Integer, String>();
		for(Entry<String, List<Set<String>>> e : sample2BarcodeSets.entrySet()){
			String _smpl = e.getKey();
			Set<String> barcodeConcatenationCombinations = generateAllBCConcatenationCombinations(e.getValue());
			for (String combin : barcodeConcatenationCombinations) {
				int hashcode = combin.hashCode();
				if(barcodehash2sample.containsKey(hashcode) && !barcodehash2sample.get(hashcode).equals(_smpl))
					throw new Jexception("Fatal error : two barcode combinations for two distinct samples ("+barcodehash2sample.get(hashcode)+" and "+_smpl+") end up with a similar hashcode ("+hashcode+"). This tool cannot be used unless you did a mistake in the barcode-sample mapping!");
				barcodehash2sample.put(hashcode, _smpl);
			}
		}
		
		
		/*
		 * Open writers for all demultiplexed FASTQ files
		 *
		 */
		fastqFactory = new JemultiplexerFastqWriterFactory();
		fastqFactory.setUseAsyncIo(asyncWrite);
		//map to hold all sample writers
		Map<String, List<FastqWriter>> sampleFastqWriters = new HashMap<String, List<FastqWriter>>(); 

		for (Entry<String, List<File>> e : sample2outputFileList.entrySet()) {
			sampleCountMap.put(e.getKey(), 0);
			List<FastqWriter> _writers = new ArrayList<FastqWriter>();
			for (File _f : e.getValue()) {
				_writers.add(fastqFactory.newWriter(_f, gzipOutput, createMD5Sum));
			}
			sampleFastqWriters.put(e.getKey(), _writers);
		}


		/*
		 * Open writers for unassigned FASTQ reads ; if requested
		 */
		if(unassignedFiles!=null){
			if(unassignedFiles.size() != fastqInFiles.length){
				throw new Jexception("The number of files to write unassigned reads must be the same as FASTQ input files but was "+unassignedFiles.size() + "and "+ fastqInFiles.length+" (respectively)");
			}
			List<FastqWriter> unassignedReadsWriters = new ArrayList<FastqWriter>();
			for (File _f : unassignedFiles) {
				unassignedReadsWriters.add(fastqFactory.newWriter(_f, gzipOutput, createMD5Sum));
			}
			sampleFastqWriters.put(UNASSIGNED, unassignedReadsWriters);
		}
		
		
		/*
		 * Open readers on all FASTQ
		 */
		List<FastqReader> fastqReaders = new ArrayList<FastqReader>(); //we need to store them to close them at the end
		List<Iterator<FastqRecord>> fastqFileIterators = new ArrayList<Iterator<FastqRecord>>();
		for (File fqFile : fastqInFiles) {
			FastqReader r = new FastqReader(fqFile);
			fastqReaders.add(r); 
			fastqFileIterators.add( r.iterator() );
		}
		
		
		/*
		 * Iterate over all records
		 */
		
		//helper
		int barcodeBlockUniqueIdNumber = JeUtils.barcodeSlotCount(readLayouts);
		
		
		//
		PrintWriter diagnosticFileWriter = null;
		if(diagnosticFile!=null){
			try {
				diagnosticFileWriter = new PrintWriter(diagnosticFile);	
			} catch (FileNotFoundException e) {
				// should just never happen
				throw new Jexception("diagnostic file was not found!", e);
			}
			
			/*
			 * init headers :
			 * - "ReadCount" : just a read counter
			 * - "Name" : read header as found in FASTQ  
			 * - "BARCODEn_readseq" : read sub-sequence to match against barcodes from BARCODE1 ; 
			 * - "BARCODEn_MM_Best", "BARCODEn_MM_Second": Mismatch number with best and second best match
			 * - "BARCODEn_identifiedbarcode" : chosen barcode for BARCODE1 according to running options
			 * - "BARCODEn_passes_cutoffs" : tells if the barcode match fullfils the cutoff (ie bm.matched)  
			 * - "sample" : sample ultimately selected 
			 */
			
			diagnosticFileWriter.print("ReadCount");
			diagnosticFileWriter.print("Name");
			for (int i = 1; i <= barcodeBlockUniqueIdNumber; i++) {
				diagnosticFileWriter.print("\t"+"BARCODE"+i+"_readseq");
				diagnosticFileWriter.print("\t"+"BARCODE"+i+"_bestbarcode");
				diagnosticFileWriter.print("\t"+"BARCODE"+i+"_MM_Best");
				diagnosticFileWriter.print("\t"+"BARCODE"+i+"_MM_Second");
				diagnosticFileWriter.print("\t"+"BARCODE"+i+"_passes_cutoffs");
			}
			diagnosticFileWriter.println("\t"+"assigned_sample"+"\t"+"notes");
		}
		
		
		//counters
		int cnt = 0; //read counter
		int unassigned = 0; //count how many pairs end up unassigned 
		int assigned = 0; //count how many pairs end up with a sample assignment
		Iterator<FastqRecord> mainIterator = fastqFileIterators.get(0);
		while(mainIterator.hasNext()){
			//reads next read from all input files
			FastqRecord[] reads = nextReads(fastqFileIterators);
			cnt ++;
			//extract the barcoding sequences from each BARCODE slot ; when List<FastqRecord> is more than 1 => means we have the SAME barcode from different slots (redundant)
			Map<Integer, List<FastqRecord>> barcodeSubsequenceBySlotIdx = FastqWriterLayout.extractBarcodeSlots(reads, readLayouts);
			
			//identify the sample matching these subsequence
			SampleMatch assignedSample = assignToSample(
					barcodeSubsequenceBySlotIdx, 
					barcodeSetBySlotId, orderedBarcodeLengths, 
					barcodeBytesBySlotId, barcodehash2sample, 
					min_base_qualities, max_mismatches, min_mismatch_deltas, strict);
			
			writeDiagnostics(reads, assignedSample, diagnosticFileWriter, cnt);
			
			//write to output(s)
			if(!assignedSample.getSample().equals(Demultiplexer.UNASSIGNED) ){
				log.debug("read set assigned to "+assignedSample.getSample());
				int c = sampleCountMap.get(assignedSample.getSample());
				c++;
				sampleCountMap.put(assignedSample.getSample(), c);
				
				assigned++;
				List<FastqWriter> writers  = sampleFastqWriters.get(assignedSample.getSample());
				for (int i = 0; i < writers.size(); i++) {
					//prepare the output according to output layout
					log.debug("Writing in output idx "+(i+1));
					FastqRecord rec = outputFastqLayout[i].assembleRecord( reads, useReadSequenceForBarcodes, assignedSample );
					writers.get(i).write(rec);
				}
			}
			else if(sampleFastqWriters.containsKey(UNASSIGNED)){
				log.debug("No assigned sample for this read set");
				unassigned++;
				//write unmodified reads when not assigned
				List<FastqWriter> writers  = sampleFastqWriters.get(UNASSIGNED);
				for (int i = 0; i < writers.size(); i++) {
					writers.get(i).write(reads[i]);
				}
			} 
			
		}
		
		
		//close readers
		for(FastqReader r : fastqReaders){
			try {
				r.close();
			} catch (Exception e) {
				// ignore
			}
		}
				
		//close all writers
		for(List<FastqWriter> _l : sampleFastqWriters.values()){
			for (FastqWriter w : _l) {
				try {
					w.close();
				} catch (Exception e) {
					// ignore
				}
			}
		}
		
		//close diag
		if(diagnosticFileWriter!=null){
			try {
				diagnosticFileWriter.close();
			} catch (Exception e) {
				// ignore
			}
		}
		
		//print counts
		jedemultiplexer.printMetricFile(sampleCountMap, cnt, unassigned, assigned);

	}

	

	/**
	 * @param reads the original reads
	 * @param assignedSample the sample match report
	 * @param diagnosticFileWriter the writer or null if no diagnostics has to be written 
	 * @param readCounter the current read iteration (starts at one) 
	 */
	private void writeDiagnostics(FastqRecord[] reads,
			SampleMatch assignedSample, 
			PrintWriter diagnosticFileWriter, 
			int readCounter) {
		
		//should we write something ?
		if(diagnosticFileWriter == null)
			return;
		
		int barcodeBlockUniqueIdNumber = assignedSample.getBarcodeMatches().size();
		
		diagnosticFileWriter.print(readCounter);
		diagnosticFileWriter.print("\t"+reads[0].getReadName().split("\\s")[0]);
		for (int i = 1; i <= barcodeBlockUniqueIdNumber; i++) {
			BarcodeMatch bm = assignedSample.getBarcodeMatches().get(i);
			diagnosticFileWriter.print("\t"+bm.readSequence);
			diagnosticFileWriter.print("\t"+bm.barcode);
			diagnosticFileWriter.print("\t"+bm.mismatches);
			diagnosticFileWriter.print("\t"+bm.mismatchesToSecondBest);
			diagnosticFileWriter.print("\t"+ (bm.matched? "yes": "no") );
		}
		diagnosticFileWriter.print("\t"+(assignedSample.getSample().equals(Demultiplexer.UNASSIGNED) ? "unassigned" : assignedSample.getSample()));
		diagnosticFileWriter.println("\t"+assignedSample.getDiagnosticNote());
		
		if(readCounter % 100 == 0)
			diagnosticFileWriter.flush();
	}







		
	/**
	 * @param barcodeSubsequenceBySlotIdx associates each BARCODE slot (keyed by its ID) with the list of (redundant) 
	 * barcodes sequences for this BARCODE slot. This is a list as a given BARCODE slot can appear more than once across
	 *  the read layouts i.e. in the case of redundant barcode
	 * @param barcodeSetBySlotIdx set of possible barcodes for a given BARCODE slot id 
	 * @param orderedBarcodeLengths ordered list of barcode length (following the concatenation ordered used for producing the hashcodes)
	 * @param barcodeBytesBySlotIdx set of possible barcodes for a given BARCODE slot id in byte format
	 * @param barcodehash2sample every single possible combination of barcode from all slots (one per slot in each combination) hash
	 * @param min_base_qualities 
	 * @param max_mismatches
	 * @param min_mismatch_deltas
	 * @param strict how to handle barcode with redundant sequence slots
	 * @return a {@link SampleMatch} in which the sample name is set to Demultiplexer.UNASSIGNED if barcode lookup failed
	 */
	private SampleMatch assignToSample(
			Map<Integer, List<FastqRecord>> barcodeSubsequenceBySlotIdx,
			Map<Integer, Set<String>> barcodeSetBySlotIdx,
			List<Integer> orderedBarcodeLengths,
			Map<Integer, byte[][]> barcodeBytesBySlotIdx,
			Map<Integer, String> barcodehash2sample,
			int [] min_base_qualities, 
			int [] max_mismatches, 
			int [] min_mismatch_deltas,
			boolean strict
			) {
		
		/*
		 * We first blindly match all barcodeSubsequence to the list of expected barcodes for this slot idx 
		 * => we convert the barcodeSubsequenceBySlotIdx map to an equivalent map of BarcodeMatch
		 */
		Map<Integer, List<BarcodeMatch>> barcodeMatchBySlotIdx = new HashMap<Integer, List<BarcodeMatch>>(); // here it is a list to cope with BARCODE with multiple redundant location across reads  
		for(Entry<Integer, List<FastqRecord>> e : barcodeSubsequenceBySlotIdx.entrySet()){
			int slotIdx = e.getKey();
			
			barcodeMatchBySlotIdx.put(slotIdx, new ArrayList<BarcodeMatch>());
			// expected barcodes for this slot
			byte[][] barcodeBytes = barcodeBytesBySlotIdx.get(slotIdx);
			Set<String> expectedBarcodes =  barcodeSetBySlotIdx.get(slotIdx);
			//now match
			for(FastqRecord _rec : e.getValue()){
				log.debug("assignToSample : looking at barcodes for slotIdx "+ slotIdx+" => current barcode sequence is \n"+_rec.toFastQString());
				BarcodeMatch bm = null;
				//optimization : we first check if the sequence is one of the expected one
				if(expectedBarcodes.contains(_rec.getReadString())){
					bm = new BarcodeMatch();
					bm.matched = true;
					bm.barcode = _rec.getReadString();
					bm.readSequence = _rec.getReadString();
					bm.mismatches = 0;
					bm.mismatchesToSecondBest = _rec.getReadLength();
				}else{
					int p = slotIdx ;
					bm = findBestBarcode(_rec, barcodeBytes, min_base_qualities[p-1], max_mismatches[p-1], min_mismatch_deltas[p-1]);
				}
				log.debug("     best barcode match is "+bm.toString());
				barcodeMatchBySlotIdx.get(slotIdx).add(bm);
			}
		}
		
		//identify the corresponding sample, if any 
		//set of all concatenated codes and the sum of their mismatch
		Map<String, Integer> concatenatedCodes = new HashMap<String, Integer>();
		concatenatedCodes.put("", 0); //init with empty string
		Map<Integer, Map<String, BarcodeMatch>> barcodeMatches = new HashMap<Integer, Map<String, BarcodeMatch>>();
		boolean hasNoSample = false;
		for(Entry<Integer, List<BarcodeMatch>> e : barcodeMatchBySlotIdx.entrySet()){
			/*
			 * for those slots with more than one BARCODE MATCH, we need to consider all possible concatenations
			 */
			int slotIdx = e.getKey();
			Map<String, BarcodeMatch> allValidAndNotRedundantBarcodeMatches = keepOnlyBestBarcodeMatches(e.getValue());  // keyed by the barcode sequence
			if(allValidAndNotRedundantBarcodeMatches.size() == 0){
				hasNoSample = true; //we can t look up a sample
				break;
			}
			//remember for later
			barcodeMatches.put(slotIdx, allValidAndNotRedundantBarcodeMatches);
			//augment the concatenated codes
			Map<String, Integer>  augmentedCodes = new HashMap<String, Integer>();
			for (Entry<String, Integer> _concat : concatenatedCodes.entrySet()) {
				for (Entry<String, BarcodeMatch> toAdd : allValidAndNotRedundantBarcodeMatches.entrySet()) {
					augmentedCodes.put(
							_concat.getKey() + toAdd.getKey() , 
							_concat.getValue() + toAdd.getValue().mismatches
							); 
				}
			}
			concatenatedCodes = augmentedCodes;
			
		}
		
		//do we have concatenated string(s) representing whole the bc slots?
		String sampleName = "";
		//a note to add to the diagnostic
		String diagNote = "";
		if( false == hasNoSample){
			//do these string resolved to the same sample ?
			Map<String, Integer> sampleNames = new HashMap<String, Integer>();
			Map<String, String> sampleName2concatenatedCode = new HashMap<String, String>();
			for (Entry<String, Integer> code : concatenatedCodes.entrySet()) {
				String sname = barcodehash2sample.get(code.getKey().hashCode());
				sampleNames.put( sname , code.getValue());
				sampleName2concatenatedCode.put(sname, code.getKey());
			}
			
			//if there is a unique sample assignment
			if(sampleNames.size() == 1){
				Entry<String, Integer> en = sampleNames.entrySet().iterator().next();
				sampleName = en.getKey();
				// pick a unique match per slot
				Map<Integer, BarcodeMatch> uniqueBCMatches = new HashMap<Integer, BarcodeMatch>();
				for (Entry<Integer, Map<String, BarcodeMatch>> e : barcodeMatches.entrySet()) {
					uniqueBCMatches.put(e.getKey(), e.getValue().values().iterator().next());
				}
				return new SampleMatch(sampleName, uniqueBCMatches);
			} 
			// OR if NOT strict
			else if (!strict) {
				// for each possible sample, compute the overall sum of mismatches
				Integer lowestMM = null;
				Map<Integer, Set<String>> mm2samples = new HashMap<Integer, Set<String>>();
				for (Entry<String, Integer> e : sampleNames.entrySet()){
					String _sample = e.getKey();
					int mm = e.getValue();
					if(!mm2samples.containsKey(mm)){
						mm2samples.put(mm, new TreeSet<String>());
					}
					mm2samples.get(mm).add(_sample);
					if(lowestMM == null || mm< lowestMM) lowestMM = mm;
				}
				
				List<Integer> orderedMMs = new ArrayList<Integer>(mm2samples.keySet());
				Collections.sort(orderedMMs);
				
				for (Integer _mm : orderedMMs) {
					for (String _sampl : mm2samples.get(_mm)) {
						diagNote += (diagNote.isEmpty() ? "" : " ; ");
						diagNote += _sampl +"("+_mm + " MMs)"; 
					}
				}
				
				// is there a better assignment ie a single sample with lowest overall MM number?
				if(mm2samples.get(lowestMM).size() == 1){
					//we have a better sample
					sampleName = mm2samples.get(lowestMM).iterator().next();
					//extract the barcodes of the concatenated barcode and indentify back the BarcodeMatch
					Map<Integer, BarcodeMatch> uniqueBCMatches = new HashMap<Integer, BarcodeMatch>();
					int from = 0;
					int _slotIdx = 0;
					String concatenatedBC = sampleName2concatenatedCode.get(sampleName);
					for(int bcLen : orderedBarcodeLengths){
						_slotIdx++;
						int end = from + bcLen;
						String _bc = concatenatedBC.substring(from, end);
						BarcodeMatch bcM = barcodeMatches.get(_slotIdx).get(_bc);
						from = end;
						uniqueBCMatches.put(_slotIdx, bcM);
					}
					log.debug( "    selecting "+sampleName+" as it has the lowest overall MM count : "+lowestMM);
					diagNote = "Selected "+sampleName+" due to lowest overall MM from : " + diagNote;
					log.debug( "    "+diagNote);
					return new SampleMatch(sampleName, uniqueBCMatches, diagNote);
				}else{
					diagNote = "Cannot select from : " + diagNote;
				}
			}else{
				for (Entry<String, Integer> e : sampleNames.entrySet()) {
					diagNote += (diagNote.isEmpty() ? "" : " ; ");
					diagNote += e.getKey() +"("+e.getValue() + " MMs)"; 
				}
				diagNote = "Cannot select from : " + diagNote;
			}
			log.debug("    barcodes'matches resolve to multiple samples :" + Arrays.toString(sampleNames.keySet().toArray()));
		}
		
		//build a fake match set for diag file
		Map<Integer, BarcodeMatch> uniqueBCMatches = new HashMap<Integer, BarcodeMatch>();
		for(Entry<Integer, List<BarcodeMatch>> e : barcodeMatchBySlotIdx.entrySet()){
			uniqueBCMatches.put(e.getKey(), buildFakeBarcodeMatchForNoSampleLookupSituation(e.getValue()));
		}
		if(diagNote.isEmpty())
			diagNote = "indicated mismatches : -1 for no match else lowest mismatch";
		return new SampleMatch(Demultiplexer.UNASSIGNED, uniqueBCMatches, diagNote);
	}

	
	
	private BarcodeMatch buildFakeBarcodeMatchForNoSampleLookupSituation(
			Collection<BarcodeMatch> values) {
		if(values.size() == 1)
			return values.iterator().next();
		
		BarcodeMatch fake = new BarcodeMatch();
		fake.readSequence = "";
		fake.barcode = "";
		fake.matched = false;
		fake.mismatches = -1;
		fake.mismatchesToSecondBest = -1;
		for (BarcodeMatch bm : values) {
			fake.readSequence += (fake.readSequence.isEmpty() ? "":"," ) + bm.readSequence;
			
			if(bm.matched){
				fake.matched = true;
				fake.barcode += (fake.barcode.isEmpty() ? "":"," ) + (bm.barcode.isEmpty() ? "NOMATCH":bm.barcode ) ;
				if(fake.mismatches < 0 || fake.mismatches > bm.mismatches ) {
					fake.mismatches = bm.mismatches;
					fake.mismatchesToSecondBest = bm.mismatchesToSecondBest;
				}
			}
		}
		return fake;
	}








	/**
	 * removes the unmatched BarcodeMatch and only keep the best match when more than one BarcodeMatch
	 * are found for the same barcode
	 * @param matches
	 * @return the best 
	 */
	private Map<String, BarcodeMatch> keepOnlyBestBarcodeMatches(List<BarcodeMatch> matches) {
		
		Map<String, BarcodeMatch> allValidAndNotRedundantBarcodeMatches = new HashMap<String, BarcodeMatch>();
		for (BarcodeMatch _bm : matches) {
			if(!_bm.matched)
				continue;

			BarcodeMatch best = allValidAndNotRedundantBarcodeMatches.get(_bm.barcode);
			if(best == null || _bm.mismatches < best.mismatches){
				best = _bm;
				allValidAndNotRedundantBarcodeMatches.put(_bm.barcode, best);
			}
		}

		return allValidAndNotRedundantBarcodeMatches;
	}


	/**
	 * Find the best barcode match for the given read sequence ; only called when no exact match was found
	 * @param readSlot the read subsequence of a BARCODE slot  
	 * @param barcodeBytes expected barcode sequences for this slot
	 * 
	 * @return a BarcodeMatch holding the matched status and extra information
	 */
	protected BarcodeMatch findBestBarcode(final FastqRecord readSlot, final byte[][] barcodeBytes, final int minQuality, final int maxMismatches, final int minMMDelta) {
		log.debug(readSlot.toFastQString());
		//turn to bytes
		byte [] subseqBytes = readSlot.getReadBases();
		byte [] qualBytes = readSlot.getBaseQualityString().getBytes();
		log.debug("   Q bytes => "+Arrays.toString(qualBytes));
		convertQuality(qualBytes, this.fastqQualityFormat);
		log.debug("   converted Q bytes => "+Arrays.toString(qualBytes));
		
		int numMismatchesInBestBarcode = readSlot.getReadLength() + 1; //init with max mismatch num + 1
		int numMismatchesInSecondBestBarcode= readSlot.getReadLength() + 1; //init with max mismatch num + 1
		
		//find best and second best barcode
		byte[] bestBarcode = null;
		
		
		for (int i = 0; i < barcodeBytes.length; i++) {
			final int numMismatches = countMismatches(barcodeBytes[i], subseqBytes, qualBytes, minQuality);
			if (numMismatches < numMismatchesInBestBarcode) {
				if (bestBarcode != null) {
					numMismatchesInSecondBestBarcode = numMismatchesInBestBarcode;
				}
				numMismatchesInBestBarcode = numMismatches;
				bestBarcode = barcodeBytes[i];
			} else if (numMismatches < numMismatchesInSecondBestBarcode) {
				numMismatchesInSecondBestBarcode = numMismatches;
			}
		}


		final boolean matched = (bestBarcode != null &&
				numMismatchesInBestBarcode <= maxMismatches &&
				numMismatchesInSecondBestBarcode - numMismatchesInBestBarcode >= minMMDelta)
				;


		final BarcodeMatch match = new BarcodeMatch();
		match.matched = matched;
		match.readSequence = readSlot.getReadString();
		match.mismatches = numMismatchesInBestBarcode;
		match.mismatchesToSecondBest = numMismatchesInSecondBestBarcode;
		match.barcode = htsjdk.samtools.util.StringUtil.bytesToString( bestBarcode );

		return match;
	}
	
	/**
	 * 
	 * Compare barcode sequence to bases from read
	 * @param barcodeBytes the barcode as array of bytes
	 * @param readSubsequence the subsequence as array of bytes
	 * @param qualities the base quality score in phred scale (+33)
	 * @param minimumBaseQuality an int in [0,40]
	 * @return how many bases did not match , 'N' or base with quality less than MINIMUM_BASE_QUALITY always considered
	 * as a mismatch 
	 * 
	 */
	protected  int countMismatches(final byte[] barcodeBytes, final byte[] readSubsequence, final byte[] qualities, final int minimumBaseQuality) {
		int numMismatches = 0;
		// Read sequence and barcode length may not be equal, so we just use the shorter of the two ; I don t see why this should happen but it is free to do this !
		final int basesToCheck = Math.min(barcodeBytes.length, readSubsequence.length);
		for (int i = 0; i < basesToCheck; ++i) {
			/*
			 * we first need to check this is a valid position in barcode (barcode can contain N to specify this position should be ignored)
			 * If the bc contains a N, we skip
			 */
			if (!SequenceUtil.isValidBase(barcodeBytes[i])) { // Returns true if the byte is in [acgtACGT].
				continue;
			}
			
			if (!SequenceUtil.isNoCall(readSubsequence[i])) { //returns true if the value of base represents a no call
				//if bases in barcode and readsubsequence are different => increase mismatch count ; else check quality
				if (!SequenceUtil.basesEqual(barcodeBytes[i], readSubsequence[i])){
					++numMismatches;
				}
				else if (qualities != null){
					final int uQual = qualities[i] & 0xff; //from FastqToSam.java ; this is kind of bit masking 
					if(uQual < minimumBaseQuality)
						++numMismatches;
				}
			}else{
				//the read subsequence contains a 'N' at this position => this is a mismatch 
				++numMismatches;
			}
		}

		return numMismatches;
	}

	
	/** 
	 * Based on the type of quality scores coming in, converts them to a numeric byte[] in phred scale. 
	 */
	protected void convertQuality(byte[] quals, final FastqQualityFormat version) {
		switch (version)  {
		case Standard:
			SAMUtils.fastqToPhred(quals);
			break ;
		case Solexa:
			SolexaQualityConverter.getSingleton().convertSolexaQualityCharsToPhredBinary(quals);
			break ;
		case Illumina:
			SolexaQualityConverter.getSingleton().convertSolexa_1_3_QualityCharsToPhredBinary(quals);
			break ;
		}
	}
	
	
	private FastqRecord[] nextReads(
			List<Iterator<FastqRecord>> fastqFileIterators) {
		FastqRecord[] reads = new FastqRecord[fastqFileIterators.size()];
		for (int j = 0; j < fastqFileIterators.size(); j++) {
			reads[j] = fastqFileIterators.get(j).next();
		}
		return reads;
	}


	/**
	 * computes all barcode combinations. Each combination has a single barcode from each set
	 * 
	 * @param barcodeSetList the list of barcode sets for a given sample. At the easiest, each set has a single entry
	 * @return the set unique barcode combinations.
	 */
	private Set<String> generateAllBCConcatenationCombinations(List<Set<String>> barcodeSetList) {
		
		//we ALWAYS have one BARCODE slot (otherwise demultiplexing is not possible!!)
		Set<String> concatenatedBarcodes = barcodeSetList.get(0); 
		//how many slots do we have in total ?
		int n = barcodeSetList.size();
		//recursively create them all 
		for(int i = 1 ; i <= n-1 ; i++){
			concatenatedBarcodes = allPairCombinations(concatenatedBarcodes, barcodeSetList.get(i)) ;	
		}
		 
		return concatenatedBarcodes;
	}
	
	private Set<String> allPairCombinations(Set<String> set1, Set<String> set2) {
		Set<String> c = new TreeSet<String>();
		for (String a : set1) {
			for (String b : set2) {
				c.add(a+b);
			}
		}
	
		return c;
	}

}
