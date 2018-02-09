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
package org.embl.gbcs.je.jeclipper;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.FastqQualityFormat;
import htsjdk.samtools.util.QualityEncodingDetector;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

import org.embl.cg.utilitytools.utils.ExceptionUtil;
import org.embl.cg.utilitytools.utils.FileUtil;
import org.embl.cg.utilitytools.utils.StringUtil;
import org.embl.cg.utilitytools.utils.parser.csv.CSVLine;
import org.embl.cg.utilitytools.utils.parser.csv.InvalidHeaderException;
import org.embl.gbcs.embase.api.EmBaseDatabaseFacade;
import org.embl.gbcs.embase.api.Queries;
import org.embl.gbcs.embase.api.exception.EmBASEConnectionException;
import org.embl.gbcs.embase.api.model.Item;
import org.embl.gbcs.embase.api.model.NGSLibrary;
import org.embl.gbcs.je.ApplicationConfiguration;
import org.embl.gbcs.je.FastqWriterLayout;
import org.embl.gbcs.je.Je;
import org.embl.gbcs.je.JeUtils;
import org.embl.gbcs.je.JemultiplexerFastqWriterFactory;
import org.embl.gbcs.je.Jexception;
import org.embl.gbcs.je.ReadLayout;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

@CommandLineProgramProperties(
		usage = "Reads records in the supplied FASTQ file(s) according to specified read layouts (RL option) and write output FASTQ file(s)"
				+ " according to supplied output layouts (OL option).\n" , 
		usageShort = "je clip F=fastq_1.txt.gz F=fastq_2.txt.gz RL=<BARCODE1:6><UMI1:8><SAMPLE1:x> RL=<BARCODE1:6><UMI2:8><SAMPLE2:x> OL=1:B1U1U2:S1 OL=2:B1U1U2:S2"
		)
public class Jeclipper extends CommandLineProgram {
	
	
	private static Logger log = LoggerFactory.getLogger(Jeclipper.class);

	/*
	 * Option defaults
	 */
	protected static final Boolean DEFAULT_GZIP_OUTPUTS = true;
	protected static final Boolean DEFAULT_WRITER_FACTORY_USE_ASYNC_IO = true;
	protected static final String DEFAULT_READ_NAME_SEPARATOR_CHAR = ":";
	
	protected static final boolean DEFAULT_ADD_SEQUENCE_LAYOUT_IN_OUTPUT_FILENAME = false;

	protected static final boolean DEFAULT_ADD_HEADER_LAYOUT_IN_OUTPUT_FILENAME = false;

	protected static final boolean DEFAULT_ADD_LAYOUT_IDX_IN_OUTPUT_FILENAME = true;

	
	
	/*
	 * Common line options
	 */
	
	// input fastq files
	@Option(shortName="F", 
			optional = false,
			printOrder=10,
			doc="Input fastq file (optionally gzipped)",
			minElements = 1	
			)
	public List<File> FASTQ;
	
	//Read layouts 
	@Option(shortName="RL", optional = false,
			printOrder=30,
			doc="Describes the read layout(s) e.g. 'RL=<BARCODE1:6><SAMPLE:x>' of input fastq file(s). "+
			"The input fastq files and read layouts are mached up by order on the command line.\n" +
			"Read layouts are only needed for complex layouts but one must provide read layouts for ALL or NONE of the input fastq files.\n"+
			"Read layouts are made of <UMIn:X>, <BARCODEn:X>, <SAMPLEn:X> blocks to describe blocks of type UMI, BARCODE or SAMPLE with : \n "+
			"   * 'n' the unique block index (an index must be unique across all read layouts for each index or each block type), use the same"+
			" index to specify redundant blocks e.g. use <BARCODE1:6> in two different layouts to specify that the barcode found in both reads are the same\n"+
			"   * 'X' : either a number indicating the length of the UMI, BARCODE or SAMPLE block or a negative number e.g. -2 to specify the last 2 bases"+
			" should be ignored/clipped) or the letter 'x' to specify to take the sequence till the end of read. Importantly, the 'x' or negative length shotcut "+
			"can only be used in the last block of a read layout (i.e. <BARCODE1:x><SAMPLE1:20> is not allowed)\n\n"
				)
	public List<String> READ_LAYOUT;

	/*
	 * Final ReadLayout created from READ_LAYOUT or defaults if READ_LAYOUT is not given 
	 */
	ReadLayout[] readLayouts = null;
	
	//output layouts
	@Option(shortName="OL", optional = false,
			printOrder=40,
			doc="Describes the output file layout(s) using the slots defined in read layouts and ':' to delimitate three parts e.g. 'OL=1:<BARCODE1><UMI1><UMI2>:<SAMPLE1>' : \n" +
					"\t"+"1.The number in the first part (i.e. from '1:' above) is the output file index and it must be unique across all 'OL' inputs. "+
					      "Inferred from order in comamnd line when not given\n"+
					"\t"+"2.The second part (i.e. '<BARCODE1><UMI1><UMI2>' above) is the read header layout ; when writing multiple UMI and BARCODE slots "+
					"in output read headers, these are always separated with the RCHAR (':' by defaults).\n"+
					"\t"+"3.The third part (i.e. '<SAMPLE1>' above) is the read sequence layout.\n"
			)
	public List<String> OUTPUT_LAYOUT;
	
	
	@Option(shortName="OWID",
			optional = true,
			printOrder=42,
			doc="Should the output layout number (output layout first slot) be injected in the filename ?\n"+
			    "Only used in absence of explicit file names in the barcode file.\n"
			)
	protected boolean ADD_LAYOUT_IDX_IN_OUTPUT_FILENAME = DEFAULT_ADD_LAYOUT_IDX_IN_OUTPUT_FILENAME;

	@Option(shortName="OWHL",
			optional = true,
			printOrder=44,
			doc="Should the output layout used for the read name (output layout second slot,in short format) be injected in the filename ? "+
			    "When true, each ouput file name contains e.g. '_B1U1' for OL='1:<BARCODE1><UMI1>:<SAMPLE1>'  \n"+
			    "Only used in absence of explicit file names in the barcode file.\n"
			)
	protected boolean ADD_HEADER_LAYOUT_IN_OUTPUT_FILENAME = DEFAULT_ADD_HEADER_LAYOUT_IN_OUTPUT_FILENAME;

	
	@Option(shortName="OWSL",
			optional = true,
			printOrder=46,
			doc="Should the output layout used for the read sequence (output layout third slot, in short format) be injected in the filename ?"+
				"When true, each ouput file name contains e.g. '_S1' for OL='1:<BARCODE1><UMI1>:<SAMPLE1>'  \n"+
			    "Only used in absence of explicit file names in the barcode file.\n"
			)
	protected boolean ADD_SEQUENCE_LAYOUT_IN_OUTPUT_FILENAME = DEFAULT_ADD_SEQUENCE_LAYOUT_IN_OUTPUT_FILENAME;

	
	/*
	 * Final FastqWriterLayout created from OUTPUT_LAYOUT or defaults if OUTPUT_LAYOUT is not given 
	 */
	FastqWriterLayout [] outLayouts = null;
	
		
	@Option(shortName = "O", optional = true,
			printOrder=90,
			doc="Output directory. By default, output files are written in running directory.\n")
	public File OUTPUT_DIR = null;

	@Option(optional = true,
			printOrder=100,
			doc="Allows to overwrite existing files (system rights still apply).\n"
			)
	public boolean FORCE = false;
	
	
	@Option(shortName="GZ", optional = true,
			printOrder=110,
			doc="Compress output files using gzip.\n"
			)
	public boolean GZIP_OUTPUTS = DEFAULT_GZIP_OUTPUTS;

	
	@Option(shortName = "SEP", optional = true,
			printOrder=170,
			doc="Separator character used to concatenate barcodes and umis in read header\n"
			)
	public String READ_NAME_SEPARATOR_CHAR = DEFAULT_READ_NAME_SEPARATOR_CHAR;
	
	
	@Option(shortName="TEST", optional = true,
			printOrder=210,
			doc="test mode ie code execution stops right before read demultiplexing starts btu after comamnd line validation"
			)
	protected static boolean TEST_MODE_STOP_AFTER_PARSING = false;

	@Option(optional=true,
			printOrder=220,
			doc="Change the default extension of created fastq files, eg 'fastqsanger'. By default uses the "
					+ "file extension from input fastq file. If result file names are given in the barcode file, "
					+ "this option is only used to adapt the unassigned file names. When using compression, a .gz is "
					+ "always appended to file names and should not be specified in FASTQ_FILE_EXTENSION i.e. \n"
					+ "use FASTQ_FILE_EXTENSION=fastq and NOT FASTQ_FILE_EXTENSION=fastq.gz\n"
			)
	public String FASTQ_FILE_EXTENSION = null;
	

	@Option(shortName="ASYNC", optional = true,
			printOrder=230,
			doc="Use one thread per Fastq Writer.\n")
	public boolean WRITER_FACTORY_USE_ASYNC_IO = DEFAULT_WRITER_FACTORY_USE_ASYNC_IO;

		
	
	/*
	 * The output files, one for each output layout
	 */
	List<File> outputFiles = new ArrayList<File>();
	
	
	
	
	@Override
	protected String[] customCommandLineValidation() {
		
		/*
		 * Check the input fastq
		 */
		Set<String> _names = new TreeSet<String>();
		for(File f : FASTQ){
			if(!f.exists() || !f.canRead())
				return new String[]{"Input FASTQ file does not exist OR cannot be read, please check: "+f.getAbsolutePath()};
			
			if(_names.contains(f.getAbsolutePath())){
				return new String[]{"Found twice the same file in FASTQ options: "+f.getAbsolutePath()};
			}
			_names.add(f.getAbsolutePath());
		}
			
		/*
		 * Validate O
		 * if not given, init to current dir else ensure the dir exists
		 * After this validation OUTDIR is SET 
		 */
		if(OUTPUT_DIR == null){
			OUTPUT_DIR = new File(System.getProperty("user.dir"));
		}

		if(!OUTPUT_DIR.exists()){
			log.info("Attempting to create output directory : "+OUTPUT_DIR.getAbsolutePath());
			try{
				FileUtil.checkWritableDir(OUTPUT_DIR, true);
			}catch(Exception e){
				return new String[]{"Failed to create output directory :"+OUTPUT_DIR.getAbsolutePath()};
			}
		}
		
		
		/*
		 * parse read layout, we must have one per input FASTQ
		 * 
		 */
		if(READ_LAYOUT.size() != FASTQ.size() )
			return new String[]{"Got "+READ_LAYOUT.size()+" read layouts for "+FASTQ.size()+" FASTQ files. You must provide as many read layouts as input FASTQ files."};

		log.debug("init ReadLayouts command line... ");

		readLayouts = new ReadLayout[READ_LAYOUT.size()];
		for (int j = 0; j < READ_LAYOUT.size(); j++) {
			String _rl = READ_LAYOUT.get(j);
			try{
				readLayouts[j] = new ReadLayout(_rl);
			}catch(Exception e){
				log.error(ExceptionUtil.getStackTrace(e));
				return new String[]{e.getMessage()};
			}
		}
		
		
		log.debug("ReadLayouts are valid.");
		
		/*
		 * Process output layout 
		 */
		if(OUTPUT_LAYOUT == null || OUTPUT_LAYOUT.isEmpty()){
			return new String[]{"Output layout(s) must be provided."};
		}

		log.debug("init output format layout from comamnd line... ");

		outLayouts = new FastqWriterLayout[OUTPUT_LAYOUT.size()];
		for (int j = 0; j < OUTPUT_LAYOUT.size(); j++) {
			String _ol = OUTPUT_LAYOUT.get(j);
			try{
				String [] parts = _ol.split(DEFAULT_READ_NAME_SEPARATOR_CHAR);
				if(parts.length > 3)
					return new String[]{"Invalid output layout. A maximum of three ':'-delimited parts are expected while "+parts.length+" are found in "+_ol};

				int olIndex = j+1;
				String headerLayout = "";
				String seqLayout = "";

				switch (parts.length) {
				case 1:
					seqLayout = parts[0];
					break;
				case 2:
					headerLayout = parts[0];
					seqLayout = parts[1];
					break;
				case 3:
					olIndex = Integer.parseInt(parts[0]);
					headerLayout = parts[1];
					seqLayout = parts[2];
					break;
				}
				outLayouts[j] = new FastqWriterLayout(seqLayout, headerLayout, readLayouts);
			}catch(Exception e){
				log.error(ExceptionUtil.getStackTrace(e));
				return new String[]{e.getMessage()};
			}
		}

		log.debug("Output format layouts are valid.");
		
		/*
		 * 
		 * Check one of  
		 * ADD_SEQUENCE_LAYOUT_IN_OUTPUT_FILENAME, 
		 * ADD_HEADER_LAYOUT_IN_OUTPUT_FILENAME, 
		 * ADD_LAYOUT_IDX_IN_OUTPUT_FILENAME 
		 * 
		 * is true 
		 * 
		 */
		
		if(!ADD_SEQUENCE_LAYOUT_IN_OUTPUT_FILENAME && !ADD_HEADER_LAYOUT_IN_OUTPUT_FILENAME && !ADD_LAYOUT_IDX_IN_OUTPUT_FILENAME){
			return new String[]{"At least one of 'OWID', 'OWHL' or 'OWSL' must be true"};
		}
		
		/*
		 * Properly init output files for samples
		 * Check if file name was provided in barcode file else create default names
		 */
		//no file name, init from output layouts
		log.debug("  initializing output file from output layouts ...");
		
		for (int i = 0; i < outLayouts.length; i++) {
			FastqWriterLayout ol = outLayouts[i];
			String fname = generateOutputFileName(i, ol);

			outputFiles.add( new File(OUTPUT_DIR , fname) );
		}

		log.debug("  output file names created, now checking if they already exist ");
		
		for (File en : outputFiles) {
			try{
				checkFilesDoNotAlreadyExist(en, FORCE);
			}catch(Jexception e){
				return new String[]{e.getMessage()};
			}
		}
		log.debug("out file path fully validated.");
		
		
		
		log.debug("Validation ended without error");
		return null; //do not return an empty array
	}


	/**
	 * @param i
	 * @param ol
	 * @return
	 */
	protected String generateOutputFileName(int i, FastqWriterLayout ol) {
		String fname = "out" + 
				( ADD_SEQUENCE_LAYOUT_IN_OUTPUT_FILENAME ? "_" +ol.getReadSequenceLayout() : "" )+
				( ADD_HEADER_LAYOUT_IN_OUTPUT_FILENAME ? "_" +ol.getReadNameLayout() : "" )+
				( ADD_LAYOUT_IDX_IN_OUTPUT_FILENAME ? "_" + (i+1) :"" )+
				"." +(FASTQ_FILE_EXTENSION == null ? "txt" : FASTQ_FILE_EXTENSION) +
				(GZIP_OUTPUTS ? ".gz" : "") ;
		return fname;
	}
	

	/**
	 * Checks if the file exists and can be overwritten. 
	 * 
	 * @param file
	 * @param overwrite
	 * @throws Jexception if the file exists and (1) overwrite is false or (2) overwrite is true and file cannot be deleted
	 */
	private void checkFilesDoNotAlreadyExist(File file, boolean overwrite) throws Jexception{
		if(file.exists() && !overwrite)
			throw new Jexception("Ouput file already exists : "+ file.getAbsolutePath()+"\nPlease delete file(s) first or use FORCE=true");
		else if(file.exists() && overwrite && !file.delete())
			throw new Jexception("Ouput file already exists but could not be deleted (file permission issue most likely): "+ file.getAbsolutePath()+"\nPlease delete file(s) first manually ro adapt file permissions.");
	}


	private void checkFilesDoNotAlreadyExist(List<File> files, boolean overwrite) {
		for (File file : files) {
			checkFilesDoNotAlreadyExist(file, overwrite);
		}
	}




	/**
	 * Checks if a given file name looks like a path
	 * @param fname
	 * @return
	 */
	protected boolean looksLikeAPath(String fname) {
		
		if(fname.contains(File.separator)){
			File f = new File(fname);
			return f.getParentFile().exists() && f.getParentFile().isDirectory();
		}
		return false;
	}

	
	
	private JemultiplexerFastqWriterFactory fastqFactory;


	@Override
	protected int doWork() {
		
		log.debug("Launching clipping...");
		
		try{
			
			/*
			 * Open writers for all output FASTQ files
			 *
			 */
			fastqFactory = new JemultiplexerFastqWriterFactory();
			fastqFactory.setUseAsyncIo(this.WRITER_FACTORY_USE_ASYNC_IO);
			//list to hold all sample writers
			List<FastqWriter> fastqWriters = new ArrayList<FastqWriter>(); 

			for (File _f : outputFiles) {
				fastqWriters.add(fastqFactory.newWriter(_f, GZIP_OUTPUTS, CREATE_MD5_FILE));
			}
				
			/*
			 * Open readers on all FASTQ
			 */
			List<FastqReader> fastqReaders = new ArrayList<FastqReader>(); //we need to store them to close them at the end
			List<Iterator<FastqRecord>> fastqFileIterators = new ArrayList<Iterator<FastqRecord>>();
			for (File fqFile : FASTQ) {
				FastqReader r = new FastqReader(fqFile);
				fastqReaders.add(r); 
				fastqFileIterators.add( r.iterator() );
			}
			
			
			/*
			 * Iterate over all records
			 */
				
			//counters
			Iterator<FastqRecord> mainIterator = fastqFileIterators.get(0);
			while(mainIterator.hasNext()){
				//reads next read from all input files
				FastqRecord[] reads = nextReads(fastqFileIterators);
			
				for (int i = 0; i < fastqWriters.size(); i++) {
					//prepare the output according to output layout
					log.debug("Writing in output idx "+(i+1));
					FastqRecord rec = outLayouts[i].assembleRecord( reads );
					fastqWriters.get(i).write(rec);
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
			for (FastqWriter w : fastqWriters) {
				try {
					w.close();
				} catch (Exception e) {
					// ignore
				}
			}

			log.debug("Clipping run ended without warning.");
			
		}catch(Exception e){
			log.error(ExceptionUtil.getStackTrace(e));
			log.error("\n\n\n\n\nAn error occurred during read clipping, please check the error message before running the process again. \n"+
			         "Error message was :\n"+e.getMessage());
			return 1;
		}
		return 0;
	}
	
	
	
	
	/**
	 * check that the indices form a series n(i+1) = n(i) + 1 
	 * @param ids
	 * @return return the max index
	 * @throws Jexception if the indices do not form a series n(i+1) = n(i) + 1 
	 */
	protected static int checkIndexSerie(String optionName, Set<Integer> ids)  {
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
			throw new Jexception("Indices for "+optionName+" option must start with number 1 (not "+min+")");

		// start at expected number according to id size ?
		int expectedEnd = min + ids.size() - 1;
		if(max != expectedEnd){
			ArrayList<Integer> _ids = new ArrayList<Integer>(ids);
			Collections.sort(_ids);
			throw new Jexception("Indices for the "+ids.size()+" "+optionName
					+" option do not form a continous sequence of integers from 1 to "+expectedEnd+" : "+StringUtil.mergeIntList(_ids,  ","));
		}
		
		return max;
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
	 * get real read length from fastq file (from first read)
	 * @param file
	 * @return the read length
	 */
	private int peekReadLength(File file) {
		
		final FastqReader r = new FastqReader(file);
		FastqRecord rec = null;
		try {
			rec = r.iterator().next();
		} finally {
			if(r!=null)
				r.close();
				
		}
		return rec.getReadLength();
	}
	
	
	
}
