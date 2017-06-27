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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.embl.cg.utilitytools.utils.ExceptionUtil;
import org.embl.cg.utilitytools.utils.FileUtil;
import org.embl.cg.utilitytools.utils.StringUtil;
import org.embl.cg.utilitytools.utils.parser.csv.CSVLine;
import org.embl.cg.utilitytools.utils.parser.csv.InvalidCSVSetUpException;
import org.embl.cg.utilitytools.utils.parser.csv.InvalidHeaderException;
import org.embl.gbcs.je.Jexception;
import org.embl.gbcs.je.ReadLayout;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;

public class Jedemultiplex extends CommandLineProgram {
	private static final String FASTQ_OPTION_SLOT_DELIMITER = ";";

	private static Logger log = LoggerFactory.getLogger(Jedemultiplex.class);

	/*
	 * We define here all common options
	 */
	
	@Option(shortName="F", optional = false,
			printOrder=10,
			doc="Input fastq file (optionally gzipped) with optional unique number and read layout " +
			"i.e. 'FASTQ=1;<BARCODE1:6><SAMPLE:x>;file.fastq.gz' ; with the ';' character as the delimiter." +
			"The first number (e.g.1:) is optional but highly recommended when more than one FASTQ input is provided (e.g. paired end situation). " +
			"When present it must be unique across all 'FASTQ' inputs. This number is used to uniquely refer to this input fastq" +
			" file and match up with the read layouts (see below). "+
			"When not present, the file index is inferred from the order on the command line."
//			"When not present, the file index is extracted from Illumina fastq read header i.e.:\n "+
//			"   @NB501764:235:HYL7KBGX2:1:11101:16625:1034 1:N:0:ATTCCT => '1' (first number after space)\n"+
//			"or inferred from the order on the command line"
			
			)
	public List<String> FASTQ;

	
	/*
	 * Barcode definition options
	 */
	@Option(shortName="BF", optional = false,  mutex={"USE_EMBASE"},
			printOrder=30,
			doc="Barcode file (tsv) matching sample names to barcode combination. \n" +
					"### SIMPLE Barcode File Format (for backward compatibility) ; please see the GENERAL format described below" + 
					"The file must have 2 columns with the sample in col1 and the corresponding barcode in col2.\n" +
					"In this format, a simple BARCODE slot is expected in the ReadLayout and NO headers are needed e.g. :\n" +
					"\t"+"\t"+"sample1\tGAGG\n" +
					"\t"+"\t"+"sample2\tCCAA\n" +
					"\t"+"The format accept the following shortcuts: \n" +
					"\t"+"1. If multiple barcodes map to the same sample, either line can be duplicated e.g.\n" +
					"\t"+"\t"+"sample1\tATAT\n" +
					"\t"+"\t"+"sample1\tGAGG\n" +
					"\t"+"\t"+"sample2\tCCAA\n" +
					"\t"+"\t"+"sample2\tTGTG\n" +
					"\t"+"Or barcodes can be combined using the OR operator '|' i.e. the file above can be re-written like\n " +
					"\t"+"\t"+"sample1\tATAT|GAGG\n" +
					"\t"+"\t"+"sample2\tCCAA|TGTG\n" +
					"\t"+"2. For the special situation of paired-end data in which barcodes differ at both ends i.e. with " +
					"BARCODE1 and BARCODE2 described for read one and two respectively, barcodes for BARCOD1 and BARCODE2 can be" +
					" distinguished using a ':' separator i.e. \n" +
					"\t"+"\t"+"sample1\tATAT:GAGG\n" +
					"\t"+"\t"+"sample2\tCCAA:TGTG\n" +
					"\t"+"This above syntax means that sample 1 is encoded with ATAT barcode from BARCODE1 slot AND GAGG barcode from BARCODE2 slot. " +
					"Note that you can still combine barcodes using | e.g. \n"+
					"\t"+"\t"+"sample1\tATAT|GAGG:CCAA|TGTG\n" +
					"\t"+"3. Extended barcode file format : 3 (single-end) or 4 (paired-end) tab-delimited colums\n"+
					"\t"+"same as the simple barcode file format but the extra columns contains the file name(s) to use to name output files." +
					" A unique extra column is expected for single-end while 2 extra columns are expected for paired-end. In case, lines are duplicated (multiple barcodes" +
					"mapping the same sample), the same file name should be indicated in the third (and fourth) column(s). \n"+
					"\t"+"\t"+"sample1\tATAT\tspl1_1.txt.gz\tspl1_2.txt.gz\n" +
					"\t"+"\t"+"sample1\tGAGG\tspl1_1.txt.gz\tspl1_2.txt.gz\n" +
					"\t"+"\t"+"sample2\tCCAA\tspl2_1.txt.gz\tspl2_2.txt.gz\n" +
					"\t"+"Or\n" +
					"\t"+"\t"+"sample1 \t ATAT|GAGG:CCAA|TGTG \t spl1_1.txt.gz \t spl1_2.txt.gz\n" +
					"### GENERAL Barcode File Format \n"+
					"In this format, the file structure is governed with headers:"+
					"\t"+"* the 'SAMPLE' column lists the sample names"+
					"\t"+"* the 'BARCODEn' columns list the matching BARCODE from the BARCODEn slot (where n is a number). "+
					"\t"+"    It is mandatory to have as many 'BARCODEn' columns as described BARCODE slots. Here again, barcodes can be combined using the OR operator '|'"+
					"\t"+"* the optional 'OUTn' columns (where n is a number) list the output file names for this sample and matching output number. "
			)
	public File BARCODE_FILE = null;
	
	

	/*
	 * This option cannot be used when using INDEX_FILE.
	 */
	@Option(shortName="EM", optional = false, mutex={"BARCODE_FILE"},
			printOrder=400,
			doc="Enables emBASE mode i.e fetch information from emBASE and place demultiplexed files directly in emBASE repository structure.\n" +
					"This option is mutually exclusive with BARCODE_FILE.\n" +
					"Note : this option forces O=null GZ=true UN=true UF1=null UF2=null STATS_ONLY=false (all other user options supported).\n" 
			)
	public boolean USE_EMBASE = false;
	
	
	@Option(shortName="RL", optional = true,
			printOrder=20,
			doc="Describes the read layout(s) i.e. 'RL=1:<BARCODE1:6><SAMPLE:x>' ; with the ':' character as the delimiter.\n" +
			"The first number (i.e. from '1:' above) is the read layout index and it must be unique across all 'RL' inputs as "+
			"it is used to match up the read layouts with the fastq files.\n"+
			"This number is optinal but highly recommended when more than one RL (or FASTQ) input are provided (e.g. paired end situation).\n" +
			"Read layouts are only needed for complex layouts and when the read layout(s) was/were not already embedded in the FASTQ files options.\n" +
			"This option must not be given when read layout were embedded in the FASTQ option"+
			"When not provided, the index is inferred from the comamnd line order"
			)
	public List<String> READ_LAYOUT;

	@Option(shortName="OL", optional = true,
			printOrder=20,
			doc="Describes the output file layout(s) using the slots defined in read layouts and ':' to delimitate three parts e.g. 'OL=1:<BARCODE1><UMI1><UMI2>:<SAMPLE1>' : \n" +
					"\t"+"1.The number in the first part (i.e. from '1:' above) is the output file index and it must be unique across all 'OL' inputs. "+
					      "Inferred from order in comamnd line when not given\n"+
					"\t"+"2.The second part (i.e. '<BARCODE1><UMI1><UMI2>' above) is the read header layout ; when writing multiple UMI and BARCODE slots "+
					"in output read headers, these are always separated with the RCHAR (':' by defaults).\n"+
					"\t"+"3.The third part (i.e. '<SAMPLE1>' above) is the read sequence layout.\n"+
			"One output file is created for each sample and each OL index. Output file names default to samplename_outputfileindex with the original extensions\n"+
			"### When no OL is described, Je considers that an output file should be created for each input FASTQ (containing a SAMPLE slot) and for each sample.\n "+
			"In this scenario:\n"+
					"\t"+"1. The output files only contain the BARCODE slots (concatenated if multiple BARCODE slots are described within the same read layout) unless CLIP is set to false\n"+
					"\t"+"2. The barcode(s) and sample names are injected in the output file names occording to the pattern 'FASTQFILENAMEn_SAMPLENAME_BARCODES.ORIGINALEXTENSIONS' ) \n"+
					"\t"+"3. All UMI slots (if any) are placed in the fastq headers following their slot index i.e. UMI1:UMI2:...:UMIn, separated with ':' (all UMIs are added to all "
							+ "file sets identically) ; unless ADD is et to false."
			
			)
	public List<String> OUTPUT_LAYOUT;

	
	@Override
	protected String[] customCommandLineValidation() {
		/*
		 * Parse and validate the barcode file
		 * we first check a valid file was given , then guess the format, convert if necessary and parse it
		 */
		if(!BARCODE_FILE.canRead()){
			return new String[]{
					"File is not readable :"+BARCODE_FILE.getAbsolutePath()
					};
		}
		
		File bcFileToParse = BARCODE_FILE;
		try {
			CSVLine l = FileUtil.readFirstValidLine(BARCODE_FILE.getAbsolutePath(), "\t", false); //hasHeaders to false to make sure to have the header line if existing 
			if(!Pattern.matches(BarcodeFileGeneralParser.headerLineRegex, l.merge("\t"))){
				log.debug("converting barcode file to general format...");
				BarcodeFileSimpleToGeneralFormatConverter conv = new BarcodeFileSimpleToGeneralFormatConverter(BARCODE_FILE);
				File tmp = File.createTempFile(
						FileUtil.removeExtension(BARCODE_FILE.getName()), 
						FileUtil.getExtension(BARCODE_FILE.getName())
						);
				conv.convertToGeneralFormat(tmp);
				bcFileToParse = tmp;
			}
		} catch (Exception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			log.error("Error while identifying the barcode file format (most likely during format convertion)");
			return new String[]{e.getMessage()};
		}
		
		/*
		 * parse the general format BC file
		 */
		
		BarcodeFileGeneralParser barcodeParser = new BarcodeFileGeneralParser(bcFileToParse);
		
		try {
			barcodeParser.parse();
		} catch (Exception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			log.error("Error while parsing the barcode file (general format) : "+bcFileToParse.getAbsolutePath());
			return new String[]{e.getMessage()};
		}
		
		//at this point we have a valid barcode file and  barcodeParser is fully initialized
		
		
		/*
		 * Check the input fastq; in particular extract associated read layouts.
		 * Also get the read length from each file and possible associated file idx e.g. from Illumina files :
		 * This is available in the second part, first number e.g.: 
		 * Read Header : 
		 * @NB501764:235:HYL7KBGX2:1:11101:16625:1034 1:N:0:ATTCCT    =====> '1'
		 * @NB501764:235:HYL7KBGX2:1:11101:16625:1034 2:N:0:ATTCCT    =====> '2'
		 */
		Map<Integer, ReadLayout> idx2layout = new HashMap<Integer, ReadLayout>();
		Map<Integer, File> idx2file = new HashMap<Integer, File>();
		for (int i = 0; i < FASTQ.size() ; i++) {
			String fileDescr = FASTQ.get(i);
			//make sure user did not give ie ";<readlayout>;file.txt" or ";;filename.gz" => make sure we dont start with a ;
			while(fileDescr.startsWith(FASTQ_OPTION_SLOT_DELIMITER))
				fileDescr = fileDescr.replace(FASTQ_OPTION_SLOT_DELIMITER, "");
			
			String [] parts = fileDescr.split(FASTQ_OPTION_SLOT_DELIMITER);
			if(parts.length > 3){
				return new String[]{"Comamnd line option FASTQ="+FASTQ.get(i)+" is not valid, a maximum of three slots delineated with '"+FASTQ_OPTION_SLOT_DELIMITER+"' is allowed)"};
			}
			//the three info to extract from the FASTQ option
			String filepath = fileDescr; //the only mandatory one => init with fileDescr in case no index nor read layout were given
			int fileindex = i; // init with position in command line, in case it is not given in the option
			ReadLayout rl = null;
			
			if(parts.length >= 2){
				//we have multiple parts, the file path is the last slot , always
				filepath = parts[parts.length - 1];
				
				//parts[0] CANNOT be empty (see above)
				try{
					fileindex = Integer.parseInt( parts[0] );
				}catch(NumberFormatException nfe){
					//no it is not; then it is a ReadLayout but only if user provided only 2 parts
					if(parts.length > 2){
						return new String[]{"Comamnd line option FASTQ="+FASTQ.get(i)+" is not valid: when three slots delineated with '"+FASTQ_OPTION_SLOT_DELIMITER+"' are used, the first one must be the file index)"};
					}
					try{
						rl = new ReadLayout(parts[0]);
					}catch(Exception e){
						log.error(ExceptionUtil.getStackTrace(e));
						return new String[]{e.getMessage()+"\nCannot make sense of command line option FASTQ="+FASTQ.get(i)+"! Please check the doc"};
					}
				}
			} 
			if(parts.length == 3){
				//then the file index was successfully parsed and parts[1] MUST be a ReadLayout or an empty slot ie '1;;file.txt'
				try{
					if(StringUtils.isNotBlank(parts[1]))
						rl = new ReadLayout(parts[1]);
				}catch(Exception e){
					log.error(ExceptionUtil.getStackTrace(e));
					return new String[]{e.getMessage()+"\nInvalid read layout in comamnd line option FASTQ="+FASTQ.get(i)+"! Please check the doc"};
				}
			}
			
			//save in tmp maps
			if(idx2file.containsKey(fileindex)){
				return new String[]{"All FASTQ index must be unique but at least 2 FASTQ options have been assigned with the index "+fileindex+"; for example: FASTQ="+FASTQ.get(i)};
			}
			idx2file.put(fileindex, new File(filepath));
			if(rl!=null)
				idx2layout.put(fileindex, rl);
		}
		//final check on indices
		int maxFastqIndex = -1;
		try{
			maxFastqIndex = checkIndexSerie("FASTQ", idx2file.keySet());
		}catch(Exception e){
			return new String[]{e.getMessage()};
		}
		
		List<ReadLayout> readLayouts = new ArrayList<ReadLayout>(FASTQ.size());
		List<File> fastqInputFiles = new ArrayList<File>(FASTQ.size());
		for (Entry<Integer, File> e : idx2file.entrySet()) {
			int indexInList = e.getKey() - 1;
			fastqInputFiles.set(indexInList, e.getValue());
			readLayouts.set(indexInList, idx2layout.get(e.getKey()));
		}
		
		/*
		 * if read layouts were given, check status with respect to above (read layout already given ?) ie
		 * - if redundant stop
		 * - if not given before, validate them (by creating ReadLayout)
		 * 
		 */
		
		/*
		 * If no read layouts were given, create defaults ones:
		 * - SE (1 fastq file given) : <BARCODE:n><SAMPLE:x> where n comes from barcode length found during barcode file parsing
		 * - PE (2 fastq given)      :   
		 *     ## option 1 : one 'barcode' column in barcode file => redundant barcode => <BARCODE:n><SAMPLE1:x> <BARCODE:n><SAMPLE2:x> 
		 *        where n is as above
		 *     ## option 2 : 2 'barcode' columns in barcode file => non redundant barcode => <BARCODE1:n1><SAMPLE1:x> <BARCODE2:n2><SAMPLE2:x>
		 *       where n1 and n2 is as n above. Only issue is how to know which barcode column matches which fastq file when n1==n2...? 
		 *       we'll assume first barcode column is for the fastq file idx 1 (from the command line or the extracted index)  
		 */
		
		
		/*
		 * Process output layout if provided. 
		 * If not provided, create default layouts according to default SE and PE situation or the number of file with a SAMPLE slot
		 */
		
		return null;
	}
	
	@Override
	protected int doWork() {
		Demultiplexer d = new Demultiplexer(fastqInFiles, readLayouts, sample2BarcodeSets, sample2outputFileList, outputFastqLayout, unassignedFiles, metricFile, fastqQualityFormat, gzipOutput, createMD5Sum);
		try{
			d.run(max_mismatches, min_mismatch_deltas, min_base_qualities, asyncWrite);
		}catch(Exception e){
			log.error(ExceptionUtil.getStackTrace(e));
			log.error("\n\n\n\n\nAn error occurred during demultiplexing, please check the error message before running the process again. \n"+
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
	

}
