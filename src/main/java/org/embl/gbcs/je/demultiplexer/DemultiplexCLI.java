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

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
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
import org.embl.gbcs.je.Jexception;
import org.embl.gbcs.je.ReadLayout;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

public class DemultiplexCLI extends CommandLineProgram {
	
	
	/*
	 * when writing all reads in same file, this is the sample name key used in the sample2outfiles map
	 */
	public static final String UNIQUE_MULTIPLEXED_SAMPLE_NAME = "MULTIPLEXED";

	private static Logger log = LoggerFactory.getLogger(DemultiplexCLI.class);

	/*
	 * Option defaults
	 */
	protected static final Integer DEFAULT_MAX_MISMATCHES = 1;
	protected static final Integer DEFAULT_MIN_MISMATCH_DELTA = 1;
	protected static final Integer DEFAULT_MIN_BASE_QUALITY = 10;
	protected static final Boolean DEFAULT_CLIP = true;
	protected static final Boolean DEFAULT_ADD = true;
	protected static final FastqQualityFormat DEFAULT_QUALITY_FORMAT = null; 
	protected static final Boolean DEFAULT_GZIP_OUTPUTS = true;
	protected static final Boolean DEFAULT_KEEP_UNASSIGNED_READ = true;
	protected static final Boolean DEFAULT_WRITER_FACTORY_USE_ASYNC_IO = true;
	protected static final String DEFAULT_UNASSIGNED_FILE_NAME_BASE = "unassigned";
	protected static final String DEFAULT_METRICS_FILE_NAME = "jemultiplexer_out_stats.txt";
	protected static final String DEFAULT_READ_NAME_SEPARATOR_CHAR = ":";
	protected static final String DEFAULT_USE_ORIGINAL_BARCODE_SEQUENCE = "0";
	protected static final boolean DEFAULT_STRICT = false;
	
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

	
	/*
	 * Barcode file
	 */
	@Option(shortName="BF", optional = false,  mutex={"USE_EMBASE"},
			printOrder=10,
			doc="Barcode file (tsv) matching sample names to barcode combination. \n\n" +
					"   ### GENERAL Barcode File Format \n"+
					"In this format, the file structure is governed with headers:\n"+
					"\t"+"* the 'SAMPLE' column lists the sample names\n"+
					"\t"+"* the 'BARCODEn' columns list the matching BARCODE from the BARCODEn slot (where n is a number, see RL option). \n"+
					"\t"+"    It is mandatory to have as many 'BARCODEn' columns as described BARCODE slots in READ LAYOUTS. Here again, barcodes can be combined using the OR operator '|'\n"+
					"\t"+"* the optional 'OUTn' columns (where n is a number) list the output file names for this sample and matching output number.\n\n"+
					"   ### SIMPLE Barcode File Format (for backward compatibility) ; please see the GENERAL format described above \n" + 
					"The file must have 2 columns with the sample in col1 and the corresponding barcode in col2.\n" +
					"In this format, a simple BARCODE slot is expected in the ReadLayout and NO headers are needed e.g. :\n" +
					"\t"+"\t"+"sample1\tGAGG\n" +
					"\t"+"\t"+"sample2\tCCAA\n" +
					"\t"+"The format accept the following shortcuts: \n" +
					"\t"+"1. If multiple barcodes map to the same sample, either lines can be duplicated e.g.\n" +
					"\t"+"\t"+"sample1\tATAT\n" +
					"\t"+"\t"+"sample1\tGAGG\n" +
					"\t"+"\t"+"sample2\tCCAA\n" +
					"\t"+"\t"+"sample2\tTGTG\n" +
					"\t"+"Or barcodes can be combined using the OR operator '|' i.e. the file above can be re-written like\n " +
					"\t"+"\t"+"sample1\tATAT|GAGG\n" +
					"\t"+"\t"+"sample2\tCCAA|TGTG\n" +
					"\t"+"2. For the special situation of paired-end data in which barcodes differ at both ends i.e. with " +
					"BARCODE1 and BARCODE2 described for read one and two respectively, barcodes for BARCODE1 and BARCODE2 can be" +
					" distinguished using a ':' separator i.e. \n" +
					"\t"+"\t"+"sample1\tATAT:GAGG\n" +
					"\t"+"\t"+"sample2\tCCAA:TGTG\n" +
					"\t"+"This above syntax means that sample 1 is encoded with ATAT barcode from BARCODE1 slot AND GAGG barcode from BARCODE2 slot. " +
					"Note that you can still combine barcodes using | e.g. \n"+
					"\t"+"\t"+"sample1\tATAT|GAGG:CCAA|TGTG\n" +
					"\t"+"3. Extended barcode file format : 3 (single-end) or 4 (paired-end) tab-delimited colums\n"+
					"\t"+"same as the simple barcode file format but the extra columns contains the file name(s) to use to name output files." +
					" A unique extra column is expected for single-end while 2 extra columns are expected for paired-end. In case lines are duplicated (multiple barcodes " +
					"mapping the same sample), the same file name should be indicated in the third (and fourth) column(s). \n"+
					"\t"+"\t"+"sample1\tATAT\tspl1_1.txt.gz\tspl1_2.txt.gz\n" +
					"\t"+"\t"+"sample1\tGAGG\tspl1_1.txt.gz\tspl1_2.txt.gz\n" +
					"\t"+"\t"+"sample2\tCCAA\tspl2_1.txt.gz\tspl2_2.txt.gz\n" +
					"\t"+"Or\n" +
					"\t"+"\t"+"sample1 \t ATAT|GAGG:CCAA|TGTG \t spl1_1.txt.gz \t spl1_2.txt.gz\n" 
					
			)
	public File BARCODE_FILE = null;
	
	
	/**
	 * Barcode file parser initialised with BARCODE_FILE or the barcode file exported from emBASE
	 */
	protected BarcodeFileGeneralParser barcodeParser ;

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
	
	
	//Read layouts 
	@Option(shortName="RL", optional = true,
			printOrder=20,
			doc="Describes the read layout(s) of input fastq file(s) e.g. RL='<BARCODE:6><SAMPLE:x>' describes a read with a barcode in the first 6 bases " +
			"followed by the sample sequence ('x' means 'till the end', see below). You MUST single quote the pattern (RL='<BARCODE:6><SAMPLE:x>') as '>' have special meaning in unix."+
			"The input fastq files and read layouts are mached up by order on the command line.\n" +
			"Read layouts are only needed for complex layouts but one must provide read layouts for ALL or NONE of the input fastq files.\n"+
			"## READ LAYOUT FORMAT DESCRIPTION:/n"+
			"Read layouts are made of <UMIn:X>, <BARCODEn:X>, <SAMPLEn:X> blocks to describe blocks of type UMI, BARCODE or SAMPLE with : \n "+
			"   * 'n' the unique block index (an index must be unique across all read layouts for each index or each block type), use the same"+
			" index to specify redundant blocks e.g. use <BARCODE1:6> in two different layouts to specify that the barcode found in both reads are the same\n"+
			"   * 'X' : either a number indicating the length of the UMI, BARCODE or SAMPLE block or a negative number e.g. -2 to specify the last 2 bases"+
			" should be ignored/clipped) or the letter 'x' to specify to take the sequence till the end of read. Importantly, the 'x' or negative length shortcuts "+
			"can *only* be used in the last block of a read layout (i.e. <BARCODE1:x><SAMPLE1:20> is not allowed)\n" +
			"In addition, layouts can contain N or fixed bases like in 'NN<BARCODE1:6>NNNN<SAMPLE1:x>' where the Ns tell Je to skip 2 and 4 bases before" +
			" extracting the barcode & sample sequence respectively.\n\n"+
			"## OMIITING READ LAYOUT IN THE COMMAND LINE:/n"+
			"When no read layout is provided, the following defaults apply :\n"+
			"   * 1 input fastq: single end with layout <BARCODE1:X><SAMPLE1:x> where X is inferred from barcode file\n"+
			"   * 2 input fastqs: \n"+
	        "       - paired end with redundant barcode if barcode file describes a single BARCODE i.e. <BARCODE1:X><SAMPLE1:x> and <BARCODE1:X><SAMPLE2:x>, where X is inferred from barcode file\n"+
	        "       - paired end with non redundant barcode if barcode file describes two BARCODE column i.e. <BARCODE1:X><SAMPLE1:x> and <BARCODE2:Y><SAMPLE2:x>, where X and Y are inferred from barcode file\n"+
	        "       - single end with index file if barcode file describes a single BARCODE and second fastq file has reads of length < 10 + barcode_length\n"+
			"   * 3 input fastqs: \n"+
	        "       - paired end with an index file i.e. <SAMPLE1:x>, <SAMPLE2:x> and <BARCODE1:X> when barcode file has a single BARCODE column (X is inferred from barcode file)\n"+
			"       - single end with two index files i.e. <SAMPLE1:x>, <BARCODE1:X> and <BARCODE2:Y> when barcode file has two BARCODE columns (X,Y is inferred from barcode file)\n"+
			"   * 4 input fastqs: paired end with either \n "+
			"       - 2 non-redundant index files i.e. <SAMPLE1:x>, <SAMPLE2:x>, <BARCODE1:X>, <BARCODE2:Y> if the barcode file has two BARCODE columns or a ATGC:GCTAGC syntax (X,Y inferred from barcode file) \n"+
			"       - 2 redundant index files <SAMPLE1:x>, <SAMPLE2:x>, <BARCODE1:X> and <BARCODE1:X> if the barcode file has a single BARCODE column (X inferred from barcode file)\n\n"
			)
	public List<String> READ_LAYOUT;

	/*
	 * Final ReadLayout created from READ_LAYOUT or defaults if READ_LAYOUT is not given 
	 */
	ReadLayout[] readLayouts = null;
	
	//output layouts
	@Option(shortName="OL", optional = true,
			printOrder=30,
			doc="Describes the output file layout(s) using the slots defined in read layouts (<UMIn>, <BARCODEn>, <SAMPLEn>) and are made of "+
			"three distinct parts separated with ':'.\n"+
			"In addition to <UMIn>, <BARCODEn>, <SAMPLEn>, <READBARn> is used as a synonym to <BARCODEn> to indicate that the real sequence "+
			" should be written (as opposed to writting the barcode when usign <BARCODEn>). \n"+
			"An output layout looks like '1:<BARCODE1><UMI1><UMI2>:<SAMPLE1>' where the three mandatory parts (':'-separated) are :\n"+
					"\t"+"- The number in the first part (i.e. from '1:' above) is the output file index and it must be unique across all 'OL' inputs. \n"+
					"\t"+"- The second part (i.e. '<BARCODE1><UMI1><UMI2>' above) is the read header layout; when writing multiple UMI and BARCODE slots "+
					"in output read headers, these are always separated with the RCHAR (':' by defaults).\n"+
					"\t"+"- The third part (i.e. '<SAMPLE1>' above) is the read sequence layout. Note that here <BARCODEn> and <READBARn> are fully synonyms as the real" +
					" sequence (i.e READBAR) is always written\n\n"+
			"Important: You MUST single quote the pattern (OL='1:<BARCODE1><UMI1><UMI2>:<SAMPLE1>') as '>' have special meaning in unix."+
			"An output file is created for each sample and each OL index. Output file names default to samplename_outputfileindex with the original fastq file extensions\n\n"+
			"## OMIITING OUTPUT LAYOUT IN THE COMMAND LINE:/n"+
			"  When no OL is described, Je considers an output file should be created for each input FASTQ (containing a SAMPLE slot) and for each sample.\n "+
			"In this scenario:\n"+
					"\t"+"1. The output files only contain the SAMPLE slot unless CLIP is set to false\n"+
					"\t"+"2. The barcode(s) and sample names are injected in the output file names according to the pattern 'FASTQFILENAMEn_SAMPLENAME_BARCODES.ORIGINALEXTENSIONS' ) \n"+
					"\t"+"3. Unless ADD is set to false, all BARCODE and UMI slots (if any) are placed in the fastq headers following their slot index i.e. BARCODE1:...:BARCODEn:UMI1:UMI2:...:UMIn and are separated "
							+ "with ':'.\n"+
			"## SHORT LAYOUT FORMAT \n"+
			"The output layout can be specified in a concise way using 'S','B', 'R' and 'U' for SAMPLE, BARCODE, READBAR and UMI, respectively. In this format,"+
			" the surounding '<>' are also omitted. For example 'OL=1:B1U1U2:S1' is a synonym of 'OL=1:<BARCODE1><UMI1><UMI2>:<SAMPLE1>'"
			
			)
	public List<String> OUTPUT_LAYOUT;
	/*
	 * Final FastqWriterLayout created from OUTPUT_LAYOUT or defaults if OUTPUT_LAYOUT is not given 
	 */
	protected FastqWriterLayout [] outLayouts = null;
	
	
	@Option(shortName="WQ", optional = true,
			printOrder=35,
			doc="Set to True to keep Phred sequence qualities in output read names. \n"+
			"This option only applies to BARCODE, READBAR and UMI described in the read name slot of output layout. "+
			"For BARCODE, the equivalent READBAR quality is used. In case of redundant slots, the best found quality is used.\n"+
			"The quality string is translated into 2 digits number representing the quality scores on the Phred scale and a e.g. UMI will look like\n"+
					"\t"+" '...:ATGCAT333023212322:...' instead of '...:ATGCAT:...'\n"+
			"This option is particularly useful with the retag module that knows how to extract quality numbers into BAM tags."
				)
	public boolean WITH_QUALITY_IN_READNAME = false;

	
	@Option(shortName = "O", optional = true,
			printOrder=40,
			doc="Output directory. By default, output files are written in running directory.\n")
	public File OUTPUT_DIR = null;

	@Option(shortName="UN", optional = true,
			printOrder=50,
			doc="Should un-assigned reads be saved in files or simply ignored. File names are automatically created"
					+ " or can be given using UF option.\n"
				)
	public boolean KEEP_UNASSIGNED_READ = DEFAULT_KEEP_UNASSIGNED_READ;


	@Option(shortName="UF", optional=true, 
			printOrder=60,
			doc="Name of unassigned files in which to write unassigned reads. When provided, Je expects as many UF files "
					+ "as input FASTQ files. UF options are matched up with FASTQ options following the order they are defined on the command line.\n" 
					+"Either a name (in which case the file will be created in the output dir) or full path.\n"
			)
	public List<String> UNASSIGNED_FILE = null;
	
	//unassignedFiles are unassignedFilePathes turned into File after validation of unassignedFilePathes
	protected List<File> unassignedFiles = null;

	
	@Option(shortName="OWID",
			optional = true,
			printOrder=70,
			doc="Should the output layout number (output layout first slot) be injected in the filename ?\n"+
			    "Only used in absence of explicit file names in the barcode file.\n"
			)
	protected boolean ADD_LAYOUT_IDX_IN_OUTPUT_FILENAME = DEFAULT_ADD_LAYOUT_IDX_IN_OUTPUT_FILENAME;

	@Option(shortName="OWHL",
			optional = true,
			printOrder=75,
			doc="Should the output layout used for the read name (output layout second slot,in short format) be injected in the filename ? "+
			    "When true, each ouput file name contains e.g. '_B1U1' for OL='1:<BARCODE1><UMI1>:<SAMPLE1>'  \n"+
			    "Only used in absence of explicit file names in the barcode file.\n"
			)
	protected boolean ADD_HEADER_LAYOUT_IN_OUTPUT_FILENAME = DEFAULT_ADD_HEADER_LAYOUT_IN_OUTPUT_FILENAME;

	
	@Option(shortName="OWSL",
			optional = true,
			printOrder=80,
			doc="Should the output layout used for the read sequence (output layout third slot, in short format) be injected in the filename ?"+
				"When true, each ouput file name contains e.g. '_S1' for OL='1:<BARCODE1><UMI1>:<SAMPLE1>'  \n"+
			    "Only used in absence of explicit file names in the barcode file.\n"
			)
	protected boolean ADD_SEQUENCE_LAYOUT_IN_OUTPUT_FILENAME = DEFAULT_ADD_SEQUENCE_LAYOUT_IN_OUTPUT_FILENAME;

	
	
	@Option(shortName="OF", optional=true, 
			printOrder=90,
			doc="Tells Je to write **all** assigned reads in the same output file(s) i.e. use this option when you do NOT want "+
			"to create per-sample demultiplexed files but rather want to keep all reads in the same file while barcode information"+
			" is gathered and injected in output formats.\n"+
			" When provided, Je expects as many 'OF=<file>' as output layouts ('OL=...') parameters or 'FASTQ=input' files when OL is not provided\n."+
			" OF options are matched up with OL/FASTQ options following the order in which they are defined on the command line.\n"+ 
			"OF expects either a name (in which case the file will be created in the output dir) or an absolute path.\n"
			)
	public List<String> OUTPUT_FILE = null;
	
	//multiplexedOutFiles are OUTPUT_FILE turned into File after validation of OUTPUT_FILE
	protected List<File> multiplexedOutFiles = null;
	protected boolean WRITE_ALL_READS_IN_SAME_OUTPUT = false; //turned to true if multiplexedOutFiles is used
	
	
	
	@Option(shortName = "MM", optional = true,
			printOrder=100,
			doc="Maximum mismatches for a barcode to be considered a match. Either exactly one or multiple values (with format MM=X:Y:Z). \n"+
					"When multiple values are provided, Je expects exactly one value for each BARCODE (with distinct indices) described in the barcode file/read layouts.\n"+
					"Values (X,Y,Z) are matched up with the sorted list of BARCODES (i.e.  X for BARCODE1, Y for BARCODE2 and Z for BARCODE3)\n"
			)
	public String MAX_MISMATCHES = DEFAULT_MAX_MISMATCHES.toString();

	// array  associated with MAX_MISMATCHES (after parsing) ; with a length that of the different barcode number found in the barcode file 
	int [] max_mismatches ; 
	
	
	@Option(shortName = "MMD", optional = true,
			printOrder=110,
			doc="Minimum difference between the number of mismatches against the best and the second best barcode. When MMD is not respected, "+
					"the read remains unassigned.\n" +
					"Either exactly one or multiple values (with format MMD=X:Y:Z). "+
					"When multiple values are provided, Je expects exactly one value for each BARCODE (with distinct indices) described in the barcode file/read layouts.\n"+
					"Values (X,Y,Z) are matched up with the sorted list of BARCODES (i.e.  X for BARCODE1, Y for BARCODE2 and Z for BARCODE3)\n"
			)
	public String MIN_MISMATCH_DELTA = DEFAULT_MIN_MISMATCH_DELTA.toString();
	
	// array  associated with MIN_MISMATCH_DELTA (after parsing) ; with a length that of the different barcode number found in the barcode file
	int [] min_mismatch_deltas;

	@Option(shortName="Q", optional = true,
			printOrder=120,
			doc="Minimum base quality during barcode matching: bases which quality is less than this cutoff are always considered as a mismatch." +
					"Either exactly one or multiple values (with format Q=X:Y:Z). "+
					"When multiple values are provided, Je expects exactly one value for each BARCODE (with distinct indices) described in the barcode file/read layouts.\n"+
					"Values (X,Y,Z) are matched up with the sorted list of BARCODES (i.e.  X for BARCODE1, Y for BARCODE2 and Z for BARCODE3)\n"
			)
	public String MIN_BASE_QUALITY = DEFAULT_MIN_BASE_QUALITY.toString();

	// array  associated with MIN_BASE_QUALITY (after parsing) ; with a length that of the different barcode number found in the barcode file
	int [] min_base_qualities;

	@Option(shortName="S", optional = true,
			printOrder=130,
			doc="When reads have redundant BARCODE slots, this option tells how to handle situation when the read sequence do not resolve to the same sample.\n"
					+" When true, the read pair is always 'unassigned'.\n"
					+" When false, the read pair is assigned to the sample with the lowest overall mismatch sum\n"
			)
	public boolean STRICT = false;
	
	
	@Option(optional = true,
			printOrder=140,
			doc="Allows to overwrite existing files (system rights still apply).\n"
			)
	public boolean FORCE = false;
	
	
	@Option(shortName="GZ", optional = true,
			printOrder=150,
			doc="Compress output files using gzip.\n"
			)
	public boolean GZIP_OUTPUTS = DEFAULT_GZIP_OUTPUTS;

	
	
	//in absence of output layout, should we clip barcodes and UMIs from sequence ? 
	@Option(optional = true,
			printOrder=160,
			doc="In absence of output layout, tell if barcode and UMI sequences should be clipped off read sequence before writing to output file.\n"
					+" If false, reads are written without modification to output file."
			)
	public boolean CLIP = DEFAULT_CLIP;

	@Option(optional = true,
			printOrder=170,
			doc="In absence of output layout, tell if barcode and UMI sequences should be added at the end of the read header.\n"
					+"BARCODE and UMI slots (in this order) are concatenated using the character defined by the SEP option\n"
			)
	public boolean ADD = DEFAULT_ADD;

	
	
	@Option(shortName = "SEP", optional = true,
			printOrder=180,
			doc="Separator character used to concatenate barcodes and umis in read header\n"
			)
	public String READ_NAME_SEPARATOR_CHAR = DEFAULT_READ_NAME_SEPARATOR_CHAR;


	
	@Option(shortName="V", optional = true,
			printOrder=190,
			doc="A value describing how the quality values are encoded in the fastq files.  Either 'Solexa' for pre-pipeline 1.3 " +
					"style scores (solexa scaling + 66), 'Illumina' for pipeline 1.3 and above (phred scaling + 64) or 'Standard' for phred scaled " +
					"scores with a character shift of 33.  If this value is not specified (or 'null' is given), the quality format is assumed to be will the 'Standard' for phred scale.\n"
			)
	public FastqQualityFormat QUALITY_FORMAT = DEFAULT_QUALITY_FORMAT;

	
	@Option(shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME, optional = true,
			printOrder=200,
			doc="File name where to write demultiplexing statistics. "
					+"Either a name (in which case the file will be created in the output dir) or an absolute path.\n"
			)
	public String METRICS_FILE_NAME = DEFAULT_METRICS_FILE_NAME;
	
	File metricsFile = null; //final File object to use with METRICS_FILE_NAME, set in cmdline validation 

	
	@Option(shortName = "DIAG", optional = true,
			printOrder=210,
			doc="Name for a barcode match reporting file (not generated by default)."
					+"Either a name (in which case the file will be created in the output dir) or full path."
					+" This file will contain a line per read set with the barcodes best matching the read subsequences "
					+"or 'null' when no match is found according to matching parameters ; and the final selected sample."
					+" This file is useful for debugging or further processing in case both ends are barcoded.\n"
					)
	String BARCODE_DIAG_FILE = null;
	File bcDiagFile = null; //final File object to use with BARCODE_DIAG_FILE, set in cmdline validation if necessary

	
	@Option(shortName="TEST", optional = true,
			printOrder=220,
			doc="test mode ie code execution stops right before read demultiplexing starts but after command line validation"
			)
	protected static boolean TEST_MODE_STOP_AFTER_PARSING = false;

	@Option(optional=true,
			printOrder=230,
			doc="Change the default extension of created fastq files, eg 'fastqsanger'. By default uses the "
					+ "file extension from input fastq file. If result file names are given in the barcode file, "
					+ "this option is only used to adapt the unassigned file names. When using compression, a .gz is "
					+ "always appended to file names and should not be specified in FASTQ_FILE_EXTENSION i.e. \n"
					+ "use FASTQ_FILE_EXTENSION=fastq and NOT FASTQ_FILE_EXTENSION=fastq.gz\n"
			)
	public String FASTQ_FILE_EXTENSION = null;
	
	@Option(optional=true,
			printOrder=240,
			doc="Indicates if the input fastq files are gzipped. Please use this option when file names are compressed but lack the typical '.gz' extension. \n"
			)
	public Boolean INPUT_FASTQ_COMPRESSION = null;
	

	@Option(shortName="ASYNC", optional = true,
			printOrder=240,
			doc="Use one thread per Fastq Writer.\n")
	public boolean WRITER_FACTORY_USE_ASYNC_IO = DEFAULT_WRITER_FACTORY_USE_ASYNC_IO;

	
	
	/**
	 * List of output fastq files for each multiplexed sample
	 */
	protected Map<String, List<File>> sample2outfiles = null;
	
	
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
		 * Check quality format
		 */
		if (QUALITY_FORMAT == null) { // we expect Standard 
			FastqReader [] readers = new FastqReader[FASTQ.size()];
			int i = 0;
			for(File f : FASTQ){
				readers[i++] = new FastqReader(f);
			}
			QUALITY_FORMAT = JeUtils.determineQualityFormat(readers, FastqQualityFormat.Standard);
			log.info( String.format("Auto-detected quality encoding format as '%s'. Please set V option explicitely if not correct.", QUALITY_FORMAT) );
		} else {
			log.info(String.format("Quality encoding format set to %s by user.", QUALITY_FORMAT));
		}
		
		
		if(USE_EMBASE){
			List<String> messages = customEmbaseCommandLineValidation(); //SET the BARCODE_FILE to the one exported from embase
			if(messages!=null && messages.size() > 0)
				return messages.toArray(new String[messages.size()]);
			
		}
		else{
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
		}
		
		
		/*
		 * Validate UN and UF
		 * We created one output per input fastq  
		 */
		if(KEEP_UNASSIGNED_READ){
			//did we have file names on command line ?
			unassignedFiles = new ArrayList<File>();

			if(UNASSIGNED_FILE==null || UNASSIGNED_FILE.isEmpty()){
				//init with default names
				for(File f : FASTQ){
					String name = 
							FileUtil.removeExtension( f.getName() ) +
							"_"+ DEFAULT_UNASSIGNED_FILE_NAME_BASE +
							"." +(FASTQ_FILE_EXTENSION == null ? "txt" : FASTQ_FILE_EXTENSION) +
							(GZIP_OUTPUTS ? ".gz" : "") ;

					//register
					unassignedFiles.add(new File(OUTPUT_DIR, name));
				}
			}
			else{
				if(UNASSIGNED_FILE.size() > 0 && UNASSIGNED_FILE.size() != FASTQ.size() )
					return new String[]{"Only "+UNASSIGNED_FILE.size()+" UF options found while "+FASTQ.size()+" FASTQ input files were give. FASTQ and UF option number must be the same."};

				for (String p : UNASSIGNED_FILE) {
					if(FASTQ_FILE_EXTENSION!=null && !p.contains("."+FASTQ_FILE_EXTENSION)){
						String newP = FileUtil.removeExtension(p) + "." +FASTQ_FILE_EXTENSION + (GZIP_OUTPUTS ? ".gz" : "");
						log.warn("Provided UF option "+p+" does not match with provided FASTQ_FILE_EXTENSION "+FASTQ_FILE_EXTENSION+". Changing file name to "+newP);
						p = newP;
					}
					File f = null;
					if(GZIP_OUTPUTS && !p.endsWith(".gz"))
						p += ".gz";

					f = (looksLikeAPath(p) ? new File(p) : new File(OUTPUT_DIR, p) );

					//register
					unassignedFiles.add(f);
				}

			}
		}

		try{
			checkFilesDoNotAlreadyExist(unassignedFiles, FORCE);
		}catch(Jexception e){
			return new String[]{e.getMessage()};
		}
		
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
						FileUtil.getExtension("."+BARCODE_FILE.getName())
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
		
		barcodeParser = new BarcodeFileGeneralParser(bcFileToParse);
		
		try {
			barcodeParser.parse();
		} catch (Exception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			log.error("Error while parsing the barcode file (general format) : "+bcFileToParse.getAbsolutePath());
			return new String[]{e.getMessage()};
		}
		
		//at this point we have a valid barcode file and  barcodeParser is fully initialized
		log.debug("valid barcode file and  barcodeParser is now initialized");
		
		/*
		 * Validate MM, MMD, Q
		 */
		int nBarcodes = barcodeParser.getBarcodeColumnNumber();
		try {
			max_mismatches = parseIntegerOptions("MAX_MISMATCHES", MAX_MISMATCHES, nBarcodes);
			min_mismatch_deltas = parseIntegerOptions("MIN_MISMATCH_DELTA", MIN_MISMATCH_DELTA, nBarcodes);
			min_base_qualities = parseIntegerOptions("MIN_BASE_QUALITY", MIN_BASE_QUALITY, nBarcodes);
		} catch (Jexception e) {
			return new String[]{e.getMessage()};
		}
		log.debug("got valid MM, MMD, Q options");
		
		/*
		 * if CLIP=true => ensure we have no barcode mix ie AGTCGA|AGGGCT from barcode file
		 * if we do , we turn it to false with a warning
		 */
		if(CLIP && barcodeParser.useBarcodeMix()){
			log.warn("CLIP cannot be set to TRUE when barcode file defines some barcode mix e.g. AGTCGA|AGGGCT (fake example) => Setting back CLIP to FALSE");
			CLIP = false;
		}
		
		/*
		 * parse read layout or initialize default ones
		 * 
		 */
		
		DefaultLayouts defaultLayouts = new DefaultLayouts(FASTQ, barcodeParser, CLIP, ADD, WITH_QUALITY_IN_READNAME, READ_NAME_SEPARATOR_CHAR, QUALITY_FORMAT);
		if(this.READ_LAYOUT == null || this.READ_LAYOUT.isEmpty()){
			if(FASTQ.size() > 4){
				return new String[]{"You must provide read layouts for each provided input FASTQ files."};
			}
			log.debug("init ReadLayouts from defaults... ");
			readLayouts = defaultLayouts.getReadLayouts();
		}
		else{
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
		}
		
		log.debug("ReadLayouts are valid.");
		
		/*
		 * Process output layout if provided. 
		 * If not provided, create default layouts according to default SE and PE situation or the number of file with a SAMPLE slot
		 */
		boolean USER_PROVIDED_OUTPUT_LAYOUT = true;
		if(this.OUTPUT_LAYOUT == null || this.OUTPUT_LAYOUT.isEmpty()){
			if(FASTQ.size() > 4){
				return new String[]{"You must provide output layouts when more than 4 FASTQ files are provided."};
			}
			log.debug("init output format layout from defaults... ");
			outLayouts = defaultLayouts.getFasqWriterLayouts();
			USER_PROVIDED_OUTPUT_LAYOUT = false;
		}
		else{
			if(READ_LAYOUT == null || READ_LAYOUT.isEmpty()){
				return new String[]{"Custom output layouts can only be provided when read layouts are also provided."};
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
					outLayouts[j] = new FastqWriterLayout(seqLayout, headerLayout, readLayouts, this.WITH_QUALITY_IN_READNAME, this.READ_NAME_SEPARATOR_CHAR, this.QUALITY_FORMAT);
				}catch(Exception e){
					log.error(ExceptionUtil.getStackTrace(e));
					return new String[]{e.getMessage()};
				}
			}
		}
		log.debug("Output format layouts are valid.");
		
		/*
		 * Check if user wants all reads in the same output file(s)
		 * 
		 */
		if(OUTPUT_FILE!=null && !OUTPUT_FILE.isEmpty()) {
			multiplexedOutFiles = new ArrayList<File>();
			//turn this global flag to true
			WRITE_ALL_READS_IN_SAME_OUTPUT = true;
			/*
			 * now check what user gave us.
			 * 1. If OL were provided, we need the same number of entries in multiplexedOutFilePathes 
			 */
			if(USER_PROVIDED_OUTPUT_LAYOUT) {
				if(OUTPUT_FILE.size() != outLayouts.length) {
					return new String[]{"You must provide as many OUTPUT_FILE (OF) options as OUTPUT_LAYOUT (OL) options."};
				}
			}else {
				if(OUTPUT_FILE.size() != FASTQ.size()) {
					return new String[]{"You must provide as many OUTPUT_FILE (OF) options as FASTQ input files."};
				}
			}
			//numbers match, initialize the file File objects for each OUTPUT_FILE
			for (String p : OUTPUT_FILE) {
				if(FASTQ_FILE_EXTENSION!=null && !p.contains("."+FASTQ_FILE_EXTENSION)){
					String newP = FileUtil.removeExtension(p) + "." +FASTQ_FILE_EXTENSION + (GZIP_OUTPUTS ? ".gz" : "");
					log.warn("Provided OF option "+p+" does not match with provided FASTQ_FILE_EXTENSION "+FASTQ_FILE_EXTENSION+". Changing file name to "+newP);
					p = newP;
				}
				File f = null;
				if(GZIP_OUTPUTS && !p.endsWith(".gz"))
					p += ".gz";

				f = (looksLikeAPath(p) ? new File(p) : new File(OUTPUT_DIR, p) );

				//register
				multiplexedOutFiles.add(f);
			}
			
		}
		
		/*
		 * Properly init output files for samples
		 * Check if file name was provided in barcode file else create default names
		 */
		log.debug("init output files for all samples...");
		sample2outfiles = new HashMap<String, List<File>>();
		if ( WRITE_ALL_READS_IN_SAME_OUTPUT ) {
			sample2outfiles.put(UNIQUE_MULTIPLEXED_SAMPLE_NAME, multiplexedOutFiles);
		} else {
			if(barcodeParser.getSample2outfiles() != null && barcodeParser.getSample2outfiles().size() != 0){
				log.debug("  file name or path provided in barcode file");
				Boolean filenamesAreAlreadyPathes = null;
				for (Entry<String, List<String>> e : barcodeParser.getSample2outfiles().entrySet()) {
					String sampleName = e.getKey();
					List<String> filenames = e.getValue(); // the parser already validated that we have a non empty value but that's all
					//init ?
					if(filenamesAreAlreadyPathes == null){
						filenamesAreAlreadyPathes = looksLikeAPath(filenames.get(0));
					}

					//convert to File
					List<File> files = new ArrayList<File>();
					for (String path : filenames) {
						//are we consistent ?
						if(filenamesAreAlreadyPathes != looksLikeAPath(path)){
							String m = null;
							if(filenamesAreAlreadyPathes)
								m = "User provided output file (from barcode) "+path+" is not a proper path while other user provided file names are fully qualified path. "; 
							else
								m = "User provided output file (from barcode) "+path+" looks like a path while other user provided file names are not. ";

							m+= "Please don't mix up and either use filename or file path for all samples.";
							return new String[]{m};
						}

						File f = null;
						if(filenamesAreAlreadyPathes)
							f = new File(path);
						else
							f = new File(this.OUTPUT_DIR, path);

						files.add(f);
					}
					sample2outfiles.put(sampleName, files);
				}
			}
			else{
				//no file name, init from output layouts
				log.debug("  file name or path NOT provided in barcode file, initializing from output layouts ...");
				for (Entry<String, List<Set<String>>> e : barcodeParser.getSample2BarcodeSets().entrySet()) {
					String sampleName = e.getKey();
					List<Set<String>> barcodeSetList = e.getValue();
					String key = generateSampleKeyForFilename(sampleName, barcodeSetList);
					List<File> files = new ArrayList<File>();
					for (int i = 0; i < outLayouts.length; i++) {
						FastqWriterLayout ol = outLayouts[i];
						String fname = key + 
								( ADD_SEQUENCE_LAYOUT_IN_OUTPUT_FILENAME ? "_" +ol.getReadSequenceLayout() : "" )+
								( ADD_HEADER_LAYOUT_IN_OUTPUT_FILENAME ? "_" +ol.getReadNameLayout() : "" )+
								( ADD_LAYOUT_IDX_IN_OUTPUT_FILENAME ? "_" + (i+1) :"" )+
								"." +(FASTQ_FILE_EXTENSION == null ? "txt" : FASTQ_FILE_EXTENSION) +
								(GZIP_OUTPUTS ? ".gz" : "") ;

						files.add( new File(OUTPUT_DIR , fname) );
					}
					sample2outfiles.put(sampleName, files);
				}
			}
		}
		log.debug("  output file names created, now checking if they already exist ");
		
		for (Entry<String, List<File>> en : sample2outfiles.entrySet()) {
			try{
				checkFilesDoNotAlreadyExist(en.getValue(), FORCE);
			}catch(Jexception e){
				return new String[]{e.getMessage()};
			}
		}
		log.debug("out file path fully validated.");
		
		
		/*
		 * init metrics and diag file
		 */
		log.debug("handling METRICS_FILE_NAME option ...");
		if(METRICS_FILE_NAME == null) 
			METRICS_FILE_NAME=DEFAULT_METRICS_FILE_NAME; //this is user mistake, reset the default name
		
		metricsFile = looksLikeAPath(METRICS_FILE_NAME) ? new File(METRICS_FILE_NAME) : new File(this.OUTPUT_DIR, METRICS_FILE_NAME) ;
		try{
			checkFilesDoNotAlreadyExist(metricsFile, FORCE);
		}catch(Jexception e){
			return new String[]{e.getMessage()};
		}
		log.debug("handling METRICS_FILE_NAME option ... OK");
		
		log.debug("handling BARCODE_DIAG_FILE option ...");
		if(BARCODE_DIAG_FILE != null){
			bcDiagFile = looksLikeAPath(BARCODE_DIAG_FILE) ? new File(BARCODE_DIAG_FILE) : new File(this.OUTPUT_DIR, BARCODE_DIAG_FILE) ;
			try{
				checkFilesDoNotAlreadyExist(bcDiagFile, FORCE);
			}catch(Jexception e){
				return new String[]{e.getMessage()};
			}
		}
		log.debug("handling BARCODE_DIAG_FILE option ... OK");	
		
		
		log.debug("Validation ended without error");
		return null; //do not return an empty array
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


	private String generateSampleKeyForFilename(
			String sampleName,
			List<Set<String>> barcodeSetList) {
		String k = sampleName.replaceAll("\\s+", "-");
		k = k.replaceAll(File.separator, "");
		String bcKey = "";
		for (Set<String> set : barcodeSetList) {
			bcKey += (bcKey.isEmpty() ? "" : "-");
			bcKey += set.iterator().next();
		}
		return k + "_" + bcKey;
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

	/** turn MM, MMD and Q options ie 'MM=2 or MM=2:3:4' into int arrays
	 * @param optionValue
	 * @param nBarcodes
	 * @return
	 * @throws NumberFormatException
	 * @throws Jexception
	 */
	private int[] parseIntegerOptions(String optionName, String optionValue, int nBarcodes) throws NumberFormatException, Jexception{
		int [] arr = new int[nBarcodes];

		if(optionValue.contains(":")){
			String [] tokens = optionValue.split(":");
			if(tokens.length != nBarcodes)
				throw new Jexception("Error while processing value of option "+optionName+" : "+tokens.length+" values found in "+optionValue+" while "+nBarcodes+" have been described.\n"+
						"You must provide either a unique value or the same number of values (separated by ':') as the number of described barcodes (BARCODE columns in barcode file).");
			for (int i = 0; i < tokens.length; i++) {
				try{
					arr[i] = Integer.parseInt(tokens[i]);
				}catch(NumberFormatException e ){
					throw new Jexception("Error while processing the "+optionName+" option '"+optionValue+"': '"+tokens[i]+"' is not a valid number");
				}
			}
		}else{
			try{
				Arrays.fill(arr, Integer.parseInt(optionValue));
			}catch(NumberFormatException e ){
				throw new Jexception("Error while processing the value of option "+optionName+" : '"+optionValue+"' is not a valid number");
			}
		}

		return arr;
	}
	
	/**
	 * parse the string 'TRUE' (case insensitive) or '1' into true ; other string values become false
	 * 
	 * @param optionName
	 * @param optionValue
	 * @param nBarcodes
	 * @return
	 * @throws Jexception
	 */
	private boolean[] parseBooleanOptions(String optionName, String optionValue, int nBarcodes) throws Jexception{
		boolean [] arr = new boolean[nBarcodes];

		if(optionValue.contains(":")){
			String [] tokens = optionValue.split(":");
			if(tokens.length != nBarcodes)
				throw new Jexception("Error while processing value of option "+optionName+" : "+tokens.length+" values found in "+optionValue+" while "+nBarcodes+" have been described.\n"+
						"You must provide either a unique value or the same number of values (separated by ':') as the number of described barcodes (BARCODE columns in barcode file).");
			for (int i = 0; i < tokens.length; i++) {
				arr[i] = Boolean.parseBoolean(tokens[i]); // give true only if token is == TRUE (case insensitive)
				if(!arr[i])
					arr[i] = tokens[i].equals("1");
			}
		}else{
			boolean b = Boolean.parseBoolean(optionValue);
			if(!b)
				b = optionValue.equals("1");
			Arrays.fill(arr, b);
		}

		return arr;
	}




	/**
	 * Fetches library names and barcodes from embase (using the fastq file path) 
	 * Export a barcode file with demultiplexed file paths automatically set to their expected 
	 * definitive locations. Also save metrics and unassigned files to their expected 
	 * definitive locations in embase NGS data structure. 
	 * Also enforces :
	 *  KEEP_UNASSIGNED_READ = true, 
	 *  GZIP_OUTPUTS=true,
	 *  CREATE_MD5_FILE = true;
	 *  
	 *  Note that the code first checks that user has sufficient permissions on item (write 
	 *  permissions is expected on NGS Assay)
	 *  
	 * @return
	 */
	private List<String> customEmbaseCommandLineValidation() {
		//error messages to be accumulated
		List<String> messages = new ArrayList<String>();
			
		OUTPUT_DIR = null;
		GZIP_OUTPUTS=true;
		CREATE_MD5_FILE = true;
		
		//we keep unassigned reads
		KEEP_UNASSIGNED_READ = true;
		//and call files according to input files, simply adding '_unassigned-reads' before extension
		for(File f : FASTQ){
			String ufname= FileUtil.removeExtension( f.getName() ) + "_unassigned-reads"+".txt.gz" ;
			UNASSIGNED_FILE.add(ufname);
		}
				
		/*
		 * Get information from emBASE 
		 */		
		
		//try load properties
		try {
			ApplicationConfiguration.init();
		} catch (IOException e) {
			log.error("Failed to read properties from property file.");
			throw new RuntimeException(e);
		}
		//Determine the caller
		String username = System.getProperty("user.name");
		//check user is not trying to trick me with a -Duser.name=another user!
		String whoami = JeUtils.whoami();
		if(!username.equalsIgnoreCase(whoami)){
			messages.add("User "+username+" does not match 'whoami' ("+whoami+") ! Are you trying to hack the system ?");
			//we cannot continue
			return messages;
		}
		
		EmBaseDatabaseFacade facade = null;
		try {
			facade = EmBaseDatabaseFacade.createEmBaseDatabaseFacade(
					username,
					ApplicationConfiguration.DB_URL,
					ApplicationConfiguration.DB_USER,
					ApplicationConfiguration.DB_PWD
					);
		} catch (EmBASEConnectionException e) {
			log.error(ExceptionUtil.getStackTrace(e));
			messages.add("Error when trying to connect emBASE DB : "+e.getMessage());
			//we cannot continue
			return messages;
		}
		
		//what s the lane number ? first get teh NGS assay
		Item ngsAssay = facade.getNGSAssay(FASTQ.get(0));
		List<CSVLine> rbaDetails = null;
		if(ngsAssay == null){
			messages.add("No NGS Assay found in emBASE for lane file : "+FASTQ.get(0).getAbsolutePath()+
					"\nAre you sure this file has been already imported in emBASE ? If so, please contact emBASE admins.");
		}else{
			//then fetch details
			ArrayList<Item> _list = new ArrayList<Item>();
			_list.add(ngsAssay);
			/*
			 * rbaDetails contains values for all org.embl.gbcs.embase.api.Queries.HEADER_*
			 */
			rbaDetails = facade.fetchRBADetailsForNGSAssays(_list).values().iterator().next();
			if(rbaDetails == null || rbaDetails.size()==0){
				messages.add("No RawBioAssay found in emBASE for NGS Assay "+ngsAssay.getName()+" (id="+ngsAssay.getId()+") "
						+ "i.e. lane file : "+FASTQ.get(0).getAbsolutePath()+
						"\nAre you sure this file has been successfully imported in emBASE ? If so, please contact emBASE admins.");
			}
		}
		
		if(messages.size()!=0){
			//we stop here
			return messages;
		}
		
		//get lane number : any RBA can be used for this
		Integer laneNumber = null;
		try {
			laneNumber = rbaDetails.get(0).getValueAsInt(Queries.HEADER_Lane);
		} catch (InvalidHeaderException e1) {
			//not possible if the code has been run at least once!
			throw new RuntimeException(e1);
		}
		
		// Check that expected directories exist
		File runDir = FASTQ.get(0).getParentFile();
		//get lane dir
		File laneDir = new File (runDir, EmBaseDatabaseFacade.getLaneDirName(laneNumber));
		//then check if expected dir exists for stats and unassigned reads
		File runStatDir = new File (laneDir, EmBaseDatabaseFacade.STAT_DIR_NAME);
		
		if(!runStatDir.exists() || !runStatDir.canWrite()){
			if(!TEST_MODE_STOP_AFTER_PARSING)
				messages.add("The expected '"+EmBaseDatabaseFacade.STAT_DIR_NAME+"' dir in the lane dir '"+laneDir.getAbsolutePath()+
						"' is either not existing or not writable");
			else{
				METRICS_FILE_NAME = new File(EmBaseDatabaseFacade.TMP_DIR, FileUtil.removeExtension( FASTQ.get(0).getName() ) +"_jemultiplexing_metrics.txt").getAbsolutePath();
				log.info("TEST MODE: METRICS_FILE_NAME set to "+METRICS_FILE_NAME);
			}
		}else{
			//set the full path so the remaining of the code leaves it like this
			METRICS_FILE_NAME = new File(runStatDir, FileUtil.removeExtension( FASTQ.get(0).getName() ) +"_jemultiplexing_metrics.txt").getAbsolutePath();
			log.info("METRICS_FILE_NAME set to "+METRICS_FILE_NAME);
		}
		
		File unassignedReadsDir = new File (laneDir, EmBaseDatabaseFacade.UNASSIGNED_DIR_NAME);
		if(!unassignedReadsDir.exists() || !unassignedReadsDir.canWrite()){
			if(!TEST_MODE_STOP_AFTER_PARSING)
				messages.add("The expected '"+EmBaseDatabaseFacade.UNASSIGNED_DIR_NAME+"' dir in the run dir '"+laneDir.getAbsolutePath()+
						"' is either not existing or not writable");
			else{
				OUTPUT_DIR = EmBaseDatabaseFacade.TMP_DIR;
				log.info("TEST MODE: OUTPUT_DIR set to "+OUTPUT_DIR.getAbsolutePath());
			}
		}else{
			 //set output dir to where the input files are located => this will be the place for unassigned reads 
			OUTPUT_DIR = unassignedReadsDir;
			log.info("OUTPUT_DIR set to "+OUTPUT_DIR.getAbsolutePath());
		}

		
		/*
		 * Check caller can actually demultiplex this dataset
		 * The user should have read access on the NGS assay
		 */
		try {
			if(!facade.canDemultiplex(FASTQ.get(0), username)){
				messages.add("According to emBASE DB, the user "+username+" is not allowed to read NGS Assay corresponding to file "+
						FASTQ.get(0).getAbsolutePath());
			}else{
				log.debug("User allowed to demultiplex");
			}
		}catch (Exception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			messages.add("Error when trying to check read rights from emBASE DB : "+e.getMessage());
		}
		
		
		/*
		 * Get all libs 
		 */
		List<NGSLibrary> libraries = null;
		try {
			libraries = facade.fetchLibraries(FASTQ.get(0));
			facade.findSampleStorageDir(laneDir, libraries);
			
			//and check dirs have been set for all libs
			for (NGSLibrary lib : libraries) {
				if(lib.getSampleDir() == null || !lib.getSampleDir().exists()){
					String mess = "No library dir for sample/library "+lib.getName()+" is available.\nThe emBASE storage architecture is not ready. " +
							"Are you sure this lane has been successfully uploaded to emBASE ?\nIf so, please contact emBASE admins.";
					log.debug(mess);
					if(TEST_MODE_STOP_AFTER_PARSING){
						log.info("TEST MODE : Setting fake dir for "+lib.getName());
						lib.setSampleDir ( EmBaseDatabaseFacade.getFakeSampleDir(lib.getName(), lib.getId()) );//does not exist, just fake for test
					}
					else{
						messages.add(mess);
					}
				}
				else {
					//update lib.getSampleDir with fastq sub-dir
					File fastqDir = new File(lib.getSampleDir(), "fastq");
					if(!fastqDir.exists() || !fastqDir.canWrite()){
						if(!TEST_MODE_STOP_AFTER_PARSING){
							String mess = "The 'fastq' dir under the library dir  "+lib.getSampleDir().getAbsolutePath()
									+" is either not existing OR not writable.\n";
							if(!fastqDir.exists()){
								mess+="The storage architecture is not ready. Are you sure this lane has been successfully uploaded to emBASE ?\nIf so, please contact emBASE admins.";
							}
							else{
								mess+="The fastq storage file has been already locked to prevent further modifications.\n" +
										"Please contact emBASE admins if you think this should not be the case.";
							}
							log.debug(mess);
							messages.add(mess);
						} else {
							lib.setSampleDir( EmBaseDatabaseFacade.TMP_DIR );
							log.info("TEST MODE : Setting fake sample dir to "+EmBaseDatabaseFacade.TMP_DIR);
						}
					}else{
						lib.setSampleDir( fastqDir );
					}
				}
			}
		} catch (Exception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			messages.add("Error when trying to fetch barcodes from emBASE DB : "+e.getMessage());
		}
		
		if(messages.size()!=0){
			//we stop here
			return messages;
		}
		
		/*
		 * All seems perfect, create barcode file 
		 * Check demultiplexed files are not already found
		 */
		PrintWriter pw = null;
		try{
			BARCODE_FILE = getBarcodeFileTmpLocation(
					FASTQ.get(0), 
					username, 
					(TEST_MODE_STOP_AFTER_PARSING ? EmBaseDatabaseFacade.TMP_DIR : runStatDir)
					);
			if(BARCODE_FILE.exists()){
				log.warn("The barcode file already exists and will be overwritten : "+BARCODE_FILE.getAbsolutePath());
			}
			pw = new PrintWriter(BARCODE_FILE);
			for (NGSLibrary lib : libraries) {
				
				File fastq1 = getDemultiplexedFile(lib, 1 , GZIP_OUTPUTS);
				File fastq2 = getDemultiplexedFile(lib, 2 , GZIP_OUTPUTS);
				
				/*
				 * check if file already exists , the test on the length is meant to avoid considering a 'fake' file as a real file
				 * we indeed create such files to prepare the architecture
				 */
				if(fastq1.exists() && fastq1.length() > 10*1000L && !FORCE){ 
					if(fastq1.canWrite()){
						messages.add("Demultiplexed file "+fastq1.getAbsolutePath()+" already exists. You need to run with FORCE=true to overwrite existing files.");
					}else{
						messages.add("Demultiplexed file "+fastq1.getAbsolutePath()+" already exists BUT you do NOT have write permission on it.\n" +
								"You first need to solve the write permission issue, then use FORCE=true to overwrite existing files.");
					}
				}
				else if(fastq1.exists() && FORCE && fastq1.canWrite()){
					log.warn("Existing demultiplexed file will be overritten (option FORCE=true) : "+fastq1.getAbsolutePath());
				}
				
				if(fastq2.exists() & fastq2.length() > 10*1000L && !FORCE){
					if(fastq2.canWrite()){
						messages.add("Demultiplexed file "+fastq2.getAbsolutePath()+" already exists. You need to run with FORCE=true to overwrite existing files.");
					}else{
						messages.add("Demultiplexed file "+fastq2.getAbsolutePath()+" already exists BUT you do NOT have write permission on it.\n" +
								"You first need to solve the write permission issue, then use FORCE=true to overwrite existing files.");
					}
				}else if( fastq2.exists() && FORCE && fastq2.canWrite()){
					log.warn("Existing demultiplexed file will be overritten (option FORCE=true) : "+fastq2.getAbsolutePath());
				}
				
				pw.print(lib.getName()+"\t"+lib.getBarcode()+"\t"+fastq1.getAbsolutePath());
				if(FASTQ.size() > 1)
					pw.print("\t"+fastq2.getAbsolutePath());
				pw.println();
			}
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}finally{
			if(pw!=null)
				pw.close();
		}
		
		//then ok just let the code go on
		
		return messages;
	
	}



	@Override
	protected int doWork() {
		
		log.debug("Launching demultiplexing...");
		
		try{
			Demultiplexer d = new Demultiplexer(
					this,
					FASTQ.toArray(new File[FASTQ.size()]), 
					readLayouts, 
					barcodeParser.getSample2BarcodeSets(), 
					sample2outfiles, 
					outLayouts, 
					this.unassignedFiles, 
					this.metricsFile, 
					this.QUALITY_FORMAT,
					this.GZIP_OUTPUTS, 
					this.CREATE_MD5_FILE);
			
			if(INPUT_FASTQ_COMPRESSION != null) {
				d.setCompressedInputFastqFiles(INPUT_FASTQ_COMPRESSION);
			}
			
			log.debug("Demultiplexer successfully initialized... now running...");
			
			d.run(
					max_mismatches, 
					min_mismatch_deltas, 
					min_base_qualities,
					STRICT,
					WRITER_FACTORY_USE_ASYNC_IO, 
					bcDiagFile);
			
			log.debug("Demultiplexer run ended without warning.");
			
		}catch(Exception e){
			log.error(ExceptionUtil.getStackTrace(e));
			log.error("\n\n\n\n\nAn error occurred during demultiplexing, please check the error message before running the process again. \n"+
			         "Error message was :\n"+e.getMessage());
			return 1;
		}
		return 0;
	}
	
	
	/**
	 * Print the final demultiplexing report
	 *   
	 * @param total_read_count 
	 * @param unassigned number of read (pair) unassigned
	 * @param assigned  number of read (pair) assigned
	 */
	protected void printMetricFile(Map<String, Integer> sampleCountMap, int total_read_count, int unassigned, int assigned) {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(this.metricsFile);
			
				
			pw.println("## Created on " + new Date() + " by je "+Je.COMMAND_DEMULTIPLEX +" version "+getVersion());
			pw.println("## Command line : "+this.getCommandLine());           
           
			pw.println("Processed Reads (pairs)"+"\t"+ total_read_count);
			pw.println("Assigned Reads (pairs)"+"\t"+ assigned );
			pw.println("Unassigned Reads (pairs)"+"\t"+ unassigned);
			pw.println("# Individual sample read (pair) counts :");
			for (Entry<String, Integer> e : sampleCountMap.entrySet()) {
				pw.println(e.getKey()+"\t"+ e.getValue() );
			}

		} catch (Exception e) {
			log.error(ExceptionUtil.getStackTrace(e));
			throw new RuntimeException(e);
		} finally {
			if(pw!=null)
				pw.close();
		}
	}
	
	/**
	 * builds a file to use for barcode file
	 * 
	 * @param fastqFile the fastq file (mate 1 in case of PE) 
	 * @param username the user how is about to demultiplex
	 * @param d the dir to use or null to use system defaults
	 * @return
	 */
	protected File getBarcodeFileTmpLocation(File fastqFile, String username, File d) {
		String tmpdir = (d == null ? System.getProperty("java.io.tmpdir") : d.getAbsolutePath() );
		if(!StringUtil.isValid(tmpdir) || !new File(tmpdir).exists())
			tmpdir = ".";
		
		File f= new File(
				tmpdir,
				FileUtil.removeExtension( fastqFile.getName() ) +"_"+username+"_barcodes.txt"
				);
		log.debug("barcode file to use is "+f.getAbsolutePath());
		return f;
	}
	
	/**
	 * 
	 * Find out what file should be used for a given library. This method does NOT fetches existing file
	 * but builds a file path to use by convention i.e. when you want to automatically demultiplex fastq file
	 *  
	 * 
	 * @param lib the lib for which we need to know the file to use to write demultiplexed reads.
	 * Note that lib should have its final sampleDir correctly set 
	 * @param readPairNumber 1 or 2 (for single end, it is always 1) 
	 * @param gzip whether output should be gzipped
	 * @return the File to use
	 */
	protected File getDemultiplexedFile(NGSLibrary lib, int readPairNumber, boolean gzip) {
		
		//for now simply build up a name with lib name
		File f =  new File(lib.getSampleDir(), lib.getName()+"_"+lib.getSampleDir().getParentFile().getName()+"_"+readPairNumber+".txt"+(gzip?".gz":"")); 
		log.debug("Demultiplexed file for "+lib.getName()+" (read pair "+readPairNumber+") should be "+f.getAbsolutePath());
		return f;
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
	
	
	/**
	 * 
	 * Init default layouts according to situation :
	 *  
	 * ########### READ LAYOUTS  
	 * 1 input fastq: single end with layout <BARCODE:X><SAMPLE:x> where X is inferred from barcode file
	 * 2 input fastqs: 
	 *       - paired end with redundant barcode if barcode file describes a single BARCODE i.e. <BARCODE1:X><SAMPLE1:x> and <BARCODE1:X><SAMPLE2:x>, where X is inferred from barcode file
	 *       - paired end with non redundant barcode if barcode file describes two BARCODE column i.e. <BARCODE1:X><SAMPLE1:x> and <BARCODE2:Y><SAMPLE2:x>, where X and Y are inferred from barcode file
	 *       - single end with index file if barcode file describes a single BARCODE and second fastq file has reads of length < 10 + barcode_length
	 * 3 input fastqs: 
	 *       - paired end with an index file i.e. <SAMPLE1:x>, <SAMPLE2:x> and <BARCODE1:X> when barcode file has a single BARCODE column (X is inferred from barcode file)
	 *       - single end with two index files i.e. <SAMPLE1:x>, <BARCODE1:X> and <BARCODE2:Y> when barcode file has two BARCODE columns (X,Y is inferred from barcode file)
	 * 4 input fastqs: paired end with either 
	 *       - 2 non-redundant index files i.e. <SAMPLE1:x>, <SAMPLE2:x>, <BARCODE1:X>, <BARCODE2:Y> if the barcode file has two BARCODE columns or a ATGC:GCTAGC syntax (X,Y inferred from barcode file) 
	 *       - 2 redundant index files <SAMPLE1:x>, <SAMPLE2:x>, <BARCODE1:X> and <BARCODE1:X> if the barcode file has a single BARCODE column (X inferred from barcode file)
	 * 
	 * 
	 * ########### OUTPUT LAYOUTS
	 * 
	 * Describes the output file layout(s) using the slots defined in read layouts and ':' to delineate three blocks e.g. 'OL=1:<BARCODE1><UMI1><UMI2>:<SAMPLE1>' : <br />
	 * <ol>
	 * <li> The first block is a integer number (i.e. from '1:' above) is the output file index and it must be unique across all 'OL' inputs. 
	 *    Inferred from order in comamnd line when not given </li>
	 * <li> The second part (i.e. '<BARCODE1><UMI1><UMI2>' above) is the read header layout ; when writing multiple UMI and BARCODE slots 
	 *   in output read headers, these are always separated with the RCHAR (':' by defaults). </li>
	 * <li> The third part (i.e. '<SAMPLE1>' above) is the read sequence layout.</li>
	 * </ol>
	 * One output file is created for each sample and each OL index. Output file names default to samplename_outputfileindex with the original extensions. <br />
	 * 
	 * When no OL is described, Je considers an output file should be created for each input FASTQ (containing a SAMPLE slot) and for each sample. <br />
	 * In this scenario:
	 * <ol> 
	 * <li> 1. The output files only contain the SAMPLE slot unless CLIP is set to false </li>
	 * <li> 2. The barcode(s) and sample names are injected in the output file names according to the pattern 'FASTQFILENAMEn_SAMPLENAME_BARCODES.ORIGINALEXTENSIONS' ) </li>
	 * <li> 3. All BARCODE and UMI slots (if any) are placed in the fastq headers following their slot index i.e. BARCODE1:...:BARCODEn:UMI1:UMI2:...:UMIn and are separated 
	 * with ':' (all UMIs are added to all file sets identically) ; unless ADD is et to false.</li>
	 * 
	 * 
	 * @param fastqFiles
	 * @param barcodeParser fully initialized parser on the validated barcode file
	 * @return
	 */
	private class DefaultLayouts{
		ReadLayout[] readLayouts;
		FastqWriterLayout[] fasqWriterLayouts;
		
		public DefaultLayouts(
				List<File> fastqFiles,
				BarcodeFileGeneralParser barcodeParser,
				boolean clip, 
				boolean add,
				boolean withQualityInReadName, 
				String readnameDelimitor, 
				FastqQualityFormat fastqQualityFormat
				) {
			
			int bcNum = barcodeParser.getBarcodeLengths().size();
			if(bcNum > 2)
				throw new Jexception("Je can't guess the layout of your FASTQ files with more thantwo BARCODE columns in barcode file. Please provide read layouts");
			
			readLayouts = new ReadLayout[fastqFiles.size()];
			fasqWriterLayouts = new FastqWriterLayout[fastqFiles.size()];
			
			switch (fastqFiles.size()) {
			case 1:
				//
				if(bcNum!=1)
					throw new Jexception("Je can't guess the layout of your FASTQ file: more than one BARCODE column found in barcode file while no read layout were provided on the command line.");
				int bcLen = barcodeParser.getBarcodeLengths().get(0);
				readLayouts[0] = new ReadLayout("<BARCODE1:"+bcLen+"><SAMPLE1:x>");
				fasqWriterLayouts[0] = new FastqWriterLayout( (clip ? "" : "B1") + "S1", (add ? "B1" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
				break;
			case 2:
				if(bcNum == 2){
					//clear case => PE with non redundant barcodes
					readLayouts[0] = new ReadLayout("<BARCODE1:"+barcodeParser.getBarcodeLengths().get(0)+"><SAMPLE1:x>");
					readLayouts[1] = new ReadLayout("<BARCODE2:"+barcodeParser.getBarcodeLengths().get(1)+"><SAMPLE2:x>");
				}else{
					bcLen = barcodeParser.getBarcodeLengths().get(0);
					if (peekReadLength(fastqFiles.get(1)) < (10 + bcLen)){
						//SE with index file
						readLayouts[0] = new ReadLayout("<SAMPLE1:x>");
						readLayouts[1] = new ReadLayout("<BARCODE1:"+bcLen+">");
						fasqWriterLayouts[0] = new FastqWriterLayout( "S1", (add ? "B1" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
						fasqWriterLayouts[1] = new FastqWriterLayout( "B1", (add ? "B1" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
					}else{
						//PE with redundant barcode
						readLayouts[0] = new ReadLayout("<BARCODE1:"+bcLen+"><SAMPLE1:x>");
						readLayouts[1] = new ReadLayout("<BARCODE1:"+bcLen+"><SAMPLE2:x>");
						fasqWriterLayouts[0] = new FastqWriterLayout( (clip ? "" : "B1") + "S1", (add ? "B1" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
						fasqWriterLayouts[1] = new FastqWriterLayout( (clip ? "" : "B1") + "S2", (add ? "B1" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
					}
				}
				break;
			case 3:
				if(bcNum == 1){
					//clear case => PE with redundant barcodes from index file
					readLayouts[0] = new ReadLayout("<SAMPLE1:x>");
					readLayouts[1] = new ReadLayout("<SAMPLE2:x>");
					readLayouts[2] = new ReadLayout("<BARCODE1:"+barcodeParser.getBarcodeLengths().get(0)+">");
					fasqWriterLayouts[0] = new FastqWriterLayout( "S1", (add ? "B1" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
					fasqWriterLayouts[1] = new FastqWriterLayout( "S2", (add ? "B1" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
					fasqWriterLayouts[2] = new FastqWriterLayout( "B1", (add ? "B1" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
				}else{
					//SE with 2 non redundant barcodes from index files
					readLayouts[0] = new ReadLayout("<SAMPLE1:x>");
					readLayouts[1] = new ReadLayout("<BARCODE1:"+barcodeParser.getBarcodeLengths().get(0)+">");
					readLayouts[2] = new ReadLayout("<BARCODE2:"+barcodeParser.getBarcodeLengths().get(1)+">");
					fasqWriterLayouts[0] = new FastqWriterLayout( "S1", (add ? "B1B2" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
					fasqWriterLayouts[1] = new FastqWriterLayout( "B1", (add ? "B1B2" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
					fasqWriterLayouts[2] = new FastqWriterLayout( "B2", (add ? "B1B2" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
				}
				break;
			case 4:
				if(bcNum == 1){
					// PE with 2 redundant barcodes from index file
					readLayouts[0] = new ReadLayout("<SAMPLE1:x>");
					readLayouts[1] = new ReadLayout("<SAMPLE2:x>");
					readLayouts[2] = new ReadLayout("<BARCODE1:"+barcodeParser.getBarcodeLengths().get(0)+">");
					readLayouts[3] = new ReadLayout("<BARCODE1:"+barcodeParser.getBarcodeLengths().get(0)+">");
					fasqWriterLayouts[0] = new FastqWriterLayout( "S1", (add ? "B1" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
					fasqWriterLayouts[1] = new FastqWriterLayout( "S2", (add ? "B1" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
					fasqWriterLayouts[2] = new FastqWriterLayout( "B1", (add ? "B1" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
				}else{
					// PE with 2 non redundant barcodes from index file
					readLayouts[0] = new ReadLayout("<SAMPLE1:x>");
					readLayouts[1] = new ReadLayout("<SAMPLE2:x>");
					readLayouts[2] = new ReadLayout("<BARCODE1:"+barcodeParser.getBarcodeLengths().get(0)+">");
					readLayouts[3] = new ReadLayout("<BARCODE2:"+barcodeParser.getBarcodeLengths().get(1)+">");
					fasqWriterLayouts[0] = new FastqWriterLayout( "S1", (add ? "B1B2" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
					fasqWriterLayouts[1] = new FastqWriterLayout( "S2", (add ? "B1B2" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
					fasqWriterLayouts[2] = new FastqWriterLayout( "B1", (add ? "B1B2" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
					fasqWriterLayouts[3] = new FastqWriterLayout( "B2", (add ? "B1B2" : ""), readLayouts, withQualityInReadName, readnameDelimitor, fastqQualityFormat);
				}
				break;
			}
			

			for (FastqWriterLayout l : fasqWriterLayouts) {
				l.setWithQualityInReadName(withQualityInReadName);
			}
			
		}

		
		/**
		 * @return the readLayouts
		 */
		public ReadLayout[] getReadLayouts() {
			return readLayouts;
		}

		/**
		 * @return the fasqWriterLayouts
		 */
		public FastqWriterLayout[] getFasqWriterLayouts() {
			return fasqWriterLayouts;
		}
		
	}
	
}
