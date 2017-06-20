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
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;

public class Demultiplexer extends CommandLineProgram {
	private static Logger log = LoggerFactory.getLogger(Demultiplexer.class);

	/*
	 * We define here all common options
	 */
	
	@Option(shortName="F", optional = false,
			printOrder=10,
			doc="Input fastq file (optionally gzipped) with optional unique number and read layout " +
			"i.e. 'FASTQ=1:<BARCODE1:6><SAMPLE:x>:file.fastq.gz' ; with the ':' character as the delimiter." +
			"The first number (e.g.1:) is mandatory when more than one FASTQ input is provided (e.g. paired end situation). " +
			"When present it must be unique across all 'FASTQ' inputs. This number is used to uniquely refer to this input fastq" +
			" file and match up with the output files (see below)"
			)
	public List<String> FASTQ;

	
	/*
	 * Barcode definition options
	 */
	@Option(shortName="BF", optional = false,  mutex={"USE_EMBASE"},
			printOrder=30,
			doc="Barcode files with optional BARCODE slot matching i.e. "
					+ "BF=<BARCODE1>:barcodes.txt or BF=<AND:BARCODE1,BARCODE2>:barcodes.txt or BF=<OR:BARCODE1,BARCODE2>:barcodes.txt (arbitrary number of BARCODE can be listed as csv list) "
					+ "The BARCODE slot matching is mandatory as soon as more than one BARCODE slot in defined across all FASTQ input. "
					+ "Barcode file describing sequence list and sample names. "
					+"Tab-delimited file with 2 columns, with the sample in col1 and the corresponding barcode in col2.\n" +
					"Simple barcode file format : 2 tab-delimited colums\n"+
					"\t"+"If multiple barcode map to the same sample, either line can be duplicated e.g.\n" +
					"\t"+"\t"+"sample1\tATAT\n" +
					"\t"+"\t"+"sample1\tGAGG\n" +
					"\t"+"\t"+"sample2\tCCAA\n" +
					"\t"+"\t"+"sample2\tTGTG\n" +
					"\t"+"Or barcodes can be combined using the OR operator '|' i.e. the file above can be re-written like\n " +
					"\t"+"\t"+"sample1\tATAT|GAGG\n" +
					"\t"+"\t"+"sample2\tCCAA|TGTG\n" +
					"\t"+"When different barcodes must be combined for the look up (e.g. paired-end data with different barcodes in fwd and rev reads), " +
					"barcodes must be separated using a ':' separator i.e. \n" +
					"\t"+"\t"+"sample1\tATAT:GAGGTG\n" +
					"\t"+"\t"+"sample2\tCCAA:TGTGCA\n" +
					"\t"+"This above syntax means that sample 1 is encoded with ATAT barcode (coming from a BARCODE slot e.g." +
					" in fwd read) AND a GAGGTG barcode (coming from a different BARCODE slot e.g. in rev read). " +
					"Note that you can still combine barcodes using | e.g. \n"+
					"\t"+"\t"+"sample1\tATAT|GAGGTG:CCAA|TGTGCA\n" +
					"Extended barcode file format : 3 (single-end) or 4 (paired-end) tab-delimited colums\n"+
					"\t"+"same as the simple barcode file format but the extra columns contains the file name(s) to use to name output files." +
					" A unique extra column is expected for single-end while 2 extra columns are expected for paired-end. In case, lines are duplicated (multiple barcodes" +
					"mapping the same sample), the same file name should be indicated in the third (and fourth) column(s). \n"+
					"\t"+"\t"+"sample1\tATAT\tspl1_1.txt.gz\tspl1_2.txt.gz\n" +
					"\t"+"\t"+"sample1\tGAGG\tspl1_1.txt.gz\tspl1_2.txt.gz\n" +
					"\t"+"\t"+"sample2\tCCAA\tspl2_1.txt.gz\tspl2_2.txt.gz\n" +
					"\t"+"Or\n" +
					"\t"+"\t"+"sample1 \t ATAT|GAGG:CCAA|TGTG \t spl1_1.txt.gz \t spl1_2.txt.gz\n" +
					"Ns in barcode sequence are allowed and are used to flag positions that should be ignored in sample matching \n"+
					"\t i.e. they will be clipped off the read sequence (like in iCLIP protocol)."
			)
	public List<String> BARCODE_FILE = null;
	

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
	
	


	@Override
	protected int doWork() {
		// TODO Auto-generated method stub
		return 0;
	}

}
