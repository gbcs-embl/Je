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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.embl.cg.utilitytools.utils.ExceptionUtil;
import org.embl.cg.utilitytools.utils.FileUtil;
import org.embl.cg.utilitytools.utils.StringUtil;
import org.embl.cg.utilitytools.utils.parser.csv.CSVLine;
import org.embl.cg.utilitytools.utils.parser.csv.CSVParser;
import org.embl.cg.utilitytools.utils.parser.csv.InvalidCSVSetUpException;
import org.embl.gbcs.je.Jexception;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;




/**
 * 
 * Validator for barcode file
 * 
 * Barcode file (tsv) matching sample names to barcode combination. 
 * 
 * SIMPLE Barcode File Format (for backward compatibility) 
 * 
 * The file must have 2 columns with the sample in col1 and the corresponding barcode in col2.
 * In this format, a simple BARCODE slot is expected in the ReadLayout and NO headers are needed e.g. :<br />
 * <pre>
 *    sample1\tGAGG
 *    sample2\tCCAA
 * </pre>
 * <br/>
 * The format accept the following shortcuts:<br />
 * <ol> 
 * <li>If multiple barcodes map to the same sample, either line can be duplicated e.g. <br />
 * <pre>
 *    sample1\tATAT
 *    sample1\tGAGG
 *    sample2\tCCAA
 *    sample2\tTGTG
 * </pre>
 * <br />
 * 
 * Or barcodes can be combined using the OR operator '|' i.e. the file above can be re-written like 
 * <pre>
 *    sample1\tATAT|GAGG
 *    sample2\tCCAA|TGTG
 * </pre>
 * </li>
 * <li> For the special situation of paired-end data in which barcodes differ at both ends i.e. with 
 * BARCODE1 and BARCODE2 described for read one and two respectively, barcodes for BARCOD1 and BARCODE2 can be
 *  distinguished using a ':' separator i.e. <br />
 * <pre>
 *    sample1\tATAT:GAGG
 *    sample2\tCCAA:TGTG
 * </pre>
 * <br />
 * This above syntax means that sample 1 is encoded with ATAT barcode from BARCODE1 slot AND GAGG barcode from BARCODE2 slot. 
 * Note that you can still combine barcodes using | e.g. <br />
 * <pre> 
 *    sample1\tATAT|GAGG:CCAA|TGTG
 * </pre>
 * </li>
 * <li> Extended barcode file format has 3 (single-end) or 4 (paired-end) tab-delimited columns i.e. like 
 * in the simple barcode file format but the extra columns contains the file name(s) to use to name output files.
 * A unique extra column is expected for single-end while 2 extra columns are expected for paired-end. 
 * In case, lines are duplicated (multiple barcodes mapping the same sample), the same file name should 
 * be indicated in the third (and fourth) column(s). <br />
 * <pre> 
 *    sample1\tATAT\tspl1_1.txt.gz\tspl1_2.txt.gz
 *    sample1\tGAGG\tspl1_1.txt.gz\tspl1_2.txt.gz
 *    sample2\tCCAA\tspl2_1.txt.gz\tspl2_2.txt.gz
 * </pre> <br/>
 * Or<br/>
 * <pre>
 *    sample1 \t ATAT|GAGG:CCAA|TGTG \t spl1_1.txt.gz \t spl1_2.txt.gz
 * </pre>
 * </li>
 * 
 * 
 * @author girardot
 *
 */
public class BarcodeFileSimpleToGeneralFormatConverter {

	private static Logger log = LoggerFactory.getLogger(BarcodeFileSimpleToGeneralFormatConverter.class);
	
	protected File file ;
	protected Integer columnNumber;
	protected Boolean splitBarcodesInTwoColumn ;
	
	public BarcodeFileSimpleToGeneralFormatConverter(File f){
		this.file = f;
	}
	
	public void convertToGeneralFormat(File convertedFile){
		
		if(convertedFile.exists() && !convertedFile.delete())
			throw new Jexception("output file "+convertedFile.getAbsolutePath()+" already exist, please delete it first.");
		
		PrintWriter pw = null;
		try {
			CSVParser p = null;
			/*
			 * first pass to init needed values
			 */
			try {
				p = new CSVParser(2,4,false);
				Iterator<CSVLine> lines = p.iterator(file.getAbsolutePath(), "\t");
				
				while (lines.hasNext()) {
					CSVLine l = lines.next();
					if(columnNumber == null){
						//it is the first line then ; we need to check the column number
						columnNumber = l.getColNum();
						continue; //make sure to ignore the first valid line in case there are headers
					}
					//second line, check if barcodes contain :
					splitBarcodesInTwoColumn = l.getValue(1).contains(":");
					break;
				}
			} finally{
				if(p!=null)
					p.terminate();
			}
			
			/*
			 * second pass to convert
			 */
			
			//writer
			pw = new PrintWriter(convertedFile);
			//print headers
			ArrayList<String> headers = new ArrayList<String>();
			headers.add(BarcodeFileGeneralParser.HEADER_SAMPLE);
			headers.add(BarcodeFileGeneralParser.HEADER_BARCODE+"1");
			if(splitBarcodesInTwoColumn)
				headers.add(BarcodeFileGeneralParser.HEADER_BARCODE+"2");
			if(columnNumber >= 3){
				headers.add(BarcodeFileGeneralParser.HEADER_OUT+"1");
			}
			if(columnNumber >= 4){
				headers.add(BarcodeFileGeneralParser.HEADER_OUT+"2");
			}
			pw.println(StringUtil.mergeList(headers, "\t"));
			
			//now write the lines
			try {
				p = new CSVParser(2,4,false);
				Iterator<CSVLine> lines = p.iterator(file.getAbsolutePath(), "\t");
				
				while (lines.hasNext()) {
					CSVLine l = lines.next();
					pw.print(l.getValue(0));
					if(splitBarcodesInTwoColumn){
						String [] twoValues = l.getValue(1).split(":");
						if(twoValues.length != 2 ){
							throw new Jexception("Expecting two barcodes separated with a ':' but got "+twoValues.length+" from barcode expression "+l.getValue(1));
						}
						pw.print("\t" + twoValues[0]);
						pw.print("\t" + twoValues[1]);
					}
					else{
						pw.print("\t" + l.getValue(1));
					}
					if(columnNumber >= 3){
						pw.print("\t" + l.getValue(2));
					}
					if(columnNumber >= 4){
						pw.print("\t" + l.getValue(3));
					}
					pw.println();
				}
			} finally{
				if(p!=null)
					p.terminate();
			}
		} 
		catch(FileNotFoundException e){
			throw new RuntimeException(e);
		}
		catch (InvalidCSVSetUpException e) {
			log.error(ExceptionUtil.getStackTrace(e));
			throw new RuntimeException(e);
		}
		catch (IOException e) {
			log.error(ExceptionUtil.getStackTrace(e));
			throw new RuntimeException(e);
		}
		finally {
			if(pw!=null)
				pw.close();
		}
	}
	

}