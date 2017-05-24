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
package org.embl.gbcs.je.jemultiplexer.functest;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.File;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.embl.cg.utilitytools.utils.StringUtil;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class BaseTesterForJemultiplexer {
	private static Logger log = LoggerFactory.getLogger(BaseTesterForJemultiplexer.class);

	

	protected boolean checkReadsInFile(File outdir, String resultFastqFilename, String expectFilename) throws URISyntaxException {

		File resFile = new File(outdir, resultFastqFilename) ;
		File expectedFile = new File(JemultiplexerFunctionalTest.class.getResource(expectFilename).toURI()) ;

		Set<String> resReads = fetchReadIdSet(resFile);
		Set<String> expectedReads = fetchReadIdSet(expectedFile);

		return resReads.equals(expectedReads);
	}

	/**
	 * extracts the set of read ids (ignore the part right the first space)
	 * @param fastq
	 * @return
	 */
	protected Set<String> fetchReadIdSet(File fastq) {
		FastqReader fqr = new FastqReader(fastq);
		Iterator<FastqRecord> it = fqr.iterator();
		Set<String> s = new TreeSet<String>();
		while(it.hasNext()){
			FastqRecord r = it.next();
			String id = r.getReadHeader().split("\\s+")[0];
			s.add(id); //remove the part after the space
		}
		return s;
	}
	
	
	protected boolean checkReadLength(File fastq, int l) {
		
		FastqReader fqr = new FastqReader(fastq);
		Iterator<FastqRecord> it = fqr.iterator();
		while(it.hasNext()){
			FastqRecord r = it.next();
			if(r.getReadString().length() != l)
				return false;
		}
		return true;
	}
	
	/**
	 * only gets the part added by jemultiplexer ie the string following the "1:N:0:". 
	 * What precedes is used as read id
	 * 
	 * @param fastq
	 * @return a map keyed by the part left to the "1:N:0:" of the header and valued with the rigth part
	 */
	protected Map<String, String> fetchBarcodesInHeader(File fastq) {
		Pattern p = Pattern.compile("^(.+)\\d\\:N\\:0\\:(.+)$");
		FastqReader fqr = new FastqReader(fastq);
		Iterator<FastqRecord> it = fqr.iterator();
		Map<String, String> m = new HashMap<String, String>();
		while(it.hasNext()){
			FastqRecord r = it.next();
			Matcher matcher = p.matcher(r.getReadHeader());
			matcher.matches();
			m.put(matcher.group(1), matcher.group(2));
		}
		fqr.close();
		return m;
	}


	/**
	 * only gets the part added by jemultiplexer :
	 * -takes the string after the space 
	 * - clip off the "1:N:0:" part of it
	 * 
	 * @param fastq
	 * @return a map keyed by the part left to the space of the header
	 */
	protected Map<String, String[]> fetchBarcodesInSameHeader(File fastq) {
		FastqReader fqr = new FastqReader(fastq);
		Iterator<FastqRecord> it = fqr.iterator();
		Map<String, String[]> m = new HashMap<String, String[]>();
		while(it.hasNext()){
			FastqRecord r = it.next();
			log.info(r.getReadHeader());
			String [] tokens = r.getReadHeader().split(":");
			String h = StringUtil.mergeArray(tokens, ":", 0, 7);
			m.put(h, StringUtil.subArray(tokens, 7, tokens.length-1));
		}
		fqr.close();
		
		return m;
	}


	
}
