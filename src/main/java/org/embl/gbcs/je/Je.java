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

import java.util.Set;
import java.util.TreeSet;

import org.embl.cg.utilitytools.utils.StringUtil;
import org.embl.gbcs.je.jeclipper.Jeclipper;
import org.embl.gbcs.je.jedropseq.Jedropseq;
import org.embl.gbcs.je.jeduplicates.MarkDuplicatesWithMolecularCode;
import org.embl.gbcs.je.jemultiplexer.Jemultiplexer;
import org.embl.gbcs.je.jemultiplexer.JemultiplexerIllumina;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;




/**
 * 
 * 
 * @author girardot
 *
 */
public class Je  {
	private static Logger log = LoggerFactory.getLogger(Je.class);

	
	public static final String COMMAND_DROPSEQ = "dropseq"; 
	public static final String COMMAND_CLIP = "clip"; 
	public static final String COMMAND_DUPES = "markdupes";
	public static final String COMMAND_MULTIPLEX = "demultiplex";
	public static final String COMMAND_MULTIPLEX_ILLUMINA = "demultiplex-illu";
	
	protected static Set<String> ALLOWED_COMMANDS = null;
	static{
		ALLOWED_COMMANDS = new TreeSet<String>();
		ALLOWED_COMMANDS.add(COMMAND_CLIP);
		ALLOWED_COMMANDS.add(COMMAND_DUPES);
		ALLOWED_COMMANDS.add(COMMAND_MULTIPLEX);
		ALLOWED_COMMANDS.add(COMMAND_MULTIPLEX_ILLUMINA);
		//ALLOWED_COMMANDS.add(COMMAND_DROPSEQ);
		
	}
	
	protected String command = null;
	
	public static void main(String[] args) {
		//we need at least one option to proceed
		if(args.length == 0 ){
			System.err.println(getUsage());
			System.exit(1); //exit with error
		}
		
		//get the first param and get rid of any prefixing '-' 
		String option = args[0];
		while(option.startsWith("-")){
			option = option.replaceFirst("-", "");
		}
		if(option.equalsIgnoreCase("h") || option.equalsIgnoreCase("help")){
			//that s ok, guy wants some help
			System.out.println(getUsage());
			System.exit(0); //exit without error
		}
		else if(option.equalsIgnoreCase("v") || option.equalsIgnoreCase("version")){
			//that s ok, guy wants some versioning
			System.out.println(getVersion());
			System.exit(0); 
		}
		else if(!ALLOWED_COMMANDS.contains(option.toLowerCase())){
			System.err.println("Unknown command name : "+option);
			System.err.println(getUsage());
			System.exit(1); //error
		}
		
		/*
		 * looks good , we delegate to proper implementation
		 */
		String [] argv = {"-h"}; // init to get help
		if(args.length > 1){
			argv = StringUtil.subArray(args, 1, args.length-1);
		}
		if(option.equalsIgnoreCase(COMMAND_CLIP)){
			new Jeclipper().instanceMainWithExit(argv); 
		}
		else if(option.equalsIgnoreCase(COMMAND_MULTIPLEX)){
			new Jemultiplexer().instanceMainWithExit(argv);
		}
		else if(option.equalsIgnoreCase(COMMAND_MULTIPLEX_ILLUMINA)){
			new JemultiplexerIllumina().instanceMainWithExit(argv);
		}
		else if(option.equalsIgnoreCase(COMMAND_DUPES)){
			new MarkDuplicatesWithMolecularCode().instanceMainWithExit(argv);
		}
		else if(option.equalsIgnoreCase(COMMAND_DROPSEQ)){
			new Jedropseq().instanceMainWithExit(argv);
		}
		else{
			System.err.println(
					"FATAL : We just reached a supposedly unreachable part of the code. Please report this bug to Je developpers indicating the options you used i.e. : \n "+
					StringUtil.mergeArray(args, " ")
					);
			System.exit(1); //error
		}
		
	}

	protected static String getUsage(){
		return "Usage:   je <command> [options] \n\n"+
				"with command in : \n"
				+"\t "+COMMAND_CLIP+"      \t\t clips molecular barcodes from fastq sequence and places them in read name headers for further use in 'dupes' module\n"
				+"\t "+COMMAND_MULTIPLEX+" \t\t demultiplex fastq file(s), with optional handling of molecular barcodes for further use in 'dupes' module\n"
				+"\t "+COMMAND_MULTIPLEX_ILLUMINA+" \t demultiplex fastq file(s) using Illumina Index files, with optional handling of molecular barcodes for further use in 'dupes' module\n"
				+"\t "+COMMAND_DUPES+"     \t\t removes read duplicates based on molecular barcodes found in read name headers (as produced by clip or plex)\n"
				//+"\t "+COMMAND_DROPSEQ+"    \t\t clips cell barcode and UMI from read 1 and adds them to header of read 2. This command is for processing drop-seq results.\n"
				+"\n"
				+"Version : "+getVersion()
				;
	}
	
	
	
	public static String getVersion(){
		return Je.class.getPackage().getImplementationVersion();
	}
	

	
}
