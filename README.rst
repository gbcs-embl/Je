Je
--

Additional documentation and support can be found at http://gbcs.embl.de/je


The Je tool suite
=================

Contains
++++++++
Je currently offers the following tools :

  **je debarcode** 
     
     demultiplexs multi-samples fastq files using user-defined input *read-layouts* and write output files following user-defined *output-layouts*.
     Replaces both **demultiplex-illu** and **demultiplex** since version 2.0.
  
  **je dropseq**
     
     to process drop-seq results: clips cell barcode and UMI from read 1 and adds them to header of read 2 (a unique output fastq is created).
  
  **je retag**
  
    extracts barcode(s) and UMI sequence(s) embedded in read names of a BAM file and migrate them to proper BAM tags.
  
  **je clip**

    removes UMIs contained in reads of fastq files that do not need sample demultiplexing
    
  **je markdupes**

     filters BAM files for read duplicates taking UMIs into account.

  **je demultiplex**

    to demultiplex multi-samples fastq files which reads contain barcodes and UMIs (or not). Deprecated since version 2.0 (use *je debarcode* instead). 

  **je demultiplex-illu**

     to demultiplex fastq files according to associated index files (contain the sample encoding barcodes).
     Reads can additionally contain UMIs (inline). Deprecated since version 2.0 (use *je debarcode* instead).

  


Distributions
++++++
dist/
    contains the different Je versions for download 

Source
++++++

src/shell/je
    is the wrapper script to call ``java -jar je_*_bundle.jar``

src/galaxy/
    contains the Je wrappers for Galaxy

src/test/
    holds the different test data