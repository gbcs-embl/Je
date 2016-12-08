Je
--

Additional documentation and support can be found at http://gbcs.embl.de/je


The Je tool suite
=================

Contains
++++++++
Je currently offers 4 tools :

  **je clip**

    to remove UMIs contained in reads of fastq files that do not need sample demultiplexing

  **je demultiplex**

    to demultiplex multi-samples fastq files which reads contain barcodes and UMIs (or not)

  **je demultiplex-illu**

     to demultiplex fastq files according to associated index files (contain the sample encoding barcodes).
     Reads can additionally contain UMIs (inline)

  **je markdupes**

     to filter BAM files for read duplicates taking UMIs into account


Distributions
++++++
dist/
    contains the different Je versions for download 

Source
++++++

src/shell/je
    is the wrapper script to call ``java -jar je_1.0_bundle.jar``

src/galaxy/
    contains the Je wrappers for Galaxy

src/test/
    holds the different test data