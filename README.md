# Je

The main public repository is at [github](https://github.com/gbcs-embl/Je/) where issues or pull request can be created.

Additional documentation and support can be found at http://gbcs.embl.de/je

## Installation

* Install from the bioconda channel with `` conda install -c bioconda je-suite ``
* Or, download the ``je_<version>.tar.gz`` from the ``dist/`` directory and unpack

## The Je tool suite

Je currently offers 4 tools:

* **je clip**

	to remove UMIs contained in reads of fastq files that do not need sample demultiplexing

* **je demultiplex**

	to demultiplex multi-samples fastq files which reads contain barcodes and UMIs (or not)

* **je demultiplex-illu**

 	to demultiplex fastq files according to associated index files (contain the sample encoding barcodes).
 	Reads can additionally contain UMIs (inline)

* **je markdupes**

 	to filter BAM files for read duplicates taking UMIs into account


### Distributions

* ``dist/``

    contains the different Je versions for download

* Bioconda

	starting from version 1.2 je-suite can be installed through conda: https://anaconda.org/bioconda/je-suite

### Source

* ``src/shell/je``

    is the wrapper script to call ``java -jar je_1.0_bundle.jar``

* ``src/galaxy/``

    contains the Je wrappers for Galaxy

* ``src/test/``

    holds the different test data
