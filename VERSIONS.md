* Version 1.0 is as published in BMC Bioinformatics 17(1) doi: 10.1186/s12859-016-1284-2
* Version 1.1 : bug correction in header writing in paired end mode with ENSURE_IDENTICAL_HEADER_NAMES=true
* Version 1.2 : change the logic of sample selection in PE when barcodes are present in both reads and redundant. 
In the rare situation where barcode matching resolve to 2 **different** samples with no mismatch in both cases, 
the code *could* assign the PE reads always to the same sample. To avoid this is the future, the PE reads are now
unassigned in this very special situation. 