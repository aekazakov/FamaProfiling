# FamaProfiling release notes
=========================================

0.0.0
-----
* Module created by kb-sdk init

0.0.1
-----
* Actual app code added. Runs functional profiling for a single ReadLibrary. Available reference library limited to nitrogen cycle enzymes.

0.0.2
-----
* Version for test purposes only, reports memory usage in the log
* Fixed bug resulted in "insufficient amount of memory" error
* Switched to DIAMOND version 0.8.38

0.0.3
-----
* Fixed report creation

0.0.4
-----
* Report in HTML format

0.0.5
-----
* Updated Fama code using smaller reference database and function-specific cutoffs
* New reference datasets

0.0.6
-----
* Reference data moved

0.1.0
-----
* Added run_FamaProteinProfiling for protein-based sequence search in KBase genomes. Results are saved as DomainAnnotation and FeatureSet objects
* Updated report in HTML format
* Krona chart generation for taxonomic profiles

0.1.1
-----
* Added proteincount metric for protein profiling
* Updated method descriptions

1.0.0
----
* First released version
* Output objects of protein profiling use gene identifiers instead of protein identifiers