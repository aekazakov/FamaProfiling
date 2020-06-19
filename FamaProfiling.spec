/*
A KBase module: FamaProfiling
*/

module FamaProfiling {

    /*
        Run functional profiling module of Fama.

        workspace_name - the name of the workspace for input/output
        read_library_ref - the name of the PE read library or SE read library
        ref_dataset - the name of Fama reference dataset
        output_read_library_ref - the name of the output filtered PE or SE read library

    */
    typedef structure {
        string workspace_name;
        string read_library_ref;
        string ref_dataset;
        string output_read_library_name;
    } FamaProfilingParams;

    /*
    Run protein functional profiling module of Fama.

    workspace_name - the name of the workspace for input/output
    genome_ref - reference to a genome object
    ref_dataset - the name of Fama reference dataset
    output_result_name - the name of the output DomainAnnotation

    */
    typedef structure {
        string workspace_name;
        list<string> genome_ref;
        string ref_dataset;
        string output_feature_set_name;
        string output_annotation_name;
    } FamaProteinProfilingParams;

    /* 
    @id ws KBaseGeneFamilies.DomainAnnotation
    */
    typedef string domain_annotation_ref;
    
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_FamaProfiling(FamaProfilingParams params) returns (ReportResults output) authentication required;
    funcdef run_FamaProteinProfiling(FamaProteinProfilingParams params) returns (ReportResults output) authentication required;

};
