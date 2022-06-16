/*
A KBase module: FamaProfiling
*/

module FamaProfiling {

    /*
    A boolean - 0 for false, 1 for true.
        @range (0, 1)
    */
    typedef int bool;

    /*
        Parameters for metagenome functional profiling.

        workspace_name - the name of the workspace for input/output
        read_library_refs - references to the name of the PE read library or SE read library
        ref_dataset - the name of Fama reference dataset
        is_paired_end - 1 for paired-end library, 0 for single-end library
        output_functional_profile_name - the name of the output functional profile
        output_read_library_ref - the name of the output filtered PE or SE read library

    */
    typedef structure {
        string workspace_name;
        list<string> read_library_refs;
        string ref_dataset;
        bool is_paired_end;
        string output_functional_profile_name;
        string output_read_library_name;
    } FamaReadProfilingParams;

    /*
    Parameters for genome functional profiling.

    workspace_name - the name of the workspace for input/output
    genome_refs - references to a genome object
    ref_dataset - the name of Fama reference dataset
    output_result_name - the name of the output DomainAnnotation

    */
    typedef structure {
        string workspace_name;
        list<string> genome_ref;
        string ref_dataset;
        string output_feature_set_name;
        string output_annotation_name;
    } FamaGenomeProfilingParams;

    /*
    Parameters for functional profile viewer.

    func_profile_ref - reference to functional profile
    workspace_name - the name of the workspace for input/output

    */
    typedef structure {
        string func_profile_ref;
        string workspace_name;
    } ViewFunctionalProfileParams;

    /* 
    Output report parameters

    report_name - the name of the report object
    report_ref - the reference to the report object
    */

    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        Run metagenome functional profiling module of Fama.
    */
    funcdef run_FamaReadProfiling(FamaReadProfilingParams params) returns (ReportResults output) authentication required;

    /*
        Run genome functional profiling module of Fama.
    */
    funcdef run_FamaGenomeProfiling(FamaGenomeProfilingParams params) returns (ReportResults output) authentication required;

    /*
		View Fama Functional Profile
    */
    funcdef view_FamaFunctionalProfile(ViewFunctionalProfileParams params) returns (ReportResults output) authentication required;

};
