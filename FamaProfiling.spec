/*
A KBase module: FamaProfiling
*/

module FamaProfiling {

	/*
		Run functional profiling module of Fama.

		workspace_name - the name of the workspace for input/output
		read_library_ref - the name of the PE read library (SE library support in the future)
		output_read_library_ref - the name of the output filtered PE read library

	*/
	typedef structure {
		string workspace_name;
		string read_library_ref;
		string output_read_library_name;
	} FamaProfilingParams;

    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_FamaProfiling(FamaProfilingParams params) returns (ReportResults output) authentication required;

};
