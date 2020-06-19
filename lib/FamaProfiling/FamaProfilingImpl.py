# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import uuid
import time
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace

from installed_clients.baseclient import ServerError 
#from Fama.fastq_pipeline import functional_profiling_pipeline
from fama.kbase_wrapper import pe_functional_profiling_pipeline as functional_profiling_pipeline
from fama.kbase_wrapper import protein_functional_profiling_pipeline
from fama.kbase_wrapper import save_domain_annotations
#END_HEADER


class FamaProfiling:
    '''
    Module Name:
    FamaProfiling

    Module Description:
    A KBase module: FamaProfiling
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.1.0"
    GIT_URL = "https://github.com/aekazakov/FamaProfiling.git"
    GIT_COMMIT_HASH = "2c7f322618c0dde0a9c1a9b0af61b5d8a297879b"

    #BEGIN_CLASS_HEADER
    def log(self, message, prefix_newline=False):
        print(('\n' if prefix_newline else '') + str(time.time()) + ': ' + message)
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.ws_url = config["workspace-url"]
        self.shared_folder = config['scratch']
        #END_CONSTRUCTOR
        pass


    def run_FamaProfiling(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of type "FamaProfilingParams" (Run functional
           profiling module of Fama. workspace_name - the name of the
           workspace for input/output read_library_ref - the name of the PE
           read library or SE read library ref_dataset - the name of Fama
           reference dataset output_read_library_ref - the name of the output
           filtered PE or SE read library) -> structure: parameter
           "workspace_name" of String, parameter "read_library_ref" of
           String, parameter "ref_dataset" of String, parameter
           "output_read_library_name" of String
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_FamaProfiling
        
        # Import Read Library and save as two paired-end FASTQ files
        input_ref = params['read_library_ref']
        fama_reference = params['ref_dataset']
        reads_params = {'read_libraries': [input_ref],
                        'interleaved': 'false',
                        'gzipped': None
                        }
        ru = ReadsUtils(self.callback_url)
        reads = ru.download_reads(reads_params)['files']

        print('Input reads files:')
        fwd_reads_file = reads[input_ref]['files']['fwd']
        rev_reads_file = reads[input_ref]['files']['rev']
        print('forward: ' + str(fwd_reads_file))
        print('reverse: ' + str(rev_reads_file))

        # Run Fama
        fama_output = functional_profiling_pipeline(fwd_reads_file, rev_reads_file, self.shared_folder, fama_reference, params['read_library_ref'])
        
        # Write filtered reads to workspace
        reads_params = {'fwd_file': fama_output['fwd_reads'],
                        'rev_file': fama_output['rev_reads'],
                        'sequencing_tech': reads[input_ref]['sequencing_tech'],
                        'wsname': params['workspace_name'],
                        'name': params['output_read_library_name'],
                        'interleaved': 'false'
                        }

        ru_ret = ru.upload_reads(reads_params)
        output_reads_ref = ru_ret['obj_ref']

        # Write HTML output to workspace
        message = 'Fama functional profiling finished successfully'
        
        dfu = DataFileUtil(self.callback_url)
        try:
            dfu_output = dfu.file_to_shock({'file_path': fama_output['html_report']})
        except ServerError as dfue:
            # not really any way to test this block
            self.log('Logging exception loading results to shock')
            self.log(str(dfue))
            raise

        html_links = [{'shock_id': dfu_output['shock_id'],
                       'description': 'HTML report for Fama App',
                       'name': 'fama_report.html',
                       'label': 'Fama_report'}
                      ]
        if 'krona_chart' in fama_output:
            try:
                dfu_output = dfu.file_to_shock({'file_path': fama_output['krona_chart']})
                html_links.append({'shock_id': dfu_output['shock_id'],
                                   'description': 'Krona chart for function taxonomy profile',
                                   'name': 'function_taxonomy_profile_chart.html',
                                   'label': 'Function taxonomy profile chart'}
                                  )
            except ServerError as dfue:
                # not really any way to test this block
                self.log('Logging exception loading results to shock')
                self.log(str(dfue))
                raise

        # Save report
        report_params = {'message': message,
                         'objects_created':[{'ref': output_reads_ref, 'description': 'Filtered Read Library'}],
                         'direct_html_link_index': 0,
                         'html_links': html_links,
                         'file_links': fama_output['report_files'],
                         'report_object_name': 'fama_profiling_report_' + str(uuid.uuid4()),
                         'workspace_name': params['workspace_name'],
                         'html_window_height': 460}
        try:
            report = KBaseReport(self.callback_url)
            report_info = report.create_extended_report(report_params)
        except ServerError as kre:
            # not really any way to test this block
            self.log('Logging exception saving report')
            self.log(str(kre))
            raise

        report_info['report_params'] = report_params
        self.log(str(report_info))
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref']
        }
        #END run_FamaProfiling

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_FamaProfiling return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def run_FamaProteinProfiling(self, ctx, params):
        """
        :param params: instance of type "FamaProteinProfilingParams" (Run
           protein functional profiling module of Fama. workspace_name - the
           name of the workspace for input/output genome_ref - reference to a
           genome object ref_dataset - the name of Fama reference dataset
           output_result_name - the name of the output DomainAnnotation) ->
           structure: parameter "workspace_name" of String, parameter
           "genome_ref" of list of String, parameter "ref_dataset" of String,
           parameter "output_feature_set_name" of String, parameter
           "output_annotation_name" of String
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_FamaProteinProfiling

        # Import protein sequences from input genome_ref
        input_genome_refs = params['genome_ref']
        fama_reference = params['ref_dataset']
        input_proteins = {}
        for input_genome_ref in input_genome_refs:
            input_proteins[input_genome_ref] = {}
            protein2fasta_params = {'genome_ref': input_genome_ref,
                            'include_functions': 'false',
                            'include_aliases': 'false'
                            }
            gfu = GenomeFileUtil(self.callback_url)
            proteins = gfu.genome_proteins_to_fasta(protein2fasta_params)['file_path']
            input_proteins[input_genome_ref]['fwd_path'] = proteins

        self.log('Input sequence files:', str(input_proteins))
        self.log('reference: ', fama_reference)
        ws_client = Workspace(self.ws_url)
        # Run Fama
        fama_params = {'input_proteins': input_proteins,
                       'work_dir': self.shared_folder,
                       'reference': fama_reference,
                       'ws_name': params['workspace_name'],
                       'ws_client': ws_client,
                       'featureset_name': params['output_feature_set_name'],
                       'annotation_prefix': params['output_annotation_name']
                       }
        fama_output = protein_functional_profiling_pipeline(fama_params)
        objects_created = fama_output['objects_created']

        dfu = DataFileUtil(self.callback_url)
        workspace_id = dfu.ws_name_to_id(params['workspace_name'])

        object_type = 'KBaseCollections.FeatureSet'
        save_object_params = {
            'id': workspace_id,
            'objects': [{'type': object_type,
                         'data': fama_output['feature_set_data'],
                         'name': params['output_feature_set_name']}]}

        try:
            dfu_oi = dfu.save_objects(save_object_params)[0]
        except ServerError as dfue:
            # not really any way to test this block
            self.log('Logging exception saving feature set')
            self.log(str(dfue))
            raise
        feature_set_obj_ref = "{}/{}/{}".format(dfu_oi[6], dfu_oi[0], dfu_oi[4])
        objects_created.append({'ref': feature_set_obj_ref, 'description': 'Filtered genome features'})

        self.log('FeatureSet saved to ' + feature_set_obj_ref)
        
        # Write HTML output to workspace
        message = 'Fama protein functional profiling finished successfully'
        
        try:
            dfu_output = dfu.file_to_shock({'file_path': fama_output['html_report']})
        except ServerError as dfue:
            # not really any way to test this block
            self.log('Logging exception loading results to shock')
            self.log(str(dfue))
            raise
        self.log('HTML report saved: ' + str(dfu_output))


        html_links = [{'shock_id': dfu_output['shock_id'],
                       'description': 'HTML report for Fama App',
                       'name': 'fama_report.html',
                       'label': 'Fama_report'}
                      ]
        for krona_file in fama_output['krona_charts']:
            try:
                dfu_output = dfu.file_to_shock({'file_path': krona_file})
                html_links.append({'shock_id': dfu_output['shock_id'],
                                   'description': 'Krona chart for function taxonomy profile',
                                   'name': krona_file.split('/')[-1],
                                   'label': 'Function taxonomy profile chart'}
                                  )
            except ServerError as dfue:
                # not really any way to test this block
                self.log('Logging exception loading results to shock')
                self.log(str(dfue))
                raise

        self.log('Krona chart saved: ' + str(dfu_output))

        # Save report
        report_params = {'message': message,
                         'objects_created':objects_created,
                         'direct_html_link_index': 0,
                         'html_links': html_links,
                         'file_links': fama_output['report_files'],
                         'report_object_name': 'fama_profiling_report_' + str(uuid.uuid4()),
                         'workspace_name': params['workspace_name'],
                         'html_window_height': 460}
        try:
            self.log('Call KBaseReport at ' + str(self.callback_url))
            report = KBaseReport(self.callback_url)
            self.log('Ready to save KBase report: ' + str(report_params))
            report_info = report.create_extended_report(report_params)
        except ServerError as kre:
            # not really any way to test this block
            self.log('Logging exception saving report')
            self.log(str(kre))
            raise

        report_info['report_params'] = report_params
        self.log('KBase report saved: ' + str(report_info))
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref']
        }

        #END run_FamaProteinProfiling

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_FamaProteinProfiling return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
