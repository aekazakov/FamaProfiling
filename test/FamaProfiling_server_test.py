# -*- coding: utf-8 -*-
import os
import time
import unittest
import shutil
from subprocess import Popen, PIPE, CalledProcessError
from configparser import ConfigParser

from FamaProfiling.FamaProfilingImpl import FamaProfiling
from FamaProfiling.FamaProfilingServer import MethodContext
from FamaProfiling.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace
from installed_clients.ReadsUtilsClient import ReadsUtils


class FamaProfilingTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('FamaProfiling'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'FamaProfiling',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = FamaProfiling(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.suffix = int(time.time() * 1000)
        cls.wsName = "test_Fama_" + str(cls.suffix)
        cls.ws_info = cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.prepare_data()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_FamaProfiling_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    @classmethod
    def prepare_data(cls):
        wd = os.getcwd()
        print('WORKING DIRECTORY', wd)
        ru = ReadsUtils(cls.callback_url)
        test_directory_name = 'fama_test_data'
        cls.test_directory_path = os.path.join(cls.scratch, test_directory_name)
        print('TEST DIRECTORY', cls.test_directory_path)
        os.makedirs(cls.test_directory_path)
        shutil.copy(os.path.join('data', 'test_fastq_pe1.fq'), cls.test_directory_path)
        shutil.copy(os.path.join('data', 'test_fastq_pe2.fq'), cls.test_directory_path)
        reads_params = {'fwd_file': os.path.join(cls.test_directory_path, 'test_fastq_pe1.fq'),
                        'rev_file': os.path.join(cls.test_directory_path, 'test_fastq_pe2.fq'),
                        'sequencing_tech': 'Illumina',
                        'wsname': cls.ws_info[1],
                        'single_genome': 0,
                        'name': 'Fama_test_pe_input',
                        'interleaved': 0
                        }        
        cls.pe_reads_ref = ru.upload_reads(reads_params)

        se_reads_params = {'fwd_file': os.path.join(cls.test_directory_path, 'test_fastq_pe1.fq'),
                        'sequencing_tech': 'Illumina',
                        'wsname': cls.ws_info[1],
                        'single_genome': 0,
                        'name': 'Fama_test_se_input'
                        }        
        cls.se_reads_ref = ru.upload_reads(se_reads_params)

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    @unittest.skip("")
    def test_pe_fama_profiling(self):

        ret = self.getImpl().run_FamaReadProfiling(self.getContext(), {'workspace_name': self.ws_info[1],
                                                                    'ref_dataset': 'nitrogen',
                                                                    'is_paired_end': "1",
#                                                                    'read_library_refs': ['22763/20/1'],  # 4701 2M susbset
                                                                    'read_library_refs': [self.pe_reads_ref['obj_ref'], '22763/2/1'], # tiny read set (20 + 20 reads), upload from local files
                                                                    'output_functional_profile_name': 'Fama_test_func_profile',
                                                                    'output_read_library_name': 'Fama_pe_test_output'})
        print ('Report name', ret[0]['report_name'])
        print ('Report reference', ret[0]['report_ref'])
        print(ret)
        
    @unittest.skip("")
    def test_se_fama_profiling(self):
        ret = self.getImpl().run_FamaReadProfiling(self.getContext(), {'workspace_name': self.ws_info[1],
                                                                    'ref_dataset': 'nitrogen',
                                                                    'is_paired_end': "0",
#                                                                    'read_library_ref': '22763/33/1',  # tiny read set (20 reads)
                                                                    'read_library_refs': [self.se_reads_ref['obj_ref']],
                                                                    'output_functional_profile_name': 'Fama_test_func_profile',
                                                                    'output_read_library_name': 'Fama_se_test_output'})
        print ('Report name', ret[0]['report_name'])
        print ('Report reference', ret[0]['report_ref'])
        print(ret)
        
    @unittest.skip("")
    def test_genome_fama_profiling(self):
        ret = self.getImpl().run_FamaGenomeProfiling(self.getContext(), {'workspace_name': self.ws_info[1],
                                                                    'ref_dataset': 'nitrogen',
                                                                    'genome_ref': ['22763/32/1', '41747/12/1'],  # S. oneidensis and S. amazonensis genomes
                                                                    'output_annotation_name': 'Fama_Ncycle.',
                                                                    'output_feature_set_name': 'Fama_protein_test_output_featureset'})
        print ('Report name', ret[0]['report_name'])
        print ('Report reference', ret[0]['report_ref'])
        print(ret)

    @unittest.skip("")
    def test_genomeset_fama_profiling(self):
        ret = self.getImpl().run_FamaGenomeProfiling(self.getContext(), {'workspace_name': self.ws_info[1],
                                                                    'ref_dataset': 'nitrogen',
                                                                    'genome_ref': ['41747/35/1'],  # two Shewanella genomes
                                                                    'output_annotation_name': 'Fama_Ncycle.',
                                                                    'output_feature_set_name': 'Fama_protein_test_output_featureset'})
        print ('Report name', ret[0]['report_name'])
        print ('Report reference', ret[0]['report_ref'])
        print(ret)

    #@unittest.skip("")
    def test_view_functional_profile(self):
        ret = self.getImpl().view_FamaFunctionalProfile(self.getContext(), {'workspace_name': self.ws_info[1],
                                                                    'func_profile_ref': '22763/77/1'})
        print ('Report name', ret[0]['report_name'])
        print ('Report reference', ret[0]['report_ref'])
        print(ret)
        
