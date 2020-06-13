# -*- coding: utf-8 -*-
import os
import time
import unittest
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
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

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

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    @unittest.skip("")
    def test_pe_fama_profiling(self):

        ret = self.getImpl().run_FamaProfiling(self.getContext(), {'workspace_name': self.getWsName(),
                                                                    'ref_dataset': 'nitrogen',
#                                                                    'read_library_ref': '22763/10/1',  # 4706 read set (LARGE!)
#                                                                   'read_library_ref': '22763/20/1',  # 4701 2M susbset
                                                                    'read_library_ref': '22763/2/1',  # tiny read set (20 + 20 reads)
#                                                                    'read_library_ref': ru_ret['obj_ref'], # tiny read set (20 + 20 reads), upload from local files
                                                                    'output_read_library_name': 'Fama_se_test_output'})
        print ('Report name', ret[0]['report_name'])
        print ('Report reference', ret[0]['report_ref'])
        
    @unittest.skip("")
    def test_se_fama_profiling(self):
        ret = self.getImpl().run_FamaProfiling(self.getContext(), {'workspace_name': self.getWsName(),
                                                                    'ref_dataset': 'nitrogen',
                                                                    'read_library_ref': '22763/33/1',  # tiny read set (20 reads)
                                                                    'output_read_library_name': 'Fama_se_test_output'})
        print ('Report name', ret[0]['report_name'])
        print ('Report reference', ret[0]['report_ref'])
        
        
    def test_protein_fama_profiling(self):
        ret = self.getImpl().run_FamaProteinProfiling(self.getContext(), {'workspace_name': self.getWsName(),
                                                                    'ref_dataset': 'nitrogen',
                                                                    'genome_ref': '22763/32/1',  # S. oneidensis genome
                                                                    'output_annotation_name': 'Fama_protein_test_output_annotation',
                                                                    'output_feature_set_name': 'Fama_protein_test_output_featureset'})
        print ('Report name', ret[0]['report_name'])
        print ('Report reference', ret[0]['report_ref'])
        

    @unittest.skip("")
    def test_create_dms(self):
        ret = self.create_dms('/data/famaprofiling/1.4/fama_nitrogen-cycle_v.10.0_functions_thresholds.tsv', 'nitrogen', '10.0', '22763')
        
        print ('Object reference', ret)


    def create_dms(self, ref_path, ref_name, ref_version, ws_name):
        domain_source = 'Fama'
        date = '2020-06-12'
        program_version = '1.0'
        model_type = 'Protein-Sequence'
        
        # make domain models
        models = {}
        acc2descr = {}
        with open(ref_path, 'r') as infile:
            for line in infile:
                row = line.rstrip('\n\r').split('\t')
                model = {'accession': row[0],
                    'name': row[0],
                    'description': row[1],
                    'length': 0,
                    'model_type': model_type
                    }
                models[row[0]] = model
                acc2descr[row[0]] = row[1]
        # make domain library
        dlib = {'id': '',
            'source': 'Fama',
            'source_url': 'https://iseq.lbl.gov',
            'version': ref_version,
            'release_date': date,
            'program': program_version,
            'domain_prefix': '',
            'dbxref_prefix': '',
            'library_files': [],
            'domains': models
            }
        ret = self.getWsClient().save_objects({'id': ws_name, #'workspace': #self.getWsName(),
                                         'objects': [{'name': domain_source + '_' + ref_name + '_v.' + ref_version + '_functions',
                                                     'type': u'KBaseGeneFamilies.DomainLibrary',
                                                     'data': dlib}]})
        print(str(ret))
        dlib_id = "{}/{}/{}".format(ret[0][6], ret[0][0], ret[0][4])
        # make domain model set
        dms = {'set_name': domain_source + '_' + ref_name + '_v.' + ref_version,
            'domain_libs': {'': dlib_id},
            'domain_prefix_to_dbxref_url': {'':'https://iseq.lbl.gov'},
            'domain_accession_to_description': acc2descr
            }
        ret = self.getWsClient().save_objects({'id': ws_name, #'workspace': #self.getWsName(),
                                         'objects': [{'name': domain_source + '_' + ref_name + '_v.' + ref_version + '_function_set',
                                                     'type': u'KBaseGeneFamilies.DomainModelSet',
                                                     'data': dms}]})
        print(str(ret))
        dms_id = "{}/{}/{}".format(ret[0][6], ret[0][0], ret[0][4])
        return dms_id
