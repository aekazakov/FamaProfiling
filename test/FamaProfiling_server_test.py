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
                                                                    'output_feature_set_name': 'Fama_protein_test_output'})
        print ('Report name', ret[0]['report_name'])
        print ('Report reference', ret[0]['report_ref'])
        

