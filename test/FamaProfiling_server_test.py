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
    def test_fama_profiling(self):

        wd = os.getcwd()
        os.chdir(self.scratch)
        with Popen(['curl','-LO','http://iseq.lbl.gov/mydocs/fama_downloads/test_fastq_pe1.fq'], stdout=PIPE, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')
        if p.returncode != 0:
            raise CalledProcessError(p.returncode, p.args)

        with Popen(['curl','-LO','http://iseq.lbl.gov/mydocs/fama_downloads/test_fastq_pe2.fq'], stdout=PIPE, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')
        if p.returncode != 0:
            raise CalledProcessError(p.returncode, p.args)
        
        #time.sleep(3.0)
        os.chdir(wd)
        
        ru = ReadsUtils(self.callback_url)
        reads_params = {'fwd_file': os.path.join(self.scratch, 'test_fastq_pe1.fq'),
                        'rev_file': os.path.join(self.scratch, 'test_fastq_pe2.fq'),
                        'sequencing_tech': 'Illumina',
                        'wsname': self.getWsName(),
                        'name': 'Fama_test_input',
                        'interleaved': 'false'
                        }

        ru_ret = ru.upload_reads(reads_params)
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                 'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        ret = self.getImpl().run_FamaProfiling(self.getContext(), {'workspace_name': self.getWsName(),
                                                                    'read_library_ref': ru_ret['obj_ref'],
                                                                    'output_read_library_name': 'Fama_test_output'})
        print ('Report name', ret[0]['report_name'])
        print ('Report reference', ret[0]['report_ref'])
        
