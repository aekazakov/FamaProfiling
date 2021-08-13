# -*- coding: utf-8 -*-
import os
import time
import unittest
import shutil
from configparser import ConfigParser

from FamaProfiling.FamaProfilingImpl import FamaProfiling
from FamaProfiling.FamaProfilingServer import MethodContext
from FamaProfiling.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.GenomeFileUtilClient import GenomeFileUtil


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
        shutil.copy(os.path.join('data', 'MR-1.gbff'), cls.test_directory_path)
        shutil.copy(os.path.join('data', 'SB2B.gbff'), cls.test_directory_path)
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
        gu = GenomeFileUtil(cls.callback_url)
        genome1_params = {'file': {'path': os.path.join(cls.test_directory_path, 'MR-1.gbff')},
                          'genome_name': 'Shewanella_oneidensis_MR1',
                          'workspace_name': cls.ws_info[1]
                          }
        cls.genome1_ref = gu.genbank_to_genome(genome1_params)['genome_ref']
        genome2_params = {'file': {'path': os.path.join(cls.test_directory_path, 'SB2B.gbff')},
                          'genome_name': 'Shewanella_amazonensis_SB2B',
                          'workspace_name': cls.ws_info[1]
                          }
        cls.genome2_ref = gu.genbank_to_genome(genome2_params)['genome_ref']
        elements = {}
        elements[cls.genome1_ref] = dict()
        elements[cls.genome1_ref]['ref'] = cls.genome1_ref
        elements[cls.genome2_ref] = dict()
        elements[cls.genome2_ref]['ref'] = cls.genome2_ref
        test_GenomeSet = {'description': 'Test GenomeSet',
                          'elements': elements}
        provenance = [{}]
        provenance[0]['input_ws_objects'] = [cls.genome1_ref, cls.genome2_ref]
        provenance[0]['service'] = 'kb_SetUtilities'
        provenance[0]['method'] = 'KButil_Build_GenomeSet'
        genome_set_info = cls.wsClient.save_objects({'workspace': cls.ws_info[1],
                                                     'objects': [{'type': 'KBaseSearch.GenomeSet',
                                                                  'data': test_GenomeSet,
                                                                  'name': 'Test_GenomeSet',
                                                                  'meta': {},
                                                                  'provenance': provenance}]})[0]
        cls.genomeset_ref = "{}/{}/{}".format(genome_set_info[6],
                                              genome_set_info[0], genome_set_info[4])
        attribute_mapping_data = {"attributes": [{"attribute": "name", "source": "Fama"},
                                                 {"attribute": "description", "source": "Fama"},
                                                 {"attribute": "category", "source": "Fama"}],
                                  "instances": {
            "AmoA_PmoA": ["AmoA_PmoA",
                          "amoA-pmoA; methane/ammonia monooxygenase subunit A [EC:1.14.18.3 1.14.99.39]",
                          "Ammonium oxidation"],
            "AmoB_PmoB": ["AmoB_PmoB",
                          "amoB-pmoB; methane/ammonia monooxygenase subunit B",
                          "Ammonium oxidation"],
            "AmoC_PmoC": ["AmoC_PmoC",
                          "amoC-pmoC; methane/ammonia monooxygenase subunit C",
                          "Ammonium oxidation"],
            "AnfG_VnfG": ["AnfG_VnfG",
                          "Nitrogenase delta subunit [EC:1.18.6.1]",
                          "Nitrogen fixation"],
            "HAO": ["HAO",
                    "hao; hydroxylamine dehydrogenase [EC:1.7.2.6]",
                    "Anaerobic ammonium oxidation"],
            "Hzo": ["Hzo",
                    "Hydrazine dehydrogenase (EC:1.7.2.8)",
                    "Anaerobic ammonium oxidation"],
            "HzsA": ["HzsA",
                     "Hydrazine synthase subunit A (EC:1.7.2.7)",
                     "Anaerobic ammonium oxidation"],
            "HzsB": ["HzsB",
                     "Hydrazine synthase subunit B (EC:1.7.2.7)",
                     "Anaerobic ammonium oxidation"],
            "HzsC": ["HzsC",
                     "Hydrazine synthase subunit C (EC:1.7.2.7)",
                     "Anaerobic ammonium oxidation"],
            "NapA": ["NapA",
                     "Periplasmic nitrate reductase precursor (EC 1.7.99.4)",
                     "Nitrate dissimilatory reduction"],
            "NapB": ["NapB",
                     "Periplasmic nitrate reductase cytochrome c550-type subunit",
                     "Nitrate dissimilatory reduction"],
            "NapC": ["NapC",
                     "Cytochrome c-type protein NapC",
                     "Nitrate dissimilatory reduction"],
            "NapD": ["NapD",
                     "Periplasmic nitrate reductase component NapD",
                     "Nitrate dissimilatory reduction"],
            "NapE": ["NapE",
                     "Periplasmic nitrate reductase component NapE",
                     "Nitrate dissimilatory reduction"],
            "NapF": ["NapF",
                     "Ferredoxin-type protein NapF (periplasmic nitrate reductase)",
                     "Nitrate dissimilatory reduction"],
            "NapG": ["NapG",
                     "Ferredoxin-type protein NapG (periplasmic nitrate reductase)",
                     "Nitrate dissimilatory reduction"],
            "NapH": ["NapH",
                     "Polyferredoxin NapH (periplasmic nitrate reductase)",
                     "Nitrate dissimilatory reduction"],
            "NapK": ["NapK",
                     "Periplasmic nitrate reductase component NapK",
                     "Nitrate dissimilatory reduction"],
            "NapL": ["NapL",
                     "Periplasmic nitrate reductase component NapL",
                     "Nitrate dissimilatory reduction"],
            "NarC": ["NarC",
                     "Respiratory nitrate reductase subunit, conjectural (EC 1.7.99.4)",
                     "Nitrate dissimilatory reduction"],
            "NarG_NxrA": ["NarG_NxrA",
                          "narG, narZ, nxrA; nitrate reductase / nitrite oxidoreductase, alpha subunit [EC:1.7.5.1 1.7.99.-]",
                          "Nitrate dissimilatory reduction"],
            "NarH_NxrB": ["NarH_NxrB",
                          "narH, narY, nxrB; nitrate reductase / nitrite oxidoreductase, beta subunit [EC:1.7.5.1 1.7.99.-]",
                          "Nitrate dissimilatory reduction"],
            "NarI": ["NarI",
                     "Respiratory nitrate reductase gamma chain (EC 1.7.99.4)",
                     "Nitrate dissimilatory reduction"],
            "NarJ": ["NarJ",
                     "narJ, narW; nitrate reductase molybdenum cofactor assembly chaperone NarJ/NarW",
                     "Nitrate dissimilatory reduction"],
            "NasA": ["NasA",
                     "Assimilatory nitrate reductase large subunit (EC:1.7.99.4)",
                     "Nitrate assimilatory reduction"],
            "NasB": ["NasB",
                     "nasB; assimilatory nitrate reductase NADH oxidase subunit [EC:1.7.99.-]",
                     "Nitrate assimilatory reduction"],
            "NasI": ["NasI",
                     "assimilatory nitrate reductase, clostridial, electron transfer subunit [EC:1.7.99.-]",
                     "Nitrite assimilation"],
            "NasJ": ["NasJ",
                     "assimilatory nitrate reductase, clostridial, NADH oxidase subunit [EC:1.7.99.-]",
                     "Nitrite assimilation"],
            "NifB": ["NifB",
                     "nifB; nitrogen fixation protein NifB",
                     "Nitrogen fixation"],
            "NifD_AnfD_VnfD": ["NifD_AnfD_VnfD",
                               "Nitrogenase alpha chain (EC 1.18.6.1) NifD/AnfD/VnfD",
                               "Nitrogen fixation"],
            "NifH_AnfH_VnfH": ["NifH_AnfH_VnfH",
                               "Nitrogenase reductase and maturation protein NifH/AnfH/VnfH",
                               "Nitrogen fixation"],
            "NifK_AnfK_VnfK": ["NifK_AnfK_VnfK",
                               "Nitrogenase beta chain (EC 1.18.6.1) NifK/AnfK/VnfK",
                               "Nitrogen fixation"],
            "NirA": ["NirA",
                     "nirA; ferredoxin-nitrite reductase [EC:1.7.7.1]",
                     "Nitrite assimilation"],
            "NirB": ["NirB",
                     "nirB; nitrite reductase (NADH) large subunit [EC:1.7.1.15]",
                     "Nitrite assimilation"],
            "NirB3": ["NirB3",
                      "Cytochrome c-552 precursor NirB",
                      "Denitrification"],
            "NirC": ["NirC",
                     "nirC; cytochrome c55X",
                     "Denitrification"],
            "NirD": ["NirD",
                     "nirD; nitrite reductase (NADH) small subunit [EC:1.7.1.15]",
                     "Nitrite assimilation"],
            "NirK": ["NirK",
                     "nirK; nitrite reductase (NO-forming) [EC:1.7.2.1]",
                     "Denitrification"],
            "NirM": ["NirM",
                     "Cytochrome c551 NirM",
                     "Denitrification"],
            "NirN": ["NirN",
                     "Nitrite reductase associated c-type cytochorome NirN",
                     "Denitrification"],
            "NirS": ["NirS",
                     "nirS; Cytochrome cd1 nitrite reductase (NO-forming) / hydroxylamine reductase [EC:1.7.2.1 1.7.99.1]",
                     "Denitrification"],
            "NirT": ["NirT",
                     "Cytochrome c-type protein NirT",
                     "Denitrification"],
            "NirU": ["NirU",
                     "assimilatory nitrite reductase, putative NADH oxidase subunit [EC:1.7.99.-]",
                     "Nitrite assimilation"],
            "NosZ": ["NosZ",
                     "nosZ; nitrous-oxide reductase [EC:1.7.2.4]",
                     "Denitrification"],
            "NrfA": ["NrfA",
                     "nrfA; nitrite reductase (cytochrome c-552) [EC:1.7.2.2]",
                     "Ammonification"],
            "NrfB": ["NrfB",
                     "nrfB; cytochrome c-type protein NrfB",
                     "Ammonification"],
            "NrfC": ["NrfC",
                     "nrfC; protein NrfC",
                     "Ammonification"],
            "NrfD": ["NrfD",
                     "nrfD; protein NrfD",
                     "Ammonification"],
            "NrfH": ["NrfH",
                     "nrfH; cytochrome c nitrite reductase small subunit",
                     "Ammonification"],
            "UreA": ["UreA",
                     "ureA; urease subunit gamma [EC:3.5.1.5]",
                     "Urease"],
            "UreB": ["UreB",
                     "ureB; urease subunit beta [EC:3.5.1.5]",
                     "Urease"],
            "UreC": ["UreC",
                     "ureC; urease subunit alpha [EC:3.5.1.5]",
                     "Urease"],
            "cNor-C": ["cNor-C",
                       "Nitric-oxide reductase subunit C (EC 1.7.99.7)",
                       "Denitrification"],
            "cNorB_qNor": ["cNorB_qNor",
                           "Nitric-oxide reductase (EC 1.7.99.7)",
                           "Denitrification"]
        },
            "ontology_mapping_method": "User curation"
        }
        am_info = cls.wsClient.save_objects({'workspace': cls.ws_info[1],
                                             'objects': [{'type': 'KBaseExperiments.AttributeMapping',
                                                          'data': attribute_mapping_data,
                                                          'name': 'Test_row_AttributeMapping',
                                                          'meta': {},
                                                          'provenance': [{}]}]})[0]
        row_attribute_mapping_ref = "{}/{}/{}".format(am_info[6], am_info[0], am_info[4])
        attribute_mapping_data = {"attributes": [{"attribute": "sample_id", "source": "KBase"}],
                                  "instances": {"Fama_test_dummy_id1": ["Fama_test_dummy_id1"], "Fama_test_dummy_id2": ["Fama_test_dummy_id2"]},
                                  "ontology_mapping_method": "User curation"}
        am_info = cls.wsClient.save_objects({'workspace': cls.ws_info[1],
                                             'objects': [{'type': 'KBaseExperiments.AttributeMapping',
                                                          'data': attribute_mapping_data,
                                                          'name': 'Test_col_AttributeMapping',
                                                          'meta': {},
                                                          'provenance': [{}]}]})[0]
        col_attribute_mapping_ref = "{}/{}/{}".format(am_info[6], am_info[0], am_info[4])
        trait_matrix_data = {"col_attributemapping_ref": col_attribute_mapping_ref,
                             "data": {"col_ids": ["Fama_test_dummy_id1", "Fama_test_dummy_id2"],
                                      "row_ids": ["AmoA_PmoA", "AmoB_PmoB", "AmoC_PmoC", "HAO", "HzsA", "NapA", "NapB", "NapC", "NapD", "NapF", "NapG", "NapH", "NapL", "NarC", "NarG_NxrA", "NarH_NxrB", "NarI", "NarJ", "NasA", "NasB", "NasI", "NasJ", "NifB", "NifD_AnfD_VnfD", "NifH_AnfH_VnfH", "NifK_AnfK_VnfK", "NirA", "NirB", "NirC", "NirD", "NirK", "NirM", "NirN", "NirS", "NirT", "NirU", "NosZ", "NrfA", "NrfB", "NrfC", "NrfD", "NrfH", "UreA", "UreB", "UreC", "cNor-C", "cNorB_qNor"],
                                      "values": [[29.0, 1862.0], [9.0, 1502.0], [20.0, 1775.0], [1.0, 6.0], [0.0, 1.0], [10.0, 335.0], [0.0, 11.0], [1.0, 16.0], [0.0, 10.0], [1.0, 27.0], [4.0, 47.0], [0.0, 26.0], [0.0, 1.0], [5.0, 424.0], [86.0, 8420.0], [46.0, 5446.0], [7.0, 128.0], [3.0, 45.0], [157.0, 5582.0], [12.0, 168.0], [0.0, 2.0], [0.0, 2.0], [1.0, 0.0], [0.0, 3.0], [1.0, 0.0], [2.0, 3.0], [129.0, 9531.0], [73.0, 994.0], [0.0, 2.0], [19.0, 1622.0], [103.0, 10225.0], [0.0, 2.0], [0.0, 9.0], [0.0, 16.0], [0.0, 3.0], [8.0, 22.0], [4.0, 149.0], [7.0, 76.0], [0.0, 1.0], [0.0, 83.0], [0.0, 18.0], [0.0, 34.0], [80.0, 2687.0], [83.0, 2689.0], [246.0, 10558.0], [4.0, 159.0], [19.0, 523.0]]
                                      },
                             "row_attributemapping_ref": row_attribute_mapping_ref,
                             "scale": "raw"}
        tm_info = cls.wsClient.save_objects({'workspace': cls.ws_info[1],
                                             'objects': [{'type': 'KBaseMatrices.TraitMatrix',
                                                          'data': trait_matrix_data,
                                                          'name': 'Test_TraitMatrix',
                                                          'meta': {},
                                                          'provenance': [{}]}]})[0]
        trait_matrix_ref = "{}/{}/{}".format(tm_info[6], tm_info[0], tm_info[4])
        func_profile_data = {"base_object_ref": trait_matrix_ref,
                             "data": {"col_ids": ["Fama_test_dummy_id1", "Fama_test_dummy_id2"],
                                      "row_ids": ["AmoA_PmoA", "AmoB_PmoB", "AmoC_PmoC", "HAO", "HzsA", "NapA", "NapB", "NapC", "NapD", "NapF", "NapG", "NapH", "NapL", "NarC", "NarG_NxrA", "NarH_NxrB", "NarI", "NarJ", "NasA", "NasB", "NasI", "NasJ", "NifB", "NifD_AnfD_VnfD", "NifH_AnfH_VnfH", "NifK_AnfK_VnfK", "NirA", "NirB", "NirC", "NirD", "NirK", "NirM", "NirN", "NirS", "NirT", "NirU", "NosZ", "NrfA", "NrfB", "NrfC", "NrfD", "NrfH", "UreA", "UreB", "UreC", "cNor-C", "cNorB_qNor"],
                                      "values": [[0.5877623222166435, 0.6194552807793001], [0.18772393783341546, 0.5158019621326551], [0.42093693109943287, 0.6389877038101813], [0.01052927042506581, 9.44285005365384E-4], [0.0, 1.163272453732729E-4], [0.07533916331087316, 0.03985967408079665], [0.0, 0.004004703811676381], [0.023525707440798897, 0.005726992499611547], [0.0, 0.005001872240847319], [0.03374102192943467, 0.012378733183834723], [0.09317287278448243, 0.01636684276395418], [0.0, 0.007156293503795473], [0.0, 2.558288744954014E-4], [0.09551482894790636, 0.11359136011797127], [0.4838004401128097, 0.7117012180568053], [0.606624219369047, 1.1214863307558443], [0.14386930938552234, 0.039556844643411124], [0.06703428893835078, 0.013033272059964287], [1.2377859668457405, 0.6541740711814014], [0.1554388223417535, 0.025585775284706288], [0.0, 4.157415511565216E-4], [0.0, 8.79384674326921E-4], [0.011357985203073509, 0.0], [0.0, 5.408864203344862E-4], [0.017785876274710746, 0.0], [0.023945607737196006, 5.321452308325661E-4], [1.3691533051493565, 1.5159715196738617], [0.5575869192146357, 0.11053408671640833], [0.0, 0.001003764424584737], [0.6020807778521705, 0.7245394400294567], [1.2760844699978113, 2.028265195299676], [0.0, 8.562698274159332E-4], [0.0, 0.0016013558654137605], [0.0, 0.002541269397417341], [0.0, 0.001148325643963228], [0.08602893305399109, 0.0036767738796032027], [0.035240090231966045, 0.020453010032479176], [0.08339398252207368, 0.013841118202571863], [0.0, 3.631862195001992E-4], [0.0, 0.02686595241534706], [0.0, 0.004760523591057187], [0.0, 0.013633043040963725], [2.4590146922690415, 1.1212330338379586], [2.1843782229213273, 1.0658240082728703], [2.5770680795310428, 1.662618906136703], [0.08799907015919287, 0.050780650063447744], [0.17980007974121962, 0.06499031332293298]]},
                             "profile_category": "community",
                             "profile_type": "sequence reads"
                             }
        fp_info = cls.wsClient.save_objects({'workspace': cls.ws_info[1],
                                             'objects': [{'type': 'KBaseProfile.FunctionalProfile',
                                                          'data': func_profile_data,
                                                          'name': 'Test_FunctionalProfile',
                                                          'meta': {},
                                                          'provenance': [{}]}]})[0]
        cls.func_profile_ref = "{}/{}/{}".format(fp_info[6], fp_info[0], fp_info[4])

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    # @unittest.skip("")
    def test_pe_fama_profiling(self):

        ret = self.getImpl().run_FamaReadProfiling(self.getContext(), {'workspace_name': self.ws_info[1],
                                                                       'ref_dataset': 'nitrogen',
                                                                       'is_paired_end': "1",
                                                                       # tiny read set (20 + 20
                                                                       # reads), upload from local
                                                                       # files
                                                                       'read_library_refs': [self.pe_reads_ref['obj_ref']],
                                                                       'output_functional_profile_name': 'Fama_test_func_profile',
                                                                       'output_read_library_name': 'Fama_pe_test_output'})
        print('Report name', ret[0]['report_name'])
        print('Report reference', ret[0]['report_ref'])
        print(ret)

    # @unittest.skip("")
    def test_se_fama_profiling(self):
        ret = self.getImpl().run_FamaReadProfiling(self.getContext(), {'workspace_name': self.ws_info[1],
                                                                       'ref_dataset': 'nitrogen',
                                                                       'is_paired_end': "0",
                                                                       'read_library_refs': [self.se_reads_ref['obj_ref']],
                                                                       'output_functional_profile_name': 'Fama_test_func_profile',
                                                                       'output_read_library_name': 'Fama_se_test_output'})
        print('Report name', ret[0]['report_name'])
        print('Report reference', ret[0]['report_ref'])
        print(ret)

    # @unittest.skip("")
    def test_genome_fama_profiling(self):
        ret = self.getImpl().run_FamaGenomeProfiling(self.getContext(), {'workspace_name': self.ws_info[1],
                                                                         'ref_dataset': 'nitrogen',
                                                                         # S. oneidensis and S.
                                                                         # amazonensis genomes
                                                                         'genome_ref': [self.genome1_ref, self.genome2_ref],
                                                                         'output_annotation_name': 'Fama_Ncycle.',
                                                                         'output_feature_set_name': 'Fama_protein_test_output_featureset'})
        print('Report name', ret[0]['report_name'])
        print('Report reference', ret[0]['report_ref'])
        print(ret)

    # @unittest.skip("")
    def test_genomeset_fama_profiling(self):
        ret = self.getImpl().run_FamaGenomeProfiling(self.getContext(), {'workspace_name': self.ws_info[1],
                                                                         'ref_dataset': 'nitrogen',
                                                                         # two Shewanella genomes
                                                                         'genome_ref': [self.genomeset_ref],
                                                                         'output_annotation_name': 'Fama_Ncycle.',
                                                                         'output_feature_set_name': 'Fama_protein_test_output_featureset'})
        print('Report name', ret[0]['report_name'])
        print('Report reference', ret[0]['report_ref'])
        print(ret)

    # @unittest.skip("")
    def test_view_functional_profile(self):
        ret = self.getImpl().view_FamaFunctionalProfile(self.getContext(), {'workspace_name': self.ws_info[1],
                                                                            'func_profile_ref': self.func_profile_ref})
        print('Report name', ret[0]['report_name'])
        print('Report reference', ret[0]['report_ref'])
        print(ret)
