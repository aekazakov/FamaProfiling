
import errno
import logging
import os
import pandas as pd
import uuid
import shutil
import math
import json

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.kb_GenericsReportClient import kb_GenericsReport


DATA_EPISTEMOLOGY = ['measured', 'asserted', 'predicted']
PROFILE_CATEGORY = ['community',  'organism']
PROFILE_TYPE = ['amplicon', 'mg', 'modelset']


class ProfileImporter:

    @staticmethod
    def _mkdir_p(path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    @staticmethod
    def _validate_params(params, expected, opt_param=set()):
        """Validates that required parameters are present. Warns if unexpected parameters appear"""
        expected = set(expected)
        opt_param = set(opt_param)
        pkeys = set(params)
        if expected - pkeys:
            raise ValueError("Required keys {} not in supplied parameters"
                             .format(", ".join(expected - pkeys)))
        defined_param = expected | opt_param
        for param in params:
            if param not in defined_param:
                logging.warning("Unexpected parameter {} supplied".format(param))

    @staticmethod
    def _convert_size(size_bytes):
        if size_bytes == 0:
            return "0B"
        size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
        i = int(math.floor(math.log(size_bytes, 1024)))
        p = math.pow(1024, i)
        s = round(size_bytes / p, 2)
        return "%s %s" % (s, size_name[i])

    def _calculate_object_size(self, func_profile_data):
        json_size = 0
        try:
            logging.info('start calculating object size')
            json_object = json.dumps(func_profile_data).encode("utf-8")
            json_size = len(json_object)
            size_str = self._convert_size(json_size)
            logging.info('serialized object JSON size: {}'.format(size_str))
        except Exception:
            logging.info('failed to calculate object size')

        return json_size

    def _generate_visualization_content(self, func_profile_ref, output_directory):
        func_profile_data = self.dfu.get_objects(
                                            {'object_refs': [func_profile_ref]})['data'][0]['data']

        data = func_profile_data.get('data')
        data_df = pd.DataFrame(data['values'],
                               index=data['row_ids'], columns=data['col_ids'])

        data_df.fillna(0, inplace=True)
        tsv_file_path = os.path.join(output_directory, 'heatmap_data_{}.tsv'.format(
                                                                    str(uuid.uuid4())))
        data_df.to_csv(tsv_file_path)
        heatmap_dir = self.report_util.build_heatmap_html({
                                            'tsv_file_path': tsv_file_path,
                                            'cluster_data': True})['html_dir']

        row_data_summary = data_df.T.describe().to_string()
        col_data_summary = data_df.describe().to_string()

        tab_def_content = ''
        tab_content = ''

        # build data summary page
        viewer_name = 'data_summary'
        tab_def_content += '''\n<div class="tab">\n'''
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += '''>Profile Statistics</button>\n'''

        tab_content += '''\n<div id="{}" class="tabcontent" style="overflow:auto">'''.format(viewer_name)
        tab_content += '''\n<h5>Profile Size: {} x {}</h5>'''.format(len(data_df.index),
                                                                     len(data_df.columns))
        tab_content += '''\n<h5>Row Aggregating Statistics</h5>'''
        html = '''\n<pre class="tab">''' + str(row_data_summary).replace("\n", "<br>") + "</pre>"
        tab_content += html
        tab_content += '''\n<br>'''
        tab_content += '''\n<hr style="height:2px;border-width:0;color:gray;background-color:gray">'''
        tab_content += '''\n<br>'''
        tab_content += '''\n<h5>Column Aggregating Statistics</h5>'''
        html = '''\n<pre class="tab">''' + str(col_data_summary).replace("\n", "<br>") + "</pre>"
        tab_content += html
        tab_content += '\n</div>\n'

        # build profile heatmap page
        viewer_name = 'ProfileHeatmapViewer'
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += ''' id="defaultOpen"'''
        tab_def_content += '''>Profile Heatmap</button>\n'''

        heatmap_report_files = os.listdir(heatmap_dir)

        heatmap_index_page = None
        for heatmap_report_file in heatmap_report_files:
            if heatmap_report_file.endswith('.html'):
                heatmap_index_page = heatmap_report_file

            shutil.copy2(os.path.join(heatmap_dir, heatmap_report_file),
                         output_directory)

        if heatmap_index_page:
            tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
            tab_content += '\n<iframe height="1300px" width="100%" '
            tab_content += 'src="{}" '.format(heatmap_index_page)
            tab_content += 'style="border:none;"></iframe>'
            tab_content += '\n</div>\n'
        else:
            tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
            tab_content += '''\n<p style="color:red;" >'''
            tab_content += '''Heatmap is too large to be displayed.</p>\n'''
            tab_content += '\n</div>\n'

        tab_def_content += '\n</div>\n'
        return tab_def_content + tab_content

    def _generate_html_report(self, func_profile_ref):

        logging.info('Start generating report page')

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start generating html report in {}'.format(output_directory))

        html_report = list()

        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'func_profile_viewer_report.html')

        visualization_content = self._generate_visualization_content(func_profile_ref,
                                                                     output_directory)

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__),
                                   'templates', 'func_profile_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Visualization_Content</p>',
                                                          visualization_content)
                result_file.write(report_template)

        report_shock_id = self.dfu.file_to_shock({'file_path': output_directory,
                                                  'pack': 'zip'})['shock_id']

        html_report.append({'shock_id': report_shock_id,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for Import Amplicon Matrix App'
                            })
        return html_report

    def gen_func_profile_report(self, func_profile_ref, workspace_name):
        logging.info('start generating report')

        objects_created = [{'ref': func_profile_ref, 'description': 'Imported FunctionalProfile'}]

        output_html_files = self._generate_html_report(func_profile_ref)

        report_params = {'message': '',
                         'objects_created': objects_created,
                         'workspace_name': workspace_name,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 1400,
                         'report_object_name': 'func_profile_viewer_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        # self.token = config['KB_AUTH_TOKEN']
        self.dfu = DataFileUtil(self.callback_url)
        self.report_util = kb_GenericsReport(self.callback_url)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)