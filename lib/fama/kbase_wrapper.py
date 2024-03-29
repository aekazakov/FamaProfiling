#!/usr/bin/python3
"""Starts Fama functional profiling pipeline from KBase environment"""
import os
import shutil
import uuid
import zipfile
from installed_clients.baseclient import ServerError

from fama.utils.const import ENDS, STATUS_GOOD
from fama.utils.utils import sanitize_file_name
from fama.project.project import Project
from fama.pe_functional_pipeline import fastq_pe_pipeline
from fama.se_functional_pipeline import fastq_pipeline as fastq_se_pipeline
from fama.protein_functional_pipeline import protein_pipeline
from fama.output.report import get_function_scores
from fama.kbase_report import generate_html_report, generate_protein_html_report

# How can I get reference data path from the server?
refdata_dir = '/data/famaprofiling/1.5/'
ref_model_set_names = {'nitrogen': 'Fama_nitrogen_v.10.0_function_set',
                       'universal': 'Fama_universal_v.1.4_function_set',
                       'rpl6': 'Fama_rpl6_v.1.2_function_set',
                       'cazy': 'Fama_cazymes_v.2._function_set'}


def pe_functional_profiling_pipeline(params):
    """Function calling functional profiling for fastq files"""
    work_dir = os.path.join(params['work_dir'], str(uuid.uuid4()))
    os.mkdir(work_dir)

    config_file = write_config_file(work_dir)
    project_file = write_project_file(params['input_reads'],
                                      params['reference'],
                                      work_dir,
                                      params['is_paired_end'])
    project = Project(config_file=config_file, project_file=project_file)

    if params['is_paired_end'] == "1":
        project = fastq_pe_pipeline(project)
    elif params['is_paired_end'] == "0":
        project = fastq_se_pipeline(project)
    else:
        raise ValueError('Wrong values of is_paired_end parameter', params['is_paired_end'])

    out_dir = os.path.join(work_dir, 'out')
    os.mkdir(out_dir)
    # export_reads
    output = {}
    output['krona_charts'] = {}
    out_fwd_fastq = os.path.join(work_dir, 'out_fwd.fastq')

    sample_id = project.list_samples()[0]
    if params['is_paired_end'] == "1":
        out_rev_fastq = os.path.join(work_dir, 'out_rev.fastq')
    else:
        out_rev_fastq = ''

    # write filtered fastq
    write_filtered_fastq(out_fwd_fastq, out_rev_fastq, project)
    output['fwd_reads'] = out_fwd_fastq
    if params['is_paired_end'] == "1":
        output['rev_reads'] = out_rev_fastq

    # Generate output
    out_report = os.path.join(out_dir, 'fama_report.html')
    generate_html_report(out_report, project, params['name2ref'])
    with zipfile.ZipFile(out_report + '.zip', 'w',
                         zipfile.ZIP_DEFLATED,
                         allowZip64=True) as zip_file:
        zip_file.write(out_report, 'fama_report.html')
    output['html_report'] = out_report + '.zip'

    # TODO: Krona charts generate_functions_chart(parser_fwd)
    report_files = {}
    if params['is_paired_end'] == "1":
        metric = 'efpkg'
        rawcount_flag = False
        for sample_id in project.list_samples():
            if project.samples[sample_id].rpkg_scaling_factor == 0.0:
                rawcount_flag = True
        if rawcount_flag:
            metric = 'fragmentcount'
    else:
        metric = 'erpkg'
        rawcount_flag = False
        for sample_id in project.list_samples():
            if project.samples[sample_id].rpkg_scaling_factor == 0.0:
                rawcount_flag = True
        if rawcount_flag:
            metric = 'readcount'

    # Create TraitMatrix object
    output['trait_matrix_ref'] = write_trait_matrix(project, params)
    output['functional_profile_ref'] = write_functional_profile(project,
                                                                params,
                                                                output['trait_matrix_ref'])
    report_files[out_report] = 'fama_report.html'
    project_xlsx_report = sanitize_file_name(os.path.join(
        project.options.work_dir,
        project.options.project_name + '_' + metric + '_functions.xlsx'
        ))
    if os.path.exists(project_xlsx_report):
        report_files[project_xlsx_report] = 'Functional_profiles_combined.xlsx'
    else:
        print('Project XLSX file not found:', project_xlsx_report)

    project_xlsx_report = sanitize_file_name(os.path.join(
        project.options.work_dir,
        project.options.project_name + '_' + metric + '_functions_taxonomy.xlsx'
        ))
    if os.path.exists(project_xlsx_report):
        report_files[project_xlsx_report] = 'Function_taxonomy_profiles_combined.xlsx'
    else:
        print('Project XLSX file not found:', project_xlsx_report)
    for sample_id in project.list_samples():
        sample_xlsx_report = sanitize_file_name(os.path.join(
            project.options.work_dir,
            sample_id + '_' + metric + '_functions_taxonomy.xlsx'
            ))
        if os.path.exists(sample_xlsx_report):
            report_files[sample_xlsx_report] = sanitize_file_name(
                sample_id + ' function taxonomy profile long.xlsx'
                )
        else:
            print('Sample XLSX file not found:', sample_xlsx_report)

        krona_file = sanitize_file_name(os.path.join(
            project.options.work_dir,
            sample_id + '_' + metric + '_functional_taxonomy_profile.xml.html'
            ))
        if os.path.exists(krona_file):
            krona_output = \
                sanitize_file_name(os.path.join(out_dir, sample_id +
                                   '_function_taxonomy_profile_chart.html'))
            shutil.copy2(krona_file, krona_output)
            with zipfile.ZipFile(krona_output + '.zip', 'w',
                                 zipfile.ZIP_DEFLATED,
                                 allowZip64=True) as zip_file:
                zip_file.write(krona_output,
                               sanitize_file_name(sample_id +
                                                  '_function_taxonomy_profile_chart.html'))
            report_files[krona_output] = \
                sanitize_file_name(sample_id + '_function_taxonomy_profile_chart.html')
            output['krona_charts'][krona_output + '.zip'] = \
                (sanitize_file_name(sample_id + '_function_taxonomy_profile_chart.html'),
                    sample_id + ' function taxonomy chart')
        else:
            print('Krona diagram file not found:', krona_file)

    output_files = list()
    result_file = os.path.join(project.options.work_dir, 'Fama_result.zip')
    with zipfile.ZipFile(result_file, 'w',
                         zipfile.ZIP_DEFLATED,
                         allowZip64=True) as zip_file:
        for filename in report_files:
            zip_file.write(filename, report_files[filename])
    output_files.append({'path': result_file,
                         'name': os.path.basename(result_file),
                         'label': os.path.basename(result_file),
                         'description': 'Files generated by Fama App'})
    output['report_files'] = output_files
    return output


def protein_functional_profiling_pipeline(params):
    """Function calling functional profiling for protein fasta file
    params = {'input_proteins': input_proteins,
               'work_dir': self.shared_folder,
               'reference': fama_reference,
               'ws_name': params['workspace_name'],
               'ws_client': ws_client,
               'featureset_name': params['output_feature_set_name'],
               'annotation_prefix': params['output_annotation_name'],
               'name2ref' : name2ref
             }
    """

    work_dir = os.path.join(params['work_dir'], str(uuid.uuid4()))
    os.mkdir(work_dir)

    config_file = write_config_file(work_dir)
    project_file = write_project_file(params['input_proteins'], params['reference'], work_dir)
    project = Project(config_file=config_file, project_file=project_file)
    # Run Fama
    project = protein_pipeline(project)

    out_dir = os.path.join(work_dir, 'out')
    os.mkdir(out_dir)
    # export_reads
    output = {}
    output['krona_charts'] = {}

    # Generate output
    report_files = {}
    metric = 'proteincount'
    project_xlsx_report = sanitize_file_name(os.path.join(
        project.options.work_dir,
        project.options.project_name + '_' + metric + '_functions.xlsx'
        ))
    if os.path.exists(project_xlsx_report):
        report_files[project_xlsx_report] = 'Functional_profiles_combined.xlsx'
    else:
        print('Project XLSX file not found:', project_xlsx_report)
    project_xlsx_report = sanitize_file_name(os.path.join(
        project.options.work_dir,
        project.options.project_name + '_' + metric + '_functions_taxonomy.xlsx'
        ))
    if os.path.exists(project_xlsx_report):
        report_files[project_xlsx_report] = 'Function_taxonomy_profiles_combined.xlsx'
    else:
        print('Project XLSX file not found:', project_xlsx_report)
    project_text_report = sanitize_file_name(os.path.join(
        project.options.work_dir,
        'all_proteins.list.txt'
        ))
    if os.path.exists(project_text_report):
        report_files[project_text_report] = 'proteins_list.txt'
    else:
        print('Proteins list not found:', project_text_report)

    featureset_elements = {}
    featureset_element_ordering = []
    objects_created = []
    genome_names = {}
    # Get Domain Model Set reference
    dms_ref = get_dms(params['reference'],
                      project.config.get_functions_file(project.collection),
                      params['ws_name'],
                      params['ws_client']
                      )

    for sample_id in project.list_samples():
        annotation_obj_ref, feature_ids, genome_name = \
            save_domain_annotations(project, dms_ref, params['ws_name'],
                                    params['ws_client'], params['annotation_prefix'],
                                    sample_id, params['name2ref'][sample_id])
        genome_names[sample_id] = genome_name
        objects_created.append({'ref': annotation_obj_ref,
                                'description': 'Functional annotations for genome '
                                + project.samples[sample_id].sample_name})
        for feature_id in feature_ids:
            if feature_id not in featureset_elements:
                featureset_elements[feature_id] = []
            featureset_elements[feature_id].append(params['name2ref'][sample_id])
            featureset_element_ordering.append(feature_id)

        sample_xlsx_report = sanitize_file_name(os.path.join(
            project.options.work_dir,
            sample_id + '_' + metric + '_functions_taxonomy.xlsx'
            ))
        if os.path.exists(sample_xlsx_report):
            report_files[sample_xlsx_report] = \
                sanitize_file_name(genome_name + '_function_taxonomy_profile_long.xlsx')
        else:
            print('Sample XLSX file not found:', sample_xlsx_report)
        krona_file = sanitize_file_name(os.path.join(
            project.options.work_dir,
            sample_id + '_' + metric + '_functional_taxonomy_profile.xml.html'
            ))

        if os.path.exists(krona_file):
            krona_output = \
                sanitize_file_name(os.path.join(out_dir, genome_name +
                                   '_function_taxonomy_profile_chart.html'))
            shutil.copy2(krona_file, krona_output)
            with zipfile.ZipFile(krona_output + '.zip', 'w',
                                 zipfile.ZIP_DEFLATED,
                                 allowZip64=True) as zip_file:
                zip_file.write(krona_output,
                               sanitize_file_name(genome_name +
                                                  '_function_taxonomy_profile_chart.html'))
            report_files[krona_output] = \
                sanitize_file_name(genome_name + '_function_taxonomy_profile_chart.html')
            output['krona_charts'][krona_output + '.zip'] = \
                (sanitize_file_name(genome_name + '_function_taxonomy_profile_chart.html'),
                    sample_id + ' function taxonomy chart')
        else:
            print('Krona diagram file not found:', krona_file)

    feature_set_data = {'description': 'FeatureSet generated by Fama protein profiling',
                        'element_ordering': featureset_element_ordering,
                        'elements': featureset_elements}

    out_report = os.path.join(out_dir, 'fama_report.html')
    generate_protein_html_report(out_report, project, params['name2ref'])
    with zipfile.ZipFile(out_report + '.zip', 'w',
                         zipfile.ZIP_DEFLATED,
                         allowZip64=True) as zip_file:
        zip_file.write(out_report, 'fama_report.html')
    output['html_report'] = out_report + '.zip'
    report_files[out_report] = 'fama_report.html'

    output_files = list()
    result_file = os.path.join(project.options.work_dir, 'Fama_result.zip')
    with zipfile.ZipFile(result_file, 'w',
                         zipfile.ZIP_DEFLATED,
                         allowZip64=True) as zip_file:
        for filename in report_files:
            zip_file.write(filename, report_files[filename])
    output_files.append({'path': result_file,
                         'name': os.path.basename(result_file),
                         'label': os.path.basename(result_file),
                         'description': 'Files generated by Fama App'})
    output['report_files'] = output_files
    output['project'] = project
    output['feature_set_data'] = feature_set_data
    output['objects_created'] = objects_created
    return output


def write_config_file(scratch):
    ret_val = os.path.join(scratch, 'config.ini')

    with open(ret_val, 'w') as of:
        # Default section
        config = """[DEFAULT]
threads = 16
identity_cutoff = 50.0
length_cutoff = 15
evalue_cutoff = 0.0001
hits_overlap_cutoff = 10
biscore_range_cutoff = 0.03
aligner_path = /kb/deployment/bin/diamond/diamond
krona_path = /kb/deployment/bin/krona/Krona/KronaTools/scripts/ImportXML.pl
taxonomy_file = {refdir}/fama_taxonomy.tsv
microbecensus_data = {refdir}

# Reference library for nitrogen cycle
[nitrogen]
functions_file = {refdir}/fama/nitrogen11/fama_nitrogen-cycle_v.11.0_functions_thresholds.tsv
proteins_list_file = {refdir}/fama/nitrogen11/fama_nitrogen-cycle_v.11.0_proteins.txt
taxonomy_file = {refdir}/fama/nitrogen11/fama_nitrogen-cycle_v.11.0_taxonomy.tsv
reference_diamond_db = {refdir}/fama/nitrogen11/preselection_db_nr100.dmnd
reference_db_size = 18388755
background_diamond_db = {refdir}/fama/nitrogen11/classification_db_nr100.dmnd
background_db_size = 64090197

# Reference library for universal markers
[universal]
functions_file = {refdir}/fama/universal1.4/fama_function_thresholds_v.1.4.txt
proteins_list_file = {refdir}/fama/universal1.4/fama_universal_v.1.4.txt
taxonomy_file = {refdir}/fama/universal1.4/fama_universal_taxonomy_v.1.4.txt
reference_diamond_db = {refdir}/fama/universal1.4/preselection_db_nr100.dmnd
reference_db_size = 31895938
background_diamond_db = {refdir}/fama/universal1.4/classification_db_nr100.dmnd
background_db_size = 49177580

# Reference library for cazymes
[cazy]
functions_file = {refdir}/fama/cazy2/cazy_v2_functions.txt
proteins_list_file = {refdir}/fama/cazy2/cazy_v2_proteins.txt
taxonomy_file = {refdir}/fama/cazy2/cazy_v2_taxonomy.txt
reference_diamond_db = {refdir}/fama/cazy2/preselection_database.dmnd
reference_db_size = 238821857
background_diamond_db = {refdir}/fama/cazy2/classification_database.dmnd
background_db_size = 1308452213

# Reference library for ribosomal protein L6
[rpl6]
functions_file = {refdir}/fama_rpl6_functions_thresholds_v.1.2.txt
proteins_list_file = {refdir}/fama_rpl6_proteins_v.1.2.txt
taxonomy_file = {refdir}/fama_rpl6_taxonomy_v.1.2.txt
reference_diamond_db = {refdir}/fama_rpl6_preselection_db_v.1.2.dmnd
reference_db_size = 2095606
background_diamond_db = {refdir}/fama_rpl6_classification_db_v.1.2.dmnd
background_db_size = 4769946
""".format(refdir=refdata_dir)
        of.write(config)
    return ret_val


def write_project_file(input_paths, ref_dataset, work_dir, is_paired_end="0"):
    ret_val = os.path.join(work_dir, 'project.ini')

    with open(ret_val, 'w') as of:
        # Default section
        of.write(('[DEFAULT]\n'
                  'project_name = \'Fama profiling, ' + ref_dataset + ' reference set\'\n'
                  'collection = ' + ref_dataset + '\n'
                  'ref_output_name = ref_tabular_output.txt\n'
                  'background_output_name = bgr_tabular_output.txt\n'
                  'ref_hits_list_name = ref_hits.txt\n'
                  'ref_hits_fastq_name = ref_hits.fq\n'
                  'reads_fastq_name = reads.fq\n'
                  'pe_reads_fastq_name = reads_pe.fq\n'
                  'output_subdir = out\n'
                  'report_name = report.txt\n'
                  'xml_name = krona.xml\n'
                  'html_name = functional_profile.html\n'
                  'reads_json_name = reads.json\n'
                  'work_dir = '))
        of.write(work_dir)
        for sample_id in input_paths:
            of.write('\n\n[' + sanitize_sample_id(sample_id) + ']\nsample_id = ' +
                     sample_id + '\nfastq_pe1 = ')
            of.write(input_paths[sample_id]['fwd'])
            if is_paired_end == "1":
                of.write('\nfastq_pe2 = ')
                of.write(input_paths[sample_id]['rev'])
            of.write('\nsample_dir = ')
            of.write(os.path.join(work_dir, sanitize_sample_id(sample_id)))
            of.write('\nreplicate = 0\n')
    return ret_val


def write_filtered_fastq(fastq_fwd, fastq_rev, project):
    reads_written = set()
    if fastq_rev == '':
        outfile1 = fastq_fwd
        with open(outfile1, 'a') as of1:
            for sample_id in project.list_samples():
                for read_id, read in project.samples[sample_id].reads['pe1'].items():
                    if read_id in reads_written:
                        continue
                    if read.status == STATUS_GOOD:
                        of1.write(read.read_id_line + '\n')
                        of1.write(read.sequence + '\n')
                        of1.write(read.line3 + '\n')
                        of1.write(read.quality + '\n')
                        reads_written.add(read_id)
    else:
        for sample_id in project.list_samples():
            for end_id in ENDS:
                outfile1 = fastq_fwd
                outfile2 = fastq_rev
                if end_id == 'pe2':
                    outfile1 = fastq_rev
                    outfile2 = fastq_fwd
                with open(outfile1, 'a') as of1:
                    with open(outfile2, 'a') as of2:
                        for read_id, read in project.samples[sample_id].reads[end_id].items():
                            if read_id in reads_written:
                                continue
                            if read.status == STATUS_GOOD:
                                of1.write(read.read_id_line + '\n')
                                of1.write(read.sequence + '\n')
                                of1.write(read.line3 + '\n')
                                of1.write(read.quality + '\n')
                                of2.write(read.pe_id + '\n')
                                of2.write(read.pe_sequence + '\n')
                                of2.write(read.pe_line3 + '\n')
                                of2.write(read.pe_quality + '\n')
                                reads_written.add(read_id)
    print('Reads exported', str(len(reads_written)*2))


def create_dms(ref_path, ref_name, ref_version, ws_name, ws_client):
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
            'source_url': 'https://iseq.lbl.gov/fama/reference/' + ref_name + '/',
            'version': ref_version,
            'release_date': date,
            'program': program_version,
            'domain_prefix': '',
            'dbxref_prefix': '',
            'library_files': [],
            'domains': models
            }
    ret = ws_client.save_objects({'workspace': ws_name,
                                  'objects': [{'name': domain_source + '_' + ref_name +
                                               '_v.' + ref_version + '_functions',
                                               'type': 'KBaseGeneFamilies.DomainLibrary',
                                               'data': dlib}]})
    print(str(ret))
    dlib_id = "{}/{}/{}".format(ret[0][6], ret[0][0], ret[0][4])
    # make domain model set
    dms = {'set_name': domain_source + '_' + ref_name + '_v.' + ref_version,
           'domain_libs': {'': dlib_id},
           'domain_prefix_to_dbxref_url': {'': 'https://iseq.lbl.gov/fama/reference/' +
                                           ref_name + '/'},
           'domain_accession_to_description': acc2descr
           }
    ret = ws_client.save_objects({'workspace': ws_name,
                                  'objects': [{'name': domain_source + '_' + ref_name + '_v.' +
                                               ref_version + '_function_set',
                                               'type': 'KBaseGeneFamilies.DomainModelSet',
                                               'data': dms}]})
    print(str(ret))
    dms_id = "{}/{}/{}".format(ret[0][6], ret[0][0], ret[0][4])
    return dms_id


def get_dms(reference_id, ref_path, ws_name, ws_client):

    # check if DomainModelSet exists
    dms_type = 'KBaseGeneFamilies.DomainModelSet'
    dms_ref = None
    try:
        dmss = ws_client.list_objects({'type': dms_type, 'workspaces': [ws_name]})
        for obj_data in dmss:
            if obj_data[1] == ref_model_set_names[reference_id]:
                dms_ref = "{}/{}/{}".format(obj_data[6], obj_data[0], obj_data[4])
                break
    except ServerError as wse:
        print('Logging exception searching for domain module set')
        print(str(wse))

    if dms_ref is None:
        dms_ref = create_dms(ref_path,
                             reference_id,
                             ref_model_set_names[reference_id].split('_')[-3],
                             ws_name,
                             ws_client)
    return dms_ref


def save_domain_annotations(project, dms_ref, ws_name, ws_client, name_prefix,
                            sample_id, genome_ref):

    ret = ws_client.get_objects2({'objects': [{'ref': genome_ref}]})
    genome = ret['data'][0]['data']
    genome_name = ret['data'][0]['info'][1]
    annotated_features = set()
    for read_id, read in project.samples[sample_id].reads['pe1'].items():
        if read.status == STATUS_GOOD:
            annotated_features.add(read_id)
    data = {}
    contig_to_size_and_feature_count = {}
    feature_to_contig_and_index = {}
    feature_count = 0
    max_start = 0
    prev_contig = None
    feature_ids = []
    for feature_index, feature in enumerate(genome['features']):
        feature_location = feature['location'][0]
        feature_count += 1
        contig = feature_location[0]
        if contig not in contig_to_size_and_feature_count and prev_contig is not None:
            contig_to_size_and_feature_count[prev_contig] = (feature_count, max_start)
            feature_count = 0
            max_start = 0
        max_start = feature_location[1]
        prev_contig = contig
        identifier = None
        if feature['id'] in annotated_features:
            identifier = feature['id']
        elif 'cdss' in feature and feature['cdss'][0] in annotated_features:
            identifier = feature['cdss'][0]
        elif 'aliases' in feature:
            for alias in feature['aliases']:
                if isinstance(alias, list):
                    if alias[1] in annotated_features:
                        identifier = alias
                        break
                else:
                    # for compatibility with KBase.Genome version 8.3 and below
                    if alias in annotated_features:
                        identifier = alias
                        break
        if identifier is None:
            continue
        feature_ids.append(identifier)
        start = feature_location[1]
        if feature_location[2] == '+':
            strand = 1
            end = start + feature_location[3] - 1
        else:
            strand = -1
            end = start - feature_location[3] + 1
        feature_mappings = {}
        gene = project.samples[sample_id].reads['pe1'][identifier]
        gene_functions = gene.functions.keys()
        for hit in gene.hit_list.hits:
            for hit_function in hit.functions:
                if hit_function in gene_functions and hit_function not in feature_mappings:
                    feature_mappings[hit_function] = [(hit.q_start,
                                                       hit.q_end,
                                                       hit.evalue,
                                                       hit.bitscore,
                                                       hit.length/len(gene.sequence))]
        annotation_element = (identifier, start, end, strand, feature_mappings)
        if contig not in data:
            data[contig] = []
        data[contig].append(annotation_element)
        feature_to_contig_and_index[identifier] = (contig, feature_index)
    contig_to_size_and_feature_count[prev_contig] = (feature_count, max_start)

    annotation_obj = {'genome_ref': genome_ref,  # project.samples[sample_id].sample_name,
                      'used_dms_ref': dms_ref,
                      'data': data,
                      'contig_to_size_and_feature_count': contig_to_size_and_feature_count,
                      'feature_to_contig_and_index': feature_to_contig_and_index
                      }

    ret = ws_client.save_objects({'workspace': ws_name,
                                  'objects': [{'name': name_prefix + genome_name,
                                               'type': u'KBaseGeneFamilies.DomainAnnotation',
                                               'data': annotation_obj}]})
    print(ret)
    result = "{}/{}/{}".format(ret[0][6], ret[0][0], ret[0][4])
    return result, feature_ids, genome_name


def get_ref_attributemapping(project, params):
    """Creates AttributeMapping objects for reference library"""
    ws_client = params['ws_client']
    row_attributemapping_ref = None
    try:
        ams = ws_client.list_objects({'type': 'KBaseExperiments.AttributeMapping',
                                      'workspaces': [params['ws_name']]})
        for obj_data in ams:
            if obj_data[1] == ref_model_set_names[params['reference']] + '_mapping':
                row_attributemapping_ref = "{}/{}/{}".format(obj_data[6], obj_data[0], obj_data[4])
                break
    except ServerError as wse:
        print('Logging exception searching for domain module set')
        print(str(wse))

    if row_attributemapping_ref is None:
        attributes = []
        attributes.append({'attribute': 'name', 'source': 'Fama'})
        attributes.append({'attribute': 'description', 'source': 'Fama'})
        attributes.append({'attribute': 'category', 'source': 'Fama'})
        instances = {}
        with open(project.config.get_functions_file(project.collection), 'r') as infile:
            for line in infile:
                row = line.rstrip('\n\r').split('\t')
                instances[row[0]] = [row[0], row[1], row[2]]
        row_attributemapping = {'instances': instances,
                                'attributes': attributes,
                                'ontology_mapping_method': 'User curation'}

        ret = ws_client.save_objects({'workspace': params['ws_name'],
                                      'objects': [{'name': ref_model_set_names[params['reference']]
                                                   + '_mapping',
                                                   'type': 'KBaseExperiments.AttributeMapping',
                                                   'data': row_attributemapping}]})
        print(str(ret))
        row_attributemapping_ref = "{}/{}/{}".format(ret[0][6], ret[0][0], ret[0][4])
    return row_attributemapping_ref


def write_trait_matrix(project, params):
    """Creates TraitMatrix object"""
    ws_client = params['ws_client']

    if params['is_paired_end'] == "1":
        metric = 'fragmentcount'
    else:
        metric = 'readcount'
    scores = get_function_scores(project, metric=metric)

    values = []
    row_ids = []
    col_ids = []
    for sample_id in project.list_samples():
        col_ids.append(sample_id)
    for function_id in sorted(scores.keys()):
        row_ids.append(function_id)
        row_values = []
        for sample_id in project.list_samples():
            if metric in scores[function_id][sample_id]:
                row_values.append(float(scores[function_id][sample_id][metric]))
            else:
                row_values.append(0.0)
        values.append(row_values)
    data = {'row_ids': row_ids, 'col_ids': col_ids, 'values': values}
    # Create AttributeMapping for columns
    attributes = []
    attributes.append({'attribute': 'sample_id', 'source': 'KBase'})
    instances = {}
    for sample_id in project.list_samples():
        instances[sample_id] = [sample_id]
    col_attributemapping = {'instances': instances,
                            'attributes': attributes,
                            'ontology_mapping_method': 'User curation'}
    ret = ws_client.save_objects({'workspace': params['ws_name'],
                                  'objects': [{'name': params['output_functional_profile_name'] +
                                               '.rawmatrix',
                                               'type': 'KBaseExperiments.AttributeMapping',
                                               'data': col_attributemapping,
                                               'provenance': [{'input_ws_objects':
                                                               params['input_read_refs']}]}]})
    print(str(ret))
    col_attributemapping_ref = "{}/{}/{}".format(ret[0][6], ret[0][0], ret[0][4])
    # find AttributeMapping for rows or create it
    row_attributemapping_ref = get_ref_attributemapping(project, params)
    trait_matrix = {'scale': 'raw',
                    'data': data,
                    'col_attributemapping_ref': col_attributemapping_ref,
                    'row_attributemapping_ref': row_attributemapping_ref}
    ret = ws_client.save_objects({'workspace': params['ws_name'],
                                  'objects': [{'name': params['output_read_library_name']
                                               + '_rawmatrix',
                                               'type': 'KBaseMatrices.TraitMatrix',
                                               'data': trait_matrix}]})
    print(str(ret))
    tm_id = "{}/{}/{}".format(ret[0][6], ret[0][0], ret[0][4])
    return tm_id


def write_functional_profile(project, params, base_obj_ref):
    """Creates FunctinalProfile object"""
    ws_client = params['ws_client']
    if params['is_paired_end'] == "1":
        metric = 'efpkg'
        rawcount_flag = False
        for sample_id in project.list_samples():
            if project.samples[sample_id].rpkg_scaling_factor == 0.0:
                rawcount_flag = True
        if rawcount_flag:
            metric = 'fragmentcount'
    else:
        metric = 'erpkg'
        rawcount_flag = False
        for sample_id in project.list_samples():
            if project.samples[sample_id].rpkg_scaling_factor == 0.0:
                rawcount_flag = True
        if rawcount_flag:
            metric = 'readcount'
    scores = get_function_scores(project, metric=metric)

    values = []
    row_ids = []
    col_ids = []
    for sample_id in project.list_samples():
        col_ids.append(sample_id)
    for function_id in sorted(scores.keys()):
        row_ids.append(function_id)
        row_values = []
        for sample_id in project.list_samples():
            if metric in scores[function_id][sample_id]:
                row_values.append(float(scores[function_id][sample_id][metric]))
            else:
                row_values.append(0.0)
        values.append(row_values)
    data = {'row_ids': row_ids, 'col_ids': col_ids, 'values': values}

    functional_profile = {'base_object_ref': base_obj_ref,
                          'data': data,
                          'profile_type': 'sequence reads',
                          'profile_category': 'community'}
    ret = ws_client.save_objects({'workspace': params['ws_name'],
                                  'objects': [{'name': params['output_functional_profile_name'],
                                               'type': 'KBaseProfile.FunctionalProfile',
                                               'data': functional_profile}]})
    print(str(ret))
    fp_ref = "{}/{}/{}".format(ret[0][6], ret[0][0], ret[0][4])
    return fp_ref


def sanitize_sample_id(sample_id):
    """Replaces slashes with underscores, so KBase refs could be used as sample IDs"""
    return sample_id.replace('/', '_')


def genome_proteins_to_fasta(genome_data, scratch):
    """
    This function creates protein FASTA file with all protein sequences in a genome object

    """
    outfasta = os.path.join(scratch, str(uuid.uuid4()) + '_protein.faa')
    with open(outfasta, 'w') as outfile:
        for cds in genome_data['cdss']:
            if 'parent_gene' in cds and 'protein_translation' in cds:
                outfile.write('>' + cds['parent_gene'] + '\n')
                outfile.write(cds['protein_translation'] + '\n')
    return outfasta
