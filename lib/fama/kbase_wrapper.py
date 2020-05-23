#!/usr/bin/python3
"""Starts Fama functional profiling pipeline from KBase environment"""
import os
import shutil
import uuid
import zipfile
from fama.utils.const import ENDS, STATUS_GOOD
from fama.utils.utils import sanitize_file_name
from fama.project.project import Project
# from fama.se_functional_pipeline import fastq_pipeline
# from fama.protein_functional_pipeline import protein_pipeline
from fama.pe_functional_pipeline import fastq_pe_pipeline
from fama.se_functional_pipeline import fastq_pipeline as fastq_se_pipeline
from fama.protein_functional_pipeline import protein_pipeline
from fama.kbase_report import generate_html_report

# How can I get reference data path from the server?
refdata_dir = '/data/famaprofiling/1.3/'


def pe_functional_profiling_pipeline(fastq_fwd, fastq_rev, scratch, ref_dataset):
    """Function calling functional profiling for fastq files"""
    work_dir = os.path.join(scratch, str(uuid.uuid4()))
    os.mkdir(work_dir)

    config_file = write_config_file(work_dir)
    project_file = write_project_file(fastq_fwd, fastq_rev, ref_dataset, work_dir)
    project = Project(config_file=config_file, project_file=project_file)

    if fastq_rev is None:
        project = fastq_se_pipeline(project)
    else:
        project = fastq_pe_pipeline(project)

    out_dir = os.path.join(work_dir, 'out')
    # export_reads
    output = {}
    out_fwd_fastq = os.path.join(work_dir, 'out_fwd.fastq')

    sample_id = project.list_samples()[0]
    if project.samples[sample_id].is_paired_end:
        out_rev_fastq = os.path.join(work_dir, 'out_rev.fastq')
    else:
        out_rev_fastq = ''

    # write filtered fastq
    write_filtered_fastq(out_fwd_fastq, out_rev_fastq, project)
    output['fwd_reads'] = out_fwd_fastq
    output['rev_reads'] = out_rev_fastq

    # Generate output
    out_report = os.path.join(out_dir, 'fama_report.html')
    generate_html_report(out_report, project)
    with zipfile.ZipFile(out_report + '.zip', 'w',
                         zipfile.ZIP_DEFLATED,
                         allowZip64=True) as zip_file:
        zip_file.write(out_report, 'fama_report.html')
    output['html_report'] = out_report + '.zip'

    # TODO: Krona charts generate_functions_chart(parser_fwd)
    report_files = {}
    if project.samples[sample_id].is_paired_end:
        metric = 'efpkg'
        if project.samples[sample_id].rpkg_scaling_factor == 0.0:
            metric = 'fragmentcount'
    else:
        metric = 'erpkg'
        if project.samples[sample_id].rpkg_scaling_factor == 0.0:
            metric = 'readcount'

    report_files[out_report + '.zip'] = 'fama_report.html'

    project_xlsx_report = sanitize_file_name(os.path.join(
        project.options.work_dir,
        project.options.project_name + '_' + metric + '_functions_taxonomy.xlsx'
        ))
    if os.path.exists(project_xlsx_report):
        report_files[project_xlsx_report] = 'function_taxonomy_profile_short.xlsx'
    else:
        print('Project XLSX file not found:', project_xlsx_report)
    sample_xlsx_report = sanitize_file_name(os.path.join(
        project.options.work_dir,
        sample_id + '_' + metric + '_functions_taxonomy.xlsx'
        ))
    if os.path.exists(sample_xlsx_report):
        report_files[sample_xlsx_report] = 'function_taxonomy_profile_full.xlsx'
    else:
        print('Sample XLSX file not found:', sample_xlsx_report)
    krona_file = sanitize_file_name(os.path.join(
        project.options.work_dir,
        sample_id + '_' + metric + '_functional_taxonomy_profile.xml.html'
        ))
    if os.path.exists(krona_file):
        krona_output = os.path.join(out_dir, 'function_taxonomy_profile_chart.html')
        shutil.copy2(krona_file, krona_output)
        with zipfile.ZipFile(krona_output + '.zip', 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:
            zip_file.write(krona_output, 'function_taxonomy_profile_chart.html')
        report_files[krona_output + '.zip'] = 'function_taxonomy_profile_chart.html'
        output['krona_chart'] = krona_output + '.zip'
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


def protein_functional_profiling_pipeline(fasta_path, scratch, ref_dataset):
    """Function calling functional profiling for protein fasta file"""
    work_dir = os.path.join(scratch, str(uuid.uuid4()))
    os.mkdir(work_dir)

    config_file = write_config_file(work_dir)
    project_file = write_project_file(fasta_path, None, ref_dataset, work_dir)
    project = Project(config_file=config_file, project_file=project_file)
    project = protein_pipeline(project)

    out_dir = os.path.join(work_dir, 'out')
    # export_reads
    output = {}

    sample_id = project.list_samples()[0]
    # Generate output
    out_report = os.path.join(out_dir, 'fama_report.html')
    generate_html_report(out_report, project)
    with zipfile.ZipFile(out_report + '.zip', 'w',
                         zipfile.ZIP_DEFLATED,
                         allowZip64=True) as zip_file:
        zip_file.write(out_report, 'fama_report.html')
    output['html_report'] = out_report + '.zip'

    # TODO: Krona charts generate_functions_chart(parser_fwd)
    report_files = {}
    metric = 'readcount'

    report_files[out_report + '.zip'] = 'fama_report.html'

    project_xlsx_report = sanitize_file_name(os.path.join(
        project.options.work_dir,
        project.options.project_name + '_' + metric + '_functions_taxonomy.xlsx'
        ))
    if os.path.exists(project_xlsx_report):
        report_files[project_xlsx_report] = 'function_taxonomy_profile_short.xlsx'
    else:
        print('Project XLSX file not found:', project_xlsx_report)
    sample_xlsx_report = sanitize_file_name(os.path.join(
        project.options.work_dir,
        sample_id + '_' + metric + '_functions_taxonomy.xlsx'
        ))
    if os.path.exists(sample_xlsx_report):
        report_files[sample_xlsx_report] = 'function_taxonomy_profile_full.xlsx'
    else:
        print('Sample XLSX file not found:', sample_xlsx_report)
    krona_file = sanitize_file_name(os.path.join(
        project.options.work_dir,
        sample_id + '_' + metric + '_functional_taxonomy_profile.xml.html'
        ))
    if os.path.exists(krona_file):
        krona_output = os.path.join(out_dir, 'function_taxonomy_profile_chart.html')
        shutil.copy2(krona_file, krona_output)
        with zipfile.ZipFile(krona_output + '.zip', 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:
            zip_file.write(krona_output, 'function_taxonomy_profile_chart.html')
        report_files[krona_output + '.zip'] = 'function_taxonomy_profile_chart.html'
        output['krona_chart'] = krona_output + '.zip'
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
    feature_ids = []
    for sample_id in project.list_samples():
        for read_id,read in project.samples[sample_id].reads['pe1'].items():
            if read.status == STATUS_GOOD:
                feature_ids.append(read_id)
    output['feature_ids'] = feature_ids
    return output


def write_config_file(scratch):
    ret_val = os.path.join(scratch, 'config.ini')

    with open(ret_val, 'w') as of:
        # Default section
        of.write(('[DEFAULT]\n'
                  'threads = 2\n'
                  'identity_cutoff = 50.0\n'
                  'length_cutoff = 15\n'
                  'evalue_cutoff = 0.0001\n'
                  'hits_overlap_cutoff = 10\n'
                  'biscore_range_cutoff = 0.03\n'
                  'aligner_path = /kb/deployment/bin/diamond/diamond\n'
                  'krona_path = /kb/deployment/bin/krona/Krona/KronaTools/scripts/ImportXML.pl\n'
                  'taxonomy_file = '))
        of.write(refdata_dir)
        of.write(('fama_taxonomy.tsv\n'
                  'microbecensus_data = '))
        of.write(refdata_dir)
        of.write('\n')

        # Reference library for nitrogen cycle
        of.write(('\n[nitrogen]\n'
                  'functions_file = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_v.10.0_functions_thresholds.tsv\n'
                  'proteins_list_file = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_v.10.0_proteins.txt\n'
                  'taxonomy_file = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_v.10.0_taxonomy.tsv\n'
                  'reference_diamond_db = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_preselection_db_v.10.0.dmnd\n'
                  'reference_db_size = 18389191\n'
                  'background_diamond_db = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_classification_db_v.10.0.dmnd\n'
                  'background_db_size = 64090633\n'))

        # Reference library for universal markers
        of.write(('\n[universal]\n'
                  'functions_file = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_v.10.0_functions_thresholds.tsv\n'
                  'proteins_list_file = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_v.10.0_proteins.txt\n'
                  'taxonomy_file = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_v.10.0_taxonomy.tsv\n'
                  'reference_diamond_db = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_preselection_db_v.10.0.dmnd\n'
                  'reference_db_size = 18389191\n'
                  'background_diamond_db = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_classification_db_v.10.0.dmnd\n'
                  'background_db_size = 64090633\n'))

        # Reference library for ribosomal protein L6
        of.write(('\n[rpl6]\n'
                  'functions_file = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_v.10.0_functions_thresholds.tsv\n'
                  'proteins_list_file = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_v.10.0_proteins.txt\n'
                  'taxonomy_file = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_v.10.0_taxonomy.tsv\n'
                  'reference_diamond_db = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_preselection_db_v.10.0.dmnd\n'
                  'reference_db_size = 18389191\n'
                  'background_diamond_db = '))
        of.write(refdata_dir)
        of.write(('fama_nitrogen-cycle_classification_db_v.10.0.dmnd\n'
                  'background_db_size = 64090633\n'))

        # Reference library for universal markers
        # of.write(('\n[universal]\n'
        #            'functions_file = /data/famaprofiling/fama_universal_functions.txt\n'
        #            'proteins_list_file = /data/famaprofiling/fama_universal_db.txt\n'
        #            'reference_diamond_db = /data/famaprofiling/fama_universal_db.dmnd\n'
        #            'reference_db_size = 18204243\n'
        #            'background_diamond_db = /data/famaprofiling/fama_background_db.dmnd\n'
        #            'background_db_size = 4156041913\n'))
        of.write('\n')

    return ret_val


def write_project_file(fastq_fwd, fastq_rev, ref_dataset, work_dir):
    ret_val = os.path.join(work_dir, 'project.ini')
    #fwd_read_count = get_read_count(fastq_fwd)
    #~ if fastq_rev is None:
        #~ fastq_rev = ''
        #~ rev_read_count = ''
    #~ else:
        #~ rev_read_count = get_read_count(fastq_rev)

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
        of.write('\n\n[sample]\nsample_id = KBase_sequences\nfastq_pe1 = ')
        of.write(fastq_fwd)
        #~ of.write('\nfastq_pe1_readcount = ')
        #~ of.write(str(fwd_read_count))
        if fastq_rev is not None:
            of.write('\nfastq_pe2 = ')
            of.write(fastq_rev)
        #~ of.write('\nfastq_pe2_readcount =')
        #~ of.write(str(rev_read_count))
        of.write('\nsample_dir = ')
        of.write(work_dir)
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


#~ def get_read_count(infile):
    #~ ret_val = 0
    #~ num_lines = sum(1 for line in open(infile))
    #~ ret_val = int(num_lines/4)
    #~ return ret_val
