import os,sys
import uuid
from subprocess import Popen, PIPE, CalledProcessError
from collections import Counter
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.Report import generate_html_report

# This program runs functional profiling for individual FASTQ file


def run_ref_search(parser):
    print ('Starting DIAMOND')
    diamond_args = ['/kb/deployment/bin/diamond/diamond',
                    'blastx',
                    '-b', '1.0', '--tmpdir', '/dev/shm',
                    '--db',
                    parser.config.get_reference_diamond_db(parser.project.get_collection(parser.sample)),
                    '--query',
                    parser.project.get_fastq_path(parser.sample,parser.end),
                    '--out',
                    os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project.get_ref_output_name()),
                    '--max-target-seqs',
                    '50',
                    '--evalue',
                    str(parser.config.get_evalue_cutoff(parser.project.get_collection(parser.sample))),
                    '--outfmt','6','qseqid','sseqid','pident','length','mismatch','slen','qstart','qend','sstart','send','evalue','bitscore'
                    ]

    with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

    print ('DIAMOND finished')

def run_bgr_search(parser):
    print ('Starting DIAMOND')
    diamond_args = ['/kb/deployment/bin/diamond/diamond',
                    'blastx',
                    '-b', '1.0', '--tmpdir', '/dev/shm',
                    '--db',
                    parser.config.get_background_diamond_db(parser.project.get_collection(parser.sample)),
                    '--query',
                    os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project.get_ref_hits_fastq_name()),
                    '--out',
                    os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project.get_background_output_name()),
                    '--max-target-seqs',
                    '100',
                    '--evalue',
                    str(parser.config.get_background_db_size(parser.project.get_collection(parser.sample)) 
                        * parser.config.get_evalue_cutoff(parser.project.get_collection(parser.sample))
                        / parser.config.get_reference_db_size(parser.project.get_collection(parser.sample))),
                    '--outfmt','6','qseqid','sseqid','pident','length','mismatch','slen','qstart','qend','sstart','send','evalue','bitscore'
                    ]

    with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

    print ('DIAMOND finished')

def write_filtered_fastq(fastq_fwd, fastq_rev, parsers):
    reads_written = {}
    for parser in parsers:
        outfile1 = fastq_fwd
        outfile2 = fastq_rev
        if parser.end == 'pe2':
            outfile1 = fastq_rev
            outfile2 = fastq_fwd
        with open(outfile1, 'a') as of1:
            with open(outfile2, 'a') as of2:
                for read_id in parser.reads.keys():
                    if not read_id in reads_written:
                        read = parser.reads[read_id]
                        if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                            of1.write(read.read_id_line + '\n')
                            of1.write(read.sequence + '\n')
                            of1.write(read.line3 + '\n')
                            of1.write(read.quality + '\n')
                            of2.write(read.pe_id + '\n')
                            of2.write(read.pe_sequence + '\n')
                            of2.write(read.pe_line3 + '\n')
                            of2.write(read.pe_quality + '\n')
                            reads_written[read_id] = 1
            of2.closed
        of1.closed
    print ('Reads exported', str(len(reads_written)*2))


def get_read_count(infile):
    ret_val = 0
    num_lines = sum(1 for line in open(infile))
    ret_val = int(num_lines/4)
    return ret_val

def write_config_file(scratch):
    ret_val = os.path.join(scratch, 'config.ini')
    
    with open (ret_val, 'w') as of:
        # Default section
        of.write(('[DEFAULT]\n'
                'identity_cutoff = 50.0\n'
                'length_cutoff = 15\n'
                'evalue_cutoff = 0.0001\n'
                'hits_overlap_cutoff = 10\n'
                'biscore_range_cutoff = 0.03\n'))
        of.write(('taxonomy_names_file = /data/famaprofiling/names.dmp\n'
                'taxonomy_nodes_file = /data/famaprofiling/nodes.dmp\n'
                'taxonomy_merged_file = /data/famaprofiling/merged.dmp\n'))
        # Reference library for nitrogen cycle
        of.write(('\n[nitrogen]\n'
                'functions_file = /data/famaprofiling/fama_nitrogen_functions.txt\n'
                'proteins_list_file = /data/famaprofiling/fama_nitrogen_db.txt\n'
                'reference_diamond_db = /data/famaprofiling/fama_nitrogen_db.dmnd\n'
                'reference_db_size = 9537414\n'
                'background_diamond_db = /data/famaprofiling/fama_background_db.dmnd\n'
                'background_db_size = 4156041913\n'))
        # Reference library for universal markers
        of.write(('\n[universal]\n'
                'functions_file = /data/famaprofiling/fama_universal_functions.txt\n'
                'proteins_list_file = /data/famaprofiling/fama_universal_db.txt\n'
                'reference_diamond_db = /data/famaprofiling/fama_universal_db.dmnd\n'
                'reference_db_size = 18204243\n'
                'background_diamond_db = /data/famaprofiling/fama_background_db.dmnd\n'
                'background_db_size = 4156041913\n'))
        of.write('\n')
        of.closed
    
    return ret_val

def write_project_file(fastq_fwd, fastq_rev, work_dir):
    ret_val = os.path.join(work_dir, 'project.ini')
    fwd_read_count = get_read_count(fastq_fwd)
    rev_read_count = get_read_count(fastq_rev)
    
    with open (ret_val, 'w') as of:
        # Default section
        of.write(('[DEFAULT]\n'
                'project_name = \'KBase job\'\n'
                'collection = nitrogen\n' #TODO; support other collections
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
        of.write('\n\n[sample]\nsample_id = KBase Read Library\nfastq_pe1 = ')
        of.write(fastq_fwd)
        of.write('\nfastq_pe1_readcount = ')
        of.write(str(fwd_read_count))
        of.write('\nfastq_pe2 = ')
        of.write(fastq_rev)
        of.write('\nfastq_pe2_readcount = ')
        of.write(str(rev_read_count))
        of.write('\nsample_dir = ')
        of.write(work_dir)
        of.write('\nreplicate = 0\n')
        of.closed

    return ret_val

def run_fastq_profiling(config_file, project_file, end):
    parser = DiamondParser(config_file=config_file, project_file=project_file, sample='sample', end=end)
    
    if not os.path.isdir(parser.project.get_project_dir(parser.sample)):
        os.mkdir(parser_fwd.project.get_project_dir(parser.sample))
    if not os.path.isdir(os.path.join(parser.project.get_project_dir(parser.sample),
                        parser.project.get_output_subdir(parser.sample))):
        os.mkdir(os.path.join(parser.project.get_project_dir(parser.sample),
                            parser.project.get_output_subdir(parser.sample)))

    # Search in reference database
    run_ref_search(parser)
    
    # Process output of reference DB search
    parser.parse_reference_output()
    
    #Import sequence data for selected sequence reads
    print ('Reading FASTQ file')
    parser.import_fastq()
    
    print ('Exporting FASTQ ')
    parser.export_hit_fastq()
    print ('Exporting hits')
    parser.export_hit_list()
    
    # Search in background database
    run_bgr_search(parser)

    # Process output of reference DB search
    parser.parse_background_output()

    parser.export_read_fastq()
    parser.export_paired_end_reads_fastq()
    
    return parser
    

def functional_profiling_pipeline(fastq_fwd, fastq_rev, scratch):
    work_dir = os.path.join(scratch, str(uuid.uuid4()))
    os.mkdir(work_dir)
    
    config_file = write_config_file(work_dir)
    project_file = write_project_file(fastq_fwd, fastq_rev, work_dir)
    
    parser_fwd = run_fastq_profiling(config_file, project_file, 'pe1')
    parser_rev = run_fastq_profiling(config_file, project_file, 'pe2')
    
    out_dir = os.path.join(work_dir, 'out')

    #export_reads
    out_fwd_fastq = os.path.join(work_dir, 'out_fwd.fastq')
    out_rev_fastq = os.path.join(work_dir, 'out_rev.fastq')

    # TODO: write filtered fastq
    write_filtered_fastq(out_fwd_fastq, out_rev_fastq, [parser_fwd, parser_rev])

    # TODO: Generate output
    out_report = os.path.join(out_dir, 'fama_report.html')
    generate_html_report(out_report, [parser_fwd, parser_rev])

    #TODO: Krona charts generate_functions_chart(parser_fwd)

    return out_fwd_fastq, out_rev_fastq, out_report
