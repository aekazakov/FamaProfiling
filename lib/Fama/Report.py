#!/usr/bin/python3
import os
from collections import defaultdict,Counter,OrderedDict

from Fama.DiamondParser.hit_utils import cleanup_protein_id,autovivify

ENDS = ['pe1','pe2']
def generate_html_report(outfile, parsers):
    
    with open (outfile, 'w') as of:
        of.write('<!doctype html><html>\n<head>\n<title>Fama report</title>\n</head>\n<body>\n')
        of.write('<h1>Fama report</h1>\n')

        for parser in parsers:


            # Write general info
            of.write('\n<h3>Run info</h3>\n')
            of.write('<p>Sample ID: ' + parser.project.get_sample_id(parser.sample) + '</p>\n')
            of.write('<p>Paired end: ' + parser.end + '</p>\n')
            of.write('<p>FASTQ file: ' + parser.project.get_fastq_path(parser.sample, parser.end) + '</p>\n')
            of.write('<p>Total number of reads:' + str(parser.project.get_fastq1_readcount(parser.sample)) + '</p>\n')

            # Write read statistics
            of.write('<h3>Read statistics</h3>\n')
            read_stats = Counter()
            for read in sorted(parser.reads.keys()):
                read_stats[parser.reads[read].get_status()] += 1
            for status in OrderedDict(read_stats.most_common()):
                if status == 'unaccounted':
                    of.write('<p>Reads missing from background DB search result: ' + str(read_stats[status]) + '</p>\n')
                elif status == 'nofunction':
                    of.write('<p>Reads not mapped to any function: ' + str(read_stats[status]) + '</p>\n')
                elif status == 'function':
                    of.write('<p>Reads mapped to a function of interest: ' + str(read_stats[status]) + '</p>\n')
                else:
                    of.write('<p>' + status + ': ' + str(read_stats[status]) + '</p>\n')
            
            of.write('<h3>Function statistics</h3>\n')
            func_stats = defaultdict(float)
            func_counts = Counter()
            func_identity = defaultdict(float)
            func_hit_counts = Counter()
            for read in parser.reads.keys():
                if parser.reads[read].get_status() == 'function,besthit' or parser.reads[read].get_status() == 'function':
                    functions = parser.reads[read].get_functions()
                    for function in functions:
                        func_stats[function] += functions[function]
                        func_counts[function] += 1/len(functions)
                    for hit in parser.reads[read].get_hit_list().get_hits():
                        for function in hit.get_functions():
                            func_identity[function] += hit.get_identity()
                            func_hit_counts[function] += 1
            for function in func_identity:
                func_identity[function] = func_identity[function]/func_hit_counts[function]
            of.write('<table>\n')
            of.write('<tr><td>Function</td><td>Definition</td><td>RPKM score</td><td>Read count</td><td>Avg. identity</td></tr>\n')
            for function in sorted(func_stats.keys()):
                of.write('<tr><td>' + function + '</td><td>' 
                        + parser.ref_data.lookup_function_name(function) + '</td><td>' 
                        + str(func_stats[function]) + '</td><td>'
                        + str(func_counts[function]) + '</td><td>'
                        + str(func_identity[function]) + '</td></tr>\n')
            of.write('</table>\n')

            # Write group scores
            of.write('<H3>Function statistics by category</H3>\n')
            func_stats = defaultdict(float)
            func_counts = Counter()
            func_identity = defaultdict(float)
            func_hit_counts = Counter()
            for read in parser.reads.keys():
                if parser.reads[read].get_status() == 'function,besthit' or parser.reads[read].get_status() == 'function':
                    functions = parser.reads[read].get_functions()
                    for function in functions:
                        func_stats[parser.ref_data.lookup_function_group(function)] += functions[function]
                        func_counts[parser.ref_data.lookup_function_group(function)] += 1/len(functions)
                    for hit in parser.reads[read].get_hit_list().get_hits():
                        for function in hit.get_functions():
                            func_identity[parser.ref_data.lookup_function_group(function)] += hit.get_identity()
                            func_hit_counts[parser.ref_data.lookup_function_group(function)] += 1
            for function in func_identity:
                func_identity[function] = func_identity[function]/func_hit_counts[function]
            of.write('<table>\n')
            of.write('<tr><td>Category</td><td>RPKM score</td><td>Read count</td><td>Avg. identity</td></tr>\n')
            for function in sorted(func_stats.keys()):
                of.write('<tr><td>' + function + '</td><td>' 
                        + str(func_stats[function]) + '</td><td>'
                        + str(func_counts[function]) + '</td><td>'
                        + str(func_identity[function]) + '</td></tr>\n')
            of.write('</table>\n')

            # Write taxonomy stats
            of.write('<h3>Taxonomy statistics for best hits</h3>\n')
            tax_stats = Counter()
            identity_stats = defaultdict(float)
            rpkm_stats = defaultdict(float)
            for read in parser.reads.keys():
    #            print (read, parser.reads[read].get_status())
                if parser.reads[read].get_status() == 'function,besthit' or parser.reads[read].get_status() == 'function':
    #            if parser.reads[read].get_status() == 'function':
    #                print ('\t',parser.reads[read].get_status())
                    hits = parser.reads[read].get_hit_list().get_hits()
                    for hit in hits:
    #                    print (hit.get_subject_id())
    #                    print (parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id())))
    #                    print (hit.get_identity())
                        protein_taxid = parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                        tax_stats[protein_taxid] += 1
                        identity_stats[protein_taxid] += hit.get_identity()
                    if len(hits) == 1:
                        read_functions = parser.reads[read].get_functions()
                        for function in read_functions:
                            rpkm_stats[parser.ref_data.lookup_protein_tax(cleanup_protein_id(hits[0].get_subject_id()))] += read_functions[function]
                    else:
                        read_functions = parser.reads[read].get_functions()
                        protein_taxids = {}
                        for hit in hits:
                            hit_taxid = parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                            hit_functions = hit.get_functions()
                            for hit_function in hit_functions:
                                protein_taxids[hit_taxid] = hit_function
                        for taxid in protein_taxids:
                            if protein_taxids[taxid] in read_functions:
                                rpkm_stats[taxid] += read_functions[protein_taxids[taxid]]
                            
            tax_data = parser.taxonomy_data
            counts_per_rank, identity_per_rank, rpkm_per_rank = tax_data.get_taxonomy_profile(counts=tax_stats, identity=identity_stats, scores = rpkm_stats)

            ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus']
            for rank in ranks:
                of.write('<h5>Taxonomy report for rank ' + rank + '</h5>\n')
                of.write('<table>\n')
                of.write('<tr><td>Taxon</td><td>Read count</td><td>RPKM score</td><td>Average identity</td></tr>\n')
                
                for tax in OrderedDict(Counter(counts_per_rank[rank]).most_common()):
                    #print (tax + '\t' + str(counts_per_rank[rank][tax]) + '\t' + str(identity_per_rank[rank][tax]))
                    of.write('<tr><td>' + rank + '</td><td>' + tax + '</td><td>' 
                            + str(counts_per_rank[rank][tax]) + '</td><td>' 
                            + str(rpkm_per_rank[rank][tax]) + '</td><td>' 
                            + str(identity_per_rank[rank][tax]) + '</td></tr>\n')
                of.write('</table>\n')
            of.write('<hr>\n')

        of.write('<p>End of report</p>\n')
        of.write('</body>\n</html>\n')
        of.closed

    
