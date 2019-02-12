#!/usr/bin/python3
import os
from collections import defaultdict,Counter,OrderedDict

from Fama.DiamondParser.hit_utils import cleanup_protein_id,autovivify

ENDS = ['pe1','pe2']
def generate_html_report(outfile, parsers):
    
    with open (outfile, 'w') as of:


        of.write('<html>\n')
        of.write('    <head>\n')
        of.write('        <title>Fama report</title>\n')
        of.write('        <!-- Compiled and minified CSS -->\n')
        of.write('        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/materialize/1.0.0/css/materialize.min.css">\n')
        of.write('        <link href="//fonts.googleapis.com/css?family=Lobster&subset=latin,latin-ext" rel="stylesheet" type="text/css">\n')
        of.write('    </head>\n')
        of.write('    <body>\n')
        of.write('        <header class="page-header">\n')
        of.write('            <h2>Fama report</h2>\n')
        of.write('        </header>\n')

        
        stats = defaultdict(dict)
        func_stats = autovivify(2,float)
        func_counts = autovivify(2,int)
        func_identity = autovivify(2,float)
        func_hit_counts = autovivify(2,int)
        function_dict = {}

        for parser in parsers:
            stats[parser.end]['reads_total'] = str(parser.project.get_fastq1_readcount(parser.sample))
            read_stats = Counter()
            for read in sorted(parser.reads.keys()):
                read_stats[parser.reads[read].get_status()] += 1
            stats[parser.end]['reads_mapped'] = str(read_stats['function'])

            for read in parser.reads.keys():
                if parser.reads[read].get_status() == 'function':
                    functions = parser.reads[read].get_functions()
                    for function in functions:
                        function_dict[function] = True
                        func_stats[parser.end][function] += functions[function]
                        func_counts[parser.end][function] += 1/len(functions)
                    for hit in parser.reads[read].get_hit_list().get_hits():
                        for function in hit.get_functions():
                            func_identity[parser.end][function] += hit.get_identity()
                            func_hit_counts[parser.end][function] += 1
            for function in func_identity[parser.end]:
                func_identity[parser.end][function] = func_identity[parser.end][function]/func_hit_counts[parser.end][function]

        of.write('\n    <h3>Run info</h3>\n')

        of.write('            <table>\n')
        of.write('                <thead>\n')
        of.write('                            <tr>\n')
        of.write('                                <th>&nbsp;</th>\n')
        for parser in parsers:
            of.write('                                <th>' + parser.end + ' </th>\n')

        of.write('                </tr>\n')
        of.write('            </thead>\n')
        of.write('            <tbody>\n')
        of.write('                    <tr>\n')
        of.write('                          <td>Reads, total</td>\n')
        for parser in parsers:
            of.write('                                <td>' + str(stats[parser.end]['reads_total']) + '</td>\n')
        of.write('                    </tr>\n')
        of.write('                    <tr>\n')
        of.write('                          <td>Reads, mapped</td>\n')
        for parser in parsers:
            of.write('                                <td>' + str(stats[parser.end]['reads_total']) + '</td>\n')
        of.write('                    </tr>\n')
        of.write('            </tbody>\n')
        of.write('        </table>\n')


        of.write('\n    <h3>Functional profile</h3>\n')

        of.write('            <table>\n')
        of.write('                <thead>\n')
        of.write('                            <tr>\n')
        of.write('                                <th>Fucntion</th>\n')
        for parser in parsers:
            of.write('                                <th>Score, ' + parser.end + '</th>\n')
            of.write('                                <th>Read count, ' + parser.end + '</th>\n')
            of.write('                                <th>Aver. % identity' + parser.end + '</th>\n')

        of.write('                </tr>\n')
        of.write('            </thead>\n')
        of.write('            <tbody>\n')
        for function in sorted(function_dict.keys()):
            of.write('                    <tr>\n')
            of.write('                          <td>' + function + '</td>\n')
            for parser in parsers:
                if function in func_stats[parser.end]:
                    of.write('                                <td>' + '{0:.5f}'.format(func_stats[parser.end][function]) + '</td>\n')
                else:
                    of.write('                                <td>N/A</td>\n')
                if function in func_counts[parser.end]:
                    of.write('                                <td>' + '{0:.0f}'.format(func_counts[parser.end][function]) + '</td>\n')
                else:
                    of.write('                                <td>0</td>\n')
                if function in func_identity[parser.end]:
                    of.write('                                <td>' + '{0:.2f}'.format(func_identity[parser.end][function]) + '%</td>\n')
                else:
                    of.write('                                <td>N/A</td>\n')
            of.write('                    </tr>\n')

        of.write('            </tbody>\n')
        of.write('        </table>\n')



        # Write group scores
        of.write('\n        <h3>Function statistics by category</h3>\n')
        func_stats = autovivify(2,float)
        func_counts = autovivify(2,int)
        func_identity = autovivify(2,float)
        func_hit_counts = autovivify(2,int)
        for parser in parsers:
            for read in parser.reads.keys():
                if parser.reads[read].get_status() == 'function':
                    functions = parser.reads[read].get_functions()
                    for function in functions:
                        func_stats[parser.end][parser.ref_data.lookup_function_group(function)] += functions[function]
                        func_counts[parser.end][parser.ref_data.lookup_function_group(function)] += 1/len(functions)
                    for hit in parser.reads[read].get_hit_list().get_hits():
                        for function in hit.get_functions():
                            func_identity[parser.end][parser.ref_data.lookup_function_group(function)] += hit.get_identity()
                            func_hit_counts[parser.end][parser.ref_data.lookup_function_group(function)] += 1
            for function in func_identity[parser.end]:
                func_identity[parser.end][function] = func_identity[parser.end][function]/func_hit_counts[parser.end][function]


        of.write('            <table>\n')
        of.write('                <thead>\n')
        of.write('                            <tr>\n')
        of.write('                                <th>Category</th>\n')
        for parser in parsers:
            of.write('                                <th>Score, ' + parser.end + '</th>\n')
            of.write('                                <th>Read count, ' + parser.end + '</th>\n')
            of.write('                                <th>Aver. % identity' + parser.end + '</th>\n')

        of.write('                </tr>\n')
        of.write('            </thead>\n')
        of.write('            <tbody>\n')
        for function in sorted(function_dict.keys()):
            of.write('                    <tr>\n')
            of.write('                          <td>' + function + '</td>\n')
            for parser in parsers:
                if function in func_stats[parser.end]:
                    of.write('                                <td>' + '{0:.5f}'.format(func_stats[parser.end][function]) + '</td>\n')
                else:
                    of.write('                                <td>N/A</td>\n')
                if function in func_counts[parser.end]:
                    of.write('                                <td>' + '{0:.0f}'.format(func_counts[parser.end][function]) + '</td>\n')
                else:
                    of.write('                                <td>0</td>\n')
                if function in func_identity[parser.end]:
                    of.write('                                <td>' + '{0:.2f}'.format(func_identity[parser.end][function]) + '%</td>\n')
                else:
                    of.write('                                <td>N/A</td>\n')
            of.write('                    </tr>\n')

        of.write('            </tbody>\n')
        of.write('        </table>\n')




        of.write('<p>End of report</p>\n')
        of.write('</body>\n</html>\n')
        of.closed

    
