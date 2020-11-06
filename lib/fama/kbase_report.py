#!/usr/bin/python3
import os
import numpy as np

from fama.utils.utils import autovivify
from fama.taxonomy.taxonomy_profile import TaxonomyProfile
from fama.output.report import get_function_scores, get_function_taxonomy_scores


def compose_run_info(project, name2ref, tab_index):
    """ Makes content for 'Run info' tab """
    result = ['<div id="tab' + tab_index + '" class="tabcontent">']
    result.append('<p>Reference data set: ' +
                  project.options.get_collection() +
                  '</p>')
    for sample_id in project.list_samples():
        if sample_id in name2ref:
            result.append('<h4>Input object: <a href="https://narrative.kbase.us/#dataview/' +
                          name2ref[sample_id] + '">' +
                          sample_id + '</a></h4>')
        else:
            result.append('<h4>Input object: <a href="https://narrative.kbase.us/#dataview/' +
                          sample_id + '">' +
                          sample_id + '</a></h4>')
        if project.samples[sample_id].is_paired_end:
            result.append('<p>Total number of forward reads: ' +
                          str(project.options.get_fastq1_readcount(sample_id)) +
                          '</p>')
            result.append('<p>Total number of reverse reads: ' +
                          str(project.options.get_fastq2_readcount(sample_id)) +
                          '</p>')
        else:
            result.append('<p>Total number of input sequences: ' +
                          str(project.options.get_fastq1_readcount(sample_id)) +
                          '</p>')
        if project.samples[sample_id].is_paired_end:
            result.append('<p>Number of mapped reads, forward: ' +
                          str(len(project.samples[sample_id].reads['pe1'])) +
                          '</p>')
            result.append('<p>Number of mapped reads, reverse: ' +
                          str(len(project.samples[sample_id].reads['pe2'])) +
                          '</p>')
        else:
            result.append('<p>Number of mapped sequences: ' +
                          str(len(project.samples[sample_id].reads['pe1'])) +
                          '</p>')
        if project.samples[sample_id].is_paired_end:
            result.append('<p>Predicted average insert size: ' +
                          '{0:.0f}'.format(project.get_insert_size(project.samples[sample_id])) +
                          '</p>')
            if (project.samples[sample_id].rpkg_scaling_factor == 0.0 or
                    project.samples[sample_id].rpkg_scaling_factor is None):
                result.append('<p>Average genome size was not calculated.</p>')
            else:
                result.append('<p>Predicted average genome size: ' +
                              '{0:.0f}'.format(
                                                project.options.get_fastq1_basecount(sample_id) *
                                                project.samples[sample_id].rpkg_scaling_factor
                                                ) +
                              '</p>')
    result.append('</div>')
    return '\n'.join(result)


def compose_functional_profile(project, sample_id, tab_index, metric=None):
    """ Makes table with functional profile """
    result = ['<div id="tab' + tab_index + '" class="tabcontent">']
    if metric is None:
        if project.samples[sample_id].is_paired_end:
            if (project.samples[sample_id].rpkg_scaling_factor != 0.0 and
                    project.samples[sample_id].rpkg_scaling_factor is not None):
                metric = 'efpkg'
            elif (project.samples[sample_id].rpkm_scaling_factor != 0.0 and
                    project.samples[sample_id].rpkm_scaling_factor is not None):
                metric = 'efpkm'
            else:
                metric = 'fragmentcount'
        else:
            if (project.samples[sample_id].rpkg_scaling_factor != 0.0 and
                    project.samples[sample_id].rpkg_scaling_factor is not None):
                metric = 'erpkg'
            elif (project.samples[sample_id].rpkm_scaling_factor != 0.0 and
                    project.samples[sample_id].rpkm_scaling_factor is not None):
                metric = 'erpkm'
            else:
                metric = 'readcount'
    scores = get_function_scores(project, sample_id=sample_id, metric=metric)
    result.append('<table><thead><tr>')
    result.append('<th>Function</th><th>' + metric + '</th>')
    if metric not in ('readcount', 'proteincount'):
        result.append('<th>Raw sequence count</th>')
    result.append('<th>Amino acid identity %, average</th>')
    result.append('<th>Description</th>')
    result.append('</thead></tr>')

    for function in sorted(scores.keys()):
        if sample_id in scores[function]:
            result.append('<tr><td>' + function + '</td>')
            if metric in ('readcount', 'fragmentcount', 'proteincount'):
                result.append('<td>' + str(int(scores[function][sample_id][metric])) + '</td>')
            else:
                result.append('<td>' + '{0:.5f}'.format(scores[function][sample_id][metric]) +
                              '</td>'
                              )
            if metric not in ('readcount', 'proteincount'):
                result.append('<td>' +
                              '{0:.0f}'.format(scores[function][sample_id]['count']) + '</td>'
                              )
            result.append('<td>' + '{0:.2f}'.format(
                scores[function][sample_id]['identity'] / scores[function][sample_id]['hit_count']
                ) + '</td>')
            result.append('<td>' + project.ref_data.lookup_function_name(function)
                          + '</td></tr>')
    result.append('</table>')
    result.append('</div>')
    return '\n'.join(result)


def compose_function_groups(project, sample_id, tab_index, metric=None):
    """ Makes table of functional groups """
    result = ['<div id="tab' + tab_index + '" class="tabcontent">']
    if metric is None:
        if project.samples[sample_id].is_paired_end:
            if (project.samples[sample_id].rpkg_scaling_factor != 0.0 and
                    project.samples[sample_id].rpkg_scaling_factor is not None):
                metric = 'efpkg'
            elif (project.samples[sample_id].rpkm_scaling_factor != 0.0 and
                    project.samples[sample_id].rpkm_scaling_factor is not None):
                metric = 'efpkm'
            else:
                metric = 'fragmentcount'
        else:
            if (project.samples[sample_id].rpkg_scaling_factor != 0.0 and
                    project.samples[sample_id].rpkg_scaling_factor is not None):
                metric = 'erpkg'
            elif (project.samples[sample_id].rpkm_scaling_factor != 0.0 and
                    project.samples[sample_id].rpkm_scaling_factor is not None):
                metric = 'erpkm'
            else:
                metric = 'readcount'
    scores = get_function_scores(project, sample_id=sample_id, metric=metric)
    categories = set()
    for function in scores:
        categories.add(project.ref_data.lookup_function_group(function))

    result.append('<table><thead><tr>')
    result.append('<th>Function category</th><th>' + metric + '</th>')
    if metric not in ('readcount', 'proteincount'):
        result.append('<th>Raw sequence count</th>')
    result.append('<th>Amino acid identity %, average</th>')
    result.append('</thead></tr>')

    for category in sorted(list(categories)):
        category_data = autovivify(2, float)
        for function in scores:
            for sample in scores[function]:
                if project.ref_data.lookup_function_group(function) != category:
                    continue
                category_data[sample][metric] += scores[function][sample][metric]
                category_data[sample]['count'] += scores[function][sample]['count']
                category_data[sample]['identity'] += scores[function][sample]['identity']
                category_data[sample]['hit_count'] += scores[function][sample]['hit_count']

        result.append('<tr><td>' + category + '</td>')
        if metric in ('readcount', 'fragmentcount', 'proteincount'):
            result.append('<td>' + str(int(category_data[sample_id][metric])) + '</td>')
        else:
            result.append('<td>' + '{0:.5f}'.format(category_data[sample_id][metric]) + '</td>')
        if metric not in ('readcount', 'proteincount'):
            result.append('<td>' + '{0:.0f}'.format(category_data[sample_id]['count']) + '</td>')
        result.append('<td>' + '{0:.2f}'.format(
            category_data[sample_id]['identity'] / category_data[sample_id]['hit_count']
            ) + '</td></tr>')

    result.append('</table>')
    result.append('</div>')
    return '\n'.join(result)


def compose_taxonomy_profile(project, sample_id, tab_index, metric=None):
    """Makes taxonomy profile """
    if metric is None:
        if project.samples[sample_id].is_paired_end:
            if (project.samples[sample_id].rpkg_scaling_factor != 0.0 and
                    project.samples[sample_id].rpkg_scaling_factor is not None):
                metric = 'efpkg'
            elif (project.samples[sample_id].rpkm_scaling_factor != 0.0 and
                    project.samples[sample_id].rpkm_scaling_factor is not None):
                metric = 'efpkm'
            else:
                metric = 'fragmentcount'
        else:
            if (project.samples[sample_id].rpkg_scaling_factor != 0.0 and
                    project.samples[sample_id].rpkg_scaling_factor is not None):
                metric = 'erpkg'
            elif (project.samples[sample_id].rpkm_scaling_factor != 0.0 and
                    project.samples[sample_id].rpkm_scaling_factor is not None):
                metric = 'erpkm'
            else:
                metric = 'readcount'
    scores = get_function_taxonomy_scores(project, sample_id=sample_id, metric=metric)

    sample_scores = autovivify(3, float)
    for taxonomy_id in scores.keys():
        for function_id in scores[taxonomy_id].keys():
            if sample_id in scores[taxonomy_id][function_id]:
                for key, val in scores[taxonomy_id][function_id][sample_id].items():
                    sample_scores[taxonomy_id][function_id][key] = val

    tax_profile = TaxonomyProfile()
    tax_profile.make_function_taxonomy_profile(project.taxonomy_data, sample_scores)
    taxonomy_df = tax_profile.convert_profile_into_score_df(metric=metric)
    taxonomy_df.replace(0.0, np.nan, inplace=True)
    result = '<div id="tab' + tab_index + '" class="tabcontent">\n'
    if metric in ('readcount', 'fragmentcount', 'proteincount'):
        result += taxonomy_df.to_html(na_rep="", float_format = '%.0f')
    else:
        result += taxonomy_df.to_html(na_rep="")
    result += '\n</div>\n'
    return result


def compose_protein_list(project, name2ref, tab_index):
    """ Makes table of proteins """
    result = ['<div id="tab' + tab_index + '" class="tabcontent">']
    result.append('<table><thead><tr>')
    result.append('<th>Input object</th>')
    result.append('<th>Protein</th>')
    result.append('<th>Function</th>')
    result.append('<th>Description</th>')
    result.append('<th>Amino acid identity %</th>')
    result.append('<th>Taxonomy</th>')
    result.append('</thead></tr>')

    for sample_id in project.list_samples():
        if 'pe1' in project.samples[sample_id].reads:

            for protein_id in sorted(project.samples[sample_id].reads['pe1'].keys()):
                protein = project.samples[sample_id].reads['pe1'][protein_id]
                if protein.status == 'function':
                    result.append('<tr><td>' + sample_id + '</td>')
                    result.append('<td>' + protein_id + '</td>')
                    fama_identity = sum(
                        [x.identity for x in protein.hit_list.hits]
                        ) / len(protein.hit_list.hits)
                    function = ','.join(sorted(protein.functions.keys()))
                    result.append('<td>' + function + '</td>')
                    description = '|'.join(
                        sorted([
                            project.ref_data.lookup_function_name(f) for f
                            in protein.functions.keys()
                            ])
                        )
                    result.append('<td>' + description + '</td>')
                    result.append('<td>' + '{0:.1f}'.format(fama_identity) + '</td>')
                    result.append('<td>' + project.taxonomy_data.data[protein.taxonomy]['name'] +
                                  '</td></tr>')
    result.append('</table>')
    result.append('</div>')
    return '\n'.join(result)


def generate_html_report(outfile, project, name2ref):
    """ Generates HTML report """
    html_template = os.path.join(os.path.dirname(__file__), 'kbase_report.template')
    report_tabs = []
    with open(outfile, 'w') as of:
        with open(html_template, 'r') as infile:
            for line in infile:
                line = line.rstrip('\n\r')
                if line == '<\InsertButtons>':
                    tab_counter = 1
                    of.write('<button class="tablinks" onclick="openTab(event, \'tab' +
                             str(tab_counter) + '\')" id="defaultOpen">Run info</button>')
                    report_tabs.append(compose_run_info(project, name2ref, str(tab_counter)))
                    for sample_id in project.list_samples():
                        tab_counter += 1
                        of.write('<button class="tablinks" onclick="openTab(event, \'tab' +
                                 str(tab_counter) + '\')" id="defaultOpen">' +
                                 sample_id + '<br>Functional profile</button>')
                        report_tabs.append(compose_functional_profile(project,
                                                                      sample_id,
                                                                      str(tab_counter)))
                        tab_counter += 1
                        of.write('<button class="tablinks" onclick="openTab(event, \'tab' +
                                 str(tab_counter) + '\')" id="defaultOpen">' +
                                 sample_id + '<br>Functional groups</button>')
                        report_tabs.append(compose_function_groups(project,
                                                                   sample_id,
                                                                   str(tab_counter)))
                        tab_counter += 1
                        of.write('<button class="tablinks" onclick="openTab(event, \'tab' +
                                 str(tab_counter) + '\')" id="defaultOpen">' +
                                 sample_id + '<br>Taxonomy profile</button>')
                        report_tabs.append(compose_taxonomy_profile(project,
                                                                    sample_id,
                                                                    str(tab_counter)))
                elif line == '<\InsertTabs>':
                    of.write('\n'.join(report_tabs))
                else:
                    of.write(line + '\n')


def generate_protein_html_report(outfile, project, name2ref):
    """ Generates HTML report for protein project """
    html_template = os.path.join(os.path.dirname(__file__), 'kbase_protein_report.template')
    report_tabs = []
    with open(outfile, 'w') as of:
        with open(html_template, 'r') as infile:
            for line in infile:
                line = line.rstrip('\n\r')
                if line == '<\InsertButtons>':
                    tab_counter = 1
                    of.write('<button class="tablinks" onclick="openTab(event, \'tab' +
                             str(tab_counter) + '\')" id="defaultOpen">Run info</button>')
                    report_tabs.append(compose_run_info(project,
                                                        name2ref,
                                                        str(tab_counter)))
                    tab_counter += 1
                    of.write('<button class="tablinks" onclick="openTab(event, \'tab' +
                             str(tab_counter) + '\')" id="defaultOpen">Protein list</button>')
                    report_tabs.append(compose_protein_list(project,
                                                            name2ref,
                                                            str(tab_counter)))

                    for sample_id in project.list_samples():
                        tab_counter += 1
                        of.write('<button class="tablinks" onclick="openTab(event, \'tab' +
                                 str(tab_counter) + '\')" id="defaultOpen">' +
                                 sample_id + '<br>Functional profile</button>')
                        report_tabs.append(compose_functional_profile(project,
                                                                      sample_id,
                                                                      str(tab_counter),
                                                                      metric='proteincount'))
                        tab_counter += 1
                        of.write('<button class="tablinks" onclick="openTab(event, \'tab' +
                                 str(tab_counter) + '\')" id="defaultOpen">' +
                                 sample_id + '<br>Functional groups</button>')
                        report_tabs.append(compose_function_groups(project,
                                                                   sample_id,
                                                                   str(tab_counter),
                                                                   metric='proteincount'))
                        tab_counter += 1
                        of.write('<button class="tablinks" onclick="openTab(event, \'tab' +
                                 str(tab_counter) + '\')" id="defaultOpen">' +
                                 sample_id + '<br>Taxonomy profile</button>')
                        report_tabs.append(compose_taxonomy_profile(project,
                                                                    sample_id,
                                                                    str(tab_counter),
                                                                    metric='proteincount'))
                elif line == '<\InsertTabs>':
                    of.write('\n'.join(report_tabs))
                else:
                    of.write(line + '\n')
