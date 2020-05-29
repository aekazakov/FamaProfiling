#!/usr/bin/python3
import os
import numpy as np

from fama.utils.utils import autovivify
from fama.taxonomy.taxonomy_profile import TaxonomyProfile
from fama.output.report import get_function_scores, get_function_taxonomy_scores


def compose_run_info(project):
    result = []
    sample_id = project.list_samples()[0]
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
    result.append('<p>Reference data set: ' +
                  project.options.get_collection(sample_id) +
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

    return '\n'.join(result)


def compose_functional_profile(project):
    result = []
    sample_id = project.list_samples()[0]
    metric = None
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
    if metric != 'readcount':
        result.append('<th>Raw sequence count</th>')
    result.append('<th>Amino acid identity %, average</th>')
    result.append('<th>Description</th>')
    result.append('</thead></tr>')

    for function in scores:
        if sample_id in scores[function]:
            result.append('<tr><td>' + function + '</td>')
            result.append('<td>' + '{0:.5f}'.format(scores[function][sample_id][metric]) + '</td>')
            if metric != 'readcount':
                result.append('<td>' +
                              '{0:.0f}'.format(scores[function][sample_id]['count']) + '</td>'
                              )
            result.append('<td>' + '{0:.2f}'.format(
                scores[function][sample_id]['identity'] / scores[function][sample_id]['hit_count']
                ) + '</td>')
            result.append('<td>' + project.ref_data.lookup_function_name(function)
                          + '</td></tr>')
    result.append('</table>')
    return '\n'.join(result)


def compose_function_groups(project):
    result = []
    sample_id = project.list_samples()[0]
    metric = None
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
    if metric != 'readcount':
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
        result.append('<td>' + '{0:.5f}'.format(category_data[sample_id][metric]) + '</td>')
        if metric != 'readcount':
            result.append('<td>' + '{0:.0f}'.format(category_data[sample_id]['count']) + '</td>')
        result.append('<td>' + '{0:.2f}'.format(
            category_data[sample_id]['identity'] / category_data[sample_id]['hit_count']
            ) + '</td></tr>')

    result.append('</table>')
    return '\n'.join(result)


def compose_taxonomy_profile(project):
    sample_id = project.list_samples()[0]
    metric = None
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
    return taxonomy_df.to_html(na_rep="")  # , float_format=lambda x: '%.2f' % x)


def compose_protein_list(project):
    result = []
    sample_id = project.list_samples()[0]
    metric = None
    if 'pe1' in project.samples[sample_id].reads:
        result.append('<table><thead><tr>')
        result.append('<th>Protein</th>')
        result.append('<th>Function</th>')
        result.append('<th>Description</th>')
        result.append('<th>Amino acid identity %</th>')
        result.append('<th>Taxonomy</th>')
        result.append('</thead></tr>')

        for protein_id in sorted(project.samples[sample_id].reads['pe1'].keys()):
            protein = project.samples[sample_id].reads['pe1'][protein_id]
            if protein.status == 'function':
                result.append('<tr><td>' + protein_id + '</td>')
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
                result.append('<td>' + project.taxonomy_data.data[protein.taxonomy]['name'] + '</td></tr>')
        result.append('</table>')
    return '\n'.join(result)


def generate_html_report(outfile, project):
    """ Generates HTML report """
    html_template = os.path.join(os.path.dirname(__file__), 'kbase_report.template')
    with open(outfile, 'w') as of:
        with open(html_template, 'r') as infile:
            for line in infile:
                line = line.rstrip('\n\r')
                if line == '<\RunInfo>':
                    of.write(compose_run_info(project))
                elif line == '<\FunctionalProfile>':
                    of.write(compose_functional_profile(project))
                elif line == '<\FunctionGroups>':
                    of.write(compose_function_groups(project))
                elif line == '<\TaxonomyProfile>':
                    taxonomy_profile = compose_taxonomy_profile(project)
                    if taxonomy_profile:
                        of.write(taxonomy_profile)
                    else:
                        of.write('<p>No taxonomy data</p>/n')
                else:
                    of.write(line + '\n')


def generate_protein_html_report(outfile, project):
    """ Generates HTML report for protein project """
    html_template = os.path.join(os.path.dirname(__file__), 'kbase_protein_report.template')
    with open(outfile, 'w') as of:
        with open(html_template, 'r') as infile:
            for line in infile:
                line = line.rstrip('\n\r')
                if line == '<\RunInfo>':
                    of.write(compose_run_info(project))
                elif line == '<\FunctionalProfile>':
                    of.write(compose_functional_profile(project))
                elif line == '<\FunctionGroups>':
                    of.write(compose_function_groups(project))
                elif line == '<\TaxonomyProfile>':
                    taxonomy_profile = compose_taxonomy_profile(project)
                    if taxonomy_profile:
                        of.write(taxonomy_profile)
                    else:
                        of.write('<p>No taxonomy data</p>/n')
                elif line == '<\ProteinList>':
                    protein_list = compose_protein_list(project)
                    if protein_list:
                        of.write(protein_list)
                    else:
                        of.write('<p>No proteins mapped</p>/n')
                else:
                    of.write(line + '\n')
