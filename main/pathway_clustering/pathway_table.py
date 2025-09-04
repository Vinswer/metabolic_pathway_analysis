#!/usr/bin/env python3

from sys import argv
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt


def parse_pathway(pathway_file):
    """Parse pathway file

    :param pathway_file: str, path to pathway file
    :return: tuple(dict, dict), [0] key: str, pathway ID
                                    value: list(str), genes
                                [1] key: str, pathway ID
                                    value: str, pathway description

    """
    pathways = {}
    descriptions = {}
    with open(pathway_file, 'r') as fo:
        for line in fo:
            line = line.strip().replace('"', '').split('\t')
            if len(line) < 3 or not line[2]:
                continue
            p_id, p_desc, p_gens = line
            pathways[p_id] = p_gens.split(';')
            descriptions[p_id] = p_desc.split(' - ')[0]
    return pathways, descriptions


def parse_deg(deg_file):
    """Parse differential expression results

    :param deg_file: str, path to differential expression results
    :return: dict, key: str, gene
                   value: dict{'logFC': float, 'PValue': float}
    """
    data = {}
    with open(deg_file, 'r') as fo:
        _ = fo.readline()  # remove header
        for line in fo:
            gene, logfc, _, _, pvalue, fdr = line.strip().split('\t')
            gene = gene.replace('"', '').split('|')[0]
            data[gene] = {'logFC': float(logfc), 'PValue': float(fdr)}
    return data


def logfc_table(deg_dict, pathway_dict):
    """Creates a table with each column being the DE file name,
    rows the pathway ID's and values the absolute logFC

    :param deg_dict: dict, key: str, file name
                           value: dict, parse_deg output
    :param pathway_dict: dict, key: str, pathway ID
                               value: list(str), genes
    :return: DataFrame(floats), columns: DE file names
                                rows: pathway ID's
                                values: absolute logFC
    """
    table = pd.DataFrame(0.0, index=pathway_dict, columns=deg_dict)

    for path, in_genes in pathway_dict.items():
        in_genes = set(in_genes)
        for file, deg in deg_dict.items():
            out_genes = set(deg) - in_genes

            in_pvalue = []
            abs_logfc = 0
            out_pvalue = [deg[g]['PValue'] for g in out_genes if g in deg]

            for g in in_genes:
                if g in deg:
                    in_pvalue.append(deg[g]['PValue'])
                    abs_logfc += abs(deg[g]['logFC'])

            pvalue = stats.ranksums(in_pvalue, out_pvalue,
                                    alternative='less')[1]

            if pvalue < .05:
                table[file][path] = abs_logfc
    return table


def clustermap(df, out, size, color_dict, color, vmax: int = 100):
    """Creates a clustermap

    :param df: DataFrame, data to cluster
    :param out: str, path to save the clustermap
    :param size: tuple(int, int), size of the clustermap
    :param color: str, cmap color
    :param color_dict: dict, key; str, label
                             value; str, color
    :param vmax: int, maximum value of color scale
    """
    sns.set(font_scale=2)
    hm = sns.clustermap(df, yticklabels=True, cmap=color, vmax=vmax,
                        method='ward', metric='euclidean',
                        dendrogram_ratio=0.1, figsize=size)

    x_labels = hm.ax_heatmap.get_xticklabels()
    labels = [label.get_text() for label in x_labels]
    for i, label in enumerate(labels):
        label = label.split('_')[0]
        if label in color_dict:
            x_labels[i].set_backgroundcolor(color_dict[label])

    plt.savefig(out, format='svg')
    plt.close()


def name_from_file(path):
    """Gets the Leaf name form a file path

    :param path: str, file path
    :return: str, Leaf name or filename excluding _EdgeR.tsv
    """
    name_dict = {'A-1': 'Axenic',
                 'A-5': 'Leaf069',
                 'A-6': 'Leaf137',
                 'A-15': 'Leaf186',
                 'A-20': 'Leaf203',
                 'A-18': 'Leaf122',
                 'A-35': 'Leaf257',
                 'A-36': 'Leaf357',
                 'A-29': 'Leaf386',
                 'A-9': 'Leaf405',
                 'A-14': 'Leaf082',
                 'A-22': 'Leaf041',
                 'A-2': 'Leaf076',
                 'A-8': 'Leaf177',
                 'A-19': 'Leaf416',
                 'A-39': 'Leaf220',
                 'A-7': 'Leaf049',
                 'A-3': 'Leaf130',
                 'A-24': 'Leaf048',
                 'A-30': 'Leaf051',
                 'A-40': 'Leaf131',
                 'A-4': 'Leaf289',
                 'A-16': 'Leaf085',
                 'A-17': 'Leaf117',
                 'A-21': 'Leaf288',
                 'A-23': 'Leaf015',
                 'A-25': 'Leaf434',
                 'A-26': 'Leaf068',
                 'A-27': 'Leaf311',
                 'A-28': 'Leaf384',
                 'A-31': 'Leaf021',
                 'A-32': 'Leaf028',
                 'A-33': 'Leaf032',
                 'A-34': 'Leaf230',
                 'A-37': 'Fr1',
                 'A-13': 'Leaf053',
                 'A-12': 'Leaf061',
                 'A-11': 'Leaf261',
                 'A-38': 'Leaf070',
                 'A-10': 'Leaf154'}
    name = path.split('/')[-1].replace('_EdgeR.tsv', '')
    if name in name_dict:
        name = name_dict[name]
    return name


def main():
    """Main Function of script"""
    # input files
    pathway_file = argv[1]
    deg_files = argv[2:]

    pathway_dict, descriptions = parse_pathway(pathway_file)
    deg_dict = {name_from_file(file): parse_deg(file) for file in deg_files}

    table = logfc_table(deg_dict, pathway_dict)

    # remove row with sum < 1
    table = table[table.sum(axis='columns') > 1]

    # add descriptions to row
    desc_list = []
    for i, path in enumerate(table.iterrows()):
        desc_list.append(descriptions[path[0]])
    table.insert(0, '', desc_list)

    # print to file
    out_name = f'output/pathway'
    table.to_csv(f'{out_name}_table.tsv', sep='\t')

    table = pd.read_csv(f'{out_name}_table.tsv', sep='\t', index_col=1)
    table = table.drop(table.columns[0], axis='columns')  # remove name column

    color_dict = {'Leaf289': '#FFCFCB',
                  'Leaf137': '#FFCFCB',
                  'Leaf069': '#FFCFCB',
                  'Leaf154': '#FFCFCB',
                  'Leaf261': '#FFCFCB',
                  'Leaf186': '#FFCFCB',
                  'Leaf203': '#FFCFCB',
                  'Leaf288': '#FFCFCB',
                  'Leaf117': '#A6E7C7',
                  'Leaf122': '#A6E7C7',
                  'Leaf085': '#A6E7C7',
                  'Leaf311': '#A6E7C7',
                  'Leaf384': '#A6E7C7',
                  'Leaf386': '#A6E7C7',
                  'Leaf068': '#A6E7C7',
                  'Leaf021': '#A6E7C7',
                  'Leaf230': '#A6E7C7',
                  'Leaf257': '#A6E7C7',
                  'Leaf028': '#A6E7C7',
                  'Leaf032': '#A6E7C7',
                  'Leaf357': '#A6E7C7',
                  'Fr1': '#A6E7C7',
                  'Leaf405': '#A6E1FF',
                  'Leaf082': '#A6E1FF',
                  'Leaf041': '#A6E1FF',
                  'Leaf076': '#DFF1CB',
                  'Leaf177': '#DFF1CB',
                  'Leaf061': '#DFF1CB',
                  'Leaf416': '#DFF1CB',
                  'Leaf220': '#DFF1CB',
                  'Leaf049': '#FFECA7',
                  'Leaf130': '#CAD8BF',
                  'Leaf053': '#CAD8BF',
                  'Leaf015': '#CAD8BF',
                  'Leaf434': '#CAD8BF',
                  'Leaf048': '#CAD8BF',
                  'Leaf051': '#CAD8BF',
                  'Leaf070': '#CAD8BF',
                  'Leaf131': '#CAD8BF'}

    clustermap(table, f'{out_name}_heatmap.svg', (30, 30), color_dict, 'RdPu')

    # col = {'Actinobacteria': 'Reds',
    #        'Alphaproteobacteria': 'BuGn',
    #        'Betaproteobacteria': 'YlGn',
    #        'Gammaproteobacteria': 'Greens',
    #        'Bacteroidetes': 'Blues',
    #        'Firmicutes': 'Yellows'}


if __name__ == '__main__':
    main()
