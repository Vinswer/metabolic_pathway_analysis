#!/usr/bin/env python3
"""
Author: Martijn Prevoo (1285572)
Description: ###
Usage: ###
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage


def hclust(df):
    """Creates a clustermap

    :param df: DataFrame, data to cluster
    """
    # df = df.transpose()
    linkage_data = linkage(df, method='ward', metric='euclidean')
    dendrogram(linkage_data, labels=df.index)

    plt.show()


def main():
    """"""
    import os
    folder = 'input/comparing_counts'  # folder with all count files
    count_files = [f'{folder}/{file}' for file in os.listdir(folder)
                   if file.endswith('_merged.csv')]
    # count_files = [f'{folder}/A-21_2_merged.csv']

    dataframes = []
    for file in count_files:
        df = pd.read_csv(file, index_col=0)
        columns = list(df.columns)
        for c in columns:
            df[c] = df[c] / sum(df[c]) * 1000000  # cpm normalization
        dataframes.append(df.transpose())

    df = pd.concat(dataframes)
    df = df.replace(np.nan, 0)

    hclust(df)


if __name__ == '__main__':
    main()
