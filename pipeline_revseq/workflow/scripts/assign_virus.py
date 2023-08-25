######
## Script name: assign_virus.py
## Date: July 2023
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq pipeline. Loads the counts and length of
##        the viral sequences computed by the pipeline, calculates RPKM,
##        detects outliers and reports the outliers in a dedicated table
######

import pandas as pd
import numpy as np
import argparse
from matplotlib import pyplot

def get_outliers(aggregated, percentile, coverage):
    perc90 = np.quantile(aggregated['rpkm_proportions'], percentile)
    q1 = np.quantile(aggregated['rpkm_proportions'], 0.25)
    iqr = perc90-q1
    upper_bound = perc90+(1.5*iqr)
    if coverage != 0:
        aggregated['outlier'] = np.where(aggregated['rpkm_proportions'] >= upper_bound, "*", "")
    else:
        aggregated['outlier'] = ''
    aggregated['percentile_threshold'] = str(percentile * 100) + " percentile: " + "{:.4f}".format(upper_bound)
    return aggregated


def rpkm(aggregated, coverage):
    aggregated['rpkm'] = aggregated['aligned']/(aggregated['length']/1000 * aggregated['aligned'].sum()/1000000)
    aggregated['rpkm_proportions'] = aggregated['rpkm'] / aggregated['rpkm'].sum() * 100
    if coverage == 0:
        aggregated['rpkm'] = 0
        aggregated['rpkm_proportions'] = 0
    return aggregated


# Script
if __name__ == '__main__':
	# Parse input args
    parser = argparse.ArgumentParser(description='fetch the primers positions on the reference', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref_table', required=True, type=str, help='the bed file listing all the reference viral genomes')
    parser.add_argument('--idxstats', required=True, type=str, help='the output of samtools idxstats')
    parser.add_argument('--counts', required=True, type=str, help='the count file output by samtools view -c -F 4')
    parser.add_argument('--out_prefix', required=True, type=str, help='the prefix for the outputs. Can include a path')
    parser.add_argument('--lookup', required=True, type=str, help='Lookup table in CSV format showing the match between strain_name and substrain_name (with column names). The table will be used to collapse substrains into strains')
    parser.add_argument('--outlier_percentile', required=True, type=float, help='percentile cutoff for outlier detection')
    parser.add_argument('--outlier_percentile_collapsed', required=True, type=float, help='percentile cutoff for outlier detection for collapsed list showing only major strains')
    parser.add_argument('--genome_res', required=True, type=str, help='The genome_results.txt file as output by qualimap bamqc')
    parser.add_argument('--dp_threshold', required=True, type=int, help='The minimum average coverage that allows for actual assignment')

    args = parser.parse_args()

    pyplot.ioff()

    refs = pd.read_table(args.ref_table, header=None)
    refs = refs.rename(columns={0: "id", 1: "ref_start", 2: "ref_end", 3: "name"})
    idxstats = pd.read_table(args.idxstats, header=None)
    idxstats = idxstats.rename(columns={0: "id", 1: "length", 2: "aligned", 3: "unaligned"})
    with open(args.counts) as f:
        counts = int(f.readline().strip())
            
    lookup = pd.read_csv(args.lookup, header=0, index_col=1).to_dict(orient="index")
    for k,v in lookup.items():
        lookup[k] = v['strain_name']
    
    with open(args.genome_res) as f:
        genome_res = f.readlines()
    if len(genome_res) == 0:
        coverage = 0
    else:
        coverage = [ float(line.split(" = ")[1].split("X")[0]) for line in genome_res if "mean coverageData = " in line ][0]

    idxstats = pd.merge(idxstats,refs, on='id', how='inner')
    aggregated_stats = idxstats.groupby(by='name')[["aligned","length"]].sum()
    aggregated_stats = rpkm(aggregated_stats, coverage)
    aggregated_stats = get_outliers(aggregated_stats, args.outlier_percentile, coverage)
    aggregated_stats = aggregated_stats.rename(columns={"name": "reference_name", "aligned": "aligned_reads", "length": "reference_length"})
    aggregated_stats["coverage"] = coverage
    aggregated_stats["coverage_threshold"] = args.dp_threshold
    if coverage < args.dp_threshold:
        aggregated_stats["qc_status"] = "failed"
    else:
        aggregated_stats["qc_status"] = "passed"
    aggregated_stats.to_csv(args.out_prefix + "substrain_count_table.tsv", sep="\t", float_format='%.5f')
    pyplot.boxplot(aggregated_stats["rpkm_proportions"])
    pyplot.savefig(fname=args.out_prefix + "substrain_proportions_boxplot.pdf", dpi=300, format="pdf")

        # Collapsing substrains in major strains using the lookup table
    aggregated_stats.index = aggregated_stats.index.map(lookup)
    aggregated_stats = aggregated_stats.groupby(level=0)[["aligned_reads","rpkm", "rpkm_proportions"]].sum()
    aggregated_stats = get_outliers(aggregated_stats, args.outlier_percentile_collapsed, coverage)
    aggregated_stats = aggregated_stats.rename(columns={"name": "reference_name", "aligned": "aligned_reads", "length": "reference_length"})
    aggregated_stats["coverage"] = coverage
    aggregated_stats["coverage_threshold"] = args.dp_threshold
    if coverage < args.dp_threshold:
        aggregated_stats["qc_status"] = "failed"
    else:
        aggregated_stats["qc_status"] = "passed"
    aggregated_stats.to_csv(args.out_prefix + "strain_count_table.tsv", sep="\t", float_format='%.5f')
    pyplot.boxplot(aggregated_stats["rpkm_proportions"])
    pyplot.savefig(fname=args.out_prefix + "strain_proportions_boxplot.pdf", dpi=300, format="pdf")

