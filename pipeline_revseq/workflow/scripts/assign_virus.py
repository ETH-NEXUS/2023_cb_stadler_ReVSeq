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
import argparse, pandas.errors
from matplotlib import pyplot

def get_outliers(aggregated, percentile):
    perc = np.percentile(aggregated['rpkm_proportions'], percentile)
    aggregated['outlier'] = np.where(np.logical_and(aggregated['rpkm_proportions'] >= perc, aggregated['rpkm_proportions'] != 0), "*", "")
    aggregated['percentile_threshold'] = str(percentile) + " percentile: " + "{:.5f}".format(perc)
    return [aggregated, perc]


def rpkm(aggregated):
    aggregated['rpkm'] = aggregated['aligned']/(aggregated['length']/1000 * aggregated['aligned'].sum()/1000000)
    aggregated['rpkm_proportions'] = aggregated['rpkm'] / aggregated['rpkm'].sum() * 100
    return aggregated


def fetch_coverage(genome_res):
    try:
        #after the string there is a blank line, so we need to move 2 indexes forward to catch the first genome
        start_index = genome_res.index(">>>>>>> Coverage per contig") + 2
    except:
        sys.exit("Error: cannot find the coverage by contig in the genome_results file. Please check that the correct file is passed and that the format has not changed")
    #end_index = start_index + number_refs - 1
    cov_df = pd.DataFrame([ [item.split('\t')[0],item.split('\t')[3] ] for item in genome_res[start_index:] if item not in [""] ], columns=["id", "coverage"], dtype='float64')
    return cov_df


def get_dp(depth, value, aggregated_stats, ref_table):
    dp_colname = "DP"
    aggregated_stats[dp_colname] = ""
    for i in aggregated_stats.index:
        dp_string = ""
        ref_id = ref_table.loc[ref_table["name"] == i, "id"].values[0]
        for j in range(1, value+1):
            if ref_id in depth[0].unique():
                depth_filtered = depth.loc[depth[0] == ref_id][2].tolist()
                depth_dp = [ v for v in depth_filtered if v >= j]
                dp_value =len(depth_dp)/len(depth_filtered)
                dp_value = '%.2f' % dp_value
            else:
                dp_value = "0.00"
            dp_string = dp_string + str(j) + ":" + dp_value + ";"
        aggregated_stats.loc[aggregated_stats.index == i, dp_colname] = dp_string
    return aggregated_stats


# Script
if __name__ == '__main__':
	# Parse input args
    parser = argparse.ArgumentParser(description='Detects, counts and reports viruses in a sample', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref_table', required=True, type=str, help='the bed file listing all the reference viral genomes')
    parser.add_argument('--idxstats', required=True, type=str, help='the output of samtools idxstats')
    parser.add_argument('--counts', required=True, type=str, help='the count file output by samtools view -c -F 4')
    parser.add_argument('--taxon_table', required=True, type=str, help='the table with all the substrain name/taxon matches')
    parser.add_argument('--depth_table', required=True, type=str, help='the table with the read depth calculated with samtools depth')
    parser.add_argument('--out_prefix', required=True, type=str, help='the prefix for the outputs. Can include a path')
    parser.add_argument('--lookup', required=True, type=str, help='Lookup table in CSV format showing the match between strain_name and substrain_name (with column names). The table will be used to collapse substrains into strains')
    parser.add_argument('--outlier_percentile', required=True, type=float, help='percentile cutoff for outlier detection')
    parser.add_argument('--outlier_percentile_collapsed', required=True, type=float, help='percentile cutoff for outlier detection for collapsed list showing only major strains')
    parser.add_argument('--genome_res', required=True, type=str, help='The genome_results.txt file as output by qualimap bamqc')
    parser.add_argument('--coverage_threshold', required=True, type=int, help='The minimum average coverage allowed for assignment')
    parser.add_argument('--readnum_threshold', required=True, type=int, help='The minimum number of reads allowed for coverage calculation')
    parser.add_argument('--dp_limit', required=True, type=int, help="The highest DP to report. All DPs from 1 to this values will be included in the table")

    args = parser.parse_args()

    pyplot.ioff()

    refs = pd.read_table(args.ref_table, header=None)
    refs = refs.rename(columns={0: "id", 1: "ref_start", 2: "ref_end", 3: "name"})
    idxstats = pd.read_table(args.idxstats, header=None)
    idxstats = idxstats.rename(columns={0: "id", 1: "length", 2: "aligned", 3: "unaligned"})
    taxon = pd.read_csv(args.taxon_table, index_col=0)
    try:
        depth = pd.read_table(args.depth_table, header=None)
    except pandas.errors.EmptyDataError:
        depth = pd.DataFrame(columns=[0, 1, 2])

    with open(args.counts) as f:
        counts = int(f.readline().strip())

    lookup = pd.read_csv(args.lookup, header=0, index_col=1)
    lookup_dict = lookup.to_dict(orient="index")
    for k,v in lookup_dict.items():
        lookup_dict[k] = v['strain_name']
    
    with open(args.genome_res) as f:
        genome_res = f.readlines()

    idxstats = pd.merge(idxstats, refs, on='id', how='inner')

    if len(genome_res) == 0:
        idxstats["coverage"] = 0
    else:
        genome_res = [ line.strip() for line in genome_res ]
        coverage = fetch_coverage(genome_res)
        idxstats = pd.merge(idxstats, coverage, on='id', how='inner')

    aggregated_stats = idxstats.groupby(by='name')[["aligned","length","coverage"]].sum()
    aggregated_stats = rpkm(aggregated_stats)
    all_outlier_info = get_outliers(aggregated_stats, args.outlier_percentile)
    aggregated_stats = all_outlier_info[0]
    outlier_threshold = all_outlier_info[1]
    aggregated_stats = aggregated_stats.rename(columns={"name": "reference_name", "aligned": "aligned_reads", "length": "reference_length"})
    aggregated_stats["readnum_threshold"] = args.readnum_threshold
    if len(genome_res) == 0:
        aggregated_stats["readnum_status"] = "FAILED"
    else:
        aggregated_stats["readnum_status"] = "PASSED"
    aggregated_stats["coverage_threshold"] = args.coverage_threshold
    aggregated_stats["coverage_status"] = "placeholder"
    aggregated_stats.loc[aggregated_stats['coverage'] >= args.coverage_threshold, "coverage_status"] = "PASSED"
    aggregated_stats.loc[aggregated_stats['coverage'] < args.coverage_threshold, "coverage_status"] = "FAILED"
    aggregated_stats = aggregated_stats.merge(taxon, right_index=True, left_index=True, how="left")

    aggregated_stats = get_dp(depth, 20, aggregated_stats, refs)

    aggregated_stats.to_csv(args.out_prefix + "substrain_count_table.tsv", sep="\t", float_format='%.5f')
    pyplot.boxplot(aggregated_stats["rpkm_proportions"])
    for row in aggregated_stats.itertuples():
        rpkm_proportions = row.rpkm_proportions
        if rpkm_proportions > outlier_threshold:
            pyplot.text(1.05, rpkm_proportions, row.Index, ha='right', va='bottom', fontsize='xx-small')
    pyplot.savefig(fname=args.out_prefix + "substrain_proportions_boxplot.pdf", dpi=300, format="pdf")

        # Collapsing substrains in major strains using the lookup table
    aggregated_stats.index = aggregated_stats.index.map(lookup_dict)
    if np.nan in aggregated_stats.index.tolist():
        sys.exit("Error: The strain/substrain lookup did not contain correct entries and produced empty names. Please check the lookup table for syntax errors or missing entries.")
    aggregated_stats = aggregated_stats.groupby(level=0)[["aligned_reads","rpkm", "rpkm_proportions"]].sum()
    all_outlier_info = get_outliers(aggregated_stats, args.outlier_percentile_collapsed)
    aggregated_stats = all_outlier_info[0]
    outlier_threshold = all_outlier_info[1]
    aggregated_stats = aggregated_stats.rename(columns={"name": "reference_name", "aligned": "aligned_reads", "length": "reference_length"})
    aggregated_stats.to_csv(args.out_prefix + "strain_count_table.tsv", sep="\t", float_format='%.5f')
    pyplot.boxplot(aggregated_stats["rpkm_proportions"])
    for row in aggregated_stats.itertuples():
        rpkm_proportions = row.rpkm_proportions
        if rpkm_proportions > outlier_threshold:
            pyplot.text(1.05, rpkm_proportions, row.Index, ha='right', va='bottom', fontsize='xx-small')
    pyplot.savefig(fname=args.out_prefix + "strain_proportions_boxplot.pdf", dpi=300, format="pdf")

