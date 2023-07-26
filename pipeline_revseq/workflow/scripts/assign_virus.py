import pandas as pd
import numpy as np
import argparse

def get_highest_alignment(aggregated, coinfection):
    maxalign = max(aggregated['aligned'])
    to_consider = aggregated [ aggregated['aligned'] >= (maxalign - (maxalign * coinfection / 100 )) ]
    return to_consider

def get_outliers(aggregated, percentile):
    perc90 = np.quantile(aggregated['rpkm_proportions'], percentile)
    q1 = np.quantile(aggregated['rpkm_proportions'], 0.25)
    iqr = perc90-q1
    upper_bound = perc90+(1.5*iqr)
    aggregated['outlier'] = np.where(aggregated['rpkm_proportions'] >= upper_bound, "*", "")
    return aggregated


# Script
if __name__ == '__main__':
	# Parse input args
    parser = argparse.ArgumentParser(description='fetch the primers positions on the reference', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref_table', required=True, type=str, help='the bed file listing all the reference viral genomes')
    parser.add_argument('--idxstats', required=True, type=str, help='the output of samtools idxstats')
    parser.add_argument('--counts', required=True, type=str, help='the count file output by samtools view -c -F 4')
    parser.add_argument('--out_prefix', required=True, type=str, help='the prefix for the outputs. Can include a path')
    parser.add_argument('-c','--collapse', nargs='+', help='generic name of viruses in the reference table to merge together in case of multiple strains. Can be called multiple times', required=True)
    parser.add_argument('--outlier_percentile', required=True, type=float, help='percentile cutoff for outlier detection')
    parser.add_argument('--outlier_percentile_collapsed', required=True, type=float, help='percentile cutoff for outlier detection for collapsed list showing only major strains')
    parser.add_argument('--aggregate_rhino_entero', required=True, type=int, help='1 if the adenovirus and enterovirus entries need to be merged. Works only if --collapse includes both adenovirus and enterovirus as separate voices')

    refs = pd.read_table(args.ref_table, header=None)
    refs = refs.rename(columns={0: "id", 1: "ref_start", 2: "ref_end", 3: "name"})
    idxstats = pd.read_table(args.idxstats, header=None)
    idxstats = idxstats.rename(columns={0: "id", 1: "length", 2: "aligned", 3: "unaligned"})
    with open(args.counts) as f:
        counts = int(f.readLine().strip())

    idxstats = pd.merge(idxstats,refs, on='id', how='inner')
    aggregated_stats = idxstats.groupby(by='name')[["aligned","length"]].sum()
    aggregated_stats['rpkm'] = aggregated_stats['aligned']/(aggregated_stats['length']/1000 * aggregated_stats['aligned'].sum()/1000000)
    aggregated_stats['rpkm_proportions'] = aggregated_stats['rpkm'] / aggregated_stats['rpkm'].sum() * 100
    aggregated_stats['normcounts'] = aggregated_stats['aligned']/aggregated_stats['length']
    aggregated_stats = get_outliers(aggregated_stats, args.outlier_percentile)
    aggregated_stats.to_csv(args.out_prefix + "count_table.tsv")


    temp = pd.DataFrame({'aligned': pd.Series(dtype='int'), 'length': pd.Series(dtype='int')})
    for virus in args.collapse:
        all_virus = aggregated_stats[[virus in s for s in aggregated_stats.index]]
        temp.loc[virus] = [all_virus['aligned'].sum(), all_virus['length'].sum()]
        aggregated_stats = aggregated_stats.drop(aggregated_stats[aggregated_stats.index.str.contains(virus)].index)

    temp = pd.concat([temp, aggregated_stats[['aligned', 'length']]])
    
    if args.aggregate_rhino_entero == 1:
        all_rhino_entero = pd.DataFrame({'aligned': pd.Series(dtype='int'), 'length': pd.Series(dtype='int')})
        rhino_entero_counts = int(temp[['rhinovirus' in s for s in temp.index]]['aligned'].iloc[0]) + int(temp[['enterovirus' in s for s in temp.index]]['aligned'].iloc[0])
        rhino_entero_length = int(temp[['rhinovirus' in s for s in temp.index]]['length'].iloc[0]) + int(temp[['enterovirus' in s for s in temp.index]]['length'].iloc[0])
        all_rhino_entero.loc["rhinovirus/enterovirus"] = [rhino_entero_counts, rhino_entero_length]
        temp = temp.drop(temp[temp.index.str.contains("rhinovirus")].index)
        temp = temp.drop(temp[temp.index.str.contains("enterovirus")].index)
        temp = pd.concat[temp, all_rhino_entero]
    
    major = temp

    major['rpkm'] = major['aligned']/(major['length']/1000 * major['aligned'].sum()/1000000)
    major['rpkm_proportions'] = major['rpkm'] / major['rpkm'].sum() * 100
    major['normcounts'] = major['aligned']/major['length']
    major = get_outliers(major, args.outlier_percentile_collapsed)
    major.to_csv(args.out_prefix + "major_count_table.csv")

