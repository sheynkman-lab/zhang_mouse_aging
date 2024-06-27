#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from Bio.Seq import translate
from collections import defaultdict
import argparse
import logging
import os

# Function to get unique accession sequences
def get_accession_seqs(seqs):
    logging.info('Getting accession sequences...')
    pb_seqs = defaultdict()  # pb_acc -> transcript_seq
    redundant_accs = []
    for entry in seqs:
        seq = str(entry.seq)
        pb_acc = entry.id.split('|')[0]
        # Skip duplicates
        if pb_acc in pb_seqs:
            redundant_accs.append(pb_acc)
        else:
            pb_seqs[pb_acc] = seq
    return pb_seqs, redundant_accs

# Function to combine by sequence
def combine_by_sequence(orfs, pb_seqs):
    logging.info('Combining by sequence...')
    orfs = orfs[['pb_acc', 'orf_start', 'orf_end', 'orf_len']]
    # Extract, translate, and aggregate protein sequences
    pb_pseqs = defaultdict(lambda: list())  # protein_seq -> list of pb acc
    for index, row in orfs.iterrows():
        pb_acc, start, end, olen = list(row)
        seq = pb_seqs[pb_acc]
        orf_seq = seq[start-1:end]
        prot_seq = translate(orf_seq, to_stop=True)
        pb_pseqs[prot_seq].append(pb_acc)
    return pb_pseqs

# Function to order IDs numerically
def order_enst_acc_numerically(accs):
    logging.info('Ordering ENST Accession Numerically...')
    accs_numerical = []
    for acc in accs:
        try:
            # Split by the dot to separate the gene and isoform parts
            gene_part, iso_idx = acc.split('.')
            # Extract the numerical part from the gene part
            gene_idx = int(gene_part[7:])  # Adjust slicing based on your specific format
            iso_idx = int(iso_idx)
            accs_numerical.append([gene_idx, iso_idx, acc])
        except ValueError as e:
            logging.error(f"Error processing accession ID '{acc}': {e}")
            continue  # Skip this accession ID if it doesn't match the expected format

    # Sort by gene_idx first, then by iso_idx
    accs_numerical.sort(key=lambda x: (x[0], x[1]))

    # Extract the sorted accession strings
    num_sorted_accs = [acc[2] for acc in accs_numerical]
    return num_sorted_accs


# Function to aggregate results
def aggregate_results(pacbio, orfs):
    logging.info('Aggregating Results')
    def get_total(accessions, orf_dict):
        total = 0
        for acc in accessions:
            total += orf_dict.get(acc, 0)
        return total
    
    pacbio['accessions'] = pacbio['pb_accs'].str.split('|')
    pacbio['base_acc'] = pacbio['accessions'].apply(lambda x: x[0])
    
    fl_dict = pd.Series(orfs.FL.values, index=orfs.pb_acc).to_dict()
    cpm_dict = pd.Series(orfs.CPM.values, index=orfs.pb_acc).to_dict()
    
    pacbio['FL'] = pacbio['accessions'].apply(lambda accs: get_total(accs, fl_dict))
    pacbio['CPM'] = pacbio['accessions'].apply(lambda accs: get_total(accs, cpm_dict))
    return pacbio

# Function to filter ORF scores
def filter_orf_scores(orfs, cutoff):
    logging.info(f'Filtering ORF coding_score to be above {cutoff}')
    orfs = orfs[orfs['coding_score'] >= cutoff]
    return orfs

# Function to convert string to boolean
def string_to_boolean(string):
    if isinstance(string, bool):
        return string
    if string.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif string.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# Main function
def main():
    parser = argparse.ArgumentParser(description='Process ORF related file locations')
    parser.add_argument('--name', action='store', dest='name', help='sample name')
    parser.add_argument('-io', '--orfs', action='store', dest='orfs', help='ORF coordinate input file location')
    parser.add_argument('-if', '--pb_fasta', action='store', dest='pb_fasta', help='PacBio fasta sequence input file location')
    parser.add_argument('-cut', '--coding_score_cutoff', dest='cutoff', type=float, default=0.0, help='CPAT coding score cutoff. Remove all below')
    results = parser.parse_args()
    name = results.name

    logging.info('Reading Fasta File...')
    seqs = SeqIO.parse(open(results.pb_fasta), 'fasta')
    logging.info('Reading ORFS...')
    # Specify dtype to avoid mixed types warning
    dtype = {'pb_acc': str, 'coding_score': float, 'orf_score': float, 'orf_calling_confidence': str, 'upstream_atgs': str, 'gene': str}
    orfs = pd.read_table(results.orfs, dtype=dtype)
    
    # Filter ORFS based on score and whether protein coding
    orfs = filter_orf_scores(orfs, results.cutoff)
    # Only keep orfs that have a stop codon
    orfs = orfs.query('has_stop_codon')
    
    pb_seqs, redundant_accs = get_accession_seqs(seqs)
    pb_pseqs = combine_by_sequence(orfs, pb_seqs)
    
    # Save combined results
    with open(f'{name}_combined.tsv', 'w') as ofile, open(f'{name}_combined.fasta', 'w') as ofile2:
        ofile.write('protein_sequence\tpb_accs\n')
        for seq, accs in pb_pseqs.items():
            accs_sorted = order_enst_acc_numerically(accs)  # Use new function for ENST IDs
            accs_str = '|'.join(accs_sorted)
            ofile.write(seq + '\t' + accs_str + '\n')
            ofile2.write('>' + accs_str + '\n' + seq + '\n')
    
    pacbio = pd.read_csv(f'{name}_combined.tsv', sep='\t')
    seqs = SeqIO.parse(open(f'{name}_combined.fasta'), 'fasta')
    os.remove(f'{name}_combined.tsv')
    os.remove(f'{name}_combined.fasta')
    
    pacbio = aggregate_results(pacbio, orfs)
    if 'gene' in orfs.columns:
        orfs = orfs[['pb_acc', 'coding_score', 'orf_score', 'orf_calling_confidence', 'upstream_atgs', 'gene']]
    else:
        logging.warning("'gene' column not found in ORFs data. Proceeding without it.")
        orfs = orfs[['pb_acc', 'coding_score', 'orf_score', 'orf_calling_confidence', 'upstream_atgs']]
    
    pacbio = pd.merge(pacbio, orfs, how='inner', left_on='base_acc', right_on='pb_acc')
    
    logging.info('Writing refined database fasta results...')
    if 'gene' in orfs.columns:
        pb_gene = pd.Series(orfs.gene.values, index=orfs.pb_acc).to_dict()
    base_map = pd.Series(pacbio.base_acc.values, index=pacbio.pb_accs).to_dict()
    with open(f'{name}_orf_refined.fasta', 'w') as ofile:
        for entry in seqs:
            seq = str(entry.seq)
            pb_acc = entry.id
            base_acc = base_map[pb_acc]
            if 'gene' in orfs.columns:
                gene = pb_gene[base_acc]
                ofile.write(f'>pb|{base_acc}|fullname GN={gene}\n{seq}\n')
            else:
                ofile.write(f'>pb|{base_acc}\n{seq}\n')
    
    logging.info('Writing refined database tsv results...')
    if 'gene' in orfs.columns:
        pacbio = pacbio[['pb_accs', 'base_acc', 'coding_score', 'orf_calling_confidence', 'upstream_atgs', 'orf_score', 'gene', 'FL', 'CPM']]
    else:
        pacbio = pacbio[['pb_accs', 'base_acc', 'coding_score', 'orf_calling_confidence', 'upstream_atgs', 'orf_score', 'FL', 'CPM']]
    pacbio.to_csv(f'{name}_orf_refined.tsv', sep='\t', index=False)
    logging.info('Refine Database Complete')
    logging.info('************************')

if __name__ == '__main__':
    main()
