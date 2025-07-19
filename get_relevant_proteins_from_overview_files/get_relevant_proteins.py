import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle as pkl
from ete3 import NCBITaxa
import numpy as np
from statsmodels.graphics.mosaicplot import mosaic
import re
import yaml
import argparse

## Utils
def replace_semicolons(var):
    if ";" in var:
        var = var.replace(';', '-')
    return var

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process protein overview file with keyword filtering.")
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="Path to YAML configuration file (default: config.yaml)",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="./test",
        help="Directory to save output files (default: ./test)",
    )
    return parser.parse_args()


# Load yaml file
def load_yaml(config_path):
    print(f'Loading config.yaml file...')
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
        overview_file_path = config['input_files']['overview_file_path']
        taxonomy_file_path = config['input_files']['taxonomy_file_path']
        cdhit_default_thresh = config['cdhit_default_thresh']
        cdhit_labelled_thresh = config['cdhit_labelled_thresh']
        cdhit_unlabelled_thresh = config['cdhit_unlabelled_thresh']
        keyword_categories = config['keyword_categories']
        interproscan_annotations = config['interproscan_annotations']
        print(f'File config.yaml loaded.')
        if interproscan_annotations:
            interpro_file_path = config['input_files']['interpro_file_path']
        else:
            interpro_file_path = ''
    
    return overview_file_path, taxonomy_file_path, interpro_file_path, interproscan_annotations, cdhit_default_thresh, cdhit_labelled_thresh, cdhit_unlabelled_thresh, keyword_categories

# In case Interproscan annotations are used
def merge_descriptions(row):
    if pd.notna(row['Interpro_accession_description_representative']):
        return f"{row['description']} | Interproscan: {row['Interpro_accession_description_representative']}"
    else:
        return row['description']


# Overview Library
def overview_table(overview_file_path='Overview_library.tsv', taxonomy_file_path='taxonomic_ranking_library.csv'):
    print('Loading overview table with all the proteins')
    overview = pd.read_csv(overview_file_path, sep = '\t')
    print(f'Table containing {len(overview)} entries')
    print('Loading taxonomy table...')
    taxonomy = pd.read_csv(taxonomy_file_path, index_col=0)
    taxonomy['domain / acellular root'] = taxonomy['domain'].combine_first(taxonomy['acellular root'])
    taxonomy = taxonomy[taxonomy['tax_id'].isna() ==False]
    taxonomy['tax_id'] = taxonomy['tax_id'].astype(int).astype(str)
    taxonomy= taxonomy.drop(columns=['domain', 'acellular root'])

    # Merge overview with taxonomy and get the domain level
    overview = create_formatted_overview_table(overview, taxonomy, cdhit_default_thresh=0.9)
    return overview, taxonomy


def create_formatted_overview_table(overview, taxonomy, cdhit_default_thresh=0.9):
    print('Table formatting ongoing...')
    overview =  pd.merge(overview[['aa_seq', 'id', 'species_exact', 'Genbank/Refseq accession', 'description']], taxonomy[['accession', 'domain / acellular root', 'tax_id']], left_on= 'Genbank/Refseq accession',
         right_on='accession').drop(columns='accession')
    # Create labels, cdhit_threshold and source columns
    overview['labels'] = np.nan
    overview['cdhit_threshold'] = cdhit_default_thresh
    overview['source'] = overview['id'].apply(lambda x: 'NCBI GenBank' if 'GCA' in x else 'NCBI RefSeq')
    # Rename columns according to pipeline guidelines
    table_formatted = overview.rename(columns={'id' : 'accession', 'species_exact' : 'source_organism', 'tax_id' : 'source_taxid', 'aa_seq' : 'ProteinSeq', 'description' : 'original_annotation', 'domain / acellular root' : 'category'})
    table_formatted.source_taxid = table_formatted.source_taxid.astype(int).astype(str)
    # Keep only important columns and reset index
    table_formatted = table_formatted[['ProteinSeq', 'accession', 'source_organism', 'Genbank/Refseq accession', 'original_annotation', 'source', 'category', 'source_taxid', 'labels', 'cdhit_threshold']].reset_index()
    table_formatted['index'] = table_formatted['index'].apply(lambda x: 'cancer_'+str(x))
    table_formatted.rename(columns={'index' : 'proteinID'}, inplace=True)
    overview = overview.reset_index()
    overview['index'] = overview['index'].apply(lambda x: 'cancer_' + str(x))
    overview.rename(columns={'index' : 'proteinID'}, inplace = True)
    # Replace the semicolons with '-' since the default csv separator in the library pipeline is the semicolon
    table_formatted['original_annotation'] = table_formatted['original_annotation'].apply(replace_semicolons)
    # Save in the current directory as table_library_complete.csv
    table_formatted.to_csv('table_library_complete.csv', sep=';', index=False)
    print('Table properly formatted.')
    return overview

#  Make directory interproscan_all_noMetacyc.tsv a variable to indicate in the script
def interpro_annotations(interpro_file_path):
    print('Interproscan annotations set to True')
    print('Merging table with additional annotations with Interproscan...')
    annotations_interpro = pd.read_csv(interpro_file_path, sep='\t', header=None)
    # Rename columns
    annotations_interpro.columns = ['Query_accession', 'MD5_digest', 'Length', 'Analysis_DB', 'Signature_accession', 'Signature_description', 'Query_start', 'Query_end', 'Score', 'Status_Match (T or F)', 'Date_run', 'InterPro_annotations_accession',
    'Interpro_accession_description', 'GO annotations']
    # Make the table slim by filtering out unnecessary ones
    annotations_slim = annotations_interpro[['Query_accession', 'Length', 'Analysis_DB', 'Signature_accession', 'Signature_description', 'InterPro_annotations_accession', 'Interpro_accession_description']].replace('-', np.nan)
    # Rename the columns and get the not nan accessions
    annotations_grouped = annotations_slim.groupby('Query_accession', as_index=False).agg({
        'Analysis_DB': lambda x: sorted(set(x.dropna())),
        'Signature_accession': lambda x: sorted(set(x.dropna())),
        'Signature_description': lambda x: sorted(set(x.dropna())),
        'InterPro_annotations_accession': lambda x: sorted(set(x.dropna())),
        'Interpro_accession_description': lambda x: sorted(set(x.dropna())),
    })

    # Representative (mode) values for Interpro accession & description
    representatives = annotations_slim.groupby('Query_accession').agg({
        'InterPro_annotations_accession': lambda x: x.dropna().mode().iloc[0] if not x.dropna().empty else None,
        'Interpro_accession_description': lambda x: x.dropna().mode().iloc[0] if not x.dropna().empty else None
    }).reset_index()

    # Rename the columns to indicate they're representative
    representatives = representatives.rename(columns={
        'InterPro_annotations_accession': 'Interpro_annotations_accession_representative',
        'Interpro_accession_description': 'Interpro_accession_description_representative'
    })

    # Merge into the grouped table
    annotations_grouped = annotations_grouped.merge(representatives, on='Query_accession', how='left')
    print('Tables merged. Interproscan annotations included.')
    return annotations_grouped

def get_selected_proteins_from_keywords(overview, keyword_categories, interproscan_annotations, cdhit_thresholds=[0.99, 0.65]):
    print('Getting relevant proteins from keywords list...')
    # This dict will hold results per category
    category_hits = {}

    # Choose which column to search in
    column_to_search = 'merged_descriptions' if interproscan_annotations else 'description'

    for category, keywords in keyword_categories.items():
        matched = overview[overview[column_to_search].apply(
            lambda text: any(re.search(kw, text, re.IGNORECASE) for kw in keywords if 'prepilin' not in text)
        )]
        matched['labels'] = category.replace('_keywords', '')  # or keep the name as-is
        category_hits[category] = matched

    df_combined = pd.concat(category_hits.values(), ignore_index=True).drop_duplicates(subset='id')

    df_grouped = (
        df_combined.groupby('id')  # ← change this!
        .agg({
            'labels': lambda x: ','.join(sorted(set(x))),
            **{col: 'first' for col in df_combined.columns if col != 'labels' and col != 'id'}
        })
        .reset_index()
    )
    # df_grouped
    print('Updating table...')
    # Update the overviw table and regenerate the standard formatted table to input in the pipeline
    overview_updated = pd.merge(overview, df_grouped[['id',  'Genbank/Refseq accession', 'labels']], on=['id', 'Genbank/Refseq accession'], how='left').drop(columns=['labels_x'])
    overview_updated.rename(columns={'labels_y': 'labels'}, inplace=True)

    overview_updated.rename(columns={'id' : 'accession', 'species_exact' : 'source_organism', 'tax_id' : 'source_taxid', 'aa_seq' : 'ProteinSeq', 'description' : 'original_annotation', 'domain / acellular root' : 'category'}, inplace=True)

    overview_updated['cdhit_threshold'] = np.where(overview_updated['labels'].notnull(), cdhit_thresholds[0], cdhit_thresholds[1])
    if interproscan_annotations:
        table_formatted = overview_updated[['proteinID', 'ProteinSeq', 'accession', 'source_organism',
            'Genbank/Refseq accession', 'original_annotation', 'source', 'category',
            'source_taxid', 'labels', 'cdhit_threshold', 'Interpro_annotations_accession_representative',
            'Interpro_accession_description_representative', 'merged_descriptions']]
    else:
         table_formatted = overview_updated[['proteinID', 'ProteinSeq', 'accession', 'source_organism',
            'Genbank/Refseq accession', 'original_annotation', 'source', 'category',
            'source_taxid', 'labels', 'cdhit_threshold']]

    print('Table updated and properly formatted.')
    return table_formatted


# Load yaml file and get all variables
args = parse_arguments()
config_path = args.config
overview_file_path, taxonomy_file_path, interpro_file_path, interproscan_annotations, cdhit_default_thresh, cdhit_labelled_thresh, cdhit_unlabelled_thresh, keyword_categories = load_yaml(config_path)

# Load tables
overview, taxonomy = overview_table(overview_file_path=overview_file_path, taxonomy_file_path=taxonomy_file_path)

# Save formatted table library
overview = create_formatted_overview_table(overview=overview, taxonomy=taxonomy)


# By default no using Interproscan annotations
if interproscan_annotations:
    # if there are interproscan annotations, add them and merge them with existing table
    annotations_interpro = interpro_annotations(interpro_file_path)
    # Merge with overview table
    overview = pd.merge(overview, annotations_interpro[['Query_accession', 'Analysis_DB', 'Signature_description', 'InterPro_annotations_accession', 'Interpro_annotations_accession_representative', 'Interpro_accession_description_representative']], left_on='proteinID', right_on='Query_accession', how='left').drop(columns=['Query_accession'])
    overview['merged_descriptions'] = overview.apply(merge_descriptions, axis=1)
#     column_annotations = 'merged_descriptions'
# else:
#     column_annotations = 'description'


# Use keywords to generate the labels to add to the final table
table_formatted = get_selected_proteins_from_keywords(overview=overview, 
                                                      keyword_categories=keyword_categories, 
                                                      interproscan_annotations=interproscan_annotations,
                                                      cdhit_thresholds=[cdhit_labelled_thresh, cdhit_unlabelled_thresh])
print('Saving to table_library_complete.csv ...')
table_formatted.to_csv('table_library_complete.csv', sep=';', index=False)
print('Table saved.')
