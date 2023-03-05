import csv
import pandas as pds
import numpy as np


###Human PrEDiCT scores
def diseases_dic(tissue):
    '''Dictionary of diseases as keys and the causal genes as values.'''
    diseases = {}
    disease_rdr = csv.reader(open("Data\\Diseases\\" + tissue + ".txt", "r"), delimiter="\t")
    next(disease_rdr)
    for row in disease_rdr:
        if row[11] == "" or row [11] == "Yes":
            if row[7] != "":
                disease_key = (row[7],row[8])
            else:
                disease_key = (row[0],row[1])
            if disease_key not in diseases:
                diseases[disease_key] = []
            if row[3] not in diseases[disease_key]:
                diseases[disease_key].append(row[3])
    return diseases

def expression_tables(tissue):
    '''Returns expression tables: normalized counts, percentage of cells expressing (>=0.5 normalized counts), Z-scores.'''
    normalized = pds.read_table("Data\\Expression\\" + tissue + "\\" + tissue + "_FinalLogAvgEx.txt", sep=" ")
    pct = pds.read_table("Data\\Expression\\" + tissue + "\\" + tissue + "_PCT.txt", sep=" ")
    zscores = pds.read_table("Data\\Expression\\" + tissue + "\\" + tissue + "_Z_scores.txt", sep=" ")
    return normalized, pct, zscores

def expressed_genes(tissue):
    '''Returns the genes that express in at least one cell type of the tissue (>=10% of a cell type and above median
    of average expression).'''
    normalized, pct, zscores = expression_tables(tissue)
    filter_indices = pct.index[pct.max(axis = 1) >= 10]
    filtered_normalized = normalized.loc[filter_indices]
    filtered_normalized = filtered_normalized.dropna(how='all')
    tis_med = np.median(filtered_normalized.values.tolist())
    genes = filtered_normalized.index[filtered_normalized.max(axis = 1) >= tis_med]
    return genes

def predict_scores(tissue):
    '''Writes files with the tissues' PReDiCT scores.'''
    diseases = diseases_dic(tissue)
    normalized, pct, zscores = expression_tables(tissue)
    exp_genes = expressed_genes(tissue)
    wrtr = open("Output\\PrEDiCT score\\" + tissue + ".tsv", "w")
    wrtr.write("Disease\Phenotypic series accession (OMIM)\tDisease\Phenotypic series name (OMIM)\tCell type\tPrEDiCT score\tNumber of causal genes\tAffected cell type (1=Yes, 0=No)\n")
    for dis in diseases:
        intersect_genes = list(set(diseases[dis]) & set(exp_genes))
        if len(intersect_genes) > 0:
            predicts = zscores.loc[intersect_genes].median()
            for i in predicts.iteritems():
                wrtr.write(dis[0] + "\t" +
                           dis[1] + "\t" +
                           i[0] + "\t" +
                           str(i[1]) + "\t" +
                           str(len(intersect_genes)) + "\t")
                if i[1] >= 2:
                    wrtr.write("1\n")
                else:
                    wrtr.write("0\n")


###Mouse PrEDiCT scores
def m_diseases_dic(tissue):
    '''Dictionary of diseases as keys and the causal genes as values.'''
    diseases = {}
    disease_rdr = csv.reader(open("Data\\Diseases\\" + tissue + ".txt", "r"), delimiter="\t")
    next(disease_rdr)
    for row in disease_rdr:
        if row[11] == "" or row [11] == "Yes" and row[4] != "":
            if row[7] != "":
                disease_key = (row[7],row[8])
            else:
                disease_key = (row[0],row[1])
            if disease_key not in diseases:
                diseases[disease_key] = []
            if row[4] not in diseases[disease_key]:
                diseases[disease_key].append(row[4])
    return diseases

def m_expression_tables(tissue):
    '''Returns expression tables: normalized counts, percentage of cells expressing (>=0.5 normalized counts), Z-scores.'''
    normalized = pds.read_table("Data\\Mouse\\Expression\\" + tissue + "\\" + tissue + "_FinalLogAvgEx.txt", sep=" ")
    pct = pds.read_table("Data\\Mouse\\Expression\\" + tissue + "\\" + tissue + "_PCT.txt", sep=" ")
    zscores = pds.read_table("Data\\Mouse\\Expression\\" + tissue + "\\" + tissue + "_Z_scores.txt", sep=" ")
    return normalized, pct, zscores

def m_expressed_genes(tissue):
    '''Returns the genes that express in at least one cell type of the tissue (>=10% of a cell type and above median
    of average expression).'''
    normalized, pct, zscores = expression_tables(tissue)
    filter_indices = pct.index[pct.max(axis = 1) >= 10]
    filtered_normalized = normalized.loc[filter_indices]
    filtered_normalized = filtered_normalized.dropna(how='all')
    tis_med = np.median(filtered_normalized.values.tolist())
    genes = filtered_normalized.index[filtered_normalized.max(axis = 1) >= tis_med]
    return genes

def m_predict_scores(tissue):
    '''Writes files with the tissues' PReDiCT scores.'''
    diseases = diseases_dic(tissue)
    normalized, pct, zscores = expression_tables(tissue)
    exp_genes = expressed_genes(tissue)
    wrtr = open("Output\\Mouse\\PrEDiCT score\\" + tissue + ".tsv", "w")
    wrtr.write("Disease\Phenotypic series accession (OMIM)\tDisease\Phenotypic series name (OMIM)\tCell type\tPrEDiCT score\tNumber of causal genes\tAffected cell type (1=Yes, 0=No)\n")
    for dis in diseases:
        intersect_genes = list(set(diseases[dis]) & set(exp_genes))
        if len(intersect_genes) > 0:
            predicts = zscores.loc[intersect_genes].median()
            for i in predicts.iteritems():
                wrtr.write(dis[0] + "\t" +
                           dis[1] + "\t" +
                           i[0] + "\t" +
                           str(i[1]) + "\t" +
                           str(len(intersect_genes)) + "\t")
                if i[1] >= 2:
                    wrtr.write("1\n")
                else:
                    wrtr.write("0\n")

tissues = ["Skeletal muscle", "Bone marrow", "Spleen", "Lung", "Trachea", "Tongue"]
