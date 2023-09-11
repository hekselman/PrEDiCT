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
    wrtr.write("Disease\Phenotypic series accession (OMIM)\tDisease\Phenotypic series name (OMIM)\tCell type\tPrEDiCT score\tNumber of causal genes\n")
    for dis in diseases:
        intersect_genes = list(set(diseases[dis]) & set(exp_genes))
        if len(intersect_genes) > 0:
            predicts = zscores.loc[intersect_genes].median()
            for i in predicts.iteritems():
                wrtr.write(dis[0] + "\t" +
                           dis[1] + "\t" +
                           i[0] + "\t" +
                           str(i[1]) + "\t" +
                           str(len(intersect_genes)) + "\n")

def all_genes_ct_table():
    '''Write a table with all cell types (columns) and genes (rows). Genes not present in a cell type are counted as 0.'''
    tissues = ["Skeletal muscle", "Bone marrow", "Spleen", "Lung", "Trachea", "Tongue"]
    full_table = pds.DataFrame()
    for tis in tissues:
        normalized, pct, zscores = expression_tables(tis)
        expressed = expressed_genes(tis)
        zscores = zscores.add_suffix("_" + tis)
        full_table = pds.concat([full_table, zscores.loc[expressed]], axis=1)
    full_table.fillna(0.0, inplace=True)
    full_table.to_csv("Output\\all_genes_expression.tsv", sep="\t")

def random_PrEDiCT_allT_total(iterations):
    '''Writes a file the same as PrEDiCT scores with X columns (=iterations) of randomly
    calculated scores with the number of genes equal to the number of causal genes of the disease.'''
    tissues = ["Skeletal muscle", "Bone marrow", "Spleen", "Lung", "Trachea", "Tongue"]
    full_exp = pds.read_csv("Output\\all_genes_expression.tsv", sep="\t")
    for tissue in tissues:
        print(tissue)
        print(datetime.now())
        exp_genes = expressed_genes(tissue)
        dis_gene_dics = diseases_dic(tissue)
        rdr = csv.reader(open("Output\\PrEDiCT score\\" + tissue + ".tsv", "r"), delimiter="\t")
        wrtr = open("Output\\PrEDiCT_Literature_Random_allT_total\\" + tissue + ".tsv", "w")
        count = 0
        for row in rdr:
            for box in row:
                wrtr.write(box + "\t")
            for i in range(0, iterations):
                if count == 0:
                    wrtr.write("Random PrEDiCT " + str(i + 1) + "\t")
                else:
                    intersect_genes = list(set(dis_gene_dics[(row[0], row[1])]) & set(exp_genes))
                    rand_z = random.sample(list(full_exp[row[2] + "_" + tissue]), k=len(intersect_genes))
                    score = float(np.median(rand_z))
                    wrtr.write(str(score) + "\t")
            wrtr.write("\n")
            count += 1

def predict_permutation_allT_total():
    '''Per disease, counts the number of random predict scores that >= true predict scores (divided by n of
    randomizations).'''
    tissues = ["Skeletal muscle", "Bone marrow", "Spleen", "Lung", "Trachea", "Tongue"]
    for tissue in tissues:
        print(tissue)
        rdr = csv.reader(open("Output\\PrEDiCT_Literature_Random_allT_total\\" + tissue + ".tsv", "r"), delimiter="\t")
        first_row = True
        p_values = []
        file = []
        wrtr = open("Output\\PrEDiCT_Literature_Permutation_allT_total\\" + tissue + ".tsv", "w")
        for row in rdr:
            if first_row:
                first_row = False
                for box in row[:-1]:
                    wrtr.write(box + "\t")
                wrtr.write("P value\n")
            else:
                randoms = row[7:-1]
                randoms = [float(i) for i in randoms]
                randoms = sum(map(lambda x: x >= float(row[3]), randoms))
                p = randoms / 1000.0
                p_values.append(p)
                file.append(row[:-1] + [str(p)])
        for n in range(len(file)):
            for box in file[n][:-1]:
                wrtr.write(box + "\t")
            wrtr.write(file[n][-1] + "\n")

def predict_permutation_dis_allT_total():
    '''Per disease, counts the number of random predict scores that >= true predict scores (divided by n of
    randomizations). FDR corrected among cell types of each disease'''
    tissues = ["Skeletal muscle", "Bone marrow", "Spleen", "Lung", "Trachea", "Tongue"]
    for tissue in tissues:
        print(tissue)
        first_row = True
        rows_dic = {}
        rdr = csv.reader(open("Output\\PrEDiCT_Literature_Permutation_allT_total\\" + tissue + ".tsv", "r"), delimiter="\t")
        wrtr = open("Output\\PrEDiCT_Literature_Permutation_allT_total_perDisease\\" + tissue + ".tsv", "w")
        for row in rdr:
            if first_row:
                first_row = False
                for box in row:
                    wrtr.write(box + "\t")
                wrtr.write("FDR_Disease\n")
            else:
                if row[0] not in rows_dic:
                    rows_dic[row[0]] = [row]
                else:
                    rows_dic[row[0]] += [row]
        for dis in rows_dic:
            p_values = []
            for row in rows_dic[dis]:
                p_values.append(float(row[-1]))
            fdrs = fdr.fdrcorrection(p_values)[1]
            for n in range(0,len(rows_dic[dis])):
                for box in rows_dic[dis][n]:
                    wrtr.write(box + "\t")
                wrtr.write(str(fdrs[n]) + "\n")

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
