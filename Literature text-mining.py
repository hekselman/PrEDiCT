from Bio import Entrez
import csv
import os
import pandas as pds
from scipy.stats import chi2_contingency
import statsmodels.api as sm
import scipy

def clean(term):
    '''Cleans phrases from ats and commas.'''
    term = term.replace(",", "")
    term = term.replace("@", "")
    return term

def search(query):
    '''Search a query in PubMed using Biopython package.'''
    Entrez.email = #Enter email
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmode='xml',
                            retmax='100000',
                            term=query)
    results = Entrez.read(handle)
    return results

def count_papers():
    '''Per tissue, creates tables with raw counts of papers that include both disease and cell type names.'''
    tissues = ["Skeletal muscle", "Bone marrow", "Spleen", "Lung", "Trachea", "Tongue"]
    for tissue in tissues:
        print("======\n" + tissue + " is now being processed")
        wrtr = open("Output\\Literature\\Raw_count\\"+tissue+".tsv", "w")
        cell_types_rdr = csv.reader(open("Data\\Human cell types\\" + tissue + ".txt", "r"), delimiter="\t")
        cell_types = []
        for row in cell_types_rdr:
            cell_types.append(row[0])
        cts_papers = {}
        print("Getting papers with cell types")
        for ct in cell_types:
            cts_papers[ct] = set(search(ct)["IdList"])
        diseases = []
        disease_rdr = csv.reader(open("Data\\Diseases\\" + tissue + ".txt", "r"), delimiter="\t")
        next(disease_rdr)
        for row in disease_rdr:
            if row[11] == "" or row [11] == "Yes":
                if row[7] != "":
                    diseases.append((row[7],row[8]))
                else:
                    diseases.append((row[0],row[1]))
        diseases = list(set(diseases))
        dis_papers = {}
        print("Getting diseases papers")
        count = 0
        percent_done = 0
        for dis in diseases:
            count += 1
            check = int(float(count)*100/len(diseases))
            if check != percent_done:
                print(str(percent_done) + "%")
                percent_done = check
            dis_papers[dis] = set(search(clean(dis[1]))["IdList"])
        print("Writing into a file")
        wrtr.write("\t")
        for ct in cell_types:
            wrtr.write("\t" + ct)
        wrtr.write("\n")
        for dis in diseases:
            wrtr.write(dis[0] + "\t" + dis[1])
            for ct in cell_types:
                overlap = len(cts_papers[ct] & dis_papers[dis])
                wrtr.write("\t" + str(overlap))
            wrtr.write("\n")

def sig_pairs():
    '''Per tissue, writes files with disease-cell-type pairs that significantly (adjusted P < 0.001) co-appear in the
    literature.'''
    count_pairs = 0
    for a, b, files in os.walk("Output\\Literature\\Raw_count\\"):
        for file in files:
            count_table = pds.read_table("Output\\Literature\\Raw_count\\" + file)
            count_table["Unnamed: 0"] = count_table[["Unnamed: 0", "Unnamed: 1"]].agg('_'.join, axis=1)
            count_table = count_table.set_index("Unnamed: 0")
            count_table = count_table.drop("Unnamed: 1", axis=1)
            count_table = count_table.loc[(count_table.sum(axis=1) > 0), (count_table.sum(axis=0) > 2)]
            table = sm.stats.Table(count_table)
            threshold = abs(scipy.stats.norm.ppf((0.001 / (count_table.size * 2)))) #Adjustment of P values
            table = table.standardized_resids > threshold
            all_pairs = table.where(table == True).stack().index.tolist()
            wrtr = open("Output\\Literature\\Pairs\\" + file, "w")
            wrtr.write("Disease\Phenotypic series accession (OMIM)\tDisease\Phenotypic series name (OMIM)\tCell type name\n")
            for pair in all_pairs:
                count_pairs += 1
                wrtr.write(pair[0].split("_")[0] + "\t")
                wrtr.write(pair[0].split("_")[1] + "\t")
                wrtr.write(pair[1] + "\n")
