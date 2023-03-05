import csv

#Read mouse orthologous genes file from MGI
rdr = csv.reader(open("Data\\Mouse\\HOM_MouseHumanSequence.rpt.txt", "r"), delimiter = "\t")
next(rdr)
genes_dic = {}
mouse_genes = []

#Create a dictionary for orthology
for row in rdr:
    if row[1] == "mouse, laboratory":
        mouse_genes.append(row[3])
    else:
        if mouse_genes[-1] not in genes_dic:
            genes_dic[mouse_genes[-1]] = [row[3]]
        else:
            genes_dic[mouse_genes[-1]].append(row[3])

#Write dictionary as file
wrtr = open("Output\\Mouse\\gene_dic.tsv", "w")
for mouse in genes_dic:
    orthologs = str(genes_dic[mouse]).replace("[","").replace("]","").replace(" ","").replace("'","")
    wrtr.write(mouse + "\t" + orthologs + "\n")
