from __future__ import print_function

from datetime import datetime
from bioservices.uniprot import UniProt
from bioservices.ncbiblast import NCBIblast
from bs4 import BeautifulSoup

import requests
import os
import time

TIMEOUT_SECONDS = 200000

with open("08-filter_ECNumbers_drugbank.txt", "r") as infile:
    data = infile.read()
my_list_uniprotid = data.splitlines()
my_list_uniprotid = list(set(my_list_uniprotid))
my_list_uniprotid.sort()

# Instancia de chamada ao UniProt
u = UniProt()
s = NCBIblast()

u.TIMEOUT = TIMEOUT_SECONDS
s.TIMEOUT = TIMEOUT_SECONDS

#my_list_uniprotid = ["P42898"]
fileResultBlast = open("11-hitsEncontradosUniprot.txt", "w")

# GET ALL JOBID FROM UNIPROT FOUND
file_jobid = open("list_jobid.txt", "w", encoding="ISO-8859-1")
start_time_first_req = time.time()
print('Tamanho da lista que vamos iterar pra ver se da problema', len(my_list_uniprotid))
for uniprotid in my_list_uniprotid:

    uniprotid = uniprotid.split(';')[3]
    print(uniprotid)
    # buscando a sequencia da proteina por ID, vindo do drugbank
    #sequence = u.retrieve(uniprotid, "fasta")
    sequence = u.get_fasta_sequence(uniprotid)
    print(sequence)
    #sequence = sequence.split("\n", 1)[1].strip("\n")
    
    # executando o job que faz o blast
    jobid = s.run(program="blastp", database="uniprotkb_reference_proteomes", sequence=sequence, stype="protein", email="thiago.merigueti@ioc.fiocruz.br", alignments='1000')
    print(jobid)
    file_jobid.write("{0};{1}\n".format(uniprotid, jobid))
    time.sleep(1)

elapsed_time = time.time() - start_time_first_req
#print("TEMPO TOTAL PRIMEIRO PASSO = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
file_jobid.close()

# AFTER GENERATE ALL JOBID FROM BLAST, WE ITERATE ALL ITENS AND GET ALL XML RETURN
with open("list_jobid.txt", "r") as infile:
    data = infile.read()
my_list_jobid = data.splitlines()
#my_list_jobid = list(set(my_list_jobid))

arquivo = open("list_jobid.txt", 'r')

for item in arquivo:
    item = item.strip()
    uniprotid = item.split(';')[0]
    jobid = item.split(';')[1]

    url_status_blast = "http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/"
    url_result_blast = "http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/"
    last_url_result_blast = "/xml"

    start_time_first_req = time.time()

    condition = True
    count_error = 0
    url_1 = url_status_blast+jobid
    #print(url_1)
    count_refresh = 0
    while condition:
        response = requests.get(url_1)
        #print(response.text)
        if response.text != "RUNNING":
            if response.text == "FINISHED":
                condition = False

            if response.text == "ERROR":
                count_error += 1
                #print("ERRO!", count_error)
                if count_error == 3:
                    condition = False

        time.sleep(3)
        if count_refresh > 60:
            condition = False

        count_refresh += 1

    elapsed_time = time.time() - start_time_first_req
    #print("TEMPO TOTAL PRIMEIRO LINK = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

    start_time_second_req = time.time()

    if response.text != 'FINISHED':
        fileResultBlast.write("{0};;{1};;{2};;{3};;{4};;{5};;{6};;{7};;{8};;{9}\n".format(
            uniprotid, "ERRO", "ERRO", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        ))
        continue

    url_2 = url_result_blast+jobid+last_url_result_blast
    #print(url_2)
    response = requests.get(url_2)
    soup = BeautifulSoup(response.text)

    elapsed_time_2 = time.time() - start_time_second_req
    #print("TEMPO TOTAL SEGUNDO LINK = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time_2)))

    # ARMAZENO EM OUTRA PASTA OS ARQUIVOS DO BLAST NA INTEGRA
    #fileHitsBlastIntegra = open(dir_blasts+"//all_hits_"+uniprotid+".xml", "w")
    #fileHitsBlastIntegra.write(str(soup))
    #fileHitsBlastIntegra.close()

    # busca por todos os hits encontrados
    listAllHits = soup.findAll('hit')

    # itera os hits e, caso encontre pseudomonas, guarda as informacoes dele
    fileListAllHist = open("//hit_organism_found_"+uniprotid+".txt", 'w')

    percentSimilarHuman = 0.0
    for hit in listAllHits:
        organism = hit['description'].split("=")[1].split("GN")[0].strip()
        if "Homo sapiens" in organism:
            percentSimilarHuman = float(hit.find('identity').string.strip())
            break

    for hit in listAllHits:
        organism = hit['description'].split("=")[1].split("GN")[0].strip()
        fileListAllHist.write("{0}\n".format(organism))

        if organism == 'Pseudomonas aeruginosa (strain ATCC 15692 / DSM 22644 / CIP 104116 / JCM 14847 / LMG 12228 / 1C / PRS 101 / PAO1)' :
        #if organismParam in organism: # se o valor da combo for encontrado no hit, armazena
            idPAO1 = hit['ac']
            percentBlast = hit.find('identity').string.strip()
            eValue = hit.find('expectation').string.strip()

            foundID = u.search(idPAO1, columns='entry name, id, genes, pathway, comment(FUNCTION), comment(CATALYTIC ACTIVITY), database(PDB), database(PSEUDOCAP), comment(SUBCELLULAR LOCATION), ec')
            splitForListFoundID = list(str(foundID).split("\t"))

            # 0-uniprotid;;1-hit pseudomonas;;2-percentidentityblast;;3-evalue;;
            # 4-Gene names;;5-Pathway;;6-Function;;7-Catalitic Activity;;
            # 8-Localization;;9-ID PDB
            if float(percentSimilarHuman) < float(percentBlast): # so grava se o percentual de similaridade for menor do que com humanos
                fileResultBlast.write("{0};;{1};;{2};;{3};;{4};;{5};;{6};;{7};;{8};;{9}\n".format(
                    uniprotid, idPAO1, percentBlast, eValue, splitForListFoundID[11], splitForListFoundID[12],
                    splitForListFoundID[13], splitForListFoundID[14], splitForListFoundID[17], splitForListFoundID[15]
                ))
            else: # Caso tenha hit com percent de humano maior, poe tudo zerado
                fileResultBlast.write("{0};;{1};;{2};;{3};;{4};;{5};;{6};;{7};;{8};;{9}\n".format(
                    uniprotid, idPAO1, percentBlast, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                ))

    fileListAllHist.close()

fileResultBlast.close()
arquivo.close()