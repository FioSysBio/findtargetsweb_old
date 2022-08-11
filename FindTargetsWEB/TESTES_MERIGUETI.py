import requests
from bs4 import BeautifulSoup
from bioservices.uniprot import UniProt
from bioservices.ncbiblast import NCBIblast
import pandas as pd
import re
from bioservices.kegg import KEGG
import cobra
import sys
import time 


with open("08-filter_ECNumbers_drugbank.txt", "r") as infile:
    data = infile.read()
my_list_uniprotid = data.splitlines()
my_list_uniprotid = list(set(my_list_uniprotid))
my_list_uniprotid.sort()

# Instancia de chamada ao UniProt
u = UniProt()
s = NCBIblast()

u.TIMEOUT = 200000
s.TIMEOUT = 200000

file_jobid = open("list_jobid.txt", "w", encoding="ISO-8859-1")
start_time_first_req = time.time()
for uniprotid in my_list_uniprotid:

    # COMENTAR ESSA LINHA QUANDO FOR FAZER SIMULACAO COM A LISTA COM 1 ITEM!
    uniprotid = uniprotid.split(';')[3]
  
    # buscando a sequencia da proteina por ID, vindo do drugbank
    sequence = u.retrieve(uniprotid, "fasta")
    sequence = sequence.split("\n", 1)[1].strip("\n")

    # executando o job que faz o blast 
    jobid = s.run(program="blastp", database="uniprotkb_reference_proteomes", sequence=sequence, stype="protein", email="thiago.merigueti@ioc.fiocruz.br", alignments='1000')
    print(jobid)
    file_jobid.write("{0};{1}\n".format(uniprotid, jobid))

elapsed_time = time.time() - start_time_first_req
print("TEMPO TOTAL PRIMEIRO PASSO = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
file_jobid.close()
    

with open("list_jobid.txt", "r") as infile:
    data = infile.read()
my_list_jobid = data.splitlines()
my_list_jobid = list(set(my_list_jobid))

fileResultBlast = open("11-hitsEncontradosUniprot.txt", "w")
organismParam = "Pseudomonas aeruginosa"

for uniprotid, jobid in zip(my_list_uniprotid, my_list_jobid):
    url_status_blast = "http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/"
    url_result_blast = "http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/"
    last_url_result_blast = "/xml"
    
    start_time_first_req = time.time()
    
    condition = True
    count_error = 0
    url_1 = url_status_blast+jobid
    print(url_1)
    count_refresh = 0
    while condition:
        response = requests.get(url_1)
        print(response.text)
        if response.text != "RUNNING":
            if response.text == "FINISHED": 
                condition = False
            
            if response.text == "ERROR":
                count_error += 1
                print("ERRO!", count_error)
                if count_error == 3:
                    condition = False
        
        time.sleep(30)
        if count_refresh > 60:
            condition = False
        
        count_refresh += 1
    
    elapsed_time = time.time() - start_time_first_req
    print("TEMPO TOTAL PRIMEIRO LINK = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    
    start_time_second_req = time.time()
    
    if response.text != 'FINISHED':
        fileResultBlast.write("{0};;{1};;{2};;{3};;{4};;{5};;{6};;{7};;{8};;{9}\n".format(
            uniprotid, "ERRO", "ERRO", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        ))
        continue
    
    url_2 = url_result_blast+jobid+last_url_result_blast
    print(url_2)
    response = requests.get(url_2)
    soup = BeautifulSoup(response.text)
    
    elapsed_time_2 = time.time() - start_time_second_req
    print("TEMPO TOTAL SEGUNDO LINK = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time_2)))
    
    # retorno do blast em formato xml
    # returnXML = s.getResult(jobid, 'xml') # <class 'bioservices.xmltools.easyXML'>
    
    # Instancia de soup para leitura do retorno XML do blast
    # soup = BeautifulStoneSoup(str(returnXML))
    # soup = returnXML.soup

    # ARMAZENO EM OUTRA PASTA OS ARQUIVOS DO BLAST NA INTEGRA
    fileHitsBlastIntegra = open("all_hits_"+uniprotid+".xml", "w")
    fileHitsBlastIntegra.write(str(soup))
    fileHitsBlastIntegra.close()

    # busca por todos os hits encontrados
    listAllHits = soup.findAll('hit')

    # itera os hits e, caso encontre pseudomonas, guarda as informacoes dele
    fileListAllHist = open("hit_organism_found_"+uniprotid+".txt", 'w')
    
    percentSimilarHuman = 0.0
    for hit in listAllHits:
        organism = hit['description'].split("=")[1].split("GN")[0].strip()
        if "Homo sapiens" in organism:
            percentSimilarHuman = float(hit.find('identity').string.strip())
            break
    
    for hit in listAllHits:
        organism = hit['description'].split("=")[1].split("GN")[0].strip()
        fileListAllHist.write("{0}\n".format(organism))
        
        #if organism == 'Pseudomonas aeruginosa (strain ATCC 15692 / DSM 22644 / CIP 104116 / JCM 14847 / LMG 12228 / 1C / PRS 101 / PAO1)' :
        if organismParam in organism: # se o valor da combo for encontrado no hit, armazena
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


class Refactoring:
    def passo7(self):
        #ec = "3.1.1.82"
        ec = "2.5.1.7"
        link01 = "https://www.drugbank.ca/unearth/q?utf8=%E2%9C%93&searcher=targets&query="+str(ec)+"&approved=1&illicit=1&investigational=1&withdrawn=1&experimental=1&us=1&canada=1&eu=1&commit=Apply+Filter"
        r = requests.get(link01)
        
        # verificar se houve retorno da pagina acessada
        if r.status_code == 200:
            
            # Instancia de bs4 do resultado da 1a tela do drugbank
            soup = BeautifulSoup(r.text) 
            
            # filtro para a lista de todos os itens encontrados na busca do requests
            list_found_ec = soup.findAll(attrs={"class": "search-result p-2 mb-4 p-sm-3 mb-sm-2"})
            
            # caso tenha item, inicia as buscas para pegar o resto da informacao
            if len(list_found_ec) > 0: 
                
                # itera na lista de itens encontrados no drugbank
                for item in list_found_ec:
                    em = item.find('em')
                    ec_encontrado = em.text
                    
                    print(ec_encontrado, ec, ec_encontrado == ec)
                    if ec_encontrado != ec:
                        continue
                    
                    # busca pelo link para acesso a 2a tela do drugbank
                    var = item.find('a')
                    
                    # Montagem do link e acesso a 2a tela do drugbank para buscar as informacoes necessarias
                    link02 = "https://www.drugbank.ca"+str(var['href'])
                    r2 = requests.get(link02)
                    if r2.status_code == 200:
                        
                        # Nova instancia de bs4 com as informacoes da 2a tela
                        # e pegando nome da proteina, organismo, uniprot id e drugbank id
                        soup2 = BeautifulSoup(r2.text)    
                        organism_data_class = soup2.find(attrs={"class": "content-container"})
                        protein_name = organism_data_class.findAll('dd')[0].text
                        organism_name = organism_data_class.findAll('dd')[2].text
                        uniprot_id = organism_data_class.find('a').text
                        drugbank_id = str(var['href']).split("/")[-1]
                    
                    else:
                        r.raise_for_status()
            
        else:
            r.raise_for_status()

    def passo8Bioservices(self):
        # Instancia de chamada ao UniProt
        u = UniProt()
        s = NCBIblast()
        
        u.TIMEOUT = 200000
        s.TIMEOUT = 200000
        
        sequence = u.retrieve("P0A749", "fasta")
        sequence = sequence.split("\n", 1)[1].strip("\n")
        
        jobid = s.run(program="blastp", database="uniprotkb_reference_proteomes", sequence=sequence, stype="protein", email="thiago.merigueti@ioc.fiocruz.br")
        
        print(jobid)
        print(type(jobid))
        print(dir(jobid))
        
        retornoOut = s.getResult(jobid, "out")[0:1000]
        print(retornoOut)
        list_retornoOut = [x for x in s.getResult(jobid, "out").split("\n") if "HUMAN" in x]
        print(list_retornoOut)
    
    def passo9Refactoring(self):
        directory = "results//testesMerigueti"

        # PEGANDO A ENTRADA DO ITEM POR MEIO DOS ARQUIVOS GERADOS E CONVERTENDO A DATAFRAME PARA FACILITAR
        data08 = pd.read_csv(directory + "//" + "08-filter_ECNumbers_drugbank.txt", sep=";", header=None)
        data08.columns = ["EC NUMBER", "PROTEIN NAME", "ORGANISM NAME", "UNIPROTID", "DRUGBANKID"]
        data08 = data08.sort_values("EC NUMBER")
        my_list_uniprotid_drugbankid = data08["UNIPROTID"].tolist()
        
        data11 = pd.read_csv(directory + "//" + "11-hitsEncontradosUniprot.txt", sep=";;", header=None)
        data11.columns = ["UNIPROTID", "HIT_UNIPROTID", "PERCENT_IDENT_BLAST", "EVALUE", "GENE_NAME",
                          "PATHWAY", "FUNCTION", "CATALYTIC ACTIVITY", "LOCALIZATION", "ID PDB"]
        data11 = data11.sort_values("UNIPROTID")
        my_list_uniprotid_hit = data11["UNIPROTID"].tolist()
        
        # GERACAO DOS ARQUIVOS DE SAIDA DESSE ITEM
        fileDataDrugs = open(directory + "//" + "13-list_inhibitors_per_target.txt", "w")
        fileInhibitorsDrugs = open(directory + "//" + "14-list_inhibitors_approved.txt", "w")
        
        for uniprot_drugbank in my_list_uniprotid_drugbankid:
        
            if uniprot_drugbank in my_list_uniprotid_hit:
                data08Aux = data08[data08["UNIPROTID"] == uniprot_drugbank]
                drugbankid = data08Aux.iloc[0][4]
        
                linkAccessForDrugsID = "https://www.drugbank.ca/biodb/bio_entities/"+str(drugbankid)
                r = requests.get(linkAccessForDrugsID)
        
                if r.status_code == 200:
                    soup = BeautifulSoup(r.text)
                    listDrugsFound = soup.find(attrs={"class": "table table-sm table-bordered datatable dt-responsive"})
                    listDrugsFound = listDrugsFound.find('tbody')
                    listDrugsFound = listDrugsFound.findAll('tr') # lista com todos os farmacos para o drugbankid informado
        
                    for item in listDrugsFound:
                        drugbank_drug_id = item('td')[0].text
                        drug_name = item('td')[1].text
                        drug_group = item('td')[2].text
                        pharma_action = item('td')[3].text
                        actions = item('td')[4].text
                        fileDataDrugs.write("{0};{1};{2};{3};{4};{5};{6}\n".format(
                            uniprot_drugbank, drugbankid, drugbank_drug_id, drug_name,
                            drug_group, pharma_action, actions
                        ))
        
                        if pharma_action == 'yes' and actions == 'inhibitor':
                            fileInhibitorsDrugs.write("{0};{1};{2};{3};{4};{5};{6}\n".format(
                            uniprot_drugbank, drugbankid, drugbank_drug_id, drug_name,
                            drug_group, pharma_action, actions
                        ))
        
                else:
                    fileDataDrugs.write("{0};{1};{2};{3};{4};{5};{6}\n".format(
                        uniprot_drugbank, drugbankid, r.raise_for_status(), r.raise_for_status(),
                        r.raise_for_status(), r.raise_for_status(), r.raise_for_status()
                    ))
        
        print("GENERATED FILE %s" % fileDataDrugs.name)
        print("GENERATED FILE %s" % fileInhibitorsDrugs.name)

    def extractListOrganismKEGG(self):
        k = KEGG(False, True)

        result_list_organism_kegg = k.list("organism")
        result_list_organism_kegg_splt = result_list_organism_kegg.split('\n')
        result_list_organism_kegg_splt = filter(None, result_list_organism_kegg_splt)
        result_list_organism_kegg_splt_filter = list(result_list_organism_kegg_splt)
        df_list_organism = pd.DataFrame(columns=['id_organism', 'acron_organism', "name_organism"])
        index_list_organism = 0
        
        for line in result_list_organism_kegg_splt_filter:
            line_splt = line.split('\t')
            id_organism = line_splt[0]
            acron_organism = line_splt[1]
            name_organism = line_splt[2]
            df_list_organism.loc[index_list_organism] = [id_organism, acron_organism, name_organism]
            index_list_organism += 1
            
        df_list_organism.to_excel("list_organism_kegg.xls")
    
    def planoAAlternativeStepGetECFromCompound(self):
        '''
        with open("react_biomass_zero_sbml_no_genes.txt", "r") as infile:
            data = infile.read()
        my_list_compound = data.splitlines()
        
        k = KEGG(False, True)
        k.TIMEOUT = 100000
        pd.options.display.max_colwidth = 1000
        
        file_vaziosno = open("vaziosno.txt", 'w', encoding="ISO-8859-1")
        file_ecs_encontrados_kegg = open("ecs_encontrados_kegg.txt", 'w', encoding="ISO-8859-1")
        
        for compound in my_list_compound:
            compound_no_stoich = re.sub("\d+\.\d+", "", compound) # retirada de todos os valores estequiometricos
            
            if " - reduced " in compound:
                compound_no_stoich = compound_no_stoich.replace(" - reduced ", "")
            
            param_splt = "" # preparo do split para separar a composicao quimica
            if "<=>" in compound_no_stoich:
                param_splt = "<=>"
            else:
                param_splt = "-->"
            
            compound_splt = compound_no_stoich.split(param_splt) # separacao da composicao quimica pelo reactante(0) e o produto(1) 
            
            # Caso seja uma reacao de insercao ou uma reacao de excrecao, ele nao busca no kegg
            if compound_splt[0] == "" or compound_splt[1] == "":
                #print("um dos lados esta vazio, entao segue o fluxo")
                #print("========================================")
                continue  
            
            # TRATAMENTO DOS DADOS DO REAGENTE DA COMPOSICAO ENCONTRADA
            reactant_sbml = compound_splt[0].strip()
            product_sbml = compound_splt[1].strip() # dados do produto da composicao quimica em questao
        
            #http://rest.kegg.jp/find/reaction/<REACTANT>
            result_reactant = k.find("reaction", str(reactant_sbml))
            if result_reactant == 400:
                #print("erro 400 no resultado do reactant")
                #print("========================================")
                continue  
            
            if len(result_reactant) == 1:
                file_vaziosno.write("{0}\n".format(compound))
                #print("vazio no resultado do reactant")
                #print("========================================")
                continue
            
            # Quebra do retorno do kegg em linhas e filtro para excluir os vazios da lista gerada
            result_reactant_splt = result_reactant.split('\n')
            result_reactant_splt = filter(None, result_reactant_splt)
            result_reactant_splt_filter = list(result_reactant_splt)
            df_reactant_kegg = pd.DataFrame(columns=['id_reaction', 'dsc_reaction_kegg'])
            
            # PREENCHER OS DADOS DO KEGG DE ID DE REACAO E A REACAO EM UM DATAFRAME
            index_reactant = 0
            for local_result in result_reactant_splt_filter:
                local_result = local_result.split('\t')
                if ";" in local_result[1]:
                    reactant_kegg = local_result[1].split(";")[-1]
                else:
                    reactant_kegg = local_result[1]
                    
                df_reactant_kegg.loc[index_reactant] = [local_result[0], reactant_kegg]
                index_reactant += 1
            
            # CONVERSAO DE TODAS AS INFORMACOES ENCONTRADAS PARA MINUSCULO
            #df_reactant_kegg = df_reactant_kegg.astype(str).apply(lambda x: x.str.lower())
            df_reactant_kegg['dsc_reaction_kegg'] = df_reactant_kegg['dsc_reaction_kegg'].str.lower() 
        
            # TRATAMENTO DOS DADOS DE PRODUTO DA COMPOSICAO ENCONTRADA    
            #http://rest.kegg.jp/find/reaction/<PRODUCT>
            result_product = k.find("reaction", str(product_sbml)) # busca no kegg por reacoes com os dados parecidos 
            if result_product == 400: # caso tenha algum problema na busca, segue o fluxo
                #print("erro 400 no resultado do product")
                #print("========================================")
                continue
            
            if len(result_product) == 1: # caso nao encontre nada, segue o fluxo tb
                file_vaziosno.write("{0}\n".format(compound))
                #print("vazio no resultado do product")
                #print("========================================")
                continue
            
            # Quebra de uma parte do retorno do kegg em linhas e filtro para excluir os vazios
            result_product_splt = result_product.split('\n')
            result_product_splt = filter(None, result_product_splt)
            result_product_splt_filter = list(result_product_splt)
            df_product_kegg = pd.DataFrame(columns=['id_reaction', 'dsc_reaction_kegg'])
            
            # Iteracao dentro da linha para separar os dados de id da reacao no kegg e o nome/composicao da reacao
            index_prod = 0
            for local_result in result_product_splt_filter:
                local_result = local_result.split('\t')
                if ";" in local_result[1]:
                    product_kegg = local_result[1].split(";")[-1]
                else:
                    product_kegg = local_result[1]
                df_product_kegg.loc[index_prod] = [local_result[0], product_kegg]
                index_prod += 1
            
            # CONVERTENDO TODOS OS RESULTADOS ENCONTRADOS NO KEGG PARA MINUSCULO
            #df_product_kegg = df_product_kegg.astype(str).apply(lambda x: x.str.lower())
            df_product_kegg['dsc_reaction_kegg'] = df_product_kegg['dsc_reaction_kegg'].str.lower() 
        
            df_merge = pd.merge(df_reactant_kegg, df_product_kegg, on=['id_reaction', 'dsc_reaction_kegg'])
            list_reactions_merge = df_merge['dsc_reaction_kegg'].tolist()
            
            if len(df_merge) > 0:
                print("composto do sbml = ", compound)
                
                for row in df_merge.itertuples():
                    #http://rest.kegg.jp/link/enzyme/rn:R01530
                    print(row.id_reaction, row.dsc_reaction_kegg)
                    
                    result_ec_number = k.link("enzyme", row.id_reaction)
                    result_ec_number_splt = result_ec_number.split('\n')
                    result_ec_number_splt = filter(None, result_ec_number_splt)
                    result_ec_number_filter = list(result_ec_number_splt)
                    
                    for result_ec in result_ec_number_filter:
                        local_result_ec = result_ec.split('\t')
                        local_result_ec[1] = local_result_ec[1].replace("ec:", "").strip()
                        file_ecs_encontrados_kegg.write("{0}\n".format(local_result_ec[1]))
                    
            else:
                file_vaziosno.write("{0}\n".format(compound))
                print("chegou no merge das listas, mas nao teve intersecao, segue o fluxo")
            
            print("========================================")
        
        file_vaziosno.close()
        file_ecs_encontrados_kegg.close()
        '''   
          
    def alternativeStepGetECNumberFromCompound(self):
        '''
        with open("react_biomass_zero_sbml_no_genes.txt", "r") as infile:
            data = infile.read()
        my_list_compound_sbml = data.splitlines()
        '''
        my_list_compound_sbml = ["--> fe2+",  "h+ + l-tyrosine --> h+ + l-tyrosine", "--> phosphate"]
        
        k = KEGG(False, True)
        k.TIMEOUT = 100000
        pd.options.display.max_colwidth = 1000
        pd.options.display.max_rows = 1000
        
        file_ecs = open("ecnumber_reactant_sbml.txt", "w", encoding="ISO-8859-1")
        file_comp_not_found_from_reactant_sbml = open("idcomp_notfound_kegg.txt", "w", encoding="ISO-8859-1")
        file_reaction_not_found_from_reactant_sbml = open("idreaction_not_found_kegg.txt", "w", encoding="ISO-8859-1")
        
        for compound in my_list_compound_sbml:
            #print(compound)
            compound_no_stoich = re.sub("\d+\.\d+", "", compound) # retirada de todos os valores estequiometricos
            
            if " - reduced " in compound:
                compound_no_stoich = compound_no_stoich.replace(" - reduced ", " ")
            
            param_splt = "" # preparo do split para separar a composicao quimica
            if "<=>" in compound_no_stoich:
                param_splt = "<=>"
            else:
                param_splt = "-->"
            
            compound_splt = compound_no_stoich.split(param_splt) # separacao da composicao quimica pelo reactante(0) e o produto(1) 
            
            # Caso seja uma reacao de insercao ou uma reacao de excrecao, ele ignora
            if compound_splt[0] == "" or compound_splt[1] == "":
                continue  
            
            # TRATAMENTO DOS DADOS DO REAGENTE DA COMPOSICAO ENCONTRADA
            reactant_sbml = compound_splt[0].strip()
        
            list_id_cpd_kegg = []
            len_product_sbml = 0
            
            # verifica os reagentes que possuem mais de um composto envolvido
            if " + " in reactant_sbml: #(ex.: acetyl-coa + atp + bicarbonate)
                reactant_sbml_splt = reactant_sbml.split(" + ")
                len_product_sbml = len(reactant_sbml_splt)
                
                # iteracao dentro dos compostos dos reactants (acetyl-coa + atp + bicarbonate)
                for cpd_product_sbml in reactant_sbml_splt:
                    #http://rest.kegg.jp/find/compound/acetyl-coa
                    result_id_cpd = k.find("compound", cpd_product_sbml)
                    result_id_cpd_splt = result_id_cpd.split('\n')
                    result_id_cpd_splt = filter(None, result_id_cpd_splt)
                    result_id_cpd_splt_filter = list(result_id_cpd_splt)
                    
                    # iteracao dentro do resultado encontrado para um dos compostos do reagente 
                    # (cpd:C00024 Acetyl-CoA; Acetyl coenzyme A)
                    if len(result_id_cpd_splt_filter) > 0:
                        for result_cpd in result_id_cpd_splt_filter:
                            local_result_cpd_splt = result_cpd.split('\t')
                            if ";" in local_result_cpd_splt[1]:
                                dsc_cpd_splt = local_result_cpd_splt[1].split(';')
                                for cpd_splt in dsc_cpd_splt:
                                    if cpd_splt.strip().lower() == cpd_product_sbml.strip().lower():
                                        list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                        break
                            else:
                                if cpd_product_sbml.lower() == local_result_cpd_splt[1].lower():
                                    list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                    break
            else:
                len_product_sbml = 1
                result_id_cpd = k.find("compound", reactant_sbml)
                result_id_cpd_splt = result_id_cpd.split('\n')
                result_id_cpd_splt = filter(None, result_id_cpd_splt)
                result_id_cpd_splt_filter = list(result_id_cpd_splt)
                
                if len(result_id_cpd_splt_filter) > 0:
                    for result_cpd in result_id_cpd_splt_filter:
                        local_result_cpd_splt = result_cpd.split('\t')
                        if ";" in local_result_cpd_splt[1]:
                            dsc_cpd_splt = local_result_cpd_splt[1].split(';')
                            for cpd_splt in dsc_cpd_splt:
                                if cpd_splt.strip().lower() == reactant_sbml.strip().lower():
                                    list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                    break
                        else:
                            if reactant_sbml.strip().lower() == local_result_cpd_splt[1].strip().lower():
                                list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                break
        
            # caso nao tenha encontrado todos os compostos, ele nao segue o fluxo
            if len(list_id_cpd_kegg) != len_product_sbml:
                file_comp_not_found_from_reactant_sbml.write("{0}\n".format(compound))
                continue
            
            # PREPARO DA APLICACAO PARA USAR OS IDS DE COMPOSTOS E BUSCAR PELAS REACOES DE CADA UM
            reactant_sbml_in_cpd = "+".join(list_id_cpd_kegg)
            result_link_reactions_cpd = k.link("reaction", str(reactant_sbml_in_cpd))
            result_link_reactions_cpd_splt = result_link_reactions_cpd.split('\n')
            result_link_reactions_cpd_splt = filter(None, result_link_reactions_cpd_splt)
            result_link_reactions_cpd_splt_filter = list(result_link_reactions_cpd_splt)
            df_link_reaction_cpd = pd.DataFrame(columns=['id_cpd', 'id_reaction'])
            
            index_link_react_cpd = 0
            if len(result_link_reactions_cpd_splt_filter) > 0:
                for item_result_link_reactions_cpd in result_link_reactions_cpd_splt_filter:
                    local_item_result_link = item_result_link_reactions_cpd.split('\t')
                    df_link_reaction_cpd.loc[index_link_react_cpd] = [local_item_result_link[0], local_item_result_link[1]]
                    index_link_react_cpd += 1
            else:
                file_reaction_not_found_from_reactant_sbml.write("{0}\n".format(compound))
                continue
            
            df_link_reaction_cpd = df_link_reaction_cpd.groupby("id_reaction").filter(lambda x: len(x) == len(reactant_sbml_splt))
            set_id_reaction_kegg = {x.id_reaction for x in df_link_reaction_cpd.itertuples()}
            
            if len(set_id_reaction_kegg) > 0:
                for item_id_react in set_id_reaction_kegg:
                    result_ec_number = k.link("enzyme", item_id_react)
                    result_ec_number_splt = result_ec_number.split('\n')
                    result_ec_number_splt = filter(None, result_ec_number_splt)
                    result_ec_number_filter = list(result_ec_number_splt)
                    
                    for result_ec in result_ec_number_filter:
                        local_result_ec = result_ec.split('\t')
                        local_result_ec[1] = local_result_ec[1].replace("ec:", "").strip()
                        file_ecs.write("{0}\n".format(local_result_ec[1]))
            
        file_comp_not_found_from_reactant_sbml.close()
        file_reaction_not_found_from_reactant_sbml.close()
        
        with open("idcomp_notfound_kegg.txt", "r") as infile:
            data = infile.read()
        my_list_compound_sbml_not_found_from_reactant = data.splitlines()
        
        file_comp_not_found_from_product_sbml = open("idcomp_notfound_2_kegg.txt", "w", encoding="ISO-8859-1")
        file_reaction_not_found_from_product_sbml = open("idreaction_notfound_2_kegg.txt", "w", encoding="ISO-8859-1")
        
        for compound in my_list_compound_sbml_not_found_from_reactant:
            print(compound)
            compound_no_stoich = re.sub("\d+\.\d+", "", compound) # retirada de todos os valores estequiometricos
            
            if "- reduced" in compound:
                compound_no_stoich = compound_no_stoich.replace("- reduced", " ")
            
            param_splt = "" # preparo do split para separar a composicao quimica
            if "<=>" in compound_no_stoich:
                param_splt = "<=>"
            else:
                param_splt = "-->"
            
            compound_splt = compound_no_stoich.split(param_splt) # separacao da composicao quimica pelo reactante(0) e o produto(1) 
            
            # Caso seja uma reacao de insercao ou uma reacao de excrecao, ele ignora
            if compound_splt[0] == "" or compound_splt[1] == "":
                continue  
            
            product_sbml = compound_splt[1].strip()
        
            list_id_cpd_kegg = []
            len_product_sbml = 0
            
            # verifica os reagentes que possuem mais de um composto envolvido
            if " + " in product_sbml: #(ex.: acetyl-coa + atp + bicarbonate)
                product_sbml_splt = product_sbml.split(" + ")
                len_product_sbml = len(product_sbml_splt)
                
                # iteracao dentro dos compostos dos reactants (acetyl-coa + atp + bicarbonate)
                for cpd_product_sbml in product_sbml_splt:
                    #http://rest.kegg.jp/find/compound/acetyl-coa
                    result_id_cpd = k.find("compound", cpd_product_sbml)
                    result_id_cpd_splt = result_id_cpd.split('\n')
                    result_id_cpd_splt = filter(None, result_id_cpd_splt)
                    result_id_cpd_splt_filter = list(result_id_cpd_splt)
                    
                    # iteracao dentro do resultado encontrado para um dos compostos do reagente 
                    # (cpd:C00024 Acetyl-CoA; Acetyl coenzyme A)
                    if len(result_id_cpd_splt_filter) > 0:
                        for result_cpd in result_id_cpd_splt_filter:
                            local_result_cpd_splt = result_cpd.split('\t')
                            if ";" in local_result_cpd_splt[1]:
                                dsc_cpd_splt = local_result_cpd_splt[1].split(';')
                                for cpd_splt in dsc_cpd_splt:
                                    if cpd_splt.strip().lower() == cpd_product_sbml.strip().lower():
                                        list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                        break
                            else:
                                if cpd_product_sbml.lower() == local_result_cpd_splt[1].lower():
                                    list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                    break
            else:
                len_product_sbml = 1
                result_id_cpd = k.find("compound", product_sbml)
                result_id_cpd_splt = result_id_cpd.split('\n')
                result_id_cpd_splt = filter(None, result_id_cpd_splt)
                result_id_cpd_splt_filter = list(result_id_cpd_splt)
                
                if len(result_id_cpd_splt_filter) > 0:
                    for result_cpd in result_id_cpd_splt_filter:
                        local_result_cpd_splt = result_cpd.split('\t')
                        if ";" in local_result_cpd_splt[1]:
                            dsc_cpd_splt = local_result_cpd_splt[1].split(';')
                            for cpd_splt in dsc_cpd_splt:
                                if cpd_splt.strip().lower() == product_sbml.strip().lower():
                                    list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                    break
                        else:
                            if product_sbml.strip().lower() == local_result_cpd_splt[1].strip().lower():
                                list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                break
        
            # caso nao tenha encontrado todos os compostos, ele nao segue o fluxo
            if len(list_id_cpd_kegg) != len_product_sbml:
                file_comp_not_found_from_product_sbml.write("{0}\n".format(compound))
                continue
            
            # PREPARO DA APLICACAO PARA USAR OS IDS DE COMPOSTOS E BUSCAR PELAS REACOES DE CADA UM
            product_sbml_in_cpd = "+".join(list_id_cpd_kegg)
            result_link_reactions_cpd = k.link("reaction", str(product_sbml_in_cpd))
            result_link_reactions_cpd_splt = result_link_reactions_cpd.split('\n')
            result_link_reactions_cpd_splt = filter(None, result_link_reactions_cpd_splt)
            result_link_reactions_cpd_splt_filter = list(result_link_reactions_cpd_splt)
            df_link_reaction_cpd = pd.DataFrame(columns=['id_cpd', 'id_reaction'])
            
            index_link_react_cpd = 0
            if len(result_link_reactions_cpd_splt_filter) > 0:
                for item_result_link_reactions_cpd in result_link_reactions_cpd_splt_filter:
                    local_item_result_link = item_result_link_reactions_cpd.split('\t')
                    df_link_reaction_cpd.loc[index_link_react_cpd] = [local_item_result_link[0], local_item_result_link[1]]
                    index_link_react_cpd += 1
            else:
                file_reaction_not_found_from_product_sbml.write("{0}\n".format(compound))
                continue
            
            df_link_reaction_cpd = df_link_reaction_cpd.groupby("id_reaction").filter(lambda x: len(x) == len(product_sbml_splt))
            set_id_reaction_kegg = {x.id_reaction for x in df_link_reaction_cpd.itertuples()}
            
            if len(set_id_reaction_kegg) > 0:
                for item_id_prod in set_id_reaction_kegg:
                    result_ec_number = k.link("enzyme", item_id_prod)
                    result_ec_number_splt = result_ec_number.split('\n')
                    result_ec_number_splt = filter(None, result_ec_number_splt)
                    result_ec_number_filter = list(result_ec_number_splt)
                    
                    for result_ec in result_ec_number_filter:
                        local_result_ec = result_ec.split('\t')
                        local_result_ec[1] = local_result_ec[1].replace("ec:", "").strip()
                        file_ecs.write("{0}\n".format(local_result_ec[1]))
            
        file_ecs.close()
        file_comp_not_found_from_product_sbml.close()
        file_reaction_not_found_from_product_sbml.close()
        
    def newPasso9_20171212(self):
        TIMEOUT_SECONDS = 200000
        
        # Instancia de chamada ao UniProt
        u = UniProt()
        s = NCBIblast()
        
        u.TIMEOUT = TIMEOUT_SECONDS
        s.TIMEOUT = TIMEOUT_SECONDS
        
        # buscando a sequencia da proteina por ID, vindo do drugbank
        sequence = u.retrieve("P0AEK2", "fasta")
        sequence = sequence.split("\n", 1)[1].strip("\n")
        
        # executando o job que faz o blast 
        jobid = s.run(program="blastp", database="uniprotkb_reference_proteomes", sequence=sequence, stype="protein", email="thiago.merigueti@ioc.fiocruz.br", alignments='1000')
        #jobid = "ncbiblast-R20171212-145943-0704-33819162-p2m"
        print(jobid)
        
        url_status_blast = "http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/"
        url_result_blast = "http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/"
        last_url_result_blast = "/xml"
        
        start_time_first_req = time.time()
        
        condition = True
        count_error = 0
        url_1 = url_status_blast+jobid
        print(url_1)
        
        while condition:
            response = requests.get(url_1)
            print(response.text)
            if response.text != "RUNNING":
                if response.text == "FINISHED": 
                    condition = False
                
                if response.text == "ERROR":
                    count_error += 1
                    print("ERRO!", count_error)
                    if count_error == 3:
                        condition = False
            
            time.sleep(2)
        
        elapsed_time = time.time() - start_time_first_req
        print("TEMPO TOTAL PRIMEIRO LINK = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        
        start_time_second_req = time.time()
        
        if response.text == 'FINISHED':
            url_2 = url_result_blast+jobid+last_url_result_blast
            print(url_2)
            response = requests.get(url_2)
            soup = BeautifulSoup(response.text)
            listAllHits = soup.findAll('hit')
            print(listAllHits)
            
        else:
            print("OCORREU ERRO!")
        
        elapsed_time_2 = time.time() - start_time_second_req
        print("TEMPO TOTAL SEGUNDO LINK = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time_2)))

'''
df_list_organism = pd.read_excel("list_organism.kegg.xls")
df_list_organism_filter = df_list_organism[df_list_organism['name_organism'].str.contains("parametro da combo da tela")]
list_acron_kegg = df_list_organism_filter['acron_organism'].tolist()
'''