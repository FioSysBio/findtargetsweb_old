import sys
#import cobra
import pandas as pd
import time
from googletrans import Translator
'''
start = time.time()

df_list_organism = pd.read_excel("c://workspace-python//findtargetsweb3//FindTargetsWEB//list_organism_kegg.xls")
df_list_organism_filter = df_list_organism[df_list_organism['name_organism'].str.contains("Salmonella typhimurium")]
list_acron_kegg = df_list_organism_filter['acron_organism'].tolist()
print(list_acron_kegg)
time.sleep(5)
end = time.time()

print(end - start, "segundos")


df_summary = pd.read_excel("sumario-result-article-v4.xls", sheet_name="List_Targets_CCBH")
#df_summary = df_summary.dropna()
print(df_summary.columns)
print(df_summary.shape)

df_result = pd.DataFrame(columns=['Via', 'Funcao', 'Localizacao']) 
translator = Translator()
count = 0
for row in df_summary.itertuples():
    print(row.EC_NUMBER)
    via = row.PATHWAY
    funcao = row.FUNCTION
    local = row.LOCALIZATION

    df_result.loc[count] = [translator.translate(via, dest='pt').text, 
                            translator.translate(funcao, dest='pt').text, 
                            translator.translate(local, dest='pt').text]
    
    count += 1
    time.sleep(2)

df_result.to_excel("result_translate_CCBH.xls")

print("ACABOU! SEJA FELIZ!")


list_gen = ['PA0025', 'PA0342', 'PA0350', 'PA0546', 'PA1758', 'PA2967',
            'PA2977', 'PA3108', 'PA3163', 'PA3296', 'PA4411', 'PA4412',
            'PA4414', 'PA4417', 'PA4450', 'PA4750', 'PA4938', 'PA5549']



model = cobra.io.read_sbml_model("c:\\temp\\iPAU1129.xml")

print([r.id for r in model.reactions])

react_biomass = model.reaction.get_by_id("PA14_Biomass")
print(react_biomass.reaction)
sys.exit()

fba = model.optimize()

print(fba.f)

for gen in list_gen:
    gen_by_id = model.genes.get_by_id(gen)
    list_lalala = [i.id for i in list(gen_by_id.reactions)]
    print(gen, ", ".join(list_lalala))

'''

param = ["CCBH", "PA14", "PAO108", "PAO117"]

list_CCBH = []
list_PA14 = []
list_PAO108 = []
list_PAO117 = []
set_ALL = []

count = 1
for p in param:

    with open("lalala_"+str(p)+".txt", "r") as infile:
        data = infile.read()
    my_list = data.splitlines()
    my_list = list(set(my_list))
    my_list.sort()
    
    for item in my_list:
        if count == 1:
            list_CCBH.append("{0}".format(item))
        elif count == 2:
            list_PA14.append("{0}".format(item))
        elif count == 3:
            list_PAO108.append("{0}".format(item))
        else:
            list_PAO117.append("{0}".format(item))
        
        set_ALL.append("{0}".format(item))
        
    count += 1

set_ALL = list(set(set_ALL))
set_ALL.sort()

print("TOTAL ENCONTRADO JUNTANDO AS 4 SEM REPETICAO", len(set_ALL))

#print('\n'.join(set_ALL))


print("EC EM TODOS OS 4 CENARIOS")
lista1 = [item for item in set_ALL if item in list_CCBH and item in list_PA14 and item in list_PAO108 and item in list_PAO117]
print(len(lista1), lista1)

print("EC APENAS NA PAO1 2008")
lista2 = [item for item in set_ALL if item not in list_CCBH and item not in list_PA14 and item in list_PAO108 and item not in list_PAO117] 
print(len(lista2), lista2)

print("EC APENAS NA PAO1 2017")
lista3 = [item for item in set_ALL if item not in list_CCBH and item not in list_PA14 and item not in list_PAO108 and item in list_PAO117] 
print(len(lista3), lista3)

print("EC APENAS NA PA14")
lista4 = [item for item in set_ALL if item not in list_CCBH and item in list_PA14 and item not in list_PAO108 and item not in list_PAO117] 
print(len(lista4), lista4)

print("EC APENAS NA CCBH")
lista5 = [item for item in set_ALL if item in list_CCBH and item not in list_PA14 and item not in list_PAO108 and item not in list_PAO117] 
print(len(lista5), lista5)

print("===========================")
print("===========================")
print("===========================")

print("EC NA PAO1 2008 E PAO1 2017")
lista6 = [item for item in set_ALL if item not in list_CCBH and item not in list_PA14 and item in list_PAO108 and item in list_PAO117] 
print(len(lista6), lista6)

print("EC NA PAO1 2008 E PA14")
lista7 = [item for item in set_ALL if item not in list_CCBH and item in list_PA14 and item in list_PAO108 and item not in list_PAO117]
print(len(lista7), lista7)

print("EC NA PAO1 2008 E CCBH")
lista8 = [item for item in set_ALL if item in list_CCBH and item not in list_PA14 and item in list_PAO108 and item not in list_PAO117]
print(len(lista8), lista8)

print("EC NA PAO1 2017 E PA14")
lista9 = [item for item in set_ALL if item not in list_CCBH and item in list_PA14 and item not in list_PAO108 and item in list_PAO117]
print(len(lista9), lista9)

print("EC NA PAO1 2017 E CCBH")
lista10 = [item for item in set_ALL if item in list_CCBH and item not in list_PA14 and item not in list_PAO108 and item in list_PAO117]
print(len(lista10), lista10)

print("EC NA PA14 E CCBH")
lista11 = [item for item in set_ALL if item in list_CCBH and item in list_PA14 and item not in list_PAO108 and item not in list_PAO117]
print(len(lista11), lista11)

print("===========================")
print("===========================")
print("===========================")

print("EC NA PAO1 2008 E PAO1 2017 E PA14")
lista12 = [item for item in set_ALL if item not in list_CCBH and item in list_PA14 and item in list_PAO108 and item in list_PAO117]
print(len(lista12), lista12)

print("EC NA PAO1 2008 E PAO1 2017 E CCBH")
lista13 = [item for item in set_ALL if item in list_CCBH and item not in list_PA14 and item in list_PAO108 and item in list_PAO117]
print(len(lista13), lista13)

print("EC NA PAO1 2017 E PA14 E CCBH")
lista14 = [item for item in set_ALL if item in list_CCBH and item in list_PA14 and item not in list_PAO108 and item in list_PAO117]
print(len(lista14), lista14)

print("EC NA PAO1 2008 E PA14 E CCBH")
lista15 = [item for item in set_ALL if item in list_CCBH and item in list_PA14 and item in list_PAO108 and item not in list_PAO117]
print(len(lista15), lista15)
