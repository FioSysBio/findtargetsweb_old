from pipelineFindTargets2 import MyThread, FindTargets
import cobra
import pandas as pd

if __name__ == '__main__':
    #modelo = FindTargets()
    #model = modelo.readModel("iMO1056_MOPS_glycerol_integrated.xml")
    #roda = MyThread("Thiago", "Pseudomonas aeruginosa", "merigueti@gmail.com", model, "1")
    #roda.run()

    ## ERRO NO MOMENTO DE GERAR O XLS DO MODELO.
    #modelo = FindTargets()
    #model = modelo.readModel("iPAO1_MOPS_acetate_integrated.xml")
    #roda = MyThread("Thiago", "Pseudomonas aeruginosa", "merigueti@gmail.com", model, "1")
    #roda.run()

    modelo = FindTargets()
    model = modelo.readModel("iPAO1_MOPS_glycerol_integrated.xml")
    roda = MyThread("Thiago", "Pseudomonas aeruginosa", "merigueti@gmail.com", model, "1")
    roda.run()

    '''
    modelo = FindTargets()
    model = modelo.readModel("MODEL1507180070_url.xml")
    roda = MyThread("Thiago", "Staphylococcus aureus", "merigueti@gmail.com", model, "1")
    roda.run()
    
    modelo = FindTargets()
    model = modelo.readModel("MODEL1507180070_url.xml")
    roda = MyThread("Thiago", "Staphylococcus aureus", "merigueti@gmail.com", model, "2")
    roda.run()
    
    modelo = FindTargets()
    model = modelo.readModel("MODEL1507180021.xml")
    roda = MyThread("Thiago", "Mycobacterium tuberculosis", "merigueti@gmail.com", model, "1")
    roda.run()
    
    modelo = FindTargets()
    model = modelo.readModel("MODEL1507180021.xml")
    roda = MyThread("Thiago", "Mycobacterium tuberculosis", "merigueti@gmail.com", model, "2")
    roda.run()
    
    modelo = FindTargets()
    model = modelo.readModel("MODEL1507180044.xml")
    roda = MyThread("Thiago", "Pseudomonas putida", "merigueti@gmail.com", model, "1")
    roda.run()
    
    modelo = FindTargets()
    model = modelo.readModel("MODEL1507180044.xml")
    roda = MyThread("Thiago", "Pseudomonas putida", "merigueti@gmail.com", model, "2")
    roda.run()
    '''