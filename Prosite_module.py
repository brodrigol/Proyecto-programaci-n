#!/usr/bin/env ipython
###################### PROSITE ######################
import re
import os
from Bio.ExPASy import Prosite,Prodoc


def Find_patterns (Querys_list): #Función de búsqueda de patrones.
    print ("\n\nSe va a proceder a la Búsqueda de patrones en sus secuencias mediante prosite.")
    input ("\nPor favor, asegurese de que su carpeta contiene un archivo denominado 'Prosite.dat antes de continuar.\nPara continuar pulse INTRO.")
    count=0
    while count<3: #Bucle para asegurar la presencia de un archivo Prosite.dat. El usuario tiene 3 oportunidades antes de que se cierre el programa.
        if os.path.isfile("Prosite.dat"):
            count=3
        elif count<2:
            count=count+1
            input("No se encuentra el archivo 'Prosite.dat'. Por favor asegurese de que está en la misma carpeta que el programa.")
        else:
            count=count+1
            print("\n\nHa excedido el número de oportunidades. Se va a proceder a salir del programa\n\n\t\tHASTA PRONTO!!!")
            return ("n","n")
    print ("\n\nParseando Prosite.dat ...")
    handle = open("prosite.dat","r")
    records = Prosite.parse(handle)
    Round_text=[]
    Search_pattern=[]
    #Prosite_pattern
    for record in records: #Se abre el archivo y se parsea. Obtenemos name, accesion, description y pattern.
        Round=[]
        Name = record.name
        Round.append(Name)
        Accesion = record.accession
        Round.append(Accesion)
        Description = record.description
        Round.append (Description)
        Pattern = record.pattern
        Round.append (Pattern)
        R="\t".join(Round)
        Round_text.append(R)
        
        ## Sustituir la sintaxis de Prosite para adaptarla al módulo RE.
        Round[3]=Round[3].replace(".","")
        Round[3]=Round[3].replace("x",".")
        Round[3]=Round[3].replace("{","[^")
        Round[3]=Round[3].replace("}","]")
        Round[3]=Round[3].replace("(","{")
        Round[3]=Round[3].replace(")","}")
        Round[3]=Round[3].replace("<","^")
        Round[3]=Round[3].replace(">","$")
        Round[3]=Round[3].replace("-","")
        Search_pattern.append(Round)
    
    ##Archivo Prosite. Genera el archivo parseado. Los patrones se mantienen en la sintaxis de prosite.
    handle.close()
    File="\n".join(Round_text)
    prosite_file = open("Parsed_Prosite","w")
    prosite_file.write("Name\tAccesion\tDescription\tPattern\n")
    prosite_file.write(File)
    prosite_file.close()
    print("\n\nSe ha generado un archivo parseado a partir de su archivo.dat y se ha guardado como 'Parsed_prosite'")
    c=True
    while c==True: #Se pregunta al usuario si quiere guardar el archivo parseado.
        Prosite_save=input('¿Desea conservar el archivo "Parsed_Prosite" al terminar de ejecutar el programa? [s/n]:\t')
        if Prosite_save.lower()=="s":
            c=False
        elif Prosite_save.lower()=="n":
            c=False
        else:
            print("La opción introducida no es válida.")
    
    ##Buscar patrones
    print("\n\nSe va a proceder a la búsqueda de patrones.\n\n\tBuscando patrones en sus secuencias ...")
    for i in Querys_list: #Abre todos los archivos procedentes del filtrado.
        Final_file=[]
        file=("blast_"+i+"_result.fasta")
        S=open(file)
        Q=S.read()
        S.close()
        c=False        
        for j in Q.split("\n"):
            if c==False: #j=Nombre del subject
                name=j[1:]
                c=True
            else: #j=Secuencia
                for k in range(len(Search_pattern)): #Compara la secuencia del query con todos los patrones de prosite.
                    positions=[]
                    if len(Search_pattern[k][3])>0:
                        if re.search(Search_pattern[k][3],j):
                            for match in re.finditer(Search_pattern[k][3],j): #Bucle para obtener la posición de inicio y final y el numero de matches de un mismo patrón.
                                start=match.start()
                                end=match.end()
                                positions.append("["+str(start)+":"+str(end)+"]")
                            p=" ".join(positions)
                            Final_file.append(name+"\t"+Round_text[k]+"\t"+p+"\t"+str(len(positions))) #Si encuentra un match lo añade al archivo final.
                c=False
        Fin="\n".join(Final_file)
        Final=open(i+"_patterns.txt","w")
        Final.write("Subject_ID\tDomain Name\tDomain Accesion\tDomain Description\tDomain Pattern\tDomain Positions\tNumber of matches\n")
        Final.write(Fin) #Se genera el archivo con todos los matches para cada query.
        Final.close()
    print("\n\nSe ha/n obtenido",len(Querys_list),"archivo/s con los patrones encontrados para sus secuencias.")
    c=True
    while c==True: #Se pregunta al usuario si quiere guardar el archivo. 
        Pattern_save=input("¿Desea conservar el/los archivo/s al terminar de ejecutar el programa? [s/n]:\t")
        if Pattern_save.lower()=="s":
            c=False
        elif Pattern_save.lower()=="n":
            c=False
        else:
            print("La opción introducida no es válida.")
    return(Prosite_save.lower(),Pattern_save.lower())
                
                
    
    

    
    