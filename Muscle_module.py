#!/usr/bin/env ipython
###################### MUSCLE ######################

import sys
from subprocess import Popen, PIPE

## Cambiar para file=querys files. Todo el bucle que lo haga en la funcion.
            

def align_and_make_tree(querys_files): #Funcion de alineamiento y generado de arboles.
    tree_files=[]
    print("\nSe va a realizar el alineamiento de sus secuencias y generar los árboles correspondientes.\n\n")
    
    for i in querys_files: #Comprueba que los archivos tengan más de una secuencia para poder alinearlos y generar arboles.
        sys.stdout.write("\r\tComprobando archivo/s ...")
        file = ("blast_"+ i + "_result.fasta")
        F=open(file)
        f=F.read()
        F.close()
        file_content=list(f.split("\n"))
        sys.stdout.flush()
        if len(file_content)>2: #Si tienen más de una seq:
            tree_files.append(i)
            
            #Alignment
            sys.stdout.write("\r\tAlineando ...            ")
            proceso = Popen(['muscle','-in',file], stdout=PIPE, stderr=PIPE) #Alineado de secuencias.
            listado = proceso.stdout.read().decode("utf-8")
            #proceso.stderr.close()
            proceso.stdout.close()
            my_output = open(i+".fasta","w") #Creación del archivo alineado
            my_output.write(listado)
            my_output.close()
            
            #Make tree
            sys.stdout.flush()
            sys.stdout.write("\r\tCreando árbol ...        ")
            proceso = Popen(['muscle','-maketree','-in',i+".fasta",'-out', i+".tre",'-cluster','neighborjoining'], stdout=PIPE, stderr=PIPE) #Creación del árbol.
            #proceso.stderr.close()
            proceso.stdout.close()
            sys.stdout.flush()
            
    if len(tree_files)>0: #Indica cuantos archivos se han obtenido.
        print("\n\n\nSe ha realizado el alineamiento y generado el árbol NJ de",len(tree_files),"querys con éxito.")
        c=True
        while c==True: #Bucle para saber si el usuario quiere guardar los archivos de alineamiento.
            Align_save=input("¿Desea conservar el/los archivo/s del alineamiento al terminar de ejecutar el programa? [s/n]:\t")
            if Align_save.lower()=="s":
                c=False
            elif Align_save.lower()=="n":
                c=False
            else:
                print("La opción introducida no es válida.")
        c=True
        while c==True: #Bucle para saber si el usuario quiere guardar los árboles.
            Tree_save=input("¿Desea conservar el/los archivo/s de los árboles al terminar de ejecutar el programa? [s/n]:\t")
            if Tree_save.lower()=="s":
                c=False
            elif Tree_save.lower()=="n":
                c=False
            else:
                print("La opción introducida no es válida.")
        return(Align_save.lower(),Tree_save.lower(),tree_files)
    else:#Si los archivos solo contenían una secuencia.
        print("\n\n\nEl/Los archivo/s query introducido/s solo contenía/n una única secuencia, por lo que no se ha podido realizar el alineamiento o el árbol correspondiente.")
        return("n","n",tree_files)