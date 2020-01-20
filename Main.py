#!/usr/bin/env ipython
###################### MAIN ######################

import os # Módulo funciones del sistema operativo.
import sys #Módulo para variables y funcionalidades.
import shutil #Módulo para mover archivos y carpetas.

from Bio import SeqIO #Modulo para parsear el GB
from subprocess import Popen, PIPE #Modulo para Blast.
from Bio.ExPASy import Prosite,Prodoc #Módulo para Prosite.
import re #Modulo para búsqueda de patrones.

#Mis módulos
import Blast_module as blast
import Muscle_module as muscle
import Prosite_module as prosite
import Intro_module as intro

#Intro
Intro_ctrl=intro.Intro()
if Intro_ctrl==True:

    #Parser
    Parser_ctrl,GB_save,Blast_save,Filtro_save,Align_save,Tree_save,Prosite_save,Pattern_save = blast.Parser() #Llama a la función parser
    if Parser_ctrl==True: #Si parser devuelve True (todo ha ido bien en el parseo), avanza a Blast. 
    
        ##Blast
        Blast_save[1],Blast_save[0],Blast_ctrl = blast.Blastp("Parsed_GB.fasta", "blast_default_name.fasta") #Llama a la función Blastp
        if Blast_ctrl==True: #Si el control devuelve True todo ha ido bien, avanza al filtrado. 
            
            #Filtro
            Filtro_save[1],Filtro_ctrl,Filtro_save[0] = blast.Filtro(Blast_save[1]) #Llamada a la función filtro.
            if Filtro_ctrl==True: #Si todo ha ido bien pasamos a muscle y prosite
    
                #Muscle 
                Align_save[0],Tree_save,Align_save[1]= muscle.align_and_make_tree(Filtro_save[1]) #Llama a Muscle
                
                ##Prosite
                Prosite_save,Pattern_save= prosite.Find_patterns(Filtro_save[1]) #Llama a prosite
                
        #Savingfiles
        intro.Save_files(GB_save,Blast_save,Filtro_save,Align_save,Tree_save,Prosite_save,Pattern_save) #Función para guardar los archivos. 
