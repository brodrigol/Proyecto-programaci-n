#!/usr/bin/env ipython
###################### INTRO ######################
import os
import shutil

def Intro(): #Mensaje de Bienvenida al usuario
    print ("""             
                BIENVENIDO AL PROGRAMA "PROYECTO FINAL" 
                  scripted by Blanca Rodrigo Lacave

OPCIONES
    a. Entrar al programa. 
    b. Ayuda para el usuario.   
    c. Salir           
""")
    #Bucle para seleccionar una de las opciones que da el menú. 
    c=True
    while c==True:
        opcion=input("Por favor, seleccione una de las opciones anteriores. [a,b,c]:\n")
        if opcion.lower()=="a":
            print ("Iniciando el programa ...\n\n")
            print("Se va a parsear su Genbank")
            c=False
            return(True)
        elif opcion.lower()=="b":
            Help()
            c=False
        elif opcion.lower()=="c":
            c=False
            print("\n\nEstá usted saliendo del programa ...\n\n\t\tHASTA PRONTO!!!\n\n")
            return (False)
        else:
            print("La opción introducida no es válida.")
    


def Help(): #Función guía para el usuario.
    print ("""
BIENVENIDO A LA AYUDA.

El funcionamiento de este programa es muy sencillo. Consta de cuatro módulos cuyo funcionamiento se explica a continuación:
    
    1. Módulo Blast.
       
       Este primer módulo cuenta con tres funciones básicas.
        a) Genbank_Parser. El usuario debe proporcionar un archivo .gbff para que el programa lo parsee. 
           Si el archivo no es válido o no existe saltará un mensaje de error.
        b) Blast-p. A partir del archivo Genbank parseado y de un query en formato fasta o multifasta 
           que el usuario debe proporcionar se realiza un Blast-p. 
           El resultado se almacena en un archivo .fasta.
        c) Filtrado. Una vez obtenido el archivo Blast, se brinda al usuario la posibilidad de filtrar 
           los resultados obtenidos por cobertura, % identidad y e-value. 
           Las secuencias subject del blast que cumplan estos requisitos se almacenarán en archivos 
           .fasta clasificados según el organismo query. 
           Si ninguno de los resultados de blast cumple con las condiciones impuestas por el usuario se 
           imprimirá un mensaje en pantalla que lo indique, dandole la posibilidad de reajustar los 
           parámetros o de salir del programa.
    
    2. Módulo Muscle.
       
       En este segundo módulo se van a alinear las secuencias para la posterior obtención de árboles 
       filógenéticos. Cuenta con una única función.
        a) Align_and_make_tree. Para el funcionamiento de este módulo el usuario no debe aportar nada. 
           Simplemente dejarlo ejecutar. 
           Por favor, no elimine o cambie de carpeta los archivos generados en el módulo anterior, ya 
           que sin ellos no podrá funcionar. 
           En primer lugar se comprueba que el archivo de las querys dadas contengan más de una secuencia,
           para garantizar que pueden alinearse y generar árboles. 
           Seguidamente se procede al alineamiento y obtención de los árboles filogenéticos.
    
    3. Módulo Prosite.
        
       Este módulo sirve para buscar patrones o dominios conocidos en sus secuencias query.
        a) Find_patterns. En primer lugar debe disponer de un archivo Prosite.dat en la carpeta en la que
           se encuentre este programa. Si no es así saltará un mensaje de error.
           A partir de un parseado del archivo Prosite.dat se generará un nuevo archivo denominado Parsed_Prosite. 
           Con este archivo y los archivos .fasta obtenidos en el primer módulo se obtienen matches en ciertos 
           dominios de cada secuencia. Estos matches se almacenarán en un archivo correspondiente a cada query, en
           los que se indica el subject, el ID del dominio,la posición en la que se encuentra...
    
    4. Main. 
       
       Es el script principal a partir del cual se llama al resto de módulos y funciones del programa. 
    
    5. Módulo Intro
       
       Simplemente contiene las funciones de Bienvenida al programa, ayuda y guardado de archivos. 
       Esta última es la que permite que se almacenen sus resultados, si es que quiere guardarlos, en carpetas organizadas. 
       Al final del programa se le indicará donde puede encontrar sus archivos.
      
Para salir de la guía pulse INTRO.
""")
    input("")
    Intro() #Vuelve a llamar a intro para que el usuario decida si usar el programa o salir del mismo. 

def Save_files (GB,blast,filtro,alineamiento,arbol,prosite,patrones): #Función para guardar o eliminar los archivos según haya indicado el usuario a lo largo del programa.
    print ("\n\nEl programa ha terminado de generar archivos.")
    ##Caso de que el usuario no quiera guardar nada
    if GB=="n" and blast[0] == "n" and filtro[0] == "n" and alineamiento[0]=="n" and arbol=="n" and prosite == "n" and patrones == "n":
        print ("\n\n\tEliminando archivos ...")
        #Comprueba la existencia de los archivos y si existen los elimina.
        if os.path.isfile("Parsed_GB.fasta"):
            os.remove("Parsed_GB.fasta")
        if os.path.isfile(blast[1]):
            os.remove("Blast_default_name.fasta")
        for i in filtro[1]:
            f="blast_"+i+"_result.fasta"
            if os.path.isfile(f):
                os.remove(f)
            g=i+"_patterns.txt"
            if os.path.isfile(g):
                os.remove(g)
        for i in alineamiento[1]:
            f=i+".fasta"
            if os.path.isfile(f):
                os.remove(f) 
            g=i+".tre"
            if os.path.isfile(g):
                os.remove(g)
        if os.path.isfile("Parsed_Prosite"):
            os.remove("Parsed_Prosite")
        
        print("\n\nSe han eliminado todos los archivos generados por este programa.\n")
    
    ##Caso de que el usuario quiera guardar al menos un archivo
    else:
        print("\n\n\tGenerando carpetas y organizando los resultados ...")
        
        #Genera una carpeta para almacenar todas las carpetas de resultados que se puedan generar.
        if os.path.isdir("Results")==False:
            os.mkdir("Results")
            os.mkdir("Results_0")
            j=0
        else: 
            j=1;ctrl=False; lista=os.listdir("Results")
            while ctrl==False: #Bucle para obtener la carpeta Results_(n).
                if "Results_"+str(j) in lista:
                    j=j+1
                else:
                    ctrl=True
                    os.mkdir("Results_"+str(j))
        
        #Si el usuario indicó que no quería guardar el archivo se elimina. 
        #Si indicó que quería conservarlo se mueve a la carpeta Results_n.
        #Si hay varios archivos del mismo tipo se crea una subcarpeta para almacenarlos dentro de Results_n
        if GB =="n" and os.path.isfile("Parsed_GB.fasta"):
            os.remove("Parsed_GB.fasta")
        elif GB=="s" and os.path.isfile("Parsed_GB.fasta"):
            shutil.move("Parsed_GB.fasta","Results_"+str(j))
            
        if blast[0] == "n" and os.path.isfile(blast[1]):
            os.remove("Blast_default_name.fasta")
        elif blast[0] == "s" and os.path.isfile(blast[1]):
            shutil.move(blast[1],"Results_"+str(j))
        
        if filtro[0] == "n":
            for i in filtro[1]:
                f="blast_"+i+"_result.fasta"
                if os.path.isfile(f):
                    os.remove(f)
        else:
            os.mkdir("Filtered_sequences")
            for i in filtro[1]:
                f="blast_"+i+"_result.fasta"
                shutil.move(f,"Filtered_sequences")
            shutil.move("Filtered_sequences","Results_"+str(j))
        
        if alineamiento[0]=="n":
            for i in alineamiento[1]:
                f=i+".fasta"
                if os.path.isfile(f):
                    os.remove(f)
        else:
            os.mkdir("Aligned_sequences")
            for i in alineamiento[1]:
                f=i+".fasta"
                shutil.move(f,"Aligned_sequences")
            shutil.move("Aligned_sequences","Results_"+str(j))
            
        if arbol=="n":
            for i in alineamiento[1]:
                f=i+".tre"
                if os.path.isfile(f):
                    os.remove(f)
        else:
            os.mkdir("Trees")
            for i in alineamiento[1]:
                f=i+".tre"
                shutil.move(f,"Trees")
            shutil.move("Trees","Results_"+str(j))
                    
        if prosite == "n" and os.path.isfile("Parsed_Prosite"):
            os.remove("Parsed_Prosite")
        elif prosite == "s":
            shutil.move("Parsed_Prosite","Results_"+str(j))
        
        if patrones == "n":
            for i in filtro[1]:
                f=i+"_patterns.txt"
                if os.path.isfile(f):
                    os.remove(f)
        else:
            os.mkdir("Patterns")
            for i in filtro[1]:
                f=i+"_patterns.txt"
                shutil.move(f,"Patterns")
            shutil.move("Patterns","Results_"+str(j))
        
        shutil.move("Results_"+str(j),"Results") #La carpeta Results_n se almacena en Results.
        print("\nSus resultados se han guardado correctamente en la carpeta","Results_"+str(j),"dentro del directorio Results.")

    print("\n\n\t\tEl programa ha finalizado. Hasta pronto!!\n\n")
