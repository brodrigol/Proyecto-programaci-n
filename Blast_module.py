#!/usr/bin/env ipython
###################### BLAST-P ######################

import os
from Bio import SeqIO
from subprocess import Popen, PIPE #Modulo para Blast

def Parser(): #Función para el parseado del Genbank
    try:
        count=0
        while count<3: #Bucle para asegurar la existencia del archivo Genbank. El usuario tiene 3 oportunidades antes de que se cierre el programa.
            Genbank = input('Por favor, introduzca el nombre del archivo que desee utilizar como Genbank:\t') ##Input del usuario con el archivo Genbank
            if os.path.isfile(Genbank):
                count=3
            elif count<2:
                count=count+1
                print("No se encuentra el archivo",Genbank, ". Por favor asegurese de que está en la misma carpeta que el programa.")
            else:
                count=count+1
                print("\n\nHa excedido el número de oportunidades. Se va a proceder a salir del programa\n\n\t\tHASTA PRONTO!!!\n\n")
                return (False,"n",["n",""],["n",""],["n",""],"n","n","n")
        if os.path.isfile("Parsed_GB.fasta")==True: #Si existe un archivo que se llame así lo elimina porque es el nombre que vamos a utilizar para guardar el genbank parseado.
            os.remove("Parsed_GB.fasta")
        with open(Genbank,"r") as input_handle: #Funcion de parseo
            print("\n\n\tParseando ...\n")
            for record in SeqIO.parse(input_handle,"genbank"):
                print()
            for feature in record.features:
                if feature.type=='CDS':
                    locus=feature.qualifiers['locus_tag'][0]
                    seq=feature.qualifiers['translation'][0]
                    fasta=str(">"+locus+"\n"+seq+"\n")
                    GB=open("Parsed_GB.fasta","a")
                    GB.write(fasta)
        input_handle.close()
        GB.close()
        print('Su Genbank parseado se ha guardado en el archivo "Parsed_GB.fasta".')
        c=True
        while c==True: #Bucle para saber si el usuario quiere guardar el archivo parseado.
            GB_save=input('¿Desea conservar el archivo "Parsed_GB.fasta" al terminar de ejecutar el programa? [s/n]:\t')
            if GB_save.lower()=="s":
                c=False
            elif GB_save.lower()=="n":
                pass
                c=False
            else:
                print("La opción introducida no es válida.")
        return(True,GB_save.lower(),["n",""],["n",""],["n",""],"n","n","n")
    except UnboundLocalError: #En caso de que el archivo no tenga el formato correcto se imprime el siguiente mensaje.
        print ('\nEl archivo "',Genbank,'", que ha introducido como Genbank no es válido.')
        print('Se va a proceder a salir del programa.\n\n\t\tHASTA PRONTO!!')
        return (False,"n",["n",""],["n",""],["n",""],"n","n","n")

### BLAST-P ###                 
def Query (): #Función para introducir el nombre de la query. Asegura que el archivo exista y sea un .fasta
    query_name=input('\nPor favor, introduzca el nombre de su archivo query:\t')
    if os.path.isfile(query_name)==False:
        print('\n\nEl archivo "',query_name,'", que ha introducido como query no existe.', 
              '\nPor favor, asegurese de que ha introducido el nombre correctamente y de que el archivo se encuentra en la misma carpeta que este programa.')
        Volver=input("Si desea volver a introducir un archivo query escriba 'volver', si desea salir del programa pulse intro.\t")
        if Volver.lower()== "volver":
            return(True)
        else:
            print ("Cerrando programa ...\n\n\t\tHasta pronto!")
            return (False)
    else: #Comprueba el formato fasta
        q=open(query_name)
        Q=q.read()
        Q=list(Q.split("\n"))
        q.close()
        b=0; fasta_ctrl=True
        while b<len(Q):
            if Q[b]=="":
                b=b+1
            elif Q[b][0]==">":
                b=len(Q)
            else:
                fasta_ctrl=False
                b=len(Q)
        if fasta_ctrl==False:
            print ('\n\nEl archivo introducido como Query no se encuentra en formato fasta.\t')
            Volver=input("Si desea volver a introducir un archivo query escriba 'volver', si desea salir del programa pulse intro.")
            if Volver.lower()== "volver":
                return(True)
            else:
                print ("Cerrando programa ...\n\n\t\tHASTA PRONTO!!")
                return (False)
        else:
            return(query_name)
        
def Blastp (subject, Blast_name):
    print ("\n\nSe va a proceder a realizar el Blast-p")
    query=True #Bucle para obtener el nombre de la query
    while query==True:
        query=Query()
    if query==False:
        return(Blast_name,"n",False)
    else:
        print ("\n\n\tBlasting ...\n\n")
        #Blast del archivo query frente al genbank.
        proceso = Popen(['blastp','-query', query ,'-subject', subject,'-outfmt',"6 qseqid sseqid qcovs pident evalue" ], stdout=PIPE, stderr=PIPE)
        error_encontrado = proceso.stderr.read().decode("utf-8")
        
        proceso.stderr.close()
        listado = proceso.stdout.read().decode("utf-8")
        proceso.stdout.close()
        
        my_output = open(Blast_name,"w")
        my_output.write(listado)
        
        if error_encontrado: 
            print("Se produjo el siguiente error:\n%s" % error_encontrado)#Mensaje de error en caso de que algo falle
            return(Blast_name,"n",False) 
        else:
            print ("Se ha realizado el Blast-p correctamente.")
            c=True #Pregunta al usuario si quiere guardar el archivo blast y que nombre quiere darle.
            while c==True:
                blast_save=input("¿Desea conservar el archivo Blast al terminar de ejecutar el programa? [s/n]:\t")
                if blast_save.lower()=="s":
                    Blast_name=input('Por favor, introduzca el nombre que quiera dar a su archivo blast:\t')
                    Blast_name=Blast_name+".fasta"
                    os.rename("blast_default_name.fasta",Blast_name)
                    c=False
                elif blast_save.lower()=="n":
                    c=False
                else:
                    print("La opción introducida no es válida.")
            return (Blast_name,blast_save.lower(),True)


### Control de errores FLOAT INPUT ### 
def floatval (mensaje, lim_min, lim_max): #Funcion para control de errores si no se introduce un float
    def error(mensaje):
        try:
            Value=float(input(mensaje))
        except:
            return ('Mal')
        else:
            return (Value)
    ctrl=True
    while ctrl==True:
        c=error(mensaje)
        if c=='Mal':
            print ('El valor introducido no es válido')
        elif lim_min=='No' and c<=lim_max:
            ctrl=False    
        elif lim_max=='No' and c>=lim_min:
            ctrl=False
        elif c>=lim_min and c<=lim_max:
            ctrl=False
        else:
            print('El valor introducido no es válido')
    return (c)
        
### FILTRO ###        
def Filtro (blast_file):
    try:
        print ("\n\nSe va a proceder a filtrar los resultados del Blast-p")
        c1=False
        while c1==False: #Bucle por si al final los valores introducidos no son superados por ninguna secuencia.
            #Introducción de los valores de filtrado.
            Cov=floatval('\nPor favor, introduzca el valor de Coverage con el que desea filtrar los resultados (se tomarán los valores mayores o iguales al introducido):\t',0,100)
            Iden=floatval('Por favor, introduzca el valor de % de identidad con el que desea filtrar los resultados (se tomarán los valores mayores o iguales al introducido):\t',0,100)
            Eval=floatval('Por favor, introduzca el valor de e value con el que desea filtrar los resultados (se tomarán los valores menores o iguales al introducido):\t', 0, 'No')
            
            #Abre el archivo blast y genera una lista.
            filtro = open(blast_file)
            F=filtro.read()
            Lista = list (F.split("\n"))
            filtro.close()
            List=[]
            for i in range(len(Lista)):
                if Lista[i]!="":     
                    L=list(Lista[i].split("\t"))
                    List.append(L)
            
            Querys=[] #Querys será una lista de los nombres de las distintas querys. 
            for i in range(len(List)):
                if List[i][0] in Querys:
                    pass
                else:    
                    Querys.append(List[i][0])
                    
            S=open("Parsed_GB.fasta")
            Subject=S.read()
            S.close()
            print("\n\n\tFiltrando...\n")
            Querys_file_names=[] #Lista las querys que han superado el filtrado y de las que se genera archivo.
            for i in range(len(Querys)):
                Sub=[]
                ctrl=True #Control para indicar cuando hay coincidencia entre ID de subject del blast y el genbank
                for j in range(len(List)):
                    if Querys[i]==List[j][0] and float(List[j][2])>=Cov and float(List[j][3])>=Iden and float(List[j][4])<=Eval: #Condición de filtro
                        Sub.append(">"+List[j][1])
                        for k in Subject.split("\n"):
                            if ctrl==False:#El ID coincide, se guarda el ID.
                                Sub.append(k)
                                ctrl=True #En la siguiente ronda se va a introducir la seq.
                            else:
                                if k[1:]==List[j][1]:
                                    ctrl=False   #Se introduce la secuencia.
                if len(Sub)>0: #Si hay al menos una secuencia se genera archivo.
                    Todo="\n".join(Sub)
                    Querys_file_names.append(Querys[i])
                    
                    N=['blast_',Querys[i],'_result.fasta']
                    Name="".join(N)
                    Filtrado=open(Name,"w")
                    Filtrado.write(Todo)
                    Filtrado.close()
            if len(Querys_file_names)>0:#Se ha generado al menos un archivo     
                print("\nEl filtrado se ha realizado con exito. Se han generado", len(Querys_file_names),"archivo/s.")
                c1=True
                c=True
                while c==True:#Bucle para saber si se quieren guardar los archivos.
                    Filtro_save=input("¿Desea conservar los archivos resultantes del filtrado al terminar de ejecutar el programa? [s/n]:\t")
                    if Filtro_save.lower()=="s":
                        c=False
                    elif Filtro_save.lower()=="n":
                        c=False
                    else:
                        print("La opción introducida no es válida.")
                return (Querys_file_names,True,Filtro_save.lower())
            else: #No se ha obtenido ningun archivo
                print("Ninguna de las secuencias introducidas cumple con las condiciones impuestas.")
                
                #Se da al usuario la oportunidad de volver a filtrar o se cierra el programa.
                option=input("Si desea introducir nuevas condiciones escriba 'filtrar', si desea salir del programa pulse INTRO.\t")
                if option.lower()=="filtrar":
                    pass
                else:
                    c1=True
                    print("\nEstá usted saliendo del programa\n\n\t\tHASTA PRONTO!!!")
                    return(Querys_file_names,False,"n")
                   
    except FileNotFoundError: #Si no se encuentra el archivo genbank o el blast salta el siquiente mensaje de error.
        if os.path.isfile("Parsed_GB.fasta"):
            print ("\n\nNo se han podido generar los archivos blast.fasta filtrados correspondientes a cada organismo.\n",
                   "Por favor, compruebe que el archivo '", blast_file,"' se encuentra en la misma carpeta que este programa")
        else:
            print ("\n\nNo se han podido generar los archivos blast.fasta filtrados correspondientes a cada organismo.",
                   "\nPor favor, compruebe que el archivo ' Parsed_GB.fasta ' se encuentra en la misma carpeta que este programa")
                   
        return([],False,"n")


                







