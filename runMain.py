#!/usr/bin/env python

import os,sys
<<<<<<< HEAD

=======
import subprocess
>>>>>>> parent of 623b2ee (upadte)
import time

import commands

from module.configure import getserver



docstring='''

runMain.py jobID

    run the full build_MSA pipeline:

'''
<<<<<<< HEAD



=======
"""
new:
    need add argv for prefunc
    need add argv for homoflag real or benchmark
"""
prefunc = 'false'
homoflag = 'real'
>>>>>>> parent of 623b2ee (upadte)
####  parse command line argument ####

if len(sys.argv)<2:

    sys.stderr.write(docstring)

    exit()

s=sys.argv[1] # protein name (jobID)

<<<<<<< HEAD
homoflag='real' # real or benchmark

=======
>>>>>>> parent of 623b2ee (upadte)
Q="casp"

heteromer_maxP=10

server="S10" # S10, amino, GL

server=getserver()

print("server=%s"%server)

#### parse directory structure ####

bindir=os.path.dirname(os.path.abspath(__file__)) #where script and programs are

datarootdir=os.path.join(os.path.dirname(bindir),"output")

datadir=os.path.join(datarootdir,s) #directory for all jobs

#datadir_amino=os.path.join(os.path.dirname(bindir),"output_amino",s)

#os.system("cp -r %s %s"%(datadir_amino,datadir))

#print(datadir)

os.system("chmod -R 777 %s &>/dev/null"%datadir)



os.chdir(datadir)



#### check input ####

if not os.path.isfile("seq.fasta") or not os.path.getsize("seq.fasta"):

    sys.stderr.write("ERROR! No sequence at %s/seq.fasta\n"%datadir)

    exit()





HHLIB=bindir+"/DeepMSAFold/DeepMSA2/bin/DeepMSA"
<<<<<<< HEAD

=======
runCOFACTOR=os.path.join(bindir,"runCOFACTOR.py")
>>>>>>> parent of 623b2ee (upadte)
AF2Mdir=os.path.join(bindir,"DeepMSAFold/alphafold_multi_2.2.0")

protein_type="multimer"

###### check configure file

msa_gene_flag="deepmsa"

config_file=open(os.path.join(datadir,'config'),'r')

lines=config_file.readlines()

config_file.close()

msa_gene_flag=lines[0].strip('\n')
<<<<<<< HEAD

# DMFOLD2

prefunc_flag=lines[1].strip('\n')

prefunc_flag = (prefunc_flag.lower() == "true")

# DMFOLD2
=======
>>>>>>> parent of 623b2ee (upadte)



######### parse sequence

seqfile=open(os.path.join(datadir,"seq.fasta"),"r")

seqtxt=seqfile.read()

seqfile.close()

seqblocks=('\n'+seqtxt).split('\n>')[1:]



if len(seqblocks)>1:

    protein_type="multimer"

else:

    protein_type="monomer"



print("protein type=%s"%protein_type)



if protein_type=="multimer":

    os.system("%s/run_alphafold_mymonomermsa_parsechain.sh -o %s -f %s -t 2099-01-01 -g false -r false -m multimer -c full_dbs -p false -l 5 "%(AF2Mdir,os.path.join(datadir,"AlphaFold2"),os.path.join(datadir,"seq.fasta")))

    os.system("%s/DeepMSAFold/parse-sequence.py %s %s %s %s"%(bindir,os.path.join(datadir,"AlphaFold2/seq/msas/chain_id_map.json"),s,os.path.join(datadir,"chain.list"),datadir))

else:

    os.system("mkdir -p %s/%s"%(datadir,s+'-A'))

    os.system("cp %s/seq.fasta %s/%s/"%(datadir,datadir,s+"-A/"))

    chainfile=open(os.path.join(datadir,'chain.list'),'w')

    chainfile.write(s+'-A\n')

    chainfile.close()



os.system("chmod -R 777 %s &>/dev/null"%datadir)



####### run DeepMSA_noIMG

qmsa="no"

if(msa_gene_flag=="deepmsa"):

    qmsa="no"

else:

    qmsa="yes"



chains=[]

chain_file=open(os.path.join(datadir,'chain.list'),'r')

lines=chain_file.readlines()

chain_file.close()

for line in lines:

    chains.append(line.strip('\n'))



print(chains)



#exit(1)



noimg_complete=0

while 1:

    noimg_complete=1

    print("run DeepMSA2 dMSA/qMSA searching with %s"%s)

    os.system("%s/DeepMSAFold/DeepMSA2/bin/DeepMSA2_noIMG.pl %s %s %s %s"%(bindir,s,datarootdir,qmsa,server))

    if qmsa=="no":

        for chain in chains:

            if not os.path.exists(os.path.join(datadir,chain,'MSA/DeepMSA.aln')):

                noimg_complete=0

    if qmsa=="yes":

        for chain in chains:

            if (not os.path.exists(os.path.join(datadir,chain,'MSA/DeepMSA.aln'))) or (not os.path.exists(os.path.join(datadir,chain,'MSA/qMSA.aln'))) or (not os.path.exists(os.path.join(datadir,chain,'MSA/AlphaFold2.aln'))):

                noimg_complete=0

    if noimg_complete:

        break

    time.sleep(1200)



#exit(1)



if msa_gene_flag=="deepmsa2img":

    #img_complete=0

    

    while 1:

        img_complete=0   

        print("run DeepMSA2 IMG searching with %s"%s)

        (status,output)=commands.getstatusoutput("%s/DeepMSAFold/DeepMSA2/bin/DeepMSA2_IMG.pl %s %s %s"%(bindir,s,datarootdir,server))

        print(output)

        if output.__contains__("JGI search still running"):

            img_complete=0

        for chain in chains:

            if output.__contains__("%s does not need additional JGI search"%chain) or output.__contains__("sJGI_"+chain+" finished"):

                img_complete+=1

        print(img_complete)

        if img_complete==len(chains):

            break

        time.sleep(1200)



#exit(1)





######## MSA ranking

while 1:

    ranking_complete=0

    os.system("%s/DeepMSAFold/DeepMSA2/bin/run_AF2_multiMSA_subM5.py %s %s %s"%(bindir,s,datarootdir,server))

    for chain in chains:

        if os.path.exists(os.path.join(datadir,chain,"seq.aln")) and os.path.exists(os.path.join(datadir,chain,"finalMSAs/MSA_ranking.info")):

            print("MSA selection for %s is complete!"%chain)

            ranking_complete+=1

    if ranking_complete==len(chains):

        break

    time.sleep(1200)



#exit(1)



######## combine MSAs and run complex modeling

if protein_type=="multimer":

    os.system("%s/DeepMSAFold/combine-MSAs_v3.py %s %s %s %s"%(bindir,s,os.path.join(datadir,"AlphaFold2/seq/msas/chain_id_map.json"),datarootdir,str(heteromer_maxP)))

    print("combine MSA complete!")

    while 1:

        (status,out)=commands.getstatusoutput("%s/DeepMSAFold/run_AF2_multiMSA_complex.py %s %s %s"%(bindir,s,datarootdir,server))

        print(out)

        if os.path.exists(os.path.join(datadir,"AF2COMP/model.info")) and (not out.__contains__("AF2COMP jobs are not complete yet!")):

            os.system("%s/DeepMSAFold/ranking_complex_allmodels_byQA.py %s %s"%(bindir,s,datarootdir))

            break

        time.sleep(1200)

else:

    os.system("mkdir -p %s/Final_ALL_AF2MODELS"%datadir)

    os.system("chmod -R 777 %s/Final_ALL_AF2MODELS"%datadir)

    os.system("cp -r %s/%s-A/AF2M5/* %s/Final_ALL_AF2MODELS/"%(datadir,s,datadir))

    

for i in range(1,6):

    os.system("cp %s %s"%(os.path.join(datadir,"Final_ALL_AF2MODELS/model_"+str(i)+".pdb"),datadir))

<<<<<<< HEAD


# DMFOLD2

# path set up

# for DMsearch

DMsearch = os.path.join(bindir, "DMsearch.py")

USalgin = os.path.join(bindir, 'DeepMSAFold', 'USalgin')

DMfold_dir = os.path.dirname(bindir)

libraray_dir = os.path.join(DMfold_dir, "library")

foldseek_db = os.path.join(libraray_dir, 'pdb_foldseek', 'PDB100')

target_dir = os.path.join(libraray_dirs, "pdb_mmcif_db")



# for function prediction

runCOFACTOR = os.path.join(bindir, "runCOFACTOR.py")

runCOFACTOR_multimer = os.path.join(bindir, "runCOFACTOR_multimer.py")

# new for function prediction and DMsearch

if protein_type == "multimer":



    # predict function

    if prefunc_flag:

        # run COFACTOR_multimer.py data_dir homoflag

        fun_cmd = [runCOFACTOR_multimer, datadir, homoflag]

        subprocess.run(fun_cmd, check=True)

    

    # run DMsearch.py

    input_model = os.path.join(datadir, "model_1.pdb")

    # set up directory

    search_dir = os.path.join(datadir, "DMsearch")

    subprocess.run(["mkdir", '-p', search_dir])

    # search 

    search_cmd = [DMsearch, 'search', input_model, target_dir, foldseek_db, search_dir]

    subprocess.run(search_cmd, check=True)



elif protein_type == "monomer":

    # monomer

    # read the 

    chain_list = open(os.path.join(datadir, "chain.list"), "r").readlines()

    chian_id = chain_list[0].strip().split("-")[-1]

    # set up directory

    cofactor_dir = os.path.join(datadir, '{}_{}_cofactor'.format(s, chian_id))

=======

# new for function prediction
if protein_type == "multimer":
    pass
else:
    # monomer
    # set up directory
    cofactor_dir = os.path.join(datadir, f"{s}-cofactor")
>>>>>>> parent of 623b2ee (upadte)
    subprocess.run(["mkdir", '-p', cofactor_dir])

    # cp seq.fasta

    subprocess.run(["cp", os.path.join(datadir, "seq.fasta"), cofactor_dir])

    # cp model_1.pdb

    model_for_cofactor = os.path.join(cofactor_dir, "model1.pdb")

    subprocess.run(["cp", os.path.join(datadir, "model_1.pdb"), model_for_cofactor])

    # start function prediction

    # runCOFACTOR.py data_dir tag homoflag
<<<<<<< HEAD

    cmd = [runCOFACTOR, cofactor_dir, s, homoflag]

    if not prefunc_flag:

        cmd.append("-s")

    subprocess.run(cmd)

# DMFOLD2
=======
    subprocess.run([runCOFACTOR, cofactor_dir, s, homoflag])


>>>>>>> parent of 623b2ee (upadte)



#exit(1)

#### [3] HTML webpage output ####

os.system("python3 %s/file2html/file2html.py %s"%(bindir,datadir))



os.system("chmod -R 777 %s &>/dev/null"%datadir)



#os.system("find -name *.pkl|xargs rm")

