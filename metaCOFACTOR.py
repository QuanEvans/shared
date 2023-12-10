#!/usr/bin/env python

docstring='''metaCOFACTOR.py datadir

    combine the output of COFACTOR, GOfreq, and PPI2GO



input files:

    GOsearchresult_*_{MF,BP}.dat - COFACTOR MF/BP final GO prediction

    COFACTOR_GOfreq_CC         - COFACTOR CC final GO prediction

    combine_gwGOfreq_{MF,BP,CC}  - GOfreq final GO prediction

    string_swGOfreq_{MF,BP,CC}   - STRING final GO prediction



output files

    dist_CofactorGOfreqPPI_{MF,BP,CC} - combine COFACTOR, GOfreq, STRING

    dist_CofactorGOfreq_{MF,BP,CC}    - combine COFACTOR, GOfreq

    dist_GOfreqPPI_{MF,BP,CC}         - combine GOfreq,   STRING

    dist_CofactorPPI_{MF,BP,CC}       - combine COFACTOR, STRING



options:

-excludeGO=GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575 

    GO term to be excluded

    (default) remove root of 3 aspect and "protein binding"

        GO:0005515 ! protein binding

        GO:0005488 ! binding

        GO:0003674 ! molecular_function

        GO:0008150 ! biological_process

        GO:0005575 ! cellular_component

'''



import sys,os

from module import obo2csv

from module.fetch import wget

from glob import glob



obo_url="http://geneontology.org/ontology/go-basic.obo"



def parse_COFACTOR_result(infile="GOsearchresult_MF.dat",obo_dict=dict(),

    excludeGO=''):

    '''parse COFACTOR format GO prediction file, back-tracing all parent 

    nodes. obo_dict is the GO hierachy object. excludeGO is a comma 

    seperated list of GO terms to be ignored'''

    if isinstance(excludeGO,str):

        excludeGO=excludeGO.split(',')

    excludeGO=set(excludeGO)



    COFACTOR_dict=dict()

    if not os.path.isfile(infile):

        return COFACTOR_dict



    fp=open(infile,'rU')

    txt=fp.read()

    fp.close()

    for line in txt.splitlines():

        line=line.split()

        GOterm=line[0]

        Cscore=float(line[1])

        if not GOterm in excludeGO:

            COFACTOR_dict[line[0]]=float(line[1])

    

    ## trace the parent terms ##

    for GOterm in COFACTOR_dict.keys():

        Cscore=COFACTOR_dict[GOterm]

        parent_term_list=obo_dict.is_a(GOterm,direct=False,

            name=False,number=False).split()

        for parent_term in list(set(parent_term_list)-excludeGO):

            if not parent_term in COFACTOR_dict or \

                Cscore>COFACTOR_dict[parent_term]:

                COFACTOR_dict[parent_term]=Cscore

    return COFACTOR_dict



def parse_GOfreq_result(infile="combine_gwGOfreq_MF",excludeGO=''):

    '''parse COFACTOR format GO prediction file, discarding all prediction 

    with 0 cscore. excludeGO is a comma seperated list of GO terms to be

    ignored'''

    if isinstance(excludeGO,str):

        excludeGO=excludeGO.split(',')



    GOfreq_dict=dict()

    if not os.path.isfile(infile):

        return GOfreq_dict



    fp=open(infile,'rU')

    txt=fp.read()

    fp.close()

    for line in txt.splitlines():

        GOterm,Aspect,Cscore=line.split()

        Cscore=float(Cscore)

        if Cscore and not GOterm in excludeGO and (not GOterm in GOfreq_dict \

            or (GOterm in GOfreq_dict and Cscore>GOfreq_dict[GOterm])):

            GOfreq_dict[GOterm]=Cscore

    return GOfreq_dict



def combine_predictor_dist(dict_list=[]):

    '''combine GO prediction by muliple predictors using distance:

        cscore=1-(1-cscore1)*(1-cscore2)*(1-cscore3)

    dict_list is a list of dict whose key is GOterm and value is Cscore'''

    combine_dict=dict()



    ## get a list of all predicted GO terms ##

    GOterm_list=[]

    for predictor_dict in dict_list:

        GOterm_list+=predictor_dict.keys()

    GOterm_list=list(set(GOterm_list))



    ## combine cscore ##

    for GOterm in GOterm_list:

        cscore=1.

        for predictor_dict in dict_list:

            if GOterm in predictor_dict:

                cscore*=(1.-predictor_dict[GOterm])

        cscore=1-cscore

        if float("%.2f"%cscore):

            combine_dict[GOterm]=cscore

    return combine_dict



def combine_seqID_weighted_predictor_dist(dict_list={

    'cofactor':dict(),'gofreq':dict(),'ppi':dict(),'seqID':dict()}):

    combine_dict=dict()



    ## get maximum seqID ##

    seqID=0.

    if dict_list['seqID']:

        seqID=max(dict_list['seqID'].values())



    ## get a list of all predicted GO terms ##

    GOterm_list=[]

    for predictor in ['cofactor','gofreq','ppi']:

        GOterm_list+=dict_list[predictor].keys()

    GOterm_list=list(set(GOterm_list))



    ## combine cscore ##

    for GOterm in GOterm_list:

        cscore=1.

        if GOterm in dict_list['cofactor']:

            cscore*=(1.-dict_list['cofactor'][GOterm])**(1-seqID)

        if GOterm in dict_list['ppi']:

            cscore*=(1.-dict_list['ppi'][GOterm])

        if GOterm in dict_list['gofreq']:

            cscore*=(1.-dict_list['gofreq'][GOterm])#**seqID

        cscore=1-cscore

        if float("%.2f"%cscore):

            combine_dict[GOterm]=cscore



    return combine_dict



def write_GOpred(combine_dict=dict(),Aspect='F',outfile="dist_GOfreqPPI_MF"):

    '''write GO predictions into outfile in GOfreq format. combine_dict is dict

    whose key is GOterm and value is Cscore'''

    txt=''

    GOterm_list=sorted([(combine_dict[GOterm],GOterm) for GOterm in \

        combine_dict if combine_dict[GOterm]>0], reverse=True)

    txt=''.join(["%s\t%s\t%.2f\n"%(GOterm,Aspect[-1],Cscore

        ) for Cscore,GOterm in GOterm_list])

    fp=open(outfile,'w')

    fp.write(txt)

    fp.close()

    return



if __name__=="__main__":    

    excludeGO="GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575"

    argv=[]

    for arg in sys.argv[1:]:

        if arg.startswith("-excludeGO="):

            excludeGO=arg[len("-excludeGO="):].upper()

        elif arg.startswith('-'):

            sys.stderr.write("ERROR! Unknown option")

        else:

            argv.append(arg)



    if len(argv)<1:

        sys.stderr.write(docstring)

        exit()



    #### parse GO hierachy ####

    fp=open(wget(obo_url,show_url=True),'rU')

    obo_txt=fp.read()

    fp.close()

    obo_dict=obo2csv.parse_obo_txt(obo_txt)



    #### combine predictors

    for datadir in argv:

        for Aspect in ["MF","BP","CC"]:

            ## parse COFACTOR result ##

            COFACTOR_dict=dict()

            if Aspect=="CC":

                COFACTOR_dict=parse_GOfreq_result(os.path.join(

                    datadir,"COFACTOR_GOfreq_CC"),excludeGO)

            GOsearchresult_list=glob(os.path.join(datadir,

                "GOsearchresult_*_%s.dat"%Aspect))

            if GOsearchresult_list and not COFACTOR_dict:

                COFACTOR_dict=parse_COFACTOR_result(

                    GOsearchresult_list[0],obo_dict,excludeGO)



            ## parse GOfreq result ##

            GOfreq_dict=parse_GOfreq_result(os.path.join(

                datadir,"combine_gwGOfreq_"+Aspect),excludeGO)



            psiblast_file=os.path.join(datadir,"psiblast_globalID_"+Aspect)

            blastp_file=os.path.join(datadir,"blastp_globalID_"+Aspect)

            seqID_dict=parse_GOfreq_result(psiblast_file if \

                os.path.isfile(psiblast_file) else blastp_file,excludeGO)

            #seqID_dict=parse_GOfreq_result(os.path.join(

                #datadir,"psiblast_globalID_"+Aspect),excludeGO)



            ## parse PPI2GO result ##

            string_dict=parse_GOfreq_result(os.path.join(

                datadir,"string_swGOfreq_"+Aspect),excludeGO)



            ## combine result ##

            dist_CofactorGOfreq_dict=combine_predictor_dist(

                [COFACTOR_dict,GOfreq_dict])

            dist_GOfreqPPI_dict=combine_predictor_dist(

                [GOfreq_dict,string_dict])

            dist_CofactorPPI_dict=combine_predictor_dist(

                [COFACTOR_dict,string_dict])

            dist_CofactorGOfreqPPI_dict=combine_predictor_dist(

                [COFACTOR_dict,GOfreq_dict,string_dict])

            swdist_CofactorGOfreqPPI_dict= \

                combine_seqID_weighted_predictor_dist({

                'cofactor':COFACTOR_dict,'gofreq':GOfreq_dict,

                'ppi':string_dict,'seqID':seqID_dict})



            ## output result ##

            write_GOpred(dist_CofactorGOfreq_dict,Aspect,

                os.path.join(datadir,"dist_CofactorGOfreq_"+Aspect))

            write_GOpred(dist_GOfreqPPI_dict,Aspect,

                os.path.join(datadir,"dist_GOfreqPPI_"+Aspect))

            write_GOpred(dist_CofactorPPI_dict,Aspect,

                os.path.join(datadir,"dist_CofactorPPI_"+Aspect))

            write_GOpred(dist_CofactorGOfreqPPI_dict,Aspect,

                os.path.join(datadir,"dist_CofactorGOfreqPPI_"+Aspect))

            write_GOpred(swdist_CofactorGOfreqPPI_dict,Aspect,

                os.path.join(datadir,"swdist_CofactorGOfreqPPI_"+Aspect))

