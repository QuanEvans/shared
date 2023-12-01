#!/usr/bin/env python
docstring='''
rescore_COFACTOR.py GOsearchresult_.dat
    re-score COFACTOR GO prediction using GOfreq score. "GOsearchresult.dat"
    is COFACTOR format GO prediction result in the following format:
    template,TMscore,RMSDa,IDEN,Cov,FCscoreGO,PDBGO

    RMSDa is RMSD of aligned region. IDEN is seqID. FCscoreGO is FC-score.
    PDBGO is comma seperated list of GO terms associated with a template. 

options:
-excludeGO=GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575 
    GO term to be excluded
    (default) remove root of 3 aspect and "protein binding"
        GO:0005515 ! protein binding
        GO:0005488 ! binding
        GO:0003674 ! molecular_function
        GO:0008150 ! biological_process
        GO:0005575 ! cellular_component
-datadir=dat
    directory of data files

output:
    COFACTOR_{GOfreq,FCwGOfreq,meanFC,maxFC}_{MF,BP,CC}

MF,BP,CC for molecular function, biological process, cellular component.
GOfreq,FCwGOfreq,maxFC are different scoring function

Template1 GO:0000001 FCscoreGO=0.01
Template2 GO:0000001 FCscoreGO=0.1
Template3 GO:0000002 FCscoreGO=0.02

Confidence Score of GO:0000001
[1] GOfreq (number of templates annotated with GO divided by template number):
    Cscore(GO:0000001)=2/3
[2] FCwGOfreq (sum of FCscoreGO for template annotated with GO divided by 
    template number)
    Cscore(GO:0000001)=(0.01+0.1)/3
[3] meanFC (mean FCscoreGO for template annotated with GO):
    Cscore(GO:0000001)=(0.01+0.1)/2
[4] maxFC (maximum FCscoreGO for template annotated with GO):
    Cscore(GO:0000001)=max{0.01,0.1}=0.1
[5] dist (combining FCscoreGO for all templates annotated with GO):
    Cscore(GO:0000001)=1-(1-0.01)*(1-0.1)=0.109
'''
import sys,os
from module import obo2csv
from module.fetch import wget

obo_url="http://geneontology.org/ontology/go-basic.obo"

def rescore_COFACTOR(GOsearchresult_dict=dict()):
    '''score each GO term by GOfreq,FCwGOfreq,meanFC,maxFC scores. return a
    dict "cscore_dict" whose key is 4 scores and values is a dict, whose key
    is 3 Aspect of GO and value is prediction'''
    cscore_dict={
        'GOfreq':   {'F':dict(),'P':dict(),'C':dict()},
        'FCwGOfreq':{'F':dict(),'P':dict(),'C':dict()},
        'meanFC':   {'F':dict(),'P':dict(),'C':dict()},
        'maxFC':    {'F':dict(),'P':dict(),'C':dict()},
        'dist':    {'F':dict(),'P':dict(),'C':dict()},
        }

    for Aspect in "FPC":
        GO_list=[]
        annotated_template=0. # template annotated with GO of this Aspect
        for template in GOsearchresult_dict:
            if GOsearchresult_dict[template]["PDBGO"][Aspect]:
                annotated_template+=1
                GO_list+=GOsearchresult_dict[template]["PDBGO"][Aspect]

        for Term in set(GO_list):
            cscore_dict["GOfreq"][Aspect][Term]=0.
            cscore_dict["FCwGOfreq"][Aspect][Term]=0.
            cscore_dict["meanFC"][Aspect][Term]=0.
            cscore_dict["maxFC"][Aspect][Term]=0.
            cscore_dict["dist"][Aspect][Term]=1.

            template_annotated_with_GO=0.
            for template in GOsearchresult_dict:
                if Term in GOsearchresult_dict[template]["PDBGO"][Aspect]:
                    template_annotated_with_GO+=1.
                    FCscoreGO=GOsearchresult_dict[template]["FCscoreGO"]

                    cscore_dict["GOfreq"][Aspect][Term]+=1.
                    cscore_dict["FCwGOfreq"][Aspect][Term]+=FCscoreGO
                    cscore_dict["meanFC"][Aspect][Term]+=FCscoreGO
                    cscore_dict["dist"][Aspect][Term]*=(1.-FCscoreGO)
                    if FCscoreGO > cscore_dict["maxFC"][Aspect][Term]:
                        cscore_dict["maxFC"][Aspect][Term]=FCscoreGO
        
            cscore_dict["GOfreq"][Aspect][Term]/=annotated_template
            cscore_dict["FCwGOfreq"][Aspect][Term]/=annotated_template
            cscore_dict["meanFC"][Aspect][Term]/=template_annotated_with_GO
            cscore_dict["dist"][Aspect][Term]=1-cscore_dict["dist"][Aspect][Term]
        
        for scoring in cscore_dict:
            cscore_dict[scoring][Aspect]=''.join(["%s\t%s\t%.2f\n"%(Term,
                Aspect,cscore_dict[scoring][Aspect][Term]) for Term in set(
                GO_list) if float("%.2f"%cscore_dict[scoring][Aspect][Term])])
    return cscore_dict

def parse_COFACTOR_GOsearchresult(GOsearchresult_file="GOsearchresult_.dat",
    obo_dict=dict(),excludeGO=""):
    '''parse COFACTOR format GO prediction result, return a dict, whose
    key is template, and value is a dict with following keys:
    'TMscore','RMSDa','IDEN','Cov','FCscoreGO','PDBGO'
    
    PDBGO is a dict with 3 keys: F,P,C

    obo_dict is obo2csv.obo class for parsing GO hierachy.

    excludeGO is a list of GO to be excluded
    '''
    if isinstance(excludeGO,str):
        excludeGO=excludeGO.split(',')
    excludeGO=set(excludeGO)
    GOsearchresult_dict=dict()
    fp=open(GOsearchresult_file,'rU')
    txt=fp.read()
    fp.close()
    for line in txt.splitlines():
        line=line.split()
        if len(line)<8:
            continue
        template,TMscore,RMSDa,IDEN,Cov,FCscoreGO,comma,PDBGO_list=line
        PDBGO_dict={'F':[],'P':[],'C':[]}
        for Term in set(PDBGO_list.split(','))-excludeGO:
            for Aspect in "FPC":
                if Term in obo_dict[Aspect]["Term"]:
                    PDBGO_dict[Aspect]+=obo_dict.is_a(
                        Term,direct=False,name=False,number=False).split()
        for Aspect in "FPC":
            PDBGO_dict[Aspect]=sorted(set(PDBGO_dict[Aspect])-excludeGO)
        GOsearchresult_dict[template]=dict(
            TMscore=float(TMscore),
            RMSDa=float(RMSDa),
            IDEN=float(IDEN),
            Cov=float(Cov),
            FCscoreGO=float(FCscoreGO),
            PDBGO=PDBGO_dict,
            )
    return GOsearchresult_dict

if __name__=="__main__":
    excludeGO="GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575"
    argv=[]
    datadir = None
    for arg in sys.argv[1:]:
        if arg.startswith("-excludeGO="):
            excludeGO=arg[len("-excludeGO="):].upper()
        elif arg.startswith('-datadir='):
            datadir=arg[len('-datadir='):]
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! Unknown option")
        else:
            argv.append(arg)

    if len(argv)!=1:
        sys.stderr.write(docstring)
        exit()

    # if data directory is specified, change to that directory
    if datadir:
        os.chdir(datadir)
    
    #### parse GO hierachy ####
    fp=open(wget(obo_url,show_url=True),'rU')
    obo_txt=fp.read()
    fp.close()
    obo_dict=obo2csv.parse_obo_txt(obo_txt)

    #### parse COFACTOR GO prediction result ####
    GOsearchresult_dict=parse_COFACTOR_GOsearchresult(argv[0],
        obo_dict,excludeGO)

    #### scoring ####
    cscore_dict=rescore_COFACTOR(GOsearchresult_dict)

    #### write output ####
    #for scoring in cscore_dict:
    for scoring in ["GOfreq"]:
        for Aspect in ["MF","BP","CC"]:
            fp=open("COFACTOR_"+scoring+'_'+Aspect,'w')
            fp.write(cscore_dict[scoring][Aspect[1]])
            fp.close()
