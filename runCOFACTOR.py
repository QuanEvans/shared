import os
import argparse
import subprocess
import re

### fixed variables related to current file path ###
bindir = os.path.dirname(os.path.abspath(__file__))
cofactor_bindir = os.path.join(bindir, 'COFACTOR') # COFACTOR bin directory
datdir = os.path.join(cofactor_bindir, 'dat') # dat file for cofactor
runGOfreq = os.path.join(cofactor_bindir, 'GOfreq', 'runGOfreq.py') # run GOfreq
runPPI2GO = os.path.join(cofactor_bindir, 'GOfreq', 'runPPI2GO.py') # run PPI2GO
COFACTORmod = os.path.join(cofactor_bindir, 'COFACTORmod') # cofactor template file
rescore_COFACTOR = os.path.join(cofactor_bindir, 'rescore_COFACTOR.pl') # rescore cofactor
metaCOFACTOR = os.path.join(cofactor_bindir, 'metaCOFACTOR.pl') # meta cofactor
qstat = f'qstat -u {os.environ["USER"]} %j' #checking submitting jobs


def create_parser()->argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run COFACTOR')
    parser.add_argument('datadir', type=str, help='data directory')
    parser.add_argument('tag', type=str, help='sequence tag (job id)')
    parse.add_argument('homoflag', type=str, help='homology flag')
    args = parser.parse_args()
    args.datadir = os.path.abspath(args.datadir)
    return args

def replace_template(template:str, datadir:str, model:str, tag:str,cofactor_tag:str,
    homoflag:str, libfile:str, ECfile:str, GOfile:str, BSfile1:str, BSfile2:str)->str:
    """
    replace template with variables
    """
    template = re.sub(r'\!S\!', tag, template)
    template = re.sub(r'\!DATADIR\!', datadir, template)
    template = re.sub(r'\!MODEL\!', model, template)
    template = re.sub(r'\!TAG\!', cofactor_tag, template)
    template = re.sub(r'\!USER\!', os.environ['USER'], template)
    template = re.sub(r'\!RUN\!', homoflag, template)
    template = re.sub(r'\!LIBFILE\!', libfile, template)
    template = re.sub(r'\!OUTPUT_EC\!', ECfile, template)
    template = re.sub(r'\!OUTPUT_GO\!', GOfile, template)
    template = re.sub(r'\!OUTPUT_BS\!', BSfile1, template)
    template = re.sub(r'\!OUTPUT_POCKET\!', BSfile2, template)
    return template

def runSequenceBasedFunctionPrediction(datadir:str, homoflag:str):
    """
    run sequence based function prediction
    """
    # make configuration file
    configfile = os.path.join(datadir, 'config.py')
    seqfile = os.path.join(datadir, 'seq.fasta')
    ss = os.path.basename(datadir) # sequence tag
    outdir = os.dirname(datadir)
    if not os.path.exists(seqfile):
        raise FileNotFoundError(f'{seqfile} does not exist')
    with open(configfile, 'w') as f:
        f.write(f'ss=[\'{ss}\']\n')
        f.write(f'outdir=\'{outdir}\'\n')
        f.write(f'run=\'{homoflag}\'\n')
    # copy gene_ontology.obo to go-basic.obo
    subprocess.run(['cp', os.path.join(datdir, 'gene_ontology.obo'), os.path.join(datadir, 'go-basic.obo')])
    # run GOfreq and PPI2GO
    print('run GOfreq and PPI2GO')
    subprocess.run([runGOfreq, configfile])
    subprocess.run([runPPI2GO, configfile])

def main(args:argparse.Namespace):
    datadir = args.datadir
    tag = args.tag
    homoflag = args.homoflag
    COFACTORmod = COFACTORmod
    # set up environment
    os.environ['PERL5LIB'] = os.environ['PERL5LIB'] + ':/nfs/amino-home/zhengwei/lib/amino-modules/2019.11-5.16.3/lib/perl5'
    # variables
    libfile:str = f'PDBsearchresult_{tag}.dat'
    ECfile:str = f'ECsearchresult_{tag}.dat'
    GOfile:str = f'GOsearchresult_{tag}.dat'
    MFfile:str = f'GOsearchresult_{tag}_MF.dat'
    BPfile:str = f'GOsearchresult_{tag}_BP.dat'
    CCfile:str = f'GOsearchresult_{tag}_CC.dat'
    BSfile1:str = f'Bsites_{tag}.dat'
    BSfile2:str = f'Bpockets_{tag}.dat'
    # first run sequence based function prediction
    runSequenceBasedFunctionPrediction(datadir, homoflag)

    # check if files exist
    if not os.path.exists(os.path.join(datadir, MFfile)) or \
        not os.path.exists(os.path.join(datadir, BPfile)) or \
        not os.path.exists(os.path.join(datadir, CCfile)) or \
        not os.path.exists(os.path.join(datadir, libfile)) or \
        not os.path.exists(os.path.join(datadir, ECfile)) or \
        not os.path.exists(os.path.join(datadir, 'BSITE_model1', BSfile1)):
        # copy model_1.pdb to model1.pdb
        model = os.path.join(datadir, 'model1.pdb')
        # create job file
        cofactor_tag = f'CFu_{tag}_{homoflag}'
        jobname = os.path.join(datadir, cofactor_tag)
        # read template file and replace variables
        with open(COFACTORmod, 'r') as f:
            template = f.read()
        template = replace_template(template, datadir, model, tag, cofactor_tag,
            homoflag, libfile, ECfile, GOfile, BSfile1, BSfile2)
        # write job file
        with open(jobname, 'w') as f:
            f.write(template)

        # run job file
        print('run COFACTOR\n')
        subprocess.run(['chmod', 'a+x', jobname])
        subprocess.run([jobname])
        # remove temporary files
        subprocess.run(['rm', '-rf', f'/tmp/{os.environ["USER"]}/{tag}'])

    # rescore COFACTOR
    print('rescore COFACTOR')
    subprocess.run(['cd', datadir])
    subprocess.run([rescore_COFACTOR, GOfile])

    # wait for output files
    while True:
        if os.path.exists(os.path.join(datadir, 'combine_gwGOfreq_MF')) and \
            os.path.exists(os.path.join(datadir, 'combine_gwGOfreq_BP')) and \
            os.path.exists(os.path.join(datadir, 'combine_gwGOfreq_CC')) and \
            os.path.exists(os.path.join(datadir, 'string_gwGOfreq_MF')) and \
            os.path.exists(os.path.join(datadir, 'string_gwGOfreq_BP')) and \
            os.path.exists(os.path.join(datadir, 'string_gwGOfreq_CC')):
            break
        qzy = subprocess.run(qstat, capture_output=True)
        if re.search(r'GOF_'+tag+'_' or re.search(r'PIG_'+tag+'_', qzy):
            time.sleep(300) # wait for 5 minutes
        else:
            print('ERROR! GOfreq or PPI2GO died. Force output generation\n')
            break
    # run metaCOFACTOR
    print('run metaCOFACTOR')
    subprocess.run(['cd', datadir])
    subprocess.run([metaCOFACTOR, '.'])

if __name__ == '__main__':
    args = create_parser()
    main(args)