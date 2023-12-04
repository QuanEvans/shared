#!/usr/bin/env python
import os
import argparse
import subprocess
import re
import time
import commands

# fixed variables related to current file path
bindir = os.path.dirname(os.path.abspath(__file__))
cofactor_bindir = os.path.join(bindir, 'COFACTOR')
datdir = os.path.join(cofactor_bindir, 'dat')
runGOfreq = os.path.join(cofactor_bindir, 'GOfreq', 'run_GOfreq.py')
runPPI2GO = os.path.join(cofactor_bindir, 'GOfreq', 'run_PPI2GO.py')
COFACTORmod = os.path.join(cofactor_bindir, 'COFACTORmod')
rescore_COFACTOR = os.path.join(cofactor_bindir, 'rescore_COFACTOR.py')
metaCOFACTOR = os.path.join(cofactor_bindir, 'metaCOFACTOR.py')
qstat = 'qstat -u ' + os.environ["USER"] + ' %j'

def create_parser():
    parser = argparse.ArgumentParser(description='Run COFACTOR')
    parser.add_argument('datadir', type=str, help='data directory')
    parser.add_argument('tag', type=str, help='sequence tag (job id)')
    parser.add_argument('homoflag', type=str, help='homology flag')
    args = parser.parse_args()
    args.datadir = os.path.abspath(args.datadir)
    return args

def getserver():
    server="S10"
    #hostname = subprocess.run("hostname", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) #  use the pipe 
#    print(hostname.stdout.decode())
    #hostname=hostname.stdout.decode().strip('\n')
    (o,hostname)=commands.getstatusoutput("hostname")
    hostname=hostname.strip('\n').strip(' ')
    if hostname.startswith("gl"):
        server="GL"
    elif hostname.startswith("lh"):
        server="S10"
    elif hostname.startswith("amino") or hostname.startswith("zhang"):
        server="amino"
    else:
        server="unknow"
    return server

def replace_template(template, datadir, tag, cofactor_tag,
                     homoflag, libfile, ECfile, GOfile, BSfile1, BSfile2, server=getserver()):
    """
    replace template with variables
    """
    template = re.sub(r'\!S\!', tag, template)
    template = re.sub(r'\!DATADIR\!', datadir, template)
    template = re.sub(r'\!MODEL\!', 'model_1.pdb', template)
    template = re.sub(r'\!TAG\!', cofactor_tag, template)
    template = re.sub(r'\!USER\!', os.environ['USER'], template)
    template = re.sub(r'\!RUN\!', homoflag, template)
    template = re.sub(r'\!LIBFILE\!', libfile, template)
    template = re.sub(r'\!OUTPUT_EC\!', ECfile, template)
    template = re.sub(r'\!OUTPUT_GO\!', GOfile, template)
    template = re.sub(r'\!OUTPUT_BS\!', BSfile1, template)
    template = re.sub(r'\!OUTPUT_POCKET\!', BSfile2, template)
    template = re.sub(r'\!SERVER\!', server, template)
    return template

def runSequenceBasedFunctionPrediction(datadir, homoflag):
    """
    run sequence based function prediction
    """
    configfile = os.path.join(datadir, 'config.py')
    seqfile = os.path.join(datadir, 'seq.fasta')
    ss = os.path.basename(datadir)
    outdir = os.path.dirname(datadir)
    if not os.path.exists(seqfile):
        raise OSError(seqfile + ' does not exist')
    with open(configfile, 'w') as f:
        f.write('ss=[\'' + ss + '\']\n')
        f.write('outdir=\'' + outdir + '\'\n')
        f.write('run=\'' + homoflag + '\'\n')
    subprocess.call(['cp', os.path.join(datdir, 'gene_ontology.obo'), os.path.join(datadir, 'go-basic.obo')])
    print('run GOfreq and PPI2GO')
    subprocess.call([runGOfreq, configfile])
    subprocess.call([runPPI2GO, configfile])

def main(args):
    datadir = args.datadir
    tag = args.tag
    homoflag = args.homoflag

    # write a start file
    with open(os.path.join(datadir, 'start'), 'w') as f:
        f.write('start\n')

    if 'PERL5LIB' in os.environ:
        os.environ['PERL5LIB'] += ':/nfs/amino-home/zhanglabs/amino-modules/2019.11-5.16.3/lib/perl5'
    else:
        os.environ['PERL5LIB'] = '/nfs/amino-home/zhanglabs/amino-modules/2019.11-5.16.3/lib/perl5'

    libfile = 'PDBsearchresult_' + tag + '.dat'
    ECfile = 'ECsearchresult_' + tag + '.dat'
    GOfile = 'GOsearchresult_' + tag + '.dat'
    MFfile = 'GOsearchresult_' + tag + '_MF.dat'
    BPfile = 'GOsearchresult_' + tag + '_BP.dat'
    CCfile = 'GOsearchresult_' + tag + '_CC.dat'
    BSfile1 = 'Bsites_' + tag + '.dat'
    BSfile2 = 'Bpockets_' + tag + '.dat'
    server = getserver()

    runSequenceBasedFunctionPrediction(datadir, homoflag)

    if not (os.path.exists(os.path.join(datadir, MFfile)) and
            os.path.exists(os.path.join(datadir, BPfile)) and
            os.path.exists(os.path.join(datadir, CCfile)) and
            os.path.exists(os.path.join(datadir, libfile)) and
            os.path.exists(os.path.join(datadir, ECfile)) and
            os.path.exists(os.path.join(datadir, 'BSITE_model1', BSfile1))):
        cofactor_tag = 'CFu_' + tag + '_' + homoflag
        jobname = os.path.join(datadir, cofactor_tag)
        with open(COFACTORmod, 'r') as f:
            template = f.read()
        template = replace_template(template, datadir, tag, cofactor_tag,
                                    homoflag, libfile, ECfile, GOfile, BSfile1, BSfile2, server)
        with open(jobname, 'w') as f:
            f.write(template)

        print('run COFACTOR\n')
        subprocess.call(['chmod', 'a+x', jobname])
        subprocess.call([jobname])
        subprocess.call(['rm', '-rf', '/tmp/' + os.environ["USER"] + '/' + tag])

    print('rescore COFACTOR')
    subprocess.call(['cd', datadir])
    subprocess.call([rescore_COFACTOR, GOfile, '-datadir=' + datadir])


    while True:
        if (os.path.exists(os.path.join(datadir, 'combine_gwGOfreq_MF')) and
            os.path.exists(os.path.join(datadir, 'combine_gwGOfreq_BP')) and
            os.path.exists(os.path.join(datadir, 'combine_gwGOfreq_CC')) and
            os.path.exists(os.path.join(datadir, 'string_gwGOfreq_MF')) and
            os.path.exists(os.path.join(datadir, 'string_gwGOfreq_BP')) and
            os.path.exists(os.path.join(datadir, 'string_gwGOfreq_CC'))):
            break
        qzy = subprocess.check_output(qstat, shell=True)
        if re.search(r'GOF_' + tag + '_', qzy) or re.search(r'PIG_' + tag + '_', qzy):
            time.sleep(300)
        else:
            print('ERROR! GOfreq or PPI2GO died. Force output generation\n')
            break

    print('run metaCOFACTOR')
    subprocess.call(['cd', datadir])
    subprocess.call([metaCOFACTOR, datadir])

    # write a finish file
    with open(os.path.join(datadir, 'finish'), 'w') as f:
        f.write('finish\n')

if __name__ == '__main__':
    args = create_parser()
    main(args)
