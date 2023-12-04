#!/usr/bin/env python
from Bio import PDB
from Bio.SeqUtils import seq1
import json, time, sys
import os, commands
import argparse
from string import Template
import subprocess

COFACTOR_template=Template("""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 10:00:00
#SBATCH --mem=7gb
#SBATCH --job-name="$TAG"
#SBATCH --output="$JOBNAME.out"
#SBATCH --error="$JOBNAME.err"
#SBATCH --partition=$PARTITION
#SBATCH --account=$ACCOUNT
###SBATCH --qos urgent
$CMD
""")

file_dir = os.path.dirname(os.path.abspath(__file__))
run_cofactor_path = os.path.join(file_dir, 'runCOFACTOR.py')

def getserver():
    server="S10"
    #hostname = subprocess.run("hostname", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) #  use the pipe 
    #print(hostname.stdout.decode())
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

def submit_job(jobname,cmd,server):
    '''write command "cmd" to PBS script "jobname", submit the script'''

    account='zhanglab'
    partition='batch'
    reservation=''
    Qos=' -q '+"casp"
    if server == "S10":
        account='sigbio_project1'
        partition='sigbio'
        Qos=''

    if server == "GL":
        account='petefred1'
        partition='standard'
        Qos=''


    fp=open(jobname,'w')
    fp.write(COFACTOR_template.substitute(dict(
        TAG=os.path.basename(jobname),
        JOBNAME=os.path.abspath(jobname),
        CMD=cmd,
        ACCOUNT=account,
        PARTITION=partition,
    )))
    fp.close()
    #os.chmod(jobname, os.stat(jobname).st_mode|0111)
    os.system("chmod 777 %s"%jobname)

    while True:
        p=subprocess.Popen("sbatch "+' '+Qos+' '+jobname,shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if stdout.strip():
            sys.stdout.write(jobname+" submitted\n")
            break
        else:
            time.sleep(5) # something wrong with qsub
    return stdout.strip() # return PBS job ID

def showq():
    '''return job queue status'''
    cmd="squeue -o %j"
    p=subprocess.Popen(cmd,shell=True,
        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=p.communicate()
    return stdout

def split_multimer(job_id, structure_pdb, seq_fasta, chain_id_map, chain_list, output_dir):
    """
    Split multimer into single chain and sequence file for a given job
    each chain and sequence file will be saved in a separate folder for running COFACTOR

    Args:
        job_id (str): job id of DMFold
        structure_pdb (str): final model structure in pdb format
        seq_fasta (str): input sequence file in fasta format
        chain_id_map (str): path to chain id map file in json format, under datadir/AlphaFold2/seq/msas/chain_id_map.json
        chain_list (str): path to non-redundant chain list file, under datadir/chain.list
        output_dir (str): output directory, also is datadir for DMFold

    Raises:
        ValueError: if sequence mismatch between pdb and fasta or chain id mismatch between pdb and fasta

    Returns:
        list of cofacor dir
    """


    parser = PDB.PDBParser()
    io = PDB.PDBIO()
    # read chain_list
    non_redundant_chain_set = set()
    with open(chain_list, 'r') as f:
        for line in f:
            chain_id = line.strip().split('-')[-1]
            non_redundant_chain_set.add(chain_id)
    # read chain_id_map
    with open(chain_id_map, 'r') as f:
        chain_id_map = json.load(f)
    # record the chain id
    sequence_id_map = {}
    for key, value in chain_id_map.items():
        sequence = value['sequence']
        if sequence not in sequence_id_map:
            sequence_id_map[sequence] = {key}
        else:
            sequence_id_map[sequence].add(key)


    cofactor_dirs = []
    # read structure
    structure = parser.get_structure('model1', structure_pdb)
    # split multimer
    # only first model
    model = structure[0]
    for chain in model:
        if chain.id in non_redundant_chain_set:

            sequence_from_pdb = ''
            for res in chain:
                if PDB.is_aa(res, standard=True):
                    sequence_from_pdb += seq1(res.get_resname())
            sequence_from_fasta = chain_id_map[chain.id]['sequence']
            sequence_name = chain_id_map[chain.id]['description']

            if sequence_from_pdb != sequence_from_fasta:
                print('sequence mismatch between pdb and fasta')
                print('sequence from pdb: {}'.format(sequence_from_pdb))
                print('sequence from fasta: {}'.format(sequence_from_fasta))
                print('sequence name: {}'.format(sequence_name))
                print('chain id: {}'.format(chain.id))
                print('job id: {}'.format(job_id))
                raise ValueError('sequence mismatch between pdb and fasta')
            elif chain.id not in sequence_id_map[sequence_from_pdb]:
                print('chain id mismatch between pdb and fasta')
                print('sequence: {}'.format(sequence_from_pdb))
                print('chain id from pdb: {}'.format(chain.id))
                print('chain id from fasta: {}'.format(sequence_id_map[sequence_from_pdb]))
                print('sequence name: {}'.format(sequence_name))
                print('job id: {}'.format(job_id))
                raise ValueError('chain id mismatch between pdb and fasta')
            else:
                chain_output_dir = os.path.join(output_dir, '{}_{}_cofactor'.format(job_id, chain.id))
                if not os.path.exists(chain_output_dir):
                    os.makedirs(chain_output_dir)
                cofactor_dirs.append(chain_output_dir)
                # save pdb and fasta
                io.set_structure(chain)
                io.save(os.path.join(chain_output_dir, 'model_1.pdb'))

                with open(os.path.join(chain_output_dir, 'seq.fasta'), 'w') as f:
                    f.write('>{}\n'.format(sequence_name))
                    f.write('{}\n'.format(sequence_from_pdb))
    return cofactor_dirs

def check_compeleted(job_id, datadir, running_jobs):
    """
    Check if COFACTOR jobs are completed

    Args:
        job_id (str): job id of DMFold
        datadir (str): output directory, also is datadir for DMFold
        running_jobs (set): set of running jobs

    Returns:
        bool: True if all COFACTOR jobs are completed, False otherwise
    """
    start_file = os.path.join(datadir, 'start')
    finish_file = os.path.join(datadir, 'finish')
    if not os.path.exists(start_file):
        return False
    if job_id in running_jobs:
        return False
    if os.path.exists(finish_file):
        return True
    if os.path.exists(start_file) and not os.path.exists(finish_file) and job_id not in running_jobs:
        print("job: {} died".format(job_id))
        return True

def run_cofactor_multimer(datadir, homoflag):
    """
    Run COFACTOR for a given job

    Args:
        datadir (str): datadir for DMFold
    """
    job_id = os.path.basename(datadir)
    structure_pdb = os.path.join(datadir, 'model_1.pdb')
    seq_fasta = os.path.join(datadir, 'seq.fasta')
    chain_id_map = os.path.join(datadir, 'AlphaFold2', 'seq', 'msas', 'chain_id_map.json')
    chain_list = os.path.join(datadir, 'chain.list')
    output_dir = datadir
    cofactor_dirs = split_multimer(job_id, structure_pdb, seq_fasta, chain_id_map, chain_list, output_dir)
    # submit jobs
    tags_datadir_pairs = []
    for chain in cofactor_dirs:
        datadir = chain
        jobname = os.path.join(datadir, 'cofactor.sh')
        tag = os.path.basename(datadir)
        tags_datadir_pairs.append((tag, datadir))
        cmd = 'python {} {} {}'.format(run_cofactor_path, chain, tag, homoflag)
        submit_job(jobname, cmd, getserver())
    
    # wait for jobs to finish
    while True:
        completed = True
        stdout = showq()
        running_jobs = stdout.split('\n')
        running_jobs = set(running_jobs)

        for tag, datadir in tags_datadir_pairs:
            if not check_compeleted(tag, datadir, running_jobs):
                completed = False
                break

        if completed:
            break
        else:
            time.sleep(300)
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split multimer into single chain and sequence file for a given job')
    parser.add_argument('datadir', help='datadir for DMFold')
    parser.add_argument('homoflag', help='homology flag')
    args = parser.parse_args()
    datadir = os.path.abspath(args.datadir)
    homoflag = args.homoflag
    run_cofactor_multimer(datadir, homoflag)
