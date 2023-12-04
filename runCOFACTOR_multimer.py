from Bio import PDB
from Bio.SeqUtils import seq1
import json
import os
import argparse
from string import Template

def split_multimer(job_id:str, structure_pdb:str, seq_fasta:str, chain_id_map:str, chain_list:str, output_dir:str):
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

    # read structure
    structure = parser.get_structure('model1', structure_pdb)
    # split multimer
    # only first model
    model = structure[0]
    for chain in model:
        if chain.id in non_redundant_chain_set:
            chain_output_dir = os.path.join(output_dir, '{}_{}_cofactor'.format(job_id, chain.id))
            if not os.path.exists(chain_output_dir):
                os.makedirs(chain_output_dir)
            io.set_structure(chain)
            io.save(os.path.join(chain_output_dir, 'model_1.pdb'))

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
                with open(os.path.join(chain_output_dir, 'seq.fasta'), 'w') as f:
                    f.write('>{}\n'.format(sequence_name))
                    f.write('{}\n'.format(sequence_from_pdb))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split multimer into single chain and sequence file for a given job')
    parser.add_argument('datadir', help='datadir for DMFold')
    args = parser.parse_args()
    datadir = os.path.abspath(args.datadir)
    job_id = os.path.basename(datadir)
    structure_pdb = os.path.join(datadir, 'model_1.pdb')
    seq_fasta = os.path.join(datadir, 'seq.fasta')
    chain_id_map = os.path.join(datadir, 'AlphaFold2', 'seq', 'msas', 'chain_id_map.json')
    chain_list = os.path.join(datadir, 'chain.list')
    output_dir = datadir
    split_multimer(job_id, structure_pdb, seq_fasta, chain_id_map, chain_list, output_dir)