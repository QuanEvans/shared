<<<<<<< HEAD
#!/nfs/amino-home/liyangum/anaconda3/bin/python3

=======
>>>>>>> parent of 623b2ee (upadte)
import os, subprocess

import argparse

from collections import defaultdict

import multiprocessing as mp

import shutil

import re

import pickle


<<<<<<< HEAD

=======
>>>>>>> parent of 623b2ee (upadte)
DMsearch_docstring= """

    DMsearch: A tool for protein complex analysis using Foldseek and MMalign.



    This program offers functionalities for searching protein complexes, creating databases, 

    and aligning protein structures. It utilizes Foldseek for searching and creating databases 

    of protein sequences and structures, and MMalign for aligning protein complexes. 



    Subcommands:

    - search: Searches protein complexes against a database and aligns the hits.

    - createdb: Creates a database from provided protein sequences or structures.

    - align: Aligns two given protein complexes.



    Each subcommand has its own set of arguments and options. Use 'DMsearch.py <subcommand> --help' 

    for more information on a specific subcommand.

    """

search_docstring = """

    Search Subcommand:

    Searches protein complexes against a specified database using Foldseek and performs alignments on the hits using USalign/MMalign.



    Arguments:

    input_dir - Directory containing the input protein complex files.

    target_dir - Directory containing the target protein complex files for alignment.

    sequenceDB - Path to the Foldseek database to be used for searching.

    output_dir - Directory where the search results and alignments will be stored.

    --min_tmscore - Minimum TM-score to consider a hit (default: 0.2).

    --threads - Number of threads to use for searching (default: Number of CPU cores).

    -v, --verbose - Enable verbose mode to output more information during the process.

    --min_hits - Minimum number of hits to consider second iteration of search (default: 100).

    -c, --continue_search - Continue search from cached common hits.

    -exclude, --exclude_list - List of PDB IDs to exclude from search.

    """

createdb_docstring = """

    Createdb Subcommand:

    Creates a database for Foldseek using provided protein sequences or structures.



    Arguments:

    input_dir - Directory containing the protein sequences or structures.

    sequenceDB - Path where the created database will be stored.

    --threads - Number of threads to use for database creation (default: Number of CPU cores).

    -v, --verbose - Enable verbose mode to output more information during the process.

    """

align_docstring = """

    Align Subcommand:

    Aligns two specified protein complexes using MMalign, providing a detailed structural comparison.



    Arguments:

    query - Path to the query protein complex file.

    target - Path to the target protein complex file.

    output - Path for the output file where alignment results will be stored.

    """



def get_absolute_path(path:str)->str:

    """

    Returns the absolute path of a file or directory.



    Args:

    path (str): Path to the file or directory.



    Returns:

    str: Absolute path of the file or directory.

    """

    if path is None:

        return None

    if shutil.which(path) is not None:

        return path

    if path.startswith('.'):

        return os.path.abspath(path)

    else:

        return os.path.abspath(os.path.join(os.getcwd(), path))

    

def create_parser()->argparse.ArgumentParser:

    file_dir = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(prog='DMsearch.py', description=DMsearch_docstring, formatter_class=argparse.RawDescriptionHelpFormatter)

    subparsers = parser.add_subparsers(dest='command', help='sub-command help')



    # create parser for search command

    parser_search = subparsers.add_parser('search', help='Search protein complexes in a database and perform alignments.', description=search_docstring,formatter_class=argparse.RawDescriptionHelpFormatter)

    parser_search.add_argument('input_path', help='input directory or file path')

    parser_search.add_argument('target_dir', help='target directory')

    parser_search.add_argument('sequenceDB', help='databases to search')

    parser_search.add_argument('output_dir', help='output directory')

    parser_search.add_argument('--min_tmscore', help='minimum TM-score to consider a hit', default=0.3, type=float)
<<<<<<< HEAD

    parser_search.add_argument('--threads', help='number of threads', default=1, type=int)

=======
    parser_search.add_argument('--threads', help='number of threads', default=mp.cpu_count(), type=int)
>>>>>>> parent of 623b2ee (upadte)
    parser_search.add_argument('-v', '--verbose', help='verbose', action='store_true')

    parser_search.add_argument('--min_hits', help='minimum number of hits to consider second iteration of search', default=100, type=int)

    parser_search.add_argument('-c', '--continue_search', help='continuous mode', action='store_true')

    parser_search.add_argument('-exclude', '--exclude_list', help='list of PDB IDs to exclude from search', default=None, type=str)
<<<<<<< HEAD

    parser_search.add_argument('--mmalign', help='mmalign executable path', default=os.path.join(file_dir,'DeepMSAFold','MMalign'))

=======
    parser_search.add_argument('--mmalign', help='mmalign executable path', default='MMalign')
>>>>>>> parent of 623b2ee (upadte)
    parser_search.add_argument('--foldseek', help='foldseek executable path', default=os.path.join(file_dir,'DeepMSAFold','foldseek'))



    # create db command

    parser_createdb = subparsers.add_parser('createdb', help='Create a Foldseek database from protein sequences or structures.', description=createdb_docstring, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser_createdb.add_argument('input_dir', help='input directory')

    parser_createdb.add_argument('sequenceDB', help='output sequenceDB')
<<<<<<< HEAD

    parser_createdb.add_argument('--threads', help='number of threads', default=1, type=int)

=======
    parser_createdb.add_argument('--threads', help='number of threads', default=mp.cpu_count(), type=int)
>>>>>>> parent of 623b2ee (upadte)
    parser_createdb.add_argument('-v', '--verbose', help='verbose', action='store_true')

    parser_createdb.add_argument('--foldseek', help='foldseek executable path', default='foldseek')



    # create parser for align command

    parser_align = subparsers.add_parser('align', help='Align two protein complexes using MMalign.', description=align_docstring, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser_align.add_argument('query', help='query file')

    parser_align.add_argument('target', help='target file')

    parser_align.add_argument('output', help='output file')

    parser_align.add_argument('--mmalign', help='mmalign executable path', default=os.path.join(file_dir,'DeepMSAFold','MMalign'))

    

    args = parser.parse_args()



    if not hasattr(args, 'command'):

        parser.print_help()

        exit(1)



    # Processing arguments based on the subcommand

    if args.command == 'search':

        args.input_path = os.path.abspath(args.input_path)

        args.target_dir = os.path.abspath(args.target_dir)

        args.sequenceDB = os.path.abspath(args.sequenceDB)

        args.output_dir = os.path.abspath(args.output_dir)

    

    elif args.command == 'createdb':

        args.input_dir = os.path.abspath(args.input_dir)

        args.sequenceDB = os.path.abspath(args.sequenceDB)

    

    elif args.command == 'align':

        args.query = os.path.abspath(args.query)

        args.target = os.path.abspath(args.target)

        args.output = os.path.abspath(args.output)



    return args



def run_foldseek(input_path: str, databases: str, output_dir: str, threads=mp.cpu_count(), verbose:bool=False, num_iterations:int=1, exhaustive:bool=False, foldseek_path:str='foldseek') -> str:

    """

    Runs the Foldseek easy-search command.



    Args:

    input_dir (str): Directory containing the input files.

    databases (str): Path to the Foldseek databases.

    output_dir (str): Directory where the output will be stored.

    threads (int, optional): Number of threads to use. Defaults to the number of available CPU cores.

    verbose (bool, optional): Enable verbose mode to output more information during the process. Defaults to False.

    num_iterations (int, optional): Number of iterations to run. Defaults to 1.

    exhaustive (bool, optional): Enable exhaustive search. Defaults to False.



    Returns:

    str: Path to the output TSV file.

    """



    # Create output and temporary directories if they do not exist

    tmp_dir = os.path.join(output_dir, 'tmp')

    os.makedirs(output_dir, exist_ok=True)

    os.makedirs(tmp_dir, exist_ok=True)



    # Define the output TSV file path

    out_tsv = os.path.join(output_dir, f'result_{num_iterations}.tsv')



    # Return the output path if the file already exists

    if os.path.exists(out_tsv):

        return out_tsv



    # Define the Foldseek search command

    search_cmd = [

        foldseek_path, 'easy-search', input_path, databases, out_tsv, tmp_dir,

        '--format-output', "query,target,evalue,bits,alntmscore,qtmscore,ttmscore,prob",

        '--threads', str(threads),  # Including threads argument in the command

        '--num-iterations', str(num_iterations), 

    ]

    if exhaustive:

        search_cmd.append('--exhaustive-search')



    if not verbose:

        search_cmd.append('-v')

        search_cmd.append('0')



    # Run the Foldseek command

    try:

        subprocess.run(search_cmd, check=True)

        print('Foldseek search finished')

    except subprocess.CalledProcessError as e:

        print('Foldseek search failed')

        print(e)

        exit(1)

    

    return out_tsv



def find_hits(result_tsv: str, min_tmscore: float = 0.3):

    """

    Processes a result TSV file to find hits above a specified TM-score without using Pandas.



    Args:

    result_tsv (str): Path to the TSV file containing Foldseek results.

    min_tmscore (float, optional): Minimum TM-score to consider as a hit. Defaults to 0.2.



    Returns:

    dict: Dictionary mapping queries to their common hits.

    """



    # Initialize dictionaries for storing complex-chain mappings and hits

    query_complex_chain_dict = defaultdict(set)

    hits_dict = defaultdict(set)



    with open(result_tsv, 'r') as file:

        for line in file:

            # Split each line into its columns

            columns = line.strip().split('\t')



            # Ensure the line has the correct number of columns

            if len(columns) >= 6:

                query, target, _, _, _, qtmscore, *_ = columns



                # Convert qtmscore to float and compare with min_tmscore

                if float(qtmscore) > min_tmscore:

                    query_complex, *_ = query.split('_')

                    target_complex, *_ = target.split('_')



                    # Skip if query and target complexes are the same

                    if query_complex == target_complex:

                        continue



                    query_complex_chain_dict[query_complex].add(query)

                    hits_dict[query].add(target_complex)



    # Compute common hits for each query complex

    common_hits = {

        query: set.union(*(hits_dict[chain] for chain in chains))

        for query, chains in query_complex_chain_dict.items()

    }



    return common_hits



def run_usalign(query:str, target:str, out:str, mmalign_path:str="MMalign")->None:

    """

    Runs USalign/MMalign command with given query and target, outputting the results to a specified file.



    Args:

    query (str): Path to the query file.

    target (str): Path to the target file.

    out (str): Path for the output file.



    If the output or error file already exists, the function will return early.

    In case of an error during execution, the output file will be renamed to an error file.

    """



    cmd = f'{mmalign_path} {query} {target} -outfmt -1 > {out}'

    err = out.replace('.out', '.err')



    # Check if the error file already exists

    if os.path.exists(err):

        return



    # Check if the output file already exists

    if os.path.exists(out):

        return



    try:

        # Execute the command

        subprocess.run(cmd, shell=True, check=True)

        # print('USalign finished')

    except subprocess.CalledProcessError as e:

        # print('USalign failed')

        # print(e)

        

        # Rename the output file to error file if execution fails

        if os.path.exists(out):

            os.rename(out, err)



def collect_results(align_out:str, output_tsv:str)->None:

    """

    Collects the results from the alignment output directory and writes them to a TSV file.



    Args:

    align_out (str): Path to the alignment output directory.

    output_tsv (str): Path to the output TSV file.

    """



    # Get the list of output files

    out_files = os.listdir(align_out)

    out_files = [f for f in out_files if f.endswith('.out') or f.endswith('.err')]

    out_files = sorted(out_files)



    # Initialize a list for storing the results

    results = []



    # Iterate over the output files

    for out_file in out_files:

        # Get the query and target names

        query, target = out_file.split('.')[0].split('_')

        # if has error file, use nan

        if out_file.endswith('.err'):

            results.append((query, target, float('nan'), float('nan')))

            continue

        # Parse the TM-scores from the output file

        qtm, ttm = parse_tmscores(os.path.join(align_out, out_file))

        # Append the results to the list

        results.append((query, target, qtm, ttm))

    

    # sort by query, qtm, ttm

    results = sorted(results, key=lambda x: (x[0], x[2], x[3]), reverse=True)



    # Write the results to a TSV file

    with open(output_tsv, 'w') as f:

        f.write('Query\tTemplate\tQtmscore\tTtmscore\n')

        for result in results:

            f.write('\t'.join(map(str, result)) + '\n')



def parse_tmscores(result_file:str)->tuple:

    """

    Parses the TM-scores from the output file of USalign/MMalign.



    Args:

    result_file (str): Path to the output file.



    Returns:

    tuple: Tuple containing the TM-scores for the query and target structures.

    """

    # Read the output file

    with open(result_file, 'r') as f:

        result_text = f.read()



    # Regex pattern to find TM-scores

    tm_score_pattern = r"TM-score= ([0-9.]+) \(normalized by length of Structure_(\d):"



    # Find all matches

    matches = re.findall(tm_score_pattern, result_text)



    # Initialize scores

    qtm, ttm = None, None



    for score, structure in matches:

        if structure == '1':

            qtm = float(score)

        elif structure == '2':

            ttm = float(score)



    return qtm, ttm





class Searcher:



    def __init__(self,

        input_dir:str=None,

        target_dir:str=None,

        sequenceDB:str=None,

        output_dir:str=None,

        min_tmscore:float=0.2,

        threads:int=mp.cpu_count(),

        verbose:bool=False,

        min_hits:int=10,

        continue_search:bool=False,

        exclude_list:str=None,

        foldseek_path:str='foldseek',

        mmalign_path:str='MMalign'

    ):

        self.target_dir = target_dir

        self.sequenceDB = sequenceDB

        self.min_tmscore = min_tmscore

        self.threads = threads

        self.verbose = verbose

        self.min_hits = min_hits

        self.continue_search = continue_search

        self.exclude_list = exclude_list

        self.prepossessing(working_dir=output_dir, input_dir=input_dir)

        self.pkl_cache = os.path.join(output_dir, 'cache.pkl')

        self.foldseek_path = foldseek_path

        self.mmalign_path = mmalign_path



        self.num_iterations = 1

        self.query_hits_dict = {}



    def prepossessing(self, working_dir:str, input_dir:str, exclude_list:str=None):

        """

        Preprocesses the input files by creating necessary directories, copying input files to the query directory,

        and storing the list of query files.



        Args:

            working_dir (str): The working directory where the preprocessing will be performed.

            input_dir (str): The directory or file containing the input files.



        Returns:

            None

        """



        self.working_dir = working_dir

        self.tmp_dir = os.path.join(working_dir, 'tmp')

        if os.path.exists(self.working_dir) and not self.continue_search:

            shutil.rmtree(self.working_dir)

        os.makedirs(self.working_dir, exist_ok=True)

        os.makedirs(self.tmp_dir, exist_ok=True)



        query_dir = os.path.join(self.working_dir, 'query_iter1')

        os.makedirs(query_dir, exist_ok=True)



        # exclude list

        if exclude_list is not None and os.path.exists(exclude_list):

            exlude_set = set()

            with open(exclude_list, 'r') as f:

                for line in f:

                    exlude_set.add(line.strip())

            self.exclude_set = exlude_set

        else:

            self.exclude_set = set()



        # copy input files to query dir

        if os.path.isdir(input_dir):

            input_files = os.listdir(input_dir)

            input_files = [f for f in input_files if f.endswith('.pdb') or f.endswith('.cif')]

            for f in input_files:

                if f in self.exclude_set:

                    continue

                shutil.copy(os.path.join(input_dir, f), query_dir)

        else:

            shutil.copy(input_dir, query_dir)



        querys = os.listdir(query_dir)

        self.querys = querys



    def search_iter(self):

        """

        Perform a searching query against a database for one iteration

        """



        search_params = self.get_search_params(self.num_iterations)

        print(f'Starting iteration {self.num_iterations}')

        out_tsv = run_foldseek(search_params['query_dir'], self.sequenceDB, self.working_dir, self.threads, \

            num_iterations=search_params['num_iterations'], exhaustive=search_params['exhaustive_search'], folseek_path=self.folseek_path)

        cur_iter_dict = find_hits(out_tsv, self.min_tmscore)

        total_hits = sum(len(v) for v in cur_iter_dict.values())

        print(f'Found {total_hits} hits in iteration {self.num_iterations}')

        self.update_query_hits_dict(cur_iter_dict)

        self.save_query_hits_dict()

        self.num_iterations += 1



    def search(self):

        """

        Perform up to 3 iterations of searching and alignments



        Returns:

            None

        """

        if self.continue_search:

            self.load_query_hits_dict()

        while not self.is_search_finished():

            self.search_iter()

            self.perpare_next_iter()

        self.run_align()



    def get_search_params(self, iteration:int)->tuple:

        """

        Get the search parameters for a given iteration.



        Args:

            iteration (int): The iteration number.



        Returns:

            tuple: A tuple containing the search parameters.

        """

    

        if iteration == 1:

            return {

                'num_iterations': 1,

                'exhaustive_search': False,

                'query_dir': os.path.join(self.working_dir, 'query_iter1'),

            }



        elif iteration == 2:

            return {

                'num_iterations': 2,

                'exhaustive_search': False,

                'query_dir': os.path.join(self.working_dir, 'query_iter2'),

            }

        elif iteration == 3:

            return {

                'num_iterations': 2,

                'exhaustive_search': True,

                'query_dir': os.path.join(self.working_dir, 'query_iter2'),

            }

        else:

            raise ValueError(f'Invalid iteration number: {iteration}')



    def is_search_finished(self):

        if len(self.querys) == len({k for k, v in self.query_hits_dict.items() if len(v) > self.min_hits}):

            return True

        elif self.num_iterations >= 3:

            return True

        else:

            return False



    def update_query_hits_dict(self, cur_iter_dict:dict):

        for k, v in cur_iter_dict.items():

            if k in self.exclude_set or v in self.exclude_set:

                continue

            if k in self.query_hits_dict:

                self.query_hits_dict[k].update(v)

            else:

                self.query_hits_dict[k] = v



    def save_query_hits_dict(self):

        save_dict = {

            'iteration': self.num_iterations,

            'query_hits_dict': self.query_hits_dict

        }

        pickle.dump(save_dict, open(self.pkl_cache, 'wb'))

    

    def load_query_hits_dict(self):

        if os.path.exists(self.pkl_cache):

            save_dict = pickle.load(open(self.pkl_cache, 'rb'))

            self.num_iterations = save_dict['iteration']

            self.query_hits_dict = save_dict['query_hits_dict']

            print(f'Loaded cached query hits dict from iteration {self.num_iterations}')

            total_hits = sum(len(v) for v in self.query_hits_dict.values())

            print(f'Found {total_hits} hits in iteration {self.num_iterations}')

        else:

            print('No cached query hits dict found')



    def perpare_next_iter(self):

        # Prepare for next iteration

        query_dir = os.path.join(self.working_dir, f'query_iter{self.num_iterations}')

        os.makedirs(query_dir, exist_ok=True)

        for query in self.querys:

            if query not in self.query_hits_dict or len(self.query_hits_dict[query]) < self.min_hits:

                shutil.copy(os.path.join(self.working_dir, 'query_iter1', query), query_dir)

                

    def run_align(self):

        """

        Run alignments for all hits in the query_hits_dict

        """

        # Prepare for alignment

        align_out = os.path.join(self.working_dir, 'align_out')

        os.makedirs(align_out, exist_ok=True)



        # Prepare arguments for alignment

        args_list = []

        for query, targets in self.query_hits_dict.items():

            for target in targets:

                # Get the paths for query and target

                query_path = os.path.join(self.working_dir, 'query_iter1', query)

                target_path = os.path.join(self.target_dir, target)

                query_name = query.split('.')[0]

                target_name = target.split('.')[0]

                if query_name == target_name:

                    continue

                out_path = os.path.join(align_out, f'{query_name}_{target_name}.out')

                err_path = os.path.join(align_out, f'{query_name}_{target_name}.err')



                if os.path.exists(out_path) and os.path.getsize(out_path) == 0:

                    os.remove(out_path)

                if os.path.exists(err_path) and os.path.getsize(out_path) == 0:

                    os.remove(err_path)

                

                if not os.path.exists(out_path) and not os.path.exists(err_path):

                    # Append the arguments to the list

                    args_list.append((query_path, target_path, out_path, self.mmalign_path))

        

        # Run alignments in parallel

        if self.threads == 1:

            for args in args_list:

                run_usalign(*args)

        

        else:

            # Run alignment in parallel

            with mp.Pool(processes=self.threads) as pool:

                pool.starmap(run_usalign, args_list)

        

        print('Finished alignments')

        # Collect results

        output_tsv = os.path.join(self.working_dir, 'search_result.tsv')

        collect_results(align_out, output_tsv)

        print(f'Finished collecting results to {output_tsv}')        



def main(args:argparse.Namespace):

    """

    Main function to handle different subcommands: search, createdb, and align.



    Depending on the subcommand specified in args.command, this function will:

    - Run a foldseek search, find hits, and perform alignments using USalign/MMalign for the 'search' command.

    - Create a database using foldseek for the 'createdb' command.

    - Align two protein complexes using MMalign for the 'align' command.



    Args:

    args: ArgumentParser namespace with command and parameters.

    """



    if args.command == 'search':

        # Get absolute paths

        query_dir = get_absolute_path(args.input_path)

        target_dir = get_absolute_path(args.target_dir)

        sequenceDB = get_absolute_path(args.sequenceDB)

        output_dir = get_absolute_path(args.output_dir)

        exlude_list = get_absolute_path(args.exlude_list)

        foldseek = get_absolute_path(args.foldseek)

        mmalign = get_absolute_path(args.mmalign)

        num_threads = args.threads

        min_tmscore = args.min_tmscore

        verbose = args.verbose

        min_hits = args.min_hits

        continue_search = args.continue_search



        # Run foldseek search

        searcher = Searcher(query_dir, target_dir, sequenceDB, output_dir, min_tmscore, num_threads, \

            verbose, min_hits, continue_search, exlude_list, foldseek, mmalign)

        searcher.search()



    elif args.command == 'createdb':

        args.input_dir = get_absolute_path(args.input_dir)

        args.sequenceDB = get_absolute_path(args.sequenceDB)

        args.foldseek = get_absolute_path(args.foldseek)

        # Run foldseek createdb

        createdb_cmd = [args.foldseek, 'createdb', args.input_dir, args.sequenceDB, '--threads', str(args.threads)]

        if args.verbose:

            createdb_cmd.append('-v')

        subprocess.run(createdb_cmd, check=True)

    

    elif args.command == 'align':

        args.query = get_absolute_path(args.query)

        args.target = get_absolute_path(args.target)

        args.output = get_absolute_path(args.output)

        args.mmalign = get_absolute_path(args.mmalign)

        # Run MMalign

        run_usalign(args.query, args.target, args.output, args.mmalign)



if __name__ == '__main__':

    args = create_parser()
<<<<<<< HEAD

    main(args)
=======
    main(args)




>>>>>>> parent of 623b2ee (upadte)
