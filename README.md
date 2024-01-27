# Messenger RNA Design via Expected Partition Function and Continuous Optimization
This repository  hosts the source code for the research paper titled "Messenger RNA Design via Expected Partition Function and Continuous Optimization." 


## Dependencies
- Clang 11.0.0 (or above) or GCC 4.8.5 (or above)
- Python 3
- Biopython

## To Compile
Run the following commands:

```
chmod +x setup.sh
./setup.sh
```

## Usage

The mRNA Design program is designed to process protein sequences provided in a FASTA file and apply our mRNA Design algorithm to each sequence. The results will be stored in `output_dir`, with each sequence having its own subfolder containing the raw outputs from all executions. The optimal solution for each sequence, determined across all runs, will be displayed to the user.

### Command Syntax

Use the following command format to run the mRNA Design program:

```
python mRNA_Design.py [--fasta <path>] [--output_dir <path>] [--beam_size <int>] [--lr <float>] [--epsilon <float>] [--num_iters <int>] [--num_runs <int>] [--num_threads <int>]
```

- `--fasta <path>`: Specifies the path to the input protein FASTA file. The default is `./examples.fasta`.
- `--output_dir <path>`: Sets the directory for saving output files. The default is `./outputs`.
- `--beam_size <int>`: Determines the beam size for beam pruning. The default is `200`. A smaller beam size can speed up the process but may result in search errors.
- `--lr <float>`: Sets the learning rate for projected gradient descent. The default is `0.03`.
- `--epsilon <float>`: Defines the epsilon parameter for soft-MFE initialization. The default is `0.5`.
- `--num_iters <int>`: Specifies the number of optimization steps. The default is `30`.
- `--num_runs <int>`: Indicates the number of execution runs per sample file. The default is `20`. More runs increase the likelihood of finding optimal solutions but require more computational resources.
- `--num_threads <int>`: Sets the number of threads in the thread pool. The default is `8`. Adjust this based on your CPU's core count to optimize parallel processing without overloading your system.

### Example Command

To run the program with a specific set of parameters, you can use a command similar to the following:

```
python mRNA_Design.py --fasta examples.fasta --output_dir ./outputs --beam_size 200 --lr 0.03 --epsilon 0.5 --num_iters 30 --num_runs 20 --num_threads 8
```
