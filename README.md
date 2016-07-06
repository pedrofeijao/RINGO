# RINGO
Ancestral gene order Reconstruction with INtermediate GenOmes.

## Synopsis
RINGO is a software for ancestral reconstruction of gene orders based on 
[References](#references)

## Installation
RINGO is implemented in Python 2.7 and [Cython 0.24](http://cython.org). 

The list of required packages is included in the `requirements.txt` file. 
If you already have Python installed, you can install all needed packages with:
```
pip install -r requirements.txt
```
### Optional Requirements

* [DeClone](https://github.com/yannponty/DeClone), to find ancestral adjacency weights.

## Running RINGO

In this section each script in RINGO will be described. As a general rule, all files with the `.py` extension are executable 
scripts. Files with `.pyx` extension are libraries compiled by Cython, and not meant to be run directly.

All python scripts can be called with a `-h` option, which shows all the basic usage for the script.

### Quick Start

The main script is `ringo.py`, that runs the ancestral reconstruction algorithm described in (CITATION?). 
Given a phylogenetic tree and a set of genomes at its leaves, RINGO reconstructs contiguous ancestral regions (CARs) 
in the internal nodes of the tree. The simplest usage is with:

```
./ringo.py -i <GENOME_FILE> -t <TREE_FILE> -o <OUTPUT_FOLDER>
```

The genome file should be in the commonly used GRIMM format (more input formats are available), and the tree in Newick format.
The output folder will be created if not existent already, and the following files are created as output:
* `ringo_genomes_custom_weight.txt`: Ancestral genomes file. Genome names correspond to the internal nodes of the tree.
* `ringo_tree.nwk`: Newick tree. Should be the same as the input tree, with new labels for the internal nodes if they were empty in the input tree.


## References

* [1] Feijão, P. (2015). Reconstruction of ancestral gene orders using intermediate genomes. *BMC Bioinformatics*, 16(Suppl 14), S3. http://doi.org/10.1186/1471-2105-16-S14-S3

* [2] Feijão, P. and Araújo, E. (2016). Fast ancestral gene order reconstruction of genomes with unequal gene content. *BMC Bioinformatics*, hopefully to appear.



