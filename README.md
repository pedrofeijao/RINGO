# RINGO
Ancestral gene order Reconstruction with INtermediate GenOmes.

## Synopsis
RINGO is a software for ancestral reconstruction of gene orders based on the concept of Intermediate Genomes `[1,2]`
([References](#references)). It also has scripts that simulate new datasets and evaluate the quality of the reconstruction on the simulated datasets, comparing with different software for ancestral reconstruction.

## Installation
RINGO is implemented in Python 2.7 and [Cython 0.24](http://cython.org). No special installation is needed, other than just cloning or downloading all files in a folder in your computer.

The list of required Python packages is included in the `requirements.txt` file. 
If you already have Python installed, you can install all needed packages with:
```
pip install -r requirements.txt
```
### Optional Requirements

* [DeClone](https://github.com/yannponty/DeClone) [3], to find ancestral adjacency weights.


## Running RINGO

In this section each script in RINGO will be described. As a general rule, all files with the `.py` extension are executable 
scripts. Files with `.pyx` extension are libraries compiled by Cython, and not meant to be run directly.

All python scripts can be called with a `-h` option, which shows all the basic usage for the script.

### Quick Start

The main script is `ringo.py`, that runs the ancestral reconstruction algorithm described in [2]. 
Given a phylogenetic tree and a set of extant genomes at the leaves, RINGO reconstructs contiguous ancestral regions (CARs) 
in the internal nodes of the tree. The simplest usage is with:

```
ringo.py -i <GENOME_FILE> -t <TREE_FILE> -o <OUTPUT_FOLDER>
```

The genome file should be in the commonly used GRIMM format (more input formats are available), and the tree in Newick format.
The output folder will be created if not existent already, and the following files are created as output:
* `ringo_genomes_custom_weight.txt`: Ancestral genomes file. Genome names correspond to the internal nodes of the tree.
* `ringo_tree.nwk`: Newick tree. Should be the same as the input tree, possibly with new labels for the internal nodes if they were empty on the input tree.

#### Examples

Ringo has some example datasets where you can test it. Inside of the `src/ringo` folder, run:

```
F='../../examples/sim5000'
./ringo.py -i $F/extant_genomes.txt -t $F/sim_tree.nwk -o $F
```
or 
```
F='../../examples/sim500'
./ringo.py -i $F/extant_genomes.txt -t $F/sim_tree.nwk -o $F
```
to run simulated datasets with 5000 or 500 genes, respectively.

#### Test pipeline

You can run all scripts on RINGO using the `test_all.sh` script, that simulates a new dataset, generates adjacency weights, using also DeClone [3] if available, runs RINGO with custom weights and with DeClone weights and runs SCJ small phylogeny [4]. It also tries to run two other ancestral reconstruction methods: [MGRA](https://github.com/ablab/mgra) [5] and [PhySca](https://github.com/nluhmann/PhySca) [6]. At the end, it outputs a table with the quality reconstruction results. Just run:
`./test_all.sh`

You can also choose the number of extant genomes and genes in the simulation:
`./test_all.sh 6 1000`
to simulate a dataset with 6 genomes and 1000 genes, and run the test pipeline.

### Configuration File

The `ringo.cfg` file has some user configuration options, such as file paths and output names. The default configuration should work without change, and some specific changes are described in more detail in the following sections.

### RINGO Usage

See the [USAGE.md](docs/USAGE.md) file for detailed documentation.

## License

The MIT License (MIT)

Copyright (c) 2016 Pedro Feijão.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## References

* [1] Feijão, P. (2015). Reconstruction of ancestral gene orders using intermediate genomes. *BMC Bioinformatics*, **16**(Suppl 14), S3. http://doi.org/10.1186/1471-2105-16-S14-S3
* [2] Feijão, P. and Araújo, E. (2016). Fast ancestral gene order reconstruction of genomes with unequal gene content. *BMC Bioinformatics*, hopefully to appear.
* [3] Chauve, C., Ponty, Y., Zanetti, J.P.P.: Evolution of genes neighborhood within reconciled phylogenies: an ensemble approach. *BMC Bioinformatics* **16**(19), 1–9 (2015)
* [4] Biller, P., Feijão, P., & Meidanis, J. (2013). Rearrangement-based phylogeny using the single-cut-or-join operation. *IEEE/ACM Transactions on Computational Biology and Bioinformatics*, **10**(1), 122–134. http://doi.org/10.1109/TCBB.2012.168
* [5] Avdeyev, P., Jiang, S., Aganezov, S., Hu, F., & Alekseyev, M. A. (2016). Reconstruction of ancestral genomes in presence of gene gain and loss. *Journal of Computational Biology*, **23**(3), 150–164. http://doi.org/10.1089/cmb.2015.0160
* [6] Luhmann, N., Thévenin, A., Ouangraoua, A., Wittler, R., & Chauve, C. (2016). The SCJ Small Parsimony Problem for Weighted Gene Adjacencies. In A. Bourgeois, P. Skums, X. Wan, & A. Zelikovsky (Eds.), *Bioinformatics Research and Applications: 12th International Symposium*, ISBRA 2016, Minsk, Belarus, June 5-8, 2016, Proceedings (pp. 200–210). Cham: Springer International Publishing. http://doi.org/10.1007/978-3-319-38782-6_17
