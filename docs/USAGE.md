# RINGO Usage

Included in RINGO there are many different scripts: 
* to calculate [adjacency weights](#adjacency-weights)
* to perform the [ancestral reconstruction](#ancestral-reconstruction-with-ringo)
* to generate [simulated datasets](#generating-simulated-datasets)
* to [evaluate quality](#evaluating-reconstruction-quality) of the reconstruction on simulated datasets, also comparing with the results of different available software (TO DO).

## Adjacency Weights

An important step in the algorithm for ancestral reconstruction used by RINGO is the calculation of *adjacency weights* for each ancestral genome in the tree. The weights are real numbers from 0 and 1 that reflect how likely it is for each adjacency (a connection between adjacent genes) to be present in a particular ancestor. 

Included in RINGO, there are two ways of generating ancestral adjacency weights. One is to use [DeClone](https://github.com/yannponty/DeClone) [3]. After installing DeClone, indicate its full path in the `ringo.cfg` file, in the section `[Paths]`. Then, running the `weighting_DeClone.py` as described below will generate a file that can be used by `ringo.py`.

Another way of determining adjacency weights is by using the algorithm proposed in [2]. This is the default way, and if no adjacency file is given to `ringo.py`, it will calculate the adjacency weights with this method. It is also possible to generate an adjacency weights file using the `ancestral_weights.py` script as described below.

##### `weighting_DeClone.py`
```
usage: weighting_DeClone.py [-h] -nf NEWICK [-i] [-sm SET_MINIMUM] -m MARKERS
                            [-kT KT] [-o OUTPUT_FOLDER]

Weights given tree in nhx-format with DeClone. Also converts tree in NEWICK-
format into tree in nhx-format

optional arguments:
  -h, --help            show this help message and exit
  -nf NEWICK, --Newick NEWICK
                        path to the file with NEWICK-tree
  -i, --ignore_weights  boolean, for either ignore or consider edge
                        length/weights, when parsing Newick Tree into nhx Tree
  -sm SET_MINIMUM, --set_minimum SET_MINIMUM
                        minimal value for any edge length, when parsing Newick
                        Tree into nhx Tree
  -m MARKERS, --markers MARKERS
                        path to marker-file
  -kT KT                deClone constant
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Folder for output files
```

##### `ancestral_weights.py`

```
usage: ancestral_weights.py [-h] -i INPUT_GENOMES -tree TREE -o OUTPUT

Generates an adjacency weight file for internal node adjacencies, given a
newick tree and the leaf genomes file.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_GENOMES, --input_genomes INPUT_GENOMES
                        Leaf genomes file.
  -tree TREE            Newick Tree file.
  -o OUTPUT, --output OUTPUT
                        Output file.
                        
```

## Ancestral Reconstruction with RINGO
##### `ringo.py`

```
usage: ringo.py [-h] [-i INPUT_GENOMES] [-t TREE] -o OUTPUT
                [-w ADJ_WEIGHTS_FILE] [-f WEIGHT_FILTER] [-a ANCESTRAL]

Solves the Small Phylogeny problem with I.G. InDel

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_GENOMES, --input_genomes INPUT_GENOMES
                        Leaf genomes file.
  -t TREE, --tree TREE  Newick Tree file.
  -o OUTPUT, --output OUTPUT
                        Output folder.
  -w ADJ_WEIGHTS_FILE, --adj_weights_file ADJ_WEIGHTS_FILE
                        internal weights file.
  -f WEIGHT_FILTER, --weight_filter WEIGHT_FILTER
                        Filter cutoff for adjacency weights, smaller weights
                        are removed.
```


## Generating simulated datasets
##### `simulation.py`

```
usage: simulation.py [-h] [-n NUM_GENES] [-c NUM_CHR] [-o OUTPUT] [-i INDEL]
                     [--indel_length INDEL_LENGTH]
                     [-r RATE | -sc SCALE | -e EVENTS_PER_EDGE]
                     (-f FILE | -s SIM) [-d DISTURB]

Simulates rearrangement evolution on a given newick tree

optional arguments:
  -h, --help            show this help message and exit
  -n NUM_GENES, --num_genes NUM_GENES
                        Number of genes in the root genome.
  -c NUM_CHR, --num_chr NUM_CHR
                        Number of chromosomes in the root genome.
  -o OUTPUT, --output OUTPUT
                        Name of the output folder.
  -i INDEL, --indel INDEL
                        Percentage of indels, from 0 to 1.0
  --indel_length INDEL_LENGTH
                        Maximum size of indel event in genes.
  -r RATE, --rate RATE  Multiplier on the input tree number of events.
  -sc SCALE, --scale SCALE
                        Scales the tree so each leaf has on average scale/2
                        *n_genes number of events. Almost the same as
                        diameter, if tree is ultrametric.
  -e EVENTS_PER_EDGE, --events_per_edge EVENTS_PER_EDGE
                        Events on edge: random [e/2, e]
  -f FILE, --file FILE  Input a Newick tree
  -s SIM, --sim SIM     Simulate a new birth_death with SIM species
  -d DISTURB, --disturb DISTURB
                        Disturb branch lengths multiplying each by e^r, where
                        r in [-d,+d].
```


## Evaluating reconstruction quality


# File Formats

## Genome formats

## Tree format

