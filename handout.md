# How to build a phylogeny
*If this isn't rendered properly, see it [online](https://github.com/DomBennett/Teaching-Phylogeny-Building/blob/master/handout.md)*

I've created a framework for generating a phylogeny-building pipeline in Python.
The pipeline takes taxonomic names and produces phylogenetic trees. It has 3
steps: download, align and estimate.

In order to get the pipeline working you must fill in the missing elements by
researching how different functions and programs work. In the process of doing
this you will come to learn about what's needed to create a phylogeny, what
programs there are and how to use them, and how to build a programmatic
pipeline. By the end, you will have a pipeline that you can call your own and
adapt to suit whatever need you require, phylogeny-building or any related
discipline.

## Pre-requisites
* Basic Python programming ability
* A computer with Python 2.7 installed.

## Aims
1. Know the major programs for phylogeny generation
2. Be confident at installing and running these programs
3. Know how to handle and run programs through Python
4. Understand the basics for creating programmatic pipelines

## 1. Setup
To get started download the git repo in a suitable directory:

`git clone https://github.com/DomBennett/Teaching-Phylogeny-Building.git`

Have a look at the contents of the cloned folder. You will find a Python script
for each stage of the pipeline: download, align and estimate. You will also find
`names.txt` which contains the names of the plant species with which we want to
make a phylogeny, also you will find `expected_tree.tre` which is a Newick tree
file of the tree your pipeline should hopefully make.

Additionally, you will find a `tools.py` for holding home-made functions,
`run.py` for running all the stages and `0_check.py`. This last script checks
that you have all the relevant packages and dependencies installed. Run it from
your terminal:

`python 0_check.py`

It will tell you if you have any missing packages. To install any Python
packages that may be missing, first ensure you have [pip](https://pip.pypa.io/en/latest/installing.html)
installed and then type:

`sudo pip install [name-of-package]`

To install any missing programs follow one of these guides:
* [RAxML installation](https://github.com/stamatak/standard-RAxML)
* [MAFFT installation](http://mafft.cbrc.jp/alignment/software/source.html)

Installation is probably the hardest and most annoying part of scientific
computing, have faith and keep trying or ask for help.

`0_check.py` also tells you the command for you RAxML installation, make note of
it as we will need it later.

Now on to the real work ....

## 2. Download
Open `1_download.py` in your favourite text-editor. Study the script and make
sure you understand each line. Every stage script will have the same structure
so your effort will pay dividends. The aim for the script is to search and
download [rbcL](http://www.uniprot.org/uniprot/O03042) gene sequences for each
name from NCBI's GenBank. If you don't understand something, either ask or
run different elements of the script in your terminal. *Always break down and
test code*.

To get the script working, you must:
* Provide your email address to Entrez (the [API](http://en.wikipedia.org/wiki/Application_programming_interface) for NCBI)
* Read in the names as a list
* Get `download()` in `tools.py` working by researching how Entrez tools work

For that last step here are some useful links:
* http://www.ncbi.nlm.nih.gov/books/NBK49540/
* Also have a play at testing your terms here: http://www.ncbi.nlm.nih.gov/nuccore
* http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc111
* http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc114

After you have got the script working, run it via terminal:

`python 1_download.py`

## 3. Align
To get this script working you need to add the MAFFT command-line function in
the `align()` in `tools.py`. It is a good a idea to try using MAFFT via terminal
to see how it works before running it through Python. There are a lot of options
for running MAFFT, as there is for lots of phylogenetic programs, but it is fine
to use the auto mode.

Links:
* http://mafft.cbrc.jp/alignment/software/
* http://biopython.org/DIST/docs/api/Bio.Align.Applications._Mafft.MafftCommandline-class.html

Obviously, there are lots of other alternatives to MAFFT for aligning sequences.
I have chosen to use it here because it is easy to use, fast and accurate. It
stands for 'Multiple Alignment Fast-Fourier Transformation' and works through
predicting the molecular structure of sequences, although details are beyond me.
Alternative mainstream multiple alignment software you can use:
* Muscle
* ClustalW
* T-Coffee

## 4. Estimate
To get this last script working, you need to recreate the same format of the
stage script like the previous two. Also, you need to repeat what you did for
`align()` for `estimate()` using the Biopython RAxML command-line. Again test
out RAxML via terminal, use '-h' to see its commands.

RAxML is probably the best ML method for creating trees, it uses in-built
heuristics to make the tree-building as fast and accurate as possible.

The main contender for this hallowed spot for tree-building software is BEAST,
which is the main Bayesian method for tree building. We're using RAxML in this
instance because it's a lot faster.

## 5. Finishing-up
Now you should have 3 different scripts which can be run individually to
generate results for each stage of your pipeline. Look at `run.py` to run them
all together. Don't worry about how `runstage()` works too much (unless you
want to), the script only needs a tiny change to get it to work.

You should be able to run your entire pipeline with `python run.py`.
Compare your tree to the expected tree using FigTree (see link in Useful
Resources) -- how does it look? What's different? What could you do to your
pipeline to improve the tree? Note, you can use FigTree to root an unrooted
tree.

To compare your answers to mine you can browse my code [here](https://github.com/DomBennett/Teaching-Phylogeny-Building/tree/answers).
Or you can download my answers with this git command:

`git clone https://github.com/DomBennett/Teaching-Phylogeny-Building.git -b answers _answers`

This will download my pipeline and put it in a folder called `_answers`.

## 6. Extras challenges
Want do do more? Here are some extra things you could do to improve your
pipeline. Implement them in whatever order you want. If you're fed-up of making
your own pipeline, check out [mine](https://github.com/DomBennett/pG-lt) --
it's an extension of what you've just made.

### 6.1 Taxonomic constraint and outgrouping
Before running the pipeline, you could generate a taxonomic tree. RAxML let's
you define an outgroup before tree generation in order to generate a rooted
tree. It also allows you to provide a taxonomic constraint to prevent your
estimated tree from being too far from reality. Or at the very least, start
your tree estimation process with the taxonomic tree to speed RAxML up.

To implement a taxonomic tree, check out my package [taxon-names-resolver](https://github.com/DomBennett/TaxonNamesResolver/wiki/How-to-use%3F).
It has a function for generating taxonomic trees. Once you've got a taxonomic
tree, you can get the most isolated tip using something like this:

```{python}
# PACKAGES
from Bio import Phylo

# INPUT
with open('taxonomic.tre', 'r') as f:
    tree = Phylo.read(f, 'newick')

# PROCESS
# get distance to nearest tip, commensurate with taxonomic isolation
distances = [tree.distance(e) for e in spp]
index = [i for i, e in enumerate(distances) if e == min(distances)]
# choose one at random
ri = random.sample(index, 1)[0]
outgroup = spp[ri]
```

### 6.2 Run things multiple times
Sometimes the tree turns out well, other times not so. Can you adapt the
pipeline so that it downloads more than 1 sequence per species, it generates
multiple alignments and multiple phylogenies?

If you do this you will create multiple phylogenetic trees. To convert these
into a single tree use Dendropy's [`sumtrees.py`](https://pythonhosted.org/DendroPy/scripts/sumtrees.html) -- this should be in your system path already.

### 6.3 Work with more than just rbcL
At the moment we're only searching for rbcL sequences, but only plants have
rbcL sequences! Why not try names from other taxonomic groups using a different
gene? cytb works quite well for mammals.

Also, how would you generalise this so that you could give your pipeline any
set of names?

## Useful Resources
* [Biopython](http://biopython.org/DIST/docs/tutorial/Tutorial.html) -- the
package for handling bioinformatic data
* [Dendropy](https://pythonhosted.org/DendroPy/)  -- the best package for
handling phylogenies in python
* [pyCogent](http://pycogent.org/) -- package for handling genomic data
* [phyloGenerator](http://willpearse.github.io/phyloGenerator/) -- Will
Pearse's python program for generating phylogenies
* [pG-lt](https://github.com/DomBennett/pG-lt) -- 'piglet' or
phyloGenerator-lite, my own one-click Python package for turning names into
phylogenies
* [SUPERSMART](http://www.supersmart-project.org/) -- another one-click
pipeline, not quite ready yet but has huge potential.
* [phylogeny.fr](http://www.phylogeny.fr/) -- French website with online
phylogeny generating resources, useful for sanity-checking.
* [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) -- the best phylogeny
visualisation software I know of.
* [GlobalNamesResolver](http://resolver.globalnames.biodinfo.org/)(GNR) --
fuzzy-search for taxonomic names
* [TaxonNamesResolver](https://github.com/DomBennett/TaxonNamesResolver/wiki)
-- my own Python package for talking to GNR.
