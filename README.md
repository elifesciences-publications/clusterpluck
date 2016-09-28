clusterpluck
=======
#####A Python package designed to efficiently process antiSMASH(1) results from large genome databases, with the ability to group homologous clusters into Operational Functional Units (OFUs).

##Installation

These instructions are written for OSX systems (tested on 10.11.6). This package requires anaconda, which is a system agnostic package and virtual environment manager. Follow the installation instructions for your system at http://conda.pydata.org/miniconda.html.

Once anaconda is installed, create a new virtual environment with python3.

```
conda create -n clusterpluck python=3
```
Now activate the environment.
```
source activate clusterpluck
```
With the clusterpluck environment activated, install the clusterpluck package and the NINJA-utils developmental tools from @knights-lab.

```
# clusterpluck
pip install git+https://github.com/RRShieldsCutler/clusterpluck.git --no-cache-dir --upgrade

# NINJA-utils
pip install git+https://github.com/knights-lab/NINJA-utils.git --no-cache-dir --upgrade
```

If you do not have pandas installed, install (or update) pandas :panda_face:
```
conda install pandas
```
In order to fully annotate the genomes, clusterpluck uses the NINJA-DOJO database package. If you want annotation of OFUs beyond just NCBI accession numbers (e.g. NC_000913.3), install DOJO and NINJA-SHOGUN in the clusterpluck environment.

```
# NINJA-DOJO
pip install git+https://github.com/knights-lab/NINJA-DOJO.git --no-cache-dir --upgrade

# NINJA-SHOGUN
pip install git+https://github.com/knights-lab/NINJA-SHOGUN.git --no-cache-dir --upgrade
```

______________

1. Weber T, Blin K, Duddela S, Krug D, Kim HU, Bruccoleri R, Lee SY, Fischbach MA, Müller R, Wohlleben W, Breitling R. (2015). antiSMASH 3.0—a comprehensive resource for the genome mining of biosynthetic gene clusters. Nucleic acids research, 43(W1), W237-W243.
