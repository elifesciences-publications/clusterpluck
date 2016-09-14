clusterpluck
=======
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
With the shogun environment activated, install the clusterpluck package and the NINJA-utils developmental tools.

```
# clusterpluck
pip install git+https://github.com/RRShieldsCutler/clusterpluck.git --no-cache-dir --upgrade

# NINJA-utils
pip install git+https://github.com/knights-lab/NINJA-utils.git --no-cache-dir --upgrade
```
If you do not have pandas installed, install (or update) pandas
```
conda install pandas
```

