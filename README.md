# type1afit
This tutorial shows how to use the __[CERN](http://home.cern/)__
minimization package Minuit (developed by __[Dr. Fred James](https://www.researchgate.net/profile/Fred_James2)__) to fit 
simple cosmological models to the [Union 2.1 Type1a](https://supernova.lbl.gov/Union/)
supernova data. 

## Dependencies
The notebook `type1afit.ipynb` depends on the following modules

| __modules__   | __description__     |
| :---          | :---        |
| cppyy         | calling of C++ code from Python |
| jupyterlab, notebook | Jupyter notebook environment
| matplotlib    | plotting module for high quality plots |
| scipy         | scientific computing    |
| pandas        | data table manipulation, often with data loaded from csv files |
| iminuit | a rewrite of the venerable CERN minimizer Minuit |


##  Installation
The simplest way to install these Python modules is first to install a software environment system. 
You could just bite the bullet and install Anaconda! However, it may be better to install
**miniconda3**, which is a very slim version of Anaconda, on your laptop. Do so by following the instructions at:

https://www.anaconda.com/docs/getting-started/miniconda/system-requirements


### Miniconda3

After installing miniconda3, it is a good idea to update conda using the command
```bash
conda update conda
```
#### Step 1 
Assuming conda is properly installed and initialized on your laptop, you can create an environment, here called *type1a* using the command
```bash
conda create --name type1a
```
and activate it by doing
```bash
conda activate type1a
```
You need create the environment only once, but you must activate the desired environment whenever you create a new terminal window.


#### Step 3
Install *jupyterlab*, *matplotlib*,  etc.
```bash
	conda install jupyterlab notebook
    conda install scipy
	conda install matplotlib
    conda install iminuit
```
Again be sure to check the exact syntax. This does change from time to time!
