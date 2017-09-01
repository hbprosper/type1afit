# type1afit
Tutorials showing how to use the __[CERN](http://home.cern/)__
minimization Minuit written by Dr. Fred James, either directly or via
RooFit, to fit simple cosmological models to the Union 2.1 Type1a
supernova data. Minuit and the probability modeling package RooFit are
released with the CERN data analysis package
__[ROOT](http://root.cern.ch)__. The tutorials use the 

## Setup

   * Install the
__[Jupyter](https://root.cern.ch/root-has-its-jupyter-kernel)__ notebook package for
your operating system.
   * Install the packages histutil and type1afit
```
cd
mkdir -p external
cd external
git clone https://github.com/hbprosper/histutil

cd
mkdir -p tutorials
cd tutorials
git clone https://github.com/hbprosper/type1afit
cd type1afit
source $HOME/external/histutil/setup.sh
```
Then, run the command
```
jupyter notebook
```
in order to run the notebook in a browser. Navigate to the
desired notebook. (NOTE: your notebook program may be called ipython.)
