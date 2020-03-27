# Graph Component Segmentation
Read in ODGI Bin output and identify co-linear components



# Developer Instructions
**Environment**: [Anaconda 3.7 ](https://www.anaconda.com/distribution/)  
Ggfapy etc. does not have an anaconda package, so it's necessary to use pip:  
`pip install -r requirements.txt`  


**IDE:**  Pycharm Professional 2019.1  
* Travis CI - automatically runs tests on master and development branches

#### Example Run Configuration
Parameters: `--json-file=data/run1.B1phi1.i1.seqwish.w100.json --out-folder=data/ --cells-per-file=5000`  
In the terminal, run: `export PYTHONPATH=$PYTHONPATH:/home/ubuntu/software/component_segmentation/matrixcomponent` but replace with the path to your matrixcomponent directory.  This solves [ModuleNotFoundError #13](https://github.com/graph-genome/component_segmentation/issues/13).
