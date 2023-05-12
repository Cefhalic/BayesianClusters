# Bayesian Clusters

<!-- ## The magic Makefile

There is a makefile to handle most boilerplate tasks for you. The makefile has a number of options:

  * `make`                - Build code
  * `make help`           - Display a help message
  * `make clean`          - Tidy all build products
  * `make all`            - Build code, generate doxygen documentation and produce PDFs of latex sources
  * `make cpp`            - Build code
  * `make verbose`        - Build code, echoing the full command
  * `make doxygen`        - Generate doxygen documentation
  * `make docs`           - Produce PDFs of latex sources

### To build the executables

Prerequisites:
 * make: https://www.gnu.org/software/make/
 * g++: https://www.gnu.org/software/gcc/
 * CERN ROOT: https://root.cern/
 * BOOST: https://www.boost.org/

To build the executables: `make` or `make cpp` or `make verbose`. The latter does a full echo of the commandline, whilst the former two display a user-friendly sanitized version.

The executables are then at `./Cluster.exe` and `./Display.exe`

### To build the doxygen documentation 

Prerequisites:
 * make: https://www.gnu.org/software/make/
 * doxygen: https://doxygen.nl/

To build the doxygen documentation: `make doxygen`

The doxygen documentation is then at `documentation/SoftwareManual.pdf`.

**However, most of the time you should never need to, as the documentation is also built as a CI pipeline and published [here](https://github.com/Cefhalic/BayesianClusters/raw/master/documentation/SoftwareManual.pdf).**

### To build the maths documentation 

Prerequisites:
 * make: https://www.gnu.org/software/make/ 
 * pdflatex: https://linux.die.net/man/1/pdflatex

To build the maths documentation: `make docs`

The maths documentation is then at `documentation/OptimizingTheMaths.pdf`

## Scan.exe
### To see help
```
./Scan.exe --help
```
### To run an RT-scan without any output (for timing)
```
./Scan.exe --cfg example-configs/config.txt -i 1_un_red.csv
```

### To run an RT-scan with JSON or XML output
```
./Scan.exe --cfg example-configs/config.txt -i 1_un_red.csv -o ScanResults.json
```
or
```
./Scan.exe --cfg example-configs/config.txt -i 1_un_red.csv -o ScanResults.xml
```
Please note - the file can have any name you please, but the extension must be `.xml` or `.json` and is case sensitive.

## Display.exe

### To run the event display
To show the whole data set
```
./Display.exe -i 1_un_red.csv
```
To show only the data in the specified region of interest
```
./Display.exe --cfg example-configs/config.txt -i 1_un_red.csv
```
-->

<!--
# Instructions for running on Imperial HPC from Python

Log into a login node on HPC or onto the HPC JupyterHub (https://jupyter.rcs.imperial.ac.uk/) and open a terminal window.  The instructions below are printed out from my command history.

## Create an anaconda environment

```
bash-4.4$ module load anaconda/personal
bash-4.4$ conda create --name bayesian python=3.8
bash-4.4$ source activate bayesian
(bayesian) bash-4.4$ conda env update --file utilities/environment.yml --prune
```

## Download the code
```
(bayesian) bash-4.4$ git clone https://github.com/Cefhalic/BayesianClusters.git Bayesian
(bayesian) bash-4.4$ cd Bayesian/
(bayesian) bash-4.4$ git checkout PythonBindings
```

## Compile the code
```
(bayesian) bash-4.4$ cd Bayesian/
(bayesian) bash-4.4$ make clean
(bayesian) bash-4.4$ make
```
-->

# Instructions for running on Imperial HPC from Python

Log into a login node on HPC or onto the HPC JupyterHub (https://jupyter.rcs.imperial.ac.uk/) and open a terminal window. 

## Download the code
```
git clone https://github.com/Cefhalic/BayesianClusters.git Bayesian
cd Bayesian/
```

## Create an anaconda environment
Should only need to be run once ever!
```
module load anaconda/personal
conda create --name bayesian python=3.8
conda activate bayesian
conda env update --file utilities/environment.yml --prune
```

## Otherwise, activate an existing anaconda environment
```
module load anaconda/personal
conda activate bayesian
```

## Compile the code
```
make clean
make -j8
```

## Test the code.  

Note that this will not produce any graphical output unless there is an X-Windows or another graphical interface.

```
(bayesian) bash-4.4$ python ./python/test.py --cfg example-configs/config.txt -i /rds/general/project/easystorm/live/bayesian/1_un_red.csv --r 20nm --t 40nm
``` 

## Running the test as a Jupyter notebook

While in a JupyterHub session first set up the conda `bayesian` environment to make the current conda environment available as `bayesian` within Jupyter (the modules should already be included with the environment.yml file above).

```
(bayesian) bash-4.4$ python -m ipykernel install --user --name bayesian --display-name "bayesian"
```

From the `python` subdirectory of the main `Bayesian` source code directory open the `test.ipynb` Jupyter notebook.

Enjoy!
