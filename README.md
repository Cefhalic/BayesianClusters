# Bayesian Clusters

## The magic Makefile

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

## Cluster.exe
### To see help
```
./Scan.exe --help
```
### To run an RT-scan without any output (for timing)
```
./Scan.exe --cfg config.txt -i 1_un_red.csv
```

### To run an RT-scan with JSON or XML output
```
./Scan.exe --cfg config.txt -i 1_un_red.csv -o ScanResults.json
```
or
```
./Scan.exe --cfg config.txt -i 1_un_red.csv -o ScanResults.xml
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
./Display.exe --cfg config.txt -i 1_un_red.csv
```
