# Bayesian Clusters

## To build the executable and generate the documentation
```
make
```
 * The executable is then at `./Cluster.exe` 
 * The doxygen documentation is at `doxygen/html/index.html`
 * The maths documentation is at `documentation/OptimizingTheMaths.pdf`

## To see help
```
./Cluster.exe --help
```
## To run an RT-scan without any output (for timing)
```
./Cluster.exe --cfg config.txt -i 1_un_red.csv
```

## To run an RT-scan with JSON or XML output
```
./Cluster.exe --cfg config.txt -i 1_un_red.csv -o ScanResults.json
```
or
```
./Cluster.exe --cfg config.txt -i 1_un_red.csv -o ScanResults.xml
```
Please note - the file can have any name you please, but the extension must be `.xml` or `.json` and is case sensitive.

## To run the event display
To show the whole data set
```
./Display.exe -i 1_un_red.csv
```
To show only the data in the specified region of interest
```
./Display.exe --cfg config.txt -i 1_un_red.csv
```
