# Assignment Problem with Diversity Constraints
This repository contains the code to reproduce the results from the paper **Assignment Problem with Diversity Constraints**(*in preparation*) and **Diversity Constraints in Public Housing Allocation**([*AAMAS 2018*](http://ifaamas.org/Proceedings/aamas2018/pdfs/p973.pdf)).

Authors: Nawal Benabbou, Mithun Chakraborty, Xuan-Vinh Ho, Jakub Sliwinski, and Yair Zick 

For any issue regarding this repository, please contact Nawal Benabbou(dcsnb@nus.edu.sg) or Xuan-Vinh Ho(hxvinh@comp.nus.edu.sg)

## Prerequisites
Our program is written in Java with available standard libraries, except for [Gurobi](http://www.gurobi.com/)(an optimization tool). Any issue in activating the license and using the library can be found [here](http://www.gurobi.com/documentation/7.5/quickstart_mac.pdf).

## Directory Structure

Our directory is organized as follows:
```bash
|-- bin
|-- Data
|   |-- chicago
|   |   `- HDB_data.xlsx
|   `-- housing
|   |   `- Chicago_data.xlsx
|-- Results
|-- src
|   |-- Agent.java
|   |-- MainTests.java
|   |-- ProcessResultsFile.java
|   `-- Solveur.java
```

All the simulated data are stored in `HDB_data.xlsx` and `Chicago_data.xlsx`, then being converted to corresponding `*.txt` files such that the program can load input from. 

To play with parameters and observe the changes in outcome, user has to modify directly in `MainTests.java` file. Specifically, you have to turn on **one** of the following boolean variables in _line 144-146_:
```java
        boolean HDB_realData = false; 
        boolean varyHDB = false;
        boolean chicago = false;
```

For each option, you are further required to modify related setting:
* HDB_realData (_line 151-160_): run experiments using the data collected from HDB website
```java
        boolean randomUtility = false;
	boolean locationBasedUtility = false;
	boolean locationEthnicityBasedUtility = false;
	boolean incomePriceBasedUtility = false;
        MainTests.variance = 10;
	int numTests = 100;
```
* varyHDB (_line 728-737_): run experiments with real data but varying the number of agents, types and blocks
```java
        boolean varyNumberAgents = false;
	boolean varyNumberBlocks = false;
	boolean varyNumberTypes = false;
        boolean locationBased = true;
        MainTests.variance = 10;
	int numTests = 100;
```
* chicago (_line 1623-1627_): run experiments using the data obtained from the Chicago public school allocation website
```java
        int numberOfChildren = 5000;
        int numTests = 100;
        double[] vars = {0,10,50,100};
```

## Result

After the simulation finishs, all results are stored in `Results` folder. Each running produces 3 types of file: `error-*.txt`, `*-var-*.txt` and `*-octave.txt`. 
* `error-*.txt` shows the standard deviation and standard error of utilitarian social welfare in optimal solution under constraint, solution in lottery unconstrained optimal solution.
* `*-var-*.txt` shows the value of unconstrained optimal solution, PoD of constrained optimal solution, PoD of lottery solution, theoretical bound based on minAlpha, betaStar and theoretical bound based on betaStar, in each iteration.
* `*-octave.txt` shows the value of the above 6 values, in expectation.

Note that if user reruns the experiment with the same setting, previous result files should be moved to somewhere else, otherwise they will be overwritten. 
