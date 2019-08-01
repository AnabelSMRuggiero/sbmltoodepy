# SBMLtoODEpy

## Overview
SBMLtoODEpy is a software package that was developed to enable conversion of Systems Biology Markup Language (SBML) models into Python classes that can be rapidly incorporated into biomedical systems modeling projects written in Python, such as the multiscale simulation platform CompuCell3D, or used simulated as ordinary equations models in Python. 

## Authors
Steve M. Ruggiero and Ashlee N. Ford Versypt, School of Chemical Engineering, Oklahoma State University

## Usage
The simpliest use of this package is to use the ParseAndGenerate function to quickly create a python implementation of an SBML model.

```python
import sbmltoodepy
ParseAndCreateModel(inputFilePath, jsonFilePath = None, outputFilePath = None, className = "SBMLmodel")
```

The parameters are:
* inputFilePath: The file path of the SBML model.
* jsonFilePath: An optional file path that if provided is where the function will create a json file containing all of the model elements. If not provided, a json file will not be created.
* outputFilePath: An optional file path of where the python model implementation will be created. If not provided, the output path will be assumed to be the same as the input but with a .py extension.
* className: The name of the class defined by the file created by this function.

This creates a new python file containing a class implementing the SBML model. 

To run the model, instantiate the class and call the RunSimulation method with the desired timestep deltaT and optional specification of tolerances.

```python
model = SBMLmodel()
model.RunSimulation(deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6)
```

## Installation Instructions

SBMLtoODEpy can be downloaded using pip:
```
pip install sbmltoodepy
```

To test the package, use the TestPackage function.

```python
import sbmltoodepy
sbmltoodepy.utilities.TestPackage()
```

The TestPackage function raises a warning about trying to set a constant species and 7 numbers.
        
1.279495335622681e-06

9.427471779552446e-10

1.330807667219949e-07

9.797580656559267e-08

4.0569605652344643e-07

2.661173190172687e-05
    
These are the average relative errors for species, parameters, and compartments between the models generated in this function and results calculated by COPASI for seven different models.

## Acknowledgments

The following python packages were used in this project and are needed to generate and run python models:
* [NumPy][1]
* [SciPy][2]
* [libSBML][3]

For more information on using this package, [click here to see the tutorial.](https://github.com/SMRuggiero/sbmltoodepy/blob/master/Tutorial.md)

For more information on SBML, including specifications and other software that supports SBML, [click here for the SBML web site](http://sbml.org/Main_Page)

For a good source of SBML models, [the BioModels repository is a great place to search.](https://www.ebi.ac.uk/biomodels/)

For this project's code documentation, [click here.](https://sbmltoodepy.readthedocs.io/en/latest/)

For this project's github repository, [click here.](https://github.com/SMRuggiero/sbmltopyode)

[1]: https://www.numpy.org/ "NumPy"
[2]: https://www.scipy.org/ "SciPy"
[3]: http://sbml.org/Software/libSBML "libSBML"
