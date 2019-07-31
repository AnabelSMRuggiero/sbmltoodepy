# SBMLtoODEpy

## Overview

The simpliest use of this package is to use the ParseAndGenerate function to quickly create a python implementation of an SBML model.

```python
import sbmltoodepy
ParseAndCreateModel(inputFilePath, jsonFilePath = None, outputFilePath = None, className = "SBMLmodel")
```

The parameters are:
* inputFilePath: The file path of the SBML model.
* jsonFilePath: An optional file path that if provided is where the function will create a json file containing all of the model elements. If not provided, a json file will not be created.
* outputFilePath: An optional file path of where the python model implementation will be created. If not provided, the output path will be assumedtobe the same as the input, but with a .py extension.
* className: The name of the class defined by the file created by this function.


This creates a new python file containing a class implementing the SBML model. 
To run the model, instantiate the class and call the RunSimulation method with the desired timestep.

```python
model = SBMLmodel()
model.RunSimulation(deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6)
```

## Installation Instructions

SBMLtoODEpy can be downloaded using pip
```
pip install sbmltoodepy
```

To test the package, use the TestPackage function.

```python
import sbmltoodepy
sbmltoodepy.utilities.TestPackage()
```

The function raises a warning about trying to set a constant species and 7 numbers.
        
1.279495335622681e-06

9.427471779552446e-10

1.330807667219949e-07

9.797580656559267e-08

4.0569605652344643e-07

2.661173190172687e-05
    
These are the average relative errors for species, parameters, and compartments between the models generated in this function and results calculated by COPASI for seven different models.

## Acknowledgements

The following python packages were used in this project and are needed to generate and run python models:
* [NumPy][1]
* [SciPy][2]
* [libSBML][3]

Models from the following papers were used to develop and test the package and are included as examples:

[Borisov, Nikolay, et al. "Systems-level interactions between insulin–EGF networks amplify mitogenic signaling." Molecular systems biology 5.1 (2009): 256.][4]

[Guyton, Arthur C., Thomas G. Coleman, and Harris J. Granger. "Circulation: overall regulation." Annual Review of Physiology 34.1 (1972): 13-44.][5]

[Kerkhoven, Eduard J., et al. "Handling uncertainty in dynamic models: the pentose phosphate pathway in Trypanosoma brucei." PLoS Computational Biology 9.12 (2013): e1003371.][6]

[Smallbone, Kieran, and Bernard M. Corfe. "A mathematical model of the colon crypt capturing compositional dynamic interactions between cell types." International Journal of Experimental Pathology 95.1 (2014): 1-7.][7]

[Waugh, Helen V., and Jonathan A. Sherratt. "Macrophage dynamics in diabetic wound dealing." Bulletin of Mathematical Biology 68.1 (2006): 197-207.][8]

[Zi, Zhike, et al. "Quantitative analysis of transient and sustained transforming growth factor‐β signaling dynamics." Molecular Systems Biology 7.1 (2011): 492.][9]

The above models are additionally located in [this projects github repository.](https://github.com/SMRuggiero/sbmltoodepy/tree/master/sbmltoodepy/sbml_files)

For more information on using this package, [click here to see the tutorial.](https://github.com/SMRuggiero/sbmltoodepy/blob/master/Tutorial.md)

For more information on SBML, including specifications and other software that supports SBML, [click here for the SBML web site](http://sbml.org/Main_Page)

For a good source of SBML models, [the BioModels repository is a great place to search.](https://www.ebi.ac.uk/biomodels/)

For this project's documentation, [click here.](https://sbmltoodepy.readthedocs.io/en/latest/)

For this project's github repository, [click here.](https://github.com/SMRuggiero/sbmltopyode)

[1]: https://www.numpy.org/ "NumPy"
[2]: https://www.scipy.org/ "SciPy"
[3]: http://sbml.org/Software/libSBML "libSBML"
[4]: https://doi.org/10.1038/msb.2009.19
[5]: https://doi.org/10.1146/annurev.ph.34.030172.000305
[6]: https://doi.org/10.1371/journal.pcbi.1003371
[7]: https://dx.doi.org/10.1111%2Fiep.12062
[8]: https://doi.org/10.1007/s11538-005-9022-3
[9]: https://doi.org/10.1038/msb.2011.22
