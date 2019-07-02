# sbmltopyode

The simpliest use of this package is to use the ParseAndGenerate function to quickly create a python implementation of an SBML model.

```python
ParseAndCreateModel(inputFilePath, jsonFilePath = None, outputFilePath = None, className = "SBMLmodel")
```

The parameters are:
* inputFilePath: The file path of the SBML model.
* jsonFilePath: An optional file path that if provided is where the function will create a json file containing all of the model elements. If not provided, a json file will not be created.
* outputFilePath: An optional file path of where the python model implementation will be created. If not provided, the output path will be assumedtobe the same as the input, but with a .py extension.
* className: The name of the class implementing the sbml model.


This creates a new python file containing a class implementing the SBML model. 
To run the model, instantiate the class and call the RunSimulation method with the desired timestep.

```python
model = SBMLmodel()
model.RunSimulation(deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6)
```

The following python packages were used in this project:
* numpy[1]
* libsbml[2]

Models from the following papers were used to develop the model and are included as examples:


For more information about this package, see **Insert link to documentation here**

For this project's github repository, [click here](https://github.com/SMRuggiero/sbmltopyode)

For the publication associated with this package, see **Insert link to publication here**

[1]: https://www.numpy.org/ "NumPy"
[2]: http://sbml.org/Software/libSBML "libSBML"