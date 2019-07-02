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

For more information about this package, see **Insert link to documentation here**
For this project's github repository, [click here](https://github.com/SMRuggiero/sbmltopyode)
For the publication associated with this package, see **Insert link to publication here**