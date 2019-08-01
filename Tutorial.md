# Tutorial

## Creating a Python model

The fastest way to create a Python model using SBMLtoODEpy is to open a Python interpreter in the same folder as the SBML file. Then, run the following lines of code using the name of your model.

```python
import sbmltoodepy
sbmltoodepy.ParseAndCreateModel("YourModelNameHere.xml")
```

If no errors occur, then a file with the same name as the SBML file but with the extension .py will be created in the same folder. To have a file created with a different name, you can use the outputFilePath keyword argument.

```python
sbmltoodepy.ParseAndCreateModel("YourModelNameHere.xml", outputFilePath = "PythonFile.py")
```

The new Python module defines a class named SBMLmodel. To customize the name of the class, use the className keyword argument.

```python
sbmltoodepy.ParseAndCreateModel("YourModelNameHere.xml", outputFilePath = "PythonFile.py", className = "ModelName")
```

The last optional keyword argument of ParseAndCreateModel() is jsonFilePath. If specified, a JSON file containing all of the model elements parsed from the SBML model will be created.

```python
sbmltoodepy.ParseAndCreateModel("YourModelNameHere.xml", jsonFilePath = "YourModelNameHere.json", outputFilePath = "PythonFile.py", className = "ModelName")
```

In reality, creating a Python model with SBMLtoODEpy is a two step process that can be broken up.

```python
dataOfModel = sbmltoodepy.parse.ParseSBMLFile("YourModelNameHere.xml")
sbmltoodepy.modulegeneration.GenerateModel(dataOfModel, "PythonFile.py", objectName = "ModelName")
```

The function ParseAndCreateModel serves as a wrapper for two other functions, ParseSBMLFile and GenerateModel.
ParseSBMLFile returns an instance of the ModelData class with all of the components pulled from the model.
While you can use that instance as input to GenerateModel, there is another feature of SBMLtoODEpy that the ModelData class provides.
The ModelData class has a method, DumpToJSON, that generates a JSON file with the model components.

```python
sbmltoodepy.dataclasses.ModelData.DumpToJSON(dataOfModel,"YourModelNameHere.json")
```

A JSON file created from the DumpToJSON method can be used a create a new instance of the ModelData class.

```python
newInstance = sbmltoodepy.dataclasses.ModelData.LoadFromJSON("YourModelNameHere.json")
```

A possible use case for this would be to generate a JSON file from an SBML model, change the JSON file, which is more human readable than an SBML file, and use the changed JSON file to create the Python model.


## Exploring a newly created Python model

After generating the Python implementation of a model, you can now import the new class and instantiate it.

```python
from PythonFile import ModelName
modelInstance = ModelName()
```

The main components that make up the state of an SBML are compartments, parameters, and species. Each of these are stored in a dictionary that is a member of the model class. These members are named c, p,and s repectively.
The key for each dictionary entry is the identification (id), as defined by the SBML specification, for the component. By printing one of these dictionaries keys, you can see the id of each compartment, parameter, or species.

```python
# get the dictionary keys for the IDs of the species in the model
print(modelInstance.s.keys())
# get the dictionary keys for the IDs of the compartments in the model
print(modelInstance.c.keys())
# get the dictionary keys for the IDs of the parameters in the model
print(modelInstance.p.keys())
```

The individual model elements are represented as instances of an appropriate class. Parameters are instances of the Parameter class, and the value member contains the parameters value.
Species are instances of the Species class, and the concentration and amount members are used in calculations. Compartments are instances of the Compartment class, and the size member contains the size of the compartment.

```python
# replace compartmentId with one of the dictionary keys returned from print(modelInstance.c.keys())
print(modelInstance.c['compartmentId'].size) 
# replace parameterId with one of the dictionary keys returned from print(modelInstance.p.keys())
print(modelInstance.p['parameterId'].value)
# replace speciesId with one of the dictionary keys returned from print(modelInstance.s.keys())
print(modelInstance.s['speciesId'].concentration)
print(modelInstance.s['speciesId'].amount)
```

Function definitions and reactions are stored in a similar manner.
Currently, rules and initial assignments are not handled the same way.
Rate rules are bound methods of the model class that are called alongside reactions, and assignment rules and initial assignments are implemented through a method of the model class, AssignmentRules().
Algebraic rules are not completely supported currently.

In SBML, each model component can have both an id and a name, but only the id is required. 
Each of the dictionaries containing compartments, parameters, species, reations, and function definitions can be searched by the component name using the appropriate method. Note that the original model must not have a blank metadata field for the corresponding entry for the search by name to yield any results. In SBML, there is no guarantee that names are unique (each component in a model does have a unique id). If a single component matches, a list with the component's id and object is returned.
If multiple components match, then a nested list is returned. Each element in the list is a list with a component's id and object.
Lastly, if no match is found, the behavior of the function depends on the keyward argument suppress. If suppress is False, an exception is raised. If suppress is True, then an empty list is returned.

```python
modelInstance.SearchParametersByName('parameter name', suppress = False)
modelInstance.SearchCompartmentsByName('compartment name', suppress = False)
modelInstance.SearchSpeciesByName('species name', suppress = False)
modelInstance.earchReactionsByName('reaction name', suppress = False)
modelInstance.SearchFunctionsByName('function name', suppress = False)
```

The last part of the model's state is simulation time.

```python
print(modelInstance.time)
```

The time member of the class is the current simulation time for the model. The units of time in the model is treated as arbitrary by the package. The units for time is dependent on the units of parameters set by the modeler. This package does not check if the units in a model are consistent.

## Simulating the model

The RunSimulation() method is used to step the ordinary differential equations (ODE) model forward in time by timeinterval. The method updates the value of each model component and the total elapsed simulation time over all calls to RunSimulation. 

```python
modelInstance.RunSimulation(timeinterval)
```

There are two optional keyword arguments, absoluteTolerance and relativeTolerance, which are for tuning the tolerances used to solve the ODE model with lsoda.

```python
modelInstance.RunSimulation(timeinterval, absoluteTolerance = 1e-12, relativeTolerance = 1e-6)
```

An important thing to note is that the model only keeps values for the current time. If you wish to graph the evolution of a species with respect to time, the concentration and time will need to be sampled at the appropriate time points.

```python
import numpy as np
times = np.zeros(101)
times[0] = modelInstance.time
concentrations = np.zeros(101)
concentrations[0] = modelInstance.s['speciesId'].concentration
timeinterval = 1
for i in range(100):
	modelInstance.RunSimulation(timeinterval)
	times[i+1] = modelInstance.time
	concentrations[i+1] = modelInstance.s['speciesId'].concentration
```

The results can then be graphed using any method you would like. For graphing directly in python, [matplotlib is a good option.](https://matplotlib.org/)

```python
import matplotlib.pyplot as plt
plt.plot(times,concentrations)
```

## Examples
Models from the following papers were used to develop and test the package and are included as example SBML files in [the sbmltoodepy/sbml_files subdirectory.](https://github.com/SMRuggiero/sbmltoodepy/tree/master/sbmltoodepy/sbml_files)

[Borisov, Nikolay, et al. "Systems-level interactions between insulin–EGF networks amplify mitogenic signaling." Molecular systems biology 5.1 (2009): 256.][1]

[Guyton, Arthur C., Thomas G. Coleman, and Harris J. Granger. "Circulation: overall regulation." Annual Review of Physiology 34.1 (1972): 13-44.][2]

[Kerkhoven, Eduard J., et al. "Handling uncertainty in dynamic models: the pentose phosphate pathway in Trypanosoma brucei." PLoS Computational Biology 9.12 (2013): e1003371.][3]

[Smallbone, Kieran, and Bernard M. Corfe. "A mathematical model of the colon crypt capturing compositional dynamic interactions between cell types." International Journal of Experimental Pathology 95.1 (2014): 1-7.][4]

[Waugh, Helen V., and Jonathan A. Sherratt. "Macrophage dynamics in diabetic wound dealing." Bulletin of Mathematical Biology 68.1 (2006): 197-207.][5]

[Zi, Zhike, et al. "Quantitative analysis of transient and sustained transforming growth factor‐β signaling dynamics." Molecular Systems Biology 7.1 (2011): 492.][6]

[1]: https://doi.org/10.1038/msb.2009.19
[2]: https://doi.org/10.1146/annurev.ph.34.030172.000305
[3]: https://doi.org/10.1371/journal.pcbi.1003371
[4]: https://dx.doi.org/10.1111%2Fiep.12062
[5]: https://doi.org/10.1007/s11538-005-9022-3
[6]: https://doi.org/10.1038/msb.2011.22
