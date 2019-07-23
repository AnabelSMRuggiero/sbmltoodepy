# Tutorial

## Creating a Python model

The fastest way to create a Python model using SBMLtoODEpy is to open a Python interpreter in the same folder as the SBML file. Then, run the following lines of code using the name of your model.

```python
import smbltoodepy
ParseAndCreateModel("YourModelNameHere.xml")
```

If no errors occur, then a file with the same name as the SBML file but with the extension .py will be created in the same folder. To have a file created with a different name, you can use the outputFilePath keyword argument.

```python
ParseAndCreateModel("YourModelNameHere.xml", outputFilePath = "PythonFile.py")
```

The new Python module defines a class named SBMLmodel. To customize the name of the class, use the className keyword argument.

```python
ParseAndCreateModel("YourModelNameHere.xml", outputFilePath = "PythonFile.py", className = "ModelName")
```

The last optional keyword argument of ParseAndCreateModel() is jsonFilePath. If specified, a JSON file containing all of the model elements parsed from the SBML model will be created.

```python
ParseAndCreateModel("YourModelNameHere.xml", jsonFilePath = "YourModelNameHere.json", outputFilePath = "PythonFile.py", className = "ModelName")
```

## Exploring a newly created Python model

After generating the Python implementation of a model, you can now import the new class and instantiate it.

```python
from PythonFile import ModelName
modelInstance = ModelName
```

The main components that make up the state of an SBML are compartments, parameters, and species. Each of these are stored in a dictionary that is a member of the model class. These members are named c, p,and s repectively.
The key for each dictionary entry is the id, as defined by the SBML specification, for the component. By printing one of these dictionaries keys, you can see the id of each compartment, parameter, or species.

```python
print(modelInstance.s.keys())
```

The individual model elements are represented as instances of an appropriate class. Parameters are instances of the Parameter class, and the value member contains the parameters value.
Species are instances of the Species class, and the concentration and amount members are used in calculations. Compartments are instances of the Compartment class, and the size member contains the size of the compartment.

```python
print(modelInstance.c['compartmentId'].size)
print(modelInstance.p['parameterId'].value)
print(modelInstance.s['speciesId'].concentration)
print(modelInstance.s['speciesId'].amount)
```

The last part of the model's state is time.

```python
print(modelInstance.time)
```

## Simulating the model

The RunSimulation() method is used to step the model forward in time.

```python
modelInstance.RunSimulation(1)
```

The method updates the value of each model component and the time member. There are two keyword arguements, absoluteTolerance and relativeTolerance, which are for tuning the tolerances used to sole the model.

```python
modelInstance.RunSimulation(1, absoluteTolerance = 1e-12, relativeTolerance = 1e-6)
```

An important thing to note is that the model only keeps values for the current time. If you wish to graph the evolution of a species with respect to time, the concentration and time will need to be sampled at the appropriate time points.

```python
times = []
concentrations = []
for i in range(100):
	modelInstance.RunSimulation(1)
	times.append(modelInstance.time)
	concentrations.append(modelInstance.s['speciesId'].concentration)
```