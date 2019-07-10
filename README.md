# sbmltopyode

The simpliest use of this package is to use the ParseAndGenerate function to quickly create a python implementation of an SBML model.

```python
import smbltoodepy
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
* [NumPy][1]
* [SciPy][2]
* [libSBML][3]

Models from the following papers were used to develop the model and are included as examples:

[Borisov, Nikolay, et al. "Systems-level interactions between insulin–EGF networks amplify mitogenic signaling." Molecular systems biology 5.1 (2009).][4]

[Cizmeci, Deniz, and Yaman Arkun. "Regulatory networks and complex interactions between the insulin and angiotensin II signalling systems: models and implications for hypertension and diabetes." PloS one 8.12 (2013): e83640.][5]

[Guyton, Arthur C., Thomas G. Coleman, and Harris J. Granger. "Circulation: overall regulation." Annual review of physiology 34.1 (1972): 13-44.][6]

[Kerkhoven, Eduard J., et al. "Handling uncertainty in dynamic models: the pentose phosphate pathway in Trypanosoma brucei." PLoS computational biology 9.12 (2013): e1003371.][7]

[Lenbury, Yongwimon, Sitipong Ruktamatakul, and Somkid Amornsamarnkul. "Modeling insulin kinetics: responses to a single oral glucose administration or ambulatory-fed conditions." Biosystems 59.1 (2001): 15-25.][8]

[Smallbone, Kieran, and Bernard M. Corfe. "A mathematical model of the colon crypt capturing compositional dynamic interactions between cell types." International journal of experimental pathology 95.1 (2014): 1-7.][9]

[Vizán, Pedro, et al. "Controlling long-term signaling: receptor dynamics determine attenuation and refractory behavior of the TGF-β pathway." Sci. Signal. 6.305 (2013): ra106-ra106.][10]

[Waugh, Helen V., and Jonathan A. Sherratt. "Macrophage dynamics in diabetic wound dealing." Bulletin of mathematical biology 68.1 (2006): 197-207.][11]

[Zi, Zhike, et al. "Quantitative analysis of transient and sustained transforming growth factor‐β signaling dynamics." Molecular systems biology 7.1 (2011).][12]

For more information about this package, see **Insert link to documentation here**

For this project's github repository, [click here](https://github.com/SMRuggiero/sbmltopyode)

For the publication associated with this package, see **Insert link to publication here**

[1]: https://www.numpy.org/ "NumPy"
[2]: https://www.scipy.org/ "SciPy"
[3]: http://sbml.org/Software/libSBML "libSBML"
[4]: https://doi.org/10.1038/msb.2009.19
[5]: https://doi.org/10.1371/journal.pone.0083640
[6]: https://doi.org/10.1146/annurev.ph.34.030172.000305
[7]: https://doi.org/10.1371/journal.pcbi.1003371
[8]: https://doi.org/10.1016/S0303-2647%2800%2900136-2
[9]: https://dx.doi.org/10.1111%2Fiep.12062
[10]: https://doi.org/10.1126/scisignal.2004416 
[11]: https://doi.org/10.1007/s11538-005-9022-3
[12]: https://doi.org/10.1038/msb.2011.22
