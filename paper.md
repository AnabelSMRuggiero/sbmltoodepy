---
title: 'SMBLtoODEpy: A Software Program for Converting SBML Models into ODE Models in Python'
tags:
  - Python

authors:
  - name: Steve M. Ruggiero
    orcid: 
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Ashlee N. Ford Versypt
    orcid:
	affiliation: 1

affiliations:
 - name: Oklahoma State University
   index: 1

date: 31 July 2019
bibliography: paper.bib
---

# Summary


Systems Biology Markup Language (SBML) is a standard intermediate language for representing models of biological systems,
particularly reaction networks [@hucka2003][@finney2003].
SBML defines models independent of the software used to create the model, which allows for easy import and export of models.
Programs for converting models into code exist for a variety of software languages but have key limitations.
For example, the Systems Biology Format Converter (SBFC), a suite of tools used to convert SBML models to various modeling and programming languages, is a tool used to convert SBML into scripts in multiple programming languages [@rodriguez2016].
The BioModels Database uses SBFC to automate conversion of SBML models into other formats that can be downloaded by users [@chelliah2014][@finney2003].
SBFC does not include a tool for generating Python scripts. For the languages SBFC supports, the codes generated are stand alone implementations of models that not easy to integrate into other projects.  
The present work describes a software program called SBMLtoODEpy taht we developed to address these limitations by enabling conversion of SBML models into Python classes that can be rapidly incorporated into biomedical systems modeling projects written in Python,
such as the multiscale simulation platform CompuCell3d.
The program aims to accelerate construction of multiscale models that import and reuse published SBML models,
many of which are available in the BioModels Database at https://www.ebi.ac.uk/biomodels/

In SMBLtoODEpy, each of the model components are extracted using libSBML,
a software library for parsing and editing SBML models [@bornstein2008].
The model components can be output to a JSON file. JSON is a format that is easier to directly edit than either the SBML or Python implementations of the model.
The extracted model components is used to create a Python file that defines a class which implements the model.
A class method is generated to solve the model using a wrapper for the lsoda algorithm in the SciPy Python package [@oliphant2007].
Also used in this package is the NumPy Python package [@van2011].


# Acknowledgements

The research results discussed in this publication were made possible in part by funding through the award for project number HR17-057, from the Oklahoma Center for the Advancement of Science and Technology.


# References