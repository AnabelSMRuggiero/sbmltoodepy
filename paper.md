---
title: 'SMBLtoODEpy: A software program for converting SBML models into ODE models in Python'
tags:
  - Python
  - Systems Biology
  - Modeling

authors:
  - name: Steve M. Ruggiero
    affiliation: 1
  - name: Ashlee N. Ford Versypt
    orcid: 0000-0001-9059-5703
    affiliation: 1

affiliations:
 - name: School of Chemical Engineering, Oklahoma State University
   index: 1

date: 1 August 2019
bibliography: paper.bib
---

# Summary


Systems Biology Markup Language (SBML) is a standard intermediate language for representing models of biological systems,
particularly reaction networks [@hucka2003; @finney2003].
SBML defines models independent of the software used to create the model, which allows for easy import and export of models.
Programs for converting models into code exist for a variety of software languages but have key limitations.
For example, the Systems Biology Format Converter (SBFC) is a suite of tools used to convert SBML models into scripts in multiple modeling and programming languages [@rodriguez2016].
The BioModels Database uses SBFC to automate conversion of a large library of SBML models into other formats that can be downloaded by users [@chelliah2015; @glont2018].
SBFC does not include a tool for generating Python scripts. For the languages SBFC supports, the codes generated are stand alone implementations of models that not easy to integrate into other projects.  
The present work describes a software program called SBMLtoODEpy that we developed to address these limitations by enabling conversion of SBML models into Python classes that can be rapidly incorporated into biomedical systems modeling projects written in Python,
such as the multiscale simulation platform CompuCell3D, or used directly in Python.
The program aims to accelerate construction of multiscale models that import and reuse published SBML models,
many of which are available in the BioModels Database at https://www.ebi.ac.uk/biomodels/

In SMBLtoODEpy, each of the model components are extracted using libSBML,
a software library for parsing and editing SBML models [@bornstein2008].
The model components can be output to a JSON file. JSON is a format that is easier for users to read and directly edit than either the SBML or Python implementations of the model.
The extracted model components are used to create a Python file that defines a class that implements the model.
A method for the class is generated to solve the model using a wrapper for the lsoda algorithm in the SciPy Python package [@oliphant2007], and the NumPy Python package [@van2011] is also used. To verify that SBMLtoODEpy properly interprets SBML files and converts them into functional differential equations models, we successfully compared the performance of SMBLtoODEpy to COPASI, a graphical user interface based platform for simulating SBML models [@hoops2006], for a set of representative SBML files downloaded from the BioModels Database that were deposited for a selection of systems biology publications [@borisov2009; @guyton1972; @kerkhoven2013; @smallbone2014; @waugh2006; @zi2011]. These files have been included in the SBMLtoODEpy package within [the sbmltoodepy/sbml_files subdirectory](https://github.com/SMRuggiero/sbmltoodepy/tree/master/sbmltoodepy/sbml_files) to serve as examples for users. In our documentation, we have provided a tutorial on how to use the SBMLtoODEpy software package.

# Acknowledgments

The research results discussed in this publication were made possible in part by funding through the award for project number HR17-057, from the Oklahoma Center for the Advancement of Science and Technology.

# References