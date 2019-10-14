# Changelog

## Version 1.0.4

Fixed logic error that arose from referencing dictionary of rules instead of initialAssignments when properly ordering rules and assignments when generating modules.

When a parameter or species does not have a value assigned in the sbml file (for instance, when the value is assigned by rule or initialAssignment), None is now used as a placeholder instead of the default of nan used by libSBML.

## Version 1.0.3

The modelclasses submodule was previously only importedwithin the module when utitilies.TestPackage() was called. This lead to the submodule only being accessible by a user after calling TestPackage(). Now the submodule is accessable imediately after importing.

Fixed mathstrings for user defined functions and reaction not being properly converted to use Python functions.

When generating the method to solve for a reaction's velocity, the code will now check if an id refers to a reaction parameter before checking global model parameters.