# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 17:28:34 2018

@author: Steve
"""

import warnings

class Compartment:

#    members:
#        size (float) - This value represents the size of the shape. 
#            Whether this is a volume or an area is determined by the compartments dimensionality.
#            If the compartment is set to be constant, assigning a value to this compartment fails
#            and a warning is raised.
#		
#        dimensionality (float or int) - The number of spatial dimensions this compartment exists in.
#            Not necessarily a whole number.
#            
#        metadata (SBMLMetadata) - Information about the model component stored in the SBML file
#        
#        _constant (boolean) - 
#            
#        
#	
#    methods:
#        __init__(self, size, dimensionality[, isConstant = False, metadata = None])
#            arguements:
#                size (float) - Sets the value of the Size member
#				
#                dimensionality (float or int) - Sets the value of the dimensionality member
#                
#                isConstant (boolean)
#				
#            returns:
#                An instance of this class
    """This class represents a compartment which is one of the model elements in SBML.
    
    Each species is associated with one compartment by referencing the appropriate 
    compartment with its own compartment member.
    
    Parameters
    ----------
    size : float
        This value represents the size of the shape. 
        Whether this is a volume or an area is determined by the compartments dimensionality.
        If the compartment is set to be constant, assigning a value to this compartment fails
        and a warning is raised.
	
    dimensionality : int or float
        The number of spatial dimensions this compartment exists in. This is not limited to whole numbers.
            
    isConstant : bool
        If true, the compartment will be set to constant and attempts to change the size attribute will raise warnings.
            
    metadata : SBMLMetadata
        An object containing information about the model component stored in the SBML file
        
    Attributes
    ----------
    dimensionality : float
        The dimensionality of the compartment
    metadata : SBMLMetadata
        Additional information about the compartment contained in the SBML model
    """
    
    def __init__(self, size, dimensionality, isConstant = False, metadata = None):

        self._size = size
        self.dimensionality = dimensionality
        self._constant = isConstant
        if metadata:
            self.metadata = metadata
            self.metadata._parentObject = self
        else:
            self.metadata = None
            
    @property
    def size(self):
        """float: The size of the comparment. This is used by species to calculate their concentrations."""
        return self._size
        
    @size.setter
    def size(self, size):
        if not self._constant:
            self._size = size
        else:
           warnings.warn("Attempted to set a constant parameter", Warning) 
        

#Currently unused
class Reaction:
    
    def __init__(self, metadata = None):
        
        if metadata:
            self.metadata = metadata
            self.metadata._parentObject = self
        else:
            self.metadata = None
        return
        
    
        
        
class Species:
    """This class represents a species as defined by SBML.
    
    Parameters
    ----------
    value : float
        Sets either the value of conentration or amount, depending on
        the concentrationIn value.

    valueType : str
        If valueType is equal to"Concentration", the amount input is assigned to the concentration member.
        Otherwise, amount is assigned to the amount member.
				
    compartment : Compartment
        Sets the compartment member to reference this input.
    
    hasOnlySubstanceUnits : bool
        Determines if the species concentration or amount is used is rate laws
				
    constant : bool
        If true, any attempts to modify the amount or concentration member fail.
    
    metadata : SBMLMetadata
        An object containing information about the model component stored in the SBML file
        
    Attributes
    ----------
    metadata : SBMLMetadata
        Additional information about the species contained in the SBML model
    """
    
    def __init__(self, value, valueType, compartment, hasOnlySubstanceUnits, constant = False, metadata = None):

        self._constant = constant
        self._hasOnlySubstanceUnits = hasOnlySubstanceUnits
        self._modifiedBy = None #If not none, this references a rule or event Id
        self._compartment = compartment
#        self._compartmentSize = compartment.size #This is needed to handle updating concentrations/amounts based on compartment resizing
        if metadata:
            self.metadata = metadata
            self.metadata._parentObject = self
        else:
            self.metadata = None
        if valueType == "Concentration":
            self._concentration = value
            self._amount = self._concentration * self._compartment.size

        else:
            self._amount = value
            self._concentration = self._amount / self._compartment.size
            
#    def UpdateCompartmentSizeMember(self):
#        self._compartmentSize = self._compartment.size
#        if (not self._modifiedBy == None) and not self._hasOnlySubstanceUnits: 
#            #Explicit check for self._modifiedBy not equalling None because who knows what id names has been picked
#            self._amount = self._concentration * self._compartmentSize
#            
#        else:
#            self._concentration = self._amount / self._compartmentSize
#                
    @property
    def concentration(self):
        """float : The concentration of the species.
        Direct changes to this member are automatically applied to amount as well.
        """
        self._concentration = self._amount/self._compartment.size
        return self._concentration
        
    @concentration.setter
    def concentration(self, concentration):
        if not self._constant:
            self._concentration =  concentration
            self._amount = concentration * self._compartment.size
#            self._amount = concentration * self._compartmentSize
        else:
            warnings.warn("Attempted to set a constant species", Warning) 
    #Todo: figure out how to make it so a user by default gets an exception, but can still get current behavior if requested
            
    @property
    def amount(self):
        """float : The total amount of the species that is present. 
        This volume is constant with respect to compartment size, unless the hasOnlySubstanceUnits
        member equals True and the species is modified by a rule or an event (events 
        are currently not supported). Changes to this member are automatically apply
        to concentration as well.
        """
        return self._amount
        
    @amount.setter
    def amount(self, amount):
        if not self._constant:
            self._amount =  amount
            self._concentration = amount / self._compartment.size
#            self._concentration = amount / self._compartmentSize
        else:
            warnings.warn("Attempted to set a constant species", Warning)
            
    @property
    def compartment(self):
        """Compartment : The compartment that this species is contained within."""
        return self._compartment
        
    @compartment.setter
    def compartment(self, compartment):
        warnings.warn("Attempted to change species compartment after species instantiation. If you know what you are doing, change species._compartment.")

            
class SBMLMetadata:
    
    def __init__(self, name = ''):
        self._parentObject = None
        self.name = name
        
class Parameter():
    """This class represents a parameter as defined by SBML.
    
    Parameters
    ----------
    value : float
        Sets either the value of the parameter.
        
    Id: str
        The Id of the parameter
				
    constant : bool
        If true, any attempts to modify the value member fail.
    
    metadata : SBMLMetadata
        An object containing information about the model component stored in the SBML file
        
    Attributes
    ----------
    metadata : SBMLMetadata
        Additional information about the species contained in the SBML model
    """
    
    def __init__(self, value, Id, constant = False, units = None, metadata = None):
        self._value = value
        self._constant = constant
        self._units = units
        self.metadata = metadata
        self.Id = Id
    
    @property
    def value(self):
        return self._value
    @value.setter  
    def value(self, value):
        if not self._constant:
            self._value = value
        else:
           warnings.warn("Attempted to set a constant parameter", Warning) 
#           
#class Test:
#    param = Parameter(1, constant = True)
#    def __init__(self, value):
#        self.param = Parameter(value, constant = True)
           
class Model():
    """This class is inhereted by model classes generated .
    
        
    Methods
    -------
    SearchParametersByName(name, suppress = False)
    
    SearchCompartmentsByName(name, suppress = False)
    
    SearchSpeciesByName(name, suppress = False)
    
    SearchReactionsByName(name, suppress = False)
    
    SearchFunctionsByName(name, suppress = False)
    
    SearchComponentsByName(componentDict, name, suppress = False)
        Searches a dictionary of SBML model components for any components that match the name argument.
        This method raises an exception if the suppress keyword argument is False and a match is not found.
        The other search methods (eg. SearchParametersByName) serve as wrappers for this method.
    """
    
#    def SearchElementsByName():
#        pass
    
    def SearchParametersByName(self, name, suppress = False):
        return self.SearchComponentsByName(self.p, name, suppress)
        
    def SearchCompartmentsByName(self, name, suppress = False):
        return self.SearchComponentsByName(self.c, name, suppress)
        
    def SearchSpeciesByName(self, name, suppress = False):
        return self.SearchComponentsByName(self.s, name, suppress)
        
    def SearchReactionsByName(self, name, suppress = False):
        return self.SearchComponentsByName(self.r, name, suppress)
        
    def SearchFunctionsByName(self, name, suppress = False):
        return self.SearchComponentsByName(self.f, name, suppress)
        
    def SearchComponentsByName(self, componentDict, name, suppress = False):
        returnList = []
        for key, component in componentDict.items():
            if component.metadata.name == name:
                returnList.append([key,component])
                
        if returnList == [] and not suppress:
            raise Exception('No matching name found')
                    
        if len(returnList) == 1:
            return returnList[0]
        else:
            return returnList
    

def Piecewise(*args):
    #    Piecewise(expression1, condition1 [, expression2, condition2 [,...]] [, otherwise])
    """
    Piecewise(expression1, condition1 [, expression2, condition2 [,...]] [, otherwise])
    
    This function implements the Piecewise function used in SBML models.
    
    Parameters
    ----------
    expressionN : float
        A numerical value calculated using the appropriate expression for a portion 
        of the functions range.
		
    conditionN : bool
        A boolean that is true if the independent variable(s) are within the portion 
        of the domain whose range is defined by the corresponding exrpession.
		
    otherwise : float
        An optional numerical value that returned if all conditions are false.
		
    Returns
    -------
    float
        The first expression passed as arguement with a true condition, read left to right.
        If all conditions are false, otherwise is returned.

    Raises
    ------
    Exception
        Raises a generic exception if no condition is true and otherwise is not specified.
		
    Notes
    -----
    This function is not intended to be used by a user, but is defined to be used by generated 
    Python models. This function is defined in a way that matches how libSBML formats piecewise
    functions used in SBML models. The inputs to this function are not actually evaluted inside
    of the function, but are evaluated before being passed to the function.
    For example, if Piecewise() is called like so:
		
    .. math:: Piecewise(x + 2, x < 3, x + 4, x > 3)
			
    and if x = 2, then the arguments will be evaluated to
		
    .. math:: Piecewise(4, True, 6, False)
			
    In this case, if x = 3, an exception will be raised, as all conditons will evaluate to False and
    otherwise is not provided.
    """
    for i in range(int(len(args)/2)):
        if args[i*2 + 1] == True:
            return args[i*2]
    
    if len(args) % 2 == 1:
        return args[-1]
    else:
        raise Exception('Piecewise called with dependant variable outside of range or Otherwise not provided.')

#def SearchReactionByName(rxnDict, name, suppress = False):
#    returnList = []
#    for key, metadata in rxnDict.items():
#        if metadata.name == name:
#            returnList.append(key)
#    
#    if returnList == [] and not suppress:
#        raise Exception('No matching name found')
#    
#    if len(returnList) == 1:
#        return returnList[0]
#    else:
#        return returnList