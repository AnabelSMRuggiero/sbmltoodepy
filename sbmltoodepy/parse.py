# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 12:12:05 2018

@author: Steve
"""
from libsbml import *
import re
import sys
from sbmltoodepy.dataclasses import *

def ParseParameterAssignment(parameterIndex, parameter):
    #if parameter.isSetValue() and parameter.isSetId():
    
#    if parameter.isSetName():
#        parameterName = parameter.getName()
#        parameterValue = parameter.getValue()
#        #print("Parameter," + str(parameterIndex + 1) + "\n" + str(parameterId) + "\nValue\n" + str(parameterValue))
#        outputFile.write(str(parameterName) + ";" + str(parameterValue) + "\n")
#        return
    newParameter = ParameterData()
    if parameter.isSetName():
        newParameter.name = parameter.getName()
    else:
        newParameter.name = ''
        
    if parameter.isSetId():
        newParameter.Id = parameter.getId()
        newParameter.value = parameter.getValue()
        if parameter.isSetConstant():
            newParameter.isConstant = parameter.getConstant()
        else:
            newParameter.isConstant = False
        #print("Parameter," + str(parameterIndex + 1) + "\n" + str(parameterId) + "\nValue\n" + str(parameterValue))
#        outputFile.write(str(parameterId) + "; " + str(parameterValue) + "; " + str(parameterConst) + "; " + parameterName + "\n")
    else:
        raise Exception('Parameter with no Id')
        

    return newParameter
        
def ParseRule(ruleIndex, rule):
    
    if rule.isAlgebraic():
        raise Exception("Algebraic rules are currently not supported")

#        outputFile.write("Rule; " + ruleId + "; Algebraic; " + ruleName + "\n")
    if rule.isAssignment():
        newRule = AssignmentRuleData()
#        outputFile.write("Rule; " + ruleId + "; Assignment; " + ruleName + "\n")
    if rule.isRate():
        newRule = RateRuleData()
#        outputFile.write("Rule; " + ruleId + "; Rate; " + ruleName + "\n")
    
    if rule.isSetName():
        newRule.name = str(rule.getName())
    else:
        newRule.name = ''
        
    if rule.isSetIdAttribute():
        newRule.Id = str(rule.getIdAttribute())
    else:
        newRule.Id = str(ruleIndex + 1)
    
    newRule.variable = rule.getVariable()
    newRule.math = formulaToL3String(rule.getMath())


        

    if rule.getVariable() == "":
        raise Exception("Algebraic rules are currently not supported")
        
    else:
        return newRule
#        outputFile.write(rule.getVariable() + ";" + formulaToL3String(rule.getMath()) + "\n")
        
def ParseSpecies(speciesIndex, species):

    newSpecies = SpeciesData()
    
    if species.isSetName():
        newSpecies.name = species.getName()
    else:
        newSpecies.name = ''
  
#    outputFile.write("Species; " + str(speciesIndex + 1) + "; " + speciesName + "\n" )
	

    newSpecies.Id = species.getId()
    newSpecies.compartment = species.getCompartment()
    newSpecies.isConstant = species.getConstant()
    newSpecies.isBoundarySpecies = species.getBoundaryCondition()
    newSpecies.hasOnlySubstanceUnits = species.getHasOnlySubstanceUnits()
    if species.isSetInitialAmount():
        newSpecies.valueType = 'Amount'
        newSpecies.value = species.getInitialAmount()
           #assert(species.isSetInitialAmount())
  
    else:
        newSpecies.valueType = 'Concentration'
        newSpecies.value = species.getInitialConcentration()
		#assert(species.isSetInitialConcentration())		
  
    return newSpecies
		
#outputFile.write(species.getId() + "; " + valueType + "; " + speciesCompartment + '; ' + str(value) + '; ' + str(speciesConstant) + '; ' + str(speciesBoundary) + '; ' + str(speciesFormulaUnits) + '\n')
    
		
    
        
def ParseReaction(reactionIndex, reaction):
    newReaction = ReactionData()
    if reaction.isSetIdAttribute():
        newReaction.Id = reaction.getIdAttribute()
    else:
        newReaction.Id = str(ruleIndex + 1)
        
    if reaction.isSetName():
        newReaction.name = reaction.getName()
    else:
        newReaction.name = ''
		
#    outputFile.write("Reaction; " + reactionId + "; " + reactionName + "\n")     
    
    numReactants = reaction.getListOfReactants().size()
    numProducts = reaction.getListOfProducts().size()
    
    newReaction.reactants = []
    i=0
    for i in range(numReactants):
        reactantStoich = -float(reaction.getListOfReactants().get(i).getStoichiometry())        
        reactantSpecies = reaction.getListOfReactants().get(i).getSpecies()

        newReaction.reactants.append([reactantStoich, reactantSpecies])

        
    i=0
    for i in range(numProducts):
        productStoich = float(reaction.getListOfProducts().get(i).getStoichiometry())
#        productsString += " , "         
        productSpecies = reaction.getListOfProducts().get(i).getSpecies()
        
        newReaction.reactants.append([productStoich, productSpecies])
        
        
    
    rateLawObject = reaction.getKineticLaw()
    
    if rateLawObject.getMath() != None:
        newReaction.rateLaw = formulaToL3String(rateLawObject.getMath())
#        outputFile.write(formulaToString(rateLaw.getMath()) + "\n")
    else:
        raise Exception("Rate law defined by plugin that is not currently supported")
#        outputFile.write("None\n")
        
    numRateLawParams = rateLawObject.getNumParameters()
    
    if numRateLawParams > 0:
        newReaction.rxnParameters = []
        for i in range(numRateLawParams):
            param = rateLawObject.getParameter(i)
#            paramsString += param.getId() + "; " + str(param.getValue())
            newReaction.rxnParameters.append([param.getId(), param.getValue()])
                
    else:
        newReaction.rxnParameters = []
        
    return newReaction
#    outputFile.write(paramsString + "\n")
        
def ParseCompartment(compartmentIndex, compartment):
    
    newCompartment = CompartmentData()
    
    newCompartment.Id = compartment.getId()
    
    if compartment.isSetName():
        newCompartment.name = compartment.getName()
    else:
        newCompartment.name = ''
        
    if compartment.isSetSize():
        newCompartment.size =  compartment.getSize()
    else:
        newCompartment.size = None
        
    if compartment.isSetSpatialDimensions():
        newCompartment.dimensionality = compartment.getSpatialDimensions()
    else:
        newCompartment.dimensionality = None
        
    if compartment.isSetConstant():
        newCompartment.isConstant = compartment.getConstant()
    else:
        newCompartment.isConstant = False
        
#    outputFile.write("Compartment; " + str(compartmentIndex + 1) + "\nSize; " + str(size) + "\nDimensionality; " + str(dimensions) + "\n")
#    outputFile.write("Compartment; " + compartment.getId() + "; " + compartmentName + "\nSize; " + str(size) + "\nDimensionality; " + str(dimensions) + "\n")
    return newCompartment
    
def ParseFunction(functionIndex, function):
    
    newFunction = FunctionData()
#    outputFile.write("Function; " + str(functionIndex + 1) + "\n")
    
    if function.isSetName():
        newFunction.name = function.getName()
    else:
        newFunction.name = ''
        
#    outputFile.write(function.getId() + ";" + functionName + "\n")
    newFunction.Id = function.getId()
    newFunction.mathString = formulaToL3String(function.getMath())
    numArguments = function.getNumArguments()
    funcStringIter = re.finditer(",", newFunction.mathString)
    newFunction.arguments = []
    for i in range(numArguments):
        match = next(funcStringIter)
        newFunction.arguments.append(formulaToString(function.getArgument(i)))


    newFunction.mathString = newFunction.mathString[match.end()+1:-1]

#    functionMath = formulaToL3String(function.getMath())
#    numArguments = function.getNumArguments()
#    funcStringIter = re.finditer(",", functionMath)
#    argumentString = ""
#    for i in range(numArguments):
#        match = next(funcStringIter)
#        argumentString += formulaToString(function.getArgument(i))
#        if i != numArguments-1:
#            argumentString += ";"
#    argumentString += "\n"
#    functionMath = functionMath[match.end()+1:-1]
#    outputFile.write(argumentString)
#    outputFile.write(functionMath + "\n")
    
    
    return newFunction
    
def ParseInitialAssignment(assignmentIndex, assignment):
    
    newAssignment = InitialAssignmentData()
    
    if assignment.isSetIdAttribute():
        newAssignment.Id = str(assignment.getIdAttribute())
    else:
        newAssignment.Id = str(assignmentIndex + 1)
        
    newAssignment.variable = assignment.getSymbol()
    newAssignment.math = formulaToL3String(assignment.getMath())
    newAssignment.name = assignment.getName()
    return newAssignment
    
def ParseSBMLFile(filePath):
    """
    Parameters
    ----------
    filePath : string
        File path of the SBML model to be parsed
		

    Returns
    -------
    ModelData
        An object containing the model's components and their properties.

    """
    doc = readSBML(filePath)
    
    assert(doc.getNumErrors() == 0)
    
    model = doc.getModel()    
    
    modelData = ModelData()
#    outputFile = open(outputPath, "w")
	
#    outputFile.write("Parameters"+ "\n")
    for i in range(model.getNumParameters()):
        newParameter = ParseParameterAssignment(i, model.getParameter(i))
        modelData.parameters[newParameter.Id] = newParameter
#    outputFile.write("Compartments" + "\n")
    for i in range(model.getNumCompartments()):
        newCompartment = ParseCompartment(i, model.getCompartment(i))        
        modelData.compartments[newCompartment.Id] = newCompartment
#    outputFile.write("Species"+ "\n")
    for i in range(model.getNumSpecies()):
        newSpecies = ParseSpecies(i, model.getSpecies(i))
        modelData.species[newSpecies.Id] = newSpecies 
#    outputFile.write("Functions" + "\n")
    for i in range(model.getNumFunctionDefinitions()):
        newFunction = ParseFunction(i, model.getFunctionDefinition(i))
        modelData.functions[newFunction.Id] = newFunction
#    outputFile.write("Rules"+ "\n")
    for i in range(model.getNumRules()):
        newRule = ParseRule(i,model.getRule(i))
        if type(newRule) == AssignmentRuleData:
            modelData.assignmentRules[newRule.Id] = newRule
        elif type(newRule) == RateRuleData:
            modelData.rateRules[newRule.Id] = newRule
#    outputFile.write("Reactions"+ "\n")
    for i in range(model.getNumReactions()):
        newReaction = ParseReaction(i, model.getReaction(i))     
        modelData.reactions[newReaction.Id] = newReaction

    for i in range(model.getNumInitialAssignments()):
        newAssignment = ParseInitialAssignment(i, model.getInitialAssignment(i))
        modelData.initialAssignments[newAssignment.Id] = newAssignment
    #print(model.getNumSpecies())
    return modelData
     
if __name__ == '__main__':
    filePath = sys.argv[1]
#    outputPath = sys.argv[2]
    modelData = ParseSBMLFile(filePath)
    