# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 15:50:45 2018

@author: Steve
"""

import re
import numpy as np
import sys


            
def GenerateModel(modelData, outputFilePath, objectName = 'SBMLmodel'):
    """
    This function takes model data, either from ParseSBMLFIle() or imported from a .json file,
    and generates a Python file containing a class that implements the SBML model.
    
    Parameters
    ----------
    modelData : ModelData
        An object containing all of the model components and values.
		
    outputFilePath : str
        The desired file path of the resulting python file.
		
    objectName : str
        The name of the class defined in the resulting python file.
		
    Returns
    -------
    None
    """
    #The library mathFuncs serves to both only allow functions supported
    #functions in SBML/user defined functions, but also the python equivalent
    
    np.set_printoptions(threshold=sys.maxsize)
    
    
    outputFile = open(outputFilePath, "w")

    parameters = modelData.parameters
    compartments = modelData.compartments
    species = modelData.species
    reactions = modelData.reactions
    functions = modelData.functions
    
    assignmentRules = modelData.assignmentRules
    rateRules = modelData.rateRules
    initialAssignments = modelData.initialAssignments
    
    mathFuncs = {'abs' : 'abs',
                 'max' : 'max',
                 'min' : 'min',
                 'pow' : 'pow',
                 'exp' : 'math.exp',
                 'floor' : 'np.floor',
                 'ceiling' : 'math.ceil',
                 'exp' : 'math.exp',
                 'ln' : 'math.log',
                 'log' : 'math.log10',
                 'factorial' : 'math.factorial',
                 'sqrt' : 'math.sqrt',
                 
                 'eq' : 'operator.eq',
                 'neq' : 'operator.ne',
                 'gt' : 'operator.gt',
                 'lt' : 'operator.lt',
                 'geq' : 'operator.ge',
                 'leq' : 'operator.le',
                 
                 'and' : 'operator.and_',
                 'or' : 'operator.or_',
                 'xor' : 'operator.xor_',
                 'not' : 'operator.not_',
                 
                 'sin' : 'np.sin',
                 'cos' : 'np.cos',
                 'tan' : 'np.tan',
                 'sec' : '1/np.cos',
                 'csc' : '1/np.sin',
                 'cot' : '1/np.tan',
                 'sinh' : 'np.sinh',
                 'cosh' : 'np.cosh',
                 'tanh' : 'np.tanh',
                 'sech' : '1/np.cosh',
                 'csch' : '1/np.sinh',
                 'coth' : '1/np.tanh',
                 'arcsin' : 'np.arcsin',
                 'arccos' : 'np.arccos',
                 'arctan' : 'np.arctan',
                 'arcsinh' : 'np.arcsinh',
                 'arccosh' : 'np.arccosh',
                 'arctanh' : 'np.arctanh',
                 
                 'true' : 'True',
                 'false' : 'False',
                 'notanumber' : 'np.nan',
                 'pi' : 'np.pi',
                 'infinity' : 'np.inf',
                 'exponentiale' : 'np.e',
                 'piecewise' : 'sbmltoodepy.modelclasses.Piecewise'
                 } 
    #Add in user defined functions
#    for function in functions:
#        mathFuncs[function] = "self." + function
		
    #Set up stoichCoeffMat, a matrix of stoichiometric coefficients for solving the reactions
    reactantCounter = 0
    reactantIndex = {}
    reactionCounter = 0
    reactionIndex = {}
    
    rateRuleVars = []
    rateParams = 0
    for specie in species:
        reactantIndex[specie] = reactantCounter
        reactantCounter += 1
    for key, rateRule in rateRules.items():
        if rateRule.variable in parameters or rateRule.variable in compartments:
            rateParams += 1
            reactantIndex[rateRule.variable] = reactantCounter
            reactantCounter += 1
            rateRuleVars.append(rateRule.variable)
        elif rateRule.variable in species:
            pass
        else:
            raise Exception("Rate Rule adjusting something other than specie amount, parameter value, or compartment size.")

    		
    stoichCoeffMat = np.zeros([len(species) + rateParams, max(len(reactions),1)])
    
    for rxnId in reactions:
        reactionIndex[rxnId] = reactionCounter
        reactionCounter += 1
        reaction = reactions[rxnId]
        for reactant in reaction.reactants:
            if reactant[1] not in reactantIndex:
                reactantIndex[reactant[1]] = reactantCounter
                reactantCounter += 1
            if not (species[reactant[1]].isBoundarySpecies == "True"):
                stoichCoeffMat[reactantIndex[reactant[1]], reactionIndex[rxnId]] += reactant[0]

    	
    # for reaction in reactions:
        # for reactant in reactions[reaction][0]:
            # if reactant[1] not in reactantIndex:
                # reactantIndex[reactant[1]] = reactantCounter
                # reactantCounter += 1
            # if not species[reactant[1]][4]:
                # stoichCoeffMat[reactantIndex[reactant[1]], reaction-1] += reactant[0]
    #print(rateParams)
    #print(stoichCoeffMat)
                
    outputFile.write("import sbmltoodepy.modelclasses\n")
    outputFile.write("from scipy.integrate import solve_ivp\n")
    outputFile.write("import numpy as np\n")
    outputFile.write("import operator\n")
    outputFile.write("import math\n\n")
    
    outputFile.write("class " + objectName +"(sbmltoodepy.modelclasses.Model):\n\n")
    
    outputFile.write("\tdef __init__(self):\n\n")
    outputFile.write("\t\tself.p = {} #Dictionary of model parameters\n")
    for paramId in parameters:
        outputFile.write("\t\tself.p[\'" + paramId + "\'] = sbmltoodepy.modelclasses.Parameter(" + 
                         str(parameters[paramId].value)+ ", \'"+ paramId + "\', " + 
                         str(parameters[paramId].isConstant) + ', metadata = sbmltoodepy.modelclasses.SBMLMetadata("' + 
                         parameters[paramId].name +'"))\n')
        
    outputFile.write("\n\t\tself.c = {} #Dictionary of compartments\n")
    for compartmentId in compartments:
        outputFile.write("\t\tself.c[\'" + compartmentId + "\'] = sbmltoodepy.modelclasses.Compartment(" + 
                         str(compartments[compartmentId].size) + ", " + str(compartments[compartmentId].dimensionality) + 
                         ", " + str(compartments[compartmentId].isConstant) + 
                         ', metadata = sbmltoodepy.modelclasses.SBMLMetadata("' + compartments[compartmentId].name +'"))\n')
        
    outputFile.write("\n\t\tself.s = {} #Dictionary of chemical species\n")
    for speciesId in species:
#        outputFile.write("\t\tspeciesMetadata = SBMLMetadata('" + species[speciesId].name +"')\n")
        outputFile.write("\t\tself.s[\'" + speciesId + "\'] = sbmltoodepy.modelclasses.Species(" + 
                         str(species[speciesId].value) + ", '" + species[speciesId].valueType + 
                         "', self.c['" + species[speciesId].compartment + "'], " + 
                         str(species[speciesId].hasOnlySubstanceUnits) + ", constant = " + 
                         str(species[speciesId].isConstant) + ', metadata = sbmltoodepy.modelclasses.SBMLMetadata("' + 
                         species[speciesId].name +'"))\n')
        for key, rule in assignmentRules.items():
            if rule.variable == speciesId:
                outputFile.write("\t\tself.s[\'" + speciesId + "\']._modifiedBy = " + rule.Id + "\n")
        for key, rule in rateRules.items():
            if rule.variable == speciesId:
                outputFile.write("\t\tself.s[\'" + speciesId + "\']._modifiedBy = " + rule.Id + "\n")
    
                
    outputFile.write("\n\t\tself.r = {} #Dictionary of reactions\n")
    for reactionId in reactions:
        outputFile.write("\t\tself.r[\'" + reactionId + "\'] = " + reactionId + "(self)\n")
    
    
    outputFile.write("\n\t\tself.f = {} #Dictionary of function definitions\n")
    for functionId in functions:
        outputFile.write("\t\tself.f[\'" + functionId + "\'] = " + functionId + "(self)\n")
    
    outputFile.write("\t\tself.time = 0\n\n")
    
#    outputFile.write("\t\tself.reactionMetadata = {")
#    commaFlag = 0
#    for reactionId in reactions:
#        if commaFlag == 0:
#            commaFlag = 1
#            outputFile.write("\n\t\t")
#        else:
#            outputFile.write(",\n\t\t")
#        outputFile.write("self.Reaction" + reactionId + ": SBMLMetadata('" + reactions[reactionId].name + "')")
#    outputFile.write("\n\t\t}\n")
#        
    outputFile.write('\t\tself.AssignmentRules()\n\n')
    

    #These functions are defined here due to reading variables in the parent function's namespace
    #These are not intended to be used elsewhere
    def ParseLHS(rawLHS):
        returnLHS = ''
        if rawLHS in parameters:
            returnLHS = "self.p[\'" + rawLHS + "\'].value = "
        elif rawLHS in species:
            if not species[rawLHS].hasOnlySubstanceUnits: 
                returnLHS = 'self.s[\'' + rawLHS + '\'].concentration = '
            else: 
                returnLHS = 'self.s[\'' + rawLHS + '\'].amount = '
        elif rawLHS in compartments:
            returnLHS = 'self.c[\'' + rawLHS + '\'].size = '
        else:
            raise(Exception("New case: rule LHS not in p: " + rawLHS))

        return returnLHS
	
    def ParseRHS(rawRHS, extendedParams = [], objectText = "self"):
        #objectText is not "self" when parsing reaction math
        
        #The main purpose of this function is to turn math strings given by libSBML into
        #code formated to properly call members of the resulting class
        #For example k_1*C_A may turn to
        
        
        rawRHS = rawRHS.replace("^", "**") #Replaces carrot notation for exponentiation with ** operator
        variables = []
        for match in re.finditer(r'\b[a-zA-Z_]\w*', rawRHS): #look for variable names
            #ToDo: check for function calls
            variables.append([rawRHS[match.start():match.end()], match.span()])
            
        returnRHS = ''
        oldSpan = None
        if variables != []:
            for variable in variables:
                if oldSpan == None and variable[1][0] != 0:
                    returnRHS += rawRHS[0:variable[1][0]]
                elif oldSpan != None:
                    returnRHS += rawRHS[oldSpan[1]:variable[1][0]]
                oldSpan = variable[1]
                if variable[0] in extendedParams:
                    if objectText == "self":
                        returnRHS += variable[0]
                    else:
                        returnRHS += "self.p[\'" + variable[0] + "\'].value"
                elif variable[0] in parameters:
                    returnRHS += objectText + '.p[\'' + variable[0] + '\'].value'
                elif variable[0] in species:
                    if not species[variable[0]].hasOnlySubstanceUnits == "True": 
                        returnRHS += objectText + '.s[\'' + variable[0] + '\'].concentration'
                    else: 
                        returnRHS += objectText + '.s[\'' + variable[0] + '\'].amount'
                elif variable[0] in compartments:
                    returnRHS += objectText + '.c[\'' + variable[0] + '\'].size'
                elif variable[0] in mathFuncs:
                    returnRHS += mathFuncs[variable[0]]
                elif variable[0] in functions:
                    returnRHS += objectText + '.f[\'' + variable[0] + '\']'
                

                elif variable[0] == "time":
                    returnRHS += objectText + '.time'
                elif variable[0] == "pi":
                    returnRHS += "np.pi"
                else:
                    raise(Exception('New case: unkown RHS variable: ' + variable[0]))
            returnRHS += rawRHS[variable[1][1]:len(rawRHS)]
    #        print(rule[1][variable[1][1]])
            #print(rule[1][-1])
        else:
            returnRHS = rawRHS
		
        return returnRHS
    
    outputFile.write("\n\n")
    outputFile.write("\tdef AssignmentRules(self):\n\n")
    #Here's where the assignment rules and initial assignments are printed. This may be changes to a class based implementation later.

    ruleDefinedVars = [rule.variable for rule in assignmentRules.values()]
    for key, assignment in initialAssignments.items():
        ruleDefinedVars.append(assignment.variable)
        
    for key, rule in assignmentRules.items():
        rule.dependents = []
        for match in re.finditer(r'\b[a-zA-Z_]\w*', rule.math): #look for variable names
            rule.dependents.append(rule.math[match.start():match.end()])
        originalLen = len(rule.dependents)
        for i in range(originalLen):
            if rule.dependents[originalLen - i -1] not in ruleDefinedVars:
                rule.dependents.pop(originalLen- i-1)
                
    for key, assignment in initialAssignments.items():
        assignment.dependents = []
        for match in re.finditer(r'\b[a-zA-Z_]\w*', assignment.math): #look for variable names
            assignment.dependents.append(assignment.math[match.start():match.end()])
        originalLen = len(assignment.dependents)
        for i in range(originalLen):
            if assignment.dependents[originalLen - i -1] not in ruleDefinedVars :
                assignment.dependents.pop(originalLen- i-1)
                
#    breakVar = False
    while True:
        continueVar = False
        breakVar = True
        varDefinedThisLoop = None
        for key, rule in assignmentRules.items():
            if rule.dependents == []:
                ruleLHS = ParseLHS(rule.variable)
                ruleRHS = ParseRHS(rule.math)
                outputFile.write("\t\t" + ruleLHS + ruleRHS + '\n\n')
                varDefinedThisLoop = rule.variable
                rule.dependents = None
                continueVar = True
                breakVar = False
                break
            elif not rule.dependents == None:
                breakVar = False
                
        if not continueVar:
            for key, assignment in initialAssignments.items():
                if assignment.dependents == []:
                    assignmentLHS = ParseLHS(assignment.variable)
                    assignmentRHS = ParseRHS(assignment.math)
                    outputFile.write("\t\tif self.time <= 0 :\n")
                    if assignment.variable in parameters:
                        outputFile.write("\t\t\tisConstantValue = self.p['" + assignment.variable + "']._constant\n")
                        outputFile.write("\t\t\tself.p['" + assignment.variable + "']._constant = False\n")
                        outputFile.write("\t\t\t" + assignmentLHS + assignmentRHS + '\n')
                        outputFile.write("\t\t\tself.p['" + assignment.variable + "']._constant = isConstantValue\n\n")
                    elif assignment.variable in species:
                        outputFile.write("\t\t\tisConstantValue = self.s['" + assignment.variable + "']._constant\n")
                        outputFile.write("\t\t\tself.s['" + assignment.variable + "']._constant = False\n")
                        outputFile.write("\t\t\t" + assignmentLHS + assignmentRHS + '\n')
                        outputFile.write("\t\t\tself.s['" + assignment.variable + "']._constant = isConstantValue\n\n")
                    elif assignment.variable in compartments:
                        outputFile.write("\t\t\tisConstantValue = self.c['" + assignment.variable + "']._constant\n")
                        outputFile.write("\t\t\tself.c['" + assignment.variable + "']._constant = False\n")
                        outputFile.write("\t\t\t" + assignmentLHS + assignmentRHS + '\n')
                        outputFile.write("\t\t\tself.c['" + assignment.variable + "']._constant = isConstantValue\n\n")
 
                    varDefinedThisLoop = assignment.variable
                    assignment.dependents = None
                    continueVar = True
                    breakVar = False
                    break
                elif not assignment.dependents == None:
                    breakVar = False
            
        for rule in assignmentRules.values():
            if not rule.dependents == None:
                originalLen = len(rule.dependents)
                for i in range(originalLen):
                    if rule.dependents[originalLen - i -1] == varDefinedThisLoop:
                        rule.dependents.pop(originalLen - i -1)
#            print(rule.variable + ':' + str(rule.dependents))

        for assignment in initialAssignments.values():
            if not assignment.dependents == None:
                originalLen = len(assignment.dependents)
                for i in range(originalLen):
                    if assignment.dependents[originalLen - i - 1] == varDefinedThisLoop:
                        assignment.dependents.pop(originalLen - i - 1)
#            print(assignment.variable + ':' + str(assignment.dependents))
        
        if continueVar:
            continue
        elif breakVar:
            break
        else:
            raise Exception('Algebraic Loop in AssignmentRules')
        
    outputFile.write("\t\treturn\n\n")
        
#    for functionId in functions:
#        arguments = functions[functionId].arguments
#        argumentString = ""
#        for i in range(len(arguments)):
#            argumentString += arguments[i]
#            if i != len(arguments) - 1:
#                argumentString += ", "
#        
#        outputFile.write("\tdef " + functionId + "(self, " + argumentString + "):\n")
#        outputFile.write("\t\treturn " + functions[functionId].mathString.replace("^", "**") + "\n")
        
#    for reactionId in reactions:
#        outputFile.write("\tdef Reaction" + str(reactionId) + "(self):\n\n")
#
#        rxnParameters = []
#        for param in reactions[reactionId].rxnParameters:
#            outputFile.write("\t\t" + param[0] + " = " + str(param[1]) + "\n")
#            rxnParameters.append(param[0])
#			
#        rateLaw = ParseRHS(reactions[reactionId].rateLaw, rxnParameters)
#           
#        outputFile.write('\t\treturn ' + rateLaw + '\n\n')

    rateRuleLHSVars = []
    for key, rateRule in rateRules.items():
        rateRuleLHSVars.append(rateRule.variable)
        outputFile.write("\tdef Rate" + rateRule.variable + "(self):\n\n")
        rateLaw = ParseRHS(rateRule.math)
        outputFile.write('\t\treturn ' + rateLaw + '\n\n')
        
    yArray = ''
    i = 0
    yArrayVars = [0 for x in range(len(species) + rateParams)]
    for variable, index in reactantIndex.items():
        yArrayVars[index] = variable
    
    for index in range(len(yArrayVars)):
        # print(yArrayVars[index])
        if index != 0:
            yArray += ', '
        
        if yArrayVars[index] in species:
            yArray += 'self.s[\'' + yArrayVars[index] + '\'].amount'
            continue
    
        if yArrayVars[index] in parameters:
            yArray += 'self.p[\'' + yArrayVars[index] + '\'].value'
            continue
        
        if yArrayVars[index] in compartments:
            yArray += 'self.c\'' + yArrayVars[index] + '\'].size'
            continue
    

            
    outputFile.write('\tdef _SolveReactions(self, t, y):\n\n')
    outputFile.write('\t\tself.time = t\n')
    outputFile.write('\t\t' + yArray + ' = y\n')
    outputFile.write('\t\tself.AssignmentRules()\n\n')
#    outputFile.write('\t\t[self.s[speciesId].UpdateCompartmentSizeMember() for speciesId in self.s]\n')
    rateArray = '[ '
    i = 0
    rateArrayVars = [0 for x in range(len(species) + rateParams)]
    
    for variable, index in reactantIndex.items():
        if variable in rateRuleLHSVars:
            rateArrayVars[index] = variable
    

    
    for variable in rateArrayVars:
        if i != 0:
            rateArray += ', '
        i += 1
        if variable == 0:
            rateArray += '0'
        else:
            rateArray += 'self.Rate' + variable + '()'
            
            
        
        
    rateArray += ']'
    outputFile.write('\t\trateRuleVector = np.array(' + str(rateArray) + ', dtype = np.float64)\n\n') 
    
    outputFile.write('\t\tstoichiometricMatrix = np.array(' + re.sub('\n,', ',\n\t\t\t\t\t', re.sub('[^[] +', ',' ,str(stoichCoeffMat))) + ', dtype = np.float64)\n\n')
    outputFile.write('\t\treactionVelocities = np.array([')
    reactionElements = ''
    if reactions:
        for reactionId in reactions:
            if reactionElements == '':
                reactionElements += ('self.r[\'' + str(reactionId) + '\']()')
            else:
                reactionElements += (', self.r[\'' + str(reactionId) + '\']()')
    else:
        reactionElements = '0'
    outputFile.write(reactionElements + '], dtype = np.float64)\n\n')
    outputFile.write('\t\trateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector\n\n')
    outputFile.write('\t\treturn rateOfSpeciesChange\n\n')
    
    outputFile.write('\tdef RunSimulation(self, deltaT, method="LSODA", args=None, absoluteTolerance=1e-12, relativeTolerance=1e-6, **options):\n\n')
    
    outputFile.write('\t\tfinalTime = self.time + deltaT\n')
    outputFile.write('\t\ty0 = np.array([' + yArray + '], dtype = np.float64)\n')
    outputFile.write('\t\t' + yArray + ' = solve_ivp(self._SolveReactions, [self.time, finalTime], y0, method=method, t_eval=[self.time, finalTime], args=args, atol=absoluteTolerance, rtol=relativeTolerance, **options).y[:, -1]\n')
    outputFile.write('\t\tself.time = finalTime\n')
    outputFile.write('\t\tself.AssignmentRules()\n')
#    outputFile.write('\t\t[self.s[speciesId].UpdateCompartmentSizeMember() for speciesId in self.s]\n')
    outputFile.write('\n')
    
    for key in reactions.keys():
        outputFile.write('class ' + key + ':\n\n')
        outputFile.write('\tdef __init__(self, parent, metadata = None):\n\n')
        outputFile.write('\t\tself.parent = parent\n')
        outputFile.write('\t\tself.p = {}\n')
        outputFile.write('\t\tif metadata:\n')
        outputFile.write('\t\t\tself.metadata = metadata\n')
        outputFile.write('\t\telse:\n')
        outputFile.write('\t\t\tself.metadata = sbmltoodepy.modelclasses.SBMLMetadata("' + reactions[key].name + '")\n')
        
        for param in reactions[key].rxnParameters:
            outputFile.write("\t\tself.p[\'" + param[0] + "\'] = sbmltoodepy.modelclasses.Parameter(" + str(param[1]) + ", '" + param[0] + "')\n")
            #"\t\tself.p[\'" + paramId + "\'] = Parameter(" + str(parameters[paramId].value)+ ", "+ paramId + ", " + str(parameters[paramId].isConstant) +")\n"
        
        outputFile.write('\n\tdef __call__(self):\n')
#        print(key)
#        print(reactions[key].rxnParameters)
        rxnParamNames = [param[0] for param in reactions[key].rxnParameters]
        rateLaw = ParseRHS(reactions[key].rateLaw, rxnParamNames, "self.parent")
        outputFile.write('\t\treturn ' + rateLaw + '\n\n')

    def ParseFunctionMath(rawRHS):
        #objectText is not "self" when parsing reaction math
        
        #The main purpose of this function is to turn math strings given by libSBML into
        #code formated to properly call members of the resulting class
        #For example k_1*C_A may turn to
        
        
        rawRHS = rawRHS.replace("^", "**") #Replaces carrot notation for exponentiation with ** operator
        variables = []
        for match in re.finditer(r'\b[a-zA-Z_]\w*', rawRHS): #look for variable names
            #ToDo: check for function calls
            variables.append([rawRHS[match.start():match.end()], match.span()])
            
        returnRHS = ''
        oldSpan = None
        if variables != []:
            for variable in variables:
                if oldSpan == None and variable[1][0] != 0:
                    returnRHS += rawRHS[0:variable[1][0]]
                elif oldSpan != None:
                    returnRHS += rawRHS[oldSpan[1]:variable[1][0]]
                oldSpan = variable[1]
                
                if variable[0] in mathFuncs:
                    returnRHS += mathFuncs[variable[0]]
                else:
                    returnRHS += variable[0]
            returnRHS += rawRHS[variable[1][1]:len(rawRHS)]
    #        print(rule[1][variable[1][1]])
            #print(rule[1][-1])
        else:
            returnRHS = rawRHS
		
        return returnRHS
    
    for key in functions.keys():
        outputFile.write('class ' + key + ':\n\n')
        outputFile.write('\tdef __init__(self, parent, metadata = None):\n\n')
        outputFile.write('\t\tself.parent = parent\n')
        outputFile.write('\t\tif metadata:\n')
        outputFile.write('\t\t\tself.metadata = metadata\n')
        outputFile.write('\t\telse:\n')
        outputFile.write('\t\t\tself.metadata = sbmltoodepy.modelclasses.SBMLMetadata("' + functions[key].name + '")\n')

        arguments = functions[key].arguments
        argumentString = ""
        for i in range(len(arguments)):
            argumentString += arguments[i]
            if i != len(arguments) - 1:
                argumentString += ", "
                
        outputFile.write('\tdef __call__(self, ' + argumentString + '):\n')
        mathString = functions[key].mathString.replace("^","**")
        mathString = ParseFunctionMath(mathString)
        outputFile.write("\t\treturn " + mathString + "\n\n")

    outputFile.close()		
		
#GenerateModel("Waugh2006_Diabetic_Wound_Healing_TGF_B_Dynamics.txt")