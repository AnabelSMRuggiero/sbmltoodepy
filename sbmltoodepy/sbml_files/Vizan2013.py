from sbmltopyode.SBMLModelClasses import *
from scipy.integrate import odeint
import numpy as np
import operator
import math

class Vizan2013:

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['parameter_1'] = Parameter(0.32, 'parameter_1', True)
		self.p['parameter_2'] = Parameter(20.0, 'parameter_2', True)
		self.p['parameter_3'] = Parameter(9.36, 'parameter_3', True)
		self.p['parameter_4'] = Parameter(0.08, 'parameter_4', True)
		self.p['parameter_5'] = Parameter(0.0, 'parameter_5', True)
		self.p['parameter_6'] = Parameter(21.3715, 'parameter_6', True)
		self.p['parameter_7'] = Parameter(24.0, 'parameter_7', True)
		self.p['parameter_8'] = Parameter(60.0, 'parameter_8', True)
		self.p['parameter_9'] = Parameter(350.877192982456, 'parameter_9', False)
		self.p['parameter_10'] = Parameter(0.171, 'parameter_10', False)
		self.p['parameter_11'] = Parameter(5.7, 'parameter_11', True)
		self.p['parameter_12'] = Parameter(4.0, 'parameter_12', True)
		self.p['parameter_13'] = Parameter(2.27, 'parameter_13', True)
		self.p['parameter_14'] = Parameter(1.0, 'parameter_14', True)
		self.p['parameter_15'] = Parameter(1.0, 'parameter_15', True)
		self.p['parameter_16'] = Parameter(0.05, 'parameter_16', True)
		self.p['parameter_17'] = Parameter(0.196565, 'parameter_17', True)
		self.p['parameter_18'] = Parameter(100.0, 'parameter_18', True)
		self.p['parameter_19'] = Parameter(9.36, 'parameter_19', False)
		self.p['parameter_20'] = Parameter(1.0, 'parameter_20', False)
		self.p['parameter_21'] = Parameter(0.710526315789474, 'parameter_21', False)
		self.p['parameter_22'] = Parameter(24.5383, 'parameter_22', True)
		self.p['parameter_23'] = Parameter(2.0, 'parameter_23', True)
		self.p['parameter_24'] = Parameter(2.0, 'parameter_24', True)
		self.p['parameter_25'] = Parameter(0.35, 'parameter_25', True)
		self.p['parameter_26'] = Parameter(0.0, 'parameter_26', True)
		self.p['parameter_27'] = Parameter(0.0, 'parameter_27', True)
		self.p['parameter_28'] = Parameter(0.0, 'parameter_28', True)
		self.p['Metabolite_9'] = Parameter(0.558933528122717, 'Metabolite_9', True)

		self.c = {} #Dictionary of compartments
		self.c['compartment_1'] = Compartment(1.0, 3, True)

		self.s = {} #Dictionary of chemical species
		speciesMetadata = SBMLMetadata('S22')
		self.s['species_1'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_1']._modifiedBy = 18
		speciesMetadata = SBMLMetadata('S24')
		self.s['species_2'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_2']._modifiedBy = 19
		speciesMetadata = SBMLMetadata('pS2tot')
		self.s['species_3'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_3']._modifiedBy = 20
		speciesMetadata = SBMLMetadata('TGF')
		self.s['species_4'] = Species(4.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_4']._modifiedBy = 21
		speciesMetadata = SBMLMetadata('R')
		self.s['species_5'] = Species(0.87962962962963, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_5']._modifiedBy = 22
		speciesMetadata = SBMLMetadata('S2c')
		self.s['species_6'] = Species(1.19430241051863, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_6']._modifiedBy = 23
		speciesMetadata = SBMLMetadata('Rcom')
		self.s['species_7'] = Species(0.12037037037037, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_7']._modifiedBy = 24
		speciesMetadata = SBMLMetadata('pS2c')
		self.s['species_8'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_8']._modifiedBy = 25
		speciesMetadata = SBMLMetadata('Rcom_S')
		self.s['species_9'] = Species(0.05, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_9']._modifiedBy = 10
		speciesMetadata = SBMLMetadata('S2n')
		self.s['species_10'] = Species(0.558933528122717, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_10']._modifiedBy = 11
		speciesMetadata = SBMLMetadata('S22n')
		self.s['species_11'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_11']._modifiedBy = 12
		speciesMetadata = SBMLMetadata('S4n')
		self.s['species_12'] = Species(1.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_12']._modifiedBy = 13
		speciesMetadata = SBMLMetadata('S22c')
		self.s['species_13'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_13']._modifiedBy = 26
		speciesMetadata = SBMLMetadata('pS2n')
		self.s['species_14'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_14']._modifiedBy = 15
		speciesMetadata = SBMLMetadata('pS2fn')
		self.s['species_15'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_15']._modifiedBy = 16
		speciesMetadata = SBMLMetadata('S24n')
		self.s['species_16'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_16']._modifiedBy = 1
		speciesMetadata = SBMLMetadata('S24c')
		self.s['species_17'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_17']._modifiedBy = 27
		speciesMetadata = SBMLMetadata('S4fc')
		self.s['species_18'] = Species(1.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_18']._modifiedBy = 2
		speciesMetadata = SBMLMetadata('S4c')
		self.s['species_19'] = Species(1.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_19']._modifiedBy = 28
		speciesMetadata = SBMLMetadata('pS2fc')
		self.s['species_20'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_20']._modifiedBy = 7
		speciesMetadata = SBMLMetadata('S4fn')
		self.s['species_21'] = Species(1.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_21']._modifiedBy = 14
		speciesMetadata = SBMLMetadata('SBI')
		self.s['species_22'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		speciesMetadata = SBMLMetadata('Rtot')
		self.s['species_23'] = Species(1.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_23']._modifiedBy = 8
		speciesMetadata = SBMLMetadata('RT')
		self.s['species_24'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_24']._modifiedBy = 29
		speciesMetadata = SBMLMetadata('Rcom_I')
		self.s['species_25'] = Species(0.0703703703703704, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_25']._modifiedBy = 9
		speciesMetadata = SBMLMetadata('Ract')
		self.s['species_26'] = Species(0.0, 'Concentration', self.c['compartment_1'], False, constant = False)
		self.s['species_26']._modifiedBy = 30

		self.r = {} #Dictionary of reactiions
		self.time = 0

		self.reactionMetadata = {
		}
		self.AssignmentRules()



	def AssignmentRules(self):

		self.s['species_16'].concentration = (self.p['parameter_13'].value + 1) * self.s['species_2'].concentration - self.p['parameter_13'].value * self.s['species_17'].concentration

		self.p['parameter_10'].value = 0.171

		self.p['parameter_9'].value = self.p['parameter_8'].value / self.p['parameter_10'].value

		self.p['parameter_19'].value = self.p['parameter_3'].value

		self.p['parameter_21'].value = (self.p['parameter_16'].value / (1 - self.p['parameter_16'].value)) * ((self.p['parameter_4'].value + 1) / self.p['parameter_4'].value)

		self.s['species_20'].concentration = self.s['species_8'].concentration - 2 * self.s['species_13'].concentration - self.s['species_17'].concentration

		self.s['species_11'].concentration = (self.p['parameter_13'].value + 1) * self.s['species_1'].concentration - self.p['parameter_13'].value * self.s['species_13'].concentration

		self.s['species_14'].concentration = (self.p['parameter_13'].value + 1) * self.s['species_3'].concentration - self.p['parameter_13'].value * self.s['species_8'].concentration

		self.s['species_15'].concentration = self.s['species_14'].concentration - 2 * self.s['species_11'].concentration - self.s['species_16'].concentration

		if self.time <= 0 :
			isConstantValue = self.s['species_4']._constant
			self.s['species_4']._constant = False
			self.s['species_4'].concentration = self.p['parameter_23'].value * self.p['parameter_24'].value
			self.s['species_4']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.s['species_5']._constant
			self.s['species_5']._constant = False
			self.s['species_5'].concentration = (1 - self.p['parameter_16'].value) / (self.p['parameter_4'].value + 1)
			self.s['species_5']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.s['species_6']._constant
			self.s['species_6']._constant = False
			self.s['species_6'].concentration = self.p['parameter_14'].value * self.p['parameter_2'].value * (1 + self.p['parameter_13'].value) / (self.p['parameter_3'].value + self.p['parameter_13'].value * self.p['parameter_2'].value)
			self.s['species_6']._constant = isConstantValue

		self.s['species_10'].concentration = (self.p['parameter_13'].value + 1) * (self.p['parameter_14'].value - self.s['species_3'].concentration) - self.p['parameter_13'].value * self.s['species_6'].concentration

		if self.time <= 0 :
			isConstantValue = self.s['species_7']._constant
			self.s['species_7']._constant = False
			self.s['species_7'].concentration = (self.p['parameter_16'].value + self.p['parameter_4'].value) / (1 + self.p['parameter_4'].value)
			self.s['species_7']._constant = isConstantValue

		self.s['species_23'].concentration = self.s['species_5'].concentration + self.s['species_7'].concentration + self.s['species_24'].concentration + self.s['species_26'].concentration

		self.s['species_25'].concentration = self.s['species_7'].concentration * (1 / (1 + self.p['parameter_21'].value))

		self.s['species_9'].concentration = self.s['species_7'].concentration * (self.p['parameter_21'].value / (1 + self.p['parameter_21'].value))

		if self.time <= 0 :
			isConstantValue = self.s['species_19']._constant
			self.s['species_19']._constant = False
			self.s['species_19'].concentration = self.p['parameter_15'].value
			self.s['species_19']._constant = isConstantValue

		self.s['species_18'].concentration = self.s['species_19'].concentration - self.s['species_17'].concentration

		self.s['species_12'].concentration = (self.p['parameter_13'].value + 1) * self.p['parameter_15'].value - self.p['parameter_13'].value * self.s['species_19'].concentration

		self.s['species_21'].concentration = self.s['species_12'].concentration - self.s['species_16'].concentration

		if self.time <= 0 :
			isConstantValue = self.p['Metabolite_9']._constant
			self.p['Metabolite_9']._constant = False
			self.p['Metabolite_9'].value = self.s['species_10'].concentration
			self.p['Metabolite_9']._constant = isConstantValue

		self.p['parameter_20'].value = (self.s['species_10'].concentration + self.s['species_14'].concentration) / self.p['Metabolite_9'].value

		return

	def Ratespecies_1(self):

		return (1 / (1 + self.p['parameter_13'].value)) * (self.p['parameter_9'].value * (self.p['parameter_13'].value * self.s['species_20'].concentration**2 + self.s['species_15'].concentration**2) - self.p['parameter_8'].value * (self.p['parameter_13'].value * self.s['species_13'].concentration + self.s['species_11'].concentration))

	def Ratespecies_2(self):

		return (1 / (self.p['parameter_13'].value + 1)) * (self.p['parameter_9'].value * (self.p['parameter_13'].value * self.s['species_18'].concentration * self.s['species_20'].concentration + self.s['species_15'].concentration * self.s['species_21'].concentration) - self.p['parameter_8'].value * (self.p['parameter_13'].value * self.s['species_17'].concentration + self.s['species_16'].concentration))

	def Ratespecies_3(self):

		return (1 / (1 + self.p['parameter_13'].value)) * (self.p['parameter_13'].value * self.p['parameter_6'].value * self.s['species_26'].concentration * (self.p['parameter_17'].value / (self.p['parameter_17'].value + self.s['species_22'].concentration)) * self.s['species_6'].concentration - self.p['parameter_7'].value * self.s['species_15'].concentration)

	def Ratespecies_4(self):

		return self.p['parameter_1'].value * (self.p['parameter_27'].value + self.p['parameter_26'].value * self.s['species_16'].concentration - (self.p['parameter_18'].value * self.s['species_9'].concentration + self.p['parameter_25'].value) * self.s['species_4'].concentration)

	def Ratespecies_5(self):

		return self.p['parameter_1'].value * ((1 - self.p['parameter_5'].value) * (1 - self.p['parameter_16'].value) - (self.p['parameter_4'].value + (1 - self.p['parameter_28'].value)) * self.s['species_5'].concentration)

	def Ratespecies_6(self):

		return self.p['parameter_2'].value * self.s['species_10'].concentration - (self.p['parameter_3'].value + self.p['parameter_6'].value * self.s['species_26'].concentration * (self.p['parameter_17'].value / (self.p['parameter_17'].value + self.s['species_22'].concentration))) * self.s['species_6'].concentration

	def Ratespecies_7(self):

		return self.p['parameter_1'].value * (self.p['parameter_4'].value * self.s['species_5'].concentration - (1 - self.p['parameter_28'].value) * self.s['species_25'].concentration - self.p['parameter_18'].value * self.s['species_4'].concentration * self.s['species_9'].concentration)

	def Ratespecies_8(self):

		return self.p['parameter_6'].value * self.s['species_26'].concentration * (self.p['parameter_17'].value / (self.p['parameter_17'].value + self.s['species_22'].concentration)) * self.s['species_6'].concentration + self.p['parameter_2'].value * self.s['species_15'].concentration - self.p['parameter_3'].value * (self.s['species_20'].concentration + self.p['parameter_11'].value * (self.s['species_17'].concentration + 2 * self.s['species_13'].concentration))

	def Ratespecies_13(self):

		return self.p['parameter_9'].value * self.s['species_20'].concentration**2 - (self.p['parameter_8'].value + self.p['parameter_3'].value * self.p['parameter_11'].value) * self.s['species_13'].concentration

	def Ratespecies_17(self):

		return self.p['parameter_9'].value * self.s['species_18'].concentration * self.s['species_20'].concentration - (self.p['parameter_8'].value + self.p['parameter_3'].value * self.p['parameter_11'].value) * self.s['species_17'].concentration

	def Ratespecies_19(self):

		return self.p['parameter_19'].value * self.s['species_21'].concentration - self.p['parameter_3'].value * (self.s['species_18'].concentration + self.p['parameter_11'].value * self.s['species_17'].concentration)

	def Ratespecies_24(self):

		return self.p['parameter_1'].value * (self.p['parameter_18'].value * self.s['species_4'].concentration * self.s['species_9'].concentration - (self.p['parameter_22'].value + self.p['parameter_12'].value * (1 - self.p['parameter_28'].value)) * self.s['species_24'].concentration)

	def Ratespecies_26(self):

		return self.p['parameter_1'].value * (self.p['parameter_22'].value * self.s['species_24'].concentration - self.p['parameter_12'].value * (1 - self.p['parameter_28'].value) * self.s['species_26'].concentration)

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['species_1'].amount, self.s['species_2'].amount, self.s['species_3'].amount, self.s['species_4'].amount, self.s['species_5'].amount, self.s['species_6'].amount, self.s['species_7'].amount, self.s['species_8'].amount, self.s['species_9'].amount, self.s['species_10'].amount, self.s['species_11'].amount, self.s['species_12'].amount, self.s['species_13'].amount, self.s['species_14'].amount, self.s['species_15'].amount, self.s['species_16'].amount, self.s['species_17'].amount, self.s['species_18'].amount, self.s['species_19'].amount, self.s['species_20'].amount, self.s['species_21'].amount, self.s['species_22'].amount, self.s['species_23'].amount, self.s['species_24'].amount, self.s['species_25'].amount, self.s['species_26'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ self.Ratespecies_1(), self.Ratespecies_2(), self.Ratespecies_3(), self.Ratespecies_4(), self.Ratespecies_5(), self.Ratespecies_6(), self.Ratespecies_7(), self.Ratespecies_8(), 0, 0, 0, 0, self.Ratespecies_13(), 0, 0, 0, self.Ratespecies_17(), 0, self.Ratespecies_19(), 0, 0, 0, 0, self.Ratespecies_24(), 0, self.Ratespecies_26()], dtype = np.float64)

		stoichiometricMatrix = np.array([[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.]], dtype = np.float64)

		reactionVelocities = np.array([0], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['species_1'].amount, self.s['species_2'].amount, self.s['species_3'].amount, self.s['species_4'].amount, self.s['species_5'].amount, self.s['species_6'].amount, self.s['species_7'].amount, self.s['species_8'].amount, self.s['species_9'].amount, self.s['species_10'].amount, self.s['species_11'].amount, self.s['species_12'].amount, self.s['species_13'].amount, self.s['species_14'].amount, self.s['species_15'].amount, self.s['species_16'].amount, self.s['species_17'].amount, self.s['species_18'].amount, self.s['species_19'].amount, self.s['species_20'].amount, self.s['species_21'].amount, self.s['species_22'].amount, self.s['species_23'].amount, self.s['species_24'].amount, self.s['species_25'].amount, self.s['species_26'].amount], dtype = np.float64)
		self.s['species_1'].amount, self.s['species_2'].amount, self.s['species_3'].amount, self.s['species_4'].amount, self.s['species_5'].amount, self.s['species_6'].amount, self.s['species_7'].amount, self.s['species_8'].amount, self.s['species_9'].amount, self.s['species_10'].amount, self.s['species_11'].amount, self.s['species_12'].amount, self.s['species_13'].amount, self.s['species_14'].amount, self.s['species_15'].amount, self.s['species_16'].amount, self.s['species_17'].amount, self.s['species_18'].amount, self.s['species_19'].amount, self.s['species_20'].amount, self.s['species_21'].amount, self.s['species_22'].amount, self.s['species_23'].amount, self.s['species_24'].amount, self.s['species_25'].amount, self.s['species_26'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

