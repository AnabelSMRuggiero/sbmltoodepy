from sbmltopyode.SBMLModelClasses import *
from scipy.integrate import odeint
import numpy as np
import operator
import math

class Waugh2006:

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['tau1'] = Parameter(-2.47, 'tau1', True)
		self.p['tau2'] = Parameter(21.94, 'tau2', True)
		self.p['tau3'] = Parameter(6.41, 'tau3', True)
		self.p['tau4'] = Parameter(1.75, 'tau4', True)
		self.p['alpha'] = Parameter(0.8, 'alpha', True)
		self.p['k1'] = Parameter(0.05, 'k1', True)
		self.p['k2'] = Parameter(0.693, 'k2', True)
		self.p['k3'] = Parameter(0.002, 'k3', True)
		self.p['k4'] = Parameter(0.07, 'k4', True)
		self.p['d1'] = Parameter(0.2, 'd1', True)
		self.p['d2'] = Parameter(9.1, 'd2', True)

		self.c = {} #Dictionary of compartments
		self.c['COMpartment'] = Compartment(1.0, 3, True)

		self.s = {} #Dictionary of chemical species
		speciesMetadata = SBMLMetadata('K_T')
		self.s['K_T'] = Species(296.53, 'Concentration', self.c['COMpartment'], False, constant = False)
		self.s['K_T']._modifiedBy = 1
		speciesMetadata = SBMLMetadata('phi_I')
		self.s['phi_I'] = Species(200.0, 'Concentration', self.c['COMpartment'], False, constant = False)
		self.s['phi_I']._modifiedBy = 2
		speciesMetadata = SBMLMetadata('phi_R')
		self.s['phi_R'] = Species(200.0, 'Concentration', self.c['COMpartment'], False, constant = False)
		self.s['phi_R']._modifiedBy = 3
		speciesMetadata = SBMLMetadata('T')
		self.s['T'] = Species(6.0, 'Concentration', self.c['COMpartment'], False, constant = False)
		self.s['T']._modifiedBy = 4

		self.r = {} #Dictionary of reactiions
		self.time = 0

		self.reactionMetadata = {
		}
		self.AssignmentRules()



	def AssignmentRules(self):

		self.s['K_T'].concentration = self.p['tau1'].value * self.s['T'].concentration**3 + self.p['tau2'].value * self.s['T'].concentration**2 + self.p['tau3'].value * self.s['T'].concentration + self.p['tau4'].value

		return

	def Ratephi_I(self):

		return self.p['alpha'].value * self.s['K_T'].concentration + self.p['k1'].value * self.p['k2'].value * self.s['phi_I'].concentration * (1 - self.p['k3'].value * (self.s['phi_I'].concentration + self.s['phi_R'].concentration)) - self.p['d1'].value * self.s['phi_I'].concentration

	def Ratephi_R(self):

		return (1 - self.p['alpha'].value) * self.s['K_T'].concentration + self.p['k1'].value * self.p['k2'].value * self.s['phi_R'].concentration * (1 - self.p['k3'].value * (self.s['phi_I'].concentration + self.s['phi_R'].concentration)) - self.p['d1'].value * self.s['phi_R'].concentration

	def RateT(self):

		return self.p['k4'].value * self.s['phi_I'].concentration - self.p['d2'].value * self.s['T'].concentration

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['K_T'].amount, self.s['phi_I'].amount, self.s['phi_R'].amount, self.s['T'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ 0, self.Ratephi_I(), self.Ratephi_R(), self.RateT()], dtype = np.float64)

		stoichiometricMatrix = np.array([[0.],[0.],[0.],[0.]], dtype = np.float64)

		reactionVelocities = np.array([0], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['K_T'].amount, self.s['phi_I'].amount, self.s['phi_R'].amount, self.s['T'].amount], dtype = np.float64)
		self.s['K_T'].amount, self.s['phi_I'].amount, self.s['phi_R'].amount, self.s['T'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

