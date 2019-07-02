from sbmltopyode.SBMLModelClasses import *
from scipy.integrate import odeint
import numpy as np
import operator
import math

class Lenbury2001:

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['time_environment'] = Parameter(0.0, 'time_environment', True)
		self.p['x'] = Parameter(5.0, 'x', False)
		self.p['r1'] = Parameter(0.15, 'r1', True)
		self.p['r2'] = Parameter(0.12, 'r2', True)
		self.p['c1'] = Parameter(0.1, 'c1', True)
		self.p['y'] = Parameter(0.0, 'y', False)
		self.p['r3'] = Parameter(0.05, 'r3', True)
		self.p['r4'] = Parameter(0.03, 'r4', True)
		self.p['c2'] = Parameter(0.1, 'c2', True)
		self.p['c3'] = Parameter(0.005, 'c3', True)
		self.p['z'] = Parameter(1.0, 'z', False)
		self.p['r5'] = Parameter(0.09, 'r5', True)
		self.p['r6'] = Parameter(0.1, 'r6', True)
		self.p['r7'] = Parameter(0.05, 'r7', True)
		self.p['z_'] = Parameter(1.01, 'z_', True)
		self.p['y_'] = Parameter(1.08, 'y_', True)
		self.p['delta'] = Parameter(0.01, 'delta', True)
		self.p['u'] = Parameter(1.0, 'u', False)
		self.p['omega'] = Parameter(0.05, 'omega', True)
		self.p['v'] = Parameter(0.0, 'v', False)
		self.p['epsilon'] = Parameter(0.1, 'epsilon', True)

		self.c = {} #Dictionary of compartments
		self.c['COMpartment'] = Compartment(1.0, 3, True)

		self.s = {} #Dictionary of chemical species

		self.r = {} #Dictionary of reactiions
		self.time = 0

		self.reactionMetadata = {
		}
		self.AssignmentRules()



	def AssignmentRules(self):

		return

	def Ratex(self):

		return self.p['z'].value * (self.p['r1'].value * self.p['y'].value + -self.p['r2'].value * self.p['x'].value + self.p['c1'].value)

	def Ratey(self):

		return self.p['epsilon'].value * (self.p['r3'].value / self.p['z'].value + -self.p['r4'].value * self.p['x'].value + self.p['c2'].value + self.p['c3'].value * self.p['u'].value)

	def Ratez(self):

		return self.p['epsilon'].value * self.p['delta'].value * (self.p['r5'].value * (self.p['y'].value - self.p['y_'].value) * (self.p['z_'].value - self.p['z'].value) + self.p['r6'].value * self.p['z'].value * (self.p['z_'].value - self.p['z'].value) - self.p['r7'].value * self.p['z'].value)

	def Rateu(self):

		return -self.p['omega'].value * self.p['v'].value

	def Ratev(self):

		return self.p['omega'].value * self.p['u'].value

	def _SolveReactions(self, y, t):

		self.time = t
		self.p['x'].value, self.p['y'].value, self.p['z'].value, self.p['u'].value, self.p['v'].value = y
		self.AssignmentRules()

		rateRuleVector = np.array([ self.Ratex(), self.Ratey(), self.Ratez(), self.Rateu(), self.Ratev()], dtype = np.float64)

		stoichiometricMatrix = np.array([[0.],[0.],[0.],[0.],[0.]], dtype = np.float64)

		reactionVelocities = np.array([0], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.p['x'].value, self.p['y'].value, self.p['z'].value, self.p['u'].value, self.p['v'].value], dtype = np.float64)
		self.p['x'].value, self.p['y'].value, self.p['z'].value, self.p['u'].value, self.p['v'].value = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

