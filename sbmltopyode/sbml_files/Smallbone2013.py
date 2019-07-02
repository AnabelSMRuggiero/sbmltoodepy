from sbmltopyode.SBMLModelClasses import *
from scipy.integrate import odeint
import numpy as np
import operator
import math

class Smallbone2013:

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['d0'] = Parameter(0.1, 'd0', True)
		self.p['d1'] = Parameter(0.420467092599869, 'd1', True)
		self.p['d2'] = Parameter(1.10138534772246, 'd2', True)
		self.p['T'] = Parameter(0.0, 'T', False)
		self.p['K'] = Parameter(120.0, 'K', True)
		self.p['r0'] = Parameter(1.83898685224665, 'r0', True)
		self.p['f0'] = Parameter(0.0, 'f0', False)
		self.p['p00'] = Parameter(0.0, 'p00', False)
		self.p['p01'] = Parameter(0.855699855699856, 'p01', True)
		self.p['r1'] = Parameter(5.88010232010212, 'r1', True)
		self.p['f1'] = Parameter(0.0, 'f1', False)
		self.p['p11'] = Parameter(0.0, 'p11', False)
		self.p['p12'] = Parameter(0.827377484810943, 'p12', True)

		self.c = {} #Dictionary of compartments
		self.c['compartment'] = Compartment(1.0, 3, True)

		self.s = {} #Dictionary of chemical species
		speciesMetadata = SBMLMetadata('N0')
		self.s['N0'] = Species(1.75444831412765, 'Amount', self.c['compartment'], False, constant = False)
		speciesMetadata = SBMLMetadata('N1')
		self.s['N1'] = Species(27.40585059, 'Amount', self.c['compartment'], False, constant = False)
		speciesMetadata = SBMLMetadata('N2')
		self.s['N2'] = Species(45.6191494109, 'Amount', self.c['compartment'], False, constant = False)

		self.r = {} #Dictionary of reactiions
		self.r['R0X'] = R0X(self, SBMLMetadata('N0 death'))
		self.r['R01'] = R01(self, SBMLMetadata('N0 differentiation'))
		self.r['R00'] = R00(self, SBMLMetadata('N0 renewal'))
		self.r['R1X'] = R1X(self, SBMLMetadata('N1 death'))
		self.r['R12'] = R12(self, SBMLMetadata('N1 differentiation'))
		self.r['R11'] = R11(self, SBMLMetadata('N1 renewal'))
		self.r['R2X'] = R2X(self, SBMLMetadata('N2 death'))
		self.time = 0

		self.reactionMetadata = {
		self.ReactionR0X: SBMLMetadata('N0 death'),
		self.ReactionR01: SBMLMetadata('N0 differentiation'),
		self.ReactionR00: SBMLMetadata('N0 renewal'),
		self.ReactionR1X: SBMLMetadata('N1 death'),
		self.ReactionR12: SBMLMetadata('N1 differentiation'),
		self.ReactionR11: SBMLMetadata('N1 renewal'),
		self.ReactionR2X: SBMLMetadata('N2 death')
		}
		self.AssignmentRules()



	def AssignmentRules(self):

		self.p['T'].value = self.s['N0'].concentration + self.s['N1'].concentration + self.s['N2'].concentration

		self.p['f0'].value = self.p['r0'].value * self.s['N0'].concentration * (1 - self.p['T'].value / self.p['K'].value)

		self.p['p00'].value = 1 - self.p['p01'].value

		self.p['f1'].value = self.p['r1'].value * self.s['N1'].concentration * (1 - self.p['T'].value / self.p['K'].value)

		self.p['p11'].value = 1 - self.p['p12'].value

		return

	def ReactionR0X(self):

		return self.p['d0'].value * self.s['N0'].concentration

	def ReactionR01(self):

		return self.p['p01'].value * self.p['f0'].value

	def ReactionR00(self):

		return self.p['p00'].value * self.p['f0'].value

	def ReactionR1X(self):

		return self.p['d1'].value * self.s['N1'].concentration

	def ReactionR12(self):

		return self.p['p12'].value * self.p['f1'].value

	def ReactionR11(self):

		return self.p['p11'].value * self.p['f1'].value

	def ReactionR2X(self):

		return self.p['d2'].value * self.s['N2'].concentration

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['N0'].amount, self.s['N1'].amount, self.s['N2'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ 0, 0, 0], dtype = np.float64)

		stoichiometricMatrix = np.array([[-1,0,1,0,0,0,0.],[ 0,1,0,-1,0,1,0.],[ 0,0,0,0,1,0,-1.]], dtype = np.float64)

		reactionVelocities = np.array([self.r['R0X'](), self.r['R01'](), self.r['R00'](), self.r['R1X'](), self.r['R12'](), self.r['R11'](), self.r['R2X']()], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['N0'].amount, self.s['N1'].amount, self.s['N2'].amount], dtype = np.float64)
		self.s['N0'].amount, self.s['N1'].amount, self.s['N2'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

class R0X:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['d0'].value * self.parent.s['N0'].concentration

class R01:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['p01'].value * self.parent.p['f0'].value

class R00:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['p00'].value * self.parent.p['f0'].value

class R1X:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['d1'].value * self.parent.s['N1'].concentration

class R12:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['p12'].value * self.parent.p['f1'].value

class R11:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['p11'].value * self.parent.p['f1'].value

class R2X:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['d2'].value * self.parent.s['N2'].concentration

