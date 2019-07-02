from sbmltopyode.SBMLModelClasses import *
from scipy.integrate import odeint
import numpy as np
import operator
import math

class Cizmeci2013:

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['K1'] = Parameter(0.05, 'K1', True)
		self.p['K2'] = Parameter(0.05, 'K2', True)
		self.p['d1'] = Parameter(1.0, 'd1', True)
		self.p['d2'] = Parameter(1.0, 'd2', True)
		self.p['k1'] = Parameter(10.0, 'k1', True)
		self.p['k2'] = Parameter(10.0, 'k2', True)
		self.p['k3'] = Parameter(0.1, 'k3', True)
		self.p['k4'] = Parameter(0.01, 'k4', True)
		self.p['k5'] = Parameter(0.1, 'k5', True)
		self.p['k6'] = Parameter(10.0, 'k6', True)
		self.p['k7'] = Parameter(0.1, 'k7', True)
		self.p['k8'] = Parameter(0.1, 'k8', True)
		self.p['delta'] = Parameter(1.0, 'delta', True)
		self.p['beta'] = Parameter(1.0, 'beta', True)
		self.p['deltaon'] = Parameter(0.01, 'deltaon', True)
		self.p['deltano'] = Parameter(0.01, 'deltano', True)
		self.p['deltaang'] = Parameter(1.0, 'deltaang', True)
		self.p['deltar'] = Parameter(0.01, 'deltar', True)
		self.p['lambda'] = Parameter(0.5, 'lambda', True)
		self.p['phi'] = Parameter(1.0, 'phi', True)
		self.p['eta'] = Parameter(0.01, 'eta', True)
		self.p['psi'] = Parameter(0.559, 'psi', True)
		self.p['E2T'] = Parameter(1.0, 'E2T', True)

		self.c = {} #Dictionary of compartments
		self.c['default'] = Compartment(1.0, 3, True)

		self.s = {} #Dictionary of chemical species
		speciesMetadata = SBMLMetadata('AKT')
		self.s['s1'] = Species(0.2, 'Amount', self.c['default'], False, constant = False)
		speciesMetadata = SBMLMetadata('pIRS1')
		self.s['s2'] = Species(1.2955, 'Amount', self.c['default'], False, constant = False)
		speciesMetadata = SBMLMetadata('NO')
		self.s['s3'] = Species(0.013, 'Amount', self.c['default'], False, constant = False)
		speciesMetadata = SBMLMetadata('ANGII')
		self.s['s4'] = Species(0.13, 'Amount', self.c['default'], True, constant = False)
		speciesMetadata = SBMLMetadata('ROS')
		self.s['s5'] = Species(1.1504, 'Amount', self.c['default'], True, constant = False)
		speciesMetadata = SBMLMetadata('ONOO')
		self.s['s6'] = Species(0.1496, 'Amount', self.c['default'], True, constant = False)
		speciesMetadata = SBMLMetadata('sa3_degraded')
		self.s['s7'] = Species(0.0, 'Amount', self.c['default'], False, constant = False)
		speciesMetadata = SBMLMetadata('sa4_degraded')
		self.s['s8'] = Species(0.0, 'Amount', self.c['default'], False, constant = False)
		speciesMetadata = SBMLMetadata('sa5_degraded')
		self.s['s9'] = Species(0.0, 'Amount', self.c['default'], False, constant = False)
		speciesMetadata = SBMLMetadata('sa6_degraded')
		self.s['s10'] = Species(0.0, 'Amount', self.c['default'], False, constant = False)

		self.r = {} #Dictionary of reactiions
		self.r['re2'] = re2(self, SBMLMetadata('activation of AKT by pIRS1'))
		self.r['re4'] = re4(self, SBMLMetadata('activation of AKT by ONOO'))
		self.r['re5'] = re5(self, SBMLMetadata('activation of pIRS1 by AKT and ANGII'))
		self.r['re6'] = re6(self, SBMLMetadata('activation of NO by AKT'))
		self.r['re7'] = re7(self, SBMLMetadata('production of ONOO from NO and ROS'))
		self.r['re8'] = re8(self, SBMLMetadata('degradation of NO'))
		self.r['re9'] = re9(self, SBMLMetadata('activation of ANG II from NO'))
		self.r['re10'] = re10(self, SBMLMetadata('degradation of ANGII'))
		self.r['re11'] = re11(self, SBMLMetadata('activation of ROS from ANGII'))
		self.r['re12'] = re12(self, SBMLMetadata('degradation of ROS'))
		self.r['re13'] = re13(self, SBMLMetadata('degradation of ONOO'))
		self.time = 0

		self.reactionMetadata = {
		self.Reactionre2: SBMLMetadata('activation of AKT by pIRS1'),
		self.Reactionre4: SBMLMetadata('activation of AKT by ONOO'),
		self.Reactionre5: SBMLMetadata('activation of pIRS1 by AKT and ANGII'),
		self.Reactionre6: SBMLMetadata('activation of NO by AKT'),
		self.Reactionre7: SBMLMetadata('production of ONOO from NO and ROS'),
		self.Reactionre8: SBMLMetadata('degradation of NO'),
		self.Reactionre9: SBMLMetadata('activation of ANG II from NO'),
		self.Reactionre10: SBMLMetadata('degradation of ANGII'),
		self.Reactionre11: SBMLMetadata('activation of ROS from ANGII'),
		self.Reactionre12: SBMLMetadata('degradation of ROS'),
		self.Reactionre13: SBMLMetadata('degradation of ONOO')
		}
		self.AssignmentRules()



	def AssignmentRules(self):

		return

	def Reactionre2(self):

		return -(((self.p['d1'].value + self.p['K1'].value) / self.p['K2'].value) * self.p['beta'].value * self.s['s1'].concentration * self.s['s2'].concentration) + (self.p['d1'].value + self.p['K1'].value) * self.p['beta'].value * self.s['s2'].concentration * self.s['s1'].concentration**2 / (self.p['K2'].value**2 + self.p['K2'].value * self.s['s1'].concentration) + self.p['d1'].value * self.p['beta'].value * self.s['s2'].concentration * self.s['s1'].concentration / (self.s['s1'].concentration + self.p['K2'].value) + self.p['K2'].value * self.p['E2T'].value * (1 - self.s['s1'].concentration) / ((1 - self.s['s1'].concentration) + self.p['K2'].value)

	def Reactionre4(self):

		return self.p['k3'].value * self.s['s6'].concentration

	def Reactionre5(self):

		return ((self.p['delta'].value / self.p['beta'].value) * self.p['K2'].value / self.p['K1'].value) * self.p['E2T'].value * self.p['lambda'].value + (self.p['phi'].value - self.p['eta'].value * (1 + self.p['k7'].value * self.s['s4'].concentration) * self.p['psi'].value) * (1 - self.s['s1'].concentration) - self.p['delta'].value * (self.p['d1'].value + self.p['K1'].value) * self.s['s2'].concentration

	def Reactionre6(self):

		return self.p['k4'].value * (1 - self.s['s1'].concentration)

	def Reactionre7(self):

		return self.p['k8'].value * self.s['s3'].concentration * self.s['s5'].concentration

	def Reactionre8(self):

		return self.p['deltano'].value * self.s['s3'].concentration

	def Reactionre9(self):

		return self.p['k6'].value * self.s['s3'].concentration

	def Reactionre10(self):

		return self.p['deltaang'].value * self.s['s4'].concentration

	def Reactionre11(self):

		return self.p['k5'].value * self.s['s4'].concentration

	def Reactionre12(self):

		return self.p['deltar'].value * self.s['s5'].concentration

	def Reactionre13(self):

		return self.p['deltaon'].value * self.s['s6'].concentration

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['s1'].amount, self.s['s2'].amount, self.s['s3'].amount, self.s['s4'].amount, self.s['s5'].amount, self.s['s6'].amount, self.s['s7'].amount, self.s['s8'].amount, self.s['s9'].amount, self.s['s10'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = np.float64)

		stoichiometricMatrix = np.array([[ 1,1,-1,-1,0,0,0,0,0,0,0.],[-1,0,1,0,0,0,0,0,0,0,0.],[ 0,0,0,1,-1,-1,-1,0,0,0,0.],[ 0,0,-1,0,0,0,1,-1,-1,0,0.],[ 0,0,0,0,-1,0,0,0,1,-1,0.],[ 0,-1,0,0,1,0,0,0,0,0,-1.],[ 0,0,0,0,0,1,0,0,0,0,0.],[ 0,0,0,0,0,0,0,1,0,0,0.],[ 0,0,0,0,0,0,0,0,0,1,0.],[ 0,0,0,0,0,0,0,0,0,0,1.]], dtype = np.float64)

		reactionVelocities = np.array([self.r['re2'](), self.r['re4'](), self.r['re5'](), self.r['re6'](), self.r['re7'](), self.r['re8'](), self.r['re9'](), self.r['re10'](), self.r['re11'](), self.r['re12'](), self.r['re13']()], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['s1'].amount, self.s['s2'].amount, self.s['s3'].amount, self.s['s4'].amount, self.s['s5'].amount, self.s['s6'].amount, self.s['s7'].amount, self.s['s8'].amount, self.s['s9'].amount, self.s['s10'].amount], dtype = np.float64)
		self.s['s1'].amount, self.s['s2'].amount, self.s['s3'].amount, self.s['s4'].amount, self.s['s5'].amount, self.s['s6'].amount, self.s['s7'].amount, self.s['s8'].amount, self.s['s9'].amount, self.s['s10'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

class re2:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return -(((self.parent.p['d1'].value + self.parent.p['K1'].value) / self.parent.p['K2'].value) * self.parent.p['beta'].value * self.parent.s['s1'].concentration * self.parent.s['s2'].concentration) + (self.parent.p['d1'].value + self.parent.p['K1'].value) * self.parent.p['beta'].value * self.parent.s['s2'].concentration * self.parent.s['s1'].concentration**2 / (self.parent.p['K2'].value**2 + self.parent.p['K2'].value * self.parent.s['s1'].concentration) + self.parent.p['d1'].value * self.parent.p['beta'].value * self.parent.s['s2'].concentration * self.parent.s['s1'].concentration / (self.parent.s['s1'].concentration + self.parent.p['K2'].value) + self.parent.p['K2'].value * self.parent.p['E2T'].value * (1 - self.parent.s['s1'].concentration) / ((1 - self.parent.s['s1'].concentration) + self.parent.p['K2'].value)

class re4:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k3'].value * self.parent.s['s6'].concentration

class re5:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return ((self.parent.p['delta'].value / self.parent.p['beta'].value) * self.parent.p['K2'].value / self.parent.p['K1'].value) * self.parent.p['E2T'].value * self.parent.p['lambda'].value + (self.parent.p['phi'].value - self.parent.p['eta'].value * (1 + self.parent.p['k7'].value * self.parent.s['s4'].concentration) * self.parent.p['psi'].value) * (1 - self.parent.s['s1'].concentration) - self.parent.p['delta'].value * (self.parent.p['d1'].value + self.parent.p['K1'].value) * self.parent.s['s2'].concentration

class re6:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k4'].value * (1 - self.parent.s['s1'].concentration)

class re7:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k8'].value * self.parent.s['s3'].concentration * self.parent.s['s5'].concentration

class re8:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['deltano'].value * self.parent.s['s3'].concentration

class re9:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k6'].value * self.parent.s['s3'].concentration

class re10:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['deltaang'].value * self.parent.s['s4'].concentration

class re11:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k5'].value * self.parent.s['s4'].concentration

class re12:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['deltar'].value * self.parent.s['s5'].concentration

class re13:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['deltaon'].value * self.parent.s['s6'].concentration

