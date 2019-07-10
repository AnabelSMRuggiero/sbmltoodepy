from sbmltopyode.SBMLModelClasses import *
from scipy.integrate import odeint
import numpy as np
import operator
import math

class Guyton1972:

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['MDFLW'] = Parameter(1.00051, 'MDFLW', True)
		self.p['ANGSCR'] = Parameter(0.0, 'ANGSCR', False)
		self.p['MDFLW3'] = Parameter(0.0, 'MDFLW3', False)
		self.p['ANX1'] = Parameter(0.0, 'ANX1', False)
		self.p['ANX'] = Parameter(0.0, 'ANX', False)
		self.p['ANPR'] = Parameter(0.0, 'ANPR', False)
		self.p['ANPRT'] = Parameter(0.0, 'ANPRT', False)
		self.p['ANPR1'] = Parameter(0.0, 'ANPR1', False)
		self.p['ANC'] = Parameter(0.859476, 'ANC', False)
		self.p['ANM'] = Parameter(0.0, 'ANM', False)
		self.p['ANU'] = Parameter(0.0, 'ANU', False)
		self.p['ANU1'] = Parameter(0.0, 'ANU1', False)
		self.p['ANUVN'] = Parameter(0.0, 'ANUVN', False)
		self.p['ANXM'] = Parameter(0.0, 'ANXM', True)
		self.p['ANV'] = Parameter(5000.0, 'ANV', True)
		self.p['REK'] = Parameter(1.0, 'REK', True)
		self.p['ANGINF'] = Parameter(0.0, 'ANGINF', True)
		self.p['ANGKNS'] = Parameter(0.0, 'ANGKNS', True)
		self.p['ANT'] = Parameter(12.0, 'ANT', True)
		self.p['Z12'] = Parameter(5.0, 'Z12', True)
		self.p['ANMUL'] = Parameter(1.8, 'ANMUL', True)
		self.p['ANMLL'] = Parameter(0.7, 'ANMLL', True)
		self.p['ANCSNS'] = Parameter(0.4, 'ANCSNS', True)
		self.p['ANULL'] = Parameter(0.8, 'ANULL', True)
		self.p['ANUM'] = Parameter(6.0, 'ANUM', True)
		self.p['ANUVM'] = Parameter(0.0, 'ANUVM', True)
		self.p['tu'] = Parameter(1.0, 'tu', True)

		self.c = {} #Dictionary of compartments
		self.c['Compartment'] = Compartment(1.0, 3, True)

		self.s = {} #Dictionary of chemical species

		self.r = {} #Dictionary of reactiions
		self.time = 0

		self.reactionMetadata = {
		}
		self.AssignmentRules()



	def AssignmentRules(self):

		self.p['MDFLW3'].value = self.p['MDFLW'].value

		self.p['ANGSCR'].value = Piecewise(1 / (1 + (self.p['MDFLW3'].value - 1) * 72), self.p['MDFLW3'].value > 1, 10 - 9 / (1 + (1 - self.p['MDFLW3'].value) * 8))

		self.p['ANX'].value = (self.p['ANGSCR'].value - 1) * self.p['ANXM'].value

		self.p['ANPRT'].value = (self.p['ANGSCR'].value + self.p['ANX1'].value) * self.p['REK'].value

		self.p['ANPR'].value = Piecewise(1e-5, self.p['ANPRT'].value < 1e-5, self.p['ANPRT'].value)

		self.p['ANPR1'].value = Piecewise(self.p['ANGKNS'].value, self.p['ANGKNS'].value > 0, self.p['ANPR'].value + self.p['ANGINF'].value)

		self.p['ANM'].value = self.p['ANMUL'].value - (self.p['ANMUL'].value - 1) / (((self.p['ANMLL'].value - 1) / (self.p['ANMLL'].value - self.p['ANMUL'].value)) * (self.p['ANC'].value - 1) * self.p['ANCSNS'].value + 1)

		self.p['ANU1'].value = (self.p['ANM'].value - 1) * self.p['ANUM'].value + 1

		self.p['ANU'].value = Piecewise(self.p['ANULL'].value, self.p['ANU1'].value < self.p['ANULL'].value, self.p['ANU1'].value)

		self.p['ANUVN'].value = (self.p['ANU'].value - 1) * self.p['ANUVM'].value + 1

		return

	def RateANX1(self):

		return (self.p['ANX'].value - self.p['ANX1'].value) / (self.p['tu'].value * self.p['ANV'].value)

	def RateANC(self):

		return (self.p['ANPR1'].value - self.p['ANC'].value) / (self.p['tu'].value * self.p['ANT'].value)

	def _SolveReactions(self, y, t):

		self.time = t
		self.p['ANX1'].value, self.p['ANC'].value = y
		self.AssignmentRules()

		rateRuleVector = np.array([ self.RateANX1(), self.RateANC()], dtype = np.float64)

		stoichiometricMatrix = np.array([[0.],[0.]], dtype = np.float64)

		reactionVelocities = np.array([0], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.p['ANX1'].value, self.p['ANC'].value], dtype = np.float64)
		self.p['ANX1'].value, self.p['ANC'].value = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

