from sbmltopyode.SBMLModelClasses import *
from scipy.integrate import odeint
import numpy as np
import operator
import math

class Zi2011:

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['stimulation_type'] = Parameter(1.0, 'stimulation_type', True)
		self.p['single_pulse_duration'] = Parameter(0.5, 'single_pulse_duration', True)
		self.p['totalNumT1R'] = Parameter(10005.0, 'totalNumT1R', False)
		self.p['totalNumT2R'] = Parameter(2272.0, 'totalNumT2R', False)
		self.p['totalNumLRC'] = Parameter(0.0, 'totalNumLRC', False)
		self.p['totalNumPSmad2'] = Parameter(0.0, 'totalNumPSmad2', False)
		self.p['totalNuclearPSmad2'] = Parameter(0.0, 'totalNuclearPSmad2', False)
		self.p['totalSmad2c'] = Parameter(60.6, 'totalSmad2c', False)
		self.p['totalSmad2n'] = Parameter(28.5, 'totalSmad2n', False)
		self.p['medium_TGF_beta_amount'] = Parameter(0.0, 'medium_TGF_beta_amount', False)
		self.p['TGF_beta_dose_mol_per_cell'] = Parameter(0.0, 'TGF_beta_dose_mol_per_cell', False)
		self.p['ki'] = Parameter(0.333, 'ki', True)
		self.p['kr'] = Parameter(0.0333, 'kr', True)
		self.p['k_T1R'] = Parameter(0.0167, 'k_T1R', True)
		self.p['k_T2R'] = Parameter(0.0190076, 'k_T2R', True)
		self.p['kdeg_T1R'] = Parameter(0.00256, 'kdeg_T1R', True)
		self.p['kdeg_T2R'] = Parameter(0.0132, 'kdeg_T2R', True)
		self.p['kdeg_LRC'] = Parameter(0.00256, 'kdeg_LRC', True)
		self.p['kdeg_TGF_beta'] = Parameter(0.347, 'kdeg_TGF_beta', True)
		self.p['klid'] = Parameter(0.0233678, 'klid', True)
		self.p['ka_LRC'] = Parameter(117.897, 'ka_LRC', True)
		self.p['kdiss_LRC'] = Parameter(0.0438111, 'kdiss_LRC', True)
		self.p['kimp_Smad2'] = Parameter(0.156, 'kimp_Smad2', True)
		self.p['kexp_Smad2'] = Parameter(0.763, 'kexp_Smad2', True)
		self.p['kimp_Smad4'] = Parameter(0.156, 'kimp_Smad4', True)
		self.p['kexp_Smad4'] = Parameter(0.359, 'kexp_Smad4', True)
		self.p['kpho_Smad2'] = Parameter(0.0488268, 'kpho_Smad2', True)
		self.p['kon_Smads'] = Parameter(0.198472, 'kon_Smads', True)
		self.p['koff_Smads'] = Parameter(1.0, 'koff_Smads', True)
		self.p['kimp_Smads'] = Parameter(0.889, 'kimp_Smads', True)
		self.p['kdepho_Smad2'] = Parameter(0.394, 'kdepho_Smad2', True)
		self.p['kon_ns'] = Parameter(0.0505413, 'kon_ns', True)
		self.p['koff_ns'] = Parameter(2.03306, 'koff_ns', False)
		self.p['KD_ns'] = Parameter(40.2257, 'KD_ns', True)

		self.c = {} #Dictionary of compartments
		self.c['default'] = Compartment(1.0, 3, True)
		self.c['Vmed'] = Compartment(2e-09, 3, False)
		self.c['Vcyt'] = Compartment(2.3e-12, 3, True)
		self.c['Vnuc'] = Compartment(1e-12, 3, True)

		self.s = {} #Dictionary of chemical species
		speciesMetadata = SBMLMetadata('TGF_beta_ex')
		self.s['TGF_beta_ex'] = Species(0.05, 'Concentration', self.c['Vmed'], False, constant = False)
		speciesMetadata = SBMLMetadata('T1R_surf')
		self.s['T1R_surf'] = Species(0.702494, 'Concentration', self.c['Vcyt'], False, constant = False)
		speciesMetadata = SBMLMetadata('T1R_endo')
		self.s['T1R_endo'] = Species(6.52344, 'Concentration', self.c['Vcyt'], False, constant = False)
		speciesMetadata = SBMLMetadata('T2R_surf')
		self.s['T2R_surf'] = Species(0.201077, 'Concentration', self.c['Vcyt'], False, constant = False)
		speciesMetadata = SBMLMetadata('T2R_endo')
		self.s['T2R_endo'] = Species(1.43997, 'Concentration', self.c['Vcyt'], False, constant = False)
		speciesMetadata = SBMLMetadata('LRC_surf')
		self.s['LRC_surf'] = Species(0.0, 'Concentration', self.c['Vcyt'], False, constant = False)
		speciesMetadata = SBMLMetadata('LRC_endo')
		self.s['LRC_endo'] = Species(0.0, 'Concentration', self.c['Vcyt'], False, constant = False)
		speciesMetadata = SBMLMetadata('Smad2c')
		self.s['Smad2c'] = Species(60.6, 'Concentration', self.c['Vcyt'], False, constant = False)
		speciesMetadata = SBMLMetadata('Smad2n')
		self.s['Smad2n'] = Species(28.5, 'Concentration', self.c['Vnuc'], False, constant = False)
		speciesMetadata = SBMLMetadata('Smad4c')
		self.s['Smad4c'] = Species(50.8, 'Concentration', self.c['Vcyt'], False, constant = False)
		speciesMetadata = SBMLMetadata('Smad4n')
		self.s['Smad4n'] = Species(50.8, 'Concentration', self.c['Vnuc'], False, constant = False)
		speciesMetadata = SBMLMetadata('PSmad2c')
		self.s['PSmad2c'] = Species(0.0, 'Concentration', self.c['Vcyt'], False, constant = False)
		speciesMetadata = SBMLMetadata('PSmad2_PSmad2_c')
		self.s['PSmad2_PSmad2_c'] = Species(0.0, 'Concentration', self.c['Vcyt'], False, constant = False)
		speciesMetadata = SBMLMetadata('PSmad2_Smad4_c')
		self.s['PSmad2_Smad4_c'] = Species(0.0, 'Concentration', self.c['Vcyt'], False, constant = False)
		speciesMetadata = SBMLMetadata('PSmad2n')
		self.s['PSmad2n'] = Species(0.0, 'Concentration', self.c['Vnuc'], False, constant = False)
		speciesMetadata = SBMLMetadata('PSmad2_PSmad2_n')
		self.s['PSmad2_PSmad2_n'] = Species(0.0, 'Concentration', self.c['Vnuc'], False, constant = False)
		speciesMetadata = SBMLMetadata('PSmad2_Smad4_n')
		self.s['PSmad2_Smad4_n'] = Species(0.0, 'Concentration', self.c['Vnuc'], False, constant = False)
		speciesMetadata = SBMLMetadata('TGF_beta_endo')
		self.s['TGF_beta_endo'] = Species(0.0, 'Concentration', self.c['Vcyt'], False, constant = False)
		speciesMetadata = SBMLMetadata('TGF_beta_ns')
		self.s['TGF_beta_ns'] = Species(0.0, 'Concentration', self.c['Vmed'], False, constant = False)
		speciesMetadata = SBMLMetadata('AA')
		self.s['AA'] = Species(0.0, 'Concentration', self.c['Vcyt'], False, constant = True)
		speciesMetadata = SBMLMetadata('empty_degraded')
		self.s['empty_degraded'] = Species(0.0, 'Concentration', self.c['Vcyt'], False, constant = True)

		self.r = {} #Dictionary of reactiions
		self.r['re1'] = re1(self, SBMLMetadata('T1R production'))
		self.r['re2'] = re2(self, SBMLMetadata('T1R internalization to early endosome'))
		self.r['re3'] = re3(self, SBMLMetadata('T1R recycling from early endosome'))
		self.r['re4'] = re4(self, SBMLMetadata('T1R constitutive degradation'))
		self.r['re5'] = re5(self, SBMLMetadata('T2R production'))
		self.r['re6'] = re6(self, SBMLMetadata('T2R internalization to early endosome'))
		self.r['re7'] = re7(self, SBMLMetadata('T2R recycling from early endosome'))
		self.r['re8'] = re8(self, SBMLMetadata('T2R constitutive degradation'))
		self.r['re9'] = re9(self, SBMLMetadata('LRC activation'))
		self.r['re10'] = re10(self, SBMLMetadata('LRC internalization to early endosome'))
		self.r['re11'] = re11(self, SBMLMetadata('LRC constitutive degradation'))
		self.r['re12'] = re12(self, SBMLMetadata('dissociation of LRC in endosome'))
		self.r['re13'] = re13(self, SBMLMetadata('TGF-beta constitutive degradation'))
		self.r['re14'] = re14(self, SBMLMetadata('Smad2 nuclear import'))
		self.r['re15'] = re15(self, SBMLMetadata('Smad2 nuclear export'))
		self.r['re16'] = re16(self, SBMLMetadata('Smad4 nuclear import'))
		self.r['re17'] = re17(self, SBMLMetadata('Smad4 nuclear export'))
		self.r['re18'] = re18(self, SBMLMetadata('Smad2 phosphorylation'))
		self.r['re19'] = re19(self, SBMLMetadata('PSmad2 nuclear import'))
		self.r['re20'] = re20(self, SBMLMetadata('PSmad2 nuclear export'))
		self.r['re21'] = re21(self, SBMLMetadata('Smad2-Smad4 complex formation'))
		self.r['re22'] = re22(self, SBMLMetadata('Smad2-Smad4 nuclear import'))
		self.r['re23'] = re23(self, SBMLMetadata('Smad2-Smad4 dissociation'))
		self.r['re24'] = re24(self, SBMLMetadata('Smad2 dephosphorylation'))
		self.r['re25'] = re25(self, SBMLMetadata('PSmad2 dimer formation'))
		self.r['re26'] = re26(self, SBMLMetadata('PSmad2 dimer nuclear  import'))
		self.r['re27'] = re27(self, SBMLMetadata('PSmad2 dimmer dissociation'))
		self.r['re28'] = re28(self, SBMLMetadata('negative feedback induced LRC degradation'))
		self.r['re29'] = re29(self, SBMLMetadata('non-specific binding of TGF-beta'))
		self.time = 0

		self.reactionMetadata = {
		self.Reactionre1: SBMLMetadata('T1R production'),
		self.Reactionre2: SBMLMetadata('T1R internalization to early endosome'),
		self.Reactionre3: SBMLMetadata('T1R recycling from early endosome'),
		self.Reactionre4: SBMLMetadata('T1R constitutive degradation'),
		self.Reactionre5: SBMLMetadata('T2R production'),
		self.Reactionre6: SBMLMetadata('T2R internalization to early endosome'),
		self.Reactionre7: SBMLMetadata('T2R recycling from early endosome'),
		self.Reactionre8: SBMLMetadata('T2R constitutive degradation'),
		self.Reactionre9: SBMLMetadata('LRC activation'),
		self.Reactionre10: SBMLMetadata('LRC internalization to early endosome'),
		self.Reactionre11: SBMLMetadata('LRC constitutive degradation'),
		self.Reactionre12: SBMLMetadata('dissociation of LRC in endosome'),
		self.Reactionre13: SBMLMetadata('TGF-beta constitutive degradation'),
		self.Reactionre14: SBMLMetadata('Smad2 nuclear import'),
		self.Reactionre15: SBMLMetadata('Smad2 nuclear export'),
		self.Reactionre16: SBMLMetadata('Smad4 nuclear import'),
		self.Reactionre17: SBMLMetadata('Smad4 nuclear export'),
		self.Reactionre18: SBMLMetadata('Smad2 phosphorylation'),
		self.Reactionre19: SBMLMetadata('PSmad2 nuclear import'),
		self.Reactionre20: SBMLMetadata('PSmad2 nuclear export'),
		self.Reactionre21: SBMLMetadata('Smad2-Smad4 complex formation'),
		self.Reactionre22: SBMLMetadata('Smad2-Smad4 nuclear import'),
		self.Reactionre23: SBMLMetadata('Smad2-Smad4 dissociation'),
		self.Reactionre24: SBMLMetadata('Smad2 dephosphorylation'),
		self.Reactionre25: SBMLMetadata('PSmad2 dimer formation'),
		self.Reactionre26: SBMLMetadata('PSmad2 dimer nuclear  import'),
		self.Reactionre27: SBMLMetadata('PSmad2 dimmer dissociation'),
		self.Reactionre28: SBMLMetadata('negative feedback induced LRC degradation'),
		self.Reactionre29: SBMLMetadata('non-specific binding of TGF-beta')
		}
		self.AssignmentRules()



	def AssignmentRules(self):

		self.c['Vmed'].size = 2e-3 / (1e6 * math.exp(math.log(1.45) * self.time / 1440))

		self.p['totalNumPSmad2'].value = (self.s['PSmad2c'].concentration + self.s['PSmad2_PSmad2_c'].concentration * 2 + self.s['PSmad2_Smad4_c'].concentration) * 2.3 * 602 + (self.s['PSmad2n'].concentration + self.s['PSmad2_PSmad2_n'].concentration * 2 + self.s['PSmad2_Smad4_n'].concentration) * 602

		self.p['totalNuclearPSmad2'].value = self.s['PSmad2n'].concentration + 2 * self.s['PSmad2_PSmad2_n'].concentration + self.s['PSmad2_Smad4_n'].concentration

		self.p['totalNumLRC'].value = (self.s['LRC_surf'].concentration + self.s['LRC_endo'].concentration) * 2.3 * 602

		self.p['totalSmad2c'].value = self.s['Smad2c'].concentration + self.s['PSmad2c'].concentration + 2 * self.s['PSmad2_PSmad2_c'].concentration + self.s['PSmad2_Smad4_c'].concentration

		self.p['totalSmad2n'].value = self.s['Smad2n'].concentration + self.s['PSmad2n'].concentration + 2 * self.s['PSmad2_PSmad2_n'].concentration + self.s['PSmad2_Smad4_n'].concentration

		self.p['koff_ns'].value = self.p['kon_ns'].value * self.p['KD_ns'].value

		self.p['medium_TGF_beta_amount'].value = self.s['TGF_beta_ex'].concentration * 1e-9 * self.c['Vmed'].size * 6e23

		if self.time <= 0 :
			isConstantValue = self.s['T1R_surf']._constant
			self.s['T1R_surf']._constant = False
			self.s['T1R_surf'].concentration = (self.p['k_T1R'].value * self.p['kdeg_T1R'].value + self.p['k_T1R'].value * self.p['kr'].value) / (self.p['kdeg_T1R'].value * self.p['ki'].value)
			self.s['T1R_surf']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.s['T1R_endo']._constant
			self.s['T1R_endo']._constant = False
			self.s['T1R_endo'].concentration = self.p['k_T1R'].value / self.p['kdeg_T1R'].value
			self.s['T1R_endo']._constant = isConstantValue

		self.p['totalNumT1R'].value = (self.s['T1R_surf'].concentration + self.s['T1R_endo'].concentration + self.s['LRC_surf'].concentration + self.s['LRC_endo'].concentration) * 2.3 * 602

		if self.time <= 0 :
			isConstantValue = self.s['T2R_surf']._constant
			self.s['T2R_surf']._constant = False
			self.s['T2R_surf'].concentration = (self.p['k_T2R'].value * self.p['kdeg_T2R'].value + self.p['k_T2R'].value * self.p['kr'].value) / (self.p['kdeg_T2R'].value * self.p['ki'].value)
			self.s['T2R_surf']._constant = isConstantValue

		if self.time <= 0 :
			isConstantValue = self.s['T2R_endo']._constant
			self.s['T2R_endo']._constant = False
			self.s['T2R_endo'].concentration = self.p['k_T2R'].value / self.p['kdeg_T2R'].value
			self.s['T2R_endo']._constant = isConstantValue

		self.p['totalNumT2R'].value = (self.s['T2R_surf'].concentration + self.s['T2R_endo'].concentration + self.s['LRC_surf'].concentration + self.s['LRC_endo'].concentration) * 2.3 * 602

		if self.time <= 0 :
			isConstantValue = self.p['TGF_beta_dose_mol_per_cell']._constant
			self.p['TGF_beta_dose_mol_per_cell']._constant = False
			self.p['TGF_beta_dose_mol_per_cell'].value = self.s['TGF_beta_ex'].concentration * 1e-9 * self.c['Vmed'].size * 6e23
			self.p['TGF_beta_dose_mol_per_cell']._constant = isConstantValue

		return

	def Reactionre1(self):

		return self.c['Vcyt'].size * self.p['k_T1R'].value

	def Reactionre2(self):

		return self.c['Vcyt'].size * self.p['ki'].value * self.s['T1R_surf'].concentration

	def Reactionre3(self):

		return self.c['Vcyt'].size * self.p['kr'].value * self.s['T1R_endo'].concentration

	def Reactionre4(self):

		return self.c['Vcyt'].size * self.p['kdeg_T1R'].value * self.s['T1R_endo'].concentration

	def Reactionre5(self):

		return self.c['Vcyt'].size * self.p['k_T2R'].value

	def Reactionre6(self):

		return self.c['Vcyt'].size * self.p['ki'].value * self.s['T2R_surf'].concentration

	def Reactionre7(self):

		return self.c['Vcyt'].size * self.p['kr'].value * self.s['T2R_endo'].concentration

	def Reactionre8(self):

		return self.c['Vcyt'].size * self.p['kdeg_T2R'].value * self.s['T2R_endo'].concentration

	def Reactionre9(self):

		return self.c['Vcyt'].size * self.p['ka_LRC'].value * self.s['TGF_beta_ex'].concentration * self.s['T2R_surf'].concentration * self.s['T1R_surf'].concentration

	def Reactionre10(self):

		return self.c['Vcyt'].size * self.p['ki'].value * self.s['LRC_surf'].concentration

	def Reactionre11(self):

		return self.c['Vcyt'].size * self.p['kdeg_LRC'].value * self.s['LRC_endo'].concentration

	def Reactionre12(self):

		return self.c['Vcyt'].size * self.p['kdiss_LRC'].value * self.s['LRC_endo'].concentration

	def Reactionre13(self):

		return self.c['Vcyt'].size * self.p['kdeg_TGF_beta'].value * self.s['TGF_beta_endo'].concentration

	def Reactionre14(self):

		return self.c['Vcyt'].size * self.p['kimp_Smad2'].value * self.s['Smad2c'].concentration

	def Reactionre15(self):

		return self.c['Vnuc'].size * self.p['kexp_Smad2'].value * self.s['Smad2n'].concentration

	def Reactionre16(self):

		return self.c['Vcyt'].size * self.p['kimp_Smad4'].value * self.s['Smad4c'].concentration

	def Reactionre17(self):

		return self.c['Vnuc'].size * self.p['kexp_Smad4'].value * self.s['Smad4n'].concentration

	def Reactionre18(self):

		return self.c['Vcyt'].size * self.p['kpho_Smad2'].value * self.s['Smad2c'].concentration * self.s['LRC_endo'].concentration

	def Reactionre19(self):

		return self.c['Vcyt'].size * self.p['kimp_Smad2'].value * self.s['PSmad2c'].concentration

	def Reactionre20(self):

		return self.c['Vnuc'].size * self.p['kexp_Smad2'].value * self.s['PSmad2n'].concentration

	def Reactionre21(self):

		return self.c['Vcyt'].size * (self.p['kon_Smads'].value * self.s['PSmad2c'].concentration * self.s['Smad4c'].concentration - self.p['koff_Smads'].value * self.s['PSmad2_Smad4_c'].concentration)

	def Reactionre22(self):

		return self.c['Vcyt'].size * self.p['kimp_Smads'].value * self.s['PSmad2_Smad4_c'].concentration

	def Reactionre23(self):

		return self.c['Vnuc'].size * (self.p['koff_Smads'].value * self.s['PSmad2_Smad4_n'].concentration - self.p['kon_Smads'].value * self.s['PSmad2n'].concentration * self.s['Smad4n'].concentration)

	def Reactionre24(self):

		return self.c['Vnuc'].size * self.p['kdepho_Smad2'].value * self.s['PSmad2n'].concentration

	def Reactionre25(self):

		return self.c['Vcyt'].size * (self.p['kon_Smads'].value * self.s['PSmad2c'].concentration**2 - self.p['koff_Smads'].value * self.s['PSmad2_PSmad2_c'].concentration)

	def Reactionre26(self):

		return self.c['Vcyt'].size * self.p['kimp_Smads'].value * self.s['PSmad2_PSmad2_c'].concentration

	def Reactionre27(self):

		return self.c['Vnuc'].size * (self.p['koff_Smads'].value * self.s['PSmad2_PSmad2_n'].concentration - self.p['kon_Smads'].value * self.s['PSmad2n'].concentration**2)

	def Reactionre28(self):

		return self.c['Vcyt'].size * self.p['klid'].value * self.s['LRC_surf'].concentration * self.p['totalNuclearPSmad2'].value

	def Reactionre29(self):

		return self.c['Vmed'].size * (self.p['kon_ns'].value * self.s['TGF_beta_ex'].concentration - self.p['koff_ns'].value * self.s['TGF_beta_ns'].concentration)

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['TGF_beta_ex'].amount, self.s['T1R_surf'].amount, self.s['T1R_endo'].amount, self.s['T2R_surf'].amount, self.s['T2R_endo'].amount, self.s['LRC_surf'].amount, self.s['LRC_endo'].amount, self.s['Smad2c'].amount, self.s['Smad2n'].amount, self.s['Smad4c'].amount, self.s['Smad4n'].amount, self.s['PSmad2c'].amount, self.s['PSmad2_PSmad2_c'].amount, self.s['PSmad2_Smad4_c'].amount, self.s['PSmad2n'].amount, self.s['PSmad2_PSmad2_n'].amount, self.s['PSmad2_Smad4_n'].amount, self.s['TGF_beta_endo'].amount, self.s['TGF_beta_ns'].amount, self.s['AA'].amount, self.s['empty_degraded'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = np.float64)

		stoichiometricMatrix = np.array([[ 0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,-1.],[ 1,-1,1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0.],[ 0,1,-1,-1,0,0,0,0,0,0,0,1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0.],[ 0,0,0,0,1,-1,1,0,-1,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,1,-1,-1,0,0,0,1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,-1,0.],[ 0,0,0,0,0,0,0,0,0,1,-1,-1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,-1.,0,0,0,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0.,0,0,0,0,0,1,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0.,0,0,-1,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0.,0,0,0,0,1,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,-1,1,-1,0,0,0,-2,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,1,-1,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,1,-1,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,1,-1,0,0,1,-1,0,0,2,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,1,-1,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,1,-1,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,1.],[-1,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0.],[ 0,0,0,1,0,0,0,1,0,0,1,0,1,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,1,0.]], dtype = np.float64)

		reactionVelocities = np.array([self.r['re1'](), self.r['re2'](), self.r['re3'](), self.r['re4'](), self.r['re5'](), self.r['re6'](), self.r['re7'](), self.r['re8'](), self.r['re9'](), self.r['re10'](), self.r['re11'](), self.r['re12'](), self.r['re13'](), self.r['re14'](), self.r['re15'](), self.r['re16'](), self.r['re17'](), self.r['re18'](), self.r['re19'](), self.r['re20'](), self.r['re21'](), self.r['re22'](), self.r['re23'](), self.r['re24'](), self.r['re25'](), self.r['re26'](), self.r['re27'](), self.r['re28'](), self.r['re29']()], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['TGF_beta_ex'].amount, self.s['T1R_surf'].amount, self.s['T1R_endo'].amount, self.s['T2R_surf'].amount, self.s['T2R_endo'].amount, self.s['LRC_surf'].amount, self.s['LRC_endo'].amount, self.s['Smad2c'].amount, self.s['Smad2n'].amount, self.s['Smad4c'].amount, self.s['Smad4n'].amount, self.s['PSmad2c'].amount, self.s['PSmad2_PSmad2_c'].amount, self.s['PSmad2_Smad4_c'].amount, self.s['PSmad2n'].amount, self.s['PSmad2_PSmad2_n'].amount, self.s['PSmad2_Smad4_n'].amount, self.s['TGF_beta_endo'].amount, self.s['TGF_beta_ns'].amount, self.s['AA'].amount, self.s['empty_degraded'].amount], dtype = np.float64)
		self.s['TGF_beta_ex'].amount, self.s['T1R_surf'].amount, self.s['T1R_endo'].amount, self.s['T2R_surf'].amount, self.s['T2R_endo'].amount, self.s['LRC_surf'].amount, self.s['LRC_endo'].amount, self.s['Smad2c'].amount, self.s['Smad2n'].amount, self.s['Smad4c'].amount, self.s['Smad4n'].amount, self.s['PSmad2c'].amount, self.s['PSmad2_PSmad2_c'].amount, self.s['PSmad2_Smad4_c'].amount, self.s['PSmad2n'].amount, self.s['PSmad2_PSmad2_n'].amount, self.s['PSmad2_Smad4_n'].amount, self.s['TGF_beta_endo'].amount, self.s['TGF_beta_ns'].amount, self.s['AA'].amount, self.s['empty_degraded'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

class re1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['k_T1R'].value

class re2:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['ki'].value * self.parent.s['T1R_surf'].concentration

class re3:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kr'].value * self.parent.s['T1R_endo'].concentration

class re4:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kdeg_T1R'].value * self.parent.s['T1R_endo'].concentration

class re5:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['k_T2R'].value

class re6:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['ki'].value * self.parent.s['T2R_surf'].concentration

class re7:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kr'].value * self.parent.s['T2R_endo'].concentration

class re8:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kdeg_T2R'].value * self.parent.s['T2R_endo'].concentration

class re9:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['ka_LRC'].value * self.parent.s['TGF_beta_ex'].concentration * self.parent.s['T2R_surf'].concentration * self.parent.s['T1R_surf'].concentration

class re10:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['ki'].value * self.parent.s['LRC_surf'].concentration

class re11:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kdeg_LRC'].value * self.parent.s['LRC_endo'].concentration

class re12:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kdiss_LRC'].value * self.parent.s['LRC_endo'].concentration

class re13:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kdeg_TGF_beta'].value * self.parent.s['TGF_beta_endo'].concentration

class re14:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kimp_Smad2'].value * self.parent.s['Smad2c'].concentration

class re15:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vnuc'].size * self.parent.p['kexp_Smad2'].value * self.parent.s['Smad2n'].concentration

class re16:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kimp_Smad4'].value * self.parent.s['Smad4c'].concentration

class re17:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vnuc'].size * self.parent.p['kexp_Smad4'].value * self.parent.s['Smad4n'].concentration

class re18:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kpho_Smad2'].value * self.parent.s['Smad2c'].concentration * self.parent.s['LRC_endo'].concentration

class re19:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kimp_Smad2'].value * self.parent.s['PSmad2c'].concentration

class re20:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vnuc'].size * self.parent.p['kexp_Smad2'].value * self.parent.s['PSmad2n'].concentration

class re21:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * (self.parent.p['kon_Smads'].value * self.parent.s['PSmad2c'].concentration * self.parent.s['Smad4c'].concentration - self.parent.p['koff_Smads'].value * self.parent.s['PSmad2_Smad4_c'].concentration)

class re22:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kimp_Smads'].value * self.parent.s['PSmad2_Smad4_c'].concentration

class re23:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vnuc'].size * (self.parent.p['koff_Smads'].value * self.parent.s['PSmad2_Smad4_n'].concentration - self.parent.p['kon_Smads'].value * self.parent.s['PSmad2n'].concentration * self.parent.s['Smad4n'].concentration)

class re24:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vnuc'].size * self.parent.p['kdepho_Smad2'].value * self.parent.s['PSmad2n'].concentration

class re25:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * (self.parent.p['kon_Smads'].value * self.parent.s['PSmad2c'].concentration**2 - self.parent.p['koff_Smads'].value * self.parent.s['PSmad2_PSmad2_c'].concentration)

class re26:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['kimp_Smads'].value * self.parent.s['PSmad2_PSmad2_c'].concentration

class re27:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vnuc'].size * (self.parent.p['koff_Smads'].value * self.parent.s['PSmad2_PSmad2_n'].concentration - self.parent.p['kon_Smads'].value * self.parent.s['PSmad2n'].concentration**2)

class re28:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vcyt'].size * self.parent.p['klid'].value * self.parent.s['LRC_surf'].concentration * self.parent.p['totalNuclearPSmad2'].value

class re29:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.c['Vmed'].size * (self.parent.p['kon_ns'].value * self.parent.s['TGF_beta_ex'].concentration - self.parent.p['koff_ns'].value * self.parent.s['TGF_beta_ns'].concentration)

