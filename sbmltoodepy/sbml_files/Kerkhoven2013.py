from sbmltopyode.SBMLModelClasses import *
from scipy.integrate import odeint
import numpy as np
import operator
import math

class Kerkhoven2013:

	def __init__(self):

		self.p = {} #Dictionary of model parameters

		self.c = {} #Dictionary of compartments
		self.c['cytosol'] = Compartment(5.4549, 3, True)
		self.c['glycosome'] = Compartment(0.2451, 3, True)
		self.c['default'] = Compartment(1.0, 3, True)

		self.s = {} #Dictionary of chemical species
		speciesMetadata = SBMLMetadata('_2PGA_c')
		self.s['_2PGA_c'] = Species(0.1, 'Concentration', self.c['cytosol'], False, constant = False)
		speciesMetadata = SBMLMetadata('DHAP_c')
		self.s['DHAP_c'] = Species(2.23132912, 'Concentration', self.c['cytosol'], False, constant = False)
		speciesMetadata = SBMLMetadata('ATP_g')
		self.s['ATP_g'] = Species(0.2405, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('DHAP_g')
		self.s['DHAP_g'] = Species(8.483130623, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('ADP_g')
		self.s['ADP_g'] = Species(1.519, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('Glc6P_g')
		self.s['Glc6P_g'] = Species(0.5, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('ADP_c')
		self.s['ADP_c'] = Species(1.3165, 'Concentration', self.c['cytosol'], False, constant = False)
		speciesMetadata = SBMLMetadata('_3PGA_c')
		self.s['_3PGA_c'] = Species(0.1, 'Concentration', self.c['cytosol'], False, constant = False)
		speciesMetadata = SBMLMetadata('Fru6P_g')
		self.s['Fru6P_g'] = Species(0.5, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('Pi_g')
		self.s['Pi_g'] = Species(0.0, 'Concentration', self.c['glycosome'], False, constant = True)
		speciesMetadata = SBMLMetadata('O2_c')
		self.s['O2_c'] = Species(1.0, 'Concentration', self.c['default'], False, constant = True)
		speciesMetadata = SBMLMetadata('Gly_e')
		self.s['Gly_e'] = Species(0.0, 'Concentration', self.c['default'], False, constant = True)
		speciesMetadata = SBMLMetadata('ATP_c')
		self.s['ATP_c'] = Species(0.3417, 'Concentration', self.c['cytosol'], False, constant = False)
		speciesMetadata = SBMLMetadata('_13BPGA_g')
		self.s['_13BPGA_g'] = Species(0.5, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('Glc_c')
		self.s['Glc_c'] = Species(0.1, 'Concentration', self.c['cytosol'], False, constant = False)
		speciesMetadata = SBMLMetadata('Glc_e')
		self.s['Glc_e'] = Species(5.0, 'Concentration', self.c['default'], False, constant = True)
		speciesMetadata = SBMLMetadata('Glc_g')
		self.s['Glc_g'] = Species(0.1, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('Pyr_c')
		self.s['Pyr_c'] = Species(10.0, 'Concentration', self.c['cytosol'], False, constant = False)
		speciesMetadata = SBMLMetadata('Pyr_e')
		self.s['Pyr_e'] = Species(0.0, 'Concentration', self.c['default'], False, constant = True)
		speciesMetadata = SBMLMetadata('NAD_g')
		self.s['NAD_g'] = Species(2.0, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('Fru16BP_g')
		self.s['Fru16BP_g'] = Species(10.0, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('GA3P_g')
		self.s['GA3P_g'] = Species(2.5, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('CO2_c')
		self.s['CO2_c'] = Species(0.0, 'Concentration', self.c['cytosol'], False, constant = True)
		speciesMetadata = SBMLMetadata('CO2_g')
		self.s['CO2_g'] = Species(0.0, 'Concentration', self.c['glycosome'], False, constant = True)
		speciesMetadata = SBMLMetadata('Gly3P_c')
		self.s['Gly3P_c'] = Species(2.76867088, 'Concentration', self.c['cytosol'], False, constant = False)
		speciesMetadata = SBMLMetadata('Gly3P_g')
		self.s['Gly3P_g'] = Species(10.51686938, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('PEP_c')
		self.s['PEP_c'] = Species(1.0, 'Concentration', self.c['cytosol'], False, constant = False)
		speciesMetadata = SBMLMetadata('AMP_g')
		self.s['AMP_g'] = Species(4.2405, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('_3PGA_g')
		self.s['_3PGA_g'] = Species(0.1, 'Concentration', self.c['glycosome'], False, constant = False)
		speciesMetadata = SBMLMetadata('AMP_c')
		self.s['AMP_c'] = Species(2.2418, 'Concentration', self.c['cytosol'], False, constant = False)
		speciesMetadata = SBMLMetadata('NADH_g')
		self.s['NADH_g'] = Species(2.0, 'Concentration', self.c['glycosome'], False, constant = False)

		self.r = {} #Dictionary of reactiions
		self.r['TPI_g'] = TPI_g(self, SBMLMetadata('TPI_g'))
		self.r['PYK_c'] = PYK_c(self, SBMLMetadata('PYK_c'))
		self.r['PFK_g'] = PFK_g(self, SBMLMetadata('PFK_g'))
		self.r['GlcT_g'] = GlcT_g(self, SBMLMetadata('GlcT_g'))
		self.r['PGAM_c'] = PGAM_c(self, SBMLMetadata('PGAM_c'))
		self.r['PyrT_c'] = PyrT_c(self, SBMLMetadata('PyrT_c'))
		self.r['GlcT_c'] = GlcT_c(self, SBMLMetadata('GlcT_c'))
		self.r['ALD_g'] = ALD_g(self, SBMLMetadata('ALD_g'))
		self.r['ENO_c'] = ENO_c(self, SBMLMetadata('ENO_c'))
		self.r['HXK_g'] = HXK_g(self, SBMLMetadata('HXK_g'))
		self.r['_3PGAT_g'] = _3PGAT_g(self, SBMLMetadata('_3PGAT_g'))
		self.r['PGK_g'] = PGK_g(self, SBMLMetadata('PGK_g'))
		self.r['G3PDH_g'] = G3PDH_g(self, SBMLMetadata('G3PDH_g'))
		self.r['GPO_c'] = GPO_c(self, SBMLMetadata('GPO_c'))
		self.r['ATPu_c'] = ATPu_c(self, SBMLMetadata('ATPu_c'))
		self.r['GK_g'] = GK_g(self, SBMLMetadata('GK_g'))
		self.r['AK_c'] = AK_c(self, SBMLMetadata('AK_c'))
		self.r['PGI_g'] = PGI_g(self, SBMLMetadata('PGI_g'))
		self.r['GAPDH_g'] = GAPDH_g(self, SBMLMetadata('GAPDH_g'))
		self.r['AK_g'] = AK_g(self, SBMLMetadata('AK_g'))
		self.r['GDA_g'] = GDA_g(self, SBMLMetadata('GDA_g'))
		self.time = 0

		self.reactionMetadata = {
		self.ReactionTPI_g: SBMLMetadata('TPI_g'),
		self.ReactionPYK_c: SBMLMetadata('PYK_c'),
		self.ReactionPFK_g: SBMLMetadata('PFK_g'),
		self.ReactionGlcT_g: SBMLMetadata('GlcT_g'),
		self.ReactionPGAM_c: SBMLMetadata('PGAM_c'),
		self.ReactionPyrT_c: SBMLMetadata('PyrT_c'),
		self.ReactionGlcT_c: SBMLMetadata('GlcT_c'),
		self.ReactionALD_g: SBMLMetadata('ALD_g'),
		self.ReactionENO_c: SBMLMetadata('ENO_c'),
		self.ReactionHXK_g: SBMLMetadata('HXK_g'),
		self.Reaction_3PGAT_g: SBMLMetadata('_3PGAT_g'),
		self.ReactionPGK_g: SBMLMetadata('PGK_g'),
		self.ReactionG3PDH_g: SBMLMetadata('G3PDH_g'),
		self.ReactionGPO_c: SBMLMetadata('GPO_c'),
		self.ReactionATPu_c: SBMLMetadata('ATPu_c'),
		self.ReactionGK_g: SBMLMetadata('GK_g'),
		self.ReactionAK_c: SBMLMetadata('AK_c'),
		self.ReactionPGI_g: SBMLMetadata('PGI_g'),
		self.ReactionGAPDH_g: SBMLMetadata('GAPDH_g'),
		self.ReactionAK_g: SBMLMetadata('AK_g'),
		self.ReactionGDA_g: SBMLMetadata('GDA_g')
		}
		self.AssignmentRules()



	def AssignmentRules(self):

		return

	def vAK(self, ADP, AMP, ATP, k1, k2):
		return k1 * ADP**2 - AMP * ATP * k2
	def v2sub2prod(self, Vfmax, Keq, S1, Ks1, S2, Ks2, P1, Kp1, P2, Kp2):
		return Vfmax * S1 * S2 * (1 - P1 * P2 / (Keq * S1 * S2)) / (Ks1 * Ks2 * (1 + S1 / Ks1 + P1 / Kp1) * (1 + S2 / Ks2 + P2 / Kp2))
	def v1sub(self, Vfmax, S, Ks):
		return Vfmax * S / (Ks * (1 + S / Ks))
	def v2sub2prod_compinhib(self, Vfmax, Keq, S1, Ks1, S2, Ks2, P1, Kp1, P2, Kp2, I1, Ki1, I2, Ki2):
		return Vfmax * S1 * S2 * (1 - P1 * P2 / (Keq * S1 * S2)) / (Ks1 * Ks2 * (1 + S1 / Ks1 + P1 / Kp1) * (1 + S2 / Ks2 + P2 / Kp2 + I1 / Ki1 + I2 / Ki2))
	def mass_action_rev(self, k1, S, k2, P):
		return k1 * S - k2 * P
	def v1sub1prod(self, Vfmax, Keq, S, Ks, P, Kp):
		return Vfmax * S * (1 - P / (Keq * S)) / (Ks * (1 + S / Ks + P / Kp))
	def mass_action_irrev(self, k, S):
		return k * S
	def ReactionTPI_g(self):

		TPI_g_Vmax = 999.3
		TPI_g_Keq = 0.046
		TPI_g_KmDHAP = 1.2
		TPI_g_KmGA3P = 0.25
		return self.v1sub1prod(TPI_g_Vmax, TPI_g_Keq, self.s['DHAP_g'].concentration, TPI_g_KmDHAP, self.s['GA3P_g'].concentration, TPI_g_KmGA3P)

	def ReactionPYK_c(self):

		PYK_c_Vmax = 1020.0
		PYK_c_KmPEP = 0.34
		PYK_c_KiATP = 0.57
		PYK_c_KiADP = 0.64
		PYK_c_n = 2.5
		PYK_c_KmADP = 0.114
		PYK_c_Keq = 10800.0
		PYK_c_KmPyr = 50.0
		PYK_c_KmATP = 15.0
		return PYK_c_Vmax * self.s['ADP_c'].concentration * (1 - self.s['Pyr_c'].concentration * self.s['ATP_c'].concentration / (PYK_c_Keq * self.s['PEP_c'].concentration * self.s['ADP_c'].concentration)) * (self.s['PEP_c'].concentration / (PYK_c_KmPEP * (1 + self.s['ADP_c'].concentration / PYK_c_KiADP + self.s['ATP_c'].concentration / PYK_c_KiATP)))**PYK_c_n / (PYK_c_KmADP * (1 + (self.s['PEP_c'].concentration / (PYK_c_KmPEP * (1 + self.s['ADP_c'].concentration / PYK_c_KiADP + self.s['ATP_c'].concentration / PYK_c_KiATP)))**PYK_c_n + self.s['Pyr_c'].concentration / PYK_c_KmPyr) * (1 + self.s['ADP_c'].concentration / PYK_c_KmADP + self.s['ATP_c'].concentration / PYK_c_KmATP))

	def ReactionPFK_g(self):

		PFK_g_Vmax = 1708.0
		PFK_g_Ki1 = 15.8
		PFK_g_KmFru6P = 0.999
		PFK_g_KmATP = 0.0648
		PFK_g_Keq = 1035.0
		PFK_g_KsATP = 0.0393
		PFK_g_KmADP = 1.0
		PFK_g_Ki2 = 10.7
		return PFK_g_Vmax * PFK_g_Ki1 * self.s['Fru6P_g'].concentration * self.s['ATP_g'].concentration * (1 - self.s['Fru16BP_g'].concentration * self.s['ADP_g'].concentration / (PFK_g_Keq * self.s['Fru6P_g'].concentration * self.s['ATP_g'].concentration)) / (PFK_g_KmFru6P * PFK_g_KmATP * (self.s['Fru16BP_g'].concentration + PFK_g_Ki1) * (PFK_g_KsATP / PFK_g_KmATP + self.s['Fru6P_g'].concentration / PFK_g_KmFru6P + self.s['ATP_g'].concentration / PFK_g_KmATP + self.s['ADP_g'].concentration / PFK_g_KmADP + self.s['Fru16BP_g'].concentration * self.s['ADP_g'].concentration / (PFK_g_KmADP * PFK_g_Ki2) + self.s['Fru6P_g'].concentration * self.s['ATP_g'].concentration / (PFK_g_KmFru6P * PFK_g_KmATP)))

	def ReactionGlcT_g(self):

		GlcT_g_k1 = 250000.0
		GlcT_g_k2 = 250000.0
		return self.mass_action_rev(GlcT_g_k1, self.s['Glc_c'].concentration, GlcT_g_k2, self.s['Glc_g'].concentration)

	def ReactionPGAM_c(self):

		PGAM_c_Vmax = 225.0
		PGAM_c_Keq = 0.17
		PGAM_c_Km3PGA = 0.15
		PGAM_c_Km2PGA = 0.16
		return self.v1sub1prod(PGAM_c_Vmax, PGAM_c_Keq, self.s['_3PGA_c'].concentration, PGAM_c_Km3PGA, self.s['_2PGA_c'].concentration, PGAM_c_Km2PGA)

	def ReactionPyrT_c(self):

		PyrT_c_Vmax = 230.0
		PyrT_c_KmPyr = 1.96
		return self.v1sub(PyrT_c_Vmax, self.s['Pyr_c'].concentration, PyrT_c_KmPyr)

	def ReactionGlcT_c(self):

		GlcT_c_Vmax = 111.7
		GlcT_c_KmGlc = 1.0
		GlcT_c_alpha = 0.75
		return GlcT_c_Vmax * (self.s['Glc_e'].concentration - self.s['Glc_c'].concentration) / (GlcT_c_KmGlc + self.s['Glc_e'].concentration + self.s['Glc_c'].concentration + GlcT_c_alpha * self.s['Glc_e'].concentration * self.s['Glc_c'].concentration / GlcT_c_KmGlc)

	def ReactionALD_g(self):

		ALD_g_Vmax = 560.0
		ALD_g_KmFru16BP = 0.009
		ALD_g_KiATP = 0.68
		ALD_g_KiADP = 1.51
		ALD_g_KiAMP = 3.65
		ALD_g_Keq = 0.084
		ALD_g_KmGA3P = 0.067
		ALD_g_KmDHAP = 0.015
		ALD_g_KiGA3P = 0.098
		return ALD_g_Vmax * self.s['Fru16BP_g'].concentration * (1 - self.s['GA3P_g'].concentration * self.s['DHAP_g'].concentration / (self.s['Fru16BP_g'].concentration * ALD_g_Keq)) / (ALD_g_KmFru16BP * (1 + self.s['ATP_g'].concentration / ALD_g_KiATP + self.s['ADP_g'].concentration / ALD_g_KiADP + self.s['AMP_g'].concentration / ALD_g_KiAMP) * (1 + self.s['GA3P_g'].concentration / ALD_g_KmGA3P + self.s['DHAP_g'].concentration / ALD_g_KmDHAP + self.s['Fru16BP_g'].concentration / (ALD_g_KmFru16BP * (1 + self.s['ATP_g'].concentration / ALD_g_KiATP + self.s['ADP_g'].concentration / ALD_g_KiADP + self.s['AMP_g'].concentration / ALD_g_KiAMP)) + self.s['GA3P_g'].concentration * self.s['DHAP_g'].concentration / (ALD_g_KmGA3P * ALD_g_KmDHAP) + self.s['Fru16BP_g'].concentration * self.s['GA3P_g'].concentration / (ALD_g_KmFru16BP * ALD_g_KiGA3P * (1 + self.s['ATP_g'].concentration / ALD_g_KiATP + self.s['ADP_g'].concentration / ALD_g_KiADP + self.s['AMP_g'].concentration / ALD_g_KiAMP))))

	def ReactionENO_c(self):

		ENO_c_Vmax = 598.0
		ENO_c_Keq = 4.17
		ENO_c_Km2PGA = 0.054
		ENO_c_KmPEP = 0.24
		return self.v1sub1prod(ENO_c_Vmax, ENO_c_Keq, self.s['_2PGA_c'].concentration, ENO_c_Km2PGA, self.s['PEP_c'].concentration, ENO_c_KmPEP)

	def ReactionHXK_g(self):

		HXK_g_Vmax = 1774.68
		HXK_g_Keq = 759.0
		HXK_g_KmGlc = 0.1
		HXK_g_KmATP = 0.116
		HXK_g_KmGlc6P = 12.0
		HXK_g_KmADP = 0.126
		return self.v2sub2prod(HXK_g_Vmax, HXK_g_Keq, self.s['Glc_g'].concentration, HXK_g_KmGlc, self.s['ATP_g'].concentration, HXK_g_KmATP, self.s['Glc6P_g'].concentration, HXK_g_KmGlc6P, self.s['ADP_g'].concentration, HXK_g_KmADP)

	def Reaction_3PGAT_g(self):

		_3PGAT_g_k = 250.0
		return self.mass_action_rev(_3PGAT_g_k, self.s['_3PGA_g'].concentration, _3PGAT_g_k, self.s['_3PGA_c'].concentration)

	def ReactionPGK_g(self):

		PGK_g_Vmax = 2862.0
		PGK_g_Keq = 3377.0
		PGK_g_Km13BPGA = 0.003
		PGK_g_KmADP = 0.1
		PGK_g_Km3PGA = 1.62
		PGK_g_KmATP = 0.29
		return self.v2sub2prod(PGK_g_Vmax, PGK_g_Keq, self.s['_13BPGA_g'].concentration, PGK_g_Km13BPGA, self.s['ADP_g'].concentration, PGK_g_KmADP, self.s['_3PGA_g'].concentration, PGK_g_Km3PGA, self.s['ATP_g'].concentration, PGK_g_KmATP)

	def ReactionG3PDH_g(self):

		G3PDH_g_Vmax = 465.0
		G3PDH_g_Keq = 17085.0
		G3PDH_g_KmDHAP = 0.1
		G3PDH_g_KmNADH = 0.01
		G3PDH_g_KmGly3P = 2.0
		G3PDH_g_KmNAD = 0.4
		return self.v2sub2prod(G3PDH_g_Vmax, G3PDH_g_Keq, self.s['DHAP_g'].concentration, G3PDH_g_KmDHAP, self.s['NADH_g'].concentration, G3PDH_g_KmNADH, self.s['Gly3P_g'].concentration, G3PDH_g_KmGly3P, self.s['NAD_g'].concentration, G3PDH_g_KmNAD)

	def ReactionGPO_c(self):

		GPO_c_Vmax = 368.0
		GPO_c_KmGly3P = 1.7
		return self.v1sub(GPO_c_Vmax, self.s['Gly3P_c'].concentration, GPO_c_KmGly3P)

	def ReactionATPu_c(self):

		ATPu_c_k = 50.0
		return ATPu_c_k * self.s['ATP_c'].concentration / self.s['ADP_c'].concentration

	def ReactionGK_g(self):

		GK_g_Vmax = 200.0
		GK_g_Keq = 0.000837
		GK_g_KmGly3P = 3.83
		GK_g_KmADP = 0.56
		GK_g_KmGly = 0.44
		GK_g_KmATP = 0.24
		return self.v2sub2prod(GK_g_Vmax, GK_g_Keq, self.s['Gly3P_g'].concentration, GK_g_KmGly3P, self.s['ADP_g'].concentration, GK_g_KmADP, self.s['Gly_e'].concentration, GK_g_KmGly, self.s['ATP_g'].concentration, GK_g_KmATP)

	def ReactionAK_c(self):

		AK_c_k1 = 480.0
		AK_c_k2 = 1000.0
		return self.vAK(self.s['ADP_c'].concentration, self.s['AMP_c'].concentration, self.s['ATP_c'].concentration, AK_c_k1, AK_c_k2)

	def ReactionPGI_g(self):

		PGI_g_Vmax = 1305.0
		PGI_g_KmGlc6P = 0.4
		PGI_g_Keq = 0.457
		PGI_g_KmFru6P = 0.12
		_6PG_g = 0.0841958
		PGI_g_Ki6PG = 0.14
		return PGI_g_Vmax * self.s['Glc6P_g'].concentration * (1 - self.s['Fru6P_g'].concentration / (PGI_g_Keq * self.s['Glc6P_g'].concentration)) / (PGI_g_KmGlc6P * (1 + self.s['Glc6P_g'].concentration / PGI_g_KmGlc6P + self.s['Fru6P_g'].concentration / PGI_g_KmFru6P + _6PG_g / PGI_g_Ki6PG))

	def ReactionGAPDH_g(self):

		GAPDH_g_Vmax = 720.9
		GAPDH_g_Keq = 0.066
		GAPDH_g_KmGA3P = 0.15
		GAPDH_g_KmNAD = 0.45
		GAPDH_g_Km13BPGA = 0.1
		GAPDH_g_KmNADH = 0.02
		return self.v2sub2prod(GAPDH_g_Vmax, GAPDH_g_Keq, self.s['GA3P_g'].concentration, GAPDH_g_KmGA3P, self.s['NAD_g'].concentration, GAPDH_g_KmNAD, self.s['_13BPGA_g'].concentration, GAPDH_g_Km13BPGA, self.s['NADH_g'].concentration, GAPDH_g_KmNADH)

	def ReactionAK_g(self):

		AK_g_k1 = 480.0
		AK_g_k2 = 1000.0
		return self.vAK(self.s['ADP_g'].concentration, self.s['AMP_g'].concentration, self.s['ATP_g'].concentration, AK_g_k1, AK_g_k2)

	def ReactionGDA_g(self):

		GDA_g_k = 600.0
		return self.s['Gly3P_g'].concentration * GDA_g_k * self.s['DHAP_c'].concentration - self.s['Gly3P_c'].concentration * GDA_g_k * self.s['DHAP_g'].concentration

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['_2PGA_c'].amount, self.s['DHAP_c'].amount, self.s['ATP_g'].amount, self.s['DHAP_g'].amount, self.s['ADP_g'].amount, self.s['Glc6P_g'].amount, self.s['ADP_c'].amount, self.s['_3PGA_c'].amount, self.s['Fru6P_g'].amount, self.s['Pi_g'].amount, self.s['O2_c'].amount, self.s['Gly_e'].amount, self.s['ATP_c'].amount, self.s['_13BPGA_g'].amount, self.s['Glc_c'].amount, self.s['Glc_e'].amount, self.s['Glc_g'].amount, self.s['Pyr_c'].amount, self.s['Pyr_e'].amount, self.s['NAD_g'].amount, self.s['Fru16BP_g'].amount, self.s['GA3P_g'].amount, self.s['CO2_c'].amount, self.s['CO2_g'].amount, self.s['Gly3P_c'].amount, self.s['Gly3P_g'].amount, self.s['PEP_c'].amount, self.s['AMP_g'].amount, self.s['_3PGA_g'].amount, self.s['AMP_c'].amount, self.s['NADH_g'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = np.float64)

		stoichiometricMatrix = np.array([[ 0,0,0,0,1,0,0,0,-1,0,0,0,0,0,0,0,0,0.,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0.,0,0,-1.],[ 0,0,-1,0,0,0,0,0,0,-1,0,1,0,0,0,1,0,0.,0,1,0.],[-1,0,0,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,0.,0,0,1.],[ 0,0,1,0,0,0,0,0,0,1,0,-1,0,0,0,-1,0,0.,0,-2,0.],[ 0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-1.,0,0,0.],[ 0,-1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-2,0.,0,0,0.],[ 0,0,0,0,-1,0,0,0,0,0,1,0,0,0,0,0,0,0.,0,0,0.],[ 0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,-1,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0.,0,0,0.],[ 0,1,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0.,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0.,1,0,0.],[ 0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0.,0,0,0.],[ 0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0.,0,0,0.],[ 0,0,0,1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0.,0,0,0.],[ 0,1,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0.],[ 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0.,-1,0,0.],[ 0,0,1,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0.,0,0,0.],[ 1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0.,-1,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0.,0,0,1.],[ 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0.,0,0,-1.],[ 0,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0.,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,1,0.],[ 0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0.,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0.,0,0,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0.,1,0,0.]], dtype = np.float64)

		reactionVelocities = np.array([self.r['TPI_g'](), self.r['PYK_c'](), self.r['PFK_g'](), self.r['GlcT_g'](), self.r['PGAM_c'](), self.r['PyrT_c'](), self.r['GlcT_c'](), self.r['ALD_g'](), self.r['ENO_c'](), self.r['HXK_g'](), self.r['_3PGAT_g'](), self.r['PGK_g'](), self.r['G3PDH_g'](), self.r['GPO_c'](), self.r['ATPu_c'](), self.r['GK_g'](), self.r['AK_c'](), self.r['PGI_g'](), self.r['GAPDH_g'](), self.r['AK_g'](), self.r['GDA_g']()], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['_2PGA_c'].amount, self.s['DHAP_c'].amount, self.s['ATP_g'].amount, self.s['DHAP_g'].amount, self.s['ADP_g'].amount, self.s['Glc6P_g'].amount, self.s['ADP_c'].amount, self.s['_3PGA_c'].amount, self.s['Fru6P_g'].amount, self.s['Pi_g'].amount, self.s['O2_c'].amount, self.s['Gly_e'].amount, self.s['ATP_c'].amount, self.s['_13BPGA_g'].amount, self.s['Glc_c'].amount, self.s['Glc_e'].amount, self.s['Glc_g'].amount, self.s['Pyr_c'].amount, self.s['Pyr_e'].amount, self.s['NAD_g'].amount, self.s['Fru16BP_g'].amount, self.s['GA3P_g'].amount, self.s['CO2_c'].amount, self.s['CO2_g'].amount, self.s['Gly3P_c'].amount, self.s['Gly3P_g'].amount, self.s['PEP_c'].amount, self.s['AMP_g'].amount, self.s['_3PGA_g'].amount, self.s['AMP_c'].amount, self.s['NADH_g'].amount], dtype = np.float64)
		self.s['_2PGA_c'].amount, self.s['DHAP_c'].amount, self.s['ATP_g'].amount, self.s['DHAP_g'].amount, self.s['ADP_g'].amount, self.s['Glc6P_g'].amount, self.s['ADP_c'].amount, self.s['_3PGA_c'].amount, self.s['Fru6P_g'].amount, self.s['Pi_g'].amount, self.s['O2_c'].amount, self.s['Gly_e'].amount, self.s['ATP_c'].amount, self.s['_13BPGA_g'].amount, self.s['Glc_c'].amount, self.s['Glc_e'].amount, self.s['Glc_g'].amount, self.s['Pyr_c'].amount, self.s['Pyr_e'].amount, self.s['NAD_g'].amount, self.s['Fru16BP_g'].amount, self.s['GA3P_g'].amount, self.s['CO2_c'].amount, self.s['CO2_g'].amount, self.s['Gly3P_c'].amount, self.s['Gly3P_g'].amount, self.s['PEP_c'].amount, self.s['AMP_g'].amount, self.s['_3PGA_g'].amount, self.s['AMP_c'].amount, self.s['NADH_g'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

class TPI_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['TPI_g_Vmax'] = Parameter(999.3, 'TPI_g_Vmax')
		self.p['TPI_g_Keq'] = Parameter(0.046, 'TPI_g_Keq')
		self.p['TPI_g_KmDHAP'] = Parameter(1.2, 'TPI_g_KmDHAP')
		self.p['TPI_g_KmGA3P'] = Parameter(0.25, 'TPI_g_KmGA3P')

	def __call__(self):
		return self.parent.v1sub1prod(self.p['TPI_g_Vmax'].value, self.p['TPI_g_Keq'].value, self.parent.s['DHAP_g'].concentration, self.p['TPI_g_KmDHAP'].value, self.parent.s['GA3P_g'].concentration, self.p['TPI_g_KmGA3P'].value)

class PYK_c:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['PYK_c_Vmax'] = Parameter(1020.0, 'PYK_c_Vmax')
		self.p['PYK_c_KmPEP'] = Parameter(0.34, 'PYK_c_KmPEP')
		self.p['PYK_c_KiATP'] = Parameter(0.57, 'PYK_c_KiATP')
		self.p['PYK_c_KiADP'] = Parameter(0.64, 'PYK_c_KiADP')
		self.p['PYK_c_n'] = Parameter(2.5, 'PYK_c_n')
		self.p['PYK_c_KmADP'] = Parameter(0.114, 'PYK_c_KmADP')
		self.p['PYK_c_Keq'] = Parameter(10800.0, 'PYK_c_Keq')
		self.p['PYK_c_KmPyr'] = Parameter(50.0, 'PYK_c_KmPyr')
		self.p['PYK_c_KmATP'] = Parameter(15.0, 'PYK_c_KmATP')

	def __call__(self):
		return self.p['PYK_c_Vmax'].value * self.parent.s['ADP_c'].concentration * (1 - self.parent.s['Pyr_c'].concentration * self.parent.s['ATP_c'].concentration / (self.p['PYK_c_Keq'].value * self.parent.s['PEP_c'].concentration * self.parent.s['ADP_c'].concentration)) * (self.parent.s['PEP_c'].concentration / (self.p['PYK_c_KmPEP'].value * (1 + self.parent.s['ADP_c'].concentration / self.p['PYK_c_KiADP'].value + self.parent.s['ATP_c'].concentration / self.p['PYK_c_KiATP'].value)))**self.p['PYK_c_n'].value / (self.p['PYK_c_KmADP'].value * (1 + (self.parent.s['PEP_c'].concentration / (self.p['PYK_c_KmPEP'].value * (1 + self.parent.s['ADP_c'].concentration / self.p['PYK_c_KiADP'].value + self.parent.s['ATP_c'].concentration / self.p['PYK_c_KiATP'].value)))**self.p['PYK_c_n'].value + self.parent.s['Pyr_c'].concentration / self.p['PYK_c_KmPyr'].value) * (1 + self.parent.s['ADP_c'].concentration / self.p['PYK_c_KmADP'].value + self.parent.s['ATP_c'].concentration / self.p['PYK_c_KmATP'].value))

class PFK_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['PFK_g_Vmax'] = Parameter(1708.0, 'PFK_g_Vmax')
		self.p['PFK_g_Ki1'] = Parameter(15.8, 'PFK_g_Ki1')
		self.p['PFK_g_KmFru6P'] = Parameter(0.999, 'PFK_g_KmFru6P')
		self.p['PFK_g_KmATP'] = Parameter(0.0648, 'PFK_g_KmATP')
		self.p['PFK_g_Keq'] = Parameter(1035.0, 'PFK_g_Keq')
		self.p['PFK_g_KsATP'] = Parameter(0.0393, 'PFK_g_KsATP')
		self.p['PFK_g_KmADP'] = Parameter(1.0, 'PFK_g_KmADP')
		self.p['PFK_g_Ki2'] = Parameter(10.7, 'PFK_g_Ki2')

	def __call__(self):
		return self.p['PFK_g_Vmax'].value * self.p['PFK_g_Ki1'].value * self.parent.s['Fru6P_g'].concentration * self.parent.s['ATP_g'].concentration * (1 - self.parent.s['Fru16BP_g'].concentration * self.parent.s['ADP_g'].concentration / (self.p['PFK_g_Keq'].value * self.parent.s['Fru6P_g'].concentration * self.parent.s['ATP_g'].concentration)) / (self.p['PFK_g_KmFru6P'].value * self.p['PFK_g_KmATP'].value * (self.parent.s['Fru16BP_g'].concentration + self.p['PFK_g_Ki1'].value) * (self.p['PFK_g_KsATP'].value / self.p['PFK_g_KmATP'].value + self.parent.s['Fru6P_g'].concentration / self.p['PFK_g_KmFru6P'].value + self.parent.s['ATP_g'].concentration / self.p['PFK_g_KmATP'].value + self.parent.s['ADP_g'].concentration / self.p['PFK_g_KmADP'].value + self.parent.s['Fru16BP_g'].concentration * self.parent.s['ADP_g'].concentration / (self.p['PFK_g_KmADP'].value * self.p['PFK_g_Ki2'].value) + self.parent.s['Fru6P_g'].concentration * self.parent.s['ATP_g'].concentration / (self.p['PFK_g_KmFru6P'].value * self.p['PFK_g_KmATP'].value)))

class GlcT_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['GlcT_g_k1'] = Parameter(250000.0, 'GlcT_g_k1')
		self.p['GlcT_g_k2'] = Parameter(250000.0, 'GlcT_g_k2')

	def __call__(self):
		return self.parent.mass_action_rev(self.p['GlcT_g_k1'].value, self.parent.s['Glc_c'].concentration, self.p['GlcT_g_k2'].value, self.parent.s['Glc_g'].concentration)

class PGAM_c:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['PGAM_c_Vmax'] = Parameter(225.0, 'PGAM_c_Vmax')
		self.p['PGAM_c_Keq'] = Parameter(0.17, 'PGAM_c_Keq')
		self.p['PGAM_c_Km3PGA'] = Parameter(0.15, 'PGAM_c_Km3PGA')
		self.p['PGAM_c_Km2PGA'] = Parameter(0.16, 'PGAM_c_Km2PGA')

	def __call__(self):
		return self.parent.v1sub1prod(self.p['PGAM_c_Vmax'].value, self.p['PGAM_c_Keq'].value, self.parent.s['_3PGA_c'].concentration, self.p['PGAM_c_Km3PGA'].value, self.parent.s['_2PGA_c'].concentration, self.p['PGAM_c_Km2PGA'].value)

class PyrT_c:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['PyrT_c_Vmax'] = Parameter(230.0, 'PyrT_c_Vmax')
		self.p['PyrT_c_KmPyr'] = Parameter(1.96, 'PyrT_c_KmPyr')

	def __call__(self):
		return self.parent.v1sub(self.p['PyrT_c_Vmax'].value, self.parent.s['Pyr_c'].concentration, self.p['PyrT_c_KmPyr'].value)

class GlcT_c:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['GlcT_c_Vmax'] = Parameter(111.7, 'GlcT_c_Vmax')
		self.p['GlcT_c_KmGlc'] = Parameter(1.0, 'GlcT_c_KmGlc')
		self.p['GlcT_c_alpha'] = Parameter(0.75, 'GlcT_c_alpha')

	def __call__(self):
		return self.p['GlcT_c_Vmax'].value * (self.parent.s['Glc_e'].concentration - self.parent.s['Glc_c'].concentration) / (self.p['GlcT_c_KmGlc'].value + self.parent.s['Glc_e'].concentration + self.parent.s['Glc_c'].concentration + self.p['GlcT_c_alpha'].value * self.parent.s['Glc_e'].concentration * self.parent.s['Glc_c'].concentration / self.p['GlcT_c_KmGlc'].value)

class ALD_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['ALD_g_Vmax'] = Parameter(560.0, 'ALD_g_Vmax')
		self.p['ALD_g_KmFru16BP'] = Parameter(0.009, 'ALD_g_KmFru16BP')
		self.p['ALD_g_KiATP'] = Parameter(0.68, 'ALD_g_KiATP')
		self.p['ALD_g_KiADP'] = Parameter(1.51, 'ALD_g_KiADP')
		self.p['ALD_g_KiAMP'] = Parameter(3.65, 'ALD_g_KiAMP')
		self.p['ALD_g_Keq'] = Parameter(0.084, 'ALD_g_Keq')
		self.p['ALD_g_KmGA3P'] = Parameter(0.067, 'ALD_g_KmGA3P')
		self.p['ALD_g_KmDHAP'] = Parameter(0.015, 'ALD_g_KmDHAP')
		self.p['ALD_g_KiGA3P'] = Parameter(0.098, 'ALD_g_KiGA3P')

	def __call__(self):
		return self.p['ALD_g_Vmax'].value * self.parent.s['Fru16BP_g'].concentration * (1 - self.parent.s['GA3P_g'].concentration * self.parent.s['DHAP_g'].concentration / (self.parent.s['Fru16BP_g'].concentration * self.p['ALD_g_Keq'].value)) / (self.p['ALD_g_KmFru16BP'].value * (1 + self.parent.s['ATP_g'].concentration / self.p['ALD_g_KiATP'].value + self.parent.s['ADP_g'].concentration / self.p['ALD_g_KiADP'].value + self.parent.s['AMP_g'].concentration / self.p['ALD_g_KiAMP'].value) * (1 + self.parent.s['GA3P_g'].concentration / self.p['ALD_g_KmGA3P'].value + self.parent.s['DHAP_g'].concentration / self.p['ALD_g_KmDHAP'].value + self.parent.s['Fru16BP_g'].concentration / (self.p['ALD_g_KmFru16BP'].value * (1 + self.parent.s['ATP_g'].concentration / self.p['ALD_g_KiATP'].value + self.parent.s['ADP_g'].concentration / self.p['ALD_g_KiADP'].value + self.parent.s['AMP_g'].concentration / self.p['ALD_g_KiAMP'].value)) + self.parent.s['GA3P_g'].concentration * self.parent.s['DHAP_g'].concentration / (self.p['ALD_g_KmGA3P'].value * self.p['ALD_g_KmDHAP'].value) + self.parent.s['Fru16BP_g'].concentration * self.parent.s['GA3P_g'].concentration / (self.p['ALD_g_KmFru16BP'].value * self.p['ALD_g_KiGA3P'].value * (1 + self.parent.s['ATP_g'].concentration / self.p['ALD_g_KiATP'].value + self.parent.s['ADP_g'].concentration / self.p['ALD_g_KiADP'].value + self.parent.s['AMP_g'].concentration / self.p['ALD_g_KiAMP'].value))))

class ENO_c:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['ENO_c_Vmax'] = Parameter(598.0, 'ENO_c_Vmax')
		self.p['ENO_c_Keq'] = Parameter(4.17, 'ENO_c_Keq')
		self.p['ENO_c_Km2PGA'] = Parameter(0.054, 'ENO_c_Km2PGA')
		self.p['ENO_c_KmPEP'] = Parameter(0.24, 'ENO_c_KmPEP')

	def __call__(self):
		return self.parent.v1sub1prod(self.p['ENO_c_Vmax'].value, self.p['ENO_c_Keq'].value, self.parent.s['_2PGA_c'].concentration, self.p['ENO_c_Km2PGA'].value, self.parent.s['PEP_c'].concentration, self.p['ENO_c_KmPEP'].value)

class HXK_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['HXK_g_Vmax'] = Parameter(1774.68, 'HXK_g_Vmax')
		self.p['HXK_g_Keq'] = Parameter(759.0, 'HXK_g_Keq')
		self.p['HXK_g_KmGlc'] = Parameter(0.1, 'HXK_g_KmGlc')
		self.p['HXK_g_KmATP'] = Parameter(0.116, 'HXK_g_KmATP')
		self.p['HXK_g_KmGlc6P'] = Parameter(12.0, 'HXK_g_KmGlc6P')
		self.p['HXK_g_KmADP'] = Parameter(0.126, 'HXK_g_KmADP')

	def __call__(self):
		return self.parent.v2sub2prod(self.p['HXK_g_Vmax'].value, self.p['HXK_g_Keq'].value, self.parent.s['Glc_g'].concentration, self.p['HXK_g_KmGlc'].value, self.parent.s['ATP_g'].concentration, self.p['HXK_g_KmATP'].value, self.parent.s['Glc6P_g'].concentration, self.p['HXK_g_KmGlc6P'].value, self.parent.s['ADP_g'].concentration, self.p['HXK_g_KmADP'].value)

class _3PGAT_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['_3PGAT_g_k'] = Parameter(250.0, '_3PGAT_g_k')

	def __call__(self):
		return self.parent.mass_action_rev(self.p['_3PGAT_g_k'].value, self.parent.s['_3PGA_g'].concentration, self.p['_3PGAT_g_k'].value, self.parent.s['_3PGA_c'].concentration)

class PGK_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['PGK_g_Vmax'] = Parameter(2862.0, 'PGK_g_Vmax')
		self.p['PGK_g_Keq'] = Parameter(3377.0, 'PGK_g_Keq')
		self.p['PGK_g_Km13BPGA'] = Parameter(0.003, 'PGK_g_Km13BPGA')
		self.p['PGK_g_KmADP'] = Parameter(0.1, 'PGK_g_KmADP')
		self.p['PGK_g_Km3PGA'] = Parameter(1.62, 'PGK_g_Km3PGA')
		self.p['PGK_g_KmATP'] = Parameter(0.29, 'PGK_g_KmATP')

	def __call__(self):
		return self.parent.v2sub2prod(self.p['PGK_g_Vmax'].value, self.p['PGK_g_Keq'].value, self.parent.s['_13BPGA_g'].concentration, self.p['PGK_g_Km13BPGA'].value, self.parent.s['ADP_g'].concentration, self.p['PGK_g_KmADP'].value, self.parent.s['_3PGA_g'].concentration, self.p['PGK_g_Km3PGA'].value, self.parent.s['ATP_g'].concentration, self.p['PGK_g_KmATP'].value)

class G3PDH_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['G3PDH_g_Vmax'] = Parameter(465.0, 'G3PDH_g_Vmax')
		self.p['G3PDH_g_Keq'] = Parameter(17085.0, 'G3PDH_g_Keq')
		self.p['G3PDH_g_KmDHAP'] = Parameter(0.1, 'G3PDH_g_KmDHAP')
		self.p['G3PDH_g_KmNADH'] = Parameter(0.01, 'G3PDH_g_KmNADH')
		self.p['G3PDH_g_KmGly3P'] = Parameter(2.0, 'G3PDH_g_KmGly3P')
		self.p['G3PDH_g_KmNAD'] = Parameter(0.4, 'G3PDH_g_KmNAD')

	def __call__(self):
		return self.parent.v2sub2prod(self.p['G3PDH_g_Vmax'].value, self.p['G3PDH_g_Keq'].value, self.parent.s['DHAP_g'].concentration, self.p['G3PDH_g_KmDHAP'].value, self.parent.s['NADH_g'].concentration, self.p['G3PDH_g_KmNADH'].value, self.parent.s['Gly3P_g'].concentration, self.p['G3PDH_g_KmGly3P'].value, self.parent.s['NAD_g'].concentration, self.p['G3PDH_g_KmNAD'].value)

class GPO_c:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['GPO_c_Vmax'] = Parameter(368.0, 'GPO_c_Vmax')
		self.p['GPO_c_KmGly3P'] = Parameter(1.7, 'GPO_c_KmGly3P')

	def __call__(self):
		return self.parent.v1sub(self.p['GPO_c_Vmax'].value, self.parent.s['Gly3P_c'].concentration, self.p['GPO_c_KmGly3P'].value)

class ATPu_c:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['ATPu_c_k'] = Parameter(50.0, 'ATPu_c_k')

	def __call__(self):
		return self.p['ATPu_c_k'].value * self.parent.s['ATP_c'].concentration / self.parent.s['ADP_c'].concentration

class GK_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['GK_g_Vmax'] = Parameter(200.0, 'GK_g_Vmax')
		self.p['GK_g_Keq'] = Parameter(0.000837, 'GK_g_Keq')
		self.p['GK_g_KmGly3P'] = Parameter(3.83, 'GK_g_KmGly3P')
		self.p['GK_g_KmADP'] = Parameter(0.56, 'GK_g_KmADP')
		self.p['GK_g_KmGly'] = Parameter(0.44, 'GK_g_KmGly')
		self.p['GK_g_KmATP'] = Parameter(0.24, 'GK_g_KmATP')

	def __call__(self):
		return self.parent.v2sub2prod(self.p['GK_g_Vmax'].value, self.p['GK_g_Keq'].value, self.parent.s['Gly3P_g'].concentration, self.p['GK_g_KmGly3P'].value, self.parent.s['ADP_g'].concentration, self.p['GK_g_KmADP'].value, self.parent.s['Gly_e'].concentration, self.p['GK_g_KmGly'].value, self.parent.s['ATP_g'].concentration, self.p['GK_g_KmATP'].value)

class AK_c:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['AK_c_k1'] = Parameter(480.0, 'AK_c_k1')
		self.p['AK_c_k2'] = Parameter(1000.0, 'AK_c_k2')

	def __call__(self):
		return self.parent.vAK(self.parent.s['ADP_c'].concentration, self.parent.s['AMP_c'].concentration, self.parent.s['ATP_c'].concentration, self.p['AK_c_k1'].value, self.p['AK_c_k2'].value)

class PGI_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['PGI_g_Vmax'] = Parameter(1305.0, 'PGI_g_Vmax')
		self.p['PGI_g_KmGlc6P'] = Parameter(0.4, 'PGI_g_KmGlc6P')
		self.p['PGI_g_Keq'] = Parameter(0.457, 'PGI_g_Keq')
		self.p['PGI_g_KmFru6P'] = Parameter(0.12, 'PGI_g_KmFru6P')
		self.p['_6PG_g'] = Parameter(0.0841958, '_6PG_g')
		self.p['PGI_g_Ki6PG'] = Parameter(0.14, 'PGI_g_Ki6PG')

	def __call__(self):
		return self.p['PGI_g_Vmax'].value * self.parent.s['Glc6P_g'].concentration * (1 - self.parent.s['Fru6P_g'].concentration / (self.p['PGI_g_Keq'].value * self.parent.s['Glc6P_g'].concentration)) / (self.p['PGI_g_KmGlc6P'].value * (1 + self.parent.s['Glc6P_g'].concentration / self.p['PGI_g_KmGlc6P'].value + self.parent.s['Fru6P_g'].concentration / self.p['PGI_g_KmFru6P'].value + self.p['_6PG_g'].value / self.p['PGI_g_Ki6PG'].value))

class GAPDH_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['GAPDH_g_Vmax'] = Parameter(720.9, 'GAPDH_g_Vmax')
		self.p['GAPDH_g_Keq'] = Parameter(0.066, 'GAPDH_g_Keq')
		self.p['GAPDH_g_KmGA3P'] = Parameter(0.15, 'GAPDH_g_KmGA3P')
		self.p['GAPDH_g_KmNAD'] = Parameter(0.45, 'GAPDH_g_KmNAD')
		self.p['GAPDH_g_Km13BPGA'] = Parameter(0.1, 'GAPDH_g_Km13BPGA')
		self.p['GAPDH_g_KmNADH'] = Parameter(0.02, 'GAPDH_g_KmNADH')

	def __call__(self):
		return self.parent.v2sub2prod(self.p['GAPDH_g_Vmax'].value, self.p['GAPDH_g_Keq'].value, self.parent.s['GA3P_g'].concentration, self.p['GAPDH_g_KmGA3P'].value, self.parent.s['NAD_g'].concentration, self.p['GAPDH_g_KmNAD'].value, self.parent.s['_13BPGA_g'].concentration, self.p['GAPDH_g_Km13BPGA'].value, self.parent.s['NADH_g'].concentration, self.p['GAPDH_g_KmNADH'].value)

class AK_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['AK_g_k1'] = Parameter(480.0, 'AK_g_k1')
		self.p['AK_g_k2'] = Parameter(1000.0, 'AK_g_k2')

	def __call__(self):
		return self.parent.vAK(self.parent.s['ADP_g'].concentration, self.parent.s['AMP_g'].concentration, self.parent.s['ATP_g'].concentration, self.p['AK_g_k1'].value, self.p['AK_g_k2'].value)

class GDA_g:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata

		self.p['GDA_g_k'] = Parameter(600.0, 'GDA_g_k')

	def __call__(self):
		return self.parent.s['Gly3P_g'].concentration * self.p['GDA_g_k'].value * self.parent.s['DHAP_c'].concentration - self.parent.s['Gly3P_c'].concentration * self.p['GDA_g_k'].value * self.parent.s['DHAP_g'].concentration

class vAK:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.metadata = metadata

	def __call__(self, ADP, AMP, ATP, k1, k2):
		return k1 * ADP**2 - AMP * ATP * k2

class v2sub2prod:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.metadata = metadata

	def __call__(self, Vfmax, Keq, S1, Ks1, S2, Ks2, P1, Kp1, P2, Kp2):
		return Vfmax * S1 * S2 * (1 - P1 * P2 / (Keq * S1 * S2)) / (Ks1 * Ks2 * (1 + S1 / Ks1 + P1 / Kp1) * (1 + S2 / Ks2 + P2 / Kp2))

class v1sub:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.metadata = metadata

	def __call__(self, Vfmax, S, Ks):
		return Vfmax * S / (Ks * (1 + S / Ks))

class v2sub2prod_compinhib:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.metadata = metadata

	def __call__(self, Vfmax, Keq, S1, Ks1, S2, Ks2, P1, Kp1, P2, Kp2, I1, Ki1, I2, Ki2):
		return Vfmax * S1 * S2 * (1 - P1 * P2 / (Keq * S1 * S2)) / (Ks1 * Ks2 * (1 + S1 / Ks1 + P1 / Kp1) * (1 + S2 / Ks2 + P2 / Kp2 + I1 / Ki1 + I2 / Ki2))

class mass_action_rev:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.metadata = metadata

	def __call__(self, k1, S, k2, P):
		return k1 * S - k2 * P

class v1sub1prod:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.metadata = metadata

	def __call__(self, Vfmax, Keq, S, Ks, P, Kp):
		return Vfmax * S * (1 - P / (Keq * S)) / (Ks * (1 + S / Ks + P / Kp))

class mass_action_irrev:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.metadata = metadata

	def __call__(self, k, S):
		return k * S

