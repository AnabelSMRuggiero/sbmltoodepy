from sbmltopyode.SBMLModelClasses import *
from scipy.integrate import odeint
import numpy as np
import operator
import math

class Borisov2009:

	def __init__(self):

		self.p = {} #Dictionary of model parameters
		self.p['EGF_tot'] = Parameter(3.4, 'EGF_tot', False)
		self.p['k1'] = Parameter(0.068, 'k1', True)
		self.p['Kd1'] = Parameter(0.58824, 'Kd1', True)
		self.p['k2'] = Parameter(0.033, 'k2', True)
		self.p['Kd2'] = Parameter(15.0, 'Kd2', True)
		self.p['k3'] = Parameter(0.4, 'k3', True)
		self.p['k4'] = Parameter(0.000666, 'k4', True)
		self.p['Kd4'] = Parameter(10.0, 'Kd4', True)
		self.p['Kd5'] = Parameter(10.0, 'Kd5', True)
		self.p['k5'] = Parameter(0.0133, 'k5', True)
		self.p['k6'] = Parameter(0.333, 'k6', True)
		self.p['k7'] = Parameter(0.000666, 'k7', True)
		self.p['Kd7'] = Parameter(400.0, 'Kd7', True)
		self.p['V8'] = Parameter(200.0, 'V8', True)
		self.p['Km8'] = Parameter(100.0, 'Km8', True)
		self.p['Kd9'] = Parameter(10.0, 'Kd9', True)
		self.p['k9'] = Parameter(0.00666, 'k9', True)
		self.p['Kd10'] = Parameter(400.0, 'Kd10', True)
		self.p['k10'] = Parameter(0.0004, 'k10', True)
		self.p['k11'] = Parameter(0.0, 'k11', False)
		self.p['k12'] = Parameter(0.00933, 'k12', True)
		self.p['Kd12'] = Parameter(12.45, 'Kd12', True)
		self.p['k13'] = Parameter(6.66e-06, 'k13', True)
		self.p['Kd13'] = Parameter(200.0, 'Kd13', True)
		self.p['k17'] = Parameter(0.000185, 'k17', True)
		self.p['k24'] = Parameter(0.011322, 'k24', True)
		self.p['Kd24'] = Parameter(0.029412, 'Kd24', True)
		self.p['k25'] = Parameter(1.66, 'k25', True)
		self.p['k26'] = Parameter(0.00933, 'k26', True)
		self.p['Kd26'] = Parameter(124.5, 'Kd26', True)
		self.p['k27'] = Parameter(6.66e-08, 'k27', True)
		self.p['Kd27'] = Parameter(2000000.0, 'Kd27', True)
		self.p['k28'] = Parameter(0.1066, 'k28', True)
		self.p['Kd28'] = Parameter(3.75, 'Kd28', True)
		self.p['k29'] = Parameter(0.66, 'k29', True)
		self.p['k30'] = Parameter(0.0066, 'k30', True)
		self.p['Kd30'] = Parameter(10.0, 'Kd30', True)
		self.p['V31'] = Parameter(333.0, 'V31', True)
		self.p['Km31'] = Parameter(143.3, 'Km31', True)
		self.p['kcat40'] = Parameter(6.6, 'kcat40', True)
		self.p['Km40'] = Parameter(110.0, 'Km40', True)
		self.p['alpha40'] = Parameter(0.00025, 'alpha40', True)
		self.p['V41'] = Parameter(6.66, 'V41', True)
		self.p['Km41'] = Parameter(50.0, 'Km41', True)
		self.p['k42'] = Parameter(0.00666, 'k42', True)
		self.p['Kd42'] = Parameter(10.0, 'Kd42', True)
		self.p['kcat43'] = Parameter(33.3, 'kcat43', True)
		self.p['Km43'] = Parameter(150.0, 'Km43', True)
		self.p['alpha43'] = Parameter(0.05, 'alpha43', True)
		self.p['Kd45'] = Parameter(100000.0, 'Kd45', True)
		self.p['k45'] = Parameter(0.000666, 'k45', True)
		self.p['k46'] = Parameter(0.00666, 'k46', True)
		self.p['Kd46'] = Parameter(1.0, 'Kd46', True)
		self.p['k47'] = Parameter(0.000666, 'k47', True)
		self.p['Kd47'] = Parameter(1000.0, 'Kd47', True)
		self.p['k48'] = Parameter(0.666, 'k48', True)
		self.p['k49'] = Parameter(0.000666, 'k49', True)
		self.p['Kd49'] = Parameter(1.0, 'Kd49', True)
		self.p['kcat50'] = Parameter(3333.0, 'kcat50', True)
		self.p['alpha50'] = Parameter(0.0001, 'alpha50', True)
		self.p['Km50'] = Parameter(150.0, 'Km50', True)
		self.p['V51'] = Parameter(333.0, 'V51', True)
		self.p['Km51'] = Parameter(130.0, 'Km51', True)
		self.p['k52'] = Parameter(0.002, 'k52', True)
		self.p['Kd52'] = Parameter(1.0, 'Kd52', True)
		self.p['k53'] = Parameter(0.0133, 'k53', True)
		self.p['Kd53'] = Parameter(2.5, 'Kd53', True)
		self.p['k54'] = Parameter(1e-05, 'k54', True)
		self.p['Kd54'] = Parameter(66666.0, 'Kd54', True)
		self.p['k55'] = Parameter(0.000666, 'k55', True)
		self.p['Kd55'] = Parameter(100.0, 'Kd55', True)
		self.p['k56'] = Parameter(0.666, 'k56', True)
		self.p['kcat57'] = Parameter(0.133, 'kcat57', True)
		self.p['Km57'] = Parameter(150.0, 'Km57', True)
		self.p['V58'] = Parameter(2.0, 'V58', True)
		self.p['Km58'] = Parameter(130.0, 'Km58', True)
		self.p['k59'] = Parameter(0.01, 'k59', True)
		self.p['Kd59'] = Parameter(20.0, 'Kd59', True)
		self.p['k60'] = Parameter(4.66, 'k60', True)
		self.p['k61'] = Parameter(3.33, 'k61', True)
		self.p['kcat62'] = Parameter(5.33, 'kcat62', True)
		self.p['Km62'] = Parameter(50.0, 'Km62', True)
		self.p['kcat63'] = Parameter(20000.0, 'kcat63', True)
		self.p['Km63'] = Parameter(50.0, 'Km63', True)
		self.p['k64'] = Parameter(0.0, 'k64', True)
		self.p['k_64'] = Parameter(2.5, 'k_64', True)
		self.p['kcat65'] = Parameter(0.1, 'kcat65', True)
		self.p['Km65'] = Parameter(400.0, 'Km65', True)
		self.p['kcat66'] = Parameter(3.33, 'kcat66', True)
		self.p['Km66'] = Parameter(10.0, 'Km66', True)
		self.p['kcat67'] = Parameter(0.666, 'kcat67', True)
		self.p['Km67'] = Parameter(10000.0, 'Km67', True)
		self.p['alpha67'] = Parameter(1e-06, 'alpha67', True)
		self.p['beta67'] = Parameter(2.0, 'beta67', True)
		self.p['kcat68'] = Parameter(0.133, 'kcat68', True)
		self.p['Km68'] = Parameter(50.0, 'Km68', True)
		self.p['V69'] = Parameter(16.6, 'V69', True)
		self.p['Km69'] = Parameter(675.299, 'Km69', True)
		self.p['kcat70'] = Parameter(0.333, 'kcat70', True)
		self.p['Km70'] = Parameter(500.0, 'Km70', True)
		self.p['kcat71'] = Parameter(0.666, 'kcat71', True)
		self.p['Km71'] = Parameter(500.0, 'Km71', True)
		self.p['V72'] = Parameter(33.3, 'V72', True)
		self.p['Km72'] = Parameter(500.0, 'Km72', True)
		self.p['V73'] = Parameter(23.33, 'V73', True)
		self.p['Km73'] = Parameter(500.0, 'Km73', True)
		self.p['k74'] = Parameter(0.00666, 'k74', True)
		self.p['Kd74'] = Parameter(100.0, 'Kd74', True)
		self.p['kcat75'] = Parameter(4.66, 'kcat75', True)
		self.p['Km75'] = Parameter(500.0, 'Km75', True)
		self.p['V76'] = Parameter(16.66, 'V76', True)
		self.p['Km76'] = Parameter(1.0, 'Km76', True)
		self.p['kcat77'] = Parameter(0.666, 'kcat77', True)
		self.p['alpha77'] = Parameter(0.5, 'alpha77', True)
		self.p['Km77'] = Parameter(100.0, 'Km77', True)
		self.p['k_77'] = Parameter(0.666, 'k_77', True)
		self.p['kcat78'] = Parameter(0.666, 'kcat78', True)
		self.p['Km78'] = Parameter(100.0, 'Km78', True)
		self.p['k_78'] = Parameter(0.666, 'k_78', True)
		self.p['kcat79'] = Parameter(0.0466, 'kcat79', True)
		self.p['Km79'] = Parameter(5000.0, 'Km79', True)
		self.p['k_79'] = Parameter(6.66e-05, 'k_79', True)
		self.p['kcat80'] = Parameter(0.04, 'kcat80', True)
		self.p['Km80'] = Parameter(700.0, 'Km80', True)
		self.p['k_80'] = Parameter(6.66e-05, 'k_80', True)
		self.p['kcat81'] = Parameter(0.166, 'kcat81', True)
		self.p['Km81'] = Parameter(300.0, 'Km81', True)
		self.p['k_81'] = Parameter(6.66e-05, 'k_81', True)
		self.p['V_82'] = Parameter(133.0, 'V_82', True)
		self.p['Km82'] = Parameter(50.0, 'Km82', True)
		self.p['k83'] = Parameter(0.0166, 'k83', True)
		self.p['V_84'] = Parameter(333.0, 'V_84', True)
		self.p['Km84'] = Parameter(266.0, 'Km84', True)
		self.p['k85'] = Parameter(0.0166, 'k85', True)
		self.p['k111'] = Parameter(0.0133, 'k111', True)
		self.p['k118'] = Parameter(0.001, 'k118', True)
		self.p['k_1'] = Parameter(0.0, 'k_1', False)
		self.p['k_2'] = Parameter(0.0, 'k_2', False)
		self.p['k_4'] = Parameter(0.0, 'k_4', False)
		self.p['k_5'] = Parameter(0.0, 'k_5', False)
		self.p['k_7'] = Parameter(0.0, 'k_7', False)
		self.p['k_9'] = Parameter(0.0, 'k_9', False)
		self.p['k_10'] = Parameter(0.0, 'k_10', False)
		self.p['k_11'] = Parameter(0.0, 'k_11', False)
		self.p['k_12'] = Parameter(0.0, 'k_12', False)
		self.p['k_13'] = Parameter(0.0, 'k_13', False)
		self.p['k_24'] = Parameter(0.0, 'k_24', False)
		self.p['k_26'] = Parameter(0.0, 'k_26', False)
		self.p['k_27'] = Parameter(0.0, 'k_27', False)
		self.p['k_28'] = Parameter(0.0, 'k_28', False)
		self.p['k_30'] = Parameter(0.0, 'k_30', False)
		self.p['k_42'] = Parameter(0.0, 'k_42', False)
		self.p['k_45'] = Parameter(0.0, 'k_45', False)
		self.p['k_46'] = Parameter(0.0, 'k_46', False)
		self.p['k_47'] = Parameter(0.0, 'k_47', False)
		self.p['k_49'] = Parameter(0.0, 'k_49', False)
		self.p['k_52'] = Parameter(0.0, 'k_52', False)
		self.p['k_53'] = Parameter(0.0, 'k_53', False)
		self.p['k_54'] = Parameter(0.0, 'k_54', False)
		self.p['k_55'] = Parameter(0.0, 'k_55', False)
		self.p['k_59'] = Parameter(0.0, 'k_59', False)
		self.p['k_74'] = Parameter(0.0, 'k_74', False)

		self.c = {} #Dictionary of compartments
		self.c['cell'] = Compartment(1.0, 3, True)
		self.c['extra'] = Compartment(34.0, 3, True)

		self.s = {} #Dictionary of chemical species
		speciesMetadata = SBMLMetadata('')
		self.s['EGF'] = Species(1.0, 'Concentration', self.c['extra'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['I'] = Species(0.0, 'Concentration', self.c['extra'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['RE'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Rd'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Rp'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['GS'] = Species(200.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Rp_GS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Shc'] = Species(270.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Rp_Shc'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Rp_pShc'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['pShc'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Rp_pShc_GS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['PI3K'] = Species(200.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Rp_PI3K'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['RasGAP'] = Species(50.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Rp_RasGAP'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IRL'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IRp'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IRp_PI3K'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IRp_RasGAP'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IRS'] = Species(300.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IRp_IRS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IRp_IRSp'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IRSp'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['iSrc'] = Species(518.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mIRS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mIRSp'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mIRSp_GS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mIRSp_PI3K'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['SHP2'] = Species(300.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mIRSp_SHP2'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['GAB'] = Species(225.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mGAB'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mGABp'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mGABp_GS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mGABp_PI3K'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mGABp_SHP2'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mGABp_pSHP2'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['PIP3'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['dRas'] = Species(150.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Raf'] = Species(100.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['aRaf'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Mek'] = Species(200.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Erk'] = Species(400.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['pErk'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['PDK1'] = Species(100.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Akt'] = Species(100.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['pAkt'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mTOR'] = Species(100.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Null'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['aaRaf'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['PKA'] = Species(100.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['pShc_GS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['ppMek'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mGABp_pSHP2_GS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['R'] = Species(100.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['ppErk'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IR'] = Species(150.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mPDK1'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['tRas'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['tRas_PI3K'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['ppAkt'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['mGABp_RasGAP'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['amTOR'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['iGS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['imGAB'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['imIRS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['aSrc'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['Ri'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IRi'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['iPX'] = Species(200.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['aPX'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['aPX_GS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IRSp_PI3K'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IRSp_GS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['IRSp_SHP2'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['GABp'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['GABp_PI3K'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['GABp_GS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['GABp_RasGAP'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['GABp_SHP2'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['GABp_pSHP2'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['GABp_pSHP2_GS'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['imGABp'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['bRasGAP'] = Species(1e-05, 'Concentration', self.c['cell'], False, constant = False)
		speciesMetadata = SBMLMetadata('')
		self.s['phosphorylated_Akt'] = Species(0.0, 'Concentration', self.c['cell'], False, constant = False)
		self.s['phosphorylated_Akt']._modifiedBy = 1

		self.r = {} #Dictionary of reactiions
		self.r['reaction_1'] = reaction_1(self, SBMLMetadata(''))
		self.r['reaction_2'] = reaction_2(self, SBMLMetadata(''))
		self.r['reaction_3'] = reaction_3(self, SBMLMetadata(''))
		self.r['reaction_4'] = reaction_4(self, SBMLMetadata(''))
		self.r['reaction_5'] = reaction_5(self, SBMLMetadata(''))
		self.r['reaction_6'] = reaction_6(self, SBMLMetadata(''))
		self.r['reaction_7'] = reaction_7(self, SBMLMetadata(''))
		self.r['reaction_8'] = reaction_8(self, SBMLMetadata(''))
		self.r['reaction_9'] = reaction_9(self, SBMLMetadata(''))
		self.r['reaction_10'] = reaction_10(self, SBMLMetadata(''))
		self.r['reaction_11'] = reaction_11(self, SBMLMetadata(''))
		self.r['reaction_12'] = reaction_12(self, SBMLMetadata(''))
		self.r['reaction_13'] = reaction_13(self, SBMLMetadata(''))
		self.r['reaction_17'] = reaction_17(self, SBMLMetadata(''))
		self.r['reaction_18'] = reaction_18(self, SBMLMetadata(''))
		self.r['reaction_19'] = reaction_19(self, SBMLMetadata(''))
		self.r['reaction_20'] = reaction_20(self, SBMLMetadata(''))
		self.r['reaction_21'] = reaction_21(self, SBMLMetadata(''))
		self.r['reaction_22'] = reaction_22(self, SBMLMetadata(''))
		self.r['reaction_23'] = reaction_23(self, SBMLMetadata(''))
		self.r['reaction_24'] = reaction_24(self, SBMLMetadata(''))
		self.r['reaction_25'] = reaction_25(self, SBMLMetadata(''))
		self.r['reaction_26'] = reaction_26(self, SBMLMetadata(''))
		self.r['reaction_27'] = reaction_27(self, SBMLMetadata(''))
		self.r['reaction_28'] = reaction_28(self, SBMLMetadata(''))
		self.r['reaction_29'] = reaction_29(self, SBMLMetadata(''))
		self.r['reaction_30'] = reaction_30(self, SBMLMetadata(''))
		self.r['reaction_31'] = reaction_31(self, SBMLMetadata(''))
		self.r['reaction_34'] = reaction_34(self, SBMLMetadata(''))
		self.r['reaction_35'] = reaction_35(self, SBMLMetadata(''))
		self.r['reaction_36'] = reaction_36(self, SBMLMetadata(''))
		self.r['reaction_37'] = reaction_37(self, SBMLMetadata(''))
		self.r['reaction_38'] = reaction_38(self, SBMLMetadata(''))
		self.r['reaction_40'] = reaction_40(self, SBMLMetadata(''))
		self.r['reaction_41'] = reaction_41(self, SBMLMetadata(''))
		self.r['reaction_42'] = reaction_42(self, SBMLMetadata(''))
		self.r['reaction_43'] = reaction_43(self, SBMLMetadata(''))
		self.r['reaction_44'] = reaction_44(self, SBMLMetadata(''))
		self.r['reaction_45'] = reaction_45(self, SBMLMetadata(''))
		self.r['reaction_46'] = reaction_46(self, SBMLMetadata(''))
		self.r['reaction_47'] = reaction_47(self, SBMLMetadata(''))
		self.r['reaction_48'] = reaction_48(self, SBMLMetadata(''))
		self.r['reaction_49'] = reaction_49(self, SBMLMetadata(''))
		self.r['reaction_50'] = reaction_50(self, SBMLMetadata(''))
		self.r['reaction_51'] = reaction_51(self, SBMLMetadata(''))
		self.r['reaction_52'] = reaction_52(self, SBMLMetadata(''))
		self.r['reaction_53'] = reaction_53(self, SBMLMetadata(''))
		self.r['reaction_54'] = reaction_54(self, SBMLMetadata(''))
		self.r['reaction_55'] = reaction_55(self, SBMLMetadata(''))
		self.r['reaction_56'] = reaction_56(self, SBMLMetadata(''))
		self.r['reaction_57'] = reaction_57(self, SBMLMetadata(''))
		self.r['reaction_58'] = reaction_58(self, SBMLMetadata(''))
		self.r['reaction_59'] = reaction_59(self, SBMLMetadata(''))
		self.r['reaction_60'] = reaction_60(self, SBMLMetadata(''))
		self.r['reaction_61'] = reaction_61(self, SBMLMetadata(''))
		self.r['reaction_62'] = reaction_62(self, SBMLMetadata(''))
		self.r['reaction_63'] = reaction_63(self, SBMLMetadata(''))
		self.r['reaction_64'] = reaction_64(self, SBMLMetadata(''))
		self.r['reaction_65'] = reaction_65(self, SBMLMetadata(''))
		self.r['reaction_66'] = reaction_66(self, SBMLMetadata(''))
		self.r['reaction_67'] = reaction_67(self, SBMLMetadata(''))
		self.r['reaction_68'] = reaction_68(self, SBMLMetadata(''))
		self.r['reaction_69'] = reaction_69(self, SBMLMetadata(''))
		self.r['reaction_70'] = reaction_70(self, SBMLMetadata(''))
		self.r['reaction_71'] = reaction_71(self, SBMLMetadata(''))
		self.r['reaction_72'] = reaction_72(self, SBMLMetadata(''))
		self.r['reaction_73'] = reaction_73(self, SBMLMetadata(''))
		self.r['reaction_74'] = reaction_74(self, SBMLMetadata(''))
		self.r['reaction_75'] = reaction_75(self, SBMLMetadata(''))
		self.r['reaction_76'] = reaction_76(self, SBMLMetadata(''))
		self.r['reaction_77'] = reaction_77(self, SBMLMetadata(''))
		self.r['reaction_78'] = reaction_78(self, SBMLMetadata(''))
		self.r['reaction_79'] = reaction_79(self, SBMLMetadata(''))
		self.r['reaction_80'] = reaction_80(self, SBMLMetadata(''))
		self.r['reaction_81'] = reaction_81(self, SBMLMetadata(''))
		self.r['reaction_82'] = reaction_82(self, SBMLMetadata(''))
		self.r['reaction_83'] = reaction_83(self, SBMLMetadata(''))
		self.r['reaction_84'] = reaction_84(self, SBMLMetadata(''))
		self.r['reaction_85'] = reaction_85(self, SBMLMetadata(''))
		self.r['reaction_88'] = reaction_88(self, SBMLMetadata(''))
		self.r['reaction_89'] = reaction_89(self, SBMLMetadata(''))
		self.r['reaction_90'] = reaction_90(self, SBMLMetadata(''))
		self.r['reaction_91'] = reaction_91(self, SBMLMetadata(''))
		self.r['reaction_92'] = reaction_92(self, SBMLMetadata(''))
		self.r['reaction_93'] = reaction_93(self, SBMLMetadata(''))
		self.r['reaction_94'] = reaction_94(self, SBMLMetadata(''))
		self.r['reaction_95'] = reaction_95(self, SBMLMetadata(''))
		self.r['reaction_96'] = reaction_96(self, SBMLMetadata(''))
		self.r['reaction_97'] = reaction_97(self, SBMLMetadata(''))
		self.r['reaction_98'] = reaction_98(self, SBMLMetadata(''))
		self.r['reaction_99'] = reaction_99(self, SBMLMetadata(''))
		self.r['reaction_100'] = reaction_100(self, SBMLMetadata(''))
		self.r['reaction_101'] = reaction_101(self, SBMLMetadata(''))
		self.r['reaction_102'] = reaction_102(self, SBMLMetadata(''))
		self.r['reaction_103'] = reaction_103(self, SBMLMetadata(''))
		self.r['reaction_104'] = reaction_104(self, SBMLMetadata(''))
		self.r['reaction_105'] = reaction_105(self, SBMLMetadata(''))
		self.r['reaction_106'] = reaction_106(self, SBMLMetadata(''))
		self.r['reaction_107'] = reaction_107(self, SBMLMetadata(''))
		self.r['reaction_108'] = reaction_108(self, SBMLMetadata(''))
		self.r['reaction_109'] = reaction_109(self, SBMLMetadata(''))
		self.r['reaction_110'] = reaction_110(self, SBMLMetadata(''))
		self.r['reaction_111'] = reaction_111(self, SBMLMetadata(''))
		self.r['reaction_112'] = reaction_112(self, SBMLMetadata(''))
		self.r['reaction_113'] = reaction_113(self, SBMLMetadata(''))
		self.r['reaction_114'] = reaction_114(self, SBMLMetadata(''))
		self.r['reaction_115'] = reaction_115(self, SBMLMetadata(''))
		self.r['reaction_117'] = reaction_117(self, SBMLMetadata(''))
		self.r['reaction_118'] = reaction_118(self, SBMLMetadata(''))
		self.time = 0

		self.reactionMetadata = {
		self.Reactionreaction_1: SBMLMetadata(''),
		self.Reactionreaction_2: SBMLMetadata(''),
		self.Reactionreaction_3: SBMLMetadata(''),
		self.Reactionreaction_4: SBMLMetadata(''),
		self.Reactionreaction_5: SBMLMetadata(''),
		self.Reactionreaction_6: SBMLMetadata(''),
		self.Reactionreaction_7: SBMLMetadata(''),
		self.Reactionreaction_8: SBMLMetadata(''),
		self.Reactionreaction_9: SBMLMetadata(''),
		self.Reactionreaction_10: SBMLMetadata(''),
		self.Reactionreaction_11: SBMLMetadata(''),
		self.Reactionreaction_12: SBMLMetadata(''),
		self.Reactionreaction_13: SBMLMetadata(''),
		self.Reactionreaction_17: SBMLMetadata(''),
		self.Reactionreaction_18: SBMLMetadata(''),
		self.Reactionreaction_19: SBMLMetadata(''),
		self.Reactionreaction_20: SBMLMetadata(''),
		self.Reactionreaction_21: SBMLMetadata(''),
		self.Reactionreaction_22: SBMLMetadata(''),
		self.Reactionreaction_23: SBMLMetadata(''),
		self.Reactionreaction_24: SBMLMetadata(''),
		self.Reactionreaction_25: SBMLMetadata(''),
		self.Reactionreaction_26: SBMLMetadata(''),
		self.Reactionreaction_27: SBMLMetadata(''),
		self.Reactionreaction_28: SBMLMetadata(''),
		self.Reactionreaction_29: SBMLMetadata(''),
		self.Reactionreaction_30: SBMLMetadata(''),
		self.Reactionreaction_31: SBMLMetadata(''),
		self.Reactionreaction_34: SBMLMetadata(''),
		self.Reactionreaction_35: SBMLMetadata(''),
		self.Reactionreaction_36: SBMLMetadata(''),
		self.Reactionreaction_37: SBMLMetadata(''),
		self.Reactionreaction_38: SBMLMetadata(''),
		self.Reactionreaction_40: SBMLMetadata(''),
		self.Reactionreaction_41: SBMLMetadata(''),
		self.Reactionreaction_42: SBMLMetadata(''),
		self.Reactionreaction_43: SBMLMetadata(''),
		self.Reactionreaction_44: SBMLMetadata(''),
		self.Reactionreaction_45: SBMLMetadata(''),
		self.Reactionreaction_46: SBMLMetadata(''),
		self.Reactionreaction_47: SBMLMetadata(''),
		self.Reactionreaction_48: SBMLMetadata(''),
		self.Reactionreaction_49: SBMLMetadata(''),
		self.Reactionreaction_50: SBMLMetadata(''),
		self.Reactionreaction_51: SBMLMetadata(''),
		self.Reactionreaction_52: SBMLMetadata(''),
		self.Reactionreaction_53: SBMLMetadata(''),
		self.Reactionreaction_54: SBMLMetadata(''),
		self.Reactionreaction_55: SBMLMetadata(''),
		self.Reactionreaction_56: SBMLMetadata(''),
		self.Reactionreaction_57: SBMLMetadata(''),
		self.Reactionreaction_58: SBMLMetadata(''),
		self.Reactionreaction_59: SBMLMetadata(''),
		self.Reactionreaction_60: SBMLMetadata(''),
		self.Reactionreaction_61: SBMLMetadata(''),
		self.Reactionreaction_62: SBMLMetadata(''),
		self.Reactionreaction_63: SBMLMetadata(''),
		self.Reactionreaction_64: SBMLMetadata(''),
		self.Reactionreaction_65: SBMLMetadata(''),
		self.Reactionreaction_66: SBMLMetadata(''),
		self.Reactionreaction_67: SBMLMetadata(''),
		self.Reactionreaction_68: SBMLMetadata(''),
		self.Reactionreaction_69: SBMLMetadata(''),
		self.Reactionreaction_70: SBMLMetadata(''),
		self.Reactionreaction_71: SBMLMetadata(''),
		self.Reactionreaction_72: SBMLMetadata(''),
		self.Reactionreaction_73: SBMLMetadata(''),
		self.Reactionreaction_74: SBMLMetadata(''),
		self.Reactionreaction_75: SBMLMetadata(''),
		self.Reactionreaction_76: SBMLMetadata(''),
		self.Reactionreaction_77: SBMLMetadata(''),
		self.Reactionreaction_78: SBMLMetadata(''),
		self.Reactionreaction_79: SBMLMetadata(''),
		self.Reactionreaction_80: SBMLMetadata(''),
		self.Reactionreaction_81: SBMLMetadata(''),
		self.Reactionreaction_82: SBMLMetadata(''),
		self.Reactionreaction_83: SBMLMetadata(''),
		self.Reactionreaction_84: SBMLMetadata(''),
		self.Reactionreaction_85: SBMLMetadata(''),
		self.Reactionreaction_88: SBMLMetadata(''),
		self.Reactionreaction_89: SBMLMetadata(''),
		self.Reactionreaction_90: SBMLMetadata(''),
		self.Reactionreaction_91: SBMLMetadata(''),
		self.Reactionreaction_92: SBMLMetadata(''),
		self.Reactionreaction_93: SBMLMetadata(''),
		self.Reactionreaction_94: SBMLMetadata(''),
		self.Reactionreaction_95: SBMLMetadata(''),
		self.Reactionreaction_96: SBMLMetadata(''),
		self.Reactionreaction_97: SBMLMetadata(''),
		self.Reactionreaction_98: SBMLMetadata(''),
		self.Reactionreaction_99: SBMLMetadata(''),
		self.Reactionreaction_100: SBMLMetadata(''),
		self.Reactionreaction_101: SBMLMetadata(''),
		self.Reactionreaction_102: SBMLMetadata(''),
		self.Reactionreaction_103: SBMLMetadata(''),
		self.Reactionreaction_104: SBMLMetadata(''),
		self.Reactionreaction_105: SBMLMetadata(''),
		self.Reactionreaction_106: SBMLMetadata(''),
		self.Reactionreaction_107: SBMLMetadata(''),
		self.Reactionreaction_108: SBMLMetadata(''),
		self.Reactionreaction_109: SBMLMetadata(''),
		self.Reactionreaction_110: SBMLMetadata(''),
		self.Reactionreaction_111: SBMLMetadata(''),
		self.Reactionreaction_112: SBMLMetadata(''),
		self.Reactionreaction_113: SBMLMetadata(''),
		self.Reactionreaction_114: SBMLMetadata(''),
		self.Reactionreaction_115: SBMLMetadata(''),
		self.Reactionreaction_117: SBMLMetadata(''),
		self.Reactionreaction_118: SBMLMetadata('')
		}
		self.AssignmentRules()



	def AssignmentRules(self):

		self.s['phosphorylated_Akt'].concentration = self.s['pAkt'].concentration + self.s['ppAkt'].concentration

		self.p['EGF_tot'].value = self.s['EGF'].concentration + (self.s['RE'].concentration + 2 * (self.s['Rd'].concentration + self.s['Rp'].concentration + self.s['Ri'].concentration + self.s['Rp_GS'].concentration + self.s['Rp_Shc'].concentration + self.s['Rp_pShc'].concentration + self.s['Rp_pShc_GS'].concentration + self.s['Rp_PI3K'].concentration + self.s['Rp_RasGAP'].concentration)) * (self.c['cell'].size / self.c['extra'].size)

		self.p['k_1'].value = self.p['Kd1'].value * self.p['k1'].value

		self.p['k_2'].value = self.p['Kd2'].value * self.p['k2'].value

		self.p['k_4'].value = self.p['Kd4'].value * self.p['k4'].value

		self.p['k_5'].value = self.p['Kd5'].value * self.p['k5'].value

		self.p['k_7'].value = self.p['Kd7'].value * self.p['k7'].value

		self.p['k_9'].value = self.p['Kd9'].value * self.p['k9'].value

		self.p['k_10'].value = self.p['Kd10'].value * self.p['k10'].value

		self.p['k_11'].value = self.p['k_9'].value

		self.p['k11'].value = self.p['k9'].value

		self.p['k_12'].value = self.p['Kd12'].value * self.p['k12'].value

		self.p['k_13'].value = self.p['Kd13'].value * self.p['k13'].value

		self.p['k_24'].value = self.p['Kd24'].value * self.p['k24'].value

		self.p['k_26'].value = self.p['Kd26'].value * self.p['k26'].value

		self.p['k_27'].value = self.p['Kd27'].value * self.p['k27'].value

		self.p['k_28'].value = self.p['Kd28'].value * self.p['k28'].value

		self.p['k_30'].value = self.p['Kd30'].value * self.p['k30'].value

		self.p['k_42'].value = self.p['Kd42'].value * self.p['k42'].value

		self.p['k_45'].value = self.p['k45'].value * self.p['Kd45'].value

		self.p['k_46'].value = self.p['Kd46'].value * self.p['k46'].value

		self.p['k_47'].value = self.p['Kd47'].value * self.p['k47'].value

		self.p['k_49'].value = self.p['Kd49'].value * self.p['k49'].value

		self.p['k_52'].value = self.p['k52'].value * self.p['Kd52'].value

		self.p['k_53'].value = self.p['Kd53'].value * self.p['k53'].value

		self.p['k_54'].value = self.p['Kd54'].value * self.p['k54'].value

		self.p['k_55'].value = self.p['Kd55'].value * self.p['k55'].value

		self.p['k_59'].value = self.p['Kd59'].value * self.p['k59'].value

		self.p['k_74'].value = self.p['k74'].value * self.p['Kd74'].value

		return

	def Reactionreaction_1(self):

		return (self.p['k1'].value * self.s['R'].concentration * self.s['EGF'].concentration - self.p['k_1'].value * self.s['RE'].concentration) * self.c['cell'].size

	def Reactionreaction_2(self):

		return (self.p['k2'].value * self.s['RE'].concentration * self.s['RE'].concentration - self.p['k_2'].value * self.s['Rd'].concentration) * self.c['cell'].size

	def Reactionreaction_3(self):

		return self.p['k3'].value * self.s['Rd'].concentration * self.c['cell'].size

	def Reactionreaction_4(self):

		return (self.p['k4'].value * self.s['Rp'].concentration * self.s['GS'].concentration - self.p['k_4'].value * self.s['Rp_GS'].concentration) * self.c['cell'].size

	def Reactionreaction_5(self):

		return (self.p['k5'].value * self.s['Rp'].concentration * self.s['Shc'].concentration - self.p['k_5'].value * self.s['Rp_Shc'].concentration) * self.c['cell'].size

	def Reactionreaction_6(self):

		return self.p['k6'].value * self.s['Rp_Shc'].concentration * self.c['cell'].size

	def Reactionreaction_7(self):

		return (self.p['k_7'].value * self.s['Rp_pShc'].concentration - self.p['k7'].value * self.s['Rp'].concentration * self.s['pShc'].concentration) * self.c['cell'].size

	def Reactionreaction_8(self):

		return (self.p['V8'].value * self.s['pShc'].concentration / (self.p['Km8'].value + self.s['pShc'].concentration)) * self.c['cell'].size

	def Reactionreaction_9(self):

		return (self.p['k9'].value * self.s['Rp_pShc'].concentration * self.s['GS'].concentration - self.p['k_9'].value * self.s['Rp_pShc_GS'].concentration) * self.c['cell'].size

	def Reactionreaction_10(self):

		return (self.p['k_10'].value * self.s['Rp_pShc_GS'].concentration - self.p['k10'].value * self.s['Rp'].concentration * self.s['pShc_GS'].concentration) * self.c['cell'].size

	def Reactionreaction_11(self):

		return (self.p['k_11'].value * self.s['pShc_GS'].concentration - self.p['k11'].value * self.s['pShc'].concentration * self.s['GS'].concentration) * self.c['cell'].size

	def Reactionreaction_12(self):

		return (self.p['k12'].value * self.s['Rp'].concentration * self.s['PI3K'].concentration - self.p['k_12'].value * self.s['Rp_PI3K'].concentration) * self.c['cell'].size

	def Reactionreaction_13(self):

		return (self.p['k13'].value * self.s['Rp'].concentration * self.s['RasGAP'].concentration - self.p['k_13'].value * self.s['Rp_RasGAP'].concentration) * self.c['cell'].size

	def Reactionreaction_17(self):

		return self.p['k17'].value * self.s['Rp'].concentration * self.c['cell'].size

	def Reactionreaction_18(self):

		return self.p['k17'].value * self.s['Rp_GS'].concentration * self.c['cell'].size

	def Reactionreaction_19(self):

		return self.p['k17'].value * self.s['Rp_Shc'].concentration * self.c['cell'].size

	def Reactionreaction_20(self):

		return self.p['k17'].value * self.s['Rp_pShc'].concentration * self.c['cell'].size

	def Reactionreaction_21(self):

		return self.p['k17'].value * self.s['Rp_pShc_GS'].concentration * self.c['cell'].size

	def Reactionreaction_22(self):

		return self.p['k17'].value * self.s['Rp_PI3K'].concentration * self.c['cell'].size

	def Reactionreaction_23(self):

		return self.p['k17'].value * self.s['Rp_RasGAP'].concentration * self.c['cell'].size

	def Reactionreaction_24(self):

		return (self.p['k24'].value * self.s['IR'].concentration * self.s['I'].concentration - self.p['k_24'].value * self.s['IRL'].concentration) * self.c['cell'].size

	def Reactionreaction_25(self):

		return self.p['k25'].value * self.s['IRL'].concentration * self.c['cell'].size

	def Reactionreaction_26(self):

		return (self.p['k26'].value * self.s['IRp'].concentration * self.s['PI3K'].concentration - self.p['k_26'].value * self.s['IRp_PI3K'].concentration) * self.c['cell'].size

	def Reactionreaction_27(self):

		return (self.p['k27'].value * self.s['IRp'].concentration * self.s['RasGAP'].concentration - self.p['k_27'].value * self.s['IRp_RasGAP'].concentration) * self.c['cell'].size

	def Reactionreaction_28(self):

		return (self.p['k28'].value * self.s['IRp'].concentration * self.s['IRS'].concentration - self.p['k_28'].value * self.s['IRp_IRS'].concentration) * self.c['cell'].size

	def Reactionreaction_29(self):

		return self.p['k29'].value * self.s['IRp_IRS'].concentration * self.c['cell'].size

	def Reactionreaction_30(self):

		return (self.p['k_30'].value * self.s['IRp_IRSp'].concentration - self.p['k30'].value * self.s['IRp'].concentration * self.s['IRSp'].concentration) * self.c['cell'].size

	def Reactionreaction_31(self):

		return (self.p['V31'].value * self.s['IRSp'].concentration / (self.p['Km31'].value + self.s['IRSp'].concentration)) * self.c['cell'].size

	def Reactionreaction_34(self):

		return self.p['k17'].value * self.s['IRp'].concentration * self.c['cell'].size

	def Reactionreaction_35(self):

		return self.p['k17'].value * self.s['IRp_PI3K'].concentration * self.c['cell'].size

	def Reactionreaction_36(self):

		return self.p['k17'].value * self.s['IRp_RasGAP'].concentration * self.c['cell'].size

	def Reactionreaction_37(self):

		return self.p['k17'].value * self.s['IRp_IRS'].concentration * self.c['cell'].size

	def Reactionreaction_38(self):

		return self.p['k17'].value * self.s['IRp_IRSp'].concentration * self.c['cell'].size

	def Reactionreaction_40(self):

		return (self.p['kcat40'].value * self.s['iSrc'].concentration * (self.s['Rp'].concentration + self.p['alpha40'].value * self.s['IRp'].concentration) / (self.p['Km40'].value + self.s['iSrc'].concentration)) * self.c['cell'].size

	def Reactionreaction_41(self):

		return (self.p['V41'].value * self.s['aSrc'].concentration / (self.p['Km41'].value + self.s['aSrc'].concentration)) * self.c['cell'].size

	def Reactionreaction_42(self):

		return (self.p['k42'].value * self.s['IRS'].concentration * self.s['PIP3'].concentration - self.p['k_42'].value * self.s['mIRS'].concentration) * self.c['cell'].size

	def Reactionreaction_43(self):

		return (self.p['kcat43'].value * self.s['mIRS'].concentration * (self.s['IRp'].concentration + self.p['alpha43'].value * self.s['Rp'].concentration) / (self.p['Km43'].value + self.s['mIRS'].concentration)) * self.c['cell'].size

	def Reactionreaction_44(self):

		return (self.p['V31'].value * self.s['mIRSp'].concentration / (self.p['Km31'].value + self.s['mIRSp'].concentration)) * self.c['cell'].size

	def Reactionreaction_45(self):

		return (self.p['k45'].value * self.s['mIRSp'].concentration * self.s['GS'].concentration - self.p['k_45'].value * self.s['mIRSp_GS'].concentration) * self.c['cell'].size

	def Reactionreaction_46(self):

		return (self.p['k46'].value * self.s['mIRSp'].concentration * self.s['PI3K'].concentration - self.p['k_46'].value * self.s['mIRSp_PI3K'].concentration) * self.c['cell'].size

	def Reactionreaction_47(self):

		return (self.p['k47'].value * self.s['mIRSp'].concentration * self.s['SHP2'].concentration - self.p['k_47'].value * self.s['mIRSp_SHP2'].concentration) * self.c['cell'].size

	def Reactionreaction_48(self):

		return self.p['k48'].value * self.s['mIRSp_SHP2'].concentration * self.c['cell'].size

	def Reactionreaction_49(self):

		return (self.p['k49'].value * self.s['GAB'].concentration * self.s['PIP3'].concentration - self.p['k_49'].value * self.s['mGAB'].concentration) * self.c['cell'].size

	def Reactionreaction_50(self):

		return (self.p['kcat50'].value * self.s['mGAB'].concentration * (self.s['Rp'].concentration + self.p['alpha50'].value * self.s['aSrc'].concentration) / (self.p['Km50'].value + self.s['mGAB'].concentration)) * self.c['cell'].size

	def Reactionreaction_51(self):

		return (self.p['V51'].value * self.s['mGABp'].concentration / (self.p['Km51'].value + self.s['mGABp'].concentration)) * self.c['cell'].size

	def Reactionreaction_52(self):

		return (self.p['k52'].value * self.s['mGABp'].concentration * self.s['GS'].concentration - self.p['k_52'].value * self.s['mGABp_GS'].concentration) * self.c['cell'].size

	def Reactionreaction_53(self):

		return (self.p['k53'].value * self.s['mGABp'].concentration * self.s['PI3K'].concentration - self.p['k_53'].value * self.s['mGABp_PI3K'].concentration) * self.c['cell'].size

	def Reactionreaction_54(self):

		return (self.p['k54'].value * self.s['mGABp'].concentration * self.s['RasGAP'].concentration - self.p['k_54'].value * self.s['mGABp_RasGAP'].concentration) * self.c['cell'].size

	def Reactionreaction_55(self):

		return (self.p['k55'].value * self.s['mGABp'].concentration * self.s['SHP2'].concentration - self.p['k_55'].value * self.s['mGABp_SHP2'].concentration) * self.c['cell'].size

	def Reactionreaction_56(self):

		return self.p['k56'].value * self.s['mGABp_SHP2'].concentration * self.c['cell'].size

	def Reactionreaction_57(self):

		return (self.p['kcat57'].value * self.s['mGABp_SHP2'].concentration * (self.s['Rp'].concentration + self.s['aSrc'].concentration) / (self.p['Km57'].value + self.s['mGABp_SHP2'].concentration)) * self.c['cell'].size

	def Reactionreaction_58(self):

		return (self.p['V58'].value * self.s['mGABp_pSHP2'].concentration / (self.p['Km58'].value + self.s['mGABp_pSHP2'].concentration)) * self.c['cell'].size

	def Reactionreaction_59(self):

		return (self.p['k59'].value * self.s['mGABp_pSHP2'].concentration * self.s['GS'].concentration - self.p['k_59'].value * self.s['mGABp_pSHP2_GS'].concentration) * self.c['cell'].size

	def Reactionreaction_60(self):

		return self.p['k60'].value * (self.s['Rp_PI3K'].concentration + self.s['IRp_PI3K'].concentration + self.s['mIRSp_PI3K'].concentration + self.s['mGABp_PI3K'].concentration + self.s['tRas_PI3K'].concentration) * self.c['cell'].size

	def Reactionreaction_61(self):

		return self.p['k61'].value * self.s['PIP3'].concentration * self.c['cell'].size

	def Reactionreaction_62(self):

		return (self.p['kcat62'].value * self.s['dRas'].concentration * (self.s['Rp_GS'].concentration + self.s['Rp_pShc_GS'].concentration + self.s['mIRSp_GS'].concentration + self.s['mGABp_GS'].concentration + self.s['mGABp_pSHP2_GS'].concentration) / (self.p['Km62'].value + self.s['dRas'].concentration)) * self.c['cell'].size

	def Reactionreaction_63(self):

		return (self.p['kcat63'].value * self.s['tRas'].concentration * (self.s['bRasGAP'].concentration + self.s['mGABp_RasGAP'].concentration + self.s['Rp_RasGAP'].concentration + self.s['IRp_RasGAP'].concentration) / (self.p['Km63'].value + self.s['tRas'].concentration)) * self.c['cell'].size

	def Reactionreaction_64(self):

		return (self.p['k64'].value * self.s['tRas'].concentration * self.s['PI3K'].concentration - self.p['k_64'].value * self.s['tRas_PI3K'].concentration) * self.c['cell'].size

	def Reactionreaction_65(self):

		return (self.p['kcat65'].value * self.s['tRas'].concentration * self.s['Raf'].concentration / (self.p['Km65'].value + self.s['Raf'].concentration)) * self.c['cell'].size

	def Reactionreaction_66(self):

		return (self.p['kcat66'].value * self.s['aSrc'].concentration * self.s['aRaf'].concentration / (self.p['Km66'].value + self.s['aRaf'].concentration)) * self.c['cell'].size

	def Reactionreaction_67(self):

		return (self.p['kcat67'].value * self.s['aaRaf'].concentration * (self.s['PKA'].concentration / (self.p['Km67'].value + self.s['aaRaf'].concentration)) + self.p['alpha67'].value * self.s['aaRaf'].concentration * (self.s['pAkt'].concentration + self.p['beta67'].value * self.s['ppAkt'].concentration)) * self.c['cell'].size

	def Reactionreaction_68(self):

		return (self.p['kcat68'].value * self.s['aaRaf'].concentration * self.s['Mek'].concentration / (self.p['Km68'].value + self.s['Mek'].concentration)) * self.c['cell'].size

	def Reactionreaction_69(self):

		return (self.p['V69'].value * self.s['ppMek'].concentration / (self.p['Km69'].value + self.s['ppMek'].concentration)) * self.c['cell'].size

	def Reactionreaction_70(self):

		return (self.p['kcat70'].value * self.s['Erk'].concentration * self.s['ppMek'].concentration / (self.p['Km70'].value + self.s['Erk'].concentration + self.s['pErk'].concentration * (self.p['Km70'].value / self.p['Km71'].value))) * self.c['cell'].size

	def Reactionreaction_71(self):

		return (self.p['kcat71'].value * self.s['pErk'].concentration * self.s['ppMek'].concentration / (self.p['Km71'].value + self.s['pErk'].concentration + self.s['Erk'].concentration * (self.p['Km71'].value / self.p['Km70'].value))) * self.c['cell'].size

	def Reactionreaction_72(self):

		return (self.p['V72'].value * self.s['ppErk'].concentration / (self.p['Km72'].value + self.s['ppErk'].concentration + self.s['pErk'].concentration * (self.p['Km72'].value / self.p['Km73'].value))) * self.c['cell'].size

	def Reactionreaction_73(self):

		return (self.p['V73'].value * self.s['pErk'].concentration / (self.p['Km73'].value + self.s['pErk'].concentration + self.s['ppErk'].concentration * (self.p['Km73'].value / self.p['Km72'].value))) * self.c['cell'].size

	def Reactionreaction_74(self):

		return (self.p['k74'].value * self.s['PDK1'].concentration * self.s['PIP3'].concentration - self.p['k_74'].value * self.s['mPDK1'].concentration) * self.c['cell'].size

	def Reactionreaction_75(self):

		return (self.p['kcat75'].value * self.s['mPDK1'].concentration * self.s['Akt'].concentration / (self.p['Km75'].value + self.s['Akt'].concentration)) * self.c['cell'].size

	def Reactionreaction_76(self):

		return (self.p['V76'].value * self.s['pAkt'].concentration / (self.p['Km76'].value + self.s['pAkt'].concentration)) * self.c['cell'].size

	def Reactionreaction_77(self):

		return (self.p['kcat77'].value * self.s['mTOR'].concentration * ((self.p['alpha77'].value * self.s['pAkt'].concentration + self.s['ppAkt'].concentration) / (self.p['Km77'].value + self.s['mTOR'].concentration)) - self.p['k_77'].value * self.s['amTOR'].concentration) * self.c['cell'].size

	def Reactionreaction_78(self):

		return (self.p['kcat78'].value * self.s['amTOR'].concentration * (self.s['pAkt'].concentration / (self.p['Km78'].value + self.s['pAkt'].concentration)) - self.p['k_78'].value * self.s['ppAkt'].concentration) * self.c['cell'].size

	def Reactionreaction_79(self):

		return (self.p['kcat79'].value * self.s['ppErk'].concentration * (self.s['GS'].concentration / (self.p['Km79'].value + self.s['GS'].concentration)) - self.p['k_79'].value * self.s['iGS'].concentration) * self.c['cell'].size

	def Reactionreaction_80(self):

		return (self.p['kcat80'].value * self.s['mGAB'].concentration * (self.s['ppErk'].concentration / (self.p['Km80'].value + self.s['mGAB'].concentration)) - self.p['k_80'].value * self.s['imGAB'].concentration) * self.c['cell'].size

	def Reactionreaction_81(self):

		return (self.p['kcat81'].value * self.s['mIRS'].concentration * (self.s['amTOR'].concentration / (self.p['Km81'].value + self.s['mIRS'].concentration)) - self.p['k_81'].value * self.s['imIRS'].concentration) * self.c['cell'].size

	def Reactionreaction_82(self):

		return (self.p['V_82'].value * self.s['Rp'].concentration / (self.p['Km82'].value + self.s['Rp'].concentration)) * self.c['cell'].size

	def Reactionreaction_83(self):

		return self.p['k83'].value * self.s['Ri'].concentration * self.c['cell'].size

	def Reactionreaction_84(self):

		return (self.p['V_84'].value * self.s['IRp'].concentration / (self.p['Km84'].value + self.s['IRp'].concentration)) * self.c['cell'].size

	def Reactionreaction_85(self):

		return self.p['k85'].value * self.s['IRi'].concentration * self.c['cell'].size

	def Reactionreaction_88(self):

		return (self.p['k_42'].value * self.s['mIRSp'].concentration - self.p['k42'].value * self.s['PIP3'].concentration * self.s['IRSp'].concentration) * self.c['cell'].size

	def Reactionreaction_89(self):

		return (self.p['k_42'].value * self.s['mIRSp_PI3K'].concentration - self.p['k42'].value * self.s['PIP3'].concentration * self.s['IRSp_PI3K'].concentration) * self.c['cell'].size

	def Reactionreaction_90(self):

		return (self.p['k_42'].value * self.s['mIRSp_GS'].concentration - self.p['k42'].value * self.s['PIP3'].concentration * self.s['IRSp_GS'].concentration) * self.c['cell'].size

	def Reactionreaction_91(self):

		return (self.p['k_42'].value * self.s['mIRSp_SHP2'].concentration - self.p['k42'].value * self.s['PIP3'].concentration * self.s['IRSp_SHP2'].concentration) * self.c['cell'].size

	def Reactionreaction_92(self):

		return (self.p['k_42'].value * self.s['mGABp'].concentration - self.p['k42'].value * self.s['PIP3'].concentration * self.s['GABp'].concentration) * self.c['cell'].size

	def Reactionreaction_93(self):

		return (self.p['k_42'].value * self.s['mGABp_PI3K'].concentration - self.p['k42'].value * self.s['PIP3'].concentration * self.s['GABp_PI3K'].concentration) * self.c['cell'].size

	def Reactionreaction_94(self):

		return (self.p['k_42'].value * self.s['mGABp_GS'].concentration - self.p['k42'].value * self.s['PIP3'].concentration * self.s['GABp_GS'].concentration) * self.c['cell'].size

	def Reactionreaction_95(self):

		return (self.p['k_42'].value * self.s['mGABp_RasGAP'].concentration - self.p['k49'].value * self.s['PIP3'].concentration * self.s['GABp_RasGAP'].concentration) * self.c['cell'].size

	def Reactionreaction_96(self):

		return (self.p['k_42'].value * self.s['mGABp_SHP2'].concentration - self.p['k42'].value * self.s['PIP3'].concentration * self.s['GABp_SHP2'].concentration) * self.c['cell'].size

	def Reactionreaction_97(self):

		return (self.p['k_42'].value * self.s['mGABp_pSHP2'].concentration - self.p['k42'].value * self.s['PIP3'].concentration * self.s['GABp_pSHP2'].concentration) * self.c['cell'].size

	def Reactionreaction_98(self):

		return (self.p['k_42'].value * self.s['mGABp_pSHP2_GS'].concentration - self.p['k42'].value * self.s['PIP3'].concentration * self.s['GABp_pSHP2_GS'].concentration) * self.c['cell'].size

	def Reactionreaction_99(self):

		return (self.p['V31'].value * self.s['IRSp_PI3K'].concentration / (self.p['Km31'].value + self.s['IRSp_PI3K'].concentration)) * self.c['cell'].size

	def Reactionreaction_100(self):

		return (self.p['V31'].value * self.s['IRSp_GS'].concentration / (self.p['Km31'].value + self.s['IRSp_GS'].concentration)) * self.c['cell'].size

	def Reactionreaction_101(self):

		return self.p['k48'].value * self.s['IRSp_SHP2'].concentration * self.c['cell'].size

	def Reactionreaction_102(self):

		return self.p['k56'].value * self.s['mGABp_pSHP2'].concentration * self.c['cell'].size

	def Reactionreaction_103(self):

		return self.p['k56'].value * self.s['mGABp_pSHP2_GS'].concentration * self.c['cell'].size

	def Reactionreaction_104(self):

		return (self.p['V51'].value * self.s['GABp'].concentration / (self.p['Km51'].value + self.s['GABp'].concentration)) * self.c['cell'].size

	def Reactionreaction_105(self):

		return (self.p['V51'].value * self.s['GABp_PI3K'].concentration / (self.p['Km51'].value + self.s['GABp_PI3K'].concentration)) * self.c['cell'].size

	def Reactionreaction_106(self):

		return (self.p['V51'].value * self.s['GABp_GS'].concentration / (self.p['Km51'].value + self.s['GABp_GS'].concentration)) * self.c['cell'].size

	def Reactionreaction_107(self):

		return (self.p['V51'].value * self.s['GABp_RasGAP'].concentration / (self.p['Km51'].value + self.s['GABp_RasGAP'].concentration)) * self.c['cell'].size

	def Reactionreaction_108(self):

		return self.p['k56'].value * self.s['GABp_SHP2'].concentration * self.c['cell'].size

	def Reactionreaction_109(self):

		return self.p['k56'].value * self.s['GABp_pSHP2'].concentration * self.c['cell'].size

	def Reactionreaction_110(self):

		return self.p['k56'].value * self.s['GABp_pSHP2_GS'].concentration * self.c['cell'].size

	def Reactionreaction_111(self):

		return self.p['k111'].value * (self.s['mGABp_SHP2'].concentration + self.s['mGABp_pSHP2'].concentration + self.s['mGABp_pSHP2_GS'].concentration + self.s['mIRSp_SHP2'].concentration) * self.s['mGABp_RasGAP'].concentration * self.c['cell'].size

	def Reactionreaction_112(self):

		return self.p['k111'].value * (self.s['mGABp_SHP2'].concentration + self.s['mGABp_pSHP2'].concentration + self.s['mGABp_pSHP2_GS'].concentration) * self.s['Rp_RasGAP'].concentration * self.c['cell'].size

	def Reactionreaction_113(self):

		return self.p['k111'].value * (self.s['mGABp_SHP2'].concentration + self.s['mGABp_pSHP2'].concentration + self.s['mGABp_pSHP2_GS'].concentration) * self.s['IRp_RasGAP'].concentration * self.c['cell'].size

	def Reactionreaction_114(self):

		return self.p['k111'].value * self.s['mIRSp_SHP2'].concentration * self.s['Rp_RasGAP'].concentration * self.c['cell'].size

	def Reactionreaction_115(self):

		return self.p['k111'].value * self.s['mIRSp_SHP2'].concentration * self.s['IRp_RasGAP'].concentration * self.c['cell'].size

	def Reactionreaction_117(self):

		return (2 * self.p['kcat80'].value * self.s['mGABp'].concentration * self.s['ppErk'].concentration / (self.p['Km80'].value + self.s['mGABp'].concentration) - self.p['k_80'].value * self.s['imGABp'].concentration) * self.c['cell'].size

	def Reactionreaction_118(self):

		return self.p['k118'].value * self.s['imGABp'].concentration * self.c['cell'].size

	def _SolveReactions(self, y, t):

		self.time = t
		self.s['EGF'].amount, self.s['I'].amount, self.s['RE'].amount, self.s['Rd'].amount, self.s['Rp'].amount, self.s['GS'].amount, self.s['Rp_GS'].amount, self.s['Shc'].amount, self.s['Rp_Shc'].amount, self.s['Rp_pShc'].amount, self.s['pShc'].amount, self.s['Rp_pShc_GS'].amount, self.s['PI3K'].amount, self.s['Rp_PI3K'].amount, self.s['RasGAP'].amount, self.s['Rp_RasGAP'].amount, self.s['IRL'].amount, self.s['IRp'].amount, self.s['IRp_PI3K'].amount, self.s['IRp_RasGAP'].amount, self.s['IRS'].amount, self.s['IRp_IRS'].amount, self.s['IRp_IRSp'].amount, self.s['IRSp'].amount, self.s['iSrc'].amount, self.s['mIRS'].amount, self.s['mIRSp'].amount, self.s['mIRSp_GS'].amount, self.s['mIRSp_PI3K'].amount, self.s['SHP2'].amount, self.s['mIRSp_SHP2'].amount, self.s['GAB'].amount, self.s['mGAB'].amount, self.s['mGABp'].amount, self.s['mGABp_GS'].amount, self.s['mGABp_PI3K'].amount, self.s['mGABp_SHP2'].amount, self.s['mGABp_pSHP2'].amount, self.s['PIP3'].amount, self.s['dRas'].amount, self.s['Raf'].amount, self.s['aRaf'].amount, self.s['Mek'].amount, self.s['Erk'].amount, self.s['pErk'].amount, self.s['PDK1'].amount, self.s['Akt'].amount, self.s['pAkt'].amount, self.s['mTOR'].amount, self.s['Null'].amount, self.s['aaRaf'].amount, self.s['PKA'].amount, self.s['pShc_GS'].amount, self.s['ppMek'].amount, self.s['mGABp_pSHP2_GS'].amount, self.s['R'].amount, self.s['ppErk'].amount, self.s['IR'].amount, self.s['mPDK1'].amount, self.s['tRas'].amount, self.s['tRas_PI3K'].amount, self.s['ppAkt'].amount, self.s['mGABp_RasGAP'].amount, self.s['amTOR'].amount, self.s['iGS'].amount, self.s['imGAB'].amount, self.s['imIRS'].amount, self.s['aSrc'].amount, self.s['Ri'].amount, self.s['IRi'].amount, self.s['iPX'].amount, self.s['aPX'].amount, self.s['aPX_GS'].amount, self.s['IRSp_PI3K'].amount, self.s['IRSp_GS'].amount, self.s['IRSp_SHP2'].amount, self.s['GABp'].amount, self.s['GABp_PI3K'].amount, self.s['GABp_GS'].amount, self.s['GABp_RasGAP'].amount, self.s['GABp_SHP2'].amount, self.s['GABp_pSHP2'].amount, self.s['GABp_pSHP2_GS'].amount, self.s['imGABp'].amount, self.s['bRasGAP'].amount, self.s['phosphorylated_Akt'].amount = y
		self.AssignmentRules()

		rateRuleVector = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype = np.float64)

		stoichiometricMatrix = np.array([[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 1,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,1,-1,-1,0,1,0,0,1,0,-1,-1,-1,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0.,0.],[ 0,0,0,-1,0,0,0,0,-1,0,1,0,0,0,1,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0.,0.],[ 0,0,0,1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0,1,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,-1,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,1,-1,0,-1,0,0,0,0,0,0,0,-1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,1,-1,0,0,1,0,0,0,0,0,1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,-1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0.,1,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,0,0,0.,0,0,0,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0.,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0.,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0.,0,1,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,1,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0.,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,-1,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,1,-1,-1,-1,0,1,0,-1,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,-1,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,-1,0,0,1,0,0,0,1,0,0,0,-1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,1,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,-1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,-1,1,0,0,0,0,0,0,-1,1,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,1,-1,1,0,0,0,0,1,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,1,-1,-1,-1,-1,-1,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0.,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,1.,-1,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0.,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,1,-1,1,-1,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,-1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,1.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0.,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0.,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0.,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0.,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0.,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0.,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0.,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0.,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0.,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.,-1.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.],[ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.,0.]], dtype = np.float64)

		reactionVelocities = np.array([self.r['reaction_1'](), self.r['reaction_2'](), self.r['reaction_3'](), self.r['reaction_4'](), self.r['reaction_5'](), self.r['reaction_6'](), self.r['reaction_7'](), self.r['reaction_8'](), self.r['reaction_9'](), self.r['reaction_10'](), self.r['reaction_11'](), self.r['reaction_12'](), self.r['reaction_13'](), self.r['reaction_17'](), self.r['reaction_18'](), self.r['reaction_19'](), self.r['reaction_20'](), self.r['reaction_21'](), self.r['reaction_22'](), self.r['reaction_23'](), self.r['reaction_24'](), self.r['reaction_25'](), self.r['reaction_26'](), self.r['reaction_27'](), self.r['reaction_28'](), self.r['reaction_29'](), self.r['reaction_30'](), self.r['reaction_31'](), self.r['reaction_34'](), self.r['reaction_35'](), self.r['reaction_36'](), self.r['reaction_37'](), self.r['reaction_38'](), self.r['reaction_40'](), self.r['reaction_41'](), self.r['reaction_42'](), self.r['reaction_43'](), self.r['reaction_44'](), self.r['reaction_45'](), self.r['reaction_46'](), self.r['reaction_47'](), self.r['reaction_48'](), self.r['reaction_49'](), self.r['reaction_50'](), self.r['reaction_51'](), self.r['reaction_52'](), self.r['reaction_53'](), self.r['reaction_54'](), self.r['reaction_55'](), self.r['reaction_56'](), self.r['reaction_57'](), self.r['reaction_58'](), self.r['reaction_59'](), self.r['reaction_60'](), self.r['reaction_61'](), self.r['reaction_62'](), self.r['reaction_63'](), self.r['reaction_64'](), self.r['reaction_65'](), self.r['reaction_66'](), self.r['reaction_67'](), self.r['reaction_68'](), self.r['reaction_69'](), self.r['reaction_70'](), self.r['reaction_71'](), self.r['reaction_72'](), self.r['reaction_73'](), self.r['reaction_74'](), self.r['reaction_75'](), self.r['reaction_76'](), self.r['reaction_77'](), self.r['reaction_78'](), self.r['reaction_79'](), self.r['reaction_80'](), self.r['reaction_81'](), self.r['reaction_82'](), self.r['reaction_83'](), self.r['reaction_84'](), self.r['reaction_85'](), self.r['reaction_88'](), self.r['reaction_89'](), self.r['reaction_90'](), self.r['reaction_91'](), self.r['reaction_92'](), self.r['reaction_93'](), self.r['reaction_94'](), self.r['reaction_95'](), self.r['reaction_96'](), self.r['reaction_97'](), self.r['reaction_98'](), self.r['reaction_99'](), self.r['reaction_100'](), self.r['reaction_101'](), self.r['reaction_102'](), self.r['reaction_103'](), self.r['reaction_104'](), self.r['reaction_105'](), self.r['reaction_106'](), self.r['reaction_107'](), self.r['reaction_108'](), self.r['reaction_109'](), self.r['reaction_110'](), self.r['reaction_111'](), self.r['reaction_112'](), self.r['reaction_113'](), self.r['reaction_114'](), self.r['reaction_115'](), self.r['reaction_117'](), self.r['reaction_118']()], dtype = np.float64)

		rateOfSpeciesChange = stoichiometricMatrix @ reactionVelocities + rateRuleVector

		return rateOfSpeciesChange

	def RunSimulation(self, deltaT, absoluteTolerance = 1e-12, relativeTolerance = 1e-6):

		finalTime = self.time + deltaT
		y0 = np.array([self.s['EGF'].amount, self.s['I'].amount, self.s['RE'].amount, self.s['Rd'].amount, self.s['Rp'].amount, self.s['GS'].amount, self.s['Rp_GS'].amount, self.s['Shc'].amount, self.s['Rp_Shc'].amount, self.s['Rp_pShc'].amount, self.s['pShc'].amount, self.s['Rp_pShc_GS'].amount, self.s['PI3K'].amount, self.s['Rp_PI3K'].amount, self.s['RasGAP'].amount, self.s['Rp_RasGAP'].amount, self.s['IRL'].amount, self.s['IRp'].amount, self.s['IRp_PI3K'].amount, self.s['IRp_RasGAP'].amount, self.s['IRS'].amount, self.s['IRp_IRS'].amount, self.s['IRp_IRSp'].amount, self.s['IRSp'].amount, self.s['iSrc'].amount, self.s['mIRS'].amount, self.s['mIRSp'].amount, self.s['mIRSp_GS'].amount, self.s['mIRSp_PI3K'].amount, self.s['SHP2'].amount, self.s['mIRSp_SHP2'].amount, self.s['GAB'].amount, self.s['mGAB'].amount, self.s['mGABp'].amount, self.s['mGABp_GS'].amount, self.s['mGABp_PI3K'].amount, self.s['mGABp_SHP2'].amount, self.s['mGABp_pSHP2'].amount, self.s['PIP3'].amount, self.s['dRas'].amount, self.s['Raf'].amount, self.s['aRaf'].amount, self.s['Mek'].amount, self.s['Erk'].amount, self.s['pErk'].amount, self.s['PDK1'].amount, self.s['Akt'].amount, self.s['pAkt'].amount, self.s['mTOR'].amount, self.s['Null'].amount, self.s['aaRaf'].amount, self.s['PKA'].amount, self.s['pShc_GS'].amount, self.s['ppMek'].amount, self.s['mGABp_pSHP2_GS'].amount, self.s['R'].amount, self.s['ppErk'].amount, self.s['IR'].amount, self.s['mPDK1'].amount, self.s['tRas'].amount, self.s['tRas_PI3K'].amount, self.s['ppAkt'].amount, self.s['mGABp_RasGAP'].amount, self.s['amTOR'].amount, self.s['iGS'].amount, self.s['imGAB'].amount, self.s['imIRS'].amount, self.s['aSrc'].amount, self.s['Ri'].amount, self.s['IRi'].amount, self.s['iPX'].amount, self.s['aPX'].amount, self.s['aPX_GS'].amount, self.s['IRSp_PI3K'].amount, self.s['IRSp_GS'].amount, self.s['IRSp_SHP2'].amount, self.s['GABp'].amount, self.s['GABp_PI3K'].amount, self.s['GABp_GS'].amount, self.s['GABp_RasGAP'].amount, self.s['GABp_SHP2'].amount, self.s['GABp_pSHP2'].amount, self.s['GABp_pSHP2_GS'].amount, self.s['imGABp'].amount, self.s['bRasGAP'].amount, self.s['phosphorylated_Akt'].amount], dtype = np.float64)
		self.s['EGF'].amount, self.s['I'].amount, self.s['RE'].amount, self.s['Rd'].amount, self.s['Rp'].amount, self.s['GS'].amount, self.s['Rp_GS'].amount, self.s['Shc'].amount, self.s['Rp_Shc'].amount, self.s['Rp_pShc'].amount, self.s['pShc'].amount, self.s['Rp_pShc_GS'].amount, self.s['PI3K'].amount, self.s['Rp_PI3K'].amount, self.s['RasGAP'].amount, self.s['Rp_RasGAP'].amount, self.s['IRL'].amount, self.s['IRp'].amount, self.s['IRp_PI3K'].amount, self.s['IRp_RasGAP'].amount, self.s['IRS'].amount, self.s['IRp_IRS'].amount, self.s['IRp_IRSp'].amount, self.s['IRSp'].amount, self.s['iSrc'].amount, self.s['mIRS'].amount, self.s['mIRSp'].amount, self.s['mIRSp_GS'].amount, self.s['mIRSp_PI3K'].amount, self.s['SHP2'].amount, self.s['mIRSp_SHP2'].amount, self.s['GAB'].amount, self.s['mGAB'].amount, self.s['mGABp'].amount, self.s['mGABp_GS'].amount, self.s['mGABp_PI3K'].amount, self.s['mGABp_SHP2'].amount, self.s['mGABp_pSHP2'].amount, self.s['PIP3'].amount, self.s['dRas'].amount, self.s['Raf'].amount, self.s['aRaf'].amount, self.s['Mek'].amount, self.s['Erk'].amount, self.s['pErk'].amount, self.s['PDK1'].amount, self.s['Akt'].amount, self.s['pAkt'].amount, self.s['mTOR'].amount, self.s['Null'].amount, self.s['aaRaf'].amount, self.s['PKA'].amount, self.s['pShc_GS'].amount, self.s['ppMek'].amount, self.s['mGABp_pSHP2_GS'].amount, self.s['R'].amount, self.s['ppErk'].amount, self.s['IR'].amount, self.s['mPDK1'].amount, self.s['tRas'].amount, self.s['tRas_PI3K'].amount, self.s['ppAkt'].amount, self.s['mGABp_RasGAP'].amount, self.s['amTOR'].amount, self.s['iGS'].amount, self.s['imGAB'].amount, self.s['imIRS'].amount, self.s['aSrc'].amount, self.s['Ri'].amount, self.s['IRi'].amount, self.s['iPX'].amount, self.s['aPX'].amount, self.s['aPX_GS'].amount, self.s['IRSp_PI3K'].amount, self.s['IRSp_GS'].amount, self.s['IRSp_SHP2'].amount, self.s['GABp'].amount, self.s['GABp_PI3K'].amount, self.s['GABp_GS'].amount, self.s['GABp_RasGAP'].amount, self.s['GABp_SHP2'].amount, self.s['GABp_pSHP2'].amount, self.s['GABp_pSHP2_GS'].amount, self.s['imGABp'].amount, self.s['bRasGAP'].amount, self.s['phosphorylated_Akt'].amount = odeint(self._SolveReactions, y0, [self.time, finalTime], atol = absoluteTolerance, rtol = relativeTolerance, mxstep=5000000)[-1]
		self.time = finalTime
		self.AssignmentRules()

class reaction_1:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k1'].value * self.parent.s['R'].concentration * self.parent.s['EGF'].concentration - self.parent.p['k_1'].value * self.parent.s['RE'].concentration) * self.parent.c['cell'].size

class reaction_2:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k2'].value * self.parent.s['RE'].concentration * self.parent.s['RE'].concentration - self.parent.p['k_2'].value * self.parent.s['Rd'].concentration) * self.parent.c['cell'].size

class reaction_3:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k3'].value * self.parent.s['Rd'].concentration * self.parent.c['cell'].size

class reaction_4:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k4'].value * self.parent.s['Rp'].concentration * self.parent.s['GS'].concentration - self.parent.p['k_4'].value * self.parent.s['Rp_GS'].concentration) * self.parent.c['cell'].size

class reaction_5:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k5'].value * self.parent.s['Rp'].concentration * self.parent.s['Shc'].concentration - self.parent.p['k_5'].value * self.parent.s['Rp_Shc'].concentration) * self.parent.c['cell'].size

class reaction_6:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k6'].value * self.parent.s['Rp_Shc'].concentration * self.parent.c['cell'].size

class reaction_7:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_7'].value * self.parent.s['Rp_pShc'].concentration - self.parent.p['k7'].value * self.parent.s['Rp'].concentration * self.parent.s['pShc'].concentration) * self.parent.c['cell'].size

class reaction_8:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V8'].value * self.parent.s['pShc'].concentration / (self.parent.p['Km8'].value + self.parent.s['pShc'].concentration)) * self.parent.c['cell'].size

class reaction_9:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k9'].value * self.parent.s['Rp_pShc'].concentration * self.parent.s['GS'].concentration - self.parent.p['k_9'].value * self.parent.s['Rp_pShc_GS'].concentration) * self.parent.c['cell'].size

class reaction_10:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_10'].value * self.parent.s['Rp_pShc_GS'].concentration - self.parent.p['k10'].value * self.parent.s['Rp'].concentration * self.parent.s['pShc_GS'].concentration) * self.parent.c['cell'].size

class reaction_11:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_11'].value * self.parent.s['pShc_GS'].concentration - self.parent.p['k11'].value * self.parent.s['pShc'].concentration * self.parent.s['GS'].concentration) * self.parent.c['cell'].size

class reaction_12:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k12'].value * self.parent.s['Rp'].concentration * self.parent.s['PI3K'].concentration - self.parent.p['k_12'].value * self.parent.s['Rp_PI3K'].concentration) * self.parent.c['cell'].size

class reaction_13:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k13'].value * self.parent.s['Rp'].concentration * self.parent.s['RasGAP'].concentration - self.parent.p['k_13'].value * self.parent.s['Rp_RasGAP'].concentration) * self.parent.c['cell'].size

class reaction_17:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k17'].value * self.parent.s['Rp'].concentration * self.parent.c['cell'].size

class reaction_18:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k17'].value * self.parent.s['Rp_GS'].concentration * self.parent.c['cell'].size

class reaction_19:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k17'].value * self.parent.s['Rp_Shc'].concentration * self.parent.c['cell'].size

class reaction_20:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k17'].value * self.parent.s['Rp_pShc'].concentration * self.parent.c['cell'].size

class reaction_21:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k17'].value * self.parent.s['Rp_pShc_GS'].concentration * self.parent.c['cell'].size

class reaction_22:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k17'].value * self.parent.s['Rp_PI3K'].concentration * self.parent.c['cell'].size

class reaction_23:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k17'].value * self.parent.s['Rp_RasGAP'].concentration * self.parent.c['cell'].size

class reaction_24:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k24'].value * self.parent.s['IR'].concentration * self.parent.s['I'].concentration - self.parent.p['k_24'].value * self.parent.s['IRL'].concentration) * self.parent.c['cell'].size

class reaction_25:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k25'].value * self.parent.s['IRL'].concentration * self.parent.c['cell'].size

class reaction_26:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k26'].value * self.parent.s['IRp'].concentration * self.parent.s['PI3K'].concentration - self.parent.p['k_26'].value * self.parent.s['IRp_PI3K'].concentration) * self.parent.c['cell'].size

class reaction_27:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k27'].value * self.parent.s['IRp'].concentration * self.parent.s['RasGAP'].concentration - self.parent.p['k_27'].value * self.parent.s['IRp_RasGAP'].concentration) * self.parent.c['cell'].size

class reaction_28:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k28'].value * self.parent.s['IRp'].concentration * self.parent.s['IRS'].concentration - self.parent.p['k_28'].value * self.parent.s['IRp_IRS'].concentration) * self.parent.c['cell'].size

class reaction_29:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k29'].value * self.parent.s['IRp_IRS'].concentration * self.parent.c['cell'].size

class reaction_30:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_30'].value * self.parent.s['IRp_IRSp'].concentration - self.parent.p['k30'].value * self.parent.s['IRp'].concentration * self.parent.s['IRSp'].concentration) * self.parent.c['cell'].size

class reaction_31:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V31'].value * self.parent.s['IRSp'].concentration / (self.parent.p['Km31'].value + self.parent.s['IRSp'].concentration)) * self.parent.c['cell'].size

class reaction_34:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k17'].value * self.parent.s['IRp'].concentration * self.parent.c['cell'].size

class reaction_35:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k17'].value * self.parent.s['IRp_PI3K'].concentration * self.parent.c['cell'].size

class reaction_36:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k17'].value * self.parent.s['IRp_RasGAP'].concentration * self.parent.c['cell'].size

class reaction_37:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k17'].value * self.parent.s['IRp_IRS'].concentration * self.parent.c['cell'].size

class reaction_38:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k17'].value * self.parent.s['IRp_IRSp'].concentration * self.parent.c['cell'].size

class reaction_40:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat40'].value * self.parent.s['iSrc'].concentration * (self.parent.s['Rp'].concentration + self.parent.p['alpha40'].value * self.parent.s['IRp'].concentration) / (self.parent.p['Km40'].value + self.parent.s['iSrc'].concentration)) * self.parent.c['cell'].size

class reaction_41:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V41'].value * self.parent.s['aSrc'].concentration / (self.parent.p['Km41'].value + self.parent.s['aSrc'].concentration)) * self.parent.c['cell'].size

class reaction_42:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k42'].value * self.parent.s['IRS'].concentration * self.parent.s['PIP3'].concentration - self.parent.p['k_42'].value * self.parent.s['mIRS'].concentration) * self.parent.c['cell'].size

class reaction_43:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat43'].value * self.parent.s['mIRS'].concentration * (self.parent.s['IRp'].concentration + self.parent.p['alpha43'].value * self.parent.s['Rp'].concentration) / (self.parent.p['Km43'].value + self.parent.s['mIRS'].concentration)) * self.parent.c['cell'].size

class reaction_44:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V31'].value * self.parent.s['mIRSp'].concentration / (self.parent.p['Km31'].value + self.parent.s['mIRSp'].concentration)) * self.parent.c['cell'].size

class reaction_45:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k45'].value * self.parent.s['mIRSp'].concentration * self.parent.s['GS'].concentration - self.parent.p['k_45'].value * self.parent.s['mIRSp_GS'].concentration) * self.parent.c['cell'].size

class reaction_46:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k46'].value * self.parent.s['mIRSp'].concentration * self.parent.s['PI3K'].concentration - self.parent.p['k_46'].value * self.parent.s['mIRSp_PI3K'].concentration) * self.parent.c['cell'].size

class reaction_47:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k47'].value * self.parent.s['mIRSp'].concentration * self.parent.s['SHP2'].concentration - self.parent.p['k_47'].value * self.parent.s['mIRSp_SHP2'].concentration) * self.parent.c['cell'].size

class reaction_48:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k48'].value * self.parent.s['mIRSp_SHP2'].concentration * self.parent.c['cell'].size

class reaction_49:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k49'].value * self.parent.s['GAB'].concentration * self.parent.s['PIP3'].concentration - self.parent.p['k_49'].value * self.parent.s['mGAB'].concentration) * self.parent.c['cell'].size

class reaction_50:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat50'].value * self.parent.s['mGAB'].concentration * (self.parent.s['Rp'].concentration + self.parent.p['alpha50'].value * self.parent.s['aSrc'].concentration) / (self.parent.p['Km50'].value + self.parent.s['mGAB'].concentration)) * self.parent.c['cell'].size

class reaction_51:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V51'].value * self.parent.s['mGABp'].concentration / (self.parent.p['Km51'].value + self.parent.s['mGABp'].concentration)) * self.parent.c['cell'].size

class reaction_52:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k52'].value * self.parent.s['mGABp'].concentration * self.parent.s['GS'].concentration - self.parent.p['k_52'].value * self.parent.s['mGABp_GS'].concentration) * self.parent.c['cell'].size

class reaction_53:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k53'].value * self.parent.s['mGABp'].concentration * self.parent.s['PI3K'].concentration - self.parent.p['k_53'].value * self.parent.s['mGABp_PI3K'].concentration) * self.parent.c['cell'].size

class reaction_54:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k54'].value * self.parent.s['mGABp'].concentration * self.parent.s['RasGAP'].concentration - self.parent.p['k_54'].value * self.parent.s['mGABp_RasGAP'].concentration) * self.parent.c['cell'].size

class reaction_55:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k55'].value * self.parent.s['mGABp'].concentration * self.parent.s['SHP2'].concentration - self.parent.p['k_55'].value * self.parent.s['mGABp_SHP2'].concentration) * self.parent.c['cell'].size

class reaction_56:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k56'].value * self.parent.s['mGABp_SHP2'].concentration * self.parent.c['cell'].size

class reaction_57:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat57'].value * self.parent.s['mGABp_SHP2'].concentration * (self.parent.s['Rp'].concentration + self.parent.s['aSrc'].concentration) / (self.parent.p['Km57'].value + self.parent.s['mGABp_SHP2'].concentration)) * self.parent.c['cell'].size

class reaction_58:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V58'].value * self.parent.s['mGABp_pSHP2'].concentration / (self.parent.p['Km58'].value + self.parent.s['mGABp_pSHP2'].concentration)) * self.parent.c['cell'].size

class reaction_59:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k59'].value * self.parent.s['mGABp_pSHP2'].concentration * self.parent.s['GS'].concentration - self.parent.p['k_59'].value * self.parent.s['mGABp_pSHP2_GS'].concentration) * self.parent.c['cell'].size

class reaction_60:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k60'].value * (self.parent.s['Rp_PI3K'].concentration + self.parent.s['IRp_PI3K'].concentration + self.parent.s['mIRSp_PI3K'].concentration + self.parent.s['mGABp_PI3K'].concentration + self.parent.s['tRas_PI3K'].concentration) * self.parent.c['cell'].size

class reaction_61:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k61'].value * self.parent.s['PIP3'].concentration * self.parent.c['cell'].size

class reaction_62:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat62'].value * self.parent.s['dRas'].concentration * (self.parent.s['Rp_GS'].concentration + self.parent.s['Rp_pShc_GS'].concentration + self.parent.s['mIRSp_GS'].concentration + self.parent.s['mGABp_GS'].concentration + self.parent.s['mGABp_pSHP2_GS'].concentration) / (self.parent.p['Km62'].value + self.parent.s['dRas'].concentration)) * self.parent.c['cell'].size

class reaction_63:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat63'].value * self.parent.s['tRas'].concentration * (self.parent.s['bRasGAP'].concentration + self.parent.s['mGABp_RasGAP'].concentration + self.parent.s['Rp_RasGAP'].concentration + self.parent.s['IRp_RasGAP'].concentration) / (self.parent.p['Km63'].value + self.parent.s['tRas'].concentration)) * self.parent.c['cell'].size

class reaction_64:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k64'].value * self.parent.s['tRas'].concentration * self.parent.s['PI3K'].concentration - self.parent.p['k_64'].value * self.parent.s['tRas_PI3K'].concentration) * self.parent.c['cell'].size

class reaction_65:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat65'].value * self.parent.s['tRas'].concentration * self.parent.s['Raf'].concentration / (self.parent.p['Km65'].value + self.parent.s['Raf'].concentration)) * self.parent.c['cell'].size

class reaction_66:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat66'].value * self.parent.s['aSrc'].concentration * self.parent.s['aRaf'].concentration / (self.parent.p['Km66'].value + self.parent.s['aRaf'].concentration)) * self.parent.c['cell'].size

class reaction_67:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat67'].value * self.parent.s['aaRaf'].concentration * (self.parent.s['PKA'].concentration / (self.parent.p['Km67'].value + self.parent.s['aaRaf'].concentration)) + self.parent.p['alpha67'].value * self.parent.s['aaRaf'].concentration * (self.parent.s['pAkt'].concentration + self.parent.p['beta67'].value * self.parent.s['ppAkt'].concentration)) * self.parent.c['cell'].size

class reaction_68:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat68'].value * self.parent.s['aaRaf'].concentration * self.parent.s['Mek'].concentration / (self.parent.p['Km68'].value + self.parent.s['Mek'].concentration)) * self.parent.c['cell'].size

class reaction_69:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V69'].value * self.parent.s['ppMek'].concentration / (self.parent.p['Km69'].value + self.parent.s['ppMek'].concentration)) * self.parent.c['cell'].size

class reaction_70:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat70'].value * self.parent.s['Erk'].concentration * self.parent.s['ppMek'].concentration / (self.parent.p['Km70'].value + self.parent.s['Erk'].concentration + self.parent.s['pErk'].concentration * (self.parent.p['Km70'].value / self.parent.p['Km71'].value))) * self.parent.c['cell'].size

class reaction_71:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat71'].value * self.parent.s['pErk'].concentration * self.parent.s['ppMek'].concentration / (self.parent.p['Km71'].value + self.parent.s['pErk'].concentration + self.parent.s['Erk'].concentration * (self.parent.p['Km71'].value / self.parent.p['Km70'].value))) * self.parent.c['cell'].size

class reaction_72:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V72'].value * self.parent.s['ppErk'].concentration / (self.parent.p['Km72'].value + self.parent.s['ppErk'].concentration + self.parent.s['pErk'].concentration * (self.parent.p['Km72'].value / self.parent.p['Km73'].value))) * self.parent.c['cell'].size

class reaction_73:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V73'].value * self.parent.s['pErk'].concentration / (self.parent.p['Km73'].value + self.parent.s['pErk'].concentration + self.parent.s['ppErk'].concentration * (self.parent.p['Km73'].value / self.parent.p['Km72'].value))) * self.parent.c['cell'].size

class reaction_74:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k74'].value * self.parent.s['PDK1'].concentration * self.parent.s['PIP3'].concentration - self.parent.p['k_74'].value * self.parent.s['mPDK1'].concentration) * self.parent.c['cell'].size

class reaction_75:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat75'].value * self.parent.s['mPDK1'].concentration * self.parent.s['Akt'].concentration / (self.parent.p['Km75'].value + self.parent.s['Akt'].concentration)) * self.parent.c['cell'].size

class reaction_76:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V76'].value * self.parent.s['pAkt'].concentration / (self.parent.p['Km76'].value + self.parent.s['pAkt'].concentration)) * self.parent.c['cell'].size

class reaction_77:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat77'].value * self.parent.s['mTOR'].concentration * ((self.parent.p['alpha77'].value * self.parent.s['pAkt'].concentration + self.parent.s['ppAkt'].concentration) / (self.parent.p['Km77'].value + self.parent.s['mTOR'].concentration)) - self.parent.p['k_77'].value * self.parent.s['amTOR'].concentration) * self.parent.c['cell'].size

class reaction_78:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat78'].value * self.parent.s['amTOR'].concentration * (self.parent.s['pAkt'].concentration / (self.parent.p['Km78'].value + self.parent.s['pAkt'].concentration)) - self.parent.p['k_78'].value * self.parent.s['ppAkt'].concentration) * self.parent.c['cell'].size

class reaction_79:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat79'].value * self.parent.s['ppErk'].concentration * (self.parent.s['GS'].concentration / (self.parent.p['Km79'].value + self.parent.s['GS'].concentration)) - self.parent.p['k_79'].value * self.parent.s['iGS'].concentration) * self.parent.c['cell'].size

class reaction_80:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat80'].value * self.parent.s['mGAB'].concentration * (self.parent.s['ppErk'].concentration / (self.parent.p['Km80'].value + self.parent.s['mGAB'].concentration)) - self.parent.p['k_80'].value * self.parent.s['imGAB'].concentration) * self.parent.c['cell'].size

class reaction_81:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['kcat81'].value * self.parent.s['mIRS'].concentration * (self.parent.s['amTOR'].concentration / (self.parent.p['Km81'].value + self.parent.s['mIRS'].concentration)) - self.parent.p['k_81'].value * self.parent.s['imIRS'].concentration) * self.parent.c['cell'].size

class reaction_82:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V_82'].value * self.parent.s['Rp'].concentration / (self.parent.p['Km82'].value + self.parent.s['Rp'].concentration)) * self.parent.c['cell'].size

class reaction_83:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k83'].value * self.parent.s['Ri'].concentration * self.parent.c['cell'].size

class reaction_84:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V_84'].value * self.parent.s['IRp'].concentration / (self.parent.p['Km84'].value + self.parent.s['IRp'].concentration)) * self.parent.c['cell'].size

class reaction_85:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k85'].value * self.parent.s['IRi'].concentration * self.parent.c['cell'].size

class reaction_88:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_42'].value * self.parent.s['mIRSp'].concentration - self.parent.p['k42'].value * self.parent.s['PIP3'].concentration * self.parent.s['IRSp'].concentration) * self.parent.c['cell'].size

class reaction_89:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_42'].value * self.parent.s['mIRSp_PI3K'].concentration - self.parent.p['k42'].value * self.parent.s['PIP3'].concentration * self.parent.s['IRSp_PI3K'].concentration) * self.parent.c['cell'].size

class reaction_90:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_42'].value * self.parent.s['mIRSp_GS'].concentration - self.parent.p['k42'].value * self.parent.s['PIP3'].concentration * self.parent.s['IRSp_GS'].concentration) * self.parent.c['cell'].size

class reaction_91:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_42'].value * self.parent.s['mIRSp_SHP2'].concentration - self.parent.p['k42'].value * self.parent.s['PIP3'].concentration * self.parent.s['IRSp_SHP2'].concentration) * self.parent.c['cell'].size

class reaction_92:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_42'].value * self.parent.s['mGABp'].concentration - self.parent.p['k42'].value * self.parent.s['PIP3'].concentration * self.parent.s['GABp'].concentration) * self.parent.c['cell'].size

class reaction_93:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_42'].value * self.parent.s['mGABp_PI3K'].concentration - self.parent.p['k42'].value * self.parent.s['PIP3'].concentration * self.parent.s['GABp_PI3K'].concentration) * self.parent.c['cell'].size

class reaction_94:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_42'].value * self.parent.s['mGABp_GS'].concentration - self.parent.p['k42'].value * self.parent.s['PIP3'].concentration * self.parent.s['GABp_GS'].concentration) * self.parent.c['cell'].size

class reaction_95:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_42'].value * self.parent.s['mGABp_RasGAP'].concentration - self.parent.p['k49'].value * self.parent.s['PIP3'].concentration * self.parent.s['GABp_RasGAP'].concentration) * self.parent.c['cell'].size

class reaction_96:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_42'].value * self.parent.s['mGABp_SHP2'].concentration - self.parent.p['k42'].value * self.parent.s['PIP3'].concentration * self.parent.s['GABp_SHP2'].concentration) * self.parent.c['cell'].size

class reaction_97:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_42'].value * self.parent.s['mGABp_pSHP2'].concentration - self.parent.p['k42'].value * self.parent.s['PIP3'].concentration * self.parent.s['GABp_pSHP2'].concentration) * self.parent.c['cell'].size

class reaction_98:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['k_42'].value * self.parent.s['mGABp_pSHP2_GS'].concentration - self.parent.p['k42'].value * self.parent.s['PIP3'].concentration * self.parent.s['GABp_pSHP2_GS'].concentration) * self.parent.c['cell'].size

class reaction_99:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V31'].value * self.parent.s['IRSp_PI3K'].concentration / (self.parent.p['Km31'].value + self.parent.s['IRSp_PI3K'].concentration)) * self.parent.c['cell'].size

class reaction_100:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V31'].value * self.parent.s['IRSp_GS'].concentration / (self.parent.p['Km31'].value + self.parent.s['IRSp_GS'].concentration)) * self.parent.c['cell'].size

class reaction_101:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k48'].value * self.parent.s['IRSp_SHP2'].concentration * self.parent.c['cell'].size

class reaction_102:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k56'].value * self.parent.s['mGABp_pSHP2'].concentration * self.parent.c['cell'].size

class reaction_103:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k56'].value * self.parent.s['mGABp_pSHP2_GS'].concentration * self.parent.c['cell'].size

class reaction_104:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V51'].value * self.parent.s['GABp'].concentration / (self.parent.p['Km51'].value + self.parent.s['GABp'].concentration)) * self.parent.c['cell'].size

class reaction_105:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V51'].value * self.parent.s['GABp_PI3K'].concentration / (self.parent.p['Km51'].value + self.parent.s['GABp_PI3K'].concentration)) * self.parent.c['cell'].size

class reaction_106:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V51'].value * self.parent.s['GABp_GS'].concentration / (self.parent.p['Km51'].value + self.parent.s['GABp_GS'].concentration)) * self.parent.c['cell'].size

class reaction_107:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (self.parent.p['V51'].value * self.parent.s['GABp_RasGAP'].concentration / (self.parent.p['Km51'].value + self.parent.s['GABp_RasGAP'].concentration)) * self.parent.c['cell'].size

class reaction_108:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k56'].value * self.parent.s['GABp_SHP2'].concentration * self.parent.c['cell'].size

class reaction_109:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k56'].value * self.parent.s['GABp_pSHP2'].concentration * self.parent.c['cell'].size

class reaction_110:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k56'].value * self.parent.s['GABp_pSHP2_GS'].concentration * self.parent.c['cell'].size

class reaction_111:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k111'].value * (self.parent.s['mGABp_SHP2'].concentration + self.parent.s['mGABp_pSHP2'].concentration + self.parent.s['mGABp_pSHP2_GS'].concentration + self.parent.s['mIRSp_SHP2'].concentration) * self.parent.s['mGABp_RasGAP'].concentration * self.parent.c['cell'].size

class reaction_112:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k111'].value * (self.parent.s['mGABp_SHP2'].concentration + self.parent.s['mGABp_pSHP2'].concentration + self.parent.s['mGABp_pSHP2_GS'].concentration) * self.parent.s['Rp_RasGAP'].concentration * self.parent.c['cell'].size

class reaction_113:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k111'].value * (self.parent.s['mGABp_SHP2'].concentration + self.parent.s['mGABp_pSHP2'].concentration + self.parent.s['mGABp_pSHP2_GS'].concentration) * self.parent.s['IRp_RasGAP'].concentration * self.parent.c['cell'].size

class reaction_114:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k111'].value * self.parent.s['mIRSp_SHP2'].concentration * self.parent.s['Rp_RasGAP'].concentration * self.parent.c['cell'].size

class reaction_115:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k111'].value * self.parent.s['mIRSp_SHP2'].concentration * self.parent.s['IRp_RasGAP'].concentration * self.parent.c['cell'].size

class reaction_117:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return (2 * self.parent.p['kcat80'].value * self.parent.s['mGABp'].concentration * self.parent.s['ppErk'].concentration / (self.parent.p['Km80'].value + self.parent.s['mGABp'].concentration) - self.parent.p['k_80'].value * self.parent.s['imGABp'].concentration) * self.parent.c['cell'].size

class reaction_118:

	def __init__(self, parent, metadata = None):

		self.parent = parent
		self.p = {}
		self.metadata = metadata


	def __call__(self):
		return self.parent.p['k118'].value * self.parent.s['imGABp'].concentration * self.parent.c['cell'].size

