# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 07:59:59 2019

@author: Steve
"""

from sbmltoodepy.parse import ParseSBMLFile
from sbmltoodepy.modulegeneration import GenerateModel
import sys
import os
import sbmltoodepy.dataclasses

def ParseAndCreateModel(inputFilePath, jsonFilePath = None, outputFilePath = None, className = "SBMLmodel"):
    """
    
    This function parses an SBML model file and generates a Python file implementing the model.
    
    Parameters
    ----------
    inputFilePath : str
        Location of the SBML model to be parsed
    jsonFilePath : str, optional
        If provided, a json file containing the data for each of the model components will be created.
        If a file exists at the provided location, the file will be overwritten.
        Otherwise, this functionality is skipped and just a Python file will be generated.
    outputFilePath : str, optional
        Location where the Python file will be created.
        If a file exists at the provided location, the file will be overwritten.
        If a location is not provided, the file path of the generated python file will be based on the inputFilePath
    className : str, optional
        Name of the class that implements the SBML model.
		
    Returns
    -------
    None

    Warnings
    --------
    This function will overwrite files located at jsonFilePath and outputFilePath
    
    Notes
    -----
    ParseAndCreateModel is intended to combine most of the functionality provided by this package.
    The intent is for this function to be suitable for most applications.
    
    """
#    if jsonFileName == None:
#        jsonFileName = inputFileName.split('.')[0] + '.json'
        
    if outputFilePath == None:
        if len(inputFilePath.split('.')) > 1:
            
            if len(inputFilePath.split('.')) == 2 and inputFilePath[0] == '.':
                #explicit relative path with no extension
                outputFilePath = inputFilePath + '.py'
            else:
                outputFilePath = ''
                for portion in inputFilePath.split('.')[0:-1]:
                    outputFilePath += portion + '.'
                outputFilePath += 'py'
        else:
            #implicit file path with no extension
            outputFilePath = inputFilePath + '.py'
        
    modelData = ParseSBMLFile(inputFilePath)
    if not jsonFilePath == None: 
        modelData.DumpToJSON(jsonFilePath)
    GenerateModel(modelData, outputFilePath, objectName = className)
    
def PrintConcentrations(model):
    [print(key + ' Concentration: ' + str(model.s[key].concentration)) for key in model.s]
     
def PrintAmounts(model):
    [print(key + ' Amount: ' + str(model.s[key].amount)) for key in model.s]
    
def TestPackage(method="LSODA"):
    
    """
    
    This function serves to test the package.
    
    
    
    Notes
    -----
    The function raises a warning about trying to set a constant species and 6 numbers.
    
    1.27949533562e-06
    
    9.42747177955e-10
    
    1.33080766722e-07
    
    9.79758065656e-08
    
    4.05696056523e-07
    
    2.66117319017e-05
    
    These are the average relative errors for species, parameters, and compartments between the model results generated using SBMLtoODEpy and results calculated by COPASI
    for six different models provided in the sbmltoodepy/sbml_files subdirectory: Smallbone2013, Borisov2009, Guyton1972, Kerkhoven2013, Waugh2006, and Zi2011, respectively.
    
    This function is intended as to serve as an example of using the sample models. 

    """

    dirname = os.path.dirname(__file__)
    ParseAndCreateModel(os.path.join(dirname, 'sbml_files/Smallbone2013_Colon_Crypt_cycle.xml'), outputFilePath = os.path.join(dirname, 'sbml_files/Smallbone2013.py'), className = 'Smallbone2013')
#    print(dirname)
    from .sbml_files.Smallbone2013 import Smallbone2013
    
    model = Smallbone2013()
    model.RunSimulation(1.0, method=method)
    
    smallBoneSpeciesConcentrations = [model.s['N0'].amount, model.s['N1'].amount, model.s['N2'].amount]
    smallBoneActualConcentrations = [1.75445, 27.4059, 45.6191]
    smallBoneIndividualRelativeError = [abs(smallBoneSpeciesConcentrations[i] - smallBoneActualConcentrations[i])/abs(smallBoneActualConcentrations[i]) for i in range(len(smallBoneActualConcentrations))]
    smallBoneAverageRelativeError = sum(smallBoneIndividualRelativeError)/len(smallBoneIndividualRelativeError)
    
    ParseAndCreateModel(os.path.join(dirname, 'sbml_files/Borisov2009_insulin_EGF.xml'), 
                        outputFilePath = os.path.join(dirname, 'sbml_files/Borisov2009.py'),
                        className = 'Borisov2009')
    
    from .sbml_files.Borisov2009 import Borisov2009
    
    model = Borisov2009()
    model.RunSimulation(1.0, method=method)
    
    borisovSpeciesConcentrations = [model.s['PIP3'].concentration, model.s['GS'].concentration, model.s['RasGAP'].concentration, model.s['Rp'].concentration,
                                    model.s['SHP2'].concentration, model.s['PI3K'].concentration, model.s['IRp'].concentration, 
                                    model.s['mGABp'].concentration, model.s['GAB'].concentration, model.s['IRS'].concentration, 
                                    model.s['mIRSp'].concentration, model.s['mGAB'].concentration, model.s['RE'].concentration, 
                                    model.s['mGABp_pSHP2'].concentration, model.s['pErk'].concentration, model.s['pShc'].concentration, 
                                    model.s['mIRS'].concentration, model.s['IRSp'].concentration, model.s['mGABp_SHP2'].concentration, 
                                    model.s['pAkt'].concentration, model.s['Rp_Shc'].concentration, model.s['tRas'].concentration, 
                                    model.s['IRL'].concentration, model.s['Rp_pShc_GS'].concentration, model.s['Null'].concentration, 
                                    model.s['IRp_IRS'].concentration, model.s['mGABp_RasGAP'].concentration, model.s['Rd'].concentration, 
                                    model.s['Raf'].concentration, model.s['Mek'].concentration, model.s['iSrc'].concentration, 
                                    model.s['Rp_pShc'].concentration, model.s['mGABp_pSHP2_GS'].concentration, model.s['imGABp'].concentration, 
                                    model.s['mIRSp_SHP2'].concentration, model.s['Rp_RasGAP'].concentration, model.s['GABp_GS'].concentration, 
                                    model.s['IRSp_PI3K'].concentration, model.s['aRaf'].concentration, model.s['mGABp_PI3K'].concentration, 
                                    model.s['mIRSp_GS'].concentration, model.s['Rp_GS'].concentration, model.s['IRi'].concentration, 
                                    model.s['GABp'].concentration, model.s['Rp_PI3K'].concentration, model.s['GABp_pSHP2'].concentration, 
                                    model.s['Ri'].concentration, model.s['IRp_PI3K'].concentration, model.s['ppErk'].concentration, 
                                    model.s['amTOR'].concentration, model.s['GABp_SHP2'].concentration, model.s['IRp_IRSp'].concentration, 
                                    model.s['mGABp_GS'].concentration, model.s['mPDK1'].concentration, model.s['iGS'].concentration, 
                                    model.s['R'].concentration, model.s['Akt'].concentration, model.s['GABp_RasGAP'].concentration, 
                                    model.s['imIRS'].concentration, model.s['pShc_GS'].concentration, model.s['IR'].concentration, 
                                    model.s['mIRSp_PI3K'].concentration, model.s['dRas'].concentration, model.s['imGAB'].concentration,
                                    model.s['GABp_pSHP2_GS'].concentration, model.s['Shc'].concentration, model.s['IRSp_GS'].concentration,
                                    model.s['GABp_PI3K'].concentration, model.s['IRp_RasGAP'].concentration, model.s['IRSp_SHP2'].concentration,
                                    model.s['ppAkt'].concentration, model.s['aaRaf'].concentration, model.s['ppMek'].concentration,
                                    model.s['aSrc'].concentration, model.s['tRas_PI3K'].concentration, model.s['I'].concentration,
                                    model.s['EGF'].concentration, model.s['PDK1'].concentration, model.s['Erk'].concentration,
                                    model.s['mTOR'].concentration, model.s['phosphorylated_Akt'].concentration, model.p['EGF_tot'].value,
                                    model.p['k11'].value, model.p['k_1'].value, model.p['k_2'].value,
                                    model.p['k_4'].value, model.p['k_5'].value, model.p['k_7'].value,
                                    model.p['k_9'].value, model.p['k_10'].value, model.p['k_11'].value,
                                    model.p['k_12'].value, model.p['k_13'].value, model.p['k_24'].value,
                                    model.p['k_26'].value, model.p['k_27'].value, model.p['k_28'].value,
                                    model.p['k_30'].value, model.p['k_42'].value, model.p['k_45'].value,
                                    model.p['k_46'].value, model.p['k_47'].value, model.p['k_49'].value,
                                    model.p['k_52'].value, model.p['k_53'].value, model.p['k_54'].value, 
                                    model.p['k_55'].value, model.p['k_59'].value, model.p['k_74'].value]
    borisovActualConcentrations = [0.00219730259, 199.999532, 49.9999991, 0.01053148607, 300, 199.9950689, 0, 9.093311137e-7, 224.9999411, 299.9992224, 6.960582322e-9, 5.768731984e-5, 5.25093293, 2.133057558e-14, 4.217439485e-14, 1.422983082e-5, 0.0007775824175, 3.722902166e-11, 1.665153817e-8, 8.224665812e-6, 0.008874610276, 0.0003445448918, 0, 0.0001031382319, 0.001309333595, 0, 4.190009525e-11, 0.2831115665, 99.99999892, 200, 517.9856616, 0.0004692834851, 2.579596703e-15, 1.256796525e-41, 1.198920022e-10, 9.000910304e-7, 1.721357033e-10, 4.35717765e-12, 1.076050031e-6, 2.352898165e-7, 1.220301636e-11, 0.0003596555963, 0, 4.834492176e-9, 0.004930853535, 8.721754773e-17, 0.007167556864, 0, 4.457520658e-30, 1.809930138e-7, 9.31127658e-11, 0, 3.547926633e-8, 0.0002379007512, 2.851107509e-34, 94.11796688, 99.99999178, 2.041802173e-13, 5.70618847e-15, 5.181092861e-6, 150, 8.918042238e-10, 149.9996555, 4.389736732e-40, 9.879025455e-18, 269.9905336, 6.443885964e-14, 1.14181038e-9, 0, 6.660741571e-13, 6.831679439e-16, 4.52159142e-10, 3.861263654e-12, 0.01433838411, 0, 0, 0.8269990259, 99.9997621, 400, 99.99999982, 8.224665812e-6, 0.9999999386, 0.00666, 0.04000032, 0.495, 0.00666, 0.133, 0.2664, 0.0666, 0.16, 0.0666, 0.1161585, 0.001332, 0.000333002664, 1.161585, 0.1332, 0.39975, 0.066, 0.0666, 66.6, 0.00666, 0.666, 0.000666, 0.002, 0.03325, 0.66666, 0.0666, 0.2, 0.666]
    borisovIndividualRelativeError = [abs(borisovSpeciesConcentrations[i] - borisovActualConcentrations[i])/abs(borisovActualConcentrations[i] + sys.float_info.epsilon) for i in range(len(borisovActualConcentrations))]
    borisovAverageRelativeError = sum(borisovIndividualRelativeError)/len(borisovIndividualRelativeError)
    
    ParseAndCreateModel(os.path.join(dirname, 'sbml_files/Guyton1972_Angiotensin.xml'), outputFilePath = os.path.join(dirname, 'sbml_files/Guyton1972.py'), className = 'Guyton1972')

    from .sbml_files.Guyton1972 import Guyton1972

    model = Guyton1972()
    model.RunSimulation(1.0, method=method)
    
    guytonSpeciesConcentrations = [model.p['ANX1'].value, model.p['ANC'].value, model.p['ANGSCR'].value, model.p['MDFLW3'].value, model.p['ANX'].value, model.p['ANPR'].value, model.p['ANPRT'].value, model.p['ANPR1'].value, model.p['ANM'].value, model.p['ANU'].value, model.p['ANU1'].value, model.p['ANUVN'].value]
    guytonActualConcentrations = [0, 0.8678797002, 0.9645806004, 1.00051, 0, 0.9645806004, 0.9645806004, 0.9645806004, 0.9883008802, 0.9298052815, 0.9298052815, 1]
    guytonIndividualRelativeError = [abs(guytonSpeciesConcentrations[i] - guytonActualConcentrations[i])/abs(guytonActualConcentrations[i] + sys.float_info.epsilon) for i in range(len(guytonActualConcentrations))]
    guytonAverageRelativeError = sum(guytonIndividualRelativeError)/len(guytonIndividualRelativeError)     

    ParseAndCreateModel(os.path.join(dirname, 'sbml_files/Kerkhoven2013_Glycolysis_in_T_brucei.xml'), outputFilePath = os.path.join(dirname, 'sbml_files/Kerkhoven2013.py'), className = 'Kerkhoven2013')

    from .sbml_files.Kerkhoven2013 import Kerkhoven2013

    model = Kerkhoven2013()
    model.RunSimulation(1.0, method=method)
    
    kerkhovenSpeciesConcentrations = [model.s['ADP_g'].concentration, model.s['ADP_c'].concentration, model.s['DHAP_g'].concentration, model.s['GA3P_g'].concentration, model.s['Glc_c'].concentration, model.s['_2PGA_c'].concentration, model.s['Glc6P_g'].concentration, model.s['_3PGA_g'].concentration, model.s['Gly3P_g'].concentration, model.s['Pyr_c'].concentration, model.s['DHAP_c'].concentration, model.s['Fru6P_g'].concentration, model.s['NAD_g'].concentration, model.s['_3PGA_c'].concentration, model.s['Glc_g'].concentration, model.s['PEP_c'].concentration, model.s['_13BPGA_g'].concentration, model.s['ATP_g'].concentration, model.s['AMP_c'].concentration, model.s['ATP_c'].concentration, model.s['AMP_g'].concentration, model.s['Fru16BP_g'].concentration, model.s['Gly3P_c'].concentration, model.s['NADH_g'].concentration]
    kerkhovenActualConcentrations = [1.520133775, 0.9624601923, 2.803270922, 0.1110184998, 0.008428687709, 0.5386705367, 0.3937007366, 11.08505023, 1.274104358, 3.088012616, 3.507803671, 0.132241083, 3.926075328, 10.37790254, 0.008059459187, 1.334530214, 0.01380450523, 4.216834927, 0.1602003523, 2.777339455, 0.2630312977, 15.15902813, 1.492196329, 0.0739246723]
    kerkhovenIndividualRelativeError = [abs(kerkhovenSpeciesConcentrations[i] - kerkhovenActualConcentrations[i])/abs(kerkhovenActualConcentrations[i]) for i in range(len(kerkhovenActualConcentrations))]
    kerkhovenAverageRelativeError = sum(kerkhovenIndividualRelativeError)/len(kerkhovenIndividualRelativeError)     
    

    ParseAndCreateModel(os.path.join(dirname, 'sbml_files/Waugh2006_Diabetic_Wound_Healing_TGF_B_Dynamics.xml'), outputFilePath = os.path.join(dirname, 'sbml_files/Waugh2006.py'), className = 'Waugh2006')

    from .sbml_files.Waugh2006 import Waugh2006

    model = Waugh2006()
    model.RunSimulation(1.0, method=method)
    
    waughSpeciesConcentrations = [model.s['phi_I'].concentration, model.s['phi_R'].concentration, model.s['T'].concentration, model.s['K_T'].concentration]
    waughActualConcentrations = [232.9131481, 181.7337441, 1.783740715, 68.97275873]
    waughIndividualRelativeError = [abs(waughSpeciesConcentrations[i] - waughActualConcentrations[i])/abs(waughActualConcentrations[i]) for i in range(len(waughActualConcentrations))]
    waughAverageRelativeError = sum(waughIndividualRelativeError)/len(waughIndividualRelativeError)     
    
    ParseAndCreateModel(os.path.join(dirname, 'sbml_files/Zi2011_TGF_beta_Pathway.xml'), outputFilePath = os.path.join(dirname, 'sbml_files/Zi2011.py'), className = 'Zi2011')

    from .sbml_files.Zi2011 import Zi2011

    model = Zi2011()
    model.RunSimulation(1.0, method=method, absoluteTolerance = 1e-16, relativeTolerance = 1e-16)
    
    ziSpeciesConcentrations = [model.s['PSmad2c'].concentration, model.s['PSmad2n'].concentration, model.s['T1R_surf'].concentration, model.s['T2R_endo'].concentration, model.s['Smad2c'].concentration, model.s['Smad4n'].concentration, model.s['T2R_surf'].concentration, model.s['LRC_endo'].concentration, model.s['T1R_endo'].concentration, model.s['LRC_surf'].concentration, model.s['PSmad2_Smad4_c'].concentration, model.s['PSmad2_PSmad2_n'].concentration, model.s['TGF_beta_ex'].concentration, model.s['TGF_beta_endo'].concentration, model.s['Smad2n'].concentration, model.s['PSmad2_Smad4_n'].concentration, model.s['TGF_beta_ns'].concentration, model.s['Smad4c'].concentration, model.s['PSmad2_PSmad2_c'].concentration, model.c['Vmed'].size, model.p['totalNumT1R'].value, model.p['totalNumT2R'].value, model.p['totalNumLRC'].value, model.p['totalNumPSmad2'].value, model.p['totalNuclearPSmad2'].value, model.p['totalSmad2c'].value, model.p['totalSmad2n'].value, model.p['medium_TGF_beta_amount'].value, model.p['koff_ns'].value]
    ziActualConcentrations = [0.01391431482, 0.001607867046, 0.5255621489, 1.395985545, 60.54721489, 50.78878884, 0.02414685438, 0.04474184166, 6.479258, 0.1763515728, 0.03276905664, 2.141642116e-6, 0.0486999112, 0.0007455529254, 28.4939252, 0.01846754609, 0.001057648159, 50.76407599, 6.327163318e-6, 1.999484006e-9, 10004.99992, 2272.441261, 306.1259416, 76.74329454, 0.02007969642, 60.59391092, 28.51400489, 58424.81612, 2.033059171]
    ziIndividualRelativeError = [abs(ziSpeciesConcentrations[i] - ziActualConcentrations[i])/abs(ziActualConcentrations[i]) for i in range(len(ziActualConcentrations))]
    ziAverageRelativeError = sum(ziIndividualRelativeError)/len(ziIndividualRelativeError)     

    print(smallBoneAverageRelativeError)
    print(borisovAverageRelativeError)
    print(guytonAverageRelativeError)
    print(kerkhovenAverageRelativeError)
    print(waughAverageRelativeError)
    print(ziAverageRelativeError)
	