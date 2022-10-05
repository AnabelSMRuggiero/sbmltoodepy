import matplotlib.pyplot as plt
import numpy as np

from sbmltoodepy.parse import ParseSBMLFile
from sbmltoodepy.modulegeneration import GenerateModel

modelData = ParseSBMLFile('sbmltoodepy/sbml_files/BIOMD0000000007_url.xml')
GenerateModel(modelData, outputFilePath = 'sbmltoodepy/sbml_files/BIOMD0000000007.py', objectName = 'BIOMD0000000007')
from sbmltoodepy.sbml_files.BIOMD0000000007 import BIOMD0000000007
model = BIOMD0000000007()

n_species = len(model.s.keys())
n_secs = 400
deltaT = 0.1
n_steps = int(n_secs/deltaT)
abs_tol = 1e-6
rel_tol = 1e-12

def simulate(modelData, model, n_steps, abs_tol=1e-3, rel_tol=1e-3, time_interval=0.01):
    """
    Performs a rollout of the model.

    Returns
    -------

    """
    reactants = []
    for k, v in modelData.reactions.items():
        for (coef, reactant) in v.reactants:
            if reactant not in reactants:
                reactants.append(reactant)
    for k, v in modelData.rateRules.items():
        if v.variable not in reactants:
            reactants.append(v.variable)
    for k in model.s.keys():
        if k not in reactants:
            reactants.append(k)


    n_reactants = len(reactants)
    times = np.zeros(n_steps)
    amounts = np.zeros((n_reactants, n_steps))

    for i in range(n_steps):
        model.RunSimulation(
            time_interval, absoluteTolerance=abs_tol, relativeTolerance=rel_tol
        )
        times[i] = model.time
        for k_idx, k in enumerate(reactants):
            if k in model.s.keys():
                amounts[k_idx, i] = model.s[k].amount

            elif k in model.p.keys():
                amounts[k_idx, i] = model.p[k].value

            elif k in model.c.keys():
                amounts[k_idx, i] = model.c[k].size

            else:
                raise ValueError

    return reactants, amounts, times

reactants, amounts, times = simulate(modelData, model, n_steps, abs_tol=abs_tol, rel_tol=rel_tol, time_interval=deltaT)


fig, ax = plt.subplots(3, 1)
for k_idx, k in enumerate(reactants):
    if k in ["G2K", "PG2", "Rum1Total", "Cdc13Total"]:
        ax[0].plot(times, amounts[k_idx], label=k)
    elif k in ["Mass"]:
        ax[1].plot(times, amounts[k_idx], label=k)
    elif k in ["Cig2Total", "G1K"]:
        ax[2].plot(times, amounts[k_idx], label=k)
ax[0].legend()
ax[1].legend()
ax[2].legend()
plt.show()