"""Density Test

This script tests the scale of the saturation pressure and acts as a primary comparison tool between the implemented
methods.
"""

import json
import sys

sys.path.append("/home/flavian/Desktop/Bachelor's End Project/adsorpyon-0.3.1/adsorpyon/")
from adsorpyon import density

import numpy as np
import matplotlib.pyplot as plt

AXES_SIZE = 19
TICK_SIZE = 13
LEGEND_SIZE = 14
FIGURE_SIZE = (7, 6)

adsorbate_data = json.load(open("../adsorpyon/Properties/N2.prop"))

temperatures = np.array(np.linspace(start=adsorbate_data["temperature_boiling"], stop=400, num=100))

densities = {'OzawaC': None, 'Empirical': None, 'HauerV': None, "OzawaV": None, "HauerC": None}

for method in densities:
    densities[method] = []
    for temperature in temperatures:
        if method == 'OzawaC':
            densities[method].append(density.ozawa(temperature=temperature,
                                                   temperature_boiling=adsorbate_data["temperature_boiling"],
                                                   density_boiling=adsorbate_data["density_boiling"]))
        elif method == 'Empirical':
            densities[method].append(density.empirical(molecular_mass=adsorbate_data["molecular_mass"],
                                                       temperature_critical=adsorbate_data["temperature_critical"],
                                                       pressure_critical=adsorbate_data["pressure_critical"]))
        elif method == 'HauerV':
            densities[method].append(density.hauer(temperature=temperature,
                                                   temperature_boiling=adsorbate_data["temperature_boiling"],
                                                   density_boiling=adsorbate_data["density_boiling"],
                                                   thermal_expansion_coefficient=adsorbate_data["thermal_expansion_coefficient"]))
        elif method == "OzawaV":
            densities[method].append(density.ozawa_modified(temperature=temperature,
                                                            temperature_boiling=adsorbate_data["temperature_boiling"],
                                                            density_boiling=adsorbate_data["density_boiling"],
                                                            thermal_expansion_coefficient=adsorbate_data["thermal_expansion_coefficient"]))
        elif method == 'HauerC':
            densities[method].append(density.hauer(temperature=temperature,
                                                   temperature_boiling=adsorbate_data["temperature_boiling"],
                                                   density_boiling=adsorbate_data["density_boiling"],
                                                   thermal_expansion_coefficient=0.00165))

plt.figure(figsize=FIGURE_SIZE)
plt.rc('axes', labelsize=AXES_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=TICK_SIZE)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=TICK_SIZE)  # fontsize of the y tick labels
plt.rc('legend', fontsize=LEGEND_SIZE)  # fontsize of the legend
plt.plot(temperatures, densities['Empirical'], label="Empirical")
plt.plot(temperatures, densities['OzawaC'], label="Ozawa (const. $\\alpha$)")
plt.plot(temperatures, densities['HauerV'], label="Hauer (var. $\\alpha$)")
plt.plot(temperatures, densities['OzawaV'], label="Ozawa (var. $\\alpha$)")
plt.plot(temperatures, densities['HauerC'], label="Hauer (const. $\\alpha$)")
plt.axvline(adsorbate_data["temperature_critical"], linestyle="--")
plt.xlabel("Temperature [K]")
plt.ylabel("Density [kg/$m^{3}$]")
plt.legend()
plt.show()
