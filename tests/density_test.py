"""Density Test

This script tests the scale of the saturation pressure and acts as a primary comparison tool between the implemented
methods.
"""

import json

from adsorpyon import density

import numpy as np
import matplotlib.pyplot as plt


adsorbate_data = json.load(open("../adsorpyon/Properties/CO2.prop"))

temperatures = np.array(np.linspace(start=240, stop=400, num=200))

densities = {'Ozawa': None, 'Empirical': None, 'Hauer': None}

for method in densities:
    densities[method] = []
    for temperature in temperatures:
        if method == 'Ozawa':
            densities[method].append(density.ozawa(temperature=temperature,
                                                   temperature_boiling=adsorbate_data["temperature_boiling"],
                                                   density_boiling=adsorbate_data["density_boiling"]))
        elif method == 'Empirical':
            densities[method].append(density.empirical(molecular_mass=adsorbate_data["molecular_mass"],
                                                       temperature_critical=adsorbate_data["temperature_critical"],
                                                       pressure_critical=adsorbate_data["pressure_critical"]))
        #elif method == 'Hauer':
            #densities[method].append(density.hauer(temperature=temperature,
            #                                       temperature_reference=adsorbate_data["temperature_reference"],
            #                                       density_boiling=adsorbate_data["density_boiling"]))


#print(f"PREOS               Dubinin             Extrapolation           PRSV")
#for p1, p2, p3 in zip(saturation_pressure['PREOS'], saturation_pressure['Dubinin'], saturation_pressure['Extrapolation']):
#    print(f"{p1}              {p2}              {p3}")

plt.figure()
plt.plot(temperatures, densities['Ozawa'], label="Ozawa")
plt.plot(temperatures, densities['Empirical'], label="Empirical")
#plt.plot(temperatures, densities['Hauer'], label="Hauer")
plt.xlabel("Temperature [K]")
plt.ylabel("Density [kg/$m^{3}$]")
plt.legend()
plt.show()
