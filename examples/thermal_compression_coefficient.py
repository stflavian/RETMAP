"""Thermal Compression Coefficient Calculator

This script is used to calculate the thermal compression coefficient of a material based on experimental data.
"""

import json
import sys
sys.path.append("/home/flavian/Desktop/Bachelor's End Project/adsorpyon-0.3.1/adsorpyon/")

from adsorpyon import density

import numpy
import matplotlib.pyplot as plt
import scipy.optimize

MATERIAL = "N2"
MOF_FILES = [f"sim-{MATERIAL}-Co-MOF-74-cCO2-saturation_T.load",
             f"sim-{MATERIAL}-IRMOF-1-saturation_T.load",
             f"sim-{MATERIAL}-MIL-47-saturation_T.load",
             f"sim-{MATERIAL}-MOF-1-saturation_T.load",
             f"sim-{MATERIAL}-ZJU-198-LP-saturation_T.load"]

NAMES = ["Co-MOF-74", "IRMOF-1", "MIL-47", "MOF-1", "ZJU-198-LP"]

temperature_range = numpy.linspace(start=253, stop=393, num=15, endpoint=True)

with open(f"../adsorpyon/Properties/N2.prop") as adsorbate_file:
    adsorbate_data = json.load(adsorbate_file)

mof_loading = [None for _ in MOF_FILES]
popt = [None for _ in MOF_FILES]

for index, file in enumerate(MOF_FILES):
    data = numpy.genfromtxt(file)[:, 1]
    mof_loading[index] = data / data[0]
    plt.scatter(x=temperature_range, y=mof_loading[index], label=NAMES[index])

    def fit_function(x, a, b):
        return a * x + b

    popt[index], pcov = scipy.optimize.curve_fit(f=fit_function, xdata=temperature_range, ydata=mof_loading[index])

alpha = - numpy.average(popt, axis=0)[0]
print(f"The thermal compression coefficient of {MATERIAL} is {alpha}.")

density_hauer = []
density_ozawa = []
density_ozawa_custom = []
for temperature in temperature_range:
    density_hauer.append(density.hauer(temperature=temperature,
                                       temperature_boiling=adsorbate_data["temperature_boiling"],
                                       density_boiling=adsorbate_data["density_boiling"],
                                       thermal_expansion_coefficient=alpha))
    density_ozawa.append(density.ozawa(temperature=temperature,
                                       temperature_boiling=adsorbate_data["temperature_boiling"],
                                       density_boiling=adsorbate_data["density_boiling"]))
    density_ozawa_custom.append(adsorbate_data["density_boiling"] 
                                * numpy.exp(-alpha * (temperature - adsorbate_data["temperature_boiling"])))

reduced_hauer = numpy.divide(density_hauer, density_hauer[0])
reduced_ozawa = numpy.divide(density_ozawa, density_ozawa[0])
reduced_ozawa_custom = numpy.divide(density_ozawa_custom, density_ozawa_custom[0])

plt.scatter(x=temperature_range, y=reduced_hauer, marker="v", label="Hauer")
plt.scatter(x=temperature_range, y=reduced_ozawa, marker="v", label="Ozawa")
plt.scatter(x=temperature_range, y=reduced_ozawa_custom, marker="v", label="Ozawa Custom")

plt.legend()
plt.show()
