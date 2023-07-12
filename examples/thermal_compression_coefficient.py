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

if MATERIAL == "CH4":
    MATERIAL_NAME = "methane"
else:
    MATERIAL_NAME = MATERIAL

MOF_FILES = [f"saturation_stuff/sim-{MATERIAL_NAME}-Co-MOF-74-cCO2-saturation_T.load",
             f"saturation_stuff/sim-{MATERIAL_NAME}-IRMOF-1-saturation_T.load",
             f"saturation_stuff/sim-{MATERIAL_NAME}-MIL-47-saturation_T.load",
             f"saturation_stuff/sim-{MATERIAL_NAME}-MOF-1-saturation_T.load",
             f"saturation_stuff/sim-{MATERIAL_NAME}-ZJU-198-LP-saturation_T.load"]

# MOF_FILES = [f"saturation_stuff/sim-CO2-Co-MOF-74-cCO2-saturation_T.load",
#              f"saturation_stuff/sim-CO2-IRMOF-1-saturation_T.load",
#              f"saturation_stuff/sim-CO2-MIL-47-saturation_T.load",
#              f"saturation_stuff/sim-CO2-MOF-1-saturation_T.load",
#              f"saturation_stuff/sim-CO2-ZJU-198-LP-saturation_T.load",
#              f"saturation_stuff/sim-methane-Co-MOF-74-cCO2-saturation_T.load",
#              f"saturation_stuff/sim-methane-IRMOF-1-saturation_T.load",
#              f"saturation_stuff/sim-methane-MIL-47-saturation_T.load",
#              f"saturation_stuff/sim-methane-MOF-1-saturation_T.load",
#              f"saturation_stuff/sim-methane-ZJU-198-LP-saturation_T.load",
#              f"saturation_stuff/sim-N2-Co-MOF-74-cCO2-saturation_T.load",
#              f"saturation_stuff/sim-N2-IRMOF-1-saturation_T.load",
#              f"saturation_stuff/sim-N2-MIL-47-saturation_T.load",
#              f"saturation_stuff/sim-N2-MOF-1-saturation_T.load",
#              f"saturation_stuff/sim-N2-ZJU-198-LP-saturation_T.load"]

NAMES = ["Co-MOF-74", "IRMOF-1", "MIL-47", "MOF-1", "ZJU-198"]

AXES_SIZE = 19
TICK_SIZE = 13
LEGEND_SIZE = 14
FIGURE_SIZE = (7, 6)

plt.figure(figsize=FIGURE_SIZE)
plt.rc('axes', labelsize=AXES_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=TICK_SIZE)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=TICK_SIZE)  # fontsize of the y tick labels
plt.rc('legend', fontsize=LEGEND_SIZE)  # fontsize of the legend

with open(f"../adsorpyon/Properties/{MATERIAL}.prop") as adsorbate_file:
    adsorbate_data = json.load(adsorbate_file)

temperature_range = numpy.linspace(start=253, stop=393, num=15, endpoint=True)

mof_loading = [None for _ in MOF_FILES]
popt = [None for _ in MOF_FILES]

for index, file in enumerate(MOF_FILES):
    data = numpy.genfromtxt(file)[:, 1]
    mof_loading[index] = data / data[0]
    plt.scatter(x=temperature_range, y=mof_loading[index], s=14, alpha=0.5)

    def fit_function(x, a, b):
        return a * x + b

    popt[index], pcov = scipy.optimize.curve_fit(f=fit_function, xdata=temperature_range, ydata=mof_loading[index])
    print(f"Slope of {file} is {-popt[index][0]}")

alpha = - numpy.average(popt, axis=0)[0]
print(f"The thermal compression coefficient of all materials is {alpha}.")

density_hauer_var = []
density_hauer_const = []
density_ozawa_const = []
density_ozawa_var = []
for temperature in temperature_range:
        density_ozawa_const.append(density.ozawa(temperature=temperature,
                                                 temperature_boiling=adsorbate_data["temperature_boiling"],
                                                 density_boiling=adsorbate_data["density_boiling"]))
        density_hauer_var.append(density.hauer(temperature=temperature,
                                               temperature_boiling=adsorbate_data["temperature_boiling"],
                                               density_boiling=adsorbate_data["density_boiling"],
                                               thermal_expansion_coefficient=adsorbate_data["thermal_expansion_coefficient"]))
        density_ozawa_var.append(density.ozawa_modified(temperature=temperature,
                                                        temperature_boiling=adsorbate_data["temperature_boiling"],
                                                        density_boiling=adsorbate_data["density_boiling"],
                                                        thermal_expansion_coefficient=adsorbate_data["thermal_expansion_coefficient"]))
        density_hauer_const.append(density.hauer(temperature=temperature,
                                                 temperature_boiling=adsorbate_data["temperature_boiling"],
                                                 density_boiling=adsorbate_data["density_boiling"],
                                                 thermal_expansion_coefficient=0.00165))

reduced_hauer_var = numpy.divide(density_hauer_var, density_hauer_var[0])
reduced_ozawa_const = numpy.divide(density_ozawa_const, density_ozawa_const[0])
reduced_ozawa_var = numpy.divide(density_ozawa_var, density_ozawa_var[0])
reduced_hauer_const = numpy.divide(density_hauer_const, density_hauer_const[0])


isobar = numpy.genfromtxt(f"saturation_stuff/{MATERIAL}_isobar_100MPa.dat")
reduced_isobar = numpy.divide(isobar[:, 2], isobar[0, 2])
plt.plot(temperature_range, reduced_isobar, label=f"{MATERIAL} 100 MPa", linewidth=3)

# isobar_CO2 = numpy.genfromtxt("saturation_stuff/CO2_isobar_100MPa.dat")
# isobar_CH4 = numpy.genfromtxt("saturation_stuff/CH4_isobar_100MPa.dat")
# isobar_N2 = numpy.genfromtxt("saturation_stuff/N2_isobar_100MPa.dat")

# reduced_isobar_CO2 = numpy.divide(isobar_CO2[:, 2], isobar_CO2[0, 2])
# reduced_isobar_CH4 = numpy.divide(isobar_CH4[:, 2], isobar_CH4[0, 2])
# reduced_isobar_N2 = numpy.divide(isobar_N2[:, 2], isobar_N2[0, 2])

# plt.plot(temperature_range, reduced_isobar_CO2, label="CO2 100 MPa", linewidth=3)
# plt.plot(temperature_range, reduced_isobar_CH4, label="CH4 100 MPa", linewidth=3)
# plt.plot(temperature_range, reduced_isobar_N2, label="N2 100 MPa", linewidth=3)

#saturation_data = numpy.genfromtxt("../adsorpyon/SaturationProperties/CO2-saturation_properties_short.dat")
#reduced_theory = numpy.divide(saturation_data[:, 2], 1032.4)
#temperature_range_theory = numpy.linspace(start=253, stop=300, num=95, endpoint=True)

plt.plot(temperature_range, reduced_ozawa_const, label="Ozawa (const. $\\alpha$)", marker="s", markersize=7)
plt.plot(temperature_range, reduced_hauer_var, label="Hauer (var. $\\alpha$)", marker="v", markersize=7)
plt.plot(temperature_range, reduced_ozawa_var, label="Ozawa (var. $\\alpha$)", marker="s", markersize=7)
plt.plot(temperature_range, reduced_hauer_const, label="Hauer (const. $\\alpha$)", marker="v", markersize=7)

#plt.scatter(x=temperature_range_theory, y=reduced_theory, marker="v", label="CO2 bulk")
plt.xlim((245, 400))
plt.ylim((0.6, 1))

plt.xlabel("Temperature [K]")
plt.ylabel("Relative density [-]")

plt.legend()
plt.show()
