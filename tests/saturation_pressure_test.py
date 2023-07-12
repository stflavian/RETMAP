"""Saturation Pressure Test

This script tests the scale of the saturation pressure and acts as a primary comparison tool between the implemented
methods.
"""

import json
import sys

sys.path.append("/home/flavian/Desktop/Bachelor's End Project/adsorpyon-0.3.1/adsorpyon/")
from adsorpyon import saturation_pressure as pressure

import numpy
import matplotlib.pyplot as plt

MOLECULE = "CH4"
SATURATION_FILE = f"../adsorpyon/IsochoricProperties/{MOLECULE}_critical_density.dat"
PROPERTIES_FILE = f"../adsorpyon/Properties/{MOLECULE}.prop"

reference_data = numpy.genfromtxt(SATURATION_FILE)
adsorbate_data = json.load(open(PROPERTIES_FILE))

temperatures = numpy.array(numpy.linspace(start=adsorbate_data["temperature_boiling"], stop=400, num=200))

saturation_pressure = {'preos': None, 'dubinin': None, 'amankwah': None, 'extrapolation': None, 'preos_extrapolation': None,
                       'prsv1_extrapolation': None, 'prsv2_extrapolation': None, 'widom-banuti': None}

AXES_SIZE = 19
TICK_SIZE = 13
LEGEND_SIZE = 14
FIGURE_SIZE = (7, 6)

for method in saturation_pressure:
    saturation_pressure[method] = []
    for temperature in temperatures:
        # if method == 'preos':
        #    saturation_pressure[method].append(pressure.pengrobinson(temperature=temperature,
        #                                                             temperature_critical=adsorbate_data["temperature_critical"],
        #                                                             pressure_critical=adsorbate_data["pressure_critical"],
        #                                                             pressure_guess=0.0001,
        #                                                             acentric_factor=adsorbate_data["acentric_factor"]))
        if method == 'dubinin':
            saturation_pressure[method].append(pressure.dubinin(temperature=temperature,
                                                                temperature_critical=adsorbate_data["temperature_critical"],
                                                                pressure_critical=adsorbate_data["pressure_critical"]))
        elif method == 'amankwah':
            saturation_pressure[method].append(pressure.amankwah(temperature=temperature,
                                                                 temperature_critical=adsorbate_data["temperature_critical"],
                                                                 pressure_critical=adsorbate_data["pressure_critical"],
                                                                 k=adsorbate_data["amankwah_exponent"]))
        elif method == 'extrapolation':
            saturation_pressure[method].append(pressure.extrapolation(temperature=temperature,
                                                                      file=SATURATION_FILE))
        elif method == 'preos_extrapolation':
            saturation_pressure[method].append(pressure.equation_extrapolation(temperature=temperature,
                                                                               temperature_critical=adsorbate_data["temperature_critical"],
                                                                               pressure_critical=adsorbate_data["pressure_critical"],
                                                                               acentric_factor=adsorbate_data["acentric_factor"],
                                                                               temperature_boiling=adsorbate_data["temperature_boiling"],
                                                                               equation="preos",
                                                                               kappa1=adsorbate_data["kappa1"],
                                                                               kappa2=adsorbate_data["kappa2"],
                                                                               kappa3=adsorbate_data["kappa3"],
                                                                               function="polynomial2"))
        # elif method == 'prsv1_extrapolation':
        #     saturation_pressure[method].append(pressure.equation_extrapolation(temperature=temperature,
        #                                                                        temperature_critical=adsorbate_data["temperature_critical"],
        #                                                                        pressure_critical=adsorbate_data["pressure_critical"],
        #                                                                        acentric_factor=adsorbate_data["acentric_factor"],
        #                                                                        temperature_boiling=adsorbate_data["temperature_boiling"],
        #                                                                        equation="prsv1",
        #                                                                        kappa1=adsorbate_data["kappa1"],
        #                                                                        kappa2=adsorbate_data["kappa2"],
        #                                                                        kappa3=adsorbate_data["kappa3"],
        #                                                                        function="polynomial2"))
        # elif method == 'prsv2_extrapolation':
        #     saturation_pressure[method].append(pressure.equation_extrapolation(temperature=temperature,
        #                                                                        temperature_critical=adsorbate_data["temperature_critical"],
        #                                                                        pressure_critical=adsorbate_data["pressure_critical"],
        #                                                                        acentric_factor=adsorbate_data["acentric_factor"],
        #                                                                        temperature_boiling=adsorbate_data["temperature_boiling"],
        #                                                                        equation="prsv2",
        #                                                                        kappa1=adsorbate_data["kappa1"],
        #                                                                        kappa2=adsorbate_data["kappa2"],
        #                                                                        kappa3=adsorbate_data["kappa3"],
        #                                                                        function="polynomial2"))
        elif method == 'widom-banuti':
            saturation_pressure[method].append((pressure.widombanuti(temperature=temperature,
                                                                     temperature_critical=adsorbate_data["temperature_critical"],
                                                                     pressure_critical=adsorbate_data["pressure_critical"],
                                                                     species_parameter=5.386,
                                                                     acentric_factor=adsorbate_data["acentric_factor"])))
    print(f"Method {method} completed!")
    print(saturation_pressure[method])

plt.figure(figsize=FIGURE_SIZE)
plt.rc('axes', labelsize=AXES_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=TICK_SIZE)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=TICK_SIZE)  # fontsize of the y tick labels
plt.rc('legend', fontsize=LEGEND_SIZE)  # fontsize of the legen
plt.ylim(-3, 100)
plt.axvline(adsorbate_data["temperature_critical"], linestyle="--")
plt.axhline(adsorbate_data["pressure_critical"], linestyle="--")
#plt.plot(temperatures, saturation_pressure['preos'], label="Peng-Robinson")
plt.plot(temperatures, saturation_pressure['dubinin'], label="Dubinin")
plt.plot(temperatures, saturation_pressure['amankwah'], label="Amankwah")
plt.plot(temperatures, saturation_pressure['extrapolation'], label="Critical isochore")
plt.plot(temperatures, saturation_pressure['preos_extrapolation'], label="PREOS Extrapolation")
#plt.plot(temperatures, saturation_pressure['prsv1_extrapolation'], label="PRSV1 Extrapolation")
#plt.plot(temperatures, saturation_pressure['prsv2_extrapolation'], label="PRSV2 Extrapolation")
plt.plot(temperatures, saturation_pressure['widom-banuti'], label="Widom-Banuti")
plt.xlabel("Temperature [K]")
plt.ylabel("Saturation pressure [MPa]")
plt.legend()
plt.show()
