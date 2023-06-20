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
SATURATION_FILE = f"../adsorpyon/SaturationProperties/{MOLECULE}-saturation_properties.dat"
PROPERTIES_FILE = f"../adsorpyon/Properties/{MOLECULE}.prop"

reference_data = numpy.genfromtxt(SATURATION_FILE)
adsorbate_data = json.load(open(PROPERTIES_FILE))

temperatures = numpy.array(numpy.linspace(start=adsorbate_data["temperature_boiling"], stop=400, num=200))

saturation_pressure = {'dubinin': None, 'amankwah': None, 'extrapolation': None, 'preos_extrapolation': None,
                       'prsv1_extrapolation': None, 'prsv2_extrapolation': None, 'widom-banuti': None}

for method in saturation_pressure:
    saturation_pressure[method] = []
    for temperature in temperatures:
        #if method == 'PREOS':
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
        #elif method == 'prsv1_extrapolation':
        #    saturation_pressure[method].append(pressure.equation_extrapolation(temperature=temperature,
        #                                                                       temperature_critical=adsorbate_data["temperature_critical"],
        #                                                                       pressure_critical=adsorbate_data["pressure_critical"],
        #                                                                       acentric_factor=adsorbate_data["acentric_factor"],
        #                                                                       temperature_boiling=adsorbate_data["temperature_boiling"],
        #                                                                       equation="prsv1",
        #                                                                       kappa1=adsorbate_data["kappa1"],
        #                                                                       kappa2=adsorbate_data["kappa2"],
        #                                                                       kappa3=adsorbate_data["kappa3"],
        #                                                                       function="polynomial2"))
        # elif method == 'prsv2_extrapolation':
        #     saturation_pressure[method].append(pressure.equation_extrapolation(temperature=temperature,
        #                                                                        temperature_critical=adsorbate_data["temperature_critical"],
        #                                                                        pressure_critical=adsorbate_data["pressure_critical"],
        #                                                                        acentric_factor=adsorbate_data["acentric_factor"],
        #                                                                        temperature_boiling=adsorbate_data["temperature_triple_point"],
        #                                                                        equation="prsv2",
        #                                                                        kappa1=adsorbate_data["kappa1"],
        #                                                                        kappa2=adsorbate_data["kappa2"],
        #                                                                        kappa3=adsorbate_data["kappa3"],
        #                                                                        function="custom"))
        elif method == 'widom-banuti':
            saturation_pressure[method].append((pressure.widombanuti(temperature=temperature,
                                                                     temperature_critical=adsorbate_data["temperature_critical"],
                                                                     pressure_critical=adsorbate_data["pressure_critical"],
                                                                     species_parameter=5.386,
                                                                     acentric_factor=adsorbate_data["acentric_factor"])))
    print(f"Method {method} completed!")

plt.figure()
plt.axvline(x=adsorbate_data["temperature_critical"])
#plt.plot(temperatures, saturation_pressure['PREOS'], label="Peng-Robinson")
#plt.plot(temperatures, saturation_pressure['Dubinin'], label="Dubinin")
plt.plot(temperatures, saturation_pressure['amankwah'], label="Amankwah")
plt.plot(temperatures, saturation_pressure['extrapolation'], label="Extrapolation")
plt.plot(temperatures, saturation_pressure['preos_extrapolation'], label="PREOS_ext")
#plt.plot(temperatures, saturation_pressure['prsv1_extrapolation'], label="PRSV1_ext")
#plt.plot(temperatures, saturation_pressure['widom-banuti'], label="WB")
plt.plot(numpy.genfromtxt(f'{MOLECULE}_critical_density.dat')[:, 0], numpy.genfromtxt(f'{MOLECULE}_critical_density.dat')[:, 1])
plt.plot(reference_data[:, 0], reference_data[:, 1], label="Validation", linestyle='-')
plt.xlabel("Temperature [K]")
plt.ylabel("Saturation pressure [MPa]")
plt.legend()
plt.show()
