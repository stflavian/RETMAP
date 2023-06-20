"""Saturation Pressure Optimizer

This script tests the scale of the saturation pressure and acts as a primary comparison tool between the implemented
methods.

"""

import json
import sys

sys.path.append("/home/flavian/Desktop/Bachelor's End Project/adsorpyon-0.3.1/adsorpyon/")
from adsorpyon import saturation_pressure as pressure
from adsorpyon import density
from adsorpyon import utils
from adsorpyon import physics

import numpy
import scipy.interpolate
import scipy.stats
import scipy.optimize
import matplotlib.pyplot as plt


MOLECULE = "CH4"
TEMPERATURES = [273, 300, 333, 373, 185]

PROPERTIES_FILE = f"../adsorpyon/Properties/{MOLECULE}.prop"
SATURATION_PROPERTIES_FILE = f"../adsorpyon/SaturationProperties/{MOLECULE}-saturation_properties.dat"

if MOLECULE == "CH4":
    MOLECULE_NAME = "methane"
else:
    MOLECULE_NAME = MOLECULE

ISOTHERMS = [f"sim-{MOLECULE_NAME}-Co-MOF-74-cCO2-{TEMPERATURES[0]}K.load",
             f"sim-{MOLECULE_NAME}-Co-MOF-74-cCO2-{TEMPERATURES[1]}K.load",
             f"sim-{MOLECULE_NAME}-Co-MOF-74-cCO2-{TEMPERATURES[2]}K.load",
             f"sim-{MOLECULE_NAME}-Co-MOF-74-cCO2-{TEMPERATURES[3]}K.load",
             f"sim-{MOLECULE_NAME}-Co-MOF-74-cCO2-{TEMPERATURES[4]}K.load"]


adsorbate_data = json.load(open(PROPERTIES_FILE))
saturation_properties = numpy.genfromtxt(SATURATION_PROPERTIES_FILE)

data = {}

for index, temperature in enumerate(TEMPERATURES):
    data[index] = {}
    data[index]["temperature"] = TEMPERATURES[index]
    # data[index]["saturation_pressure"] = pressure.equation_extrapolation(temperature=data[index]["temperature"],
    #                                                                      temperature_critical=adsorbate_data["temperature_critical"],
    #                                                                      pressure_critical=adsorbate_data["pressure_critical"],
    #                                                                      acentric_factor=adsorbate_data["acentric_factor"],
    #                                                                      temperature_boiling=adsorbate_data["temperature_triple_point"],
    #                                                                      equation="preos",
    #                                                                      kappa1=adsorbate_data["kappa1"],
    #                                                                      kappa2=adsorbate_data["kappa2"],
    #                                                                      kappa3=adsorbate_data["kappa3"],
    #                                                                      function="polynomial2")

    data[index]["saturation_pressure"] = pressure.amankwah(temperature=data[index]["temperature"],
                                                           temperature_critical=adsorbate_data["temperature_critical"],
                                                           pressure_critical=adsorbate_data["pressure_critical"],
                                                           k=adsorbate_data["amankwah_exponent"])
    data[index]["density"] = density.hauer(temperature=data[index]["temperature"],
                                           temperature_boiling=adsorbate_data["temperature_boiling"],
                                           density_boiling=adsorbate_data["density_boiling"],
                                           thermal_expansion_coefficient=adsorbate_data["thermal_expansion_coefficient"])

    isotherm_data = numpy.genfromtxt(ISOTHERMS[index])
    data[index]["pressure"] = isotherm_data[:, 0] * utils.convert_input(unit="kPa", adsorbate_data=adsorbate_data)
    data[index]["adsorbed_amount"] = isotherm_data[:, 1] * utils.convert_input(unit="mol/kg", adsorbate_data=adsorbate_data)
    data[index]["adsorption_potential"] = physics.get_adsorption_potential(temperature=data[index]["temperature"],
                                                                           saturation_pressure=data[index]["saturation_pressure"],
                                                                           pressure=data[index]["pressure"])
    data[index]["adsorption_volume"] = physics.get_adsorption_volume(adsorbed_amount=data[index]["adsorbed_amount"],
                                                                     adsorbate_density=data[index]["density"])


interpolation_function = scipy.interpolate.interp1d(x=data[0]["adsorption_potential"],
                                                    y=data[0]["adsorption_volume"],
                                                    fill_value="extrapolate")

array = []
for index in data:
    data[index]["saturation_pressure_optimized"] = data[index]["saturation_pressure"]
    if index >= 1:
        def characteristic_curve_overlap(offset):
            new_pressure = data[index]["saturation_pressure"] - offset
            new_potential = physics.get_adsorption_potential(temperature=data[index]["temperature"],
                                                             saturation_pressure=new_pressure,
                                                             pressure=data[index]["pressure"])
            error = 1 - scipy.stats.linregress(x=data[index]["adsorption_volume"], y=interpolation_function(new_potential)).rvalue**2
            print(index, offset, error)
            return error

        optimized_offset = scipy.optimize.fsolve(func=characteristic_curve_overlap, x0=0)
        data[index]["saturation_pressure_optimized"] = data[index]["saturation_pressure"] - optimized_offset
        print(optimized_offset, data[index]["saturation_pressure"], data[index]["saturation_pressure_optimized"])

        data[index]["adsorption_potential_optimized"] = physics.get_adsorption_potential(temperature=data[index]["temperature"],
                                                                                         saturation_pressure=data[index]["saturation_pressure_optimized"],
                                                                                         pressure=data[index]["pressure"])
    array.append(data[index]["saturation_pressure_optimized"])


plt.figure()
for index in data:
    if index < 1:
        plt.scatter(data[index]["adsorption_potential"], data[index]["adsorption_volume"], label=f'{data[index]["temperature"]}K')
    else:
        plt.scatter(data[index]["adsorption_potential_optimized"], data[index]["adsorption_volume"], label=f'{data[index]["temperature"]}K')
plt.gca().set_xlim(left=0)
plt.xlabel("Adsorption potential [kJ/mol]")
plt.ylabel("Adsorption volume [ml/g]")
plt.legend()
plt.show()

plt.figure()
for index in data:
    plt.scatter(data[index]["adsorption_potential"], data[index]["adsorption_volume"], label=f'{data[index]["temperature"]}K')
plt.gca().set_xlim(left=0)
plt.xlabel("Adsorption potential [kJ/mol]")
plt.ylabel("Adsorption volume [ml/g]")
plt.legend()
plt.show()

plt.figure()
plt.axvline(x=adsorbate_data["temperature_critical"])
plt.scatter(TEMPERATURES, array)
plt.plot(saturation_properties[:, 0], saturation_properties[:, 1])
plt.show()
