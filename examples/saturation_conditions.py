"""Saturation Conditions Calculator

This script returns the saturation pressure for a temperature range chosen by the user.
"""

import json
import sys
sys.path.append("/home/flavian/Desktop/Bachelor's End Project/adsorpyon-0.3.1/adsorpyon/")

from adsorpyon import saturation_pressure as pressure

import numpy
import matplotlib.pyplot as plt


MATERIAL = "CO2"
START_TEMPERATURE = 253
STOP_TEMPERATURE = 393
STEP_TEMPERATURE = 10


reference_data = numpy.genfromtxt(f"../adsorpyon/SaturationProperties/{MATERIAL}-saturation_properties.dat")
adsorbate_data = json.load(open(f"../adsorpyon/Properties/{MATERIAL}.prop"))

temperatures = numpy.array(numpy.linspace(start=START_TEMPERATURE, stop=STOP_TEMPERATURE,
                                          num=int((STOP_TEMPERATURE - START_TEMPERATURE)/STEP_TEMPERATURE + 1),
                                          endpoint=True))

saturation_pressure = []
print(f"Temperature [K]             Saturation pressure [Pa]")
for index, temperature in enumerate(temperatures):
    saturation_pressure.append(pressure.preos_extrapolation(temperature=temperature,
                                                            temperature_critical=adsorbate_data["temperature_critical"],
                                                            pressure_critical=adsorbate_data["pressure_critical"],
                                                            acentric_factor=adsorbate_data["acentric_factor"],
                                                            temperature_boiling=adsorbate_data["temperature_boiling"]))
    print(f"{temperature}               {saturation_pressure[index]*1_000_000}")


plt.figure()
plt.axvline(x=adsorbate_data["temperature_critical"])
plt.plot(temperatures, saturation_pressure, label="Experimental")
plt.plot(reference_data[:, 0], reference_data[:, 1], label="Validation", linestyle='-')
plt.xlabel("Temperature [K]")
plt.ylabel("Saturation pressure [MPa]")
plt.legend()
plt.show()
