"""Saturation Conditions Calculator

This script returns the saturation pressure for a temperature range chosen by the user.
"""

import json

import numpy

from adsorpyon import saturation_pressure as pressure

import numpy as np
import matplotlib.pyplot as plt


MATERIAL = "CH4"
START_TEMPERATURE = 253
STOP_TEMPERATURE = 393
STEP_TEMPERATURE = 10


reference_data = numpy.genfromtxt(f"../SaturationProperties/{MATERIAL}-saturation_properties.dat")
adsorbate_data = json.load(open(f"../adsorpyon/Properties/{MATERIAL}.prop"))

temperatures = np.array(np.linspace(start=START_TEMPERATURE, stop=STOP_TEMPERATURE,
                                    num=int((STOP_TEMPERATURE - START_TEMPERATURE)/STEP_TEMPERATURE + 1)))

saturation_pressure = []
print(f"Temperature [K]             Saturation pressure [Pa]")
for index, temperature in enumerate(temperatures):
    saturation_pressure.append(pressure.extrapolation_experimental(temperature=temperature,
                                                                   temperature_critical=adsorbate_data["temperature_critical"],
                                                                   pressure_critical=adsorbate_data["pressure_critical"],
                                                                   pressure_guess=0.0001,
                                                                   acentric_factor=adsorbate_data["acentric_factor"]))
    print(f"{temperature}               {saturation_pressure[index]*1_000_000}")


plt.figure()
plt.axvline(x=adsorbate_data["temperature_critical"])
plt.plot(temperatures, saturation_pressure, label="Experimental")
plt.plot(reference_data[:, 0], reference_data[:, 1], label="Validation", linestyle='-')
plt.xlabel("Temperature [K]")
plt.ylabel("Saturation pressure [MPa]")
plt.legend()
plt.show()
