
import json

from adsorpyon import physics

import numpy


adsorbate_data = json.load(open("../adsorpyon/Properties/CO2.prop"))

temperatures = numpy.array(numpy.linspace(start=200, stop=350, num=300))

compressibilities = []

for temperature in temperatures:
    compressibilities.append(physics.get_compressibility(pressure_critical=adsorbate_data["pressure_critical"],
                                                         temperature_critical=adsorbate_data["temperature_critical"],
                                                         temperature=temperature, pressure=))
