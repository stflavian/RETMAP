"""Adsorbate Density Calculator

This script contains the implemented methods of calculating the temperature dependent density of the adsorbate. All
subsequent method of calculating the aforementioned value should be added in this file.

The methods supported are:
    1. Ozawa's method
    2. Empirical method
    3. Hauer's method

All methods take in a set of material dependent and environmental parameters and return a float value representing the
density of the adsorbate in kg/m3. In case of different units needed, the output value can be converted using an
external function. Changing the units in the functions is not recommended, as it may impact the functionality of the
other modules using this code.

Usage:
"""

# Standard libraries
import math
import importlib.resources

# Local libraries
from adsorpyon import constants
from adsorpyon import input_reader

# Third-party libraries
import numpy
import scipy.interpolate


def empirical(pressure_critical: float, temperature_critical: float, molecular_mass: float) -> float:
    """
    Calculates the adsorbate density based on an empirical formula, which is not temperature dependent.
    :param pressure_critical: Critical pressure of the adsorbate in MPa
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param molecular_mass: Molecular mass of the adsorbate in g/mol.
    :return: Density in kg/m3.
    """
    return 8 * pressure_critical * molecular_mass / (constants.GAS_CONSTANT * temperature_critical) * 1000


def hauer(temperature: float, temperature_boiling: float, density_boiling: float,
          thermal_expansion_coefficient: float) -> float:
    """
    Calculates the temperature dependent adsorbate density based on
    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_boiling: Boiling temperature of the adsorbate in K.
    :param thermal_expansion_coefficient: Thermal expansion coefficient in the adsorbed phase in 1/K.
    :param density_boiling: Density of the adsorbate at the boiling point in kg/m3.
    :return: Density in kg/m3.
    """
    return density_boiling * (1 - thermal_expansion_coefficient * (temperature - temperature_boiling))


def ozawa(temperature: float, temperature_boiling: float, density_boiling: float,
          thermal_expansion_coefficient: float) -> float:
    """
    Calculates the temperature dependent adsorbate density based on Ozawa's method, represented by an exponential
    formula.
    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_boiling: Boiling temperature of the adsorbate in K.
    :param density_boiling: Density of the adsorbate at the boiling point in kg/m3.
    :param thermal_expansion_coefficient: Thermal expansion coefficient in 1/K.
    :return: Density in kg/m3.
    """
    return density_boiling * math.exp(-thermal_expansion_coefficient * (temperature - temperature_boiling))


def extrapolation(temperature: float, file: str, adsorbate_name: str) -> float:
    """
    Calculates the density by extrapolating the data found in a data file.

    The target should be a two-column file, where the first column contains the temperature and the second one the
    density. For the extrapolation, a second order polynomial is used. If the input temperature is found in
    the temperature range covered by the data file, interpolation is used to determine the value.

    :param temperature: Temperature at which the experiment is conducted in K.
    :param file: Path to file containing reference data.
    :return: Density in the same units as the input file.
    """

    if file == "local":
        file = importlib.resources.files("adsorpyon").joinpath(f"library/density/{adsorbate_name}.dat")

    data = input_reader.create_data_list(file)
    data = numpy.array(data)

    def fit_function(x, a, b, c):
        return a * x**2 + b * x + c

    if temperature <= numpy.max(data[:, 0]):
        interpolation_function = scipy.interpolate.CubicSpline(data[:, 0], data[:, 1], extrapolate=True)
        return interpolation_function(temperature)
    else:
        # noinspection PyTupleAssignmentBalance
        popt, pcov = scipy.optimize.curve_fit(fit_function, data[:, 0], data[:, 1])
        return fit_function(temperature, *popt)