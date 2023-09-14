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

# Local libraries
import constants


def ozawa(temperature: float, temperature_boiling: float, density_boiling: float) -> float:
    """
    Calculates the temperature dependent adsorbate density based on Ozawa's method, represented by an exponential
    formula.
    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_boiling: Boiling temperature of the adsorbate in K.
    :param density_boiling: Density of the adsorbate at the boiling point in kg/m3.
    :return: Density in kg/m3.
    """
    return density_boiling * math.exp(-0.0025 * (temperature - temperature_boiling))


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


def ozawa_modified(temperature: float, temperature_boiling: float, density_boiling: float,
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
