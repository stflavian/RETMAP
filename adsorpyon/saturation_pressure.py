"""Adsorbate Saturation Pressure Calculator

This script contains the implemented methods of calculating the temperature dependent saturation pressure of the
adsorbate. All subsequent method of calculating the aforementioned value should be added in this file.

The methods supported are:
    1. Dubinin's method
    2. Amankwah's method
    3. Extrapolation of experimental data
    4. Polynomial expression for water
    5. Fugacity equilibration using the Peng-Robinson equation of state

All methods take in a set of material dependent and environmental parameters and return a float value representing the
saturation pressure of the adsorbate in MPa. In case of different units needed, the output value can be converted using
an external function. Changing the units in the functions is not recommended, as it may impact the functionality of the
other modules using this code.

Usage:
"""

import warnings

from . import physics

import numpy
import scipy.optimize


def dubinin(temperature: float, temperature_critical: float, pressure_critical: float) -> float:
    """
    Calculates the temperature dependent saturation pressure based on Dubinin's method.
    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :return: Saturation pressure in MPa.
    """
    return pressure_critical * (temperature / temperature_critical) ** 2


def amankwah(temperature: float, temperature_critical: float, pressure_critical: float, k: float) -> float:
    """
    Calculates the temperature dependent saturation pressure based on Amankwah's method, a modified version of Dubinin's
    method where the exponent is unique for each individual adsorbate-adsorbent pair.
    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :param k: Exponent of the reduced temperature. k=2 results in Dubinin's method.
    :return: Saturation pressure in MPa.
    """
    return pressure_critical * (temperature / temperature_critical) ** k


def extrapolation(temperature: float, file: str) -> float:
    """
    Calculates the temperature dependent saturation pressure by extrapolating experimental data. The file should contain
    two columns, where the first one is the temperature and the second is the saturation pressure. For the extrapolation
    a second degree polynomial is used.
    :param temperature: Temperature at which the experiment is conducted in K.
    :param file: Path to file containing the experimental data.
    :return: Saturation pressure in MPa.
    """
    data = numpy.genfromtxt(file)

    def fit_function(x, a, b, c):
        return a * x**2 + b * x + c

    # noinspection PyTupleAssignmentBalance
    popt, pcov = scipy.optimize.curve_fit(fit_function, data[:, 0], data[:, 1])
    return fit_function(temperature, *popt)


def polynomial_water(temperature: float) -> float:
    """
    Calculates the temperature dependent saturation pressure based on empirical formula represented by a seventh degree
    polynomial.
    :param temperature: Temperature at which the experiment is conducted in K.
    :return: Saturation pressure in MPa.
    """
    return (- 1.14798e-11*temperature**7 + 2.23756e-8*temperature**6 - 1.54376e-5*temperature**5
            + 0.00443279*temperature**4 - 0.177671*temperature**3 - 193.14*temperature**2
            + 42890.6*temperature - 2.87726e+6) / 1000 / 1000


def pengrobinson(temperature: float, temperature_critical: float, pressure_critical: float, pressure_guess: float,
                 acentric_factor: float) -> float:
    """
    Calculates the temperature dependent saturation pressure by equilibrating the fugacities of the vapor and liquid
    phases according to the Peng-Robinson equation of state. It calculates the roots of the compressibility polynomial
    form of the EoS and uses them to determine the fugacity coefficients in liquid and vapor phase. It then solves for
    the pressure at which the two coefficients are equal to each other, thus satisfying saturation conditions.
    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :param pressure_guess: Saturation pressure in MPa.
    :param acentric_factor: The acentric factor of the adsorbate.
    :return: Saturation pressure in MPa.
    """

    # Ignore warning regarding Deprecation since Python3.7 is used
    warnings.filterwarnings("ignore", category=numpy.VisibleDeprecationWarning)

    # Create a function for the solver to determine the saturation pressure
    def fugacity_ratio(p_guess):
        compressibility_vapor = physics.get_compressibility(pressure_critical=pressure_critical, method='pr',
                                                            temperature_critical=temperature_critical,
                                                            temperature=temperature, pressure=p_guess,
                                                            acentric_factor=acentric_factor, state='vapor')

        compressibility_liquid = physics.get_compressibility(pressure_critical=pressure_critical, method='pr',
                                                             temperature_critical=temperature_critical,
                                                             temperature=temperature, pressure=p_guess,
                                                             acentric_factor=acentric_factor, state='liquid')

        fugacity_vapor = physics.get_fugacity_coefficient(compressibility=compressibility_vapor, method='pr',
                                                          pressure_critical=pressure_critical,
                                                          temperature_critical=temperature_critical,
                                                          temperature=temperature, pressure=p_guess,
                                                          acentric_factor=acentric_factor)

        fugacity_liquid = physics.get_fugacity_coefficient(compressibility=compressibility_liquid, method='pr',
                                                           pressure_critical=pressure_critical,
                                                           temperature_critical=temperature_critical,
                                                           temperature=temperature, pressure=p_guess,
                                                           acentric_factor=acentric_factor)

        return fugacity_vapor / fugacity_liquid - 1

    return scipy.optimize.fsolve(func=fugacity_ratio, x0=numpy.array(pressure_guess))[0]


def prsv1(temperature: float, temperature_critical: float, pressure_critical: float, pressure_guess: float,
          acentric_factor: float) -> float:
    """
    Calculates the temperature dependent saturation pressure by equilibrating the fugacities of the vapor and liquid
    phases according to the Peng-Robinson equation of state. It calculates the roots of the compressibility polynomial
    form of the EoS and uses them to determine the fugacity coefficients in liquid and vapor phase. It then solves for
    the pressure at which the two coefficients are equal to each other, thus satisfying saturation conditions.
    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :param pressure_guess: Saturation pressure in MPa.
    :param acentric_factor: The acentric factor of the adsorbate.
    :return: Saturation pressure in MPa.
    """

    # Ignore warning regarding Deprecation since Python3.7 is used
    warnings.filterwarnings("ignore", category=numpy.VisibleDeprecationWarning)

    # Create a function for the solver to determine the saturation pressure
    def fugacity_ratio(p_guess):
        compressibility_vapor = physics.get_compressibility(pressure_critical=pressure_critical, method='prsv1',
                                                            temperature_critical=temperature_critical,
                                                            temperature=temperature, pressure=p_guess,
                                                            acentric_factor=acentric_factor, state='vapor')

        compressibility_liquid = physics.get_compressibility(pressure_critical=pressure_critical, method='prsv1',
                                                             temperature_critical=temperature_critical,
                                                             temperature=temperature, pressure=p_guess,
                                                             acentric_factor=acentric_factor, state='liquid')

        fugacity_vapor = physics.get_fugacity_coefficient(compressibility=compressibility_vapor, method='prsv1',
                                                          pressure_critical=pressure_critical,
                                                          temperature_critical=temperature_critical,
                                                          temperature=temperature, pressure=p_guess,
                                                          acentric_factor=acentric_factor)

        fugacity_liquid = physics.get_fugacity_coefficient(compressibility=compressibility_liquid, method='prsv1',
                                                           pressure_critical=pressure_critical,
                                                           temperature_critical=temperature_critical,
                                                           temperature=temperature, pressure=p_guess,
                                                           acentric_factor=acentric_factor)

        return fugacity_vapor / fugacity_liquid - 1

    return scipy.optimize.fsolve(func=fugacity_ratio, x0=numpy.array(pressure_guess))[0]


def extrapolation_experimental(temperature: float, temperature_critical: float, pressure_critical: float,
                               acentric_factor: float) -> float:

    temp_range = numpy.linspace(start=50, stop=temperature_critical, num=200)
    temp_range = numpy.flipud(temp_range)

    pressure_guess = 1

    subcritical_pressures = []
    for index, temp in enumerate(temp_range):
        subcritical_pressures.append(pengrobinson(temperature=temp, temperature_critical=temperature_critical,
                                                  pressure_critical=pressure_critical, pressure_guess=pressure_guess,
                                                  acentric_factor=acentric_factor))
        pressure_guess = subcritical_pressures[index]

    def fit_function(x, a, b, c, d):
        return a * x ** 3 + b * x ** 2 + c * x + d

    # noinspection PyTupleAssignmentBalance
    popt, pcov = scipy.optimize.curve_fit(fit_function, temp_range, subcritical_pressures)
    return fit_function(temperature, *popt)

