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

import numpy
import scipy.optimize


UNIVERSAL_GAS_CONSTANT = 8.31446  # [J/(mol*K)] or [m3*Pa/(mol*K)]


def dubinin(temperature: float, temperature_critical: float, pressure_critical: float) -> float:
    """
    Calculates the temperature dependent saturation pressure based on Dubinin's method.
    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :return: Saturation pressure in MPa.
    """
    return pressure_critical * (temperature / temperature_critical)**2


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
    return pressure_critical * (temperature / temperature_critical)**k


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

    # Create the function for the root finder
    def peng_robinson_polynomial(p_critical, t_critical, t, p_guess):
        a = 0.45724 * (UNIVERSAL_GAS_CONSTANT * t_critical) ** 2 / p_critical
        b = 0.07780 * UNIVERSAL_GAS_CONSTANT * t_critical / p_critical
        kappa = 0.37464 + 1.54226 * acentric_factor - 0.26992 * acentric_factor ** 2
        alpha = (1 + kappa * (1 - (t / t_critical) ** 0.5)) ** 2
        A = a * alpha * p_guess / (UNIVERSAL_GAS_CONSTANT * t) ** 2
        B = b * p_guess / (UNIVERSAL_GAS_CONSTANT * t)
        return [1, B - 1, A - 3 * B ** 2 - 2 * B, B ** 3 + B ** 2 - A * B]

    # Create the function to calculate the fugacity coefficient
    def fugacity_coefficient(Z, p_critical, t_critical, t, p_guess):
        a = 0.45724 * (UNIVERSAL_GAS_CONSTANT * t_critical) ** 2 / p_critical
        b = 0.07780 * UNIVERSAL_GAS_CONSTANT * t_critical / p_critical
        kappa = 0.37464 + 1.54226 * acentric_factor - 0.26992 * acentric_factor ** 2
        alpha = (1 + kappa * (1 - (t / t_critical) ** 0.5)) ** 2
        A = a * alpha * p_guess / (UNIVERSAL_GAS_CONSTANT * t) ** 2
        B = b * p_guess / (UNIVERSAL_GAS_CONSTANT * t)
        return np.exp(Z - 1 - np.log(Z - B) - A / (B * 2 * 2 ** 0.5) * np.log((Z + 2.414 * B) / (Z - 0.414 * B)))

    # Create a function for the solver to determine the saturation pressure
    def fugacity_ratio(p_guess):
        roots = np.roots(peng_robinson_polynomial(p_critical=pressure_critical,
                                                  t_critical=temperature_critical, t=temperature,
                                                  p_guess=p_guess))

        abs_roots = np.absolute(roots)

        fugacity_vapor = fugacity_coefficient(Z=np.max(abs_roots), p_critical=pressure_critical,
                                              t_critical=temperature_critical, t=temperature,
                                              p_guess=p_guess)

        fugacity_liquid = fugacity_coefficient(Z=np.min(abs_roots), p_critical=pressure_critical,
                                               t_critical=temperature_critical, t=temperature,
                                               p_guess=p_guess)

        return fugacity_vapor / fugacity_liquid - 1

    return fsolve(func=fugacity_ratio, x0=np.array(pressure_guess))[0]