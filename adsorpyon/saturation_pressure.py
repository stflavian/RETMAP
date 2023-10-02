"""Saturation Pressure Methods

This file contains a library of methods that can be used for calculating the temperature dependent saturation pressure
of an adsorbate. Further methods of determining this value should be added in this file.

Currently supported methods are:
    1. Dubinin's method;
    2. Amankwah's method;
    3. Extrapolation of experimental data;
    4. Polynomial expression for water;
    5. Fugacity equilibration using the polynomial form of the Peng-Robinson equation of state;
    6. Fugacity equilibration using the PRSV1 equation;
    7. Fugacity equilibration using the PRSV2 equation;
    8. Extrapolation of data obtained using Peng-Robinson, PRSV1, and PRSV2 below the critical temperature;
    9. Banuti's equation for the Widom line.

All methods take as input the temperature and a set of material dependent parameters and return a float value
representing the saturation pressure of the adsorbate in megapascals (MPa). In case different pressure units are needed,
the output value can be converted from MPa using the function utils.convert_output(). Changing the units inside the
functions is not recommended, as it may break other modules that make use of this file.
"""

# Standard libraries
import warnings

# Local libraries
import physics

# Third-party libraries
import numpy
import scipy.optimize
import scipy.interpolate


def dubinin(temperature: float, temperature_critical: float, pressure_critical: float) -> float:
    """
    Calculates the saturation pressure based on Dubinin's method.

    Dubinin's method matches the reduced saturation pressure at a given temperature with the square of the
    reduced temperature. Source material: https://doi.org/10.1021/cr60204a006.

    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :return: Saturation pressure in MPa.
    """
    return pressure_critical * (temperature / temperature_critical) ** 2


def amankwah(temperature: float, temperature_critical: float, pressure_critical: float, k: float) -> float:
    """
    Calculates the saturation pressure based on Amankwah's method.

    Amankwah's method is a modified version of Dubinin's method where the exponent 2 is replaced by a unique value for
    each individual adsorbate-adsorbent pair. Source material: https://doi.org/10.1016/0008-6223(95)00079-S.

    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :param k: Exponent of the reduced temperature; k=2 results in Dubinin's method.
    :return: Saturation pressure in MPa.
    """
    return pressure_critical * (temperature / temperature_critical) ** k


def extrapolation(temperature: float, file: str) -> float:
    """
    Calculates the saturation pressure by extrapolating the data found in a data file.

    The target should be a two-column file, where the first column contains the temperature and the second one the
    saturation pressure. For the extrapolation, a second order polynomial is used. If the input temperature is found in
    the temperature range covered by the data file, interpolation is used to determine the value.

    :param temperature: Temperature at which the experiment is conducted in K.
    :param file: Path to file containing reference data.
    :return: Saturation pressure in the same units as the input file.
    """
    data = numpy.genfromtxt(file)

    def fit_function(x, a, b, c):
        return a * x**2 + b * x + c

    if temperature <= numpy.max(data[:, 0]):
        interpolation_function = scipy.interpolate.interp1d(data[:, 0], data[:, 1], fill_value="extrapolate")
        return interpolation_function(temperature)
    else:
        # noinspection PyTupleAssignmentBalance
        popt, pcov = scipy.optimize.curve_fit(fit_function, data[:, 0], data[:, 1])
        return fit_function(temperature, *popt)


def polynomial_water(temperature: float) -> float:
    """
    Calculates the saturation pressure of water based on an empirically determined seventh order polynomial.

    This function SHOULD NOT be used for any adsorbent other than water, NOR should it be used for
    temperatures above the critical temperature of water.

    :param temperature: Temperature at which the experiment is conducted in K.
    :return: Saturation pressure in MPa.
    """
    return (- 1.14798e-11 * temperature**7 + 2.23756e-8 * temperature**6 - 1.54376e-5 * temperature**5
            + 0.00443279 * temperature**4 - 0.177671 * temperature**3 - 193.14 * temperature**2
            + 42890.6 * temperature - 2.87726e+6) / 1_000_000


def pengrobinson(temperature: float, temperature_critical: float, pressure_critical: float, pressure_guess: float,
                 acentric_factor: float) -> float:
    """
    Calculates the saturation pressure using the polynomial form of Peng-Robinson's equation of state.

    Calculates the saturation pressure by equilibrating the fugacities of the vapor and liquid phases of the adsorbate
    according to the Peng-Robinson equation of state. It calculates the roots of the polynomial, which represent the
    compressibilities of the vapor and liquid phases, and uses them to determine the fugacity coefficients for the
    two phases. It then solves for the pressure at which the two coefficients are equal to each other, thus satisfying
    saturation conditions. Source material: https://doi.org/10.1021/i160057a011.

    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :param pressure_guess: Initial guess of the saturation pressure in MPa.
    :param acentric_factor: The acentric factor of the adsorbate.
    :return: Saturation pressure in MPa.
    """

    # Ignore warning regarding Deprecation since Python3.7 is used
    warnings.filterwarnings("ignore", category=numpy.VisibleDeprecationWarning)

    # Create a function for the solver to determine the saturation pressure
    def fugacity_ratio(p_guess):
        compressibility_vapor = physics.get_compressibility(pressure_critical=pressure_critical, equation='preos',
                                                            temperature_critical=temperature_critical,
                                                            temperature=temperature, pressure=p_guess,
                                                            acentric_factor=acentric_factor, state='vapor', kappa1=0,
                                                            kappa2=0, kappa3=0)

        compressibility_liquid = physics.get_compressibility(pressure_critical=pressure_critical, equation='preos',
                                                             temperature_critical=temperature_critical,
                                                             temperature=temperature, pressure=p_guess,
                                                             acentric_factor=acentric_factor, state='liquid', kappa1=0,
                                                             kappa2=0, kappa3=0)

        fugacity_vapor = physics.get_fugacity_coefficient(compressibility=compressibility_vapor, equation='preos',
                                                          pressure_critical=pressure_critical,
                                                          temperature_critical=temperature_critical,
                                                          temperature=temperature, pressure=p_guess,
                                                          acentric_factor=acentric_factor, kappa1=0, kappa2=0, kappa3=0)

        fugacity_liquid = physics.get_fugacity_coefficient(compressibility=compressibility_liquid, equation='preos',
                                                           pressure_critical=pressure_critical,
                                                           temperature_critical=temperature_critical,
                                                           temperature=temperature, pressure=p_guess,
                                                           acentric_factor=acentric_factor, kappa1=0, kappa2=0,
                                                           kappa3=0)

        return fugacity_vapor / fugacity_liquid - 1

    return scipy.optimize.fsolve(func=fugacity_ratio, x0=numpy.array(pressure_guess))[0]


def prsv1(temperature: float, temperature_critical: float, pressure_critical: float, pressure_guess: float,
          acentric_factor: float, kappa1: float) -> float:
    """
    Calculates the saturation pressure using the PRSV1 equation.

    Calculates the saturation pressure by equilibrating the fugacities of the vapor and liquid phases of the adsorbate
    according to the PRSV1 equation, a modified version of the Peng-Robinson equation of state. It calculates the roots
    of the polynomial, which represent the compressibilities of the vapor and liquid phases, and uses them to determine
    the fugacity coefficients for the two phases. It then solves for the pressure at which the two coefficients are
    equal to each other, thus satisfying saturation conditions. Source material: https://doi.org/10.1002/cjce.5450640224.

    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :param pressure_guess: Saturation pressure in MPa.
    :param acentric_factor: The acentric factor of the adsorbate.
    :param kappa1: First molecule specific constant, given in the PRSV1 paper.
    :return: Saturation pressure in MPa.
    """

    # Ignore warning regarding Deprecation since Python3.7 is used
    warnings.filterwarnings("ignore", category=numpy.VisibleDeprecationWarning)

    # Create a function for the solver to determine the saturation pressure
    def fugacity_ratio(p_guess):
        compressibility_vapor = physics.get_compressibility(pressure_critical=pressure_critical, equation='prsv1',
                                                            temperature_critical=temperature_critical,
                                                            temperature=temperature, pressure=p_guess,
                                                            acentric_factor=acentric_factor, state='vapor',
                                                            kappa1=kappa1, kappa2=0, kappa3=0)

        compressibility_liquid = physics.get_compressibility(pressure_critical=pressure_critical, equation='prsv1',
                                                             temperature_critical=temperature_critical,
                                                             temperature=temperature, pressure=p_guess,
                                                             acentric_factor=acentric_factor, state='liquid',
                                                             kappa1=kappa1, kappa2=0, kappa3=0)

        fugacity_vapor = physics.get_fugacity_coefficient(compressibility=compressibility_vapor, equation='prsv1',
                                                          pressure_critical=pressure_critical,
                                                          temperature_critical=temperature_critical,
                                                          temperature=temperature, pressure=p_guess,
                                                          acentric_factor=acentric_factor, kappa1=kappa1, kappa2=0,
                                                          kappa3=0)

        fugacity_liquid = physics.get_fugacity_coefficient(compressibility=compressibility_liquid, equation='prsv1',
                                                           pressure_critical=pressure_critical,
                                                           temperature_critical=temperature_critical,
                                                           temperature=temperature, pressure=p_guess,
                                                           acentric_factor=acentric_factor, kappa1=kappa1, kappa2=0,
                                                           kappa3=0)

        return fugacity_vapor / fugacity_liquid - 1

    return scipy.optimize.fsolve(func=fugacity_ratio, x0=numpy.array(pressure_guess))[0]


def prsv2(temperature: float, temperature_critical: float, pressure_critical: float, pressure_guess: float,
          acentric_factor: float, kappa1: float, kappa2: float, kappa3: float) -> float:
    """
    Calculates the saturation pressure using the PRSV2 equation.

    Calculates the saturation pressure by equilibrating the fugacities of the vapor and liquid phases of the adsorbate
    according to the PRSV2 equation, a modified version of the PRSV1 equation. It calculates the roots of the
    polynomial, which represent the compressibilities of the vapor and liquid phases, and uses them to determine
    the fugacity coefficients for the two phases. It then solves for the pressure at which the two coefficients are
    equal to each other, thus satisfying saturation conditions. Source material: https://doi.org/10.1002/cjce.5450640516.

    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :param pressure_guess: Saturation pressure in MPa.
    :param acentric_factor: The acentric factor of the adsorbate.
    :param kappa1: First molecule specific constant, given in the PRSV1 paper.
    :param kappa2: Second molecule specific constant, given in the PRSV2 paper.
    :param kappa3: Third molecule specific constant, given in the PRSV2 paper.
    :return: Saturation pressure in MPa.
    """

    # Ignore warning regarding Deprecation since Python3.7 is used
    warnings.filterwarnings("ignore", category=numpy.VisibleDeprecationWarning)

    # Create a function for the solver to determine the saturation pressure
    def fugacity_ratio(p_guess):
        compressibility_vapor = physics.get_compressibility(pressure_critical=pressure_critical, equation='prsv2',
                                                            temperature_critical=temperature_critical,
                                                            temperature=temperature, pressure=p_guess,
                                                            acentric_factor=acentric_factor, state='vapor',
                                                            kappa1=kappa1, kappa2=kappa2, kappa3=kappa3)

        compressibility_liquid = physics.get_compressibility(pressure_critical=pressure_critical, equation='prsv2',
                                                             temperature_critical=temperature_critical,
                                                             temperature=temperature, pressure=p_guess,
                                                             acentric_factor=acentric_factor, state='liquid',
                                                             kappa1=kappa1, kappa2=kappa2, kappa3=kappa3)

        fugacity_vapor = physics.get_fugacity_coefficient(compressibility=compressibility_vapor, equation='prsv2',
                                                          pressure_critical=pressure_critical,
                                                          temperature_critical=temperature_critical,
                                                          temperature=temperature, pressure=p_guess,
                                                          acentric_factor=acentric_factor, kappa1=kappa1, kappa2=kappa2,
                                                          kappa3=kappa3)

        fugacity_liquid = physics.get_fugacity_coefficient(compressibility=compressibility_liquid, equation='prsv2',
                                                           pressure_critical=pressure_critical,
                                                           temperature_critical=temperature_critical,
                                                           temperature=temperature, pressure=p_guess,
                                                           acentric_factor=acentric_factor, kappa1=kappa1,
                                                           kappa2=kappa2, kappa3=kappa3)

        return fugacity_vapor / fugacity_liquid - 1

    return scipy.optimize.fsolve(func=fugacity_ratio, x0=numpy.array(pressure_guess))[0]


def equation_extrapolation(temperature: float, temperature_critical: float, pressure_critical: float,
                           acentric_factor: float, temperature_boiling: float, equation: str, kappa1: float,
                           kappa2: float, kappa3: float, function: str) -> float:
    """
    Calculates the saturation pressure above the critical point by extrapolating the results of Peng-Robinson's
    equation of state, PRSV1 equation or PRSV2 equation.

    For the extrapolation, the user can choose between a second order polynomial, Amankwah's equation, and a custom
    equation. If the input temperature is found in the temperature range covered by the data file, interpolation is used
    to determine the value.

    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :param acentric_factor: The acentric factor of the adsorbate.
    :param temperature_boiling: Boiling temperature of the adsorbate in K.
    :param equation: Equation used below the critical point; preos, prsv1, or prsv2.
    :param kappa1: First molecule specific constant, given in the PRSV1 paper.
    :param kappa2: Second molecule specific constant, given in the PRSV2 paper.
    :param kappa3: Third molecule specific constant, given in the PRSV2 paper.
    :param function: Function used for the extrapolation of the results above the critical temperature.
    :return: Saturation pressure in MPa.
    """

    temp_range = numpy.linspace(start=temperature_boiling, stop=temperature_critical, num=50)
    temp_range = numpy.flipud(temp_range)

    pressure_guess = 1

    subcritical_pressures = []
    if equation == "preos":
        for index, temp in enumerate(temp_range):
            subcritical_pressures.append(pengrobinson(temperature=temp, temperature_critical=temperature_critical,
                                                      pressure_critical=pressure_critical, pressure_guess=pressure_guess,
                                                      acentric_factor=acentric_factor))
            pressure_guess = subcritical_pressures[index]
    elif equation == "prsv1":
        for index, temp in enumerate(temp_range):
            subcritical_pressures.append(prsv1(temperature=temp, temperature_critical=temperature_critical,
                                               pressure_critical=pressure_critical, pressure_guess=pressure_guess,
                                               acentric_factor=acentric_factor, kappa1=kappa1))
            pressure_guess = subcritical_pressures[index]
    elif equation == "prsv2":
        for index, temp in enumerate(temp_range):
            subcritical_pressures.append(prsv2(temperature=temp, temperature_critical=temperature_critical,
                                               pressure_critical=pressure_critical, pressure_guess=pressure_guess,
                                               acentric_factor=acentric_factor, kappa1=kappa1, kappa2=kappa2,
                                               kappa3=kappa3))
            pressure_guess = subcritical_pressures[index]
    else:
        raise ValueError(f"Equation type {equation} is not 'preos', 'prsv1' or 'prsv2'. Check the string!")

    if function == "polynomial2":
        def fit_function(x, a, b, c):
            return a * x ** 2 + b * x + c
    elif function == "amankwah":
        def fit_function(x, k):
            return pressure_critical * (x / temperature_critical)**k
    elif function == "custom":
        def fit_function(x, a, b, c):
            return a * x ** 1.1 + b * x ** 0.5 + c
    else:
        raise ValueError(f"No known function {function}!")

    if temperature <= temperature_critical:
        interpolation_function = scipy.interpolate.interp1d(temp_range, subcritical_pressures)
        return interpolation_function(temperature)
    else:
        # noinspection PyTupleAssignmentBalance
        popt, pcov = scipy.optimize.curve_fit(fit_function, temp_range, subcritical_pressures)
        return fit_function(temperature, *popt)


def widombanuti(temperature: float, temperature_critical: float, pressure_critical: float,
                species_parameter: float, acentric_factor: float) -> float:
    """
    Calculates the saturation pressure using Banuti's empirical equation.

    Source material: https://doi.org/10.1103/PhysRevE.95.052120.

    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :param species_parameter: Molecule specific constant, given in the source material.
    :param acentric_factor: The acentric factor of the adsorbate.
    :return: Saturation pressure in MPa.
    """
    if temperature >= temperature_critical:
        return numpy.exp(species_parameter*(temperature/temperature_critical - 1)) * pressure_critical
    else:
        return pengrobinson(temperature=temperature, temperature_critical=temperature_critical,
                            pressure_critical=pressure_critical, pressure_guess=0.001, acentric_factor=acentric_factor)


def critical_isochore_model(temperature: float, temperature_critical: float, pressure_critical: float) -> float:
    """
    Calculate the pressure on the critical isochore using an empirical model.
    :param temperature: Temperature at which the experiment is conducted in K.
    :param temperature_critical: Critical temperature of the adsorbate in K.
    :param pressure_critical: Critical pressure of the adsorbate in MPa.
    :return: Saturation pressure in MPa.
    """
    return temperature * 5.65 * pressure_critical / temperature_critical

