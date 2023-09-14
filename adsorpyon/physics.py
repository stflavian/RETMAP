"""Adsorption Physics

This script contains the functions used to calculate various thermodynamic properties in the project, such as pressure,
adsorption volume, and much more.

This script makes use of the local script constants and the third-party library numpy. Mathematical functions employed
are part of the numpy package to facilitate the usage of arrays.

Usage:
"""

# Local libraries
import constants

# Third-party libraries
import numpy


def get_adsorption_volume(adsorbed_amount: float, adsorbate_density: float) -> float:
    """
    Calculates the adsorption volume in terms of the adsorbed amount and adsorbate density.

    :param adsorbed_amount: Amount of adsorbate adsorbed, measured in mg/g or g/kg.
    :param adsorbate_density: Density of the adsorbate measured in kg/m3, g/l or mg/ml.
    :return: Adsorption volume in ml/g.
    """
    return adsorbed_amount / adsorbate_density


def get_adsorption_potential(temperature: float, saturation_pressure: float, pressure: float) -> float:
    """
    Calculates the adsorption potential in terms of the temperature, pressure, and saturation pressure.

    :param temperature: Temperature at which the experiment is conducted in K.
    :param saturation_pressure: Saturation pressure at given temperature in MPa.
    :param pressure: Pressure at which the experiment is conducted in MPa.
    :return: Adsorption potential in kJ/mol.
    """
    return constants.UNIVERSAL_GAS_CONSTANT * temperature * numpy.log(saturation_pressure / pressure) * 0.001


def get_adsorbed_amount(adsorption_volume: float, adsorbate_density: float) -> float:
    """
    The inverse function of get_adsorption_volume. Calculates the adsorbed amount given the adsorption volume and
    adsorbate density.

    :param adsorption_volume: Adsorption volume in ml/g.
    :param adsorbate_density: Density of the adsorbate in kg/m3, g/l or mg/ml.
    :return: Amount of adsorbate adsorbed in mg/g or g/kg.
    """
    return adsorption_volume * adsorbate_density


def get_pressure(adsorption_potential: float, saturation_pressure: float, temperature: float) -> float:
    """
    The inverse function of get_adsorption_potential. Calculates the pressure at which the experiment is conducted
    given the adsorption potential, saturation pressure of the adsorbate, and temperature of the environment.

    :param adsorption_potential: The adsorption potential in kJ/mol.
    :param saturation_pressure: The saturation pressure of the adsorbate in MPa.
    :param temperature: The temperature at which the experiment is conducted in K.
    :return: Pressure in MPa.
    """
    return saturation_pressure * numpy.exp(-adsorption_potential * 1000
                                           / (constants.UNIVERSAL_GAS_CONSTANT * temperature))


def get_temperature(adsorption_potential: float, saturation_pressure: float, pressure: float) -> float:
    """
    The inverse function of get_adsorption_potential. Calculates the temperature at which the experiment is conducted
    given the adsorption potential, saturation pressure of the adsorbate, and pressure of the environment.

    :param adsorption_potential: The adsorption potential in kJ/mol.
    :param saturation_pressure: The saturation pressure of the adsorbate in MPa.
    :param pressure: The pressure at which the experiment is conducted in MPa.
    :return: Temperature in K.
    """
    return 1000 * adsorption_potential / constants.UNIVERSAL_GAS_CONSTANT / numpy.log(saturation_pressure / pressure)


def _peng_robinson_coefficients(temperature_critical: float, pressure_critical: float, acentric_factor: float,
                                temperature: float, pressure: float) -> list:
    a = 0.45724 * (constants.UNIVERSAL_GAS_CONSTANT * temperature_critical) ** 2 / pressure_critical
    b = 0.07780 * constants.UNIVERSAL_GAS_CONSTANT * temperature_critical / pressure_critical
    kappa = 0.37464 + 1.54226 * acentric_factor - 0.26992 * acentric_factor ** 2
    alpha = (1 + kappa * (1 - (temperature / temperature_critical) ** 0.5)) ** 2
    A = a * alpha * pressure / (constants.UNIVERSAL_GAS_CONSTANT * temperature) ** 2
    B = b * pressure / (constants.UNIVERSAL_GAS_CONSTANT * temperature)
    return [a, b, A, B]


def _prsv1_coefficients(temperature_critical: float, pressure_critical: float, acentric_factor: float,
                        temperature: float, pressure: float, kappa1: float) -> list:
    a = 0.45724 * (constants.UNIVERSAL_GAS_CONSTANT * temperature_critical) ** 2 / pressure_critical
    b = 0.07780 * constants.UNIVERSAL_GAS_CONSTANT * temperature_critical / pressure_critical
    kappa0 = 0.378893 + 1.4897153 * acentric_factor - 0.17131848 * acentric_factor**2 + 0.0196554 * acentric_factor**3
    reduced_temperature = temperature/temperature_critical
    if reduced_temperature <= 0.7:
        kappa = kappa0 + kappa1 * (1 + reduced_temperature**0.5) * (0.7 - reduced_temperature)
    else:
        kappa = kappa0
    alpha = (1 + kappa * (1 - reduced_temperature**0.5)) ** 2
    A = a * alpha * pressure / (constants.UNIVERSAL_GAS_CONSTANT * temperature) ** 2
    B = b * pressure / (constants.UNIVERSAL_GAS_CONSTANT * temperature)
    return [a, b, A, B]


def _prsv2_coefficients(temperature_critical: float, pressure_critical: float, acentric_factor: float,
                        temperature: float, pressure: float, kappa1: float, kappa2: float, kappa3: float) -> list:
    a = 0.45724 * (constants.UNIVERSAL_GAS_CONSTANT * temperature_critical) ** 2 / pressure_critical
    b = 0.07780 * constants.UNIVERSAL_GAS_CONSTANT * temperature_critical / pressure_critical
    kappa0 = 0.378893 + 1.4897153 * acentric_factor - 0.17131848 * acentric_factor**2 + 0.0196554 * acentric_factor**3
    reduced_temperature = temperature / temperature_critical
    kappa = kappa0 + (kappa1 + kappa2 * (kappa3 - reduced_temperature) * (1 - reduced_temperature**0.5)) * (1 + reduced_temperature**0.5) * (0.7 - reduced_temperature)
    alpha = (1 + kappa * (1 - reduced_temperature**0.5)) ** 2
    A = a * alpha * pressure / (constants.UNIVERSAL_GAS_CONSTANT * temperature) ** 2
    B = b * pressure / (constants.UNIVERSAL_GAS_CONSTANT * temperature)
    return [a, b, A, B]


def get_fugacity_coefficient(compressibility: float, pressure_critical: float, temperature_critical: float,
                             temperature: float, pressure: float, acentric_factor: float, kappa1: float,
                             kappa2: float, kappa3: float, equation: str) -> float:

    if equation == "preos":
        coefficients = _peng_robinson_coefficients(temperature_critical=temperature_critical,
                                                   pressure_critical=pressure_critical, acentric_factor=acentric_factor,
                                                   temperature=temperature, pressure=pressure)
    elif equation == "prsv1":
        coefficients = _prsv1_coefficients(temperature_critical=temperature_critical,
                                           pressure_critical=pressure_critical, acentric_factor=acentric_factor,
                                           temperature=temperature, pressure=pressure, kappa1=kappa1)
    elif equation == "prsv2":
        coefficients = _prsv2_coefficients(temperature_critical=temperature_critical,
                                           pressure_critical=pressure_critical, acentric_factor=acentric_factor,
                                           temperature=temperature, pressure=pressure, kappa1=kappa1, kappa2=kappa2,
                                           kappa3=kappa3)
    else:
        raise ValueError(f"Equation {equation} is not a known equation! Check the string for typos!")
    A = coefficients[2]
    B = coefficients[3]
    return numpy.exp(compressibility - 1 - numpy.log(compressibility - B)
                     - A / (B * 2 * 2**0.5) * numpy.log((compressibility + 2.414 * B) / (compressibility - 0.414 * B)))


def get_compressibility(pressure_critical: float, temperature_critical: float, temperature: float, pressure: float,
                        acentric_factor: float, kappa1: float, kappa2: float, kappa3: float, equation: str,
                        state: str) -> float:
    def peng_robinson_polynomial(p):
        if equation == "preos":
            coefficients = _peng_robinson_coefficients(temperature_critical=temperature_critical,
                                                       pressure_critical=pressure_critical,
                                                       acentric_factor=acentric_factor,
                                                       temperature=temperature, pressure=p)
        elif equation == "prsv1":
            coefficients = _prsv1_coefficients(temperature_critical=temperature_critical,
                                               pressure_critical=pressure_critical, acentric_factor=acentric_factor,
                                               temperature=temperature, pressure=pressure, kappa1=kappa1)
        elif equation == "prsv2":
            coefficients = _prsv2_coefficients(temperature_critical=temperature_critical,
                                               pressure_critical=pressure_critical, acentric_factor=acentric_factor,
                                               temperature=temperature, pressure=pressure, kappa1=kappa1, kappa2=kappa2,
                                               kappa3=kappa3)
        else:
            raise ValueError(f"Method {equation} is not supported! Check the string for typos!")
        A = coefficients[2]
        B = coefficients[3]
        return [1, B - 1, A - 3 * B ** 2 - 2 * B, B ** 3 + B ** 2 - A * B]

    compressibility_roots = numpy.absolute(numpy.roots(peng_robinson_polynomial(p=pressure)))

    if state == 'liquid':
        return numpy.min(compressibility_roots)
    elif state == 'vapor':
        return numpy.max(compressibility_roots)
    else:
        raise ValueError(f"Selected state {state} is not valid! Supported states are liquid and vapor!")


def get_adsorption_enthalpy(enthalpy_vaporization: float, adsorption_potential: numpy.ndarray, temperature: float,
                            adsorption_volume: numpy.ndarray, thermal_expansion_coefficient: float) -> float:
     entropy = thermal_expansion_coefficient * adsorption_volume * numpy.gradient(adsorption_potential, adsorption_volume)
     return enthalpy_vaporization + adsorption_potential - temperature * entropy


