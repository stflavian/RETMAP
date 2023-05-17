import numpy as np
from scipy.optimize import minimize
import json

adsorbate_data = json.load(open("../Properties/CO2_carbon_dioxide.prop"))

# Adsorbate properties and general constants
UNIVERSAL_GAS_CONSTANT = 8.31446  # [J/(mol*K)] or [m3*Pa/(mol*K)]
ACENTRIC_FACTOR = adsorbate_data["acentric_factor"]  # [unitless]
TEMPERATURE_CRITICAL = adsorbate_data["temperature_critical"]  # [K]
PRESSURE_CRITICAL = adsorbate_data["pressure_critical"] * 10**6  # [Pa]

# Experimental variables
TEMPERATURE = 300
PRESSURE_GUESS = 6.7*10**6


def peng_robinson_polynomial(pressure_critical, temperature_critical, temperature, pressure_guess):
    a = 0.45724 * (UNIVERSAL_GAS_CONSTANT * temperature_critical)**2 / pressure_critical
    b = 0.07780 * UNIVERSAL_GAS_CONSTANT * temperature_critical / pressure_critical
    kappa = 0.37464 + 1.54226 * ACENTRIC_FACTOR - 0.26992 * ACENTRIC_FACTOR**2
    alpha = (1 + kappa * (1 - (temperature/temperature_critical)**0.5))**2
    A = a * alpha * pressure_guess / (UNIVERSAL_GAS_CONSTANT * temperature)**2
    B = b * pressure_guess / (UNIVERSAL_GAS_CONSTANT * temperature)
    return [1, B - 1, A - 3*B**2 - 2*B, B**3 + B**2 - A*B]


def fugacity_coefficient(Z, pressure_critical, temperature_critical, temperature, pressure_guess):
    a = 0.45724 * (UNIVERSAL_GAS_CONSTANT * temperature_critical)**2 / pressure_critical
    b = 0.07780 * UNIVERSAL_GAS_CONSTANT * temperature_critical / pressure_critical
    kappa = 0.37464 + 1.54226 * ACENTRIC_FACTOR - 0.26992 * ACENTRIC_FACTOR**2
    alpha = (1 + kappa * (1 - (temperature/temperature_critical)**0.5))**2
    A = a * alpha * pressure_guess / (UNIVERSAL_GAS_CONSTANT * temperature)**2
    B = b * pressure_guess / (UNIVERSAL_GAS_CONSTANT * temperature)
    return np.exp(Z - 1 - np.log(Z - B) - A / (B * 2 * 2**0.5) * np.log((Z + 2.414 * B) / (Z - 0.414 * B)))


res = np.roots(peng_robinson_polynomial(pressure_critical=PRESSURE_CRITICAL, temperature_critical=TEMPERATURE_CRITICAL,
                                        temperature=TEMPERATURE, pressure_guess=PRESSURE_GUESS))

vapor_fc = fugacity_coefficient(Z=np.max(res), pressure_critical=PRESSURE_CRITICAL,
                                temperature_critical=TEMPERATURE_CRITICAL, temperature=TEMPERATURE,
                                pressure_guess=PRESSURE_GUESS)

liquid_fc = fugacity_coefficient(Z=np.min(res), pressure_critical=PRESSURE_CRITICAL,
                                 temperature_critical=TEMPERATURE_CRITICAL, temperature=TEMPERATURE,
                                 pressure_guess=PRESSURE_GUESS)

print(res, vapor_fc/liquid_fc)
