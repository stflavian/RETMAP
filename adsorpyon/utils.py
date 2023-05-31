"""Utility functions

This script contains various utility functions for plotting, converting and writing data to files.

Usage:
"""

import os

import numpy
import scipy.optimize
import matplotlib.pyplot as plt


def plot_isotherm(data: dict, logarithmic: str, save: str):
    """
    Plots the adsorption isotherm from the experimental data read.
    :param data: Data dictionary containing all necessary information.
    :param logarithmic: Indicates if the pressure should be plotted using a logarithmic scale.
    :param save: Indicates if the resulting plot should be saved; saves for "yes", does not save for any other string.
    It is case-insensitive.
    :return: A plot with the adsorption isotherms from the experimental data.
    """
    plt.figure()
    for index in data:
        plt.scatter(data[index]["pressure"], data[index]["adsorbed_amount"], label=data[index]["name"])
    plt.xlabel("Pressure [MPa]")
    plt.ylabel("Adsorbed amount [mg/g]")
    if logarithmic.lower() == "yes":
        plt.xscale('log')
    plt.legend()
    if save.lower() == "yes":
        os.makedirs(name="Plots", exist_ok=True)
        plt.savefig("Plots/isotherm.png")


def plot_characteristic_curve(data: dict, save: str):
    """
    Plots the characteristic curve from the experimental data read.
    :param data: Data dictionary containing all necessary information.
    :param save: Indicates if the resulting plot should be saved; saves for "yes", does not save for any other string.
    It is case-insensitive.
    :return: A plot with the characteristic curves from the experimental data.
    """
    plt.figure()
    for index in data:
        plt.scatter(data[index]["adsorption_potential"], data[index]["adsorption_volume"], label=data[index]["name"])
    plt.gca().set_xlim(left=0)
    plt.xlabel("Adsorption potential [kJ/mol]")
    plt.ylabel("Adsorption volume [ml/g]")
    plt.legend()
    if save.lower() == "yes":
        os.makedirs(name="Plots", exist_ok=True)
        plt.savefig("Plots/characteristic_curve.png")


def show_plots():
    plt.show()


def write_data(data: dict, index: int):
    os.makedirs(name="Output", exist_ok=True)
    if data[index]["data_type"] == "isotherm":
        with open(file=f"Output/char_curve_{data[index]['name']}.dat", mode="w") as file:
            file.write("#Adsorption_potential [kJ/mol]    Adsorption_volume [ml/g] \n")
            for potential, volume in zip(data[index]["adsorption_potential"], data[index]["adsorption_volume"]):
                file.write(f"{potential}       {volume} \n")
    elif data[index]["data_type"] == "characteristic curve":
        with open(file=f"Output/isotherm_{data[index]['name']}.dat", mode="w") as file:
            file.write("#Pressure [MPa]       Adsorbed_amount [mg/g] \n")
            for pressure, amount in zip(data[index]["pressure"], data[index]["adsorbed_amount"]):
                file.write(f"{pressure}       {amount} \n")


def convert_input(unit: str, adsorbate_data: dict) -> float:
    """
    Returns a conversion factor for the input units to the standard ones: MPa, mg/g, kJ/mol, ml/g.
    :param unit: The unit of the input data.
    :param adsorbate_data: The dictionary containing the adsorbate data.
    :return: A number that the input is multiplied with to be converted to the intended unit.
    """
    # Pressure
    if unit == "MPa":
        conversion_factor = 1
    elif unit == "kPa":
        conversion_factor = 0.001
    elif unit == "Pa":
        conversion_factor = 0.000001
    elif unit == "bar":
        conversion_factor = 0.1

    # Adsorbed amount
    elif unit in ["mg/g", "g/kg"]:
        conversion_factor = 1
    elif unit in ["mol/kg", "mmol/g"]:
        conversion_factor = adsorbate_data["molecular_mass"]

    # Adsorption potential
    elif unit == "kJ/mol":
        conversion_factor = 1
    elif unit == "J/mol":
        conversion_factor = 0.001

    # Adsorption volume
    elif unit in ["ml/g", "l/kg", "cm3/g", "dm3/kg"]:
        conversion_factor = 1

    # Not a recognized unit
    else:
        raise ValueError(f"Input unit {unit} is not a recognized unit!")
    return conversion_factor


def convert_output(unit: str, adsorbate_data: dict) -> float:
    """
    Returns a conversion factor for the standard units to the user defined ones.
    :param unit: The unit of the input data.
    :param adsorbate_data: The dictionary containing the adsorbate data.
    :return: A number that the input is multiplied with to be converted to the intended unit.
    """
    return 1 / convert_input(unit=unit, adsorbate_data=adsorbate_data)


def evaluate(data: dict, adsorbate_data: dict) -> float:

    def fit_function(x, a, b, c):
        return a / (1 + numpy.exp(-b * (x - c)))

    for key in data:
        scipy.optimize.curve_fit()
    return 0

