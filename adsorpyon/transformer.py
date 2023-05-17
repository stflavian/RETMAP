import numpy as np
import matplotlib.pyplot as plt
import os


UNIVERSAL_GAS_CONSTANT = 8.3144  # [J/(mol*K)]


def get_adsorption_volume(adsorbed_amount: float, adsorbate_density: float) -> float:
    """
    Calculates the adsorption volume in terms of the adsorbed amount and adsorbate density.
    :param adsorbed_amount: Amount of adsorbate adsorbed, measured in mg/g or g/kg.
    :param adsorbate_density: Density of the adsorbate measured in kg/m3, g/l or mg/ml.
    :return: Adsorption volume in ml/g.
    """
    return adsorbed_amount / adsorbate_density


def get_adsorption_potential(temperature: float, saturation_pressure: float,
                             pressure: float) -> float:
    """
    Calculates the adsorption potential in terms of the temperature, pressure, and saturation pressure.
    :param temperature: Temperature at which the experiment is conducted in K.
    :param saturation_pressure: Saturation pressure at given temperature in MPa.
    :param pressure: Pressure at which the experiment is conducted in MPa.
    :return: Adsorption potential in kJ/mol.
    """
    return UNIVERSAL_GAS_CONSTANT * temperature * np.log(saturation_pressure / pressure) * 0.001


def get_adsorbed_amount(adsorption_volume: float, adsorbate_density: float) -> float:
    """
    The inverse function of get_adsorption_volume. Calculates the adsorbed amount given the adsorption volume and
    adsorbate density.
    :param adsorption_volume: Adsorption volume in ml/g.
    :param adsorbate_density: Density of the adsorbate in kg/m3, g/l or mg/ml.
    :return: Amount of adsorbate adsorbed in mg/g or g/kg.
    """
    return adsorption_volume * adsorbate_density


def get_pressure(adsorption_potential: float, saturation_pressure: float,
                 temperature: float) -> float:
    """
    The inverse function of get_adsorption_potential. Calculates the pressure at which the experiment is conducted
    given the adsorption potential, saturation pressure of the adsorbate, and temperature of the environment.
    :param adsorption_potential: The adsorption potential in kJ/mol.
    :param saturation_pressure: The saturation pressure of the adsorbate in MPa.
    :param temperature: The temperature at which the experiment is conducted in K.
    :return: Pressure in MPa.
    """
    return saturation_pressure * np.exp(-adsorption_potential * 1000 / (UNIVERSAL_GAS_CONSTANT * temperature))


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


def convert(unit: str, adsorbate_data: dict) -> float:
    """
    Returns a conversion factor based to convert the input units to the standard ones: MPa, mg/g, kJ/mol, ml/g
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
    elif unit == "mol/kg":
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
        conversion_factor = None
    return conversion_factor
