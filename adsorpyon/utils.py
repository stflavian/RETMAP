"""Utility functions

This script contains various utility functions for plotting, converting and writing data to files.

Usage:
"""

import os

import physics

import numpy
import scipy.interpolate
import scipy.stats
import matplotlib.pyplot as plt

AXES_SIZE = 19
TICK_SIZE = 13
LEGEND_SIZE = 14
FIGURE_SIZE = (7, 6)
LEFT_LIMIT = None
RIGHT_LIMIT = None


def plot_isotherm(data: dict, logarithmic: str, save: str):
    """
    Plots the adsorption isotherm from the experimental data read.
    :param data: Data dictionary containing all necessary information.
    :param logarithmic: Indicates if the pressure should be plotted using a logarithmic scale.
    :param save: Indicates if the resulting plot should be saved; saves for "yes", does not save for any other string.
    It is case-insensitive.
    :return: A plot with the adsorption isotherms from the experimental data.
    """
    plt.figure(figsize=FIGURE_SIZE)
    plt.rc('axes', labelsize=AXES_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=TICK_SIZE)  # fontsize of the x tick labels
    plt.rc('ytick', labelsize=TICK_SIZE)  # fontsize of the y tick labels
    plt.rc('legend', fontsize=LEGEND_SIZE)  # fontsize of the legend
    # plt.xlim(left=LEFT_LIMIT)
    # plt.xlim(right=RIGHT_LIMIT)
    for index in data:
        plt.scatter(data[index]["pressure"], data[index]["adsorbed_amount"], label=data[index]["name"])
    plt.xlabel("Pressure [MPa]")
    plt.ylabel("Adsorbed amount [mg/g]")
    if logarithmic.lower() == "yes":
        plt.xscale('log')
    plt.legend()
    if save.lower() == "yes":
        os.makedirs(name="Plots", exist_ok=True)
        plt.savefig(f"Plots/{data[0]['adsorbate']}_in_{data[0]['adsorbent']}_isotherms.png")


def plot_enthalpy(data: dict, save: str):
    plt.figure(figsize=FIGURE_SIZE)
    plt.rc('axes', labelsize=AXES_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=TICK_SIZE)  # fontsize of the x tick labels
    plt.rc('ytick', labelsize=TICK_SIZE)  # fontsize of the y tick labels
    plt.rc('legend', fontsize=LEGEND_SIZE)  # fontsize of the legend
    for index in data:
        plt.scatter(data[index]["adsorbed_amount"], data[index]["adsorption_enthalpy"], label=data[index]["name"])
    plt.xlabel("Adsorbed amount [mg/g]")
    plt.ylabel("Adsorption enthalpy [kJ/mol]")
    plt.legend()
    if save.lower() == "yes":
        os.makedirs(name="Plots", exist_ok=True)
        plt.savefig(f"Plots/{data[0]['adsorbate']}_in_{data[0]['adsorbent']}_enthalpy.png")


def plot_characteristic_curve(data: dict, save: str):
    """
    Plots the characteristic curve from the experimental data read.
    :param data: Data dictionary containing all necessary information.
    :param save: Indicates if the resulting plot should be saved; saves for "yes", does not save for any other string.
    It is case-insensitive.
    :return: A plot with the characteristic curves from the experimental data.
    """
    plt.figure(figsize=FIGURE_SIZE)
    plt.rc('axes', labelsize=AXES_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=TICK_SIZE)  # fontsize of the x tick labels
    plt.rc('ytick', labelsize=TICK_SIZE)  # fontsize of the y tick labels
    plt.rc('legend', fontsize=LEGEND_SIZE)  # fontsize of the legend
    for index in data:
        plt.scatter(data[index]["adsorption_potential"], data[index]["adsorption_volume"], label=data[index]["name"])
    plt.gca().set_xlim(left=0)
    plt.xlabel("Adsorption potential [kJ/mol]")
    plt.ylabel("Adsorption volume [ml/g]")
    plt.legend()
    if save.lower() == "yes":
        os.makedirs(name="Plots", exist_ok=True)
        plt.savefig(f"Plots/{data[0]['adsorbate']}_in_{data[0]['adsorbent']}_characteristic_curve.png")


def show_plots():
    plt.show()


def write_data(data: dict, index: int):
    os.makedirs(name="Output", exist_ok=True)
    if data[index]["data_type"] == "isotherm":
        with open(file=f"Output/{data[index]['adsorbate']}_in_{data[index]['adsorbent']}_char_curve_{data[index]['temperature']}.dat", mode="w") as file:
            file.write("#Adsorption_potential [kJ/mol]    Adsorption_volume [ml/g] \n")
            for potential, volume in zip(data[index]["adsorption_potential"], data[index]["adsorption_volume"]):
                file.write(f"{potential}       {volume} \n")
    elif data[index]["data_type"] == "characteristic curve":
        with open(file=f"Output/{data[index]['adsorbate']}_in_{data[index]['adsorbent']}_isotherm_{data[index]['temperature']}.dat", mode="w") as file:
            file.write("#Pressure [MPa]       Adsorbed_amount [mg/g] \n")
            for pressure, amount in zip(data[index]["pressure"], data[index]["adsorbed_amount"]):
                file.write(f"{pressure}       {amount} \n")


def convert_input(unit: str, molecular_mass: float) -> float:
    """
    Returns a conversion factor for the input units to the standard ones: MPa, mg/g, kJ/mol, ml/g.
    :param unit: The unit of the input data.
    :param molecular_mass: The molecular mass of the molecule.
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

    # Temperature
    elif unit == "K":
        conversion_factor = 1

    # Adsorbed amount
    elif unit in ["mg/g", "g/kg"]:
        conversion_factor = 1
    elif unit in ["mol/kg", "mmol/g"]:
        conversion_factor = molecular_mass

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


def convert_output(unit: str, molecular_mass: float) -> float:
    """
    Returns a conversion factor for the standard units to the user defined ones.
    :param unit: The unit of the input data.
    :param molecular_mass: The molecular mass of the molecule.
    :return: A number that the input is multiplied with to be converted to the intended unit.
    """
    return 1 / convert_input(unit=unit, molecular_mass=molecular_mass)


def evaluate_characteristic_curve(data: dict, temperature_reference_isotherm: float, save: str):
    """
    Calculate the accuracy of the characteristic curves.
    :param data:
    :param temperature_reference_isotherm:
    :param save:
    :return:
    """
    plt.figure(figsize=FIGURE_SIZE)
    plt.rc('axes', labelsize=AXES_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=TICK_SIZE)  # fontsize of the x tick labels
    plt.rc('ytick', labelsize=TICK_SIZE)  # fontsize of the y tick labels
    plt.rc('legend', fontsize=LEGEND_SIZE)  # fontsize of the legend
    for index in data:
        if temperature_reference_isotherm == data[index]["temperature"]:
            index_reference = index

    interpolation_function = scipy.interpolate.interp1d(x=data[index_reference]["adsorption_potential"],
                                                        y=data[index_reference]["adsorption_volume"],
                                                        fill_value="extrapolate")

    amount_reference = numpy.array([])
    amount_predicted = numpy.array([])
    for index in data:
        if index is not index_reference:
            amount_reference = numpy.append(amount_reference, data[index]["adsorption_volume"])
            amount_predicted = numpy.append(amount_predicted, interpolation_function(data[index]["adsorption_potential"]))

    r_score = scipy.stats.linregress(x=amount_reference, y=amount_predicted).rvalue
    rss_score = numpy.sum(numpy.square(numpy.subtract(amount_reference, amount_predicted)))

    os.makedirs(name="Output", exist_ok=True)
    with open(file=f"Output/{data[index_reference]['adsorbate']}_in_{data[index_reference]['adsorbent']}_evaluation.log", mode="w") as file:
        file.write(f"Characteristic curve evaluation metrics \n")
        file.write("\n")
        file.write(f"Adsorbate: {data[index_reference]['adsorbate']}            Adsorbent: {data[index_reference]['adsorbent']} \n")
        file.write(f"Temperature of the reference isotherm: {data[index_reference]['temperature']}K \n")
        file.write("========================================================================= \n")
        file.write(f"Correlation coefficient (r): {r_score} \n")
        file.write(f"Residual sum of squares (RSS): {rss_score} \n")
        file.write("\n")

    plt.scatter(x=amount_predicted, y=amount_reference)
    plt.plot(data[index_reference]["adsorption_volume"], data[index_reference]["adsorption_volume"],
             linestyle="dashed", color="red")

    plt.xlabel("Predicted adsorption volume [ml/g]")
    plt.ylabel(f"Reference adsorption volume [ml/g]")
    if save.lower() == "yes":
        os.makedirs(name="Plots", exist_ok=True)
        plt.savefig(f"Plots/{data[index_reference]['adsorbate']}_in_{data[index_reference]['adsorbent']}_evaluation.png")


def predict_isotherms(data: dict, temperature_reference_isotherm: float, logarithmic: str, save: str):
    """
    Calculate the accuracy of the characteristic curves.
    :param data:
    :param temperature_reference_isotherm:
    :param logarithmic:
    :param save:
    :return:
    """
    plt.figure(figsize=FIGURE_SIZE)
    plt.rc('axes', labelsize=AXES_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=TICK_SIZE)  # fontsize of the x tick labels
    plt.rc('ytick', labelsize=TICK_SIZE)  # fontsize of the y tick labels
    plt.rc('legend', fontsize=LEGEND_SIZE)  # fontsize of the legend
    #plt.xlim(left=9*10**-4)
    #plt.xlim(right=3*10**1)
    #plt.ylim((0, 150))
    for index in data:
        if temperature_reference_isotherm == data[index]["temperature"]:
            index_reference = index

    interpolation_function = scipy.interpolate.interp1d(x=data[index_reference]["adsorption_potential"],
                                                        y=data[index_reference]["adsorption_volume"],
                                                        fill_value="extrapolate")

    for index in data:
        max_potential = numpy.max(data[index]["adsorption_potential"])
        min_potential = numpy.min(data[index]["adsorption_potential"])
        adsorption_potential = numpy.linspace(start=min_potential, stop=max_potential, num=50)

        pressure = physics.get_pressure(adsorption_potential=adsorption_potential,
                                        saturation_pressure=data[index]["saturation_pressure"],
                                        temperature=data[index]["temperature"])

        amount = physics.get_adsorbed_amount(adsorption_volume=interpolation_function(adsorption_potential),
                                             adsorbate_density=data[index]["density"])

        predicted_amount = physics.get_adsorbed_amount(adsorption_volume=interpolation_function(data[index]["adsorption_potential"]),
                                                       adsorbate_density=data[index]["density"])

        plt.scatter(x=data[index]["pressure"], y=data[index]["adsorbed_amount"],
                    label=f"{data[index]['name']} Simulated")

        plt.plot(pressure, amount, label=f"{data[index]['name']} Predicted", linestyle="--")

        predicted_amount[predicted_amount < 0] = 0
        relative_error = numpy.absolute(numpy.divide(predicted_amount, data[index]["adsorbed_amount"]) - 1)
        relative_error[numpy.isnan(relative_error)] = numpy.nanmean(relative_error)
        mape_score = numpy.sum(relative_error) / len(relative_error) * 100

        absolute_error = numpy.absolute(numpy.subtract(predicted_amount, data[index]["adsorbed_amount"]))
        absolute_error[numpy.isnan(absolute_error)] = numpy.nanmean(absolute_error)
        mae_score = numpy.sum(absolute_error) / len(absolute_error)

        mae_score_scaled = mae_score / numpy.max(data[index]["adsorbed_amount"])

        rmse_score = numpy.sqrt(numpy.sum(numpy.square(absolute_error)) / len(absolute_error))

        os.makedirs(name="Output", exist_ok=True)
        with open(file=f"Output/{data[index_reference]['adsorbate']}_in_{data[index_reference]['adsorbent']}_prediction_{data[index]['temperature']}K.log", mode="w") as file:
            file.write(f"Isotherm prediction metrics \n")
            file.write("\n")
            file.write(f"Adsorbate: {data[index_reference]['adsorbate']}            Adsorbent: {data[index_reference]['adsorbent']} \n")
            file.write(f"Temperature of the reference isotherm: {data[index_reference]['temperature']}K \n")
            file.write("========================================================================= \n")
            file.write(f"Mean absolute percentage error (MAPE): {mape_score} \n")
            file.write(f"Mean absolute error (MAE): {mae_score} \n")
            file.write(f"Root mean squared error (RMSE): {rmse_score} \n")
            file.write(f"Mean absolute error scaled (MAE / maximum_loading): {mae_score_scaled} \n")
            file.write(f"Maximum recorded loading: {numpy.max(data[index]['adsorbed_amount'])} \n")
            file.write(f"Maximum absolute error loading: {numpy.max(absolute_error)} \n")
            file.write("\n")


        with open(file=f"Output/{data[index_reference]['adsorbate']}_in_{data[index_reference]['adsorbent']}_prediction_{data[index]['temperature']}K.dat", mode="w") as file:
            file.write("#Pressure [MPa]       Adsorbed_amount [mg/g] \n")
            for amount_point, pressure_point in zip(amount, pressure):
                file.write(f"{pressure_point}           {amount_point} \n")


    if logarithmic.lower() == "yes":
        plt.xscale('log')
    plt.legend()
    plt.xlabel("Pressure [MPa]")
    plt.ylabel("Adsorbed amount [mg/g]")
    if save.lower() == "yes":
        os.makedirs(name="Plots", exist_ok=True)
        plt.savefig(f"Plots/{data[index_reference]['adsorbate']}_in_{data[index_reference]['adsorbent']}_predictions.png")
