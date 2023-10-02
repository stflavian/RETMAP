
# Standard libraries
import os

# Local libraries
import density
import saturation_pressure
import physics

# Third-party libraries
import matplotlib.pyplot as plt
import scipy.interpolate
import scipy.optimize
import numpy


AXES_SIZE = 19
TICK_SIZE = 13
LEGEND_SIZE = 14
FIGURE_SIZE = (7, 6)


def compute_density_from_method(method: str, temperature: float, properties_dictionary: dict) -> float:
    """
    Compute the adsorbate density using the method specified in the input file and the respective molecular properties
    and environmental conditions.

    :param method: Name of the method used to compute the adsorbate density.
    :param temperature: Temperature at which the adsorbate density is computed in K.
    :param properties_dictionary: Dictionary containing the properties of the molecule used.
    :return: Adsorbate density in kg/m3.
    """

    def density_empirical() -> float:
        return density.empirical(pressure_critical=properties_dictionary["PRESSURE_CRITICAL"],
                                 temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
                                 molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def density_hauer() -> float:
        return density.hauer(temperature=temperature,
                             temperature_boiling=properties_dictionary["TEMPERATURE_BOILING"],
                             density_boiling=properties_dictionary["DENSITY_BOILING"],
                             thermal_expansion_coefficient=properties_dictionary["THERMAL_EXPANSION_COEFFICIENT"])

    def density_ozawa() -> float:
        return density.ozawa_modified(temperature=temperature,
                                      temperature_boiling=properties_dictionary["TEMPERATURE_BOILING"],
                                      density_boiling=properties_dictionary["DENSITY_BOILING"],
                                      thermal_expansion_coefficient=properties_dictionary["THERMAL_EXPANSION_COEFFICIENT"])

    density_methods = {
        "empirical": density_empirical,
        "hauer": density_hauer,
        "ozawa": density_ozawa
    }

    if method in density_methods.keys():
        adsorbate_density = density_methods[method]()
    else:
        raise ValueError(f"{method} is not a valid adsorbate density computation method."
                         f" Change the method or check for spelling errors!")

    return adsorbate_density


def compute_saturation_pressure_from_method(method: str, temperature: float, properties_dictionary: dict,
                                            saturation_pressure_file: str) -> float:
    """
    Compute the adsorbate saturation pressure using the method specified in the input file and the respective molecular
    properties and environmental conditions.

    :param method: Name of the method used to compute the adsorbate saturation pressure.
    :param temperature: Temperature at which the adsorbate saturation pressure is computed in K.
    :param properties_dictionary: Dictionary containing the properties of the molecule used.
    :param saturation_pressure_file: Path to the file containing saturation pressure data.
    :return: Adsorbate saturation pressure in MPa.
    """

    def saturation_pressure_dubinin() -> float:
        return saturation_pressure.dubinin(temperature=temperature,
                                           temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
                                           pressure_critical=properties_dictionary["PRESSURE_CRITICAL"])

    def saturation_pressure_amankwah() -> float:
        return saturation_pressure.amankwah(temperature=temperature,
                                            temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
                                            pressure_critical=properties_dictionary["PRESSURE_CRITICAL"],
                                            k=properties_dictionary["AMANKWAH_EXPONENT"])

    def saturation_pressure_extrapolation() -> float:
        return saturation_pressure.extrapolation(temperature=temperature,
                                                 file=saturation_pressure_file)

    def saturation_pressure_polynomial_water() -> float:
        return saturation_pressure.polynomial_water(temperature=temperature)

    def saturation_pressure_peng_robinson() -> float:
        return saturation_pressure.pengrobinson(temperature=temperature,
                                                temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
                                                pressure_critical=properties_dictionary["PRESSURE_CRITICAL"],
                                                pressure_guess=1,
                                                acentric_factor=properties_dictionary["ACENTRIC_FACTOR"])

    def saturation_pressure_preos_extrapolation() -> float:
        return saturation_pressure.equation_extrapolation(temperature=temperature,
                                                          temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
                                                          pressure_critical=properties_dictionary["PRESSURE_CRITICAL"],
                                                          acentric_factor=properties_dictionary["ACENTRIC_FACTOR"],
                                                          temperature_boiling=properties_dictionary["TEMPERATURE_BOILING"],
                                                          equation="preos",
                                                          kappa1=properties_dictionary["PRSV_KAPPA1"],
                                                          kappa2=properties_dictionary["PRSV_KAPPA2"],
                                                          kappa3=properties_dictionary["PRSV_KAPPA3"],
                                                          function="polynomial2")

    def saturation_pressure_prsv1_extrapolation() -> float:
        return saturation_pressure.equation_extrapolation(temperature=temperature,
                                                          temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
                                                          pressure_critical=properties_dictionary["PRESSURE_CRITICAL"],
                                                          acentric_factor=properties_dictionary["ACENTRIC_FACTOR"],
                                                          temperature_boiling=properties_dictionary["TEMPERATURE_BOILING"],
                                                          equation="prsv1",
                                                          kappa1=properties_dictionary["PRSV_KAPPA1"],
                                                          kappa2=properties_dictionary["PRSV_KAPPA2"],
                                                          kappa3=properties_dictionary["PRSV_KAPPA3"],
                                                          function="polynomial2")

    def saturation_pressure_prsv2_extrapolation() -> float:
        return saturation_pressure.equation_extrapolation(temperature=temperature,
                                                          temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
                                                          pressure_critical=properties_dictionary["PRESSURE_CRITICAL"],
                                                          acentric_factor=properties_dictionary["ACENTRIC_FACTOR"],
                                                          temperature_boiling=properties_dictionary["TEMPERATURE_BOILING"],
                                                          equation="prsv2",
                                                          kappa1=properties_dictionary["PRSV_KAPPA1"],
                                                          kappa2=properties_dictionary["PRSV_KAPPA2"],
                                                          kappa3=properties_dictionary["PRSV_KAPPA3"],
                                                          function="polynomial2")

    def saturation_pressure_widom_banuti() -> float:
        return saturation_pressure.widombanuti(temperature=temperature,
                                               temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
                                               pressure_critical=properties_dictionary["PRESSURE_CRITICAL"],
                                               species_parameter=5.589,
                                               acentric_factor=properties_dictionary["ACENTRIC_FACTOR"])

    def satuartion_pressure_isochore() -> float:
        return saturation_pressure.critical_isochore_model(temperature=temperature,
                                                           temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
                                                           pressure_critical=properties_dictionary["PRESSURE_CRITICAL"])

    saturation_pressure_methods = {
        "dubinin": saturation_pressure_dubinin,
        "amankwah": saturation_pressure_amankwah,
        "extrapolation": saturation_pressure_extrapolation,
        "polynomial_water": saturation_pressure_polynomial_water,
        "peng_robinson": saturation_pressure_peng_robinson,
        "peng_robinson_extrapolation": saturation_pressure_preos_extrapolation,
        "prsv1_extrapolation": saturation_pressure_prsv1_extrapolation,
        "prsv2_extrapolation": saturation_pressure_prsv2_extrapolation,
        "widom_banuti": saturation_pressure_widom_banuti,
        "critical_isochore": satuartion_pressure_isochore
    }

    if method in saturation_pressure_methods.keys():
        adsorbate_saturation_pressure = saturation_pressure_methods[method]()
    else:
        raise ValueError(f"{method} is not a valid adsorbate saturation "
                         f"pressure computation method. Change the method or check for spelling errors!")

    return adsorbate_saturation_pressure


def plot_data(source_dictionary: dict, input_dictionary: dict, plot_format: str) -> None:
    """

    :param source_dictionary:
    :param input_dictionary:
    :param plot_format:
    :return:
    """

    def plot_isotherm(index):
        label = f"{source_dictionary[index]['temperature']}K"
        plt.scatter(source_dictionary[index]["pressure"], source_dictionary[index]["loading"], label=label)
        plt.xlabel("Pressure [MPa]")
        plt.ylabel("Adsorbed amount [mg/g]")

    def plot_isobar(index):
        label = f"{source_dictionary[index]['pressure']}MPa"
        plt.scatter(source_dictionary[index]["temperature"], source_dictionary[index]["loading"], label=label)
        plt.xlabel("Temperature [K]")
        plt.ylabel("Adsorbed amount [mg/g]")

    def plot_isosurface(index):
        plt.scatter(source_dictionary[index]["pressure"], source_dictionary[index]["loading"])
        plt.xlabel("Pressure [MPa]")
        plt.ylabel("Adsorbed amount [mg/g]")

    def plot_characteristic(index):
        label = f"{input_dictionary[index]['TEMPERATURES']}K"
        plt.scatter(source_dictionary[index]["potential"], source_dictionary[index]["volume"], label=label)
        plt.xlabel("Adsorption potential [kJ/mol]")
        plt.ylabel("Adsorption volume [ml/g]")

    plot_formats = {
        "isotherm": plot_isotherm,
        "isobar": plot_isobar,
        "characteristic": plot_characteristic
    }

    plt.figure(figsize=FIGURE_SIZE)
    plt.rc('axes', labelsize=AXES_SIZE)
    plt.rc('xtick', labelsize=TICK_SIZE)
    plt.rc('ytick', labelsize=TICK_SIZE)
    plt.rc('legend', fontsize=LEGEND_SIZE)

    for index in source_dictionary:
        if plot_format in plot_formats:

            plot_formats[plot_format](index)

            if plot_format == "isotherm" and input_dictionary[0]["LOGARITHMIC_PLOT"] == "yes":
                plt.xscale('log')

        elif plot_format == "isosurface":
            plt.figure(figsize=FIGURE_SIZE).add_subplot(projection='3d')

            plot_isosurface(index)

    plt.legend()



def write_data(source_dictionary: dict, input_dictionary: dict, write_format: str) -> None:
    """

    :param source_dictionary:
    :param input_dictionary:
    :param write_format:
    :return:
    """

    def write_isotherm(index, base_name) -> None:
        file_name = f"{base_name}_isotherm_{source_dictionary[index]['temperature']}K.dat"
        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write("#Pressure [MPa]         Loading [mg/g] \n")
            for pressure, loading in zip(source_dictionary[index]["pressure"], source_dictionary[index]["loading"]):
                file.write(f"{pressure}         {loading} \n")

    def write_isobar(index, base_name) -> None:
        file_name = f"{base_name}_isobar_{source_dictionary[index]['pressure']}MPa.dat"
        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write("#Temperature [K]         Loading [mg/g] \n")
            for temperature, loading in zip(source_dictionary[index]["temperature"], source_dictionary[index]["loading"]):
                file.write(f"{temperature}         {loading} \n")

    def write_isosurface(index, base_name) -> None:
        file_name = f"{base_name}_isotherm_{input_dictionary[index]['TEMPERATURES']}K.dat"
        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write("#Pressure [MPa]         Loading [mg/g] \n")
            for pressure, loading in zip(source_dictionary[index]["loading"], source_dictionary[index]["pressure"]):
                file.write(f"{pressure}         {loading} \n")

    def write_characteristic(index, base_name) -> None:
        file_name = f"{base_name}_characteristic_{input_dictionary[index]['TEMPERATURES']}K.dat"
        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write("#Potential [kJ/mol]         Volume [ml/g] \n")
            for potential, volume in zip(source_dictionary[index]["potential"], source_dictionary[index]["volume"]):
                file.write(f"{potential}         {volume} \n")

    file_write_formats = {
        "isotherm": write_isotherm,
        "isobar": write_isobar,
        "characteristic": write_characteristic
    }

    os.makedirs(name="Output", exist_ok=True)

    for index in source_dictionary:
        base_name = f"{input_dictionary[0]['ADSORBATE']}_in_{input_dictionary[0]['ADSORBENT']}"
        if write_format in file_write_formats:
            file_write_formats[write_format](index, base_name)
        elif write_format == "isosurface":
            write_isosurface(index, base_name)


def predict_data(data_dictionary: dict, input_dictionary: dict, prediction_type: str, properties_dictionary: dict) -> dict:

    def _get_pressure_boundaries(temperature: float, potential: numpy.ndarray) -> list:
        sat_pres = compute_saturation_pressure_from_method(
                                            method=input_dictionary[0]["ADSORBATE_SATURATION_PRESSURE"],
                                            temperature=temperature,
                                            properties_dictionary=properties_dictionary,
                                            saturation_pressure_file=input_dictionary[0]["SATURATION_PRESSURE_FILE"])

        minimum_pressure = physics.get_pressure(adsorption_potential=numpy.max(potential),
                                                saturation_pressure=sat_pres,
                                                temperature=temperature)

        maximum_pressure = physics.get_pressure(adsorption_potential=numpy.min(potential),
                                                saturation_pressure=sat_pres,
                                                temperature=temperature)

        return [minimum_pressure, maximum_pressure]

    def _get_temperature_boundaries(pressure: float, potential: numpy.ndarray) -> list:
        def minimum_temperature(temperature_guess: float) -> float:
            sat_pres = compute_saturation_pressure_from_method(
                method=input_dictionary[0]["ADSORBATE_SATURATION_PRESSURE"],
                temperature=temperature_guess,
                properties_dictionary=properties_dictionary,
                saturation_pressure_file=input_dictionary[0]["SATURATION_PRESSURE_FILE"])

            potential_computed = physics.get_adsorption_potential(temperature=temperature_guess,
                                                                  saturation_pressure=sat_pres,
                                                                  pressure=pressure)
            return numpy.min(potential) - potential_computed

        def maximum_temperature(temperature_guess: float) -> float:
            sat_pres = compute_saturation_pressure_from_method(
                method=input_dictionary[0]["ADSORBATE_SATURATION_PRESSURE"],
                temperature=temperature_guess,
                properties_dictionary=properties_dictionary,
                saturation_pressure_file=input_dictionary[0]["SATURATION_PRESSURE_FILE"])

            potential_computed = physics.get_adsorption_potential(temperature=temperature_guess,
                                                                  saturation_pressure=sat_pres,
                                                                  pressure=pressure)
            return numpy.max(potential) - potential_computed

        return [scipy.optimize.fsolve(minimum_temperature, x0=273)[0], scipy.optimize.fsolve(maximum_temperature, x0=273)[0]]

    def predict_isotherm():

        prediction_dictionary = {}
        for index, temperature in enumerate(input_dictionary[0]["PREDICTION_TEMPERATURES"]):
            prediction_dictionary[index] = {}
            prediction_dictionary[index]["temperature"] = temperature
            prediction_dictionary[index]["saturation_pressure"] = compute_saturation_pressure_from_method(
                method=input_dictionary[0]["ADSORBATE_SATURATION_PRESSURE"],
                temperature=temperature,
                properties_dictionary=properties_dictionary,
                saturation_pressure_file=input_dictionary[0]["SATURATION_PRESSURE_FILE"])

            prediction_dictionary[index]["density"] = compute_density_from_method(
                method=input_dictionary[0]["ADSORBATE_DENSITY"],
                temperature=temperature,
                properties_dictionary=properties_dictionary)

            pressure_boundaries = _get_pressure_boundaries(
                temperature=temperature,
                potential=data_dictionary[0]["potential"])

            prediction_dictionary[index]["pressure"] = numpy.linspace(
                start=pressure_boundaries[0],
                stop=pressure_boundaries[1])

            potential_range = physics.get_adsorption_potential(
                temperature=prediction_dictionary[index]["temperature"],
                saturation_pressure=prediction_dictionary[index]["saturation_pressure"],
                pressure=prediction_dictionary[index]["pressure"])

            prediction_dictionary[index]["loading"] = physics.get_adsorbed_amount(
                adsorption_volume=interpolation_function(potential_range),
                adsorbate_density=prediction_dictionary[index]["density"])

        return prediction_dictionary

    def predict_isobar():

        prediction_dictionary = {}
        for index, pressure in enumerate(input_dictionary[0]["PREDICTION_PRESSURES"]):
            prediction_dictionary[index] = {}
            prediction_dictionary[index]["pressure"] = pressure
            temperature_boundaries = _get_temperature_boundaries(
                pressure=pressure,
                potential=data_dictionary[0]["potential"])

            prediction_dictionary[index]["temperature"] = numpy.linspace(
                start=temperature_boundaries[0],
                stop=temperature_boundaries[1])

            prediction_dictionary[index]["saturation_pressure"] = []
            prediction_dictionary[index]["density"] = []

            for temperature in prediction_dictionary[index]["temperature"]:
                prediction_dictionary[index]["saturation_pressure"].append(compute_saturation_pressure_from_method(
                    method=input_dictionary[0]["ADSORBATE_SATURATION_PRESSURE"],
                    temperature=temperature,
                    properties_dictionary=properties_dictionary,
                    saturation_pressure_file=input_dictionary[0]["SATURATION_PRESSURE_FILE"]))

                prediction_dictionary[index]["density"].append(compute_density_from_method(
                    method=input_dictionary[0]["ADSORBATE_DENSITY"],
                    temperature=temperature,
                    properties_dictionary=properties_dictionary))

            prediction_dictionary[index]["saturation_pressure"] = numpy.array(prediction_dictionary[index]["saturation_pressure"])
            prediction_dictionary[index]["density"] = numpy.array(prediction_dictionary[index]["density"])

            potential_range = physics.get_adsorption_potential(
                temperature=prediction_dictionary[index]["temperature"],
                saturation_pressure=prediction_dictionary[index]["saturation_pressure"],
                pressure=prediction_dictionary[index]["pressure"])

            prediction_dictionary[index]["loading"] = physics.get_adsorbed_amount(
                adsorption_volume=interpolation_function(potential_range),
                adsorbate_density=prediction_dictionary[index]["density"])

        return prediction_dictionary

    prediction_formats = {
        "isotherm": predict_isotherm,
        "isobar": predict_isobar,
    }

    interpolation_function = scipy.interpolate.interp1d(x=data_dictionary[0]["potential"],
                                                        y=data_dictionary[0]["volume"],
                                                        fill_value="extrapolate")

    if prediction_type in prediction_formats:
        prediction_dictionary = prediction_formats[prediction_type]()
    else:
        raise ValueError(f"{prediction_type} is not a valid prediction type. Change the method or check "
                         f"for spelling errors!")

    return prediction_dictionary


def show_plots() -> None:
    """
    Show all created plots simultaneously in separate windows.
    """
    plt.show()
