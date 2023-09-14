
# Local libraries
import density
import saturation_pressure

# Third-party libraries
import matplotlib.pyplot as plt


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

    saturation_pressure_methods = {
        "dubinin": saturation_pressure_dubinin,
        "amankwah": saturation_pressure_amankwah,
        "extrapolation": saturation_pressure_extrapolation,
        "polynomial_water": saturation_pressure_polynomial_water,
        "peng_robinson": saturation_pressure_peng_robinson,
        "peng_robinson_extrapolation": saturation_pressure_preos_extrapolation,
        "prsv1_extrapolation": saturation_pressure_prsv1_extrapolation,
        "prsv2_extrapolation": saturation_pressure_prsv2_extrapolation,
        "widom_banuti": saturation_pressure_widom_banuti
    }

    if method in saturation_pressure_methods.keys():
        adsorbate_saturation_pressure = saturation_pressure_methods[method]()
    else:
        raise ValueError(f"{method} is not a valid adsorbate saturation "
                         f"pressure computation method. Change the method or check for spelling errors!")

    return adsorbate_saturation_pressure


def plot_data(index: int, data: dict, input_dictionary: dict):
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

    if input_dictionary[index]["TEMPERATURES"]:
        label = f"{input_dictionary[index]['TEMPERATURES']}K"
    elif input_dictionary[index]["PRESSURES"]:
        label = f"{input_dictionary[index]['PRESSURES']} MPa"
    else:
        label = input_dictionary[index]["DATA_FILES"]

    if input_dictionary[index]["DATA_TYPES"] is "isotherm":
        plt.scatter(data[index]["pressure"], data[index]["loading"], label=label)
        plt.xlabel("Pressure [MPa]")
        plt.ylabel("Adsorbed amount [mg/g]")
    elif input_dictionary[index]["DATA_TYPES"] is "isobar":
        plt.scatter(data[index]["temperature"], data[index]["loading"], label=label)
        plt.xlabel("Temperature [K]")
        plt.ylabel("Adsorbed amount [mg/g]")
    elif input_dictionary[index]["DATA_TYPES"] is "characteristic":
        plt.scatter(data[index]["potential"], data[index]["volume"], label=label)
        plt.xlabel("Adsorption potential [kJ/mol]")
        plt.ylabel("Adsorption volume [g/ml]")
    else:
        raise ValueError(f"{input_dictionary[index]['DATA_TYPES']} is not a recognized argument!")

    # if logarithmic.lower() == "yes":
    #     plt.xscale('log')
    # plt.legend()
    # if save.lower() == "yes":
    #     os.makedirs(name="Plots", exist_ok=True)
    #     plt.savefig(f"Plots/{data[0]['adsorbate']}_in_{data[0]['adsorbent']}_isotherms.png")
