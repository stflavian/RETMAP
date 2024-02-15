
# Standard libraries
import os

# Local libraries
import density
import saturation_pressure
import physics
import input_reader

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
        return density.empirical(
            pressure_critical=properties_dictionary["PRESSURE_CRITICAL"],
            temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def density_hauer() -> float:
        return density.hauer(
            temperature=temperature,
            temperature_boiling=properties_dictionary["TEMPERATURE_BOILING"],
            density_boiling=properties_dictionary["DENSITY_BOILING"],
            thermal_expansion_coefficient=properties_dictionary["THERMAL_EXPANSION_COEFFICIENT"])

    def density_ozawa() -> float:
        return density.ozawa_modified(
            temperature=temperature,
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
        return saturation_pressure.dubinin(
            temperature=temperature,
            temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
            pressure_critical=properties_dictionary["PRESSURE_CRITICAL"])

    def saturation_pressure_amankwah() -> float:
        return saturation_pressure.amankwah(
            temperature=temperature,
            temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
            pressure_critical=properties_dictionary["PRESSURE_CRITICAL"],
            k=properties_dictionary["AMANKWAH_EXPONENT"])

    def saturation_pressure_extrapolation() -> float:
        return saturation_pressure.extrapolation(temperature=temperature, file=saturation_pressure_file)

    def saturation_pressure_polynomial_water() -> float:
        return saturation_pressure.polynomial_water(temperature=temperature)

    def saturation_pressure_peng_robinson() -> float:
        return saturation_pressure.pengrobinson(
            temperature=temperature,
            temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
            pressure_critical=properties_dictionary["PRESSURE_CRITICAL"],
            pressure_guess=1,
            acentric_factor=properties_dictionary["ACENTRIC_FACTOR"])

    def saturation_pressure_preos_extrapolation() -> float:
        return saturation_pressure.equation_extrapolation(
            temperature=temperature,
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
        return saturation_pressure.equation_extrapolation(
            temperature=temperature,
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
        return saturation_pressure.equation_extrapolation(
            temperature=temperature,
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
        return saturation_pressure.widombanuti(
            temperature=temperature,
            temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
            pressure_critical=properties_dictionary["PRESSURE_CRITICAL"],
            species_parameter=5.589,
            acentric_factor=properties_dictionary["ACENTRIC_FACTOR"])

    def satuartion_pressure_isochore() -> float:
        return saturation_pressure.critical_isochore_model(
            temperature=temperature,
            temperature_critical=properties_dictionary["TEMPERATURE_CRITICAL"],
            pressure_critical=properties_dictionary["PRESSURE_CRITICAL"],
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
        "widom_banuti": saturation_pressure_widom_banuti,
        "critical_isochore": satuartion_pressure_isochore
    }

    if method in saturation_pressure_methods.keys():
        adsorbate_saturation_pressure = saturation_pressure_methods[method]()
    else:
        raise ValueError(f"{method} is not a valid adsorbate saturation "
                         f"pressure computation method. Change the method or check for spelling errors!")

    return adsorbate_saturation_pressure


def read_data(source_dictionary: dict, properties_dictionary: dict, input_dictionary: dict) -> None:

    def read_isotherm(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]

        pressure = numpy.array([row[0] for row in file_data])
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        loading = numpy.array([row[1] for row in file_data])
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_isobar(index: int, file_data: list) -> None:
        source_dictionary[index]["pressure"] = input_dictionary[index]["PRESSURES"]

        temperature = numpy.array([row[0] for row in file_data])
        source_dictionary[index]["temperature"] = temperature * convert_input(
            unit=source_dictionary[index]["TEMPERATURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        loading = numpy.array([row[1] for row in file_data])
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_isostere(index: int, file_data: list) -> None:
        source_dictionary[index]["loading"] = input_dictionary[index]["LOADINGS"]

        temperature = numpy.array([row[0] for row in file_data])
        source_dictionary[index]["temperature"] = temperature * convert_input(
            unit=source_dictionary[index]["TEMPERATURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        pressure = numpy.array([row[1] for row in file_data])
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_characteristic(index: int, file_data: list) -> None:

        potential = numpy.array([row[0] for row in file_data])
        source_dictionary[index]["potential"] = potential * convert_input(
            unit=input_dictionary[index]["POTENTIAL_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        volume = numpy.array([row[1] for row in file_data])
        source_dictionary[index]["volume"] = volume * convert_input(
            unit=input_dictionary[index]["VOLUME_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_langmuir(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        bp = param2 * pressure
        loading = param1 * bp / (1 + bp)
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])


    def read_n_langmuir(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        params1 = file_data[1][0]
        params2 = file_data[2][0]
        if len(params1) != len(params2):
            raise ValueError("The number of parameters for the n-Langmuir isotherm does not match!")

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        loading = numpy.zeros(shape=num)
        for param1, param2 in zip(params1, params2):
            bp = param2 * pressure
            loading += param1 * bp / (1 + bp)

        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_bet(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]
        param4 = file_data[4][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        pr = pressure / param1
        bpr = param3 * pr
        loading = param2 * bpr / ((1 - param4 * pr) * (1 - param4 + bpr))
        source_dictionary[index]["loading"] = loading* convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_anti_langmuir(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        loading = (param1 * pressure) / (1 - param2 * pressure)
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_henry(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        loading = param1 * pressure
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_freundlich(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        loading = param1 * numpy.power(pressure, 1/param2)
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_sips(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        bpiv = numpy.power(param2 * pressure, 1/param3)
        loading = param1 * bpiv / (1 + bpiv)
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_n_sips(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        params1 = file_data[1]
        params2 = file_data[2]
        params3 = file_data[3]
        if len(params1) != len(params2) or len(params1) != len(params3):
            raise ValueError("The number of parameters for the n-Langmuir isotherm does not match!")

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        loading = numpy.zeros(shape=num)
        for param1, param2, param3 in zip(params1, params2, params3):
            bpiv = numpy.power(param2 * pressure, 1/param3)
            loading += param1 * bpiv / (1 + bpiv)

        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_langmuir_freundlich(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        bpv = param2 * numpy.power(pressure, param3)
        loading = param1 * bpv / (1 + bpv)
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_n_langmuir_freundlich(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        params1 = file_data[1]
        params2 = file_data[2]
        params3 = file_data[3]
        if len(params1) != len(params2) or len(params1) != len(params3):
            raise ValueError("The number of parameters for the n-Langmuir isotherm does not match!")

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        loading = numpy.zeros(shape=num)
        for param1, param2, param3 in zip(params1, params2, params3):
            bpv = param2 * numpy.power(pressure, param3)
            loading += param1 * bpv / (1 + bpv)

        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])


    def read_redlich_peterson(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        loading = param1 * pressure / (1 + param2 * numpy.power(pressure, param3))
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])


    def read_toth(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        bp = param2 * pressure
        loading = param1 * bp / numpy.power((1 + numpy.power(bp, param3)), 1/param3)
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_unilan(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        bp = param2 * pressure
        loading = (param1 / (2 * param3) * numpy.log((1 + bp * numpy.exp(param3)) / (1 + bp * numpy.exp(-param3))))
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_obrien_myers(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]
        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        bp = param2 * pressure
        loading = param1 * (bp / (1 + bp) + param3**2 * bp * (1 - bp) / (2 * numpy.power(1 + bp, 3)))
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_quadratic(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]
        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        bp = param2 * pressure
        cp2 = param3 * numpy.power(pressure, 2)
        loading = param1 * (bp + 2 * cp2 / (1 + bp + cp2))
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_asymptotic_temkin(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]
        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        bp = param2 * pressure
        frac = bp / (1 + bp)
        loading = param1 * frac + param1 * param3 * numpy.power(frac, 2) * (frac - 1)
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    def read_bingel_walton(index: int, file_data: list) -> None:
        source_dictionary[index]["temperature"] = input_dictionary[index]["TEMPERATURES"]
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]
        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]["pressure"] = pressure * convert_input(
            unit=input_dictionary[index]["PRESSURE_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        term = -(param2 + param3) * pressure
        loading = param1 * (1 - numpy.exp(term) / (1 + param3/param2 * numpy.exp(term)))
        source_dictionary[index]["loading"] = loading * convert_input(
            unit=input_dictionary[index]["LOADING_UNITS"],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

    data_types = {
        "isotherm": read_isotherm,
        "isobar": read_isobar,
        "isostere": read_isostere,
        "characteristic": read_characteristic,
        "langmuir": read_langmuir,
        "n-langmuir": read_n_langmuir,
        "bet": read_bet,
        "anti-langmuir": read_anti_langmuir,
        "henry": read_henry,
        "freundlich": read_freundlich,
        "sips": read_sips,
        "n-sips": read_n_sips,
        "langmuir-freundlich": read_langmuir_freundlich,
        "n-langmuir-freundlich": read_n_langmuir_freundlich,
        "redlich-peterson": read_redlich_peterson,
        "toth": read_toth,
        "unilan": read_unilan,
        "obrien-myers": read_obrien_myers,
        "quadratic": read_quadratic,
        "asymptotic-temkin": read_asymptotic_temkin,
        "bingel-walton": read_bingel_walton
    }

    for index in input_dictionary:
        file_data = input_reader.create_data_list(path=input_dictionary[index]["DATA_FILES"])
        if input_dictionary[index]["DATA_TYPES"] in data_types:
            source_dictionary[index] = {}
            data_types[input_dictionary[index]["DATA_TYPES"]](index, file_data)


def compute_characteristic(source_dictionary: dict, input_dictionary: dict, properties_dictionary: dict):

    def from_isotherm(index):
        source_dictionary[index]["saturation_pressure"] = compute_saturation_pressure_from_method(
            method=input_dictionary[index]["ADSORBATE_SATURATION_PRESSURE"],
            temperature=source_dictionary[index]["temperature"],
            properties_dictionary=properties_dictionary,
            saturation_pressure_file=input_dictionary[index]["SATURATION_PRESSURE_FILE"])

        source_dictionary[index]["density"] = compute_density_from_method(
            method=input_dictionary[index]["ADSORBATE_DENSITY"],
            temperature=source_dictionary[index]["temperature"],
            properties_dictionary=properties_dictionary)

        source_dictionary[index]["potential"] = physics.get_adsorption_potential(
            temperature=source_dictionary[index]["temperature"],
            saturation_pressure=source_dictionary[index]["saturation_pressure"],
            pressure=source_dictionary[index]["pressure"])

        source_dictionary[index]["volume"] = physics.get_adsorption_volume(
            adsorbed_amount=source_dictionary[index]["loading"],
            adsorbate_density=source_dictionary[index]["density"])

    def from_isobar(index):
        saturation_pressure_array = []
        density_array = []
        for temperature in source_dictionary[index]["temperature"]:
            saturation_pressure_array.append(compute_saturation_pressure_from_method(
                method=input_dictionary[index]["ADSORBATE_SATURATION_PRESSURE"],
                temperature=temperature,
                properties_dictionary=properties_dictionary,
                saturation_pressure_file=input_dictionary[index]["SATURATION_PRESSURE_FILE"]))

            density_array.append(compute_density_from_method(
                method=input_dictionary[index]["ADSORBATE_DENSITY"],
                temperature=temperature,
                properties_dictionary=properties_dictionary))

        source_dictionary[index]["saturation_pressure"] = numpy.array(saturation_pressure_array)
        source_dictionary[index]["density"] = numpy.array(density_array)

        source_dictionary[index]["potential"] = physics.get_adsorption_potential(
            temperature=source_dictionary[index]["temperature"],
            saturation_pressure=source_dictionary[index]["saturation_pressure"],
            pressure=source_dictionary[index]["pressure"])

        source_dictionary[index]["volume"] = physics.get_adsorption_volume(
            adsorbed_amount=source_dictionary[index]["loading"],
            adsorbate_density=source_dictionary[index]["density"])

    def from_characteristic(index):
        pass

    data_types = {
        "isotherm": from_isotherm,
        "isobar": from_isobar,
        "isostere": from_characteristic,  # CHANGE LATER
        "characteristic": from_characteristic,
        "langmuir": from_isotherm,
        "n-langmuir": from_isotherm,
        "bet": from_isotherm,
        "anti-langmuir": from_isotherm,
        "henry": from_isotherm,
        "freundlich": from_isotherm,
        "sips": from_isotherm,
        "n-sips": from_isotherm,
        "langmuir-freundlich": from_isotherm,
        "n-langmuir-freundlich": from_isotherm,
        "redlich-peterson": from_isotherm,
        "toth": from_isotherm,
        "unilan": from_isotherm,
        "obrien-myers": from_isotherm,
        "quadratic": from_isotherm,
        "asymptotic-temkin": from_isotherm,
        "bingel-walton": from_isotherm
    }

    for index in source_dictionary:
        if input_dictionary[index]["DATA_TYPES"] in data_types:
            data_types[input_dictionary[index]["DATA_TYPES"]](index)


def plot_data(source_dictionary: dict, input_dictionary: dict, plot_format: str, save: str) -> None:
    """
    Create plot based on the input data type. Supports between isotherm, isobar, and characteristic curve.

    :param source_dictionary: Dictionary containing the source data for the plotting.
    :param input_dictionary: Dictionary containing the arguments found in the input file.
    :param plot_format: Format of the plot, given by the source data format. Can be isobar, isotherm or characteristic
    curve.
    :param save: Dictates if the plot is saved. Saved if "yes", otherwise do not save.
    """

    def plot_isotherm(index):
        label = f"{source_dictionary[index]['temperature']}K"
        plt.scatter(source_dictionary[index]["pressure"], source_dictionary[index]["loading"], label=label)
        plt.xlabel("Pressure [MPa]")
        plt.ylabel("Adsorbed amount [mg/g]")

    def plot_isobar(index):
        label = f"{source_dictionary[index]['pressure']} MPa"
        plt.scatter(source_dictionary[index]["temperature"], source_dictionary[index]["loading"], label=label)
        plt.xlabel("Temperature [K]")
        plt.ylabel("Adsorbed amount [mg/g]")

    def plot_isostere(index):
        label = f"{source_dictionary[index]['loading']} mg/g"
        plt.scatter(source_dictionary[index]["temperature"], source_dictionary[index]["pressure"], label=label)
        plt.xlabel("Temperature [K]")
        plt.ylabel("Pressure [MPa]")

    def plot_enthalpy(index):
        label = source_dictionary[index]["source"]
        plt.scatter(source_dictionary[index]["loading"], source_dictionary[index]["enthalpy"], label=label)
        plt.xlabel("Loading [mg/g]")
        plt.ylabel("Enthalpy of adsorption [kJ/mol]")

    def plot_characteristic(index):
        if type(source_dictionary[index]["temperature"]) is not numpy.ndarray:
            label = f"{source_dictionary[index]['temperature']}K"
        else:
            label = f"{source_dictionary[index]['pressure']} MPa"

        plt.scatter(source_dictionary[index]["potential"], source_dictionary[index]["volume"], label=label)
        plt.xlabel("Adsorption potential [kJ/mol]")
        plt.ylabel("Adsorption volume [ml/g]")

    plot_formats = {
        "isotherm": plot_isotherm,
        "isobar": plot_isobar,
        "isostere": plot_isostere,
        "enthalpy": plot_enthalpy,
        "characteristic": plot_characteristic,
        "langmuir": plot_isotherm,
        "n-langmuir": plot_isotherm,
        "bet": plot_isotherm,
        "anti-langmuir": plot_isotherm,
        "henry": plot_isotherm,
        "freundlich": plot_isotherm,
        "sips": plot_isotherm,
        "n-sips": plot_isotherm,
        "langmuir-freundlich": plot_isotherm,
        "n-langmuir-freundlich": plot_isotherm,
        "redlich-peterson": plot_isotherm,
        "toth": plot_isotherm,
        "unilan": plot_isotherm,
        "obrien-myers": plot_isotherm,
        "quadratic": plot_isotherm,
        "asymptotic-temkin": plot_isotherm,
        "bingel-walton": plot_isotherm
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

    plt.legend()

    if save == "yes":
        os.makedirs(name="Plots", exist_ok=True)
        figure_name = f"{input_dictionary[0]['ADSORBATE']}_in_{input_dictionary[0]['ADSORBENT']}_{plot_format}"
        plt.savefig(f"Plots/{figure_name}.png")


def write_data(source_dictionary: dict, properties_dictionary: dict, input_dictionary: dict, write_format: str) -> None:
    """
    Create data files based on the input data type. Supports between isotherm, isobar, and characteristic curve.

    :param source_dictionary: Dictionary containing the source data for the plotting.
    :param properties_dictionary: Dictionary containing the properties of the source data.
    :param input_dictionary: Dictionary containing the arguments found in the input file.
    :param write_format: Format of the file, given by the source data format. Can be isobar, isotherm or characteristic
    curve.
    """

    def write_isotherm(index, base_name) -> None:
        file_name = f"{base_name}_isotherm_{source_dictionary[index]['temperature']}K.dat"
        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write(f"#Pressure [{input_dictionary[0]['OUTPUT_PRESSURE_UNITS']}]         Loading [{input_dictionary[0]['OUTPUT_LOADING_UNITS']}] \n")
            cf_pressure = convert_output(
                input_dictionary[0]['OUTPUT_PRESSURE_UNITS'],
                molecular_mass=properties_dictionary["MOLECULAR_MASS"])

            cf_loading = convert_output(
                input_dictionary[0]['OUTPUT_LOADING_UNITS'],
                molecular_mass=properties_dictionary["MOLECULAR_MASS"])

            for pressure, loading in zip(source_dictionary[index]["pressure"], source_dictionary[index]["loading"]):
                file.write(f"{pressure * cf_pressure}         {loading * cf_loading} \n")

    def write_isobar(index, base_name) -> None:
        file_name = f"{base_name}_isobar_{source_dictionary[index]['pressure']}MPa.dat"
        cf_temperature = convert_output(
            input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS'],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        cf_loading = convert_output(
            input_dictionary[0]['OUTPUT_LOADING_UNITS'],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])
        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write(f"#Temperature [{input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS']}]         Loading [{input_dictionary[0]['OUTPUT_LOADING_UNITS']}] \n")
            for temperature, loading in zip(source_dictionary[index]["temperature"], source_dictionary[index]["loading"]):
                file.write(f"{temperature * cf_temperature}         {loading * cf_loading} \n")

    def write_isostere(index, base_name) -> None:
        file_name = f"{base_name}_isostere_{source_dictionary[index]['loading']}mgpg.dat"
        cf_temperature = convert_output(
            input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS'],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        cf_pressure = convert_output(
            input_dictionary[0]['OUTPUT_PRESSURE_UNITS'],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])
        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write(f"#Temperature [{input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS']}]         Pressure [{input_dictionary[0]['OUTPUT_PRESSURE_UNITS']}] \n")
            for temperature, pressure in zip(source_dictionary[index]["temperature"], source_dictionary[index]["pressure"]):
                file.write(f"{temperature * cf_temperature}         {pressure * cf_pressure} \n")

    def write_characteristic(index, base_name) -> None:
        if type(source_dictionary[index]["temperature"]) is not numpy.ndarray:
            condition = f"{source_dictionary[index]['temperature']}K"
        else:
            condition = f"{source_dictionary[index]['pressure']}MPa"

        file_name = f"{base_name}_characteristic_{condition}.dat"
        cf_potential = convert_output(
            input_dictionary[0]['OUTPUT_POTENTIAL_UNITS'],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])

        cf_volume = convert_output(
            input_dictionary[0]['OUTPUT_VOLUME_UNITS'],
            molecular_mass=properties_dictionary["MOLECULAR_MASS"])
        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write(f"#Potential [{input_dictionary[0]['OUTPUT_POTENTIAL_UNITS']}]         Volume [{input_dictionary[0]['OUTPUT_VOLUME_UNITS']}] \n")
            for potential, volume in zip(source_dictionary[index]["potential"], source_dictionary[index]["volume"]):
                file.write(f"{potential * cf_potential}         {volume * cf_volume} \n")

    def write_enthalpy(index, base_name) -> None:
        file_name = f"{base_name}_enthalpy_{source_dictionary[index]['source']}.dat"
        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write("#Loading [mg/g]         Enthalpy of adsorption [kJ/mol] \n")
            for loading, enthalpy in zip(source_dictionary[index]["loading"], source_dictionary[index]["enthalpy"]):
                file.write(f"{loading}         {enthalpy} \n")

    file_write_formats = {
        "isotherm": write_isotherm,
        "isobar": write_isobar,
        "isostere": write_isostere,
        "characteristic": write_characteristic,
        "enthalpy": write_enthalpy,
        "langmuir": write_isotherm,
        "n-langmuir": write_isotherm,
        "bet": write_isotherm,
        "anti-langmuir": write_isotherm,
        "henry": write_isotherm,
        "freundlich": write_isotherm,
        "sips": write_isotherm,
        "n-sips": write_isotherm,
        "langmuir-freundlich": write_isotherm,
        "n-langmuir-freundlich": write_isotherm,
        "redlich-peterson": write_isotherm,
        "toth": write_isotherm,
        "unilan": write_isotherm,
        "obrien-myers": write_isotherm,
        "quadratic": write_isotherm,
        "asymptotic-temkin": write_isotherm,
        "bingel-walton": write_isotherm
    }

    os.makedirs(name="Output", exist_ok=True)

    for index in source_dictionary:
        base_name = f"{input_dictionary[0]['ADSORBATE']}_in_{input_dictionary[0]['ADSORBENT']}"
        if write_format in file_write_formats:
            file_write_formats[write_format](index, base_name)


def predict_data(data_dictionary: dict, input_dictionary: dict, prediction_type: str,
                 properties_dictionary: dict) -> dict:
    """
    Predict adsorption properties based on the data from the first file entered in the input file.

    :param data_dictionary: Dictionary containing the data from the files.
    :param input_dictionary: Dictionary containing the arguments found in the input file.
    :param prediction_type: Type of predicted data. Can be isobar or isotherm
    :param properties_dictionary: Dictionary containing the properties of the molecule used.
    :return: Dictionary containing the predicted data.
    """

    def _get_pressure_boundaries(temperature: float, potential: numpy.ndarray) -> list:
        sat_pres = compute_saturation_pressure_from_method(
            method=input_dictionary[0]["ADSORBATE_SATURATION_PRESSURE"],
            temperature=temperature,
            properties_dictionary=properties_dictionary,
            saturation_pressure_file=input_dictionary[0]["SATURATION_PRESSURE_FILE"])

        minimum_pressure = physics.get_pressure(
            adsorption_potential=numpy.max(potential),
            saturation_pressure=sat_pres,
            temperature=temperature)

        maximum_pressure = physics.get_pressure(
            adsorption_potential=numpy.min(potential),
            saturation_pressure=sat_pres,
            temperature=temperature)

        return [minimum_pressure, maximum_pressure]

    def _get_temperature_boundaries(pressure: float, potential: numpy.ndarray) -> list:
        def minimum_temperature_function(temperature_guess: float) -> float:
            sat_pres = compute_saturation_pressure_from_method(
                method=input_dictionary[0]["ADSORBATE_SATURATION_PRESSURE"],
                temperature=temperature_guess,
                properties_dictionary=properties_dictionary,
                saturation_pressure_file=input_dictionary[0]["SATURATION_PRESSURE_FILE"])

            potential_computed = physics.get_adsorption_potential(
                temperature=temperature_guess,
                saturation_pressure=sat_pres,
                pressure=pressure)

            return numpy.min(potential) - potential_computed

        def maximum_temperature_function(temperature_guess: float) -> float:
            sat_pres = compute_saturation_pressure_from_method(
                method=input_dictionary[0]["ADSORBATE_SATURATION_PRESSURE"],
                temperature=temperature_guess,
                properties_dictionary=properties_dictionary,
                saturation_pressure_file=input_dictionary[0]["SATURATION_PRESSURE_FILE"])

            potential_computed = physics.get_adsorption_potential(
                temperature=temperature_guess,
                saturation_pressure=sat_pres,
                pressure=pressure)

            return numpy.max(potential) - potential_computed

        minimum_temperature = scipy.optimize.fsolve(minimum_temperature_function, x0=273)[0]
        maximum_temperature = scipy.optimize.fsolve(maximum_temperature_function, x0=273)[0]

        return [minimum_temperature, maximum_temperature]

    def _get_isostere_boundaries(loading: float, volume: numpy.ndarray) -> list:

        def minimum_temperature_function(temperature_guess: float) -> float:
            ads_dens = compute_density_from_method(
                method=input_dictionary[0]["ADSORBATE_DENSITY"],
                temperature=temperature_guess,
                properties_dictionary=properties_dictionary)

            volume_computed = physics.get_adsorption_volume(
                adsorbed_amount=loading,
                adsorbate_density=ads_dens)

            return numpy.max(volume) - volume_computed

        def maximum_temperature_function(temperature_guess: float) -> float:
            ads_dens = compute_density_from_method(
                method=input_dictionary[0]["ADSORBATE_DENSITY"],
                temperature=temperature_guess,
                properties_dictionary=properties_dictionary)

            volume_computed = physics.get_adsorption_volume(
                adsorbed_amount=loading,
                adsorbate_density=ads_dens)

            return numpy.min(volume) - volume_computed

        minimum_temperature = scipy.optimize.fsolve(minimum_temperature_function, x0=273)[0]
        maximum_temperature = scipy.optimize.fsolve(maximum_temperature_function, x0=273)[0]

        return [minimum_temperature, maximum_temperature]

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

            boundaries = _get_pressure_boundaries(
                temperature=temperature,
                potential=data_dictionary[0]["potential"])

            if (input_dictionary[0]["PREDICTION_PRESSURE_RANGE"] is not None and
                    boundaries[0] <= input_dictionary[0]["PREDICTION_PRESSURE_RANGE"][0] <= boundaries[1]):
                start_pressure = input_dictionary[0]["PREDICTION_PRESSURE_RANGE"][0]
            else:
                start_pressure = boundaries[0]

            if (input_dictionary[0]["PREDICTION_PRESSURE_RANGE"] is not None and
                    start_pressure <= input_dictionary[0]["PREDICTION_PRESSURE_RANGE"][1] <= boundaries[1]):
                end_pressure = input_dictionary[0]["PREDICTION_PRESSURE_RANGE"][1]
            else:
                end_pressure = boundaries[1]

            prediction_dictionary[index]["pressure"] = numpy.geomspace(
                start=start_pressure,
                stop=end_pressure,
                num=int(input_dictionary[0]["NUMBER_PRESSURE_POINTS"]))

            potential_range = physics.get_adsorption_potential(
                temperature=prediction_dictionary[index]["temperature"],
                saturation_pressure=prediction_dictionary[index]["saturation_pressure"],
                pressure=prediction_dictionary[index]["pressure"])

            prediction_dictionary[index]["loading"] = physics.get_adsorbed_amount(
                adsorption_volume=volume_interpolation_function(potential_range),
                adsorbate_density=prediction_dictionary[index]["density"])

        return prediction_dictionary

    def predict_isobar():

        prediction_dictionary = {}
        for index, pressure in enumerate(input_dictionary[0]["PREDICTION_PRESSURES"]):
            prediction_dictionary[index] = {}
            prediction_dictionary[index]["pressure"] = pressure
            boundaries = _get_temperature_boundaries(
                pressure=pressure,
                potential=data_dictionary[0]["potential"])

            if (input_dictionary[0]["PREDICTION_TEMPERATURE_RANGE"] is not None and
                    boundaries[0] <= input_dictionary[0]["PREDICTION_TEMPERATURE_RANGE"][0] <= boundaries[1]):
                start_temperature = input_dictionary[0]["PREDICTION_TEMPERATURE_RANGE"][0]
            else:
                start_temperature = boundaries[0]

            if (input_dictionary[0]["PREDICTION_TEMPERATURE_RANGE"] is not None and
                    start_temperature <= input_dictionary[0]["PREDICTION_TEMPERATURE_RANGE"][1] <= boundaries[1]):
                end_temperature = input_dictionary[0]["PREDICTION_TEMPERATURE_RANGE"][1]
            else:
                end_temperature = boundaries[1]

            prediction_dictionary[index]["temperature"] = numpy.linspace(
                start=start_temperature,
                stop=end_temperature,
                num=int(input_dictionary[0]["NUMBER_TEMPERATURE_POINTS"]))

            saturation_pressure_list = []
            density_list = []

            for temperature in prediction_dictionary[index]["temperature"]:
                saturation_pressure_list.append(compute_saturation_pressure_from_method(
                    method=input_dictionary[0]["ADSORBATE_SATURATION_PRESSURE"],
                    temperature=temperature,
                    properties_dictionary=properties_dictionary,
                    saturation_pressure_file=input_dictionary[0]["SATURATION_PRESSURE_FILE"]))

                density_list.append(compute_density_from_method(
                    method=input_dictionary[0]["ADSORBATE_DENSITY"],
                    temperature=temperature,
                    properties_dictionary=properties_dictionary))

            prediction_dictionary[index]["saturation_pressure"] = numpy.array(saturation_pressure_list)
            prediction_dictionary[index]["density"] = numpy.array(density_list)

            potential_range = physics.get_adsorption_potential(
                temperature=prediction_dictionary[index]["temperature"],
                saturation_pressure=prediction_dictionary[index]["saturation_pressure"],
                pressure=prediction_dictionary[index]["pressure"])

            prediction_dictionary[index]["loading"] = physics.get_adsorbed_amount(
                adsorption_volume=volume_interpolation_function(potential_range),
                adsorbate_density=prediction_dictionary[index]["density"])

        return prediction_dictionary

    def predict_isostere():

        prediction_dictionary = {}
        for index, loading in enumerate(input_dictionary[0]["PREDICTION_LOADINGS"]):
            prediction_dictionary[index] = {}
            prediction_dictionary[index]["loading"] = loading

            boundaries = _get_isostere_boundaries(
                loading=loading,
                volume=data_dictionary[0]["volume"])

            boundaries.sort()

            if (input_dictionary[0]["PREDICTION_ISOSTERE_RANGE"] is not None and
                    boundaries[0] <= input_dictionary[0]["PREDICTION_ISOSTERE_RANGE"][0] <= boundaries[1]):
                start_temperature = input_dictionary[0]["PREDICTION_ISOSTERE_RANGE"][0]
            else:
                start_temperature = boundaries[0]

            if (input_dictionary[0]["PREDICTION_ISOSTERE_RANGE"] is not None and
                    start_temperature <= input_dictionary[0]["PREDICTION_ISOSTERE_RANGE"][1] <= boundaries[1]):
                end_temperature = input_dictionary[0]["PREDICTION_ISOSTERE_RANGE"][1]
            else:
                end_temperature = boundaries[1]

            prediction_dictionary[index]["temperature"] = numpy.linspace(
                start=start_temperature,
                stop=end_temperature,
                num=int(input_dictionary[0]["NUMBER_ISOSTERE_POINTS"]))

            saturation_pressure_list = []
            density_list = []

            for temperature in prediction_dictionary[index]["temperature"]:
                saturation_pressure_list.append(compute_saturation_pressure_from_method(
                    method=input_dictionary[0]["ADSORBATE_SATURATION_PRESSURE"],
                    temperature=temperature,
                    properties_dictionary=properties_dictionary,
                    saturation_pressure_file=input_dictionary[0]["SATURATION_PRESSURE_FILE"]))

                density_list.append(compute_density_from_method(
                    method=input_dictionary[0]["ADSORBATE_DENSITY"],
                    temperature=temperature,
                    properties_dictionary=properties_dictionary))

            prediction_dictionary[index]["saturation_pressure"] = numpy.array(saturation_pressure_list)
            prediction_dictionary[index]["density"] = numpy.array(density_list)

            volume_range = physics.get_adsorption_volume(
                adsorbed_amount=loading,
                adsorbate_density=prediction_dictionary[index]["density"])

            prediction_dictionary[index]["pressure"] = physics.get_pressure(
                adsorption_potential=potential_interpolation_function(volume_range),
                saturation_pressure=prediction_dictionary[index]["saturation_pressure"],
                temperature=prediction_dictionary[index]["temperature"])

        return prediction_dictionary

    prediction_formats = {
        "isotherm": predict_isotherm,
        "isobar": predict_isobar,
        "isostere": predict_isostere
    }

    volume_interpolation_function = scipy.interpolate.interp1d(
        x=data_dictionary[0]["potential"],
        y=data_dictionary[0]["volume"],
        fill_value="extrapolate")

    potential_interpolation_function = scipy.interpolate.interp1d(
        x=data_dictionary[0]["volume"],
        y=data_dictionary[0]["potential"],
        fill_value="extrapolate")

    if prediction_type in prediction_formats:
        prediction_dictionary = prediction_formats[prediction_type]()
    else:
        raise ValueError(f"{prediction_type} is not a valid prediction type. Change the method or check "
                         f"for spelling errors!")

    return prediction_dictionary


def compute_adsorption_enthalpy(data_dictionary: dict, input_dictionary: dict) -> dict:

    def isotherm_enthalpy():
        enthalpy_dictionary = {}
        old_interpolation_function = None
        for index in data_dictionary:
            new_interpolation_function = scipy.interpolate.interp1d(
                x=data_dictionary[index]["loading"],
                y=data_dictionary[index]["pressure"],
                fill_value="extrapolate")

            if index > 0:
                enthalpy_dictionary[index-1] = {}
                new_temperature = data_dictionary[index]["temperature"]
                old_pressure = old_interpolation_function(data_dictionary[index]["loading"])
                new_pressure = data_dictionary[index]["pressure"]

                enthalpy_dictionary[index-1]["source"] = f"{old_temperature}K - {new_temperature}K"
                enthalpy_dictionary[index-1]["loading"] = data_dictionary[index]["loading"]
                enthalpy_dictionary[index-1]["enthalpy"] = physics.get_adsorption_enthalpy(
                    temperature_1=old_temperature,
                    temperature_2=new_temperature,
                    pressure_1=old_pressure,
                    pressure_2=new_pressure)

            old_temperature = data_dictionary[index]["temperature"]
            old_interpolation_function = new_interpolation_function

        return enthalpy_dictionary

    def isobar_enthalpy():
        enthalpy_dictionary = {}
        old_interpolation_function = None
        for index in data_dictionary:
            new_interpolation_function = scipy.interpolate.interp1d(
                x=data_dictionary[index]["loading"],
                y=data_dictionary[index]["temperature"],
                fill_value="extrapolate")

            enthalpy_dictionary[index] = {}
            if index > 0:
                new_pressure = data_dictionary[index]["pressure"]
                old_temperature = old_interpolation_function(data_dictionary[index]["loading"])
                new_temperature = data_dictionary[index]["temperature"]

                enthalpy_dictionary[index-1]["source"] = f"{old_pressure}MPa - {new_pressure}MPa"
                enthalpy_dictionary[index-1]["loading"] = data_dictionary[index]["loading"]
                enthalpy_dictionary[index-1]["enthalpy"] = physics.get_adsorption_enthalpy(
                    temperature_1=old_temperature,
                    temperature_2=new_temperature,
                    pressure_1=old_pressure,
                    pressure_2=new_pressure)

            old_pressure = data_dictionary[index]["pressure"]
            old_interpolation_function = new_interpolation_function

        return enthalpy_dictionary

    data_types = {
        "isotherm": isotherm_enthalpy,
        "isobar": isobar_enthalpy,
        "enthalpy": isotherm_enthalpy,
        "langmuir": isotherm_enthalpy,
        "n-langmuir": isotherm_enthalpy,
        "bet": isotherm_enthalpy,
        "anti-langmuir": isotherm_enthalpy,
        "henry": isotherm_enthalpy,
        "freundlich": isotherm_enthalpy,
        "sips": isotherm_enthalpy,
        "n-sips": isotherm_enthalpy,
        "langmuir-freundlich": isotherm_enthalpy,
        "n-langmuir-freundlich": isotherm_enthalpy,
        "redlich-peterson": isotherm_enthalpy,
        "toth": isotherm_enthalpy,
        "unilan": isotherm_enthalpy,
        "obrien-myers": isotherm_enthalpy,
        "quadratic": isotherm_enthalpy,
        "asymptotic-temkin": isotherm_enthalpy,
        "bingel-walton": isotherm_enthalpy
    }

    if input_dictionary[0]["DATA_TYPES"] in data_types:
        enthalpy_dictionary = data_types[input_dictionary[0]["DATA_TYPES"]]()
    else:
        raise ValueError(f"{input_dictionary[0]['DATA_TYPES']} is not a valid data type for enthalpy calculations. "
                         f"Change the method or check for spelling errors!")

    return enthalpy_dictionary


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


def show_plots() -> None:
    """
    Show all created plots simultaneously in separate windows.
    """
    plt.show()
