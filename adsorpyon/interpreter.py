
# Standard libraries
import os
import logging

# Local libraries
from adsorpyon import density
from adsorpyon import saturation_pressure
from adsorpyon import physics
from adsorpyon import input_reader
from adsorpyon import constants

# Third-party libraries
import matplotlib.pyplot as plt
import scipy.interpolate
import scipy.optimize
import numpy


logger = logging.getLogger(__name__)
logging.basicConfig(filename="cappa.log", filemode="w+", level=logging.INFO, datefmt="%d-%m-%y %H:%M:%S",
                    format="%(asctime)s %(levelname)s -> %(message)s")


def compute_density_from_method(method: str, temperature: float, properties_dictionary: dict,
                                input_dictionary: dict) -> float:
    """
    Compute the adsorbate density using the method specified in the input file and the respective molecular properties
    and environmental conditions.

    :param method: Name of the method used to compute the adsorbate density.
    :param temperature: Temperature at which the adsorbate density is computed in K.
    :param properties_dictionary: Dictionary containing the properties of the molecule used.
    :param input_dictionary: Dictionary containing the arguments found in the input file.
    :return: Adsorbate density in kg/m3.
    """

    logger.info(f"Computing density at {temperature} K using method {method}.")

    def density_empirical() -> float:
        return density.empirical(
            pressure_critical=properties_dictionary['PRESSURE_CRITICAL'],
            temperature_critical=properties_dictionary['TEMPERATURE_CRITICAL'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def density_hauer() -> float:
        return density.hauer(
            temperature=temperature,
            temperature_boiling=properties_dictionary['TEMPERATURE_BOILING'],
            density_boiling=properties_dictionary['DENSITY_BOILING'],
            thermal_expansion_coefficient=input_dictionary[0]['THERMAL_EXPANSION_COEFFICIENT'])

    def density_ozawa() -> float:
        return density.ozawa(
            temperature=temperature,
            temperature_boiling=properties_dictionary['TEMPERATURE_BOILING'],
            density_boiling=properties_dictionary['DENSITY_BOILING'],
            thermal_expansion_coefficient=input_dictionary[0]['THERMAL_EXPANSION_COEFFICIENT'])

    density_methods = {
        "empirical": density_empirical,
        "hauer": density_hauer,
        "ozawa": density_ozawa
    }

    if method in density_methods.keys():
        adsorbate_density = density_methods[method]()
        logger.info(f"Obtained density {adsorbate_density} kg/m3.")
    else:
        logger.error(f"{method} is not a valid adsorbate density computation method.")
        raise ValueError(f"{method} is not a valid adsorbate density computation method."
                         f" Change the method or check for spelling errors!")

    return adsorbate_density


def compute_saturation_pressure_from_method(method: str, temperature: float, properties_dictionary: dict,
                                            saturation_pressure_file: str, input_dictionary: dict) -> float:
    """
    Compute the adsorbate saturation pressure using the method specified in the input file and the respective molecular
    properties and environmental conditions.

    :param method: Name of the method used to compute the adsorbate saturation pressure.
    :param temperature: Temperature at which the adsorbate saturation pressure is computed in K.
    :param properties_dictionary: Dictionary containing the properties of the molecule used.
    :param saturation_pressure_file: Path to the file containing saturation pressure data.
    :param input_dictionary: Dictionary containing the arguments found in the input file.
    :return: Adsorbate saturation pressure in MPa.
    """

    logger.info(f"Computing saturation pressure at {temperature} K using method {method}.")

    def saturation_pressure_dubinin() -> float:
        return saturation_pressure.dubinin(
            temperature=temperature,
            temperature_critical=properties_dictionary['TEMPERATURE_CRITICAL'],
            pressure_critical=properties_dictionary['PRESSURE_CRITICAL'])

    def saturation_pressure_amankwah() -> float:
        return saturation_pressure.amankwah(
            temperature=temperature,
            temperature_critical=properties_dictionary['TEMPERATURE_CRITICAL'],
            pressure_critical=properties_dictionary['PRESSURE_CRITICAL'],
            k=input_dictionary[0]['AMANKWAH_EXPONENT'])

    def saturation_pressure_extrapolation() -> float:
        return saturation_pressure.extrapolation(temperature=temperature, file=saturation_pressure_file)

    def saturation_pressure_polynomial_water() -> float:
        return saturation_pressure.polynomial_water(temperature=temperature)

    def saturation_pressure_peng_robinson() -> float:
        return saturation_pressure.pengrobinson(
            temperature=temperature,
            temperature_critical=properties_dictionary['TEMPERATURE_CRITICAL'],
            pressure_critical=properties_dictionary['PRESSURE_CRITICAL'],
            pressure_guess=1,
            acentric_factor=properties_dictionary['ACENTRIC_FACTOR'])

    def saturation_pressure_preos_extrapolation() -> float:
        return saturation_pressure.equation_extrapolation(
            temperature=temperature,
            temperature_critical=properties_dictionary['TEMPERATURE_CRITICAL'],
            pressure_critical=properties_dictionary['PRESSURE_CRITICAL'],
            acentric_factor=properties_dictionary['ACENTRIC_FACTOR'],
            temperature_boiling=properties_dictionary['TEMPERATURE_BOILING'],
            equation="preos",
            kappa1=properties_dictionary['PRSV_KAPPA1'],
            kappa2=properties_dictionary['PRSV_KAPPA2'],
            kappa3=properties_dictionary['PRSV_KAPPA3'],
            function="polynomial2")

    def saturation_pressure_prsv1_extrapolation() -> float:
        return saturation_pressure.equation_extrapolation(
            temperature=temperature,
            temperature_critical=properties_dictionary['TEMPERATURE_CRITICAL'],
            pressure_critical=properties_dictionary['PRESSURE_CRITICAL'],
            acentric_factor=properties_dictionary['ACENTRIC_FACTOR'],
            temperature_boiling=properties_dictionary['TEMPERATURE_BOILING'],
            equation="prsv1",
            kappa1=properties_dictionary['PRSV_KAPPA1'],
            kappa2=properties_dictionary['PRSV_KAPPA2'],
            kappa3=properties_dictionary['PRSV_KAPPA3'],
            function="polynomial2")

    def saturation_pressure_prsv2_extrapolation() -> float:
        return saturation_pressure.equation_extrapolation(
            temperature=temperature,
            temperature_critical=properties_dictionary['TEMPERATURE_CRITICAL'],
            pressure_critical=properties_dictionary['PRESSURE_CRITICAL'],
            acentric_factor=properties_dictionary['ACENTRIC_FACTOR'],
            temperature_boiling=properties_dictionary['TEMPERATURE_BOILING'],
            equation="prsv2",
            kappa1=properties_dictionary['PRSV_KAPPA1'],
            kappa2=properties_dictionary['PRSV_KAPPA2'],
            kappa3=properties_dictionary['PRSV_KAPPA3'],
            function="polynomial2")

    def saturation_pressure_widom_banuti() -> float:
        return saturation_pressure.widombanuti(
            temperature=temperature,
            temperature_critical=properties_dictionary['TEMPERATURE_CRITICAL'],
            pressure_critical=properties_dictionary['PRESSURE_CRITICAL'],
            species_parameter=5.589,
            acentric_factor=properties_dictionary['ACENTRIC_FACTOR'])

    def saturation_pressure_isochore() -> float:
        return saturation_pressure.critical_isochore_model(
            temperature=temperature,
            temperature_critical=properties_dictionary['TEMPERATURE_CRITICAL'],
            pressure_critical=properties_dictionary['PRESSURE_CRITICAL'],
            acentric_factor=properties_dictionary['ACENTRIC_FACTOR'])

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
        "critical_isochore": saturation_pressure_isochore
    }

    if method in saturation_pressure_methods.keys():
        adsorbate_saturation_pressure = saturation_pressure_methods[method]()
        logger.info(f"Obtained saturation pressure {adsorbate_saturation_pressure} MPa.")
    else:
        logger.error(f"{method} is not a valid adsorbate saturation pressure computation method.")
        raise ValueError(f"{method} is not a valid adsorbate saturation "
                         f"pressure computation method. Change the method or check for spelling errors!")

    return adsorbate_saturation_pressure


def read_data(source_dictionary: dict, properties_dictionary: dict, input_dictionary: dict) -> None:
    """
    Read adsorption data from the provided files.

    Takes the list of data files provided in the input dictionary and uses input_reader.create_data_list() to parse it.
    The data files can either be two-column files or equation parameter files. For the former the data is stored as it,
    while for the latter the data is generated using the provided parameters.

    :param source_dictionary: Dictionary used to store the parsed data.
    :param properties_dictionary: Dictionary containing the properties of the molecule used.
    :param input_dictionary: Dictionary containing the arguments found in the input file.
    """

    logger.info(f"Starting reading data procedure.")

    def read_isotherm(index: int, file_data: list) -> None:
        """
        From a two-column file, store the first column as pressure and second column as loading.
        """
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']

        pressure = numpy.array([row[0] for row in file_data])
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        loading = numpy.array([row[1] for row in file_data])
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_isobar(index: int, file_data: list) -> None:
        """
        From a two-column file, store the first column as temperature and second column as loading.
        """
        source_dictionary[index]['pressure'] = input_dictionary[index]['PRESSURES']

        temperature = numpy.array([row[0] for row in file_data])
        source_dictionary[index]['temperature'] = temperature * convert_input(
            unit=source_dictionary[index]['TEMPERATURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        loading = numpy.array([row[1] for row in file_data])
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_isostere(index: int, file_data: list) -> None:
        """
        From a two-column file, store the first column as temperature and second column as pressure.
        """
        source_dictionary[index]['loading'] = input_dictionary[index]['LOADINGS']

        temperature = numpy.array([row[0] for row in file_data])
        source_dictionary[index]['temperature'] = temperature * convert_input(
            unit=source_dictionary[index]['TEMPERATURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        pressure = numpy.array([row[1] for row in file_data])
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_characteristic(index: int, file_data: list) -> None:
        """
        From a two-column file, store the first column as adsorption potential and second column as volume filling.
        """
        potential = numpy.array([row[0] for row in file_data])
        source_dictionary[index]['potential'] = potential * convert_input(
            unit=input_dictionary[index]['POTENTIAL_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        volume = numpy.array([row[1] for row in file_data])
        source_dictionary[index]['volume'] = volume * convert_input(
            unit=input_dictionary[index]['VOLUME_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_langmuir(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        bp = param2 * pressure
        loading = param1 * bp / (1 + bp)
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])


    def read_n_langmuir(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        params1 = file_data[1][0]
        params2 = file_data[2][0]
        if len(params1) != len(params2):
            raise ValueError("The number of parameters for the n-Langmuir isotherm does not match!")

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        loading = numpy.zeros(shape=num)
        for param1, param2 in zip(params1, params2):
            bp = param2 * pressure
            loading += param1 * bp / (1 + bp)

        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_bet(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]
        param4 = file_data[4][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        pr = pressure / param1
        bpr = param3 * pr
        loading = param2 * bpr / ((1 - param4 * pr) * (1 - param4 + bpr))
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_anti_langmuir(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        loading = (param1 * pressure) / (1 - param2 * pressure)
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_henry(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        loading = param1 * pressure
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_freundlich(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        loading = param1 * numpy.power(pressure, 1/param2)
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_sips(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        bpiv = numpy.power(param2 * pressure, 1/param3)
        loading = param1 * bpiv / (1 + bpiv)
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_n_sips(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        params1 = file_data[1]
        params2 = file_data[2]
        params3 = file_data[3]
        if len(params1) != len(params2) or len(params1) != len(params3):
            raise ValueError("The number of parameters for the n-Langmuir isotherm does not match!")

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        loading = numpy.zeros(shape=num)
        for param1, param2, param3 in zip(params1, params2, params3):
            bpiv = numpy.power(param2 * pressure, 1/param3)
            loading += param1 * bpiv / (1 + bpiv)

        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_langmuir_freundlich(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        bpv = param2 * numpy.power(pressure, param3)
        loading = param1 * bpv / (1 + bpv)

        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_n_langmuir_freundlich(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        params1 = file_data[1]
        params2 = file_data[2]
        params3 = file_data[3]
        if len(params1) != len(params2) or len(params1) != len(params3):
            raise ValueError("The number of parameters for the n-Langmuir isotherm does not match!")

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        loading = numpy.zeros(shape=num)
        for param1, param2, param3 in zip(params1, params2, params3):
            bpv = param2 * numpy.power(pressure, param3)
            loading += param1 * bpv / (1 + bpv)

        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_redlich_peterson(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        loading = param1 * pressure / (1 + param2 * numpy.power(pressure, param3))
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])


    def read_toth(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        bp = param2 * pressure
        loading = param1 * bp / numpy.power((1 + numpy.power(bp, param3)), 1/param3)
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_unilan(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]

        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        bp = param2 * pressure
        loading = (param1 / (2 * param3) * numpy.log((1 + bp * numpy.exp(param3)) / (1 + bp * numpy.exp(-param3))))
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_obrien_myers(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]
        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        bp = param2 * pressure
        loading = param1 * (bp / (1 + bp) + param3**2 * bp * (1 - bp) / (2 * numpy.power(1 + bp, 3)))
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_quadratic(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]
        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        bp = param2 * pressure
        cp2 = param3 * numpy.power(pressure, 2)
        loading = param1 * (bp + 2 * cp2 / (1 + bp + cp2))
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_asymptotic_temkin(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]
        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        bp = param2 * pressure
        frac = bp / (1 + bp)
        loading = param1 * frac + param1 * param3 * numpy.power(frac, 2) * (frac - 1)
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

    def read_bingel_walton(index: int, file_data: list) -> None:
        source_dictionary[index]['temperature'] = input_dictionary[index]['TEMPERATURES']
        start_pressure, stop_pressure, num = file_data[0][0], file_data[0][1], int(file_data[0][2])
        param1 = file_data[1][0]
        param2 = file_data[2][0]
        param3 = file_data[3][0]
        pressure = numpy.linspace(start=start_pressure, stop=stop_pressure, num=num)
        source_dictionary[index]['pressure'] = pressure * convert_input(
            unit=input_dictionary[index]['PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        term = -(param2 + param3) * pressure
        loading = param1 * (1 - numpy.exp(term) / (1 + param3/param2 * numpy.exp(term)))
        source_dictionary[index]['loading'] = loading * convert_input(
            unit=input_dictionary[index]['LOADING_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

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
        logger.info(f"Attempting to parse data from {input_dictionary[index]['DATA_FILES']}.")
        file_data = input_reader.create_data_list(path=input_dictionary[index]['DATA_FILES'])
        logger.info(f"Finished parsing data from {input_dictionary[index]['DATA_FILES']}.")
        if input_dictionary[index]["DATA_TYPES"] in data_types:
            logger.info(f"Found {input_dictionary[index]['DATA_FILES']} is {input_dictionary[index]['DATA_TYPES']}.")
            source_dictionary[index] = {}
            data_types[input_dictionary[index]['DATA_TYPES']](index, file_data)
            logger.info(f"Successfully stored data from {input_dictionary[index]['DATA_FILES']} at index {index}.")
        else:
            logger.error(f"{input_dictionary[index]['DATA_TYPES']} is not a valid data type.")
            raise ValueError(f"{input_dictionary[index]['DATA_TYPES']} is not a valid data type. Change the method or "
                             f"check for spelling errors!")


def compute_characteristic(source_dictionary: dict, input_dictionary: dict, properties_dictionary: dict) -> None:
    """
    Compute the characteristic curve from the input data.

    If the input data type is either isobar or isotherm, use the appropriate transformation to convert it to
    a characteristic curve. If the input data is already a characteristic curve, store the data as is. All data is
    stored in source_dictionary.

    :param source_dictionary: Dictionary used to store the parsed data.
    :param properties_dictionary: Dictionary containing the properties of the molecule used.
    :param input_dictionary: Dictionary containing the arguments found in the input file.
    """

    logger.info(f"Starting characteristic calculation procedure.")

    def from_isotherm(index):
        source_dictionary[index]['saturation_pressure'] = compute_saturation_pressure_from_method(
            method=input_dictionary[index]['ADSORBATE_SATURATION_PRESSURE'],
            temperature=source_dictionary[index]['temperature'],
            properties_dictionary=properties_dictionary,
            saturation_pressure_file=input_dictionary[index]['SATURATION_PRESSURE_FILE'],
            input_dictionary=input_dictionary)

        source_dictionary[index]['density'] = compute_density_from_method(
            method=input_dictionary[index]['ADSORBATE_DENSITY'],
            temperature=source_dictionary[index]['temperature'],
            properties_dictionary=properties_dictionary,
            input_dictionary=input_dictionary)

        source_dictionary[index]['potential'] = physics.get_adsorption_potential(
            temperature=source_dictionary[index]['temperature'],
            saturation_pressure=source_dictionary[index]['saturation_pressure'],
            pressure=source_dictionary[index]['pressure'])

        source_dictionary[index]['volume'] = physics.get_adsorption_volume(
            adsorbed_amount=source_dictionary[index]['loading'],
            adsorbate_density=source_dictionary[index]['density'])

    def from_isobar(index):
        saturation_pressure_array = []
        density_array = []
        for temperature in source_dictionary[index]['temperature']:
            saturation_pressure_array.append(compute_saturation_pressure_from_method(
                method=input_dictionary[index]['ADSORBATE_SATURATION_PRESSURE'],
                temperature=temperature,
                properties_dictionary=properties_dictionary,
                saturation_pressure_file=input_dictionary[index]['SATURATION_PRESSURE_FILE'],
                input_dictionary=input_dictionary))

            density_array.append(compute_density_from_method(
                method=input_dictionary[index]['ADSORBATE_DENSITY'],
                temperature=temperature,
                properties_dictionary=properties_dictionary,
                input_dictionary=input_dictionary))

        source_dictionary[index]['saturation_pressure'] = numpy.array(saturation_pressure_array)
        source_dictionary[index]['density'] = numpy.array(density_array)

        source_dictionary[index]['potential'] = physics.get_adsorption_potential(
            temperature=source_dictionary[index]['temperature'],
            saturation_pressure=source_dictionary[index]['saturation_pressure'],
            pressure=source_dictionary[index]['pressure'])

        source_dictionary[index]['volume'] = physics.get_adsorption_volume(
            adsorbed_amount=source_dictionary[index]['loading'],
            adsorbate_density=source_dictionary[index]['density'])

    def from_characteristic(index):
        """
        Do nothing since the data is already stored properly in source_dictionary.
        """
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
        if input_dictionary[index]['DATA_TYPES'] in data_types:
            logger.info(f"Attempting to compute characteristic for {input_dictionary[index]['DATA_FILES']}.")
            data_types[input_dictionary[index]['DATA_TYPES']](index)
            logger.info(f"Finished computing characteristic for {input_dictionary[index]['DATA_FILES']}.")
        else:
            logger.error(f"{input_dictionary[index]['DATA_TYPE']} is not a valid data type for characteristic "
                         f"calculations")
            raise ValueError(f"{input_dictionary[index]['DATA_TYPE']} is not a valid data type for characteristic "
                             f"calculations. Earlier checks should have raised this error. Something went "
                             f"terribly wrong!")


def compute_saturation_pressure_curve(input_dictionary: dict, properties_dictionary: dict) -> None:
    """
    Compute the saturation pressure curve for a given range and store it in a file in the Output folder.

    :param input_dictionary: Dictionary containing the arguments found in the input file.
    :param properties_dictionary: Dictionary containing the properties of the molecule used.
    """

    decimals = 4

    logger.info(f"Starting saturation pressure curve calculations.")

    start_temperature = input_dictionary[0]['SATURATION_PRESSURE_RANGE'][0]
    end_temperature = input_dictionary[0]['SATURATION_PRESSURE_RANGE'][1]
    num = input_dictionary[0]['NUMBER_SATURATION_PRESSURE_POINTS']
    logger.info(f"Found temperature interval {start_temperature}K - {end_temperature}K with {num} points in between.")

    temperatures = numpy.linspace(start_temperature, end_temperature, num)
    saturation_pressures = numpy.zeros(num)
    logger.info(f"Successfully generated temperature interval and saturation pressure variable.")

    for index, temperature in enumerate(temperatures):
        saturation_pressures[index] = compute_saturation_pressure_from_method(
            method=input_dictionary[0]['ADSORBATE_SATURATION_PRESSURE'],
            temperature=temperature,
            properties_dictionary=properties_dictionary,
            saturation_pressure_file=input_dictionary[0]['SATURATION_PRESSURE_FILE'],
            input_dictionary=input_dictionary)
        logger.info(f"For temperature {temperature}K got saturation pressure {saturation_pressures[index]} MPa.")

    molecule = input_dictionary[0]['ADSORBATE']

    file_name = f"{molecule}_saturation_pressure.dat"
    unit_pressure = input_dictionary[0]['OUTPUT_PRESSURE_UNITS']
    unit_temperature = input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS']

    logger.info(f"Starting file writing procedure for saturation pressure curve.")
    os.makedirs(name="Output", exist_ok=True)
    with open(file=f"Output/{file_name}", mode="w") as file:
        logger.info(f"Successfully created file Output/{file_name}.")

        file.write(f"# Temperature [{unit_temperature}] \t Saturation pressure [{unit_pressure}] \n")
        cf_pressure = convert_output(
            input_dictionary[0]['OUTPUT_PRESSURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        cf_temperature = convert_output(
            input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        for temperature, pressure in zip(temperatures, saturation_pressures):
            pressure = numpy.round(pressure * cf_pressure, decimals=decimals)
            temperature = numpy.round(temperature * cf_temperature, decimals=decimals)

            pressure = str(pressure).rjust(11)
            temperature = str(temperature).rjust(11)

            file.write(f"{temperature} \t {pressure} \n")

    logger.info(f"Finished writing to file Output/{file_name}.")


def compute_density_curve(input_dictionary: dict, properties_dictionary: dict):

    decimals = 4

    logger.info(f"Starting density curve calculations.")

    start_temperature = input_dictionary[0]['DENSITY_RANGE'][0]
    end_temperature = input_dictionary[0]['DENSITY_RANGE'][1]
    num = input_dictionary[0]['NUMBER_DENSITY_POINTS']
    logger.info(f"Found temperature interval {start_temperature}K - {end_temperature}K with {num} points in between.")

    temperatures = numpy.linspace(start_temperature, end_temperature, num)
    densities = numpy.zeros(num)
    logger.info(f"Successfully generated temperature interval and saturation pressure variable.")

    for index, temperature in enumerate(temperatures):
        densities[index] = compute_density_from_method(
            method=input_dictionary[0]['ADSORBATE_DENSITY'],
            temperature=temperature,
            properties_dictionary=properties_dictionary,
            input_dictionary=input_dictionary)
        logger.info(f"For temperature {temperature}K got density {densities[index]} kg/m3.")

    molecule = input_dictionary[0]['ADSORBATE']

    file_name = f"{molecule}_density.dat"
    unit_temperature = input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS']

    logger.info(f"Starting file writing procedure for density curve.")
    os.makedirs(name="Output", exist_ok=True)
    with open(file=f"Output/{file_name}", mode="w") as file:
        logger.info(f"Successfully created file Output/{file_name}.")

        file.write(f"# Temperature [{unit_temperature}] \t Density [kg/m3] \n")
        cf_density = convert_output(
            input_dictionary[0]['OUTPUT_DENSITY_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        cf_temperature = convert_output(
            input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS'],
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        for temperature, density in zip(temperatures, densities):
            density = numpy.round(density * cf_density, decimals=decimals)
            temperature = numpy.round(temperature * cf_temperature, decimals=decimals)

            density = str(density).rjust(11)
            temperature = str(temperature).rjust(11)

            file.write(f"{temperature} \t {density} \n")

    logger.info(f"Finished writing to file Output/{file_name}.")


def plot_data(source_dictionary: dict, input_dictionary: dict, properties_dictionary: dict, plot_format: str,
              save: str, from_input: bool) -> None:
    """
    Create plot based on the input data type. Supports between isotherm, isobar, and characteristic curve.

    :param source_dictionary: Dictionary containing the source data for the plotting.
    :param input_dictionary: Dictionary containing the arguments found in the input file.
    :param plot_format: Format of the plot, given by the source data format. Can be isobar, isotherm or characteristic
    curve.
    :param save: Dictates if the plot is saved. Saved if "yes", otherwise do not save.
    :param from_input: Dictates if the data comes from the input files, and sets a separate name for the plot.
    """

    logger.info(f"Starting plotting procedure.")

    def plot_isotherm(index):

        unit_temperature = input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS']
        unit_pressure = input_dictionary[0]['OUTPUT_PRESSURE_UNITS']
        unit_loading = input_dictionary[0]['OUTPUT_LOADING_UNITS']

        cf_temperature = convert_output(
            unit_temperature,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        cf_pressure = convert_output(
            unit_pressure,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        cf_loading = convert_output(
            unit_loading,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        temperature = source_dictionary[index]['temperature'] * cf_temperature
        pressure = source_dictionary[index]['pressure'] * cf_pressure
        loading = source_dictionary[index]['loading'] * cf_loading

        label = f"{temperature:.2f}{unit_temperature}"
        plt.scatter(pressure, loading, label=label)
        plt.xlabel(f"Pressure [{unit_pressure}]")
        plt.ylabel(f"Adsorbed amount [{unit_loading}]")

    def plot_isobar(index):

        unit_temperature = input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS']
        unit_pressure = input_dictionary[0]['OUTPUT_PRESSURE_UNITS']
        unit_loading = input_dictionary[0]['OUTPUT_LOADING_UNITS']

        cf_temperature = convert_output(
            unit_temperature,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        cf_pressure = convert_output(
            unit_pressure,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        cf_loading = convert_output(
            unit_loading,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        temperature = source_dictionary[index]['temperature'] * cf_temperature
        pressure = source_dictionary[index]['pressure'] * cf_pressure
        loading = source_dictionary[index]['loading'] * cf_loading

        label = f"{pressure:.2f} {unit_pressure}"
        plt.scatter(temperature, loading, label=label)
        plt.xlabel(f"Temperature [{unit_temperature}]")
        plt.ylabel(f"Adsorbed amount [{unit_loading}]")

    def plot_isostere(index):

        unit_temperature = input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS']
        unit_pressure = input_dictionary[0]['OUTPUT_PRESSURE_UNITS']
        unit_loading = input_dictionary[0]['OUTPUT_LOADING_UNITS']

        cf_temperature = convert_output(
            unit_temperature,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        cf_pressure = convert_output(
            unit_pressure,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        cf_loading = convert_output(
            unit_loading,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        temperature = source_dictionary[index]['temperature'] * cf_temperature
        pressure = source_dictionary[index]['pressure'] * cf_pressure
        loading = source_dictionary[index]['loading'] * cf_loading

        label = f"{loading:.2f} {unit_loading}"
        plt.scatter(temperature, pressure, label=label)
        plt.xlabel(f"Temperature [{unit_temperature}]")
        plt.ylabel(f"Pressure [{unit_pressure}]")

    def plot_enthalpy(index):

        unit_loading = input_dictionary[0]['OUTPUT_LOADING_UNITS']

        cf_loading = convert_output(
            unit_loading,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        loading = source_dictionary['loading'] * cf_loading
        enthalpy = source_dictionary['enthalpy']

        plt.scatter(loading, enthalpy)
        plt.xlabel(f"Loading [{unit_loading}]")
        plt.ylabel(f"Enthalpy of adsorption [kJ/mol]")

    def plot_characteristic(index):

        unit_temperature = input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS']
        unit_pressure = input_dictionary[0]['OUTPUT_PRESSURE_UNITS']
        unit_potential = input_dictionary[0]['OUTPUT_POTENTIAL_UNITS']
        unit_volume = input_dictionary[0]['OUTPUT_VOLUME_UNITS']

        cf_temperature = convert_output(
            unit_temperature,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        cf_pressure = convert_output(
            unit_pressure,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        cf_potential = convert_output(
            unit_potential,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        cf_volume = convert_output(
            unit_volume,
            molecular_mass=properties_dictionary['MOLECULAR_MASS'])

        temperature = source_dictionary[index]['temperature'] * cf_temperature
        pressure = source_dictionary[index]['pressure'] * cf_pressure
        potential = source_dictionary[index]['potential'] * cf_potential
        volume = source_dictionary[index]['volume'] * cf_volume

        if type(source_dictionary[index]['temperature']) is not numpy.ndarray:
            label = f"{temperature:.2f}{unit_temperature}"
        else:
            label = f"{pressure:2.f} {unit_pressure}"

        plt.scatter(potential, volume, label=label)
        plt.xlabel(f"Adsorption potential [{unit_potential}]")
        plt.ylabel(f"Adsorption volume [{unit_volume}]")

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

    plt.figure(figsize=(7, 6))
    plt.rc('axes', labelsize="xx-large")
    plt.rc('xtick', labelsize="xx-large")
    plt.rc('ytick', labelsize="xx-large")
    plt.rc('legend', fontsize="x-large")

    for index in source_dictionary:
        if plot_format in plot_formats:
            logger.info(f"Attempting to plot {plot_format} {index}.")
            plot_formats[plot_format](index)
            logger.info(f"Finished plotting {plot_format} {index}.")
        else:
            logger.error(f"{plot_format} at index {index} is not a valid data type for plotting!")
            raise ValueError(f"{plot_format} at index {index} is not a valid data type for plotting!")

    if plot_format == "isotherm" and input_dictionary[0]['LOGARITHMIC_PLOT'] == "yes":
        plt.xscale('log')

    plt.legend()
    plt.tight_layout()

    if save.lower() == "yes":
        os.makedirs(name="Plots", exist_ok=True)

        if from_input is True:
            figure_name = f"{input_dictionary[0]['ADSORBATE']}_in_{input_dictionary[0]['ADSORBENT']}_input"
        else:
            figure_name = f"{input_dictionary[0]['ADSORBATE']}_in_{input_dictionary[0]['ADSORBENT']}_{plot_format}"

        plt.savefig(f"Plots/{figure_name}.png")
        logger.info(f"Successfully saved plot at Plots/{figure_name}.png.")


def write_data(source_dictionary: dict, properties_dictionary: dict, input_dictionary: dict, write_format: str) -> None:
    """
    Create data files based on the input data type. Supports between isotherm, isobar, and characteristic curve.

    :param source_dictionary: Dictionary containing the source data for the plotting.
    :param properties_dictionary: Dictionary containing the properties of the source data.
    :param input_dictionary: Dictionary containing the arguments found in the input file.
    :param write_format: Format of the file, given by the source data format. Can be isobar, isotherm or characteristic
    curve.
    """

    decimals = 4
    logger.info(f"Starting writing procedure.")

    def write_isotherm(index, base_name) -> None:

        file_name = f"{base_name}_isotherm_{source_dictionary[index]['temperature']}K.dat"
        unit_pressure = input_dictionary[0]['OUTPUT_PRESSURE_UNITS']
        unit_loading = input_dictionary[0]['OUTPUT_LOADING_UNITS']

        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write(f"# Pressure [{unit_pressure}] \t Loading [{unit_loading}] \n")
            cf_pressure = convert_output(
                unit_pressure,
                molecular_mass=properties_dictionary['MOLECULAR_MASS'])

            cf_loading = convert_output(
                unit_loading,
                molecular_mass=properties_dictionary['MOLECULAR_MASS'])

            for pressure, loading in zip(source_dictionary[index]['pressure'], source_dictionary[index]['loading']):
                pressure = numpy.round(pressure * cf_pressure, decimals=decimals)
                loading = numpy.round(loading * cf_loading, decimals=decimals)

                pressure = str(pressure).rjust(11)
                loading = str(loading).rjust(11)

                file.write(f"{pressure} \t {loading} \n")

    def write_isobar(index, base_name) -> None:

        file_name = f"{base_name}_isobar_{source_dictionary[index]['pressure']}MPa.dat"
        unit_temperature = input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS']
        unit_loading = input_dictionary[0]['OUTPUT_LOADING_UNITS']

        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write(f"#Temperature [{unit_temperature}] \t Loading [{unit_loading}] \n")
            cf_temperature = convert_output(
                unit_temperature,
                molecular_mass=properties_dictionary['MOLECULAR_MASS'])

            cf_loading = convert_output(
                unit_loading,
                molecular_mass=properties_dictionary['MOLECULAR_MASS'])

            for temp, loading in zip(source_dictionary[index]['temperature'], source_dictionary[index]['loading']):
                temp = numpy.round(temp * cf_temperature, decimals=decimals)
                loading = numpy.round(loading * cf_loading, decimals=decimals)

                temp = str(temp).rjust(11)
                loading = str(loading).rjust(11)

                file.write(f"{temp} \t {loading} \n")


    def write_isostere(index, base_name) -> None:

        file_name = f"{base_name}_isostere_{source_dictionary[index]['loading']}mgpg.dat"
        unit_temperature = input_dictionary[0]['OUTPUT_TEMPERATURE_UNITS']
        unit_pressure = input_dictionary[0]['OUTPUT_PRESSURE_UNITS']

        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write(f"#Temperature [{unit_temperature}] \t Pressure [{unit_pressure}] \n")
            cf_temperature = convert_output(
                unit_temperature,
                molecular_mass=properties_dictionary['MOLECULAR_MASS'])

            cf_pressure = convert_output(
                unit_pressure,
                molecular_mass=properties_dictionary['MOLECULAR_MASS'])

            for temp, pressure in zip(source_dictionary[index]['temperature'], source_dictionary[index]['pressure']):
                temp = numpy.round(temp * cf_temperature, decimals=decimals)
                pressure = numpy.round(pressure * cf_pressure, decimals=decimals)

                temp = str(temp).rjust(11)
                pressure = str(pressure).rjust(11)

                file.write(f"{temp} \t {pressure} \n")


    def write_characteristic(index, base_name) -> None:

        if type(source_dictionary[index]['temperature']) is not numpy.ndarray:
            condition = f"{source_dictionary[index]['temperature']}K"
        else:
            condition = f"{source_dictionary[index]['pressure']}MPa"

        file_name = f"{base_name}_characteristic_{condition}.dat"
        unit_potential = input_dictionary[0]['OUTPUT_POTENTIAL_UNITS']
        unit_volume = input_dictionary[0]['OUTPUT_VOLUME_UNITS']

        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write(f"#Potential [{unit_potential}] \t Volume [{unit_volume}] \n")

            cf_potential = convert_output(
                unit_potential,
                molecular_mass=properties_dictionary['MOLECULAR_MASS'])

            cf_volume = convert_output(
                unit_volume,
                molecular_mass=properties_dictionary['MOLECULAR_MASS'])

            for potential, volume in zip(source_dictionary[index]['potential'], source_dictionary[index]['volume']):
                potential = numpy.round(potential * cf_potential, decimals=decimals)
                volume = numpy.round(volume * cf_volume, decimals=decimals)

                potential = str(potential).rjust(11)
                volume = str(volume).rjust(11)

                file.write(f"{potential} \t {volume} \n")

    def write_enthalpy(index, base_name) -> None:

        file_name = f"{base_name}_enthalpy.dat"
        unit_loading = input_dictionary[0]['OUTPUT_LOADING_UNITS']

        with open(file=f"Output/{file_name}", mode="w") as file:
            file.write(f"#Loading [{unit_loading}] \t Enthalpy of adsorption [kJ/mol] \n")

            cf_loading = convert_output(
                unit_loading,
                molecular_mass=properties_dictionary['MOLECULAR_MASS'])

            for loading, enthalpy in zip(source_dictionary['loading'], source_dictionary['enthalpy']):
                loading = numpy.round(loading * cf_loading, decimals=decimals)
                enthalpy = numpy.round(enthalpy, decimals=decimals)

                loading = str(loading).rjust(11)
                enthalpy = str(enthalpy).rjust(11)

                file.write(f"{loading} \t {enthalpy} \n")

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
            logger.info(f"Attempting to write {write_format} {index}.")
            file_write_formats[write_format](index, base_name)
            logger.info(f"Finished writing {write_format} {index}.")
        else:
            logger.error(f"{write_format} at index {index} is not a valid data type for writing!")
            raise ValueError(f"{write_format} at index {index} is not a valid data type for writing!")


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

    logger.info(f"Starting predicting procedure.")

    def _get_pressure_boundaries(temperature: float, potential: numpy.ndarray) -> list:
        logger.info(f"Computing pressure boundaries procedure.")
        sat_pres = compute_saturation_pressure_from_method(
            method=input_dictionary[0]['ADSORBATE_SATURATION_PRESSURE'],
            temperature=temperature,
            properties_dictionary=properties_dictionary,
            saturation_pressure_file=input_dictionary[0]['SATURATION_PRESSURE_FILE'],
            input_dictionary=input_dictionary)

        minimum_pressure = physics.get_pressure(
            adsorption_potential=numpy.max(potential),
            saturation_pressure=sat_pres,
            temperature=temperature)

        maximum_pressure = physics.get_pressure(
            adsorption_potential=numpy.min(potential),
            saturation_pressure=sat_pres,
            temperature=temperature)

        logger.info(f"Obtained pressure interval between {minimum_pressure} MPa and {maximum_pressure} MPa.")
        return [minimum_pressure, maximum_pressure]

    def _get_temperature_boundaries(pressure: float, potential: numpy.ndarray) -> list:
        logger.info(f"Computing temperature boundaries procedure.")
        def minimum_temperature_function(temperature_guess: float) -> float:
            sat_pres = compute_saturation_pressure_from_method(
                method=input_dictionary[0]['ADSORBATE_SATURATION_PRESSURE'],
                temperature=temperature_guess,
                properties_dictionary=properties_dictionary,
                saturation_pressure_file=input_dictionary[0]['SATURATION_PRESSURE_FILE'],
                input_dictionary=input_dictionary)

            potential_computed = physics.get_adsorption_potential(
                temperature=temperature_guess,
                saturation_pressure=sat_pres,
                pressure=pressure)

            return numpy.min(potential) - potential_computed

        def maximum_temperature_function(temperature_guess: float) -> float:
            sat_pres = compute_saturation_pressure_from_method(
                method=input_dictionary[0]['ADSORBATE_SATURATION_PRESSURE'],
                temperature=temperature_guess,
                properties_dictionary=properties_dictionary,
                saturation_pressure_file=input_dictionary[0]['SATURATION_PRESSURE_FILE'],
                input_dictionary=input_dictionary)

            potential_computed = physics.get_adsorption_potential(
                temperature=temperature_guess,
                saturation_pressure=sat_pres,
                pressure=pressure)

            return numpy.max(potential) - potential_computed

        minimum_temperature = scipy.optimize.fsolve(minimum_temperature_function, x0=273)[0]
        maximum_temperature = scipy.optimize.fsolve(maximum_temperature_function, x0=273)[0]

        logger.info(f"Obtained temperature interval between {minimum_temperature} MPa and {maximum_temperature} MPa.")
        return [minimum_temperature, maximum_temperature]

    def _get_isostere_boundaries(loading: float, volume: numpy.ndarray) -> list:
        logger.info(f"Computing isostere boundaries procedure.")
        def minimum_temperature_function(temperature_guess: float) -> float:
            ads_dens = compute_density_from_method(
                method=input_dictionary[0]['ADSORBATE_DENSITY'],
                temperature=temperature_guess,
                properties_dictionary=properties_dictionary,
                input_dictionary=input_dictionary)

            volume_computed = physics.get_adsorption_volume(
                adsorbed_amount=loading,
                adsorbate_density=ads_dens)

            return numpy.max(volume) - volume_computed

        def maximum_temperature_function(temperature_guess: float) -> float:
            ads_dens = compute_density_from_method(
                method=input_dictionary[0]['ADSORBATE_DENSITY'],
                temperature=temperature_guess,
                properties_dictionary=properties_dictionary,
                input_dictionary=input_dictionary)

            volume_computed = physics.get_adsorption_volume(
                adsorbed_amount=loading,
                adsorbate_density=ads_dens)

            return numpy.min(volume) - volume_computed

        minimum_temperature = scipy.optimize.fsolve(minimum_temperature_function, x0=273)[0]
        maximum_temperature = scipy.optimize.fsolve(maximum_temperature_function, x0=273)[0]

        logger.info(f"Obtained isostere interval between {minimum_temperature} MPa and {maximum_temperature} MPa.")
        return [minimum_temperature, maximum_temperature]

    def predict_isotherm():
        logger.info(f"Starting isotherm prediction procedure.")

        prediction_dictionary = {}
        for index, temperature in enumerate(input_dictionary[0]['PREDICTION_TEMPERATURES']):
            logger.info(f"Predicting isotherm at {temperature} K.")

            prediction_dictionary[index] = {}
            prediction_dictionary[index]['temperature'] = temperature
            prediction_dictionary[index]['saturation_pressure'] = compute_saturation_pressure_from_method(
                method=input_dictionary[0]['ADSORBATE_SATURATION_PRESSURE'],
                temperature=temperature,
                properties_dictionary=properties_dictionary,
                saturation_pressure_file=input_dictionary[0]['SATURATION_PRESSURE_FILE'],
                input_dictionary=input_dictionary)

            prediction_dictionary[index]['density'] = compute_density_from_method(
                method=input_dictionary[0]['ADSORBATE_DENSITY'],
                temperature=temperature,
                properties_dictionary=properties_dictionary,
                input_dictionary=input_dictionary)

            boundaries = _get_pressure_boundaries(
                temperature=temperature,
                potential=data_dictionary[0]['potential'])

            if (input_dictionary[0]['PREDICTION_PRESSURE_RANGE'] is not None and
                    boundaries[0] <= input_dictionary[0]['PREDICTION_PRESSURE_RANGE'][0] <= boundaries[1]):
                start_pressure = input_dictionary[0]['PREDICTION_PRESSURE_RANGE'][0]
            else:
                start_pressure = boundaries[0]

            if (input_dictionary[0]['PREDICTION_PRESSURE_RANGE'] is not None and
                    start_pressure <= input_dictionary[0]['PREDICTION_PRESSURE_RANGE'][1] <= boundaries[1]):
                end_pressure = input_dictionary[0]['PREDICTION_PRESSURE_RANGE'][1]
            else:
                end_pressure = boundaries[1]

            prediction_dictionary[index]['pressure'] = numpy.geomspace(
                start=start_pressure,
                stop=end_pressure,
                num=int(input_dictionary[0]['NUMBER_PRESSURE_POINTS']))

            potential_range = physics.get_adsorption_potential(
                temperature=prediction_dictionary[index]['temperature'],
                saturation_pressure=prediction_dictionary[index]['saturation_pressure'],
                pressure=prediction_dictionary[index]['pressure'])

            prediction_dictionary[index]['loading'] = physics.get_adsorbed_amount(
                adsorption_volume=volume_interpolation_function(potential_range),
                adsorbate_density=prediction_dictionary[index]['density'])
        return prediction_dictionary

    def predict_isobar():
        logger.info(f"Starting isobar prediction procedure.")

        prediction_dictionary = {}
        for index, pressure in enumerate(input_dictionary[0]['PREDICTION_PRESSURES']):
            logger.info(f"Predicting isobar at {pressure} MPa.")

            prediction_dictionary[index] = {}
            prediction_dictionary[index]['pressure'] = pressure
            boundaries = _get_temperature_boundaries(
                pressure=pressure,
                potential=data_dictionary[0]['potential'])

            if (input_dictionary[0]['PREDICTION_TEMPERATURE_RANGE'] is not None and
                    boundaries[0] <= input_dictionary[0]['PREDICTION_TEMPERATURE_RANGE'][0] <= boundaries[1]):
                start_temperature = input_dictionary[0]['PREDICTION_TEMPERATURE_RANGE'][0]
            else:
                start_temperature = boundaries[0]

            if (input_dictionary[0]['PREDICTION_TEMPERATURE_RANGE'] is not None and
                    start_temperature <= input_dictionary[0]['PREDICTION_TEMPERATURE_RANGE'][1] <= boundaries[1]):
                end_temperature = input_dictionary[0]['PREDICTION_TEMPERATURE_RANGE'][1]
            else:
                end_temperature = boundaries[1]

            prediction_dictionary[index]['temperature'] = numpy.linspace(
                start=start_temperature,
                stop=end_temperature,
                num=int(input_dictionary[0]['NUMBER_TEMPERATURE_POINTS']))

            saturation_pressure_list = []
            density_list = []

            for temperature in prediction_dictionary[index]['temperature']:
                saturation_pressure_list.append(compute_saturation_pressure_from_method(
                    method=input_dictionary[0]['ADSORBATE_SATURATION_PRESSURE'],
                    temperature=temperature,
                    properties_dictionary=properties_dictionary,
                    saturation_pressure_file=input_dictionary[0]['SATURATION_PRESSURE_FILE'],
                    input_dictionary=input_dictionary))

                density_list.append(compute_density_from_method(
                    method=input_dictionary[0]['ADSORBATE_DENSITY'],
                    temperature=temperature,
                    properties_dictionary=properties_dictionary,
                    input_dictionary=input_dictionary))

            prediction_dictionary[index]['saturation_pressure'] = numpy.array(saturation_pressure_list)
            prediction_dictionary[index]['density'] = numpy.array(density_list)

            potential_range = physics.get_adsorption_potential(
                temperature=prediction_dictionary[index]['temperature'],
                saturation_pressure=prediction_dictionary[index]['saturation_pressure'],
                pressure=prediction_dictionary[index]['pressure'])

            prediction_dictionary[index]['loading'] = physics.get_adsorbed_amount(
                adsorption_volume=volume_interpolation_function(potential_range),
                adsorbate_density=prediction_dictionary[index]['density'])

        return prediction_dictionary

    def predict_isostere():
        logger.info(f"Starting isostere prediction procedure.")

        prediction_dictionary = {}
        for index, loading in enumerate(input_dictionary[0]['PREDICTION_LOADINGS']):
            logger.info(f"Predicting isostere at {loading} mg/g.")

            prediction_dictionary[index] = {}
            prediction_dictionary[index]['loading'] = loading

            boundaries = _get_isostere_boundaries(
                loading=loading,
                volume=data_dictionary[0]['volume'])

            boundaries.sort()

            if (input_dictionary[0]['PREDICTION_ISOSTERE_RANGE'] is not None and
                    boundaries[0] <= input_dictionary[0]['PREDICTION_ISOSTERE_RANGE'][0] <= boundaries[1]):
                start_temperature = input_dictionary[0]['PREDICTION_ISOSTERE_RANGE'][0]
            else:
                start_temperature = boundaries[0]

            if (input_dictionary[0]['PREDICTION_ISOSTERE_RANGE'] is not None and
                    start_temperature <= input_dictionary[0]['PREDICTION_ISOSTERE_RANGE'][1] <= boundaries[1]):
                end_temperature = input_dictionary[0]['PREDICTION_ISOSTERE_RANGE'][1]
            else:
                end_temperature = boundaries[1]

            prediction_dictionary[index]['temperature'] = numpy.linspace(
                start=start_temperature,
                stop=end_temperature,
                num=int(input_dictionary[0]['NUMBER_ISOSTERE_POINTS']))

            saturation_pressure_list = []
            density_list = []

            for temperature in prediction_dictionary[index]['temperature']:
                saturation_pressure_list.append(compute_saturation_pressure_from_method(
                    method=input_dictionary[0]['ADSORBATE_SATURATION_PRESSURE'],
                    temperature=temperature,
                    properties_dictionary=properties_dictionary,
                    saturation_pressure_file=input_dictionary[0]['SATURATION_PRESSURE_FILE'],
                    input_dictionary=input_dictionary))

                density_list.append(compute_density_from_method(
                    method=input_dictionary[0]['ADSORBATE_DENSITY'],
                    temperature=temperature,
                    properties_dictionary=properties_dictionary,
                    input_dictionary=input_dictionary))

            prediction_dictionary[index]['saturation_pressure'] = numpy.array(saturation_pressure_list)
            prediction_dictionary[index]['density'] = numpy.array(density_list)

            volume_range = physics.get_adsorption_volume(
                adsorbed_amount=loading,
                adsorbate_density=prediction_dictionary[index]['density'])

            prediction_dictionary[index]['pressure'] = physics.get_pressure(
                adsorption_potential=potential_interpolation_function(volume_range),
                saturation_pressure=prediction_dictionary[index]['saturation_pressure'],
                temperature=prediction_dictionary[index]['temperature'])

        return prediction_dictionary

    prediction_formats = {
        "isotherm": predict_isotherm,
        "isobar": predict_isobar,
        "isostere": predict_isostere
    }

    volume_interpolation_function = scipy.interpolate.interp1d(
        x=data_dictionary[0]['potential'],
        y=data_dictionary[0]['volume'],
        fill_value="extrapolate")

    potential_interpolation_function = scipy.interpolate.CubicSpline(
        x=data_dictionary[0]['volume'],
        y=data_dictionary[0]['potential'],
        extrapolate=True)

    if prediction_type in prediction_formats:
        prediction_dictionary = prediction_formats[prediction_type]()
    else:
        raise ValueError(f"{prediction_type} is not a valid prediction type. Change the method or check "
                         f"for spelling errors!")

    return prediction_dictionary


def compute_adsorption_enthalpy(data_dictionary: dict, input_dictionary: dict, properties_dictionary: dict) -> dict:

    def _get_isostere_boundaries(loading: float, volume: numpy.ndarray) -> list:

        def minimum_temperature_function(temperature_guess: float) -> float:
            ads_dens = compute_density_from_method(
                method=input_dictionary[0]['ADSORBATE_DENSITY'],
                temperature=temperature_guess,
                properties_dictionary=properties_dictionary,
                input_dictionary=input_dictionary)

            volume_computed = physics.get_adsorption_volume(
                adsorbed_amount=loading,
                adsorbate_density=ads_dens)

            return numpy.max(volume) - volume_computed

        def maximum_temperature_function(temperature_guess: float) -> float:
            ads_dens = compute_density_from_method(
                method=input_dictionary[0]['ADSORBATE_DENSITY'],
                temperature=temperature_guess,
                properties_dictionary=properties_dictionary,
                input_dictionary=input_dictionary)

            volume_computed = physics.get_adsorption_volume(
                adsorbed_amount=loading,
                adsorbate_density=ads_dens)

            return numpy.min(volume) - volume_computed

        minimum_temperature = scipy.optimize.fsolve(minimum_temperature_function, x0=273)[0]
        maximum_temperature = scipy.optimize.fsolve(maximum_temperature_function, x0=273)[0]

        return [minimum_temperature, maximum_temperature]
    
    potential_interpolation_function = scipy.interpolate.PchipInterpolator(
        x=data_dictionary[0]['volume'],
        y=data_dictionary[0]['potential'])

    loadings = numpy.linspace(
        start=input_dictionary[0]['ENTHALPY_RANGE'][0],
        stop=input_dictionary[0]['ENTHALPY_RANGE'][1],
        num=int(input_dictionary[0]['NUMBER_ENTHALPY_POINTS']))

    enthalpy_dictionary = {"loading": loadings, "enthalpy": []}

    prediction_dictionary = {}
    plt.figure()
    for index, loading in enumerate(loadings):
        prediction_dictionary[index] = {}
        prediction_dictionary[index]['loading'] = loading

        boundaries = _get_isostere_boundaries(
            loading=loading,
            volume=data_dictionary[0]['volume'])

        boundaries.sort()
        boundaries[0] = max(boundaries[0], 0)

        prediction_dictionary[index]['temperature'] = numpy.linspace(
            start=input_dictionary[0]['ENTHALPY_TEMPERATURE_RANGE'][0],
            stop=input_dictionary[0]['ENTHALPY_TEMPERATURE_RANGE'][1],
            num=3)

        saturation_pressure_list = []
        density_list = []

        for temperature in prediction_dictionary[index]['temperature']:
            saturation_pressure_list.append(compute_saturation_pressure_from_method(
                method=input_dictionary[0]['ADSORBATE_SATURATION_PRESSURE'],
                temperature=temperature,
                properties_dictionary=properties_dictionary,
                saturation_pressure_file=input_dictionary[0]['SATURATION_PRESSURE_FILE'],
                input_dictionary=input_dictionary))

            density_list.append(compute_density_from_method(
                method=input_dictionary[0]['ADSORBATE_DENSITY'],
                temperature=temperature,
                properties_dictionary=properties_dictionary,
                input_dictionary=input_dictionary))

        prediction_dictionary[index]['saturation_pressure'] = numpy.array(saturation_pressure_list)
        prediction_dictionary[index]['density'] = numpy.array(density_list)

        volume_range = physics.get_adsorption_volume(
            adsorbed_amount=loading,
            adsorbate_density=prediction_dictionary[index]['density'])

        prediction_dictionary[index]['pressure'] = physics.get_pressure(
            adsorption_potential=potential_interpolation_function(volume_range),
            saturation_pressure=prediction_dictionary[index]['saturation_pressure'],
            temperature=prediction_dictionary[index]['temperature'])

        pressures = []
        temperatures = []
        for press, temp in zip(prediction_dictionary[index]['pressure'], prediction_dictionary[index]['temperature']):
            if press <= 0 or numpy.isnan(press) or numpy.isnan(temp):
                continue
            else:
                pressures.append(press)
                temperatures.append(temp)

        prediction_dictionary[index]['pressure'] = numpy.log(pressures)
        prediction_dictionary[index]['temperature'] = numpy.divide(1, temperatures)

        plt.scatter(prediction_dictionary[index]['temperature'], prediction_dictionary[index]['pressure']/prediction_dictionary[index]['pressure'][0])

        def fit_function(itemperature, heat, offset):
            return heat * 1000 / constants.GAS_CONSTANT * itemperature + offset

        opt, cvt = scipy.optimize.curve_fit(
            fit_function,
            prediction_dictionary[index]['temperature'],
            prediction_dictionary[index]['pressure'])

        plt.plot(prediction_dictionary[index]['temperature'],
                 fit_function(prediction_dictionary[index]['temperature'], *opt)/prediction_dictionary[index]['pressure'][0])

        enthalpy_dictionary['enthalpy'].append(-opt[0])

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
    elif unit == "atm":
        conversion_factor = 0.09869232667160
    elif unit == "Torr":
        conversion_factor = 0.000133322
    elif unit == "mmHg":
        conversion_factor = 133.322 * 0.000001

    # Temperature
    elif unit == "K":
        conversion_factor = 1
    elif unit == "R":
        conversion_factor = 1.8

    # Adsorbed amount
    elif unit in ["mg/g", "g/kg"]:
        conversion_factor = 1
    elif unit in ["mol/kg", "mmol/g"]:
        conversion_factor = molecular_mass
    elif unit in ["cm3/g", "mm3/mg", "dm3/kg", "l/kg", "ml/g"]:
        conversion_factor = molecular_mass / 22.4139757476

    # Adsorption potential
    elif unit == "kJ/mol":
        conversion_factor = 1
    elif unit == "J/mol":
        conversion_factor = 0.001

    # Adsorption volume
    elif unit in ["ml/g", "l/kg", "cm3/g", "dm3/kg"]:
        conversion_factor = 1

    # Density
    elif unit == "kg/m3":
        conversion_factor = 1

    # Not a recognized unit
    else:
        logger.error(f"Input unit {unit} is not a recognized unit!")
        raise ValueError(f"Input unit {unit} is not a recognized unit!")

    logger.info(f"Used conversion factor {conversion_factor} for input {unit}.")
    return conversion_factor


def convert_output(unit: str, molecular_mass: float) -> float:
    """
    Returns a conversion factor for the standard units to the user defined ones.
    :param unit: The unit of the input data.
    :param molecular_mass: The molecular mass of the molecule.
    :return: A number that the input is multiplied with to be converted to the intended unit.
    """

    logging.disable(level=logging.INFO)
    conversion_factor = 1 / convert_input(unit, molecular_mass)
    logging.disable(level=logging.NOTSET)

    logger.info(f"Used conversion factor {conversion_factor} for output {unit}.")
    return conversion_factor


def show_plots() -> None:
    """
    Show all created plots simultaneously in separate windows.
    """
    logger.info(f"Displaying the plots.")
    plt.show()
