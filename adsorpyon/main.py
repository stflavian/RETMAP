"""Main

This represents the main script.

Usage: run
"""

import json

import saturation_pressure
import density
import utils
import physics

from numpy import genfromtxt, nan


def get_saturation_pressure(data: dict, index: int, adsorbate_data: dict) -> float:
    if data[index]["saturation_pressure_function"] == "dubinin":
        return saturation_pressure.dubinin(temperature=data[index]["temperature"],
                                           temperature_critical=adsorbate_data["temperature_critical"],
                                           pressure_critical=adsorbate_data["pressure_critical"])
    elif data[index]["saturation_pressure_function"] == "amankwah":
        return saturation_pressure.amankwah(temperature=data[index]["temperature"],
                                            temperature_critical=adsorbate_data["temperature_critical"],
                                            pressure_critical=adsorbate_data["pressure_critical"],
                                            k=adsorbate_data["amankwah_exponent"])
    elif data[index]["saturation_pressure_function"] == "extrapolation":
        return saturation_pressure.extrapolation(temperature=data[index]["temperature"],
                                                 file=data[index]["saturation_pressure_file"])
    elif data[index]["saturation_pressure_function"] == "polynomial":
        return saturation_pressure.polynomial_water(temperature=data[index]["temperature"])
    elif data[index]["saturation_pressure_function"] == "peng-robinson":
        return saturation_pressure.pengrobinson(temperature=data[index]["temperature"],
                                                temperature_critical=adsorbate_data["temperature_critical"],
                                                pressure_critical=adsorbate_data["pressure_critical"],
                                                pressure_guess=1,
                                                acentric_factor=adsorbate_data["acentric_factor"])
    elif data[index]["saturation_pressure_function"] == "preos-extrapolation":
        return saturation_pressure.preos_extrapolation(temperature=data[index]["temperature"],
                                                       temperature_critical=adsorbate_data["temperature_critical"],
                                                       pressure_critical=adsorbate_data["pressure_critical"],
                                                       acentric_factor=adsorbate_data["acentric_factor"],
                                                       temperature_boiling=adsorbate_data["temperature_boiling"])
    else:
        raise ValueError(f"Input string {data[index]['saturation_pressure_function']} does not correspond to a valid"
                         f"method of determining the saturation pressure")


def get_density(data: dict, index: int, adsorbate_data: dict) -> float:
    if data[index]["density_function"] == "ozawa":
        return density.ozawa(temperature=data[index]["temperature"],
                             temperature_boiling=adsorbate_data["temperature_boiling"],
                             density_boiling=adsorbate_data["density_boiling"])
    elif data[index]["density_function"] == "empirical":
        return density.empirical(pressure_critical=adsorbate_data["pressure_critical"],
                                 temperature_critical=adsorbate_data["temperature_critical"],
                                 molecular_mass=adsorbate_data["molecular_mass"])
    elif data[index]["density_function"] == "hauer":
        return density.hauer(temperature=data[index]["temperature"],
                             temperature_boiling=adsorbate_data["temperature_boiling"],
                             density_boiling=adsorbate_data["density_boiling"],
                             thermal_expansion_coefficient=adsorbate_data["thermal_expansion_coefficient"])
    elif data[index]["density_function"] == "ozawa-modified":
        return density.ozawa_modified(temperature=data[index]["temperature"],
                                      temperature_boiling=adsorbate_data["temperature_boiling"],
                                      density_boiling=adsorbate_data["density_boiling"],
                                      thermal_expansion_coefficient=adsorbate_data["thermal_expansion_coefficient"])
    else:
        raise ValueError(f"Input string {data[index]['density_function']} does not correspond to a valid"
                         f"method of determining the density")


def main():

    # Read the config and adsorbate data files
    config = json.load(open("config.json"))
    adsorbate_data = json.load(open(config["adsorbate_data_file"]))

    # Create the data array and assign values
    data = {}
    for index, file in enumerate(config["data_files"]):
        data[index] = {}

        # Assign data from the config file
        data[index]["data_file"] = file
        data[index]["name"] = config["data_names"][index]
        data[index]["temperature"] = config["temperatures"][index]
        data[index]["data_type"] = config["data_types"][index]
        data[index]["saturation_pressure_function"] = config["adsorbate_saturation_pressure"]
        if data[index]["saturation_pressure_function"] == "extrapolation":
            data[index]["saturation_pressure_file"] = config["saturation_pressure_file"]
        data[index]["density_function"] = config["adsorbate_density"]
        data[index]["adsorbate_data_file"] = config["adsorbate_data_file"]

        # Calculate the saturation pressure and density based on the desired functions
        data[index]["saturation_pressure"] = get_saturation_pressure(data=data, index=index, adsorbate_data=adsorbate_data)
        data[index]["density"] = get_density(data=data, index=index, adsorbate_data=adsorbate_data)

        # Read the input data and calculate the isotherm or characteristic curve depending on data type
        file_data = genfromtxt(fname=file, filling_values=nan)
        if data[index]["data_type"] == "isotherm":
            data[index]["pressure"] = file_data[:, 0] * utils.convert_input(unit=config["column1_units"][index],
                                                                            adsorbate_data=adsorbate_data)  # [MPa]

            data[index]["adsorbed_amount"] = file_data[:, 1] * utils.convert_input(unit=config["column2_units"][index],
                                                                                   adsorbate_data=adsorbate_data)  # [mg/g]

            data[index]["adsorption_potential"] = physics.get_adsorption_potential(temperature=data[index]["temperature"],
                                                                                   saturation_pressure=data[index]["saturation_pressure"],
                                                                                   pressure=data[index]["pressure"])  # [kJ/mol]

            data[index]["adsorption_volume"] = physics.get_adsorption_volume(adsorbed_amount=data[index]["adsorbed_amount"],
                                                                             adsorbate_density=data[index]["density"])  # [ml/g]
        if data[index]["data_type"] == "characteristic curve":
            data[index]["adsorption_potential"] = file_data[:, 0] * utils.convert_input(unit=config["column1_units"][index],
                                                                                        adsorbate_data=adsorbate_data)  # [kJ/mol]

            data[index]["adsorption_volume"] = file_data[:, 1] * utils.convert_input(unit=config["column2_units"][index],
                                                                                     adsorbate_data=adsorbate_data)  # [ml/g]

            data[index]["pressure"] = physics.get_pressure(adsorption_potential=data[index]["adsorption_potential"],
                                                           saturation_pressure=data[index]["saturation_pressure"],
                                                           temperature=data[index]["temperature"])  # [MPa]

            data[index]["adsorbed_amount"] = physics.get_adsorbed_amount(adsorption_volume=data[index]["adsorption_volume"],
                                                                         adsorbate_density=data[index]["density"])  # [mg/g]
            
        if config["write_results"][index].lower() == "yes":
            utils.write_data(data=data, index=index)

    plot_presence = False
    if config["show_isotherm"].lower() == "yes":
        plot_presence = True
        utils.plot_isotherm(data=data, logarithmic=config["logarithmic_isotherm"], save=config["save_isotherm"])

    if config["show_characteristic_curve"].lower() == "yes":
        plot_presence = True
        utils.plot_characteristic_curve(data=data, save=config["save_characteristic_curve"])

    if plot_presence is True:
        utils.show_plots()


if __name__ == "__main__":
    main()
