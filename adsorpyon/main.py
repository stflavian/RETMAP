"""Main

This represents the main script.

Usage: run
"""

import utils
import interpreter
import input_reader
import physics

from numpy import genfromtxt, nan


def main():

    # Read the config and adsorbate data files
    config = input_reader.create_input_dictionary("config.in")
    adsorbate_data = input_reader.create_properties_dictionary(config[0]["ADSORBATE_DATA_FILE"])

    # Create the data array and assign values
    data = {}
    for index in config:
        data[index] = {}

        # Assign data from the config file
        data[index]["data_file"] = config[index]["DATA_FILES"]
        data[index]["name"] = config[index]["DATA_NAMES"]
        data[index]["temperature"] = config[index]["TEMPERATURES"]
        data[index]["data_type"] = config[index]["DATA_TYPES"]
        data[index]["saturation_pressure_function"] = config[index]["ADSORBATE_SATURATION_PRESSURE"]
        data[index]["saturation_pressure_file"] = config[index]["SATURATION_PRESSURE_FILE"]
        data[index]["density_function"] = config[index]["ADSORBATE_DENSITY"]
        data[index]["adsorbate_data_file"] = config[index]["ADSORBATE_DATA_FILE"]
        data[index]["adsorbent"] = config[index]["ADSORBENT"]
        data[index]["adsorbate"] = config[index]["ADSORBATE"]

        # Calculate the saturation pressure and density based on the desired functions
        data[index]["saturation_pressure"] = interpreter.compute_saturation_pressure_from_method(
                                                    method=data[index]["saturation_pressure_function"],
                                                    temperature=data[index]["temperature"],
                                                    properties_dictionary=adsorbate_data,
                                                    saturation_pressure_file=data[index]["saturation_pressure_file"])

        data[index]["density"] = interpreter.compute_density_from_method(
                                                    method=data[index]["density_function"],
                                                    temperature=data[index]["temperature"],
                                                    properties_dictionary=adsorbate_data)

        # Read the input data and calculate the isotherm or characteristic curve depending on data type
        file_data = genfromtxt(fname=data[index]["data_file"], filling_values=nan)
        if data[index]["data_type"] == "isotherm":
            data[index]["pressure"] = file_data[:, 0] * utils.convert_input(
                                                                        unit=config[index]["COLUMN_1_UNITS"],
                                                                        molecular_mass=adsorbate_data["MOLECULAR_MASS"])

            data[index]["adsorbed_amount"] = file_data[:, 1] * utils.convert_input(
                                                                        unit=config[index]["COLUMN_2_UNITS"],
                                                                        molecular_mass=adsorbate_data["MOLECULAR_MASS"])

            data[index]["adsorption_potential"] = physics.get_adsorption_potential(
                                                                temperature=data[index]["temperature"],
                                                                saturation_pressure=data[index]["saturation_pressure"],
                                                                pressure=data[index]["pressure"])

            data[index]["adsorption_volume"] = physics.get_adsorption_volume(
                                                                        adsorbed_amount=data[index]["adsorbed_amount"],
                                                                        adsorbate_density=data[index]["density"])

            data[index]["adsorption_enthalpy"] = physics.get_adsorption_enthalpy(
                                        enthalpy_vaporization=15,
                                        temperature=data[index]["temperature"],
                                        adsorption_volume=data[index]["adsorption_volume"],
                                        adsorption_potential=data[index]["adsorption_potential"],
                                        thermal_expansion_coefficient=adsorbate_data["THERMAL_EXPANSION_COEFFICIENT"])

        if data[index]["data_type"] == "characteristic":
            data[index]["adsorption_potential"] = file_data[:, 0] * utils.convert_input(
                                                                        unit=config[index]["COLUMN_1_UNITS"],
                                                                        molecular_mass=adsorbate_data["MOLECULAR_MASS"])

            data[index]["adsorption_volume"] = file_data[:, 1] * utils.convert_input(
                                                                        unit=config[index]["COLUMN_2_UNITS"],
                                                                        molecular_mass=adsorbate_data["MOLECULAR_MASS"])

            data[index]["pressure"] = physics.get_pressure(
                                                            adsorption_potential=data[index]["adsorption_potential"],
                                                            saturation_pressure=data[index]["saturation_pressure"],
                                                            temperature=data[index]["temperature"])

            data[index]["adsorbed_amount"] = physics.get_adsorbed_amount(
                                                                    adsorption_volume=data[index]["adsorption_volume"],
                                                                    adsorbate_density=data[index]["density"])

            data[index]["adsorption_enthalpy"] = physics.get_adsorption_enthalpy(
                                        enthalpy_vaporization=15,
                                        temperature=data[index]["temperature"],
                                        adsorption_volume=data[index]["adsorption_volume"],
                                        adsorption_potential=data[index]["adsorption_potential"],
                                        thermal_expansion_coefficient=adsorbate_data["thermal_expansion_coefficient"])
            
        if config[0]["WRITE_RESULTS"].lower() == "yes":
            utils.write_data(data=data, index=index)

    plot_presence = False
    if config[0]["SHOW_ISOTHERM"].lower() == "yes":
        plot_presence = True
        utils.plot_isotherm(data=data, logarithmic=config[0]["LOGARITHMIC_ISOTHERM"], save=config[0]["SAVE_ISOTHERM"])

    if config[0]["SHOW_SPECIFIC_ENTHALPY"].lower() == "yes":
        plot_presence = True
        utils.plot_enthalpy(data=data, save=config[0]["SAVE_SPECIFIC_ENTHALPY"])

    if config[0]["SHOW_CHARACTERISTIC_CURVE"].lower() == "yes":
        plot_presence = True
        utils.plot_characteristic_curve(data=data, save=config[0]["SAVE_CHARACTERISTIC_CURVE"])

    if config[0]["EVALUATE_CHARACTERISTIC_CURVE"].lower() == "yes":
        plot_presence = True
        utils.evaluate_characteristic_curve(data=data, save=config[0]["SAVE_EVALUATION"],
                                            temperature_reference_isotherm=config[0]["TEMPERATURE_REFERENCE_ISOTHERM"])

    if config[0]["PREDICT_ISOTHERMS"].lower() == "yes":
        plot_presence = True
        utils.predict_isotherms(data=data, logarithmic=config[0]["LOGARITHMIC_ISOTHERM"],
                                save=config[0]["SAVE_PREDICTIONS"],
                                temperature_reference_isotherm=config[0]["TEMPERATURE_REFERENCE_ISOTHERM"])

    if plot_presence is True:
        utils.show_plots()


if __name__ == "__main__":
    main()
