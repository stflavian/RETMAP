
# Local libraries
import input_reader
import interpreter
import utils
import physics

# Third-party libraries
from numpy import genfromtxt, nan


INPUT_FILE_NAME = "input.config"


def main():
    input_dict = input_reader.create_input_dictionary(INPUT_FILE_NAME)
    properties_dict = input_reader.create_properties_dictionary(input_dict[0]["ADSORBATE_DATA_FILE"])
    data_dict = {}

    for index in input_dict:
        data_dict[index] = {}
        properties_dict = input_reader.create_properties_dictionary(input_dict[index]["ADSORBATE_DATA_FILE"])

        file_data = genfromtxt(fname=input_dict[index]["DATA_FILES"], filling_values=nan)
        if input_dict[index]["DATA_TYPES"] is "isotherm":
            data_dict[index]["pressure"] = file_data[:, 0] * utils.convert_input(
                                                                    unit=input_dict[index]["PRESSURE_UNITS"],
                                                                    molecular_mass=properties_dict["MOLECULAR_MASS"])

            data_dict[index]["loading"] = file_data[:, 1] * utils.convert_input(
                                                                    unit=input_dict[index]["LOADING_UNITS"],
                                                                    molecular_mass=properties_dict["MOLECULAR_MASS"])
        elif input_dict[index]["DATA_TYPES"] is "isobar":
            data_dict[index]["temperature"] = file_data[:, 0] * utils.convert_input(
                                                                    unit=input_dict[index]["TEMPERATURE_UNITS"],
                                                                    molecular_mass=properties_dict["MOLECULAR_MASS"])

            data_dict[index]["loading"] = file_data[:, 1] * utils.convert_input(
                                                                    unit=input_dict[index]["LOADING_UNITS"],
                                                                    molecular_mass=properties_dict["MOLECULAR_MASS"])
        elif input_dict[index]["DATA_TYPES"] is "characteristic":
            data_dict[index]["potential"] = file_data[:, 0] * utils.convert_input(
                                                                    unit=input_dict[index]["POTENTIAL_UNITS"],
                                                                    molecular_mass=properties_dict["MOLECULAR_MASS"])

            data_dict[index]["volume"] = file_data[:, 1] * utils.convert_input(
                                                                    unit=input_dict[index]["VOLUME_UNITS"],
                                                                    molecular_mass=properties_dict["MOLECULAR_MASS"])
        else:
            raise ValueError(f"{input_dict[index]['DATA_TYPES']} is not a recognized argument!")

    if input_dict[0]["PLOT_DATA"] is "yes":
        print("Plot data")

    if input_dict[0]["COMPUTE_SPECIFIC_ENTHALPY"] is "yes":
        print("This function is not implemented yet!")

    if input_dict[0]["COMPUTE_CHARACTERISTIC_CURVE"] is "yes":

        for index in data_dict:
            if input_dict[index]["DATA_TYPES"] is "isotherm":
                data_dict[index]["saturation_pressure"] = interpreter.compute_saturation_pressure_from_method(
                                                method=input_dict[index]["ADSORBATE_SATURATION_PRESSURE"],
                                                temperature=input_dict[index]["TEMPERATURES"],
                                                properties_dictionary=properties_dict,
                                                saturation_pressure_file=input_dict[index]["SATURATION_PRESSURE_FILE"])

                data_dict[index]["density"] = interpreter.compute_density_from_method(
                                                                        method=input_dict[index]["ADSORBATE_DENSITY"],
                                                                        temperature=input_dict[index]["TEMPERATURES"],
                                                                        properties_dictionary=properties_dict)

            elif input_dict[index]["DATA_TYPES"] is "isobar":
                saturation_pressure_array = []
                density_array = []
                for temperature in data_dict[index]["temperature"]:
                    saturation_pressure = interpreter.compute_saturation_pressure_from_method(
                                                method=input_dict[index]["ADSORBATE_SATURATION_PRESSURE"],
                                                temperature=temperature,
                                                properties_dictionary=properties_dict,
                                                saturation_pressure_file=input_dict[index]["SATURATION_PRESSURE_FILE"])
                    saturation_pressure_array.append(saturation_pressure)

                    density = interpreter.compute_density_from_method(
                                                method=input_dict[index]["ADSORBATE_DENSITY"],
                                                temperature=temperature,
                                                properties_dictionary=properties_dict)
                    density_array.append(density)

                data_dict[index]["saturation_pressure"] = saturation_pressure_array
                data_dict[index]["density"] = density_array

            elif input_dict[index]["DATA_TYPE"] is "characteristic":
                continue

            data_dict[index]["potential"] = physics.get_adsorption_potential(
                                                            temperature=input_dict[index]["TEMPERATURES"],
                                                            saturation_pressure=data_dict[index]["saturation_pressure"],
                                                            pressure=data_dict[index]["pressure"])

            data_dict[index]["volume"] = physics.get_adsorption_volume(
                                                            adsorbed_amount=data_dict[index]["loading"],
                                                            adsorbate_density=data_dict[index]["density"])

            if input_dict[0]["SAVE_CHARACTERISTIC_CURVE_DATA"] is "yes":


            if input_dict[0]["PLOT_CHARACTERISTIC_CURVE"] is "yes":


            if input_dict[0]["SAVE_CHARACTERISTIC_CURVE_PLOT"] is "yes":




    if input_dict[0]["EVALUATE_CHARACTERISTIC_CURVE"] is "yes":
        print("Evaluate characteristic curve")

    if input_dict[0]["PREDICT_ISOTHERMS"] is "yes":
        print("Predict isotherms")

    if input_dict[0]["PREDICT_ISOBARS"] is "yes":
        print("Predict isobars")

    if input_dict[0]["PREDICT_ISOSURFACE"] is "yes":
        print("Predict isosurface")



    if config["write_results"][index].lower() == "yes":
        utils.write_data(data=data, index=index)

    plot_presence = False
    if config["show_isotherm"].lower() == "yes":
        plot_presence = True
        utils.plot_isotherm(data=data, logarithmic=config["logarithmic_isotherm"], save=config["save_isotherm"])

    if config["show_specific_enthalpy"].lower() == "yes":
        plot_presence = True
        utils.plot_enthalpy(data=data, save=config["save_specific_enthalpy"])

    if config["show_characteristic_curve"].lower() == "yes":
        plot_presence = True
        utils.plot_characteristic_curve(data=data, save=config["save_characteristic_curve"])

    if config["evaluate_characteristic_curve"].lower() == "yes":
        plot_presence = True
        utils.evaluate_characteristic_curve(data=data, save=config["save_evaluation"],
                                            temperature_reference_isotherm=config["temperature_reference_isotherm"])

    if config["predict_isotherms"].lower() == "yes":
        plot_presence = True
        utils.predict_isotherms(data=data, logarithmic=config["logarithmic_isotherm"], save=config["save_predictions"],
                                temperature_reference_isotherm=config["temperature_reference_isotherm"])

    if input_dict[0]["SHOW_PLOTS"] is True:
        utils.show_plots()


if __name__ == "__main__":
    main()
