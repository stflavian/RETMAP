import matplotlib.pyplot as plt

# Local libraries
import input_reader
import interpreter
import utils
import physics

# Third-party libraries
from numpy import genfromtxt, nan


INPUT_FILE_NAME = "config.in"


def main():
    input_dict = input_reader.create_input_dictionary(INPUT_FILE_NAME)
    properties_dict = input_reader.create_properties_dictionary(input_dict[0]["ADSORBATE_DATA_FILE"])
    data_dict = {}

    for index in input_dict:
        data_dict[index] = {}
        properties_dict = input_reader.create_properties_dictionary(input_dict[index]["ADSORBATE_DATA_FILE"])

        file_data = genfromtxt(fname=input_dict[index]["DATA_FILES"], filling_values=nan)
        if input_dict[index]["DATA_TYPES"] == "isotherm":
            data_dict[index]["temperature"] = input_dict[index]["TEMPERATURES"]
            data_dict[index]["pressure"] = file_data[:, 0] * utils.convert_input(
                                                                    unit=input_dict[index]["PRESSURE_UNITS"],
                                                                    molecular_mass=properties_dict["MOLECULAR_MASS"])

            data_dict[index]["loading"] = file_data[:, 1] * utils.convert_input(
                                                                    unit=input_dict[index]["LOADING_UNITS"],
                                                                    molecular_mass=properties_dict["MOLECULAR_MASS"])
        elif input_dict[index]["DATA_TYPES"] == "isobar":
            data_dict[index]["pressure"] = input_dict[index]["PRESSURES"]
            data_dict[index]["temperature"] = file_data[:, 0] * utils.convert_input(
                                                                    unit=input_dict[index]["TEMPERATURE_UNITS"],
                                                                    molecular_mass=properties_dict["MOLECULAR_MASS"])

            data_dict[index]["loading"] = file_data[:, 1] * utils.convert_input(
                                                                    unit=input_dict[index]["LOADING_UNITS"],
                                                                    molecular_mass=properties_dict["MOLECULAR_MASS"])
        elif input_dict[index]["DATA_TYPES"] == "characteristic":
            data_dict[index]["potential"] = file_data[:, 0] * utils.convert_input(
                                                                    unit=input_dict[index]["POTENTIAL_UNITS"],
                                                                    molecular_mass=properties_dict["MOLECULAR_MASS"])

            data_dict[index]["volume"] = file_data[:, 1] * utils.convert_input(
                                                                    unit=input_dict[index]["VOLUME_UNITS"],
                                                                    molecular_mass=properties_dict["MOLECULAR_MASS"])
        else:
            raise ValueError(f"{input_dict[index]['DATA_TYPES']} is not a recognized argument!")

    if input_dict[0]["PLOT_DATA"] == "yes":
        interpreter.plot_data(source_dictionary=data_dict, input_dictionary=input_dict,
                              plot_format=input_dict[0]["DATA_TYPES"], save=input_dict[0]["SAVE_DATA_PLOT"])

    if input_dict[0]["COMPUTE_SPECIFIC_ENTHALPY"] is "yes":
        print("This function is not implemented yet!")

    if input_dict[0]["COMPUTE_CHARACTERISTIC_CURVE"] == "yes":

        for index in data_dict:
            if input_dict[index]["DATA_TYPES"] == "isotherm":
                data_dict[index]["saturation_pressure"] = interpreter.compute_saturation_pressure_from_method(
                                                method=input_dict[index]["ADSORBATE_SATURATION_PRESSURE"],
                                                temperature=input_dict[index]["TEMPERATURES"],
                                                properties_dictionary=properties_dict,
                                                saturation_pressure_file=input_dict[index]["SATURATION_PRESSURE_FILE"])

                data_dict[index]["density"] = interpreter.compute_density_from_method(
                                                                        method=input_dict[index]["ADSORBATE_DENSITY"],
                                                                        temperature=input_dict[index]["TEMPERATURES"],
                                                                        properties_dictionary=properties_dict)

            elif input_dict[index]["DATA_TYPES"] == "isobar":
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

            elif input_dict[index]["DATA_TYPE"] == "characteristic":
                continue

            data_dict[index]["potential"] = physics.get_adsorption_potential(
                                                            temperature=input_dict[index]["TEMPERATURES"],
                                                            saturation_pressure=data_dict[index]["saturation_pressure"],
                                                            pressure=data_dict[index]["pressure"])

            data_dict[index]["volume"] = physics.get_adsorption_volume(
                                                            adsorbed_amount=data_dict[index]["loading"],
                                                            adsorbate_density=data_dict[index]["density"])

        if input_dict[0]["SAVE_CHARACTERISTIC_CURVE_DATA"] == "yes":
            interpreter.write_data(source_dictionary=data_dict, input_dictionary=input_dict,
                                   write_format="characteristic")

        if input_dict[0]["PLOT_CHARACTERISTIC_CURVE"] == "yes":
            interpreter.plot_data(source_dictionary=data_dict, input_dictionary=input_dict,
                                  plot_format="characteristic", save=input_dict[0]["SAVE_CHARACTERISTIC_CURVE_PLOT"])

    if input_dict[0]["EVALUATE_CHARACTERISTIC_CURVE"] == "yes":
        print("Evaluate characteristic curve")

    if input_dict[0]["PREDICT_ISOTHERMS"] == "yes":
        predicted_isotherms = interpreter.predict_data(data_dictionary=data_dict,
                                                       input_dictionary=input_dict,
                                                       prediction_type="isotherm",
                                                       properties_dictionary=properties_dict)

        if input_dict[0]["PLOT_PREDICTED_ISOTHERMS"] == "yes":
            interpreter.plot_data(source_dictionary=predicted_isotherms, input_dictionary=input_dict,
                                  plot_format="isotherm", save=input_dict[0]["SAVE_PREDICTED_ISOTHERMS_PLOT"])

        if input_dict[0]["SAVE_PREDICTED_ISOTHERMS_DATA"] == "yes":
            interpreter.write_data(source_dictionary=predicted_isotherms, input_dictionary=input_dict,
                                   write_format="isotherm")

    if input_dict[0]["PREDICT_ISOBARS"] == "yes":
        predicted_isobars = interpreter.predict_data(data_dictionary=data_dict,
                                                     input_dictionary=input_dict,
                                                     prediction_type="isobar",
                                                     properties_dictionary=properties_dict)

        if input_dict[0]["PLOT_PREDICTED_ISOBARS"] == "yes":
            interpreter.plot_data(source_dictionary=predicted_isobars, input_dictionary=input_dict,
                                  plot_format="isobar", save=input_dict[0]["SAVE_PREDICTED_ISOBARS_PLOT"])

        if input_dict[0]["SAVE_PREDICTED_ISOBARS_DATA"] == "yes":
            interpreter.write_data(source_dictionary=predicted_isobars, input_dictionary=input_dict,
                                   write_format="isobar")

    if input_dict[0]["PREDICT_ISOSURFACE"] == "yes":
        print("Predict isosurface")

    interpreter.show_plots()


if __name__ == "__main__":
    main()
