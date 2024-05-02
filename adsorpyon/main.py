
# Local libraries
import input_reader
import interpreter


INPUT_FILE_NAME = "config.in"


def main():
    input_dict = input_reader.create_input_dictionary(INPUT_FILE_NAME)
    properties_dict = input_reader.create_properties_dictionary(input_dict[0]["ADSORBATE_DATA_FILE"])
    data_dict = {}

    interpreter.read_data(
        source_dictionary=data_dict,
        properties_dictionary=properties_dict,
        input_dictionary=input_dict)

    if input_dict[0]["PLOT_DATA"] == "yes":
        interpreter.plot_data(
            source_dictionary=data_dict,
            input_dictionary=input_dict,
            properties_dictionary=properties_dict,
            plot_format=input_dict[0]["DATA_TYPES"],
            save=input_dict[0]["SAVE_DATA_PLOT"],
            from_input=True)

    if input_dict[0]["COMPUTE_CHARACTERISTIC_CURVE"] == "yes":

        interpreter.compute_characteristic(
            source_dictionary=data_dict,
            input_dictionary=input_dict,
            properties_dictionary=properties_dict)

        if input_dict[0]["COMPUTE_SATURATION_PRESSURE_CURVE"] == "yes":
            interpreter.compute_saturation_pressure_curve(
                input_dictionary=input_dict,
                properties_dictionary=properties_dict)

        if input_dict[0]["COMPUTE_DENSITY_CURVE"] == "yes":
            interpreter.compute_density_curve(
                input_dictionary=input_dict,
                properties_dictionary=properties_dict)

        if input_dict[0]["SAVE_CHARACTERISTIC_CURVE_DATA"] == "yes":
            interpreter.write_data(
                source_dictionary=data_dict,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                write_format="characteristic")

        if input_dict[0]["PLOT_CHARACTERISTIC_CURVE"] == "yes":
            interpreter.plot_data(
                source_dictionary=data_dict,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                plot_format="characteristic",
                save=input_dict[0]["SAVE_CHARACTERISTIC_CURVE_PLOT"],
                from_input=False)

    if input_dict[0]["COMPUTE_ENTHALPY"] == "yes":
        enthalpy = interpreter.compute_adsorption_enthalpy(
            data_dictionary=data_dict,
            input_dictionary=input_dict,
            properties_dictionary=properties_dict)

        if input_dict[0]["PLOT_ENTHALPY"] == "yes":
            interpreter.plot_data(
                source_dictionary=enthalpy,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                plot_format="enthalpy",
                save=input_dict[0]["SAVE_ENTHALPY_PLOT"],
                from_input=False)

        if input_dict[0]["SAVE_ENTHALPY_DATA"] == "yes":
            interpreter.write_data(
                source_dictionary=enthalpy,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                write_format="enthalpy")

    if input_dict[0]["PREDICT_ISOTHERMS"] == "yes":
        predicted_isotherms = interpreter.predict_data(
            data_dictionary=data_dict,
            input_dictionary=input_dict,
            prediction_type="isotherm",
            properties_dictionary=properties_dict)

        if input_dict[0]["PLOT_PREDICTED_ISOTHERMS"] == "yes":
            interpreter.plot_data(
                source_dictionary=predicted_isotherms,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                plot_format="isotherm",
                save=input_dict[0]["SAVE_PREDICTED_ISOTHERMS_PLOT"],
                from_input=False)

        if input_dict[0]["SAVE_PREDICTED_ISOTHERMS_DATA"] == "yes":
            interpreter.write_data(
                source_dictionary=predicted_isotherms,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                write_format="isotherm")

    if input_dict[0]["PREDICT_ISOBARS"] == "yes":
        predicted_isobars = interpreter.predict_data(
            data_dictionary=data_dict,
            input_dictionary=input_dict,
            prediction_type="isobar",
            properties_dictionary=properties_dict)

        if input_dict[0]["PLOT_PREDICTED_ISOBARS"] == "yes":
            interpreter.plot_data(
                source_dictionary=predicted_isobars,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                plot_format="isobar",
                save=input_dict[0]["SAVE_PREDICTED_ISOBARS_PLOT"],
                from_input=False)

        if input_dict[0]["SAVE_PREDICTED_ISOBARS_DATA"] == "yes":
            interpreter.write_data(
                source_dictionary=predicted_isobars,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                write_format="isobar")

    if input_dict[0]["PREDICT_ISOSTERES"] == "yes":
        predicted_isosteres = interpreter.predict_data(
            data_dictionary=data_dict,
            input_dictionary=input_dict,
            prediction_type="isostere",
            properties_dictionary=properties_dict)

        if input_dict[0]["PLOT_PREDICTED_ISOSTERES"] == "yes":
            interpreter.plot_data(
                source_dictionary=predicted_isosteres,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                plot_format="isostere",
                save=input_dict[0]["SAVE_PREDICTED_ISOSTERES_PLOT"],
                from_input=False)

        if input_dict[0]["SAVE_PREDICTED_ISOSTERES_DATA"] == "yes":
            interpreter.write_data(
                source_dictionary=predicted_isosteres,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                write_format="isostere")

    if input_dict[0]["SHOW_PLOTS"] == "yes":
        interpreter.show_plots()


if __name__ == "__main__":
    main()
