#!/home/flavian/.venv/envScience/bin/python

# Local libraries
from src import input_reader
from src import interpreter


INPUT_FILE_NAME = "config.in"
output_file = open("cappa.out", "w+")


def main():
    
    output_file.write("\n")
    input_dict = input_reader.create_input_dictionary(INPUT_FILE_NAME)
    message = " Input file data "
    message = message.center(80, "=")
    output_file.write(f"{message}\n")
    for entry in input_dict:
        output_file.write("\n")
        message = f" Entry {entry} "
        message = message.center(80, "-")
        output_file.write(f"{message}\n")
        for key in input_dict[entry]:
            key_message = f" {key} "
            key_message = key_message.ljust(40)
            value_message = f" {input_dict[entry][key]} "
            value_message = value_message.rjust(40)
            output_file.write(f"{key_message} {value_message}\n")

    output_file.write("\n")
    properties_dict = input_reader.create_properties_dictionary(input_dict[0]['ADSORBATE_DATA_FILE'])
    message = " Molecular properties "
    message = message.center(80, "=")
    output_file.write(f"{message}\n")
    for key in properties_dict:
        key_message = f" {key} "
        key_message = key_message.ljust(40)
        value_message = f" {properties_dict[key]} "
        value_message = value_message.rjust(40)
        output_file.write(f"{key_message} {value_message}\n")

    data_dict = {}

    interpreter.read_data(
        source_dictionary=data_dict,
        properties_dictionary=properties_dict,
        input_dictionary=input_dict)

    output_file.write("\n")
    message = " Input data "
    message = message.center(80, "=")
    output_file.write(f"{message}\n")
    for entry in data_dict:
        output_file.write("\n")
        message = f" Entry {entry} "
        message = message.center(80, "-")
        output_file.write(f"{message}\n")
        for key in data_dict[entry]:
            output_file.write(f" {key} \n")
            output_file.write(f" {data_dict[entry][key]} \n")
            output_file.write("\n")

    if input_dict[0]['PLOT_DATA'].lower() == "yes":
        interpreter.plot_data(
            source_dictionary=data_dict,
            input_dictionary=input_dict,
            properties_dictionary=properties_dict,
            plot_format=input_dict[0]['DATA_TYPES'],
            save=input_dict[0]['SAVE_DATA_PLOT'],
            from_input=True)

    if input_dict[0]['COMPUTE_CHARACTERISTIC_CURVE'].lower() == "yes":

        interpreter.compute_characteristic(
            source_dictionary=data_dict,
            input_dictionary=input_dict,
            properties_dictionary=properties_dict)

        output_file.write("\n")
        message = " Characteristic curve calculations "
        message = message.center(80, "=")
        output_file.write(f"{message}\n")
        for entry in data_dict:
            output_file.write("\n")
            message = f" Entry {entry} "
            message = message.center(80, "-")
            output_file.write(f"{message}\n")
            for key in ["potential", "volume"]:
                output_file.write(f" {key} \n")
                output_file.write(f" {data_dict[entry][key]} \n")
                output_file.write("\n")

        if input_dict[0]['COMPUTE_SATURATION_PRESSURE_CURVE'].lower() == "yes":
            interpreter.compute_saturation_pressure_curve(
                input_dictionary=input_dict,
                properties_dictionary=properties_dict)

        if input_dict[0]['COMPUTE_DENSITY_CURVE'].lower() == "yes":
            interpreter.compute_density_curve(
                input_dictionary=input_dict,
                properties_dictionary=properties_dict)

        if input_dict[0]['SAVE_CHARACTERISTIC_CURVE_DATA'].lower() == "yes":
            interpreter.write_data(
                source_dictionary=data_dict,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                write_format="characteristic")

        if input_dict[0]['PLOT_CHARACTERISTIC_CURVE'].lower() == "yes":
            interpreter.plot_data(
                source_dictionary=data_dict,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                plot_format="characteristic",
                save=input_dict[0]['SAVE_CHARACTERISTIC_CURVE_PLOT'],
                from_input=False)

    if input_dict[0]['COMPUTE_ENTHALPY'].lower() == "yes":
        enthalpy = interpreter.compute_adsorption_enthalpy(
            data_dictionary=data_dict,
            input_dictionary=input_dict,
            properties_dictionary=properties_dict)

        output_file.write("\n")
        message = " Enthalpy of adsorption calculations "
        message = message.center(80, "=")
        output_file.write(f"{message}\n")
        for key in enthalpy:
            output_file.write(f" {key} \n")
            output_file.write(f" {enthalpy[key]} \n")

        if input_dict[0]['PLOT_ENTHALPY'].lower() == "yes":
            interpreter.plot_data(
                source_dictionary=enthalpy,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                plot_format="enthalpy",
                save=input_dict[0]['SAVE_ENTHALPY_PLOT'],
                from_input=False)

        if input_dict[0]['SAVE_ENTHALPY_DATA'].lower() == "yes":
            interpreter.write_data(
                source_dictionary=enthalpy,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                write_format="enthalpy")

    if input_dict[0]['PREDICT_ISOTHERMS'].lower() == "yes":
        predicted_isotherms = interpreter.predict_data(
            data_dictionary=data_dict,
            input_dictionary=input_dict,
            prediction_type="isotherm",
            properties_dictionary=properties_dict)

        output_file.write("\n")
        message = " Isotherm predictions "
        message = message.center(80, "=")
        output_file.write(f"{message}\n")
        for entry in predicted_isotherms:
            output_file.write("\n")
            message = f" Entry {entry} "
            message = message.center(80, "-")
            output_file.write(f"{message}\n")
            for key in predicted_isotherms[entry]:
                output_file.write(f" {key} \n")
                output_file.write(f" {predicted_isotherms[entry][key]} \n")
                output_file.write("\n")

        if input_dict[0]['PLOT_PREDICTED_ISOTHERMS'].lower() == "yes":
            interpreter.plot_data(
                source_dictionary=predicted_isotherms,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                plot_format="isotherm",
                save=input_dict[0]['SAVE_PREDICTED_ISOTHERMS_PLOT'],
                from_input=False)

        if input_dict[0]['SAVE_PREDICTED_ISOTHERMS_DATA'].lower() == "yes":
            interpreter.write_data(
                source_dictionary=predicted_isotherms,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                write_format="isotherm")

    if input_dict[0]['PREDICT_ISOBARS'].lower() == "yes":
        predicted_isobars = interpreter.predict_data(
            data_dictionary=data_dict,
            input_dictionary=input_dict,
            prediction_type="isobar",
            properties_dictionary=properties_dict)
        
        output_file.write("\n")
        message = " Isobar predictions "
        message = message.center(80, "=")
        output_file.write(f"{message}\n")
        for entry in predicted_isobars:
            output_file.write("\n")
            message = f" Entry {entry} "
            message = message.center(80, "-")
            output_file.write(f"{message}\n")
            for key in predicted_isobars[entry]:
                output_file.write(f" {key} \n")
                output_file.write(f" {predicted_isobars[entry][key]} \n")
                output_file.write("\n")

        if input_dict[0]['PLOT_PREDICTED_ISOBARS'].lower() == "yes":
            interpreter.plot_data(
                source_dictionary=predicted_isobars,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                plot_format="isobar",
                save=input_dict[0]['SAVE_PREDICTED_ISOBARS_PLOT'],
                from_input=False)

        if input_dict[0]['SAVE_PREDICTED_ISOBARS_DATA'].lower() == "yes":
            interpreter.write_data(
                source_dictionary=predicted_isobars,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                write_format="isobar")

    if input_dict[0]['PREDICT_ISOSTERES'].lower() == "yes":
        predicted_isosteres = interpreter.predict_data(
            data_dictionary=data_dict,
            input_dictionary=input_dict,
            prediction_type="isostere",
            properties_dictionary=properties_dict)
        
        output_file.write("\n")
        message = " Isostere predictions "
        message = message.center(80, "=")
        output_file.write(f"{message}\n")
        for entry in predicted_isosteres:
            output_file.write("\n")
            message = f" Entry {entry} "
            message = message.center(80, "-")
            output_file.write(f"{message}\n")
            for key in predicted_isosteres[entry]:
                output_file.write(f" {key} \n")
                output_file.write(f" {predicted_isosteres[entry][key]} \n")
                output_file.write("\n")

        if input_dict[0]['PLOT_PREDICTED_ISOSTERES'].lower() == "yes":
            interpreter.plot_data(
                source_dictionary=predicted_isosteres,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                plot_format="isostere",
                save=input_dict[0]['SAVE_PREDICTED_ISOSTERES_PLOT'],
                from_input=False)

        if input_dict[0]['SAVE_PREDICTED_ISOSTERES_DATA'].lower() == "yes":
            interpreter.write_data(
                source_dictionary=predicted_isosteres,
                input_dictionary=input_dict,
                properties_dictionary=properties_dict,
                write_format="isostere")

    if input_dict[0]['SHOW_PLOTS'].lower() == "yes":
        interpreter.show_plots()


if __name__ == "__main__":
    main()
