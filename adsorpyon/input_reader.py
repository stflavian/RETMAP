"""
Scripts needed to read the input file and molecular properties file and convert them to structured dictionaries that can
be interpreted by the main functions.
"""

# Standard libraries
import warnings

DEFAULT_INPUT_DICTIONARY = {
    "DATA_FILES": None,
    "DATA_TYPES": None,
    "ADSORBATE_DATA_FILE": None,
    "PRESSURES": None,
    "TEMPERATURES": None,
    "ADSORBATE": None,
    "ADSORBENT": None,
    "SOURCE": None,
    "TEMPERATURE_UNITS": None,
    "PRESSURE_UNITS": None,
    "LOADING_UNITS": None,
    "POTENTIAL_UNITS": None,
    "VOLUME_UNITS": None,
    "OUTPUT_TEMPERATURE_UNITS": "K",
    "OUTPUT_PRESSURE_UNITS": "MPa",
    "OUTPUT_LOADING_UNITS": "mg/g",
    "OUTPUT_POTENTIAL_UNITS": "kJ/mol",
    "OUTPUT_VOLUME_UNITS": "ml/g",
    "PLOT_DATA": "no",
        "LOGARITHMIC_PLOT": "no",
        "SHOW_DATA_PLOT": "yes",
        "SAVE_DATA_PLOT": "no",
    "COMPUTE_SPECIFIC_ENTHALPY": "no",
        "PLOT_SPECIFIC_ENTHALPY": "yes",
        "SAVE_SPECIFIC_ENTHALPY": "no",
        "SAVE_SPECIFIC_ENTHALPY_PLOT": "no",
    "COMPUTE_CHARACTERISTIC_CURVE": "no",
        "ADSORBATE_SATURATION_PRESSURE": None,
        "SATURATION_PRESSURE_FILE": None,
        "ADSORBATE_DENSITY": None,
        "PLOT_CHARACTERISTIC_CURVE": "yes",
        "SAVE_CHARACTERISTIC_CURVE_DATA": "no",
        "SAVE_CHARACTERISTIC_CURVE_PLOT": "no",
    "EVALUATE_CHARACTERISTIC_CURVE": "no",
        "PLOT_EVALUATION": "yes",
        "SAVE_EVALUATION_DATA": "no",
        "SAVE_EVALUATION_PLOT": "no",
    "PREDICT_ISOTHERMS": "no",
        "PREDICTION_TEMPERATURES": None,
        "PREDICTION_PRESSURE_RANGE": None,
        "REFERENCE_ISOTHERMS": None,
        "PLOT_PREDICTED_ISOTHERMS": "yes",
        "SAVE_PREDICTED_ISOTHERMS_DATA": "no",
        "SAVE_PREDICTED_ISOTHERMS_PLOT": "no",
    "PREDICT_ISOBARS": "no",
        "PREDICTION_PRESSURES": None,
        "PREDICTION_TEMPERATURE_RANGE": None,
        "REFERENCE_ISOBARS": None,
        "PLOT_PREDICTED_ISOBARS": "yes",
        "SAVE_PREDICTED_ISOBARS_DATA": "no",
        "SAVE_PREDICTED_ISOBARS_PLOT": "no",
    "PREDICT_ISOSURFACE": "no",
        "REFERENCE_ISOSURFACES": None,
        "REFERENCE_PRESSURE_UNITS": None,
        "REFERENCE_TEMPERATURE_UNITS": None,
        "REFERENCE_LOADING_UNITS": None,
        "PLOT_PREDICTED_ISOSURFACES": "yes",
        "SAVE_PREDICTED_ISOSURFACES_DATA": "no",
        "SAVE_PREDICTED_ISOSURFACES_PLOT": "no",
    "SHOW_PLOTS": "no",
    "DATA_NAMES": None,
    "WRITE_RESULTS": "yes",
    "SHOW_ISOTHERM": "no",
    "LOGARITHMIC_ISOTHERM": "no",
    "SAVE_ISOTHERM": "yes",
    "SHOW_SPECIFIC_ENTHALPY": "no",
    "SHOW_CHARACTERISTIC_CURVE": "no",
    "SAVE_CHARACTERISTIC_CURVE": "yes",
    "SAVE_EVALUATION": "yes",
    "TEMPERATURE_REFERENCE_ISOTHERM": None,
    "SAVE_PREDICTIONS": "no",
    "COLUMN_1_UNITS": None,
    "COLUMN_2_UNITS": None
}

DEFAULT_PROPERTIES_DICTIONARY = {
    "NAME": None,
    "MOLECULAR_MASS": None,
    "PRESSURE_CRITICAL": None,
    "TEMPERATURE_CRITICAL": None,
    "TEMPERATURE_BOILING": None,
    "TEMPERATURE_TRIPLE_POINT": None,
    "DENSITY_BOILING": None,
    "ACENTRIC_FACTOR": None,
    "THERMAL_EXPANSION_COEFFICIENT": None,
    "AMANKWAH_EXPONENT": None,
    "PRSV_KAPPA1": None,
    "PRSV_KAPPA2": None,
    "PRSV_KAPPA3": None,
    "__SOURCE__": None
}


def create_input_dictionary(path: str) -> dict:
    """
    Convert the input file to structured dictionary that can be read by the rest of the application.

    Open the input file and parse each line separately. If the first word on a line is part of the recognised key-words
    from DEFAULT_DICTIONARY then assign the following words as arguments, otherwise ignore the line. If there are
    multiple data files and only one unique value assigned, use that value as argument for all files.
    :param path: The path of the input file.
    :return: A structured dictionary.
    """
    input_dictionary = {}

    with open(path, "rt") as input_file:
        f = input_file.read()
        lines = f.splitlines()
        for line in lines:
            if line:
                line_input = line.split()
                key_word = line_input[0]
                line_input.pop(0)

                if key_word == "DATA_FILES":
                    for index, _ in enumerate(line_input):
                        input_dictionary[index] = DEFAULT_INPUT_DICTIONARY.copy()

                for index in input_dictionary.keys():
                    if key_word in DEFAULT_INPUT_DICTIONARY and len(line_input) == len(input_dictionary):
                        try:
                            converted_type = float(line_input[index])
                        except ValueError:
                            input_dictionary[index][key_word] = line_input[index]
                        else:
                            input_dictionary[index][key_word] = converted_type
                    elif key_word in DEFAULT_INPUT_DICTIONARY and len(line_input) == 1:
                        try:
                            converted_type = float(line_input[0])
                        except ValueError:
                            input_dictionary[index][key_word] = line_input[0]
                        else:
                            input_dictionary[index][key_word] = converted_type
    return input_dictionary


def create_properties_dictionary(path: str) -> dict:
    """
    Convert the properties file to structured dictionary that can be read by the rest of the application.

    Open the properties file and parse each line separately. If the first word on a line is part of the recognised
    key-words from DEFAULT_DICTIONARY then assign the following words as arguments, otherwise raise error. If there are
    key-words from the dictionary that have no value assigned, raise a warning and continue.
    :param path: The path of the properties file.
    :return: A structured dictionary.
    """
    properties_dictionary = DEFAULT_PROPERTIES_DICTIONARY.copy()

    with open(path, "rt") as properties_file:
        f = properties_file.read()
        lines = f.splitlines()
        for line in lines:
            if line:
                line_input = line.split()
                key_word = line_input[0]
                line_input.pop(0)

                if key_word in properties_dictionary:
                    try:
                        converted_type = float(line_input[0])
                    except ValueError:
                        properties_dictionary[key_word] = line_input[0]
                    else:
                        properties_dictionary[key_word] = converted_type
                else:
                    raise ValueError(f"{key_word} in {path} is not a recognised tag for a properties file. Remove the "
                                     f"tag or check for spelling errors!")

    for key_word in properties_dictionary.keys():
        if not properties_dictionary[key_word]:
            warnings.warn(f"{key_word} was not found in the properties file at {path}, which may be a requirement for "
                          f"some methods. In case of errors please add the tag to the properties file!")

    return properties_dictionary
