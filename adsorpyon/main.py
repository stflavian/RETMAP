import json
from adsorpyon import literature_functions as lf, transformer as trf
from numpy import genfromtxt, nan


def saturation_pressure(data: dict, index: int, adsorbate_data: dict) -> float:
    if data[index]["saturation_pressure_function"] == "dubinin":
        return lf.pressure_dubinin(temperature=data[index]["temperature"],
                                   temperature_critical=adsorbate_data["temperature_critical"],
                                   pressure_critical=adsorbate_data["pressure_critical"])
    elif data[index]["saturation_pressure_function"] == "amankwah":
        return lf.pressure_amankwah(temperature=data[index]["temperature"],
                                    temperature_critical=adsorbate_data["temperature_critical"],
                                    pressure_critical=adsorbate_data["pressure_critical"],
                                    k=1.42)
    elif data[index]["saturation_pressure_function"] == "extrapolation":
        return lf.pressure_extrapolation(temperature=data[index]["temperature"],
                                         file=data[index]["saturation_pressure_file"])
    elif data[index]["saturation_pressure_function"] == "polynomial":
        return lf.pressure_polynomial(temperature=data[index]["temperature"])
    elif data[index]["saturation_pressure_function"] == "peng-robinson":
        return lf.pressure_pengrobinson(temperature=data[index]["temperature"],
                                        temperature_critical=adsorbate_data["temperature_critical"],
                                        pressure_critical=adsorbate_data["pressure_critical"],
                                        pressure_guess=1,
                                        acentric_factor=adsorbate_data["acentric_factor"])
    else:
        raise ValueError(f"Input string {data[index]['saturation_pressure_function']} does not coresspond to a valid"
                         f"method of determining the saturation pressure")


def density(data: dict, index: int, adsorbate_data: dict) -> float:
    if data[index]["density_function"] == "ozawa":
        return lf.density_ozawa(temperature=data[index]["temperature"],
                                temperature_boiling=adsorbate_data["temperature_boiling"],
                                density_boiling=adsorbate_data["density_boiling"])
    elif data[index]["density_function"] == "empirical":
        return lf.density_empirical(pressure_critical=adsorbate_data["pressure_critical"],
                                    temperature_critical=adsorbate_data["temperature_critical"],
                                    molecular_mass=adsorbate_data["molecular_mass"])
    elif data[index]["density_function"] == "hauer":
        return lf.density_hauer(temperature=data[index]["temperature"],
                                temperature_reference=adsorbate_data["temperature_reference"],
                                density_boiling=adsorbate_data["density_boiling"])
    else:
        raise ValueError(f"Input string {data[index]['density_function']} does not coresspond to a valid"
                         f"method of determining the density")


def main():

    # Read the config and adsorbate data files
    config = json.load(open("../config.json"))
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
        data[index]["saturation_pressure"] = saturation_pressure(data=data, index=index, adsorbate_data=adsorbate_data)
        data[index]["density"] = density(data=data, index=index, adsorbate_data=adsorbate_data)

        # Read the input data and calculate the isotherm or characteristic curve depending on data type
        file_data = genfromtxt(fname=file, filling_values=nan)
        if data[index]["data_type"] == "isotherm":
            data[index]["pressure"] = file_data[:, 0] * trf.convert(unit=config["column1_units"][index],
                                                                    adsorbate_data=adsorbate_data)  # [MPa]

            data[index]["adsorbed_amount"] = file_data[:, 1] * trf.convert(unit=config["column2_units"][index],
                                                                           adsorbate_data=adsorbate_data)  # [mg/g]

            data[index]["adsorption_potential"] = trf.get_adsorption_potential(temperature=data[index]["temperature"],
                                                                               saturation_pressure=data[index]["saturation_pressure"],
                                                                               pressure=data[index]["pressure"])  # [kJ/mol]

            data[index]["adsorption_volume"] = trf.get_adsorption_volume(adsorbed_amount=data[index]["adsorbed_amount"],
                                                                         adsorbate_density=data[index]["density"])  # [ml/g]
        if data[index]["data_type"] == "characteristic curve":
            data[index]["adsorption_potential"] = file_data[:, 0] * trf.convert(unit=config["column1_units"][index],
                                                                                adsorbate_data=adsorbate_data)  # [kJ/mol]

            data[index]["adsorption_volume"] = file_data[:, 1] * trf.convert(unit=config["column2_units"][index],
                                                                             adsorbate_data=adsorbate_data)  # [ml/g]

            data[index]["pressure"] = trf.get_pressure(adsorption_potential=data[index]["adsorption_potential"],
                                                       saturation_pressure=data[index]["saturation_pressure"],
                                                       temperature=data[index]["temperature"])  # [MPa]

            data[index]["adsorbed_amount"] = trf.get_adsorbed_amount(adsorption_volume=data[index]["adsorption_volume"],
                                                                     adsorbate_density=data[index]["density"])  # [mg/g]
            
        if config["write_results"][index].lower() == "yes":
            trf.write_data(data=data, index=index)

    plot_presence = False
    if config["show_isotherm"].lower() == "yes":
        plot_presence = True
        trf.plot_isotherm(data=data, logarithmic=config["logarithmic_isotherm"], save=config["save_isotherm"])

    if config["show_characteristic_curve"].lower() == "yes":
        plot_presence = True
        trf.plot_characteristic_curve(data=data, save=config["save_characteristic_curve"])

    if plot_presence is True:
        trf.show_plots()


if __name__ == "__main__":
    main()
