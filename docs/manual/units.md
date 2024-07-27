---
icon: material/ruler
---

# Measurements units

CAPPA works with six main metrics: pressure, temperature, loading, potential, volume,
and density. These metrics are stored and processed using the following units: MPa, K,
mg/g, kJ/mol, ml/g, and kg/m3. As such, all values provided inside the configuration 
file (e.g. `PREDICTION_TEMPERATURES`) should be provided using the aforementioned units.
However, the user can specify inside the configuration file if the data read from the
input files is provided in different units, in which case conversion factors are used 
to reach the default units. Similarly, the user can specify the units used for writing
the results. 

Currently, the following measurement units are supported:


!!! info "`TEMPERATURE_UNITS`, `OUTPUT_TEMPERATURE_UNITS` (default: K, type: String)"
    
    These keywords can be used to specify the measurement unit for temperature for the
    input and output files, respectively. The keyword `TEMPERATURE_UNITS` only applies
    to the data read from isobars and isosteres, and not to the input provided inside the 
    configuration file. Moreover, it accepts an entry for each input data file. The 
    keyword `OUTPUT_TEMPERATURE_UNITS` is used for both writing the results to files
    and plotting the results. However, the `cappa.out` file will still use the standard
    units. Only one output unit can be provided.

    __Options__: _K_, _R_


!!! info "`PRESSURE_UNITS`, `OUTPUT_PRESSURE_UNITS` (default: MPa, type: String)"
    
    These keywords can be used to specify the measurement unit for pressure for the
    input and output files, respectively. The keyword `PRESSURE_UNITS` only applies
    to the data read from isotherms and isosteres, and not to the input provided inside 
    the configuration file. Moreover, it accepts an entry for each input data file. The 
    keyword `OUTPUT_PRESSURE_UNITS` is used for both writing the results to files
    and plotting the results. However, the `cappa.out` file will still use the standard
    units. Only one output unit can be provided.

    __Options__: _MPa_, _kPa_, _Pa_, _bar_, _atm_, _Torr_, _mmHg_


!!! info "`LOADING_UNITS`, `OUTPUT_LOADING_UNITS` (default: mg/g, type: String)" 

    These keywords can be used to specify the measurement unit for loading for the
    input and output files, respectively. The keyword `LOADING_UNITS` only applies
    to the data read from isotherms and isobars, and not to the input provided inside 
    the configuration file. Moreover, it accepts an entry for each input data file. The 
    keyword `OUTPUT_LOADING_UNITS` is used for both writing the results to files
    and plotting the results. However, the `cappa.out` file will still use the standard
    units. Only one output unit can be provided.

    __Options__: mg/g, _g/kg_, _mol/kg_, _mmol/g_, _cm3/g_, _mm3/mg_, _dm3/kg_, _l/kg_, 
    _ml/g_


!!! info "`POTENTIAL_UNITS`, `OUTPUT_POTENTIAL_UNITS` (default: kJ/mol, type: String)"

    These keywords can be used to specify the measurement unit for potential energy 
    for the input and output files, respectively. The keyword `POTENTIAL_UNITS` only 
    applies to the data read from characteristics. Moreover, it accepts an entry for each 
    input data file. The keyword `OUTPUT_POTENTIAL_UNITS` is used for both writing the 
    results to files and plotting the results. However, the `cappa.out` file will still 
    use the standard units. Only one output unit can be provided.

    __Options__: _kJ/mol_, _J/mol_


!!! info "`VOLUME_UNITS`, `OUTPUT_VOLUME_UNITS` (default: ml/g, type: String)"

    These keywords can be used to specify the measurement unit for volume for the
    input and output files, respectively. The keyword `VOLUME_UNITS` only applies
    to the data read from characteristics. Moreover, it accepts an entry for each input 
    data file. The keyword `OUTPUT_VOLUME_UNITS` is used for both writing the results to 
    files and plotting the results. However, the `cappa.out` file will still use the 
    standard units. Only one output unit can be provided.

    __Options__: _ml/g_, _l/kg_, _cm3/g_, _dm3/kg_


!!! info "`OUTPUT_DENSITY_UNITS` (default: kg/m3, type: String)"

    The keyword `OUTPUT_DENSITY_UNITS` is used for both writing the results to files
    and plotting the results. However, the `cappa.out` file will still use the standard
    units. Only one output unit can be provided.

    __Options__: _kg/m3_
