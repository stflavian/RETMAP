# Units

CAPPA works with six main metrics: pressure, temperature, loading, potential, volume,
and density. These metrics are stored and processed using the following units: MPa, K,
mg/g, kJ/mol, ml/g, and kg/m3. As such, all values provided inside the configuration 
file (e.g. `PREDICTION_TEMPERATURES`) should be provided using the aforementioned units.
However, the user can specify inside the configuration file if the data read from the
input files is provided in different units, in which case conversion factors are used 
to reach the default units. Similarly, the user can specify the units used for writing
the results. 

The following units are allowed:

PRESSURE
: MPa, kPa, Pa, bar, atm, Torr, mmHg

LOADING
: mg/g, g/kg, mol/kg, mmol/g, cm3/g, mm3/mg, dm3/kg, l/kg, ml/g

POTENTIAL
: kJ/mol, J/mol

VOLUME
: ml/g, l/kg, cm3/g, dm3/kg
