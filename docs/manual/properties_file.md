# Properties file

When processing adsorption data, molecular properties are required for preforming most
calculations. These molecular properties are usually stored in a properties file,
which is then provided to the configuration file using the key word 
`ADSORBATE_DATA_FILE`. Properties files follow the same structure as configuration
files, and are usually saved using the `.prop` file extension.

> Warning: The interpreter allows for one or more keywords to be skipped inside the
> properties files. Should that happen, individual warnings are raised during the
> run for each individual missing keyword. Moreover, if one of the missing keywords
> is required for calculations, the program will fail.

The following keywords should be provided, in no particular order, inside the 
properties files:

NAME (default: None, type: String)
: The name of the molecule for which the properties are provided. This string is not
used in the code and is only provided for readability reasons.

MOLECULAR_MASS (default: None, type: Float)
: The molecular mass of the molecule in g/mol.
    
PRESSURE_CRITICAL (default: None, type: Float)
: The pressure at the critical point in MPa.
    
TEMPERATURE_CRITICAL (default: None, type: Float)
: The temperature at the critical point in K.
    
TEMPERATURE_BOILING (default: None, type: Float)
: The temperature at the boiling point in K.
    
TEMPERATURE_TRIPLE_POINT (default: None, type: Float)
: The temperature at the triple point in K.
    
DENSITY_BOILING (default: None, type: Float)
: The density at the boiling point in kg/m3.
    
ACENTRIC_FACTOR (default: None, type: Float)
: The acentric factor of the molecule.
    
PRSV_KAPPA1 (default: None, type: Float)
: The first fitting parameter of the PRSV and PRSV2 equations. A list of given values 
can be found in [PRSV: An improved peng—Robinson equation of state for pure compounds 
and mixtures](https://onlinelibrary.wiley.com/doi/10.1002/cjce.5450640224) by R. 
Stryjek and J. H. Vera.
    
PRSV_KAPPA2 (default: None, type: Float)
: The second fitting parameter of the PRSV2 equation. A list of given values can be 
found in [PRSV2: A cubic equation of state for accurate vapor—liquid equilibria 
calculations](https://onlinelibrary.wiley.com/doi/10.1002/cjce.5450640516) by R. Stryjek 
and J. H. Vera.
    
PRSV_KAPPA3 (default: None, type: Float)
: The third fitting parameter of the PRSV2 equation. A list of given values can be 
found in [PRSV2: A cubic equation of state for accurate vapor—liquid equilibria 
calculations](https://onlinelibrary.wiley.com/doi/10.1002/cjce.5450640516) by R. Stryjek 
and J. H. Vera.
    
__SOURCE__ (default: None, type: String)
: The source from which the molecular properties are taken. This string is not
used in the code and is only provided for transparency reasons. The format used for
this is left to the user.
