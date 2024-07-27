---
icon: material/scale
---

# Density

!!! info "`ADSORBATE_DENSITY` (default: None, type: String)"

    Specifies the density model used for the calculations. This keyword can take as many 
    arguments as the number of input files provided. However, since the keyword 
    `DENSITY_FILE` can only take one argument, using the _extrapolation_ method multiple 
    times would produce the same results. 

    __Options__: _empirical_, _hauer_, _ozawa_, _extrapolation_


!!! info "`DENSITY_FILE` (default: None, type: String)"
 
    The path towards the file that contains the density data. This keyword accepts only 
    one argument. If the argument provided is _local_, the program will search its 
    internal library for pre-made density data file based on the argument given to the 
    `ADSORBATE` keyword. 


!!! info "`THERMAL_EXPANSION_COEFFICIENT` (default: 0.00165, type: Float)"

    The thermal expansion coefficient used for the _hauer_ and _ozawa_ models.


### Computing density curves

!!! info "`COMPUTE_DENSITY_CURVE` (default: no, type: String)"
    
    Specify if a density curve should be computed and saved using the first calculation 
    method provided to `ADSORBATE_DENSITY`.

    __Options__: _yes_, _no_


!!! info "`DENSITY_RANGE` (default: None, type: Tuple{Float})"

    The range in which the density is computed.


!!! info "`NUMBER_DENSITY_POINTS` (default: 50, type: Integer)"
    
    The number of points used to for computing the density curve.
