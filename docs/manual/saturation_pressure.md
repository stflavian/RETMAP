---
icon: material/water-thermometer-outline
---

# Saturation pressure


!!! info "`ADSORBATE_SATURATION_PRESSURE` (default: None, type: String)"

    Specifies the saturation pressure model used for the calculations. This keyword
    can take as many arguments as the number of input files provided. However, since
    the keyword `SATURATION_PRESSURE_FILE` can only take one argument, using the
    _extrapolation_ method multiple times would produce the same results. 
        
    __Options__: _dubinin_, _amankwah_, _extrapolation_, _polynomial_water_,
    _peng_robinson_, _peng_robinson_extrapolation_, _prsv1_extrapolation_,
    _prsv2_extrapolation_, _widom_banuti_, _critical_isochore_,


!!! info "`SATURATION_PRESSURE_FILE` (default: None, type: String)"

    The path towards the file that contains the saturation pressure data. This keyword 
    accepts only one argument. If the argument provided is _local_, the program will
    search its internal library for pre-made saturation pressure data file based on the 
    argument given to the `ADSORBATE` keyword. 


!!! info "`AMANKWAH_EXPONENT` (default: 3.0, type: Float)"

    The exponent used for the _amankwah_ method.


### Computing saturation pressure curves



!!! info "`COMPUTE_SATURATION_PRESSURE_CURVE` (default: no, type: String)"

    Specify if a stuartion pressure curve should be computed and saved using the
    first calculation method provided to `ADSORBATE_SATURATION_PRESSURE`.

    __Options__: _yes_, _no_


!!! info "`SATURATION_PRESSURE_RANGE` (default: None, type: Tuple{Float})"

    The range in which the saturation pressure is computed.


!!! info "`NUMBER_SATURATION_PRESSURE_POINTS` (default: 50, type: Integer)"

    The number of points used to for computing the saturation pressure curve.
