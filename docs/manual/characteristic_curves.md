---
icon: material/cached
---

# Computing characteristic curves

The most important feature of CAPPA is the ability to compute adsorption characteristic 
curves based on the adsorption properties provided. 

!!! info "`COMPUTE_CHARACTERISTIC_CURVE` (default: no, type: String)"
    
    Specify if characteristic curves should be computed for all input files provided.
    For the computing process, saturation pressure models and density models should be
    specified using the keywords `ADSORBATE_SATURATION_PRESSURE` and `ADSORBATE_DENSITY`.

    __Options__: _yes_, _no_


!!! info "`SAVE_CHARACTERISTIC_CURVE_DATA` (default: no, type: String)"

    Specify if the resulting characteristic curves should be written to files.

    __Options__: _yes_, _no_


!!! info "`PLOT_CHARACTERISTIC_CURVE` (default: yes, type: String)"

    Specify if the resulting characteristic curves should be plotted.    

    __Options__: _yes_, _no_


!!! info "`SAVE_CHARACTERISTIC_CURVE_PLOT` (default: no, type: String)"

    Specify if the characteristic cruve plot should be saved inside the __Plots__ 
    directory. 

    __Options__: _yes_, _no_