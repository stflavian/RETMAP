---
icon: material/chart-bell-curve-cumulative
---

# Prediction-making


### Isotherm prediction

!!! info "`PREDICT_ISOTHERMS` (default: no, type: String)"

    Specify if isotherms should be predicted using the first characteristic curve 
    provided or computed.

    __Options__: _yes_, _no_


!!! info "`PREDICTION_TEMPERATURES` (default: None, type: Array{Float})"
    
    The temperatures at which the isotherms should be computed.    

         
!!! info "`PREDICTION_PRESSURE_RANGE` (default: None, type: Tuple{Float})"

    The pressure range covered by the isotherms. If the values provided are out of
    reach, the program will fallback to the minimum and maximum pressures that can
    be reached for the given adsorption potential range.
     
   
!!! info "`NUMBER_PRESSURE_POINTS` (default: 50, type: Integer)"
    
    The number of points used for the isotherm.
     

!!! info "`SAVE_PREDICTED_ISOTHERMS_DATA` (default: no, type: String)"

    Specify if the predicted isotherms should be written to files.

    __Options__: _yes_, _no_

   
!!! info "`PLOT_PREDICTED_ISOTHERMS` (default: yes, type: String)"

    Specify if the predicted isotherms should be plotted.  

    __Options__: _yes_, _no_

    
    
!!! info "`SAVE_PREDICTED_ISOTHERMS_PLOT` (default: no, type: String)"

    Specify if the predicted isotherms plot should be saved inside the __Plots__ 
    directory. 

    __Options__: _yes_, _no_


### Isobar prediction

!!! info "`PREDICT_ISOBARS` (default: no, type: String)"

    Specify if isobars should be predicted using the first characteristic curve 
    provided or computed.

    __Options__: _yes_, _no_


!!! info "`PREDICTION_PRESSURES` (default: None, type: Array{Float})"

    The pressures at which the isobars should be computed.    
        
 
!!! info "`PREDICTION_TEMPERATURE_RANGE` (default: None, type: Tuple{Float})"
    
    The temperature range covered by the isobars. If the values provided are out of
    reach, the program will fallback to the minimum and maximum temperatures that can
    be reached for the given adsorption potential range.
        

!!! info "`NUMBER_TEMPERATURE_POINTS` (default: 50, type: Integer)"

    The number of points used for the isobars.


!!! info "`SAVE_PREDICTED_ISOBARS_DATA` (default: no, type: String)"

    Specify if the predicted isobars should be written to files.

    __Options__: _yes_, _no_


!!! info "`PLOT_PREDICTED_ISOBARS` (default: yes, type: String)"

    Specify if the predicted isobars should be plotted.  

    __Options__: _yes_, _no_       
 

!!! info "`SAVE_PREDICTED_ISOBARS_PLOT` (default: no, type: String)"

    Specify if the predicted isobars plot should be saved inside the __Plots__ 
    directory. 

    __Options__: _yes_, _no_


### Isostere prediction

!!! info "`PREDICT_ISOSTERES` (default: no, type: String)"

    Specify if isosteres should be predicted using the first characteristic curve 
    provided or computed.

    __Options__: _yes_, _no_


!!! info "`PREDICTION_LOADINGS` (default: None, type: Array{Float})"

    The loadings at which the isosteres should be computed.    
       
  
!!! info "`PREDICTION_ISOSTERE_RANGE` (default: None, type: Tuple{Float})"

    The temperature range covered by the isosteres. If the values provided are out of
    reach, the program will fallback to the minimum and maximum temperatures that can
    be reached for the given adsorbed volume range.
      
  
!!! info "`NUMBER_ISOSTERE_POINTS` (default: 50, type: Integer)"

    The number of points used for the isosteres.
       

!!! info "`SAVE_PREDICTED_ISOSTERES_DATA` (default: no, type: String)"

    Specify if the predicted isosteres should be written to files.

    __Options__: _yes_, _no_ 


!!! info "`PLOT_PREDICTED_ISOSTERES` (default: yes, type: String)"

    Specify if the predicted isosteres should be plotted.  

    __Options__: _yes_, _no_
      
  
!!! info "`SAVE_PREDICTED_ISOSTERES_PLOT` (default: no, type: String)"

    Specify if the predicted isosteres plot should be saved inside the __Plots__ 
    directory. 

    __Options__: _yes_, _no_

