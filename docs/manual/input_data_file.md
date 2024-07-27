---
icon: material/database-outline
---

# Input data files

Inside the configuration file, the first keyword used should always be `DATA_FILES`, as
the number of arguments provided here dictates the number of data structures crated for
storage. Data files contain information that is used to recreate adsorption properties, 
such as isotherms, isobars, isosteres, and characteristic curves. While the main format 
used for data files is the two-column format, which works for all the aforementioned 
data types, isotherms can also be provided in a variety of formats depending on the
fitting equation used. 

!!! info "`DATA_FILES` (default: None, type: String)"
    
    The paths towards the input files used for calculations. This keyword accepts
    an arbitrary number of arguments.

### Supported data types


!!! info "`DATA_TYPES` (default: None, type: String)"
    
    The format of the input file. If only one argument is provided, it is assumed that 
    all input files have the same format. Alternatively, mixed input types can be
    provided.

    __Options__: _isotherm_, _isobar_, _isostere_, _characteristic_, _langmuir_, 
    _n-langmuir_, _bet_, _anti-langmuir_, _henry_, _freundlich_, _sips_, _n-sips_, 
    _langmuir-freundlich_, _n-langmuir-freundlich_, _redlich-peterson_, _toth_, _unilan_, 
    _obrien-myers_, _quadratic_, _asymptotic-temkin_, _bingel-walton_


#### Two-column files

Two-column files are the simplest and most common way to provide data, as they are easy
to read and work with. The first line of a two-column file should always contain a
comment, most often used to specify the values recorded and the measurement units. The
data is read by parsing each line after the comment. Only the first two entries of each 
line are stored, with the rest being ignored.

For isotherms, the first column should contain the pressure and the second column should 
contain the loading. For isobars, the first column should contain the temperature and the 
second column should contain the loading. For isosteres, the first column should contain 
the temperature and the second column should contain the pressure. For characteristic 
curves, the first column should contain the adsorption potential and the second column 
should contain the adsorbed volume.

```text title="two_column_file.dat"
# first line is a comment line
prop_1  prop_2  (skipped entries)
prop_1  prop_2  (skipped entries)
```

#### Langmuir isotherm

$$
q(p) = q_{sat} \frac{b p}{1 + b p}
$$

```text title="langmuir_isotherm.dat"
p_min   p_max   num
q_sat
b
```

#### n-site Langmuir isotherm

$$
q(p) = \sum^{n}_{i=1} q_{sat, i} \frac{b_{i} p}{1 + b_{i} p}
$$

```text title="n_langmuir_isotherm.dat"
p_min   p_max   num
q_sat_1 q_sat_2 ...
b_1     b_2     ...
```

#### anti-Langmuir isotherm

$$
q(p) = \frac{a p}{1 - b p}
$$

```text title="anti_langmuir_isotherm.dat"
p_min   p_max   num
a
b
```


#### BET isotherm

$$
p_{r} = p / p_{0}
$$

$$
q(p_{r}) = q_{sat} \frac{b p_{r}}{(1 - c p_{r})(1 - c + b p_{r})}
$$

```text title="bet_isotherm.dat"
p_min   p_max   num
p_0
q_sat
b
c
```


#### Henry isotherm

$$
q(p) = a p
$$

```text title="henry_isotherm.dat"
p_min   p_max   num
a
```


#### Freundlich isotherm

$$
q(p) = a p^{1/\nu}
$$

```text title="freundlich_isotherm.dat"
p_min   p_max   num
a
nu
```

#### Sips isotherm

$$
q(p) = q_{sat} \frac{(b p)^{1/\nu}}{1 + (b p)^{1/\nu}}
$$

```text title="sips_isotherm.dat"
p_min   p_max   num
q_sat
b
nu
```


#### n-site Sips isotherm

$$
q(p) = \sum^{n}_{i=1} q_{sat, i} \frac{(b_{i} p)^{1/\nu_{i}}}{1 + (b_{i} p)^{1/\nu_{i}}}
$$

```text title="n_sips_isotherm.dat"
p_min   p_max   num
q_sat_1 q_sat_2 ...   
b_1     b_2     ...
nu_1    nu_2    ...
```


#### Langmuir-Freundlich isotherm

$$
q(p) = q_{sat} \frac{b p^{\nu}}{1 + b p^{\nu}}
$$

```text title="langmuir_freundlich_isotherm.dat"
p_min   p_max   num
q_sat 
b
nu
```


#### n-site Langmuir-Freundlich isotherm

$$
q(p) = \sum^{n}_{i=1} q_{sat, i} \frac{b_{i} p^{\nu_{i}}}{1 + b_{i} p^{\nu_{i}}}
$$

```text title="n_langmuir_freundlich_isotherm.dat"
p_min   p_max   num
q_sat_1 q_sat_2 ... 
b_1     b_2     ...
nu_1    nu_2    ...
```


#### Redlich-Peterson isotherm

$$
q(p) = \frac{a p}{1 + b p^{\nu}}
$$

```text title="redlich_peterson_isotherm.dat"
p_min   p_max   num
a
b
nu
```


#### Toth isotherm

$$
q(p) = q_{sat} \frac{b p}{[1 + (b p)^{\nu}]^{1/\nu}}
$$

```text title="toth_isotherm.dat"
p_min   p_max   num
q_sat
b
nu
```


#### Unilan isotherm

$$
q(p) = q_{sat} \frac{1}{2 \eta} \ln{\left[\frac{1 + b e^{\eta} p}{1 + b e^{-\eta} p}\right]}
$$

```text title="unilan_isotherm.dat"
p_min   p_max   num
q_sat
b
eta
```


#### O'Brien & Myers isotherm

$$
q(p) = q_{sat} \left[ \frac{b p}{1 + b p} + \sigma^{2} \frac{b p(1 - b p)}{2(1 + b p)^{3}} \right]
$$

```text title="obrien_myers_isotherm.dat"
p_min   p_max   num
q_sat
b
sigma
```


#### Quadratic isotherm

$$
q(p) = q_{sat} \frac{b p + 2 c p^{2}}{1 + b p + c p^{2}} 
$$

```text title="quadratic_isotherm.dat"
p_min   p_max   num
q_sat
b
c
```


#### Asymptotic Temkin isotherm

$$
f(p) = \frac{b p}{1 + b p} 
$$

$$
q(p) = q_{sat} f(p) + g_{sat} \theta f(p)^{2} (f(p) - 1)
$$

```text title="asymptotic_temkin_isotherm.dat"
p_min   p_max   num
q_sat
b
theta
```


#### Bingel & Walton isotherm

$$
q(p) = q_{sat} \frac{1 - \exp{[-(a+b)p]}}{1+ b/a \exp{[-(a+b)p]}}
$$

```text title="bingel_walton_isotherm.dat"
p_min   p_max   num
q_sat
a
b
```

### Identifiers

Depending on the type of the data files. additional information is required to process
the data (e.g. the temperature at which the isotherm is measured). This is done by 
using the following keywords, which can be referred to as identifiers. Identifiers are
mandatory and serve a variety of purposes, from being used in the calculations to 
helping with labeling in plots or written files. 

!!! info "`ADSORBATE` (default: None, type: String)"

    The name of the adsorbate studied. This keyword can only take one value, which is
    primarely used for naming output files and plots. Morover, it serves as a serching
    term inside the local library in case the keywords `ADSORBATE_DATA_FILE`.
    `SATURATION_PRESSURE_FILE`, `DENSITY_FILE` or  receive the argument `local`. 


!!! info "`ADSORBENT` (default: None, type: String)"

    The name of the adsorbent studied. This keyword can only take one value, which is
    only used for naming output files and plots.


!!! info "`LOADINGS` (default: None, type: Float)"

    The loadings at which the isosteres are measured. This keyword should take as many
    arguments as the number of data files provided. If mixed data types are used,
    placeholder values should be provided for the non-isostere entries.


!!! info "`PRESSURES` (default: None, type: Float)"

    The pressures at which the isobars are measured. This keyword should take as many
    arguments as the number of data files provided. If mixed data types are used,
    placeholder values should be provided for the non-isobar entries.


!!! info "`TEMPERATURES` (default: None, type: Float)"

    The temperatures at which the isosteres are measured. This keyword should take as 
    many arguments as the number of data files provided. If mixed data types are used,
    placeholder values should be provided for the non-isotherm entries.



### Plotting data

In some cases, it is helpful to plot the input data to use it as a reference for the 
results computed. Below, a list of keywords that affect plotting is provided. 

!!! info "`PLOT_DATA` (default: no, type: String)"
    
    Specifies if the input data should be plotted. Keep in mind that plotting the data 
    does not entail that the plot is immediately displayed or saved. It is technically
    possible to plot the data and then ignore the resulting plot, so the user has to 
    pay attention to the keywords used.
    
    __Options__: _yes_, _no_


!!! info "`LOGARITHMIC_PLOT` (default: no, type: String)"
    
    Specifies if the x-axis of the plot should be in logarithmic scale. This only
    affects isotherms. 
    
    __Options__: _yes_, _no_


!!! info "`SAVE_DATA_PLOT` (default: no, type: String)"
    
    Specifies if the plot should be saved inside the __Plots__ directory. 
    
    __Options__: _yes_, _no_


!!! info "`SHOW_PLOTS` (default: no, type: String)"
    
    Specifies if a pop up with all of the created plots is displayed. This keyword
    affects all of the plots created throughout the run. 
    
    __Options__: _yes_, _no_