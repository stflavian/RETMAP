# Configuration file

When running `cappa`, the only required argument is the name of the configuration file, 
usually called `config.in`. This contains all the information needed by the application 
to run, such as the location of the needed files and the expected output.

The configuration file is a plain text file structured as:

> <pre>
> KEY_WORD           argument
> </pre>

where `KEY_WORD` represents the option needed to be edited and `argument` the respective
option. The distance between the two entries should be at least one empty space. In case
multiple entries are needed (e.g. when multiple data files should be processed at the 
same time), the arguments should be added on the same line with at least one space in
between them.

> <pre>
> KEY_WORD           argument1  argument2  ...  argumentN
> </pre>

In case multiple data files are processed at the same time, the user can choose to use
either the same number of arguments for all settings or only one argument. In the later
case it is assumed that all the arguments should take the same value.

> <pre>
> KEY_WORD1          argument1  argument2  ...  argumentN
> KEY_WORD2          argument1 (argument1) ... (argument1) 
> </pre>

There are however a few special keys that do not abide by this rule (e.g. those that
denote ranges, prediction parameters, etc.), which will be discussed in their respective
parts.

To allow for easy organization of the `config.in` file, empty lines can be used to 
separate different settings blocks. Also, if the first word of the line is not in the 
list of recognized key-words, the line will be ignored. This behaviour allows for 
comments to be used inside the configuration file, with little limitations. Still, it 
is recommended that comments start with some sequence of symbols (e.g. `//`) to ensure
that the first word does not match an implemented key-word. This is also good for
future proofing, as no key-words starting with symbols are planned to be
implemented.

> <pre>
> // This is the first block of settings
> KEY_WORD1          argument1  argument2  ...  argumentN
> KEY_WORD2          argument1
>
> // This is the second block of settings
> KEY_WORD3          argument1
> </pre>

Both key-words and arguments are case-sensitive, with the former written in uppercase
and the latter in lower case. For the boolean clauses (a.i. those that evaluate for
either `True` or `False`), `yes` represents `True` while any other string represents
`False`. However, it is recommended that `no` is used for false to future-proof the
file.