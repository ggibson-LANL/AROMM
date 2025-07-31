This folder contains the following folders and scripts

[output](https://drive.google.com/open?id=1weENhS7XOsUsnHUfVpwx97HDrWEX-xs6&usp=drive_copy)

[Porcupine_river.py](https://drive.google.com/open?id=159tiD4Dl9p4UgBvfs4wEhEdrq13T9yBC&usp=drive_copy)

[output_1_4000](https://drive.google.com/open?id=10MdPJwYUjlzQI-CTaS6HkogbcDUG-smA&usp=drive_copy)

[Data_Ana_2.py](https://drive.google.com/open?id=1o7boLka0LPjI_4o2TJNS6NHsNrGhu7J0&usp=drive_copy)

[Data_Ana_3.py](https://drive.google.com/open?id=1sCYugXsY4o5vsfafb6aPVCzeoyUDhwLQ&usp=drive_copy)

[Koyukuk_river.py](https://drive.google.com/open?id=1ah1DMEI0XzXD_AzNJX5tGddwzq52hSAr&usp=drive_copy)

[Pelly_river.py](https://drive.google.com/open?id=1f0WLSkTmAh0DfmaqqC98Oqly1vmx00_a&usp=drive_copy)

[Yukon_river.py](https://drive.google.com/open?id=1eXYNDPGl30xynB3LfLI90FcZUu43wGYP&usp=drive_copy)

[READ](https://drive.google.com/open?id=1hW_pxPA8AA-iIF2p9JRBKDkxgQgevl-XPiAU2otG8ww&usp=drive_copy)ME

[river_initial_monac_code.py](https://drive.google.com/open?id=1p4xGzHb1JxoInmQC--59CrlipCWyhbUb&usp=drive_copy)

[Stewart_river.py](https://drive.google.com/open?id=1Ti9PO7bvwwSmvVvD3eP0KekaBmBUnQ8_&usp=drive_copy)

[Tanana_river.py](https://drive.google.com/open?id=1852ufT2-fHT33a6VREAfEfOoZbS4X4T2&usp=drive_copy)

[Teslin_river_plots.py](https://drive.google.com/open?id=1tilWPCdmHuYwFyK0jqKBQq5VWUVH0-Ex&usp=drive_copy)

[Teslin_river.py](https://drive.google.com/open?id=1iD6flY3aNLN3pR4a52PAIL-OyY4a9KCZ&usp=drive_copy)

[White_Donjec_river.py](https://drive.google.com/open?id=1aGeyoAZgVCyL2XzAb-xQ7BKjHyuRSHi3&usp=drive_copy)

[Yukon_river_plot.py](https://drive.google.com/open?id=1L1w0OZTZDraJHCNU75_Ya3WL95c3RUK3&usp=drive_copy)

\*\*\*\*\*\*\*\*\*\*\*\*

This folder contains the following: a README file, folders, and Python scripts

Readme file

[READ](https://drive.google.com/open?id=1hW_pxPA8AA-iIF2p9JRBKDkxgQgevl-XPiAU2otG8ww&usp=drive_copy)ME

Folders

[output](https://drive.google.com/open?id=1weENhS7XOsUsnHUfVpwx97HDrWEX-xs6&usp=drive_copy)

The output folder contains an analysis folder and separate folders for the Yukon River tributaries and the main Yukon River stem.

[output_1_4000](https://drive.google.com/open?id=10MdPJwYUjlzQI-CTaS6HkogbcDUG-smA&usp=drive_copy)

This output folder specifically has 4000 simulation runs outputs, contains an analysis folder, and separate folders for the Yukon River tributaries and the main Yukon River stem.

Python Codes

[river_initial_monac_code.py](https://drive.google.com/open?id=1p4xGzHb1JxoInmQC--59CrlipCWyhbUb&usp=drive_copy)

This code allows the user to manually adjust the average values and boundaries (minimum and maximum) of the input parameters and also model parameters. It generates random values for each parameter within the specified boundaries. The number of simulations can be controlled by changing the num_samples parameter. Also, it contains codes that create distribution plots for each input parameter that show their distribution types, max and min boundary values that are distributed around their mean values. For some distribution functions, code has been written to truncate negative input parameter values.

The outputs of this script are,

[CDOM_component_fractions.csv](https://drive.google.com/open?id=1e-25TrWO72cwrdX1J8fnobwZYm--_sUf&usp=drive_copy)

[Chemical_fraction.csv](https://drive.google.com/open?id=1QlpA8AgiG62V_nM_A1RRJSf8VGUFHhda&usp=drive_copy)

[Dilution_fractions.csv](https://drive.google.com/open?id=1JGUnzcTTAi6hMbuh2tDD-FcpxbqoJxx1&usp=drive_copy)

[Initial_DOC_Values.csv](https://drive.google.com/open?id=1SoFvxgtjyeqqCdAAjSEd74HogNNqdQNQ&usp=drive_copy)

[Production_Values.csv](https://drive.google.com/open?id=14Dr89eDtCHnkr61wR5A-uLxwuVD4lQub&usp=drive_copy)

[Tau_values.csv](https://drive.google.com/open?id=1ZSi2diNgDxTtqplbr8Hi3DrExtzfGifQ&usp=drive_copy)

[Velocity_Values.csv](https://drive.google.com/open?id=1Vc1BbRt_DLGwmZiBQO4QrXMKPSjLuKn2&usp=drive_copy) , which are saved as CSV files in the output/random_initial_val/ directory ([random_initial_val](https://drive.google.com/open?id=1Dm2LEWdGxyU8KJuFnaA1zm_YHySKtGdD&usp=drive_copy)).

This also generates **distribution plots of each input parameter**.

[Teslin_river.py](https://drive.google.com/open?id=1iD6flY3aNLN3pR4a52PAIL-OyY4a9KCZ&usp=drive_copy)

This is the Teslin(first) tributary module of the Yukon River. It contains Python code that includes:

- A function for chemical decay,
- A loop for Monte Carlo analysis,
- Code to generate values for imax (maximum number of iterations), tmax (maximum number of time steps), and dmax (maximum distance of the tributary),
- Time step calculations (time_calculation),
- Distance calculations (distance_calculation),
- Creation of arrays for C_Conc, Tau, and Prod,
- Initialization of C_Conc, Tau, and Prod values for the chemical
- The main river code
- Creating CSV files that have the tributary river mouth concentrations of each macromolecule

_Output of the code_

The output of this script is the river mouth concentration of the tributary, which is saved as a CSV file in the output/Teslin/ directory([teslin](https://drive.google.com/open?id=1klacGE0QGnZlPzFaFFloTXzwdm9EzAih&usp=drive_copy)).

[Teslin_river_plots.py](https://drive.google.com/open?id=1tilWPCdmHuYwFyK0jqKBQq5VWUVH0-Ex&usp=drive_copy)

This code is similar to the “[Teslin_river.py](https://drive.google.com/open?id=1yZtmOod1TNqkFkuEH37wTpxt1AmiQutp&usp=drive_copy)”, but it contains an additional code part that generates plots of Chemical concentration vs distance.

[Koyukuk_river.py](https://drive.google.com/open?id=1ah1DMEI0XzXD_AzNJX5tGddwzq52hSAr&usp=drive_copy)

This is the Koyukuk tributary module of the Yukon River. It contains Python code that includes:

- A function for chemical decay,
- A loop for Monte Carlo analysis,
- Code to generate values for imax (maximum number of iterations), tmax (maximum number of time steps), and dmax (maximum distance of the tributary),
- Time step calculations (time_calculation),
- Distance calculations (distance_calculation),
- Creation of arrays for C_Conc, Tau, and Prod,
- Initialization of C_Conc, Tau, and Prod values for the chemical
- The main river code
- Creating CSV files that have the tributary river mouth concentrations of each macromolecule

_Output of the code_

The output of this script is the river mouth concentration of the tributary, which is saved as a CSV file in the output/Koyukuk/ directory ([koyukuk](https://drive.google.com/open?id=1pAD2YLfUoc0xTJMDK8G8Ha5wMKTD3C2P&usp=drive_copy))

[Pelly_river.py](https://drive.google.com/open?id=1f0WLSkTmAh0DfmaqqC98Oqly1vmx00_a&usp=drive_copy)

This is the Pelly tributary module of the Yukon River. It contains Python code that includes:

- A function for chemical decay,
- A loop for Monte Carlo analysis,
- Code to generate values for imax (maximum number of iterations), tmax (maximum number of time steps), and dmax (maximum distance of the tributary),
- Time step calculations (time_calculation),
- Distance calculations (distance_calculation),
- Creation of arrays for C_Conc, Tau, and Prod,
- Initialization of C_Conc, Tau, and Prod values for the chemical
- The main river code
- Creating CSV files that have the tributary river mouth concentrations of each macromolecule

_Output of the code_

The output of this script is the river mouth concentration of the tributary, which is saved as a CSV file in the output/Pelly/ directory ([pelly](https://drive.google.com/open?id=1Ra3cSI4s0t1xR6fuheNtgKivDy-PKR8Q&usp=drive_copy))

[Porcupine_river.py](https://drive.google.com/open?id=159tiD4Dl9p4UgBvfs4wEhEdrq13T9yBC&usp=drive_copy)

This is the Porcupine tributary module of the Yukon River. It contains Python code that includes:

- A function for chemical decay,
- A loop for Monte Carlo analysis,
- Code to generate values for imax (maximum number of iterations), tmax (maximum number of time steps), and dmax (maximum distance of the tributary),
- Time step calculations (time_calculation),
- Distance calculations (distance_calculation),
- Creation of arrays for C_Conc, Tau, and Prod,
- Initialization of C_Conc, Tau, and Prod values for the chemical
- The main river code
- Creating CSV files that have the tributary river mouth concentrations of each macromolecule

_Output of the code_

The output of this script is the river mouth concentration of the tributary, which is saved as a CSV file in the output/Porcupine/ directory ([porcupine](https://drive.google.com/open?id=1TJvoK7JuJujnUUrnUW01sAQ5aCb6R5w4&usp=drive_copy))

[Stewart_river.py](https://drive.google.com/open?id=1Ti9PO7bvwwSmvVvD3eP0KekaBmBUnQ8_&usp=drive_copy)

This is the Stewart tributary module of the Yukon River. It contains Python code that includes:

- A function for chemical decay,
- A loop for Monte Carlo analysis,
- Code to generate values for imax (maximum number of iterations), tmax (maximum number of time steps), and dmax (maximum distance of the tributary),
- Time step calculations (time_calculation),
- Distance calculations (distance_calculation),
- Creation of arrays for C_Conc, Tau, and Prod,
- Initialization of C_Conc, Tau, and Prod values for the chemical
- The main river code
- Creating CSV files that have the tributary river mouth concentrations of each macromolecule

_Output of the code_

The output of this script is the river mouth concentration of the tributary, which is saved as a CSV file in the output/Stewart/ directory ([stewart](https://drive.google.com/open?id=1rQiTqlf9RXNgyu-MTKQGYqmeB7-J5pQt&usp=drive_copy))

[Tanana_river.py](https://drive.google.com/open?id=1852ufT2-fHT33a6VREAfEfOoZbS4X4T2&usp=drive_copy)

This is the Tanana tributary module of the Yukon River. It contains Python code that includes:

- A function for chemical decay,
- A loop for Monte Carlo analysis,
- Code to generate values for imax (maximum number of iterations), tmax (maximum number of time steps), and dmax (maximum distance of the tributary),
- Time step calculations (time_calculation),
- Distance calculations (distance_calculation),
- Creation of arrays for C_Conc, Tau, and Prod,
- Initialization of C_Conc, Tau, and Prod values for the chemical
- The main river code
- Creating CSV files that have the tributary river mouth concentrations of each macromolecule

_Output of the code_

The output of this script is the river mouth concentration of the tributary, which is saved as a CSV file in the output/Tanana/ directory ([tanana](https://drive.google.com/open?id=18hAkW-uXODzh79xkcndLdQLKVFXp2dU4&usp=drive_copy))

[White_Donjec_river.py](https://drive.google.com/open?id=1aGeyoAZgVCyL2XzAb-xQ7BKjHyuRSHi3&usp=drive_copy)

This is the White_Donjec tributary module of the Yukon River. It contains Python code that includes:

- A function for chemical decay,
- A loop for Monte Carlo analysis,
- Code to generate values for imax (maximum number of iterations), tmax (maximum number of time steps), and dmax (maximum distance of the tributary),
- Time step calculations (time_calculation),
- Distance calculations (distance_calculation),
- Creation of arrays for C_Conc, Tau, and Prod,
- Initialization of C_Conc, Tau, and Prod values for the chemical
- The main river code
- Creating CSV files that have the tributary river mouth concentrations of each macromolecule

_Output of the code_

The output of this script is the river mouth concentration of the tributary, which is saved as a CSV file in the output/White_Donjec/ directory ([white+donjec](https://drive.google.com/open?id=1OO5jz7pSOOUK0Fzd45SOeGE_XjG7pl33&usp=drive_copy))

[Yukon_river.py](https://drive.google.com/open?id=1eXYNDPGl30xynB3LfLI90FcZUu43wGYP&usp=drive_copy)

This is the main module of the Yukon River, which represents the Yukon River stem. This is the script that calls all the tributary modules to the main river code in this script. It contains Python code that includes:

- A function for chemical decay,
- A loop for Monte Carlo analysis,
- Code to generate values for imax (maximum number of iterations), tmax (maximum number of time steps), and dmax (maximum distance of the tributary),
- Time step calculations (time_calculation),
- Distance calculations (distance_calculation),
- Creation of arrays for C_Conc, Tau, and Prod,
- Initialization of C_Conc, Tau, and Prod values for the chemical
- The main river code
- Creating CSV files that have the tributary river mouth concentrations of each macromolecule

But since this is the main river code of the stem, it contains the dilution function at the node points where tributaries are connecting to the main stem.

_Output of the code_

The output of this script is the river mouth concentration of the tributary, which is saved as a CSV file in the output/White_Donjec/ directory ([yukon](https://drive.google.com/open?id=1fWe2rrwjifOMPZq7bDnswTYYqGr7fHW4&usp=drive_copy))

[Yukon_river_plot.py](https://drive.google.com/open?id=1L1w0OZTZDraJHCNU75_Ya3WL95c3RUK3&usp=drive_copy)

This is the main module of the Yukon River, which represents the Yukon River stem. It contains Python code that includes:

- A function for chemical decay,
- A loop for Monte Carlo analysis,
- Code to generate values for imax (maximum number of iterations), tmax (maximum number of time steps), and dmax (maximum distance of the tributary),
- Time step calculations (time_calculation),
- Distance calculations (distance_calculation),
- Creation of arrays for C_Conc, Tau, and Prod,
- Initialization of C_Conc, Tau, and Prod values for the chemical
- The main river code
- Creating CSV files that have the tributary river mouth concentrations of each macromolecule

_Output of the code_

In addition to the Yukon_river.py script, this Python code generates a CSV file that stores output data from all iterations across all time steps for n sample runs (e.g., from a Monte Carlo simulation).

For example, in a single sample run, the concentration is recorded at each time step of the simulation. This data can then be used to plot distance vs. concentration, illustrating how concentration changes along the length of the tributary or river over time.

[Data_Ana_2.py](https://drive.google.com/open?id=1o7boLka0LPjI_4o2TJNS6NHsNrGhu7J0&usp=drive_copy)

This is one of the main codes. Once you set the initial average values and parameter boundaries in River_initial_monac_code.py, there is no need to manually run the separate tributary or Yukon_river.py scripts.

_Output of the code_

The code automatically calls the necessary modules and generates the following:

- **Plots of DOC concentration vs. input variables  
    **

[Data_Ana_3.py](https://drive.google.com/open?id=1sCYugXsY4o5vsfafb6aPVCzeoyUDhwLQ&usp=drive_copy)

This is the main (ultimate) code. Once you set the initial average values and parameter boundaries in River_initial_monac_code.py, there is no need to manually run the separate tributary or Yukon_river.py scripts. The difference between [Data_Ana_2.py](https://drive.google.com/open?id=1nBc7Ej53XkLhNfvfU1aTow3re9QpVhyt&usp=drive_copy) and [Data_Ana_3.py](https://drive.google.com/open?id=1kiEoYr5F1G7YmkmSdB-ibfbcrwM-GuTv&usp=drive_copy) is that [Data_Ana_3.py](https://drive.google.com/open?id=1kiEoYr5F1G7YmkmSdB-ibfbcrwM-GuTv&usp=drive_copy) generates **Plots of macromolecular concentration vs. input variables** for each macromolecule addition to the DOC, and also all the figures that are generated during the run are going to be stored in the figure folder ([Figures](https://drive.google.com/open?id=1Uve_PV1XoI7c2BsflqX19tOSpmziOcOQ&usp=drive_copy)) in the /output/analysis/ directory

_Output of the code_

The code automatically calls the necessary modules and generates the following:

- **Plots of DOC concentration vs. input variables**
- **Plots of macromolecular concentration vs. input variables** for each macromolecule

(This is the additional feature of the [Data_Ana_3.py](https://drive.google.com/open?id=1kiEoYr5F1G7YmkmSdB-ibfbcrwM-GuTv&usp=drive_copy) concerning [Data_Ana_2.py](https://drive.google.com/open?id=1nBc7Ej53XkLhNfvfU1aTow3re9QpVhyt&usp=drive_copy) )

- **Plots of macromolecular concentration vs. input variables** for each m

This code also generates an analysis folder, located at /output/analysis/, which combines all the CSV files needed for further analysis. ([analysis](https://drive.google.com/open?id=1He76INuI5yqH-0WZWyHgOAkQjcv4WO5H&usp=drive_copy)). In the analysis folder, in addition to these CSV files, a folder containing figures is created.

[Figures](https://drive.google.com/open?id=1Uve_PV1XoI7c2BsflqX19tOSpmziOcOQ&usp=drive_copy)

[CDOM_component_fractions.csv](https://drive.google.com/open?id=1b8C0S76dJDRNl4O0Nhy6UvRcKoNwwQru&usp=drive_copy)

[Chemical_fraction.csv](https://drive.google.com/open?id=1t-32pAAgJ6JYj939_ZfQuDHapuBuvtAV&usp=drive_copy)

[Dilution_fractions.csv](https://drive.google.com/open?id=1PJAx8wG5jYHKPc2k4h800d9NlT66i1HW&usp=drive_copy)

[Initial_DOC_Values.csv](https://drive.google.com/open?id=1YXfbgyWi_WFQAEXivgaeYBveUntd983S&usp=drive_copy)

[InputData.csv](https://drive.google.com/open?id=1xwXesAoDGHGIBQpc7SOk-avMxL0_022X&usp=drive_copy)

[Koyukuk_River_Mouth_Values.csv](https://drive.google.com/open?id=1Djn0cvSmKWVAmzVsNEBxQfTkmNfZugqp&usp=drive_copy)

[OutputData.csv](https://drive.google.com/open?id=1jeUCmj3rAozSE3NEWu04xwZW1Wrz_3bS&usp=drive_copy)

[Pelly_River_Mouth_Values.csv](https://drive.google.com/open?id=1nyE278ou3FxDzN6h1hkAaq_AIVlXMCyv&usp=drive_copy)

[Porcupine_River_Mouth_Values.csv](https://drive.google.com/open?id=1fghg38NIjkXDBT-Pr0tzgQYbxaxKqV5h&usp=drive_copy)

[Production_Values.csv](https://drive.google.com/open?id=1lJS2D8zihxFaDP7148b-iJe8UAQt2HwM&usp=drive_copy)

[Stewart_River_Mouth_Values.csv](https://drive.google.com/open?id=1d8zcNCF9hqh0HQ4RPS6YILKl9tQkquWb&usp=drive_copy)

[Tanana_River_Mouth_Values.csv](https://drive.google.com/open?id=1sT7QHn8TwaneN6hEQbB2Bp6vA0nSqUwk&usp=drive_copy)

[Tau_Values.csv](https://drive.google.com/open?id=1CF-lK23HHMcr2GXhIZhCs5njhqJpZeyS&usp=drive_copy)

[Teslin_River_Mouth_Values.csv](https://drive.google.com/open?id=10ApeleN9VTEKvG2QuWdvCZsVsBn8kv9g&usp=drive_copy)

[Velocity_Values.csv](https://drive.google.com/open?id=1o9aSdiBKISLSwRLSKUFUiV_4ggCUsquS&usp=drive_copy)

[White+Donjec_River_Mouth_Values.csv](https://drive.google.com/open?id=1wH6C0bv_Qa-g4oBAdYuEUF-rQSeERnz0&usp=drive_copy)

[Yukon_River_Mouth_Values.csv](https://drive.google.com/open?id=1U89FDyuZVO5NjDq_iFidp3nOCqaZiHeD&usp=drive_copy)