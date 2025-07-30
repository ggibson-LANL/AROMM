# AROMM
OPEN SOURCE ASSERTION: The ‘Yukon River Chemistry Model’ was approved for Open-Source Assertion on 4/9/2025 and assigned #O4870 

## Model Overview:
This project focused on modeling and analyzing the transport and transformation of Dissolved Organic Carbon (DOC) and its macromolecular composition in Arctic rivers, with a particular emphasis on the Yukon River. A central goal was to develop a process-based river model capable of simulating the reactivity and fate of dissolved organic components under conditions characteristic of Arctic river systems. The specific objectives were to: (1) model the chemical reactivity of DOC components; (2) parameterize the model to reflect the physical and biogeochemical characteristics of the Yukon River; (3) estimate DOC and it’s macromolecular component concentrations at the river mouth; and (4) evaluate spatial and temporal changes in chemical concentrations along the river continuum.

The AROMM model developed under this project simulates the evolution of dissolved organic chemicals along the river, from the headwater to the river mouth. 

Model code is written in Python and runs on Mac OS. The model should be compatible with any operating system. Python ≥ 3.9 is required, along with the following Python packages:
Numpy - Numerical operations
Pandas - Data handling and input/output
Satplotlib - Plotting and visualization
Scipy - Statistical distributions and optimization
pyDOE -Design of experiments (Latin Hypercube Sampling, content, and the structure of the model)

## Points of Contact: 
Repository manager: Dr. Georgina Anne Gibson  
Affiliation: Scientist, CCS-2 Division, Los Alamos National Laboratory (LANL)  
Email: ggibson@lanl.gov 

Model Developer: Dr. Amadini Mendis Jayasinghe  
Affiliation*: Consultant Senior Lecturer, CINEC Campus, Malabe, Sri Lanka  
Email: amadini.mj@gmail.com
* Amadini Mendis Jayasinghe developed the model while serving as a Postdoctoral Research Associate at Los Alamos National Laboratory (LANL), CCS-2 Division, Los Alamos, NM, USA. 

Funding Source: This research was supported by the Regional and Global Model Analysis (RGMA) component of the Earth and Environmental System Modeling (EESM) program of the U.S. Department of Energy's Office of Science, as a contribution to the HiLAT-RASM project. 

###Study Location: Yukon River Basin, Alaska, USA

#### Bounding Coordinates:
-Northernmost Latitude: 68.0°N
-Southernmost Latitude: 60.0°N
-Westernmost Longitude: -165.0°W
-Easternmost Longitude: -130.0°W

Geographic Coverage Description: Although no new field samples were collected directly by the authors for this work, all data used in this file pertain to the Yukon River Basin in Alaska, USA. The data were compiled from previously published literature and publicly available experimental portals, including observational and biogeochemical datasets relevant to this region. Although the data originated from various sources, they all describe environmental or geochemical conditions within the specified region.

## Github link:
https://github.com/ggibson-LANL/AROMM
Project Summary:



## Folder Organization
The model code is located in:

### Yukon_River_Model_Code

The main model folder contains four subdirectories, each of which includes model code designed to simulate the transformation of Dissolved Organic Carbon (DOC) components in the Yukon River. Each subdirectory serves a specific purpose:

There are detailed readme files in each subdirectory with additional information on file contents

#### /Original_Yukon_River_Model_Code_plot
Purpose:
 This directory contains the original version of the river model, which is used as a reference baseline.
Key Features:
Simulates DOC transformation along the Yukon River.
Generates plots of DOC concentration versus river distance.
Usage:
 Use this model when comparing new developments against the original implementation or when visual outputs are needed for DOC trends along the river.

#### /Original_Yukon_River_Model_Code_for_No_of_Simulations
Purpose:
 Designed to perform Monte Carlo simulations with user-defined numbers of iterations.
Key Features:
Outputs simulation results for varying numbers of Monte Carlo runs.
The number of simulations must be manually adjusted before each run.
Usage:
 Use this model when performing sensitivity analyses or uncertainty assessments requiring different simulation counts.

#### /Original_Yukon_River_Model_Code_Experiment_I
This model code creates outputs for a  fixed number of Monte Carlo simulations and initial values that are related to Experiment I
Purpose:
 Implements Experiment I, using a fixed number of Monte Carlo simulations and a specific set of initial conditions.
Key Features:
Tailored for reproducible runs of Experiment I.
Inputs and settings are predefined and consistent.
Usage:
 Use this model to generate results specific to Experiment I or when reproducing previously reported findings from that experiment.

#### /Original_Yukon_River_Model_Code_Experiment_II
This model code creates outputs for a  fixed number of Monte Carlo simulations and initial values that are related to Experiment II
Purpose:
 Implements Experiment II, similar in structure to Experiment I but using a different set of initial conditions.
Key Features:
Runs Monte Carlo simulations with fixed settings related to Experiment II.
Usage:
 Use this version for analyses specific to Experiment II or when comparing outcomes between experimental configurations.
Lineage
The Yukon River model code was developed from the Idealized_Arctic_River_Model, which has been archived for prosperity. 

### Idealized_Arctic_River_Model

#### /Idealized_Arctic_River_Model/Idealized_Yukon_Multi_Node_Model
This directory contains the model code for a multi-node, idealized river chemistry model, which represents the Yukon River as a series of interconnected spatial compartments.

#### /Idealized_Arctic_River_Model/One_Node_Model
This folder includes the single-node version of the idealized river model.

#### /Idealized_Arctic_River_Model/Two_Node_Model
Contains the two-node version of the idealized river chemistry model.

### How to run the model
#### 1. Navigate to the Model Code Folder
Start by navigating to the folder that contains the model scripts. The model is organized as a set of modular scripts, each serving a specific function in the simulation pipeline.

#### 2. Customize Initial Conditions
To set or modify the initial values for the model (e.g., initial DOC concentrations, reaction rate constants, flow conditions), and the no of simulations, use the script: 
river_initial_monac_code.py. This script defines the boundary and initial conditions used by the main simulation and should be edited before each new run if parameter changes are needed.
This script contains the main parameters for the model simulation: <br>
-Main stem and tributary distance
-No of tributary points
-Number of distance points for each tributary and main stem.
-Velocity values for the main river and the tributaries.
-Dilution factor values for confluence points.
-Initial DOC and macromolecule concentration values.
-Chemical turnover time values.
-Initial production values for macromolecules.
-This file will be called by each tributary file and the main Yukon River file to gain the necessary data.
-This script also performs the Monte Carlo Sensitivity Analysis to generate the random initial values for parameters within the boundaries we have gives for each run.

#### 3. Run the Main Simulation
After customizing the initial conditions, execute the main model script:
python Yukon_river.py    
This script calls all of the tributary scripts (i.e. Teslin_river.py, Pelly_river.py, White_Donjec_river.py, Porcupine_river.py, Tanana_river.py, Stewart_river.py, and Koyukuk_river.py) as well as the initial river condition (river_initial_monac_code.py) script
This script runs the simulation of DOC transport and transformation along the Yukon River, using the inputs defined in river_initial_monac_code.py

#### 4. Analyze and Store Model Outputs
To analyze simulation results and automatically save the output data, run the script:
Data_Ana_3.py
This script will:
This will run the model script (yukon_river.py)
Extract and analyze the output data from the simulation.
Generate summary statistics, plots, or visualizations as needed.
Save the processed output files for further review.

#### 5. Optional: Run Sub-basin Models
If needed, you can also simulate individual tributaries or sub-basins of the Yukon River by running any of the following scripts:
Teslin_river.py
Pelly_river.py
White_Donjec_river.py
Porcupine_river.py
Tanana_river.py
Stewart_river.py
Koyukuk_river.py

These can be useful for localized sensitivity analyses or when validating sub-regional hydrological and biogeochemical behavior.

### Workflow Summary:
Edit river_initial_monac_code.py to define your scenario<br>
-Run Yukon_river.py to execute the simulation<br>
-Analyze the output with Data_Ana_3.py<br>
-Optionally simulate individual tributaries with scripts like Teslin_river.py or Tanana_river.py<br>
-Customize the model using river_initial_monac_code.py, which defines the fundamental structure and initial conditions of the Yukon River system. This script must be configured prior to each simulation run and serves as the central place to define key model parameters.<br>

The script allows the user to specify:<br>
-Distances for the main stem and tributaries<br>
-Number of tributary points<br>
-Number of spatial distance points for each tributary and the main stem<br>
-Flow velocity values for the main river and each tributary<br>
-Dilution factors at tributary confluence points<br>
-Initial dissolved organic carbon (DOC) concentrations and macromolecule concentrations<br>
-Chemical turnover times of macromolecular compounds<br>

The script also sets the number of runs for the Monte Carlo simulation<br>
-Defines boundary ranges for sensitivity analysis<br>
-Performs Monte Carlo Sensitivity Analysis, generating randomized parameter values within the defined boundaries for each simulation run<br>
-Defining the River Structure<br>

The river structure defined in river_initial_monac_code.py includes:<br>
-Number of tributaries<br>
-Specifies how many major tributaries feed into the main stem (e.g., Pelly, Tanana, Koyukuk).<br>
-Distances<br>
-The total length of the main river stem, the length of each tributary, and the length to the each confluence point from the starting point of the main river stem.<br>

### Confluence points
Indicate the spatial location (along the main stem) where each tributary joins. The length These are critical for simulating dilution and mixing processes.

### Points of interest
Specific locations along the main stem or tributaries where output data may be recorded or where particular biogeochemical events (e.g., monitoring stations, hotspots) are expected to occur.
We use published literature to obtain the geographical coordinates of key locations in the Yukon River system, including the starting points of the main river and tributary headwaters, the confluence points where tributaries meet the Yukon main stem, and the endpoint of the Yukon main stem.

### Flow Velocity Data: 
The river velocities of the main Yukon stem and the input is customized in river_initial_monac_code.py, using observational data by calculating the average velocity separately for the main stem and the tributaries. 

### Dilution factor Data:
In river_initial_monac_code.py, the dilution factors at the river confluence points are customized using observational data of river fluxes from main stem and tributaries.

### Initial DOC Data:
In river_initial_monac_code.py, the initial DOC concentration at each tributary headwater is set based on the average soil organic carbon obtained from observational data, which is then converted to dissolved organic carbon through soil processing mechanisms.

### Chemical fraction:
The chemical fractions of DOC components in river_initial_monac_code.py are defined using information from the literature and expert judgment. These fractions establish the initial macromolecular composition of DOC at the headwaters.

### Chemical turnover time: 
The chemical turnover of macromolecular components in river_initial_monac_code.py are defined using information from the literature and expert judgment. These fractions establish the initial macromolecular composition of DOC at the headwaters.
​​
