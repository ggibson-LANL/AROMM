# AROMM
Arctic Riverine Organic Macromolecular Model 
Project title: Arctic Riverine Organic Macromolecular Model (AROMM)                                                                                                                      

Keywords: Yukon River, Alaska, Arctic Rivers, Subarctic Regions, Yukon Watersheds

Model Developer: Dr. Amadini Mendis Jayasinghe  
Affiliation*: Consultant Senior Lecturer, CINEC Campus, Malabe, Sri Lanka  
Email: amadini.mj@gmail.com
* Amadini Mendis Jayasinghe developed the model while serving as a Postdoctoral Research Associate at Los Alamos National Laboratory (LANL), CCS-2 Division, Los Alamos, NM, USA. 

Repository manager: Dr. Georgina Anne Gibson  
Affiliation: Scientist, CCS-2 Division, Los Alamos National Laboratory (LANL)  
Email: ggibson@lanl.gov 

Funding Source: This research was supported by the Regional and Global Model Analysis (RGMA) component of the Earth and Environmental System Modeling (EESM) program of the U.S. Department of Energy's Office of Science, as a contribution to the HiLAT-RASM project. 

Link to manuscript.  https://www.authorea.com/users/929544/articles/1300830-modeling-the-sensitivity-of-yukon-river-biogeochemical-dynamics-to-environmental-and-chemical-drivers-implications-for-dissolved-organic-carbon?commit=adc20a7ca10f4b81bb41d96f543d7390dedbe324
DOI: 10.22541/au.174861009.92610692/v1
Hydrological Processes

Abstract
Riverine dissolved organic carbon (DOC) is a critical biogeochemical component that transmits information from Arctic soils to the Arctic Ocean, significantly influencing carbon dynamics in this unique ecosystem. As DOC travels downstream, it undergoes transformations that alter its composition and fate. The Yukon River serves as an effective testbed for modeling these dynamics, offering sufficient scale to capture key biogeochemical processes, simpler hydrology than other major Arctic rivers, and access to long-term DOC data for model validation. To investigate DOC transformations during transit, we adapted our Arctic Riverine Organic Macromolecular Model by applying regional-specific parameterizations. Our model simulates the transport and transformation of 15 organic macromolecules, including CDOM (Coloured Dissolved Organic Matter), proteins, polysaccharides, lipids, lignin phenols, and humic substances. Initial DOC concentrations were derived from surrounding observed soil organic carbon stocks, while chemical transformations and hydrological dynamics were modeled along the river’s course. Sensitivity and uncertainty analyses were conducted using a Monte Carlo approach under two experimental setups. Results revealed that variability in DOC and CDOM concentrations at the river mouth were predominantly driven by initial DOC concentration (~70% of variability explained) and dilution at confluence points (~10%). The refractory fraction of DOC explained 21-88% of the variability in 14 macromolecular concentrations. River velocity, which determines residence time, explained 8-47% of the variability in protein, polysaccharide, lipid, pigments, and lignin phenols at the river mouth. Incontrast, chemical turnover times contributed only 1–5% to output variability. Our findings underscore the need for improved land-specific headwater observations, including seasonal soil moisture and lateral transport dynamics that control the initial tributary-specific DOC inputs. With accelerated permafrost thaw and increasing river discharge, extending our model to other Arctic River systems and seasons will enhance understanding of Arctic riverine carbon fluxes and their contributions to the Arctic Ocean.                                                                                                           

Study Location Name: 
Yukon River Basin, Alaska, USA
 
  Bounding Coordinates:
  Northernmost Latitude: 68.0°N
  Southernmost Latitude: 60.0°N
  Westernmost Longitude: -165.0°W
  Easternmost Longitude: -130.0°W
  Geographic Coverage Description: 
A
lthough no new field samples were collected directly by the authors for this work, all data used in this file pertain to the Yukon River Basin in Alaska, USA. The data were compiled from previously published literature and publicly available experimental portals, including observational and biogeochemical datasets relevant to this region. Even though the data came from other sources, they all describe environmental or geochemical conditions within this specified region.

Github link:
https://github.com/ggibson-LANL/AROMM

Project Summary:

This project focused on modeling and analyzing the transport and transformation of Dissolved Organic Carbon (DOC) and its macromolecular composition in Arctic rivers, with a particular emphasis on the Yukon River. A central goal was to develop a process-based river model capable of simulating the reactivity and fate of dissolved organic components under conditions characteristic of Arctic river systems. The specific objectives were to: (1) model the chemical reactivity of DOC components; (2) parameterize the model to reflect the physical and biogeochemical characteristics of the Yukon River; (3) estimate DOC and it’s macromolecular component concentrations at the river mouth; and (4) evaluate spatial and temporal changes in chemical concentrations along the river continuum.

The AROMM model developed under this project simulates the evolution of dissolved organic chemicals along the river, from the headwater to the river mouth. 

Model code is written in Python and runs on Mac OS. The model should be compatible with any operating system. Python ≥ 3.9 is required, along with the following Python packages:
Numpy - Numerical operations
Pandas - Data handling and input/output
Satplotlib - Plotting and visualization
Scipy - Statistical distributions and optimization
pyDOE -Design of experiments (Latin Hypercube Sampling, content, and the structure of the model)
