# The GEUS Surface Energy Balance and Firn Model
Baptiste Vandecrux
bav@geus.dk
All the scripts were developed on Matlab 2015b.

# Table of Contents
1. [Overview](#overview)
2. [Structure of the program](#structure-of-the-program)
3. [Getting started](#start)
4. [Coding and running tips](#tips)
5. [Frequently asked questions](#faq)
6. [Supported studies and references](#ref)
7. [Funding](#funding)

## Overview
These are the scripts for the Surface Energy Balance model from [van As et al. (2005)](https://doi.org/10.1007/s10546-004-4631-1) as used in Vandecrux et al. ([2018](https://doi.org/10.1029/2017JF004597), [2020a](https://doi.org/10.1017/jog.2020.30), [2020b](https://doi.org/10.5194/tc-2019-331)) and from the GEUS firn model also described in Vandecrux et al. ([2018](https://doi.org/10.1029/2017JF004597), [2020a](https://doi.org/10.1017/jog.2020.30), [2020b](https://doi.org/10.5194/tc-2019-331)). They take as input gap-free time series of air temperature, humidity, wind speed, downward and upward shortwave radiation and downward longwave radiation, along with instrument heights and snowfall. These data can be prepared at PROMICE and GC-Net sites using this [AWS data processing suite](https://github.com/BaptisteVandecrux/AWS_processing). 

The surface energy balance model finds iteratively the subfreezing surface temperature that closes the energy budget defined as the sum of:
- upward shortwave radiation (given by input file)
- downward shortwave radiation (given by input file)
- downward longwave radiation (given by input file)
- upward shortwave radiation (calculated from surface temperature)
- latent and sensible heat fluxes (calculated from surface temperature and input meteorology using the [Monin-Obukhov similarity theory](https://en.wikipedia.org/wiki/Monin%E2%80%93Obukhov_similarity_theory))
- conductive heat flux to the underlying snow and firn (calculated from the surface temperature and the firn model)

If it is not possible to find a subfreezing surface temperature that nullify the sum of energy fluxes, then surface temperature is set to 0°C and the sum of all energy fluxes is then used to melt surface material, snow or ice.

The GEUS firn model is a multilayer snow and firn model that tracks the temperature, density, grain size and water content for each model layer. At each time step, the model column is updated for temperature diffusion, firn compaction, grain growth and meltwater infiltration. More information about each of these routines can be found in Vandecrux et al. ([2018](https://doi.org/10.1029/2017JF004597), [2020a](https://doi.org/10.1017/jog.2020.30), [2020b](https://doi.org/10.5194/tc-2019-331)) and in the references therein.

## Structure of the program:
The working directory should contain:
* “Input” folder. It should contain the three folders:
    * Constants. Containing different csv files where the value of the constants used by the model are defined. If the user wants to change these constant values, it can either be done in the csv files or in the main script of the model.
    * Initial state. Contains the initial density and temperature profiles used for the various stations. They can be generated using the InitialDensityProfileGenerator.m script.
    * Weather data. Contains all the weather station data. The weather data files that I used as input contain too much GCnet data to be distributed.  Yet I provide one year of my input files to give an outlook of the model.  Please contact me if you need the original input files. An alternative is to feed the model with the standardized weather data that is freely available on https://arcticdata.io (needs some formatting).
* “Output” folder: Where all output plots and documents will be placed. For each new run a new folder is created.
* “lib” folder: Collection of various functions used by the model. Each of them should contain a short description.
* “main.m” script: Main script for running the model. It is possible to modify some options within that script. The output of the model is saved in nc format in the output folder.
* “Analysis.m” script: Generates the plots used in the publication.

## Getting started <a name="start" />

- Download the scripts from Github
- Place the input data file in the Input/Weather data folder. For the example, let's continue with the data_KAN_M_combined_hour.txt previously generated using the [AWS data treatment suite](https://github.com/BaptisteVandecrux/AWS_processing). A short version of that file is already present in the Weather data folder.

- Open Main.m in Matlab and set the name of the station you want to process:
```
station_list =   {'KAN_M'};
```
and RCM that was used for gap-filling (only used for naming of folders)
```
RCM_list = {'RACMO'};
```
add the station and its location to the list of station data location (~line 85):
```
        case 'KAN_M'
            param{kk}.InputAWSFile = './Input/Weather data/data_KAN_M_combined_hour.txt';
``` 

- We now need to prescribe the initial density and temperature profile. At KAN_M the surface is composed of bare ice seasonally covered by snow. In the folder Input/Initial state/density, there is no file for KAN_M so we can make a copy of *DensityProfile_NUK_K.csv* (NUK_K is also composed of bare ice overlaid by snow) and name it *DensityProfile_KAN_M.csv*. You can see that there are some density profiles in that folder that can be used to initialize the model. More info about how the model initialize density, temperature and liquid water content are available in [InitializationSubsurface.m](/lib/HHsubsurf/Initialization/InitializationSubsurface.m).

- Make sure that the station is listed in [StationInfo.csv](/Input/Constants/StationInfo.csv) along with the appropriate information about the locaiton. Default values are currently used for KAN_M and updated info should be inserted in that file.

- You can now run the Main.m script. It will create a folder in *Output* where outputs will be saved in netcdf format. Consider the [Panoply tool](https://www.giss.nasa.gov/tools/panoply/) for easy visualization.

- Many ideas for plots and analysis are available in [Analysis.m](Analysis.m) or in the [Study-specific scripts](/Study-specific%20scripts)

## Coding and running tips <a name="tips" />
- I have divided every scripts into sections using the command ‘%%’. Matlab allow to expend or collapse each section using the ‘+’ or ‘-‘ sign next to ‘%%’. This works also for ‘for’ loops and ‘if’ loops. For better visibility, you can always start by collapsing all sections (‘View’ tab> collapse all) and open only the section that you are interested in. Each section should have a clear name to help you navigate into the script. This should avoid you to scroll up and down all day.
- ‘ctr + d’ when you have your cursor on a function will open the function script in Matlab. A description of the function is available at the beginning of every function script. When used on a variable, will display its content in the variable tab of Matlab.
- ‘F9’ when you have selected a command, one or multiple lines will run only what was selected.
- ‘ctr + enter’ will run the section in which your cursor is located (from previous ‘%%’ to next ‘%%’ command). This tip and the previous ones are useful for running your script step by step and debug it.
- Stop points (red dot that you insert/remove by clicking on the left side of a line) makes Matlab stop every time it comes across that point. Useful to debug a loop or see what happens in a function.
- I use Matlab structure objects, for example to simplify the input/output of functions or group some variables together. I also used time series objects and table objects. They all allow manipulating complex data in an easier way. Read the online documentation about these different objects to know how to use them.
- ‘ctr + c’ when the cursor is in the command prompt to interrupt whatever Matlab is doing.

## Frequently asked questions:<a name="faq" />
* What is the difference between variables labeled Origin e.g. AirPressurehPa and AirPressurehPa_Origin?
The weather data provided are a composite of various sources (see the JGR paper and figures in the supplementary material). Each weather variable has its own “Origin” field where:
0 = from main station;  1 = adjusted from CP2;  2 = adjusted from Swiss Camp; 3 = adjusted from KAN_U;4 = adjusted from HIRHAM5; 5 = calculated using MODIS albedo.

* Which variable is surface meltwater production in m of w.e. ?
For the melt, the melt flux variable (W m-2) is named meltflux and is generated in SurfEnergyBudget.m. It is later on converted to melt amount (m w.eq.) in the variable dH_melt_weq in the function MassBudget. Later on, it is renamed snmel (name for melt in HIRHAM’s subsurface). The sign convention is that meltflux  and dH_melt_weq  are both negative (because they take energy and mass away from the surface) while  is positive (since it is an amount to be injected into the subsurface). It is then printed in the file “surface-bin1.nc” under the variable “H_melt”.

## Supported studies and references<a name="ref" />
This repository contains the surface energy and mass balance as well as the firn model used for the following studies. However I do not own the input files (station data, firn cores...) and therefore cannot put them online. If you need these input files, please contact me.

1. [Vandecrux et al. (2018) Drivers of Firn Densification on the Greenland Ice Sheet Revealed by Weather Station Observations and Modelling.](https://doi.org/10.1029/2017JF004597)
In that study we process data from four weather stations in the accumulation area of the Greenland ice sheet, run and evaluate the firn model and use the model output to document the processes behind the seasonal and long term changes in firn density. The scripts, as they were at the time of the sutdy are available [here](https://github.com/BaptisteVandecrux/SEB_Firn_model/releases/tag/0.1) and the data is available for free [here](https://arcticdata.io/catalog/view/doi:10.18739/A2TH8BM6X). 

2. [Vanderux et al. (2020a) Firn cold content evolution at nine sites on the Greenland ice sheet between 1998 and 2017](https://doi.org/10.1017/jog.2020.30). In that study, we use the firn model at nine weather station sites to document the seasonal and long term changes in the firn's cold content, which controls the firn's ability to retain meltwater. The data from this study is available [here](https://arcticdata.io/catalog/view/doi%3A10.18739%2FA2QF8JJ9B), [here](https://arcticdata.io/catalog/view/doi%3A10.18739%2FA2V698C33), and [there](https://arcticdata.io/catalog/view/doi%3A10.18739%2FA2GX44V1P).

3. [Vandecrux et al. (2020b) The firn meltwater Retention Model Intercomparison Project (RetMIP): Evaluation of nine firn models at four weather station sites on the Greenland ice sheet](https://doi.org/10.5194/tc-2019-331). In that study, we compare the GEUS firn model to other firn models currently used in Greenland and evaluate its output against a set of firn observations. All model runs will be made available on the [PROMICE data portal](https://www.promice.org/PromiceDataPortal/) and the plotting scripts are available [here](https://github.com/BaptisteVandecrux/RetMIP).

4. [Vandecrux et al. (2021) Firn evolution at Camp Century, Greenland: 1966-2100](https://doi.org/10.3389/feart.2021.578978). In that study we use the GEUS model to evaluate the probability of meltwater interacting with military waste buried in the firn given three different climate scenarios.

When using any piece of script, cite the articles above and give the link to the repository in the acknowledgment.

If using specifically the Surface Energy Balance Model the following reference should be added:

van As, D., van den Broeke, M., Reijmer, C., & van de Wal, R. (2005). [The summer surface energy balance of the high antartic plateau](https://doi.org/10.1007/s10546-004-4631-1). Boundary Layer Meteorol., 115, 289-317. 

If using specifically the Firn Model the following reference should be added:

Langen, P., Fausto, R. S., Vandecrux, B., Mottram, R. H., & Box, J. E. (2017). [Liquid Water Flow and Retention on the Greenland Ice Sheet in the Regional Climate Model HIRHAM5: Local and Large-Scale Impacts](https://doi.org/10.3389/feart.2016.00110). Front. Earth Sci., 4(110). 

## Funding <a name="funding" />

This work is part of the Retain project funded by the Danish Council for Independent research (Grant no. 4002-00234). The GC-Net data is available at http://cires1.colorado.edu/steffen/gcnet. The data produced by this study (processed gap-filled hourly standardized weather data and SEB and firn model output) are available for download at https://arcticdata.io . The scripts used for this study are available at at https://github.com/BaptisteVandecrux/SEB_Firn_model. The PARCA cores are available at http://research.bpcrc.osu.edu/Icecore/data/ and their collection was supported by NASA grants NAG5-5032, 6817, 5031, 6779, NAGW-4248 and NSF/OPP grant 9423530. KAN_U weather station data is funded by the Greenland Analogue Project (GAP), and made available through the Programme for Monitoring of the Greenland Ice Sheet (PROMICE). PROMICE data is freely accessible at http://promice.org. HIRHAM5 output is available at http://prudence.dmi.dk/data/temp/RUM/HIRHAM/.
