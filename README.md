# Surface Energy Balance and Firn Model v1.3

This repository contains the surface energy and mass balance as well as the firn model used for the study: 

Vandecrux, B.; Fausto, R. S.; Langen, P. L.; As, D. Van; MacFerrin, M.; Colgan, W.T.; Ingeman-Nielsen, T.; Steffen, K.; Jensen, N.S.; Møller, M. T.; Box, J.E.. Drivers of Firn Densification on the Greenland Ice Sheet Revealed by Weather Station Observations and Modelling. Journal of Geophysical Research. 2018. https://doi.org/10.1029/2017JF004597

It is available for free at https://github.com/BaptisteVandecrux/SEB_Firn_model
When using any piece of script, cite the article above and give the link to the repository in the acknowledgment.

If using specifically the Surface Energy Balance Model the following reference should be added:

van As, D., van den Broeke, M., Reijmer, C., & van de Wal, R. (2005). The summer surface energy balance of the high antartic plateau. Boundary Layer Meteorol., 115, 289-317. https://doi.org/10.1007/s10546-004-4631-1

If using specifically the Firn Model the following reference should be added:

Langen, P., Fausto, R. S., Vandecrux, B., Mottram, R. H., & Box, J. E. (2017). Liquid Water Flow and Retention on the Greenland Ice Sheet in the Regional Climate Model HIRHAM5: Local and Large-Scale Impacts. Front. Earth Sci., 4(110). https://doi.org/10.3389/feart.2016.00110

The repository contains all the scripts needed to run the model and produce the plots and data used in the Vandecrux et al. (2018). However I do not own the input files (station data, firn cores...) and therefore cannot put them online. If you need these input files, please contat me.

The ouput of the model for Vandecrux et al. (2018) are available here:

Baptiste Vandecrux. 2018. Surface energy balance and modeled firn density, Greenland ice sheet, 1998-2015. Arctic Data Center. doi:10.18739/A2TH8BM6X. 

Baptiste Vandecrux
b.vandecrux@gmail.com



This work is part of the Retain project funded by the Danish Council for Independent research (Grant no. 4002-00234). The GC-Net data is available at http://cires1.colorado.edu/steffen/gcnet. The data produced by this study (processed gap-filled hourly standardized weather data and SEB and firn model output) are available for download at https://arcticdata.io . The scripts used for this study are available at at https://github.com/BaptisteVandecrux/SEB_Firn_model. The PARCA cores are available at http://research.bpcrc.osu.edu/Icecore/data/ and their collection was supported by NASA grants NAG5-5032, 6817, 5031, 6779, NAGW-4248 and NSF/OPP grant 9423530. KAN_U weather station data is funded by the Greenland Analogue Project (GAP), and made available through the Programme for Monitoring of the Greenland Ice Sheet (PROMICE). PROMICE data is freely accessible at http://promice.org. HIRHAM5 output is available at http://prudence.dmi.dk/data/temp/RUM/HIRHAM/.
All the scripts were developed on Matlab 2015b.
Baptiste Vandecrux
b.vandecrux@gmail.com
Last update : 09/09/2017

Structure of the program:
The working folder should contain:
“Input” folder. It should contain the three foldes:
Constants. Containing different csv files where the value of the constants used by the model are defined. If the user wants to change these constant values, it can either be done in the csv files or in the main script of the model.
Initial state. Contains the initial density and temperature profiles used for the various stations. They can be generated using the InitialDensityProfileGenerator.m script.
Weather data. Contains all the weather station data. The weather data files that I used as input contain too much GCnet data to be distributed.  Yet I provide one year of my input files to give an outlook of the model.  Please contact me if you need the original input files. An alternative is to feed the model with the standardized weather data that is freely available on https://arcticdata.io (needs some formatting).
     “Output” folder
Where all output plots and documents will be placed. For each new run a new folder is created.
“lib” folder
Collection of various functions used by the model. Each of them should contain a short description.
“main.m” script
Main script for running the model. It is possible to modify some options within that script. The output of the model is saved in nc format in the output folder.
“Analysis.m” script
Generate the plots used in the publication.

Coding and running tips
I have divided every scripts into sections using the command ‘%%’. Matlab allow to expend or collapse each section using the ‘+’ or ‘-‘ sign next to ‘%%’. This works also for ‘for’ loops and ‘if’ loops. For better visibility, you can always start by collapsing all sections (‘View’ tab> collapse all) and open only the section that you are interested in. Each section should have a clear name to help you navigate into the script. This should avoid you to scroll up and down all day.
‘ctr + d’ when you have your cursor on a function will open the function script in Matlab. A description of the function is available at the beginning of every function script. When used on a variable, will display its content in the variable tab of Matlab.
‘F9’ when you have selected a command, one or multiple lines will run only what was selected.
‘ctr + enter’ will run the section in which your cursor is located (from previous ‘%%’ to next ‘%%’ command). This tip and the previous ones are useful for running your script step by step and debug it.
Stop points (red dot that you insert/remove by clicking on the left side of a line) makes Matlab stop every time it comes across that point. Useful to debug a loop or see what happens in a function.
I use Matlab structure objects, for example to simplify the input/output of functions or group some variables together. I also used time series objects and table objects. They all allow an to manipulate complex data in an easier way. Read the online documentation about these different objects to know how to use them.
‘ctr + c’ when the cursor is in the command prompt to interrupt whatever Matlab is doing.

Frequently asked questions:
What is the difference between variables labeled Origin e.g. AirPressurehPa and AirPressurehPa_Origin?
The weather data provided are a composite of various sources (see the JGR paper and figures in the supplementary material). Each weather variable has its own “Origin” field where:
0 = from main station;  1 = adjusted from CP2;  2 = adjusted from Swiss Camp; 3 = adjusted from KAN_U;4 = adjusted from HIRHAM5; 5 = calculated using MODIS albedo.

Which variable is surface meltwater production in m of w.e. ?
For the melt, the melt flux variable (W m-2) is named meltflux and is generated in SurfEnergyBudget.m. It is later on converted to melt amount (m w.eq.) in the variable dH_melt_weq in the function MassBudget. Later on, it is renamed snmel (name for melt in HIRHAM’s subsurface). The sign convention is that meltflux  and dH_melt_weq  are both negative (because they take energy and mass away from the surface) while  is positive (since it is an amount to be injected into the subsurface). It is then printed in the file “surface-bin1.nc” under the variable “H_melt”.






