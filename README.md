# dust_extinction_estimator
This is a simple GUI code built for estimating the extinction caused by dust in QSO-DLAs

GUI code to estimate the extinction caused by dust at the rest frame of the QSO as well as any absorber.

This is a GUI code that takes in a spectra from a source file and estimates the extinction caused by dust. The GUI features are currently available only for one rest frame absorber. The code further lets you choose absorption free data points in the source file that can be used for fitting for dust extinction. The code can be run by simply downloading all the files in the folder and compiling the 'main_GUI.py' through a linux terminal (if all required modules are available, see below). The code is linux terminal friendly and has not been checked for other operating systems. 

Required dependancies - This is a simple python code (version - 2.7). Due to issues in running matplotlib widget in version 3, unfortuantely the code might not work in python 3. However, experts are welcome to change the code and try. Here is the list of python modules required to run this code -

numpy

matplotlib

sys

itertools

os

scipy

pathlib

subprocess

tabulate



Other essential non python files include -

initial_parameters.dat - This is the primary parameter file and is required to give initial values for many paramters used by the GUI. The file is self descriptive. Please have a look to change the default.

data_file_J2140-0321.dat - This is an example of an acceptable format of the source file. It has to be a binary with wavelength, flux and flux error in that order. Any comments should start with '#'.

selsing_new.dat - This is an example of the optical quasar template spectra. The fitting function takes this spectra (flux calibrated) as the model for the spectra without reddening. The change in the instrinsic slope of this quasar template can be controlled by the delta-beta paramter (see below). Any other similar file can be used. The name of the template file should be updated in the GUI code, if changed. 



Here is a short description of interactive functions that one can perform and all the buttons that you see in the GUI -


Selecting data points for fitting purpose - The quasar spectra might be contaminated with all sorts of absorption from the metal lines and sky residuals. For the purpose of fitting the dust continuum, all these can be treated as contaminants. Hence, for fitting, I have created a point selection method. One has to specify the quasar continuum area that are free of absorption line contaminants before one starts the fit. To do the same, one has to press the buttons 't' -> 'y' -> 'h' in that sequence. This activates the selection function and now one can draw a rectangle around data points to select them for fitting. In case you would like to remove some data points that have been selected by mistake, you can activate the remove function by pressing 'j' (or 't' -> 'y' -> 'j', if the function is inactive). After this, the points inside any selected rectangle will disappear. Once all the relevant areas for the fit are selected, one can press the 'save_selec' button to save the list of selected data points.

All the matplotlib basic plotting buttons and keys work as their defaults.



All custom button functionalities (clockwise) -

Bar: 'Flux_reduction' - The quasar template has a specific continuum flux that may or may not be on the same scale as that of the quasar being analysed. This bar is used to scale the continuum up and down to match the flux of the quasar. We have used the template of a quasar for this code. However, the template can easily be changed to other type of bright source. The scaling done here goes as an initial guess to the fit, as it is one of the free parameters.

Radio button: 'Fit' - The button is for executing the fitting function. It will iteratively fit for all different extinction classes given - 'Galactic', 'SMC_Bar', 'LMC_Supershell', 'LMC_Average'. The best fit class is decided by the minimal chi-squared obtained and this solution (with all other parameters) are updated in the GUI.   


Radio button: 'Fit Selected' - The button is for executing the fitting function for the extinction classes specified in the GUI. If, say the GUI is set to 'LMC_Average' extinction class, the code will give best fit results with 'LMC_Average' extinction class. The solution (with all other parameters) are updated in the GUI.  


Radio buttons: 'save_selec/get_selec' - The 'save_selec' button can be used to save (in a file) all the data points chosen for performing the fit. The 'get_selec' button will retrieve the selected data point information from a file, if previously saved.


Radio buttons: 'save_params/get_params' - The 'save_params' button can be used to save (in a file) all the fitted or updated values for all paramters. The 'get_params' button will retrieve the best fitted parameters from a file, if previously saved.


Bar: 'Flux_smoothing_model' - This bar can be used to smooth the quasar model template to get rid of any small scale variations that might be understood as relevant features by the fit. The fitting function takes this smoothing of data into account before proceeding with the fit. 


Bar: 'Flux_smoothing_data' - This bar can be used to smooth the source data. The fitting function does not take this smoothing of data into account before proceeding with the fit.


Bar: 'Delta_Beta' - This bar controls the value of Delta-beta, which is a parameter that controls the power-law slope relative to the intrinsic shape of the quasar. This bar will be used to obtain the initial guess on delta-beta before starting the fit.


Bar: 'E(B-V) DLA' - This bar controls the value of reddening, E(B-V), at the redshift of a potential absorber. This bar will be used to obtain the initial guess on E(B-V) at the absorber redshift before starting the fit.


Bar: 'E(B-V) QSO' - This bar controls the value of reddening, E(B-V), at the redshift of the quasar. This bar will be used to obtain the initial guess on E(B-V) at the quasar redshift before starting the fit.


Radio buttons: 'Galactic', 'SMC_Bar', 'LMC_Supershell' and 'LMC_Average' - The radio buttons are used to select an extinction class for the dust extinction at the absorber redshift. The dust extinction at the quasar redshift is fixed to 'SMC_Bar' based on the previous studies and observations.


Check buttons: 'Reddened_Flux', 'Flux' - The check buttons are used to display the reddened and unreddened version of the quasar template in the GUI screen.


Button: 'Save Fig' - Saves the dust extinction diagram with the updated extinction parameters in a file.


Button: 'Plot' - Plots the dust extinction diagram with the updated extinction parameters on screen.


Button: 'print/load_latex' - Loads the previously saved fit parameters from a file and also prints it in the terminal screen.


Button: 'print/save_latex' - Saves the fit parameters from the GUI screen and also prints it in the terminal screen.


Bar: 'DLA-Redshift' - Bar that can be used to vary the absorption redshift

Bar: 'QSO-Redshift' - Bar that can be used to vary the QSO emission redshift




The code is free to use for all. Thank you for using the code and please feel free to contact me at - 'ranjan_adarsh@yahoo.com' for any comments, suggestions and bugs.
