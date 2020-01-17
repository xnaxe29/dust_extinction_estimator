
#Imports

import numpy
import numpy as np
import matplotlib.pyplot as plt
import sys
import itertools
import os.path
import os
import matplotlib.patches as mpatches
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from pathlib import Path
from matplotlib.widgets import Slider, Button, RadioButtons
from scipy.signal import savgol_filter
from matplotlib.widgets import CheckButtons
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import fmin
import subprocess
import os
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons, RectangleSelector
from tabulate import tabulate
import matplotlib.ticker as ticker
import matplotlib






#Getting source filenames and initial guess of the required paramters 

initial_guesses = np.genfromtxt('initial_parameters.dat', dtype=[('mystring','S30')], comments='#')

z = float(initial_guesses[0][0])
file_name1 = str(initial_guesses[1][0])
z_qso = float(initial_guesses[2][0])
size_of_font = float(initial_guesses[3][0])
quasar_template_file = str(initial_guesses[4][0])
quasar_template_smooth_factor = int(initial_guesses[5][0])
E_bv_dla = float(initial_guesses[6][0])
flux_red = float(initial_guesses[7][0])
extinction_class_string_initial = str(initial_guesses[8][0])
E_bv_qso = float(initial_guesses[9][0])
del_beta = float(initial_guesses[10][0])

font = {'family' : 'times new roman',
        #'weight' : 'bold',
        'size'   : size_of_font}

matplotlib.rc('font', **font)





########################################FLIP########################################

def flip(m, axis):
    if not hasattr(m, 'ndim'):
        m = asarray(m)
    indexer = [slice(None)] * m.ndim
    try:
        indexer[axis] = slice(None, None, -1)
    except IndexError:
        raise ValueError("axis=%i is invalid for the %i-dimensional input array"
                         % (axis, m.ndim))
    return m[tuple(indexer)]

########################################FLIP########################################


####################ARRAY_SMOOTHING_FUNCTION####################
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
####################ARRAY_SMOOTHING_FUNCTION####################







#######################################DEFINING_FIXED_PARAMETERS#######################################
#Note: The values here are taken from Gordon et al. 2003. Please feel free to change them if required
 
def extinction(extinction_class):
	if (extinction_class=='Galactic'):
		x0 = 4.592; gamma = 0.922; c1 = -0.175;	c2 = 0.807; c3 = 2.991; c4 = 0.319; c5 = 6.097; O1 = 2.055; O2 = 1.322; O3 = 0.000
		R_v = 3.001; k_IR = 1.057; x0_err = 0.00; gamma_err = 0.00; c1_err = 0.00; c2_err = 0.00; c3_err = 0.00; c4_err = 0.00
		c5_err = 0.00; O1_err = 0.00; O2_err = 0.00; O3_err = 0.00; R_v_err = 0.00; k_IR_err = 0.00

	elif (extinction_class=='SMC_Bar'):
		x0 = 4.600; gamma = 1.000; c1 = -4.959; c2 = 2.264; c3 = 0.389; c4 = 0.461; c5 = 6.097; O1 = 2.055; O2 = 1.322; O3 = 0.000
		R_v = 2.74; k_IR = 1.057; x0_err = 0.00; gamma_err = 0.00; c1_err = 0.197; c2_err = 0.040; c3_err = 0.110; c4_err = 0.079
		c5_err = 0.00; O1_err = 0.00; O2_err = 0.00; O3_err = 0.00; R_v_err = 0.00; k_IR_err = 0.00

	elif (extinction_class=='SMC_Wing'):
		x0 = 4.703; gamma = 1.212; c1 = -0.856; c2 = 1.038; c3 = 3.215;	c4 = 0.107; c5 = 6.097; O1 = 2.055; O2 = 1.322; O3 = 0.000
		R_v = 2.05; k_IR = 1.057; x0_err = 0.018; gamma_err = 0.019; c1_err = 0.246; c2_err = 0.074; c3_err = 0.439; c4_err = 0.038
		c5_err = 0.00; O1_err = 0.00; O2_err = 0.00; O3_err = 0.00; R_v_err = 0.00; k_IR_err = 0.00

	elif (extinction_class=='LMC_Supershell'):
		x0 = 4.558; gamma = 0.945; c1 = -1.475; c2 = 1.132; c3 = 1.463; c4 = 0.294; c5 = 6.097; O1 = 2.055; O2 = 1.322; O3 = 0.000
		R_v = 2.76; k_IR = 1.057; x0_err = 0.021; gamma_err = 0.026; c1_err = 0.152; c2_err = 0.029; c3_err = 0.121; c4_err = 0.057
		c5_err = 0.00; O1_err = 0.00; O2_err = 0.00; O3_err = 0.00; R_v_err = 0.00; k_IR_err = 0.00

	elif (extinction_class=='LMC_Average'):
		x0 = 4.579; gamma = 0.934; c1 = -0.890; c2 = 0.998; c3 = 2.719; c4 = 0.400; c5 = 6.097; O1 = 2.055; O2 = 1.322; O3 = 0.000
		R_v = 3.41; k_IR = 1.057; x0_err = 0.007; gamma_err = 0.016; c1_err = 0.142; c2_err = 0.027; c3_err = 0.137; c4_err = 0.036
		c5_err = 0.00; O1_err = 0.00; O2_err = 0.00; O3_err = 0.00; R_v_err = 0.00; k_IR_err = 0.00


	return (x0, gamma, c1, c2, c3, c4, c5, O1, O2, O3, R_v, k_IR)

#######################################DEFINING_FIXED_PARAMETERS#######################################


styles = mpatches.ArrowStyle.get_styles()
c = 299792.458 # Speed in Light in Km/s
k = 1.38064852e-23 #Boltzmann's constant in m^2 kg s^(-2) K^(-1)


#Data and file names
file_sel = str(file_name1)[:-4] + "_selected_data.txt"
file9 = str(file_name1)[:-4] + "_fit_parameters.txt"
file10 = str(file_name1)[:-4] + "_fit_parameters_latex.txt"


wave, flux, flux_err = np.loadtxt(str(file_name1), comments='#', unpack=True)




#MODEL_SPECTRUM
trish = np.loadtxt(str(quasar_template_file), comments='#')
x_trish = trish[:,0]*(1+z_qso)
y_trish = trish[:,1]

if (min(wave) < min(x_trish)):
	print ("Lower end (Shorter Wavelength) of the Model does not cover the lower end of the Data. Hence using simple extension.")
	x_trish = np.insert(x_trish, 0, min(wave))
	y_trish = np.insert(y_trish, 0, y_trish[0])

if (max(wave) > max(x_trish)):
	print ("Higher end (Longer Wavelength) of the Model does not cover the higher end of the Data. Hence using simple extension.")
	x_trish = np.append(x_trish, max(wave))
	y_trish = np.append(y_trish, y_trish[-1])





#SMOOTHING BY A FACTOR OF 20 (IDEAL FOR SELSING SPECTRUM (selsing_new.dat)). Please change if required.
y_trish = smooth(y_trish, quasar_template_smooth_factor)




#SYNCHRONISING SELSING MODEL AND DATA
idx_lower = np.searchsorted(wave, min(wave), side='left')
idx_upper = np.searchsorted(wave, max(wave), side='right')

f2 = interpolate.interp1d(x_trish[idx_lower:idx_upper], y_trish[idx_lower:idx_upper])
ynew = f2(wave)


#'SMC Wing' is not ideal for usage in dust extinction modelling. So removed here. Taking the rest 4 templates as defined in Gordon et al. 2003
raj = np.array(['Galactic', 'SMC_Bar', 'LMC_Supershell', 'LMC_Average'])










####################################REDDENING_FUNCTION####################################

def reddening_func2(F_rest_lambda, wave, E_bv_dla, R_v_dla, initial_result_dla, E_bv_qso, R_v_qso, initial_result_qso, del_beta):
	constant = np.zeros([len(F_rest_lambda)])
	F_lambda = np.zeros([len(F_rest_lambda)])
	for i in range(len(F_rest_lambda)):
		constant[i] = -0.4 * ((E_bv_dla * (initial_result_dla[i] + R_v_dla)) + (E_bv_qso * (initial_result_qso[i] + R_v_qso)))
		F_lambda[i] = ((F_rest_lambda[i]*((wave[i]/5510)**(del_beta))) * (10**(constant[i])))

	return (F_lambda)


def k_lambda_V(x, x0, gamma, c1, c2, c3, c4, c5):
	D_func = np.zeros([len(x)])
	result = np.zeros([len(x)])
	for i in range(len(x)):
		D_func[i] = x[i]**2 / ( ( (x[i]**2) - (x0**2) )**2 + (x[i] * gamma)**2 )
		if (x[i] <= c5):
			result[i] = c1 + (c2*x[i]) + c3*D_func[i]
		else:
			result[i] = c1 + (c2*x[i]) + c3*D_func[i] + c4*((x[i]-c5)**2)
		
	return (result)


#GLOBAL PARAMETER DEFINITIONS
U1_pos = 0.27; U2_pos = 0.26; O1_pos = 0.33; O2_pos = 0.4; O3_pos = 0.553; O4_pos = 0.7; I1_pos = 0; I2_pos = 1/0.25; I3_pos = 1/0.50; I4_pos = 1/0.75; I5_pos = 1/1.00

####################################REDDENING_FUNCTION####################################


#########################################REDDENING_BACKBONE_DLA#########################################


def spline_part_new_fit(extinction_class_string, c3_fit_new):
	x0, gamma, c1, c2, c3, c4, c5, O1, O2, O3, R_v, k_IR = extinction(extinction_class_string)
	c3=c3_fit_new

	I1_pos_new = 0.00000000001
	I2_pos_new = 1/I2_pos	
	I3_pos_new = 1/I3_pos
	I4_pos_new = 1/I4_pos
	I5_pos_new = 1/I5_pos

	O1_val = O1
	O2_val = O2
	O3_val = O3


	def I_n(k_IR, pos, R_v):
		return (((k_IR * (pos**(-1.84))) - R_v))


	I1_val = I_n(k_IR, I1_pos_new, R_v)
	I2_val = I_n(k_IR, I2_pos_new, R_v)
	I3_val = I_n(k_IR, I3_pos_new, R_v)
	I4_val = I_n(k_IR, I4_pos_new, R_v)
	I5_val = I_n(k_IR, I5_pos_new, R_v)

	x_new = np.array([1/O1_pos, 1/O2_pos, 1/O3_pos, 1/I4_pos, 1/I3_pos, 1/I2_pos, 0])
	y_new = np.array([O1_val, O2_val, O3_val, I4_val, I3_val, I2_val, I1_val])

	func = interpolate.interp1d(x_new, y_new)
	
	wave_new = (wave*1e-4)/(1+z)
	x = np.sort(1/(wave_new))
	x1 = np.array([])
	x2 = np.array([])
	for i in range(len(x)):
		if (x[i] < (1/O1_pos)):
			x1 = np.append(x1, x[i])
		else:
			x2 = np.append(x2, x[i])


	y1 = func(x1)
	position = np.searchsorted(x, (1/O1_pos))
	y = np.zeros([len(x)])

	y2 = k_lambda_V(x2, x0, gamma, c1, c2, c3, c4, c5)
	for i in range(len(x)):
		if (x[i] < (1/O1_pos)):
			y[i] = y1[i]
		else:
			y[i] = y2[i-position]


	y_k_lambda_V = y
	y_k_lambda_V = flip(y_k_lambda_V, 0)

	return (y_k_lambda_V)	


#########################################REDDENING_BACKBONE_DLA#########################################


#########################################REDDENING_BACKBONE_QSO#########################################

def spline_part_new_fit_qso():
	extinction_class_string = 'SMC_Bar'
	x0, gamma, c1, c2, c3, c4, c5, O1, O2, O3, R_v, k_IR = extinction(extinction_class_string)

	I1_pos_new = 0.00000000001
	I2_pos_new = 1/I2_pos	
	I3_pos_new = 1/I3_pos
	I4_pos_new = 1/I4_pos
	I5_pos_new = 1/I5_pos

	O1_val = O1
	O2_val = O2
	O3_val = O3

	def I_n(k_IR, pos, R_v):
		return (((k_IR * (pos**(-1.84))) - R_v))

	I1_val = I_n(k_IR, I1_pos_new, R_v)
	I2_val = I_n(k_IR, I2_pos_new, R_v)
	I3_val = I_n(k_IR, I3_pos_new, R_v)
	I4_val = I_n(k_IR, I4_pos_new, R_v)
	I5_val = I_n(k_IR, I5_pos_new, R_v)

	x_new = np.array([1/O1_pos, 1/O2_pos, 1/O3_pos, 1/I4_pos, 1/I3_pos, 1/I2_pos, 0])
	y_new = np.array([O1_val, O2_val, O3_val, I4_val, I3_val, I2_val, I1_val])

	func = interpolate.interp1d(x_new, y_new)
	
	wave_new = (wave*1e-4)/(1+z)
	x = np.sort(1/(wave_new))
	x1 = np.array([])
	x2 = np.array([])
	for i in range(len(x)):
		if (x[i] < (1/O1_pos)):
			x1 = np.append(x1, x[i])
		else:
			x2 = np.append(x2, x[i])


	y1 = func(x1)
	position = np.searchsorted(x, (1/O1_pos))
	y = np.zeros([len(x)])

	y2 = k_lambda_V(x2, x0, gamma, c1, c2, c3, c4, c5)
	for i in range(len(x)):
		if (x[i] < (1/O1_pos)):
			y[i] = y1[i]
		else:
			y[i] = y2[i-position]


	y_k_lambda_V = y
	y_k_lambda_V = flip(y_k_lambda_V, 0)

	return (y_k_lambda_V)	


#########################################REDDENING_BACKBONE_QSO#########################################












###################################################################################################
###################################################################################################
##############################THIS_IS_THE_DYNAMIC_PART_OF_THE_CODE#################################
###################################################################################################
###################################################################################################

###############################VALUE_AND_PLOT_INITIALISATION################################

E_bv = E_bv_dla
x0_dla, gamma_dla, c1_dla, c2_dla, c3_dla, c4_dla, c5_dla, O1_dla, O2_dla, O3_dla, R_v_dla, k_IR_dla = extinction(extinction_class_string_initial)

#Default extinction class for quasars set to - 'SMC_Bar'
x0_qso, gamma_qso, c1_qso, c2_qso, c3_qso, c4_qso, c5_qso, O1_qso, O2_qso, O3_qso, R_v_qso, k_IR_qso = extinction('SMC_Bar')

c3_fit_new = c3_dla
y_k_lambda_V = spline_part_new_fit(extinction_class_string_initial, c3_fit_new)
y_k_lambda_V_dla = spline_part_new_fit(extinction_class_string_initial, c3_fit_new)
y_k_lambda_V_qso = spline_part_new_fit_qso()
y_trish_new = reddening_func2(ynew, wave, E_bv_dla, R_v_dla, y_k_lambda_V_dla, E_bv_qso, R_v_qso, y_k_lambda_V_qso, del_beta)











#Making the matplotlib window and plotting initial lines
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax2 = fig.add_subplot(111)

line2, = ax2.plot([], [], 'g.')
line1, = ax.plot(wave,flux, 'k-', alpha=0.4, drawstyle='steps')
line5, = ax.plot(wave, (y_trish_new/flux_red), 'r-', alpha=0.6, visible=True)
line6, = ax.plot(x_trish, (y_trish/flux_red), 'b-', alpha=0.5, visible=True)






#Creating Objects for the Matplotlib Window
axcolor = 'lightgoldenrodyellow'
axqso = plt.axes([0.065, 0.96, 0.3, 0.02])
sqso = Slider(axqso, 'QSO-Redshift', (z_qso-1), (z_qso+1), valinit=z_qso)

axdla = plt.axes([0.065, 0.92, 0.3, 0.02])
sdla = Slider(axdla, 'DLA-Redshift', (z-1), (z+1), valinit=z)

axflux_red = plt.axes([0.60, 0.95, 0.3, 0.02])
sflux_red = Slider(axflux_red, 'Flux_reduction', (1e-7), (flux_red+20), valinit=flux_red)

axdel_beta = plt.axes([0.085, 0.005, 0.35, 0.015])
sdel_beta = Slider(axdel_beta, 'Delta Beta', (del_beta-1), (del_beta+1), valinit=del_beta)

axred_func = plt.axes([0.085, 0.030, 0.35, 0.015])
sred_func = Slider(axred_func, 'E(B-V) DLA', (0.0), (E_bv_dla+1), valinit=E_bv_dla)

axred_qso_func = plt.axes([0.085, 0.055, 0.35, 0.015])
sred_qso_func = Slider(axred_qso_func, 'E(B-V) QSO', (0.0), (E_bv_qso+1), valinit=E_bv_qso)




#Setting smoothing buttons for data and quasar template
global data_smooth
data_smooth = 1.0
axdata_smooth = plt.axes([0.86, 0.04, 0.1, 0.02])
sdata_smooth = Slider(axdata_smooth, 'Flux_smoothing_model', (data_smooth+0), (data_smooth+101), valinit=data_smooth+19)

global data_smooth2
data_smooth2 = 1.0
axdata_smooth2 = plt.axes([0.86, 0.01, 0.1, 0.02])
sdata_smooth2 = Slider(axdata_smooth2, 'Flux_smoothing_data', (data_smooth2+0), (data_smooth2+101), valinit=data_smooth2)



#Setting check buttons to display redenned and original lines
rax = plt.axes([0.005, 0.4, 0.1, 0.15])
check = CheckButtons(rax, ('Reddened_Flux', 'Flux'), (True, True))

# 1 - Tick Marks
def func(label):
	if label == 'Reddened_Flux':
		line5.set_visible(not line5.get_visible())
	elif label == 'Flux':
		line6.set_visible(not line6.get_visible())
	
	fig.canvas.draw()
check.on_clicked(func)

###############################VALUE_AND_PLOT_INITIALISATION################################






############################################SELECTOR############################################
#Used to selecting data points that would be used for fitting 



xdata = wave
ydata = flux
ydata_err = flux_err

text = ax.text(0, 0, "")
xnew_sel = wave
ynew_sel = np.full([len(flux)], np.nan)
ynew_sel_err = np.full([len(flux_err)], np.nan)
str1 = ''

def line_select_callback(eclick, erelease):
	global str1
	global xnew_sel, ynew_sel, ynew_sel_err
	x1, y1 = eclick.xdata, eclick.ydata
	x2, y2 = erelease.xdata, erelease.ydata
	idx5 = np.searchsorted(xdata, min(x1, x2))
	idx6 = np.searchsorted(xdata, max(x1, x2))
	min_y = min(y1, y2)
	max_y = max(y1, y2)

	if (str1 == 'add'):
		for i in range(idx5, idx6):
        		if (min_y < ydata[i] < max_y):
                		xnew_sel[i] = xdata[i]
                		ynew_sel[i] = ydata[i]
                		ynew_sel_err[i] = ydata_err[i]


	elif (str1 == 'rem'):
        	for i in range(idx5, idx6):
        		for j in range(len(xnew)):
                		if (min_y < ydata[i] < max_y) and (xdata[i] == xnew_sel[j]) and (ydata[i] == ynew_sel[j]):
                    			xnew_sel[j] = np.NAN
                    			ynew_sel[j] = np.NAN
                    			ynew_sel_err[j] = np.NAN


	else:
        	print ("Function Inactive.....")


	line2.set_data(xnew_sel, ynew_sel)
	fig.canvas.draw()


def toggle_selector(event):
	global str1
	print(' Key pressed.')
	if event.key in ['T', 't'] and toggle_selector.RS.active:
        	print(' RectangleSelector deactivated.')
        	toggle_selector.RS.set_active(False)
	if event.key in ['Y', 'y'] and not toggle_selector.RS.active:
        	print(' RectangleSelector activated.')
        	toggle_selector.RS.set_active(True)
	if event.key in ['H', 'h'] and toggle_selector.RS.active:
        	print('Add function activated')
        	str1 = 'add'
        	toggle_selector.RS.set_active(True)
	if event.key in ['J', 'j'] and toggle_selector.RS.active:
        	print('Remove function activated')
        	str1 = 'rem'
        	toggle_selector.RS.set_active(True)


toggle_selector.RS = RectangleSelector(ax2, line_select_callback, drawtype='box', useblit=False, button=[1], minspanx=5, minspany=5, spancoords='pixels', interactive=True)
plt.connect('key_press_event', toggle_selector)


############################################SELECTOR############################################









################LOAD_AND_SAVE_SELECTION_POINTS####################

save_selec_ax = plt.axes([0.91, 0.525, 0.08, 0.02])
button6 = Button(save_selec_ax, 'save_selec')
def save_selec_ax(event):
	x2 = np.nan_to_num(xnew_sel)
	y2 = np.nan_to_num(ynew_sel)
	err2 = np.nan_to_num(ynew_sel_err)
	data2 = np.transpose(np.array([x2, y2, err2]))
	if (os.path.isfile(file_sel)):
		save_prompt2=raw_input("File already exists. Overwrite?(y/n) : ")
        	if (save_prompt2=='y'):
        		np.savetxt(file_sel, data2)
        		print ('Selected data points loaded from file.')
        	else:
        		print ('data file not saved')
	else:
       		np.savetxt(file_sel, data2)
       		print ('saved')
button6.on_clicked(save_selec_ax)




get_selec_ax = plt.axes([0.91, 0.50, 0.08, 0.02])
button4 = Button(get_selec_ax, 'get_selc')
def get_selec_ax(event):
	x, y, err2 = np.loadtxt(file_sel, unpack=True)
	x[x==0] = np.nan
	y[y==0] = np.nan
	err2[err2==0] = np.nan
	global xnew_sel, ynew_sel, ynew_sel_err
	xnew_sel = x
	ynew_sel = y
	ynew_sel_err = err2
	line2.set_data(xnew_sel, ynew_sel)
	fig.canvas.draw()   
	print ('Selected data points loaded from file.')
button4.on_clicked(get_selec_ax)

################LOAD_AND_SAVE_SELECTION_POINTS####################










################################SAVE_AND_LOAD_PARAMETERS####################################

save_params_ax = plt.axes([0.91, 0.40, 0.08, 0.02])
button9 = Button(save_params_ax, 'save_params')
def save_params_ax(event):
	params_fitted = []
	z_qso = sqso.val
	z = sdla.val
	flux_red = sflux_red.val
	red_func = sred_func.val
	E_bv_qso = sred_qso_func.val
	box_pts = int(sdata_smooth.val)
	box_pts2 = int(sdata_smooth2.val)
	del_beta_new = sdel_beta.val
	raj_custom = np.array([str(radio.value_selected)])
	idx_best_class = int(np.where(raj == str(raj_custom[0]))[0][0])

	params_fitted = np.append(params_fitted, [z_qso, z, flux_red, red_func, E_bv_qso, box_pts, box_pts2, del_beta_new, idx_best_class]) 
	np.savetxt(file9, np.transpose(params_fitted))
	print ("Fitted Parameters Saved")
button9.on_clicked(save_params_ax)


get_params_ax = plt.axes([0.91, 0.370, 0.08, 0.02])
button10 = Button(get_params_ax, 'get_params')
def get_params_ax(event):
	params_fitted_new = np.loadtxt(file9)
	z_qso, z, flux_red, red_func, E_bv_qso, box_pts, box_pts2, del_beta_new, idx_best_class = params_fitted_new

	sqso.set_val(z_qso)
	sdla.set_val(z)
	sflux_red.set_val(flux_red)
	sred_func.set_val(red_func)
	sred_qso_func.set_val(E_bv_qso)
	sdata_smooth.set_val(box_pts)
	sdata_smooth2.set_val(box_pts2)
	sdel_beta.set_val(del_beta_new)
	radio.set_active(int(idx_best_class))

	print ("Fitted Parameters Loaded")
button10.on_clicked(get_params_ax)









############################################MAKE_LATEX_TABLE_FOR_PARAMETERS_AND_SAVE_PARAMETERS############################################

make_params_table_ax = plt.axes([0.005, 0.7, 0.1, 0.05])
button_make_params_table = Button(make_params_table_ax, 'print/save_latex')
def make_params_table_ax(event):

	global best_class, best_E_bv_dla_new, best_E_bv_dla_err_new, best_flux_redu_new, best_flux_redu_err_new, best_E_bv_qso_new, best_E_bv_qso_err_new, best_del_beta_fit_new, best_del_beta_fit_err_new, best_A_v_dla_new, best_A_v_dla_err_new, best_A_v_qso_new, best_A_v_qso_err_new, A_dla_bump, A_dla_bump_err, ratio, A_UV_dla, A_UV_qso

	fit_parameter_array = np.array([str(best_class), float(best_E_bv_dla_new), float(best_E_bv_dla_err_new), float(best_flux_redu_new), float(best_flux_redu_err_new), float(best_E_bv_qso_new), float(best_E_bv_qso_err_new), float(best_del_beta_fit_new), float(best_del_beta_fit_err_new), float(best_A_v_dla_new), float(best_A_v_dla_err_new), float(best_A_v_qso_new), float(best_A_v_qso_err_new), float(A_dla_bump), float(A_dla_bump_err), float(ratio), float(A_UV_dla), float(A_UV_qso)])

	np.savetxt(file10, fit_parameter_array, fmt='%s', delimiter=',')
	print ("\n" "\n")
	print ("Latex Parameter File saved")

	print ("Creating Parameter Table for Latex")
	print ("##################################")
	print ("\n" "\n")

	for i in range(len(fit_parameter_array)-1):
		fit_parameter_array[i+1] = np.round((np.float(fit_parameter_array[i+1])), 2)
	

	fitted_var_names = np.chararray([], itemsize=35)
	fitted_var_names = np.delete(fitted_var_names, 0)
	fitted_var_names = np.append(fitted_var_names, ['Extinction Class', 'E(B-V)(DLA)', 'E(B-V)(DLA)(err)', 'Flux Reduced', 'Flux Reduced Err', 'E(B-V)(QSO)', 'E(B-V)(QSO)(err)', 'Intrinsic QSO Shape', 'Intrinsic QSO Shape Err', 'A(V)(DLA)', 'A(V)(DLA)(err)', 'A(V)(QSO)', 'A(V)(QSO)(err)', 'A(bump)', 'A(bump) Err', 'A(UV) to A(V) ratio', 'A(UV)(DLA)', 'A(UV)(QSO)'])

	#print (len(fit_parameter_array), len(fitted_var_names))

	table_array = np.chararray([len(fit_parameter_array), 2], itemsize=35)
	table_array[:,0] = fitted_var_names
	table_array[:,1] = fit_parameter_array

	print(tabulate(table_array, tablefmt="latex", floatfmt="2.2f"))

	print ("\n" "\n")
	print ("##################################")
	print ("Parameter Table for Latex Created")

	
	
button_make_params_table.on_clicked(make_params_table_ax)

############################################MAKE_LATEX_TABLE_FOR_PARAMETERS_AND_SAVE_PARAMETERS############################################








############################################LOAD_PREVIOUS_PARAMETERS_FROM_A_SAVED_FILE############################################

load_params_table_ax = plt.axes([0.005, 0.65, 0.1, 0.05])
button_load_params_table = Button(load_params_table_ax, 'print/load_latex')
def load_params_table_ax(event):

	global best_class, best_E_bv_dla_new, best_E_bv_dla_err_new, best_flux_redu_new, best_flux_redu_err_new, best_E_bv_qso_new, best_E_bv_qso_err_new, best_del_beta_fit_new, best_del_beta_fit_err_new, best_A_v_dla_new, best_A_v_dla_err_new, best_A_v_qso_new, best_A_v_qso_err_new, A_dla_bump, A_dla_bump_err, ratio, A_UV_dla, A_UV_qso

	best_class, best_E_bv_dla_new, best_E_bv_dla_err_new, best_flux_redu_new, best_flux_redu_err_new, best_E_bv_qso_new, best_E_bv_qso_err_new, best_del_beta_fit_new, best_del_beta_fit_err_new, best_A_v_dla_new, best_A_v_dla_err_new, best_A_v_qso_new, best_A_v_qso_err_new, A_dla_bump, A_dla_bump_err, ratio, A_UV_dla, A_UV_qso = np.loadtxt(file10, unpack=True, dtype='S35')
	print ("\n" "\n")
	print ("Latex Parameter File loaded")


	print ("Creating Parameter Table for Latex")
	print ("##################################")
	print ("\n" "\n")


	fit_parameter_array = np.array([str(best_class), float(best_E_bv_dla_new), float(best_E_bv_dla_err_new), float(best_flux_redu_new), float(best_flux_redu_err_new), float(best_E_bv_qso_new), float(best_E_bv_qso_err_new), float(best_del_beta_fit_new), float(best_del_beta_fit_err_new), float(best_A_v_dla_new), float(best_A_v_dla_err_new), float(best_A_v_qso_new), float(best_A_v_qso_err_new), float(A_dla_bump), float(A_dla_bump_err), float(ratio), float(A_UV_dla), float(A_UV_qso)])

	for i in range(len(fit_parameter_array)-1):
		fit_parameter_array[i+1] = np.round((np.float(fit_parameter_array[i+1])), 2)
	

	fitted_var_names = np.chararray([], itemsize=35)
	fitted_var_names = np.delete(fitted_var_names, 0)
	fitted_var_names = np.append(fitted_var_names, ['Extinction Class', 'E(B-V)(DLA)', 'E(B-V)(DLA)(err)', 'Flux Reduced', 'Flux Reduced Err', 'E(B-V)(QSO)', 'E(B-V)(QSO)(err)', 'Intrinsic QSO Shape', 'Intrinsic QSO Shape Err', 'A(V)(DLA)', 'A(V)(DLA)(err)', 'A(V)(QSO)', 'A(V)(QSO)(err)', 'A(bump)', 'A(bump) Err', 'A(UV) to A(V) ratio', 'A(UV)(DLA)', 'A(UV)(QSO)'])
	table_array = np.chararray([len(fit_parameter_array), 2], itemsize=35)
	table_array[:,0] = fitted_var_names
	table_array[:,1] = fit_parameter_array

	print(tabulate(table_array, tablefmt="latex", floatfmt="2.2f"))

	print ("\n" "\n")
	print ("##################################")
	print ("Parameter Table for Latex Created")
button_load_params_table.on_clicked(load_params_table_ax)

############################################LOAD_PREVIOUS_PARAMETERS_FROM_A_SAVED_FILE############################################


################################SAVE_AND_LOAD_PARAMETERS####################################













##########################Creating_Update_for_the_Dynamic_actions##########################

def update(val):
	z_qso = sqso.val
	z = sdla.val
	flux_red = sflux_red.val
	red_func = sred_func.val
	E_bv_qso = sred_qso_func.val
	box_pts = int(sdata_smooth.val)
	box_pts2 = int(sdata_smooth2.val)
	del_beta_new = sdel_beta.val

	line5.set_ydata(smooth((reddening_func2(ynew, wave, red_func, R_v_dla, y_k_lambda_V_dla, E_bv_qso, R_v_qso, y_k_lambda_V_qso, del_beta_new)/flux_red),box_pts))
	line6.set_ydata(smooth((y_trish/flux_red),box_pts))
	line1.set_ydata(smooth(flux,box_pts2))


fig.canvas.draw_idle()
sqso.on_changed(update)
sdla.on_changed(update)
sflux_red.on_changed(update)
sred_func.on_changed(update)
sdata_smooth.on_changed(update)
sdata_smooth2.on_changed(update)
sdel_beta.on_changed(update)
sred_qso_func.on_changed(update)

##########################Creating_Update_for_the_Dynamic_actions##########################









##########################RADIO_BUTTON_FOR_SELECTING_EXTINCTION_CLASS##########################

rax2 = plt.axes([0.005, 0.2, 0.1, 0.15])
radio = RadioButtons(rax2, ('Galactic', 'SMC_Bar', 'LMC_Supershell', 'LMC_Average'), active=3)


def hzfunc(label2):
	hzdict = {'Galactic': 'Galactic', 'SMC_Bar': 'SMC_Bar', 'LMC_Supershell': 'LMC_Supershell', 'LMC_Average': 'LMC_Average'}

	extinction_class_string_initial = str(hzdict[label2])
	z_qso = sqso.val
	z = sdla.val
	flux_red = sflux_red.val
	E_bv_dla = sred_func.val
	E_bv_qso = sred_qso_func.val
	box_pts = int(sdata_smooth.val)
	box_pts2 = int(sdata_smooth2.val)
	del_beta_new = sdel_beta.val
	E_bv = E_bv_dla
	x0_dla, gamma_dla, c1_dla, c2_dla, c3_dla, c4_dla, c5_dla, O1_dla, O2_dla, O3_dla, R_v_dla, k_IR_dla = extinction(extinction_class_string_initial)
	x0_qso, gamma_qso, c1_qso, c2_qso, c3_qso, c4_qso, c5_qso, O1_qso, O2_qso, O3_qso, R_v_qso, k_IR_qso = extinction('SMC_Bar')
	c3_fit_new = c3_dla
	del_beta = del_beta_new
	y_k_lambda_V = spline_part_new_fit(extinction_class_string_initial, c3_fit_new)
	y_k_lambda_V_dla = spline_part_new_fit(extinction_class_string_initial, c3_fit_new)
	y_k_lambda_V_qso = spline_part_new_fit_qso()
	y_trish_new = reddening_func2(ynew, wave, E_bv_dla, R_v_dla, y_k_lambda_V_dla, E_bv_qso, R_v_qso, y_k_lambda_V_qso, del_beta)

	line5.set_ydata(y_trish_new/flux_red)
	plt.draw()
radio.on_clicked(hzfunc)

##########################RADIO_BUTTON_FOR_SELECTING_EXTINCTION_CLASS##########################









##########################MAKE_A_NICE_FIGURE##########################


plot_fig_ax = plt.axes([0.005, 0.6, 0.1, 0.05])
button_fig = Button(plot_fig_ax, 'Plot')
def plot_fig_ax(event):

	print ('Plotting.....')
	
	global x_wave, y_wave, y_wave_err, ynew1
	z_qso = sqso.val
	z = sdla.val
	best_flux_redu_new = sflux_red.val
	best_E_bv_dla_new = sred_func.val
	best_E_bv_qso_new = sred_qso_func.val
	box_pts = int(sdata_smooth.val)
	box_pts2 = int(sdata_smooth2.val)
	best_del_beta_fit_new = sdel_beta.val
	raj_custom = np.array([str(radio.value_selected)])
	idx_best_class = int(np.where(raj == str(raj_custom[0]))[0][0])
	best_class = str(raj_custom[0])

	x0_dla, gamma_dla, c1_dla, c2_dla, c3_dla, c4_dla, c5_dla, O1_dla, O2_dla, O3_dla, R_v_dla, k_IR_dla = extinction(best_class)
	x0_qso, gamma_qso, c1_qso, c2_qso, c3_qso, c4_qso, c5_qso, O1_qso, O2_qso, O3_qso, R_v_qso, k_IR_qso = extinction('SMC_Bar')
	
	best_A_v_dla_new = R_v_dla * best_E_bv_dla_new
	best_A_v_qso_new = R_v_qso * best_E_bv_qso_new
	c3_without_bump = 0.0
	A_dla_bump = (((np.pi)*c3_dla)/(2*gamma_dla*R_v_dla))*best_A_v_dla_new

	idx_V_band1 = np.searchsorted(wave, 5510)
	idx_UV_band1 = np.searchsorted(wave, 3650)
	idx_V_band2 = np.searchsorted(x_trish, 5510)
	idx_UV_band2 = np.searchsorted(x_trish, 3650)
	ratio = (y_trish[idx_UV_band2] - flux[idx_UV_band1]) / (y_trish[idx_V_band2] - flux[idx_V_band1]) 

	A_UV_dla = ratio * best_A_v_dla_new
	A_UV_qso = ratio * best_A_v_qso_new

	model_new_unreddened = (reddening_func2(ynew, wave, best_E_bv_dla_new, R_v_dla, spline_part_new_fit(best_class, c3_dla), best_E_bv_qso_new, R_v_qso, spline_part_new_fit_qso(), best_del_beta_fit_new)/1.0)
	model_new_reddened = (reddening_func2(ynew, wave, best_E_bv_dla_new, R_v_dla, spline_part_new_fit(best_class, c3_dla), best_E_bv_qso_new, R_v_qso, spline_part_new_fit_qso(), best_del_beta_fit_new)/best_flux_redu_new)
	model_new_reddened_without_bump = (reddening_func2(ynew, wave, best_E_bv_dla_new, R_v_dla, spline_part_new_fit(best_class, c3_without_bump), best_E_bv_qso_new, R_v_qso, spline_part_new_fit_qso(), best_del_beta_fit_new)/best_flux_redu_new)

	y_wave_plot = y_wave
	y_wave_plot[y_wave==0] = np.nan
	flux_plot = flux
	

	fig_new = plt.figure(2, figsize=(25, 16), dpi=300)
	ax_new = fig_new.add_subplot(111)

	ax_new.plot(wave,(smooth(flux_plot,box_pts2)), 'k-', alpha=0.4, drawstyle='steps', zorder=1, label='Data')	
	ax_new.plot(x_wave, (smooth(y_wave_plot,box_pts2)), 'k.', alpha=0.3, label='Selection', zorder=2)
	ax_new.plot(wave, (smooth(model_new_reddened,box_pts)), 'r-', alpha=0.6, zorder=7, label= '\n' ' Fit' '\n' r' A$\rm_{V}$(DLA)=%2.2f' '\n' r' A$\rm_{V}$(QSO)=%2.2f' '\n' r' $\rm \Delta \beta$(QSO)=%2.2f' %(float(best_A_v_dla_new), float(best_A_v_qso_new), float(best_del_beta_fit_new)))
	#ax_new.plot(wave, (smooth(model_new_unreddened,box_pts)), 'b-', alpha=0.5, zorder=3, label=' Selsing et al. 2016''\n' r' Extinction Class = %s' %str(best_class))
	ax_new.plot(x_trish, (y_trish/best_flux_redu_new), 'b-', alpha=0.5, zorder=3, label=' Selsing et al. 2016''\n' r' Extinction Class = %s' %str(best_class))
	ax_new.fill_between(wave, (smooth(model_new_reddened_without_bump,box_pts)), (smooth(model_new_reddened,box_pts)), facecolor='red', alpha=0.3, zorder=4, label= '\n' r' 2175 $\rm \AA$ bump' '\n' r' A$\rm_{bump}$=%2.2f' '\n' %(float(A_dla_bump)))
	ax_new.axvline(float(2175*(1+z)), linewidth=1, color='r', linestyle='--', label=r' 2175 $\rm \AA$ bump center', zorder=5)
	ax_new.axhline(0.0, linewidth=1, color='g', linestyle='--', zorder=6)
	ax_new.set_xlim((wave.min()-(0.1*wave.mean())), (wave.max()+(0.1*wave.mean())))
	ax_new.set_ylim((flux.min()-(0.1*flux.mean())), (flux.max()+(2*flux.mean())))
	ax_new.legend()

	plt.show(2)
	
button_fig.on_clicked(plot_fig_ax)






plot_save_fig_ax = plt.axes([0.005, 0.55, 0.1, 0.05])
button2_fig = Button(plot_save_fig_ax, 'Save Fig')
def plot_save_fig_ax(event):

	print ('Plotting.....')
	
	global x_wave, y_wave, y_wave_err, ynew1
	z_qso = sqso.val
	z = sdla.val
	best_flux_redu_new = sflux_red.val
	best_E_bv_dla_new = sred_func.val
	best_E_bv_qso_new = sred_qso_func.val
	box_pts = int(sdata_smooth.val)
	box_pts2 = int(sdata_smooth2.val)
	best_del_beta_fit_new = sdel_beta.val
	raj_custom = np.array([str(radio.value_selected)])
	idx_best_class = int(np.where(raj == str(raj_custom[0]))[0][0])
	best_class = str(raj_custom[0])

	x0_dla, gamma_dla, c1_dla, c2_dla, c3_dla, c4_dla, c5_dla, O1_dla, O2_dla, O3_dla, R_v_dla, k_IR_dla = extinction(best_class)
	x0_qso, gamma_qso, c1_qso, c2_qso, c3_qso, c4_qso, c5_qso, O1_qso, O2_qso, O3_qso, R_v_qso, k_IR_qso = extinction('SMC_Bar')
	
	best_A_v_dla_new = R_v_dla * best_E_bv_dla_new
	best_A_v_qso_new = R_v_qso * best_E_bv_qso_new
	c3_without_bump = 0.0
	A_dla_bump = (((np.pi)*c3_dla)/(2*gamma_dla*R_v_dla))*best_A_v_dla_new

	idx_V_band1 = np.searchsorted(wave, 5510)
	idx_UV_band1 = np.searchsorted(wave, 3650)
	idx_V_band2 = np.searchsorted(x_trish, 5510)
	idx_UV_band2 = np.searchsorted(x_trish, 3650)
	ratio = (y_trish[idx_UV_band2] - flux[idx_UV_band1]) / (y_trish[idx_V_band2] - flux[idx_V_band1]) 

	A_UV_dla = ratio * best_A_v_dla_new
	A_UV_qso = ratio * best_A_v_qso_new

	model_new_unreddened = (reddening_func2(ynew, wave, best_E_bv_dla_new, R_v_dla, spline_part_new_fit(best_class, c3_dla), best_E_bv_qso_new, R_v_qso, spline_part_new_fit_qso(), best_del_beta_fit_new)/1.0)
	model_new_reddened = (reddening_func2(ynew, wave, best_E_bv_dla_new, R_v_dla, spline_part_new_fit(best_class, c3_dla), best_E_bv_qso_new, R_v_qso, spline_part_new_fit_qso(), best_del_beta_fit_new)/best_flux_redu_new)
	model_new_reddened_without_bump = (reddening_func2(ynew, wave, best_E_bv_dla_new, R_v_dla, spline_part_new_fit(best_class, c3_without_bump), best_E_bv_qso_new, R_v_qso, spline_part_new_fit_qso(), best_del_beta_fit_new)/best_flux_redu_new)

	y_wave_plot = y_wave
	y_wave_plot[y_wave==0] = np.nan
	flux_plot = flux
	

	fig_new = plt.figure(3, figsize=(25, 16), dpi=300)
	ax_new = fig_new.add_subplot(111, rasterized=False)

	ax_new.plot(wave,(smooth(flux_plot,box_pts2)), 'k-', alpha=0.4, drawstyle='steps', zorder=1, label='Data')	
	ax_new.plot(x_wave, (smooth(y_wave_plot,box_pts2)), 'k.', alpha=0.3, label='Selection', zorder=2)
	ax_new.plot(wave, (smooth(model_new_reddened,box_pts)), 'r-', alpha=0.6, zorder=7, label= '\n' ' Fit' '\n' r' A$\rm_{V}$(DLA)=%2.2f'  %(float(best_A_v_dla_new)))
	#ax_new.plot(wave, (smooth(model_new_unreddened,box_pts)), 'b-', alpha=0.5, zorder=3, label=' Selsing et al. 2016''\n' r' Extinction Class = %s' %str(best_class))
	ax_new.plot(x_trish, (y_trish/best_flux_redu_new), 'b-', alpha=0.5, zorder=3, label=' Selsing et al. 2016''\n' r' Extinction Class = %s' %str(best_class))
	ax_new.fill_between(wave, (smooth(model_new_reddened_without_bump,box_pts)), (smooth(model_new_reddened,box_pts)), facecolor='red', alpha=0.3, zorder=4, label= '\n' r' 2175 $\rm \AA$ bump' '\n' r' A$\rm_{bump}$=%2.2f' '\n' %(float(A_dla_bump)))
	ax_new.axvline(float(2175*(1+z)), linewidth=1, color='r', linestyle='--', label=r' 2175 $\rm \AA$ bump center', zorder=5)
	ax_new.axhline(0.0, linewidth=1, color='g', linestyle='--', zorder=6)

	#ax_new.set_xlim((wave.min()-(0.1*wave.mean())), (wave.max()+(0.1*wave.mean())))
	#ax_new.set_ylim((flux.min()-(0.1*flux.mean())), (flux.max()+(2*flux.mean())))

	ax_new.set_xlim((wave.min()-(0.0*wave.mean())), (wave.max()+(0.0*wave.mean())))
	ax_new.set_ylim(-0.2, 8.)

	ax_new.tick_params(axis = 'both', which = 'major', direction='in', length=size_of_font, width=2, colors='k', labelsize=size_of_font)
	ax_new.tick_params(axis = 'both', which = 'minor', direction='in', length=0.5*size_of_font, width=1, colors='k', labelsize=size_of_font)

	ax_new.yaxis.set_major_locator(ticker.MultipleLocator(2))
	ax_new.yaxis.set_minor_locator(ticker.MultipleLocator(0.4))

	ax_new.xaxis.set_major_locator(ticker.MultipleLocator(5000))
	ax_new.xaxis.set_minor_locator(ticker.MultipleLocator(1000))


	ax_new.legend(fontsize=0.75*size_of_font)

	ax_new.set_xlabel(r'Wavelength $\AA$', fontsize=size_of_font)
	ax_new.set_ylabel(r'Flux ($10^{-17}$ $erg$ $cm^{-2}$ $s^{-1}$ $AA^{-1}$)', fontsize=size_of_font)
	plt.title(str(file_name1)[:-4], fontsize=1.5*size_of_font)
	#plt.legend(fontsize=size_of_font)

	file_name = str(file_name1)[:-4] + '_dust_extinction_v.pdf'
	print (file_name)
	plt.savefig(file_name)
	print ("Figure saved")
	#plt.show(3)
	
button2_fig.on_clicked(plot_save_fig_ax)


##########################MAKE_A_NICE_ASS_FIGURE##########################









##########################FITTING_ESSENTIALS##########################

loadtxt_file_name = str(file_name1)[:-4] + '_selected_data.txt'


try:
	f = open(str(loadtxt_file_name))
	x_wave, y_wave, y_wave_err = np.loadtxt(loadtxt_file_name, unpack=True)
	ynew1 = np.zeros([len(ynew)])
	y_wave_err[y_wave==0.0] = 0.00000000001

	for i in range(len(y_wave)):
		if y_wave[i]!=0.0:
			ynew1[i] = ynew[i]
		else:
			ynew1[i] = 0.0

except IOError:
	print("Selected data file not found. Setting all values to zero")

	x_wave, y_wave, y_wave_err = np.loadtxt(str(file_name1), unpack=True)
	y_wave[:] = 0.0
	y_wave_err[:] = 0.0






def func_5(xdata, *t):
	global y_wave
	E_bv_dla, flux_redu, E_bv_qso, del_beta=t
	x0_dla, gamma_dla, c1_dla, c2_dla, c3_dla, c4_dla, c5_dla, O1_dla, O2_dla, O3_dla, R_v_dla, k_IR_dla = extinction(extinction_class_string)
	x0_qso, gamma_qso, c1_qso, c2_qso, c3_qso, c4_qso, c5_qso, O1_qso, O2_qso, O3_qso, R_v_qso, k_IR_qso = extinction('SMC_Bar')
	model = ((reddening_func2(ynew1, wave, E_bv_dla, R_v_dla, spline_part_new_fit(extinction_class_string, c3_dla), E_bv_qso, R_v_qso, spline_part_new_fit_qso(), del_beta))/flux_redu)
	model[y_wave==0] = 0.0
	return (model)



def func_6(xdata, *t):
	global y_wave
	E_bv_dla, flux_redu, E_bv_qso, del_beta=t
	x0_dla, gamma_dla, c1_dla, c2_dla, c3_dla, c4_dla, c5_dla, O1_dla, O2_dla, O3_dla, R_v_dla, k_IR_dla = extinction(extinction_class_string)
	x0_qso, gamma_qso, c1_qso, c2_qso, c3_qso, c4_qso, c5_qso, O1_qso, O2_qso, O3_qso, R_v_qso, k_IR_qso = extinction('SMC_Bar')
	model = ((reddening_func2(ynew1, wave, E_bv_dla, R_v_dla, spline_part_new_fit(extinction_class_string, c3_dla), 0., R_v_qso, spline_part_new_fit_qso(), 0.))/flux_redu)
	model[y_wave==0] = 0.0
	return (model)





def fit_curvefit_complete_red(p0, datax, datay, yerror, function, method_str, **kwargs):

	#Choosing reasonable default bounds. Please change if required.
	bounds_lower = np.zeros([len(p0)])
	bounds_upper = np.zeros([len(p0)])
	bounds_lower[:] = 0.0
	bounds_upper[:] = np.inf
	bounds_lower[1] = 0.01
	#bounds_lower[3] = -np.inf
	bounds_lower[2] = 0.
	bounds_upper[2] = 1.
	bounds_lower[3] = -1.
	bounds_upper[3] = 1.


	pfit, pcov = curve_fit(function,datax,datay,p0=p0, bounds=((bounds_lower), (bounds_upper)), sigma=yerror, method=method_str, **kwargs)

	error = [] 
	for i in range(len(pfit)):
		try:
			error.append(np.absolute(pcov[i][i])**0.5)
		except:
			error.append(0.00)
	pfit_curvefit = pfit
	perr_curvefit = np.array(error)

	return pfit_curvefit, perr_curvefit




def fitting_function(raj_custom):
	print ('Fitting in progress......')
	global x_wave, y_wave, y_wave_err, ynew1
	choose_class = np.zeros([len(raj_custom)])
	E_bv_dla_new = np.zeros([len(raj_custom)])
	E_bv_qso_new = np.zeros([len(raj_custom)])
	flux_redu_new = np.zeros([len(raj_custom)])
	E_bv_dla_new_err = np.zeros([len(raj_custom)])
	E_bv_qso_new_err = np.zeros([len(raj_custom)])
	flux_redu_new_err = np.zeros([len(raj_custom)])
	c3_fit_new = np.zeros([len(raj_custom)])
	c3_fit_new_err = np.zeros([len(raj_custom)])
	del_beta_fit_new = np.zeros([len(raj_custom)])
	del_beta_fit_new_err = np.zeros([len(raj_custom)])

	E_bv_dla = sred_func.val
	flux_redu = sflux_red.val
	E_bv_qso = sred_qso_func.val
	del_beta = sdel_beta.val
	data_smoothing_parameter = int(sdata_smooth2.val)


	for i in range(len(raj_custom)):
		global extinction_class_string
		extinction_class_string = str(raj_custom[i])
		x0_dla, gamma_dla, c1_dla, c2_dla, c3_dla, c4_dla, c5_dla, O1_dla, O2_dla, O3_dla, R_v_dla, k_IR_dla = extinction(extinction_class_string)
		x0_qso, gamma_qso, c1_qso, c2_qso, c3_qso, c4_qso, c5_qso, O1_qso, O2_qso, O3_qso, R_v_qso, k_IR_qso = extinction('SMC_Bar')
		
		y_wave = smooth(y_wave,data_smoothing_parameter)
		x0_fit = np.array([E_bv_dla, flux_redu, E_bv_qso, del_beta])

		method_str = 'trf'
		#method_str = 'dogbox'

		#popt, perr = fit_curvefit_complete_red(x0_fit, wave, y_wave, y_wave_err, func_5, method_str)
		popt, perr = fit_curvefit_complete_red(x0_fit, wave, y_wave, y_wave_err, func_6, method_str)


		E_bv_dla_new[i] = popt[0]
		flux_redu_new[i] = popt[1]
		E_bv_qso_new[i] = popt[2]
		del_beta_fit_new[i] = popt[3]
		E_bv_dla_new_err[i] = perr[0]
		flux_redu_new_err[i] = perr[1]
		E_bv_qso_new_err[i] = perr[2]
		del_beta_fit_new_err[i] = perr[3]
		model = func_5(wave, *popt)
		chi_squared = np.zeros([len(wave)])

		for j in range(len(wave)):
			if (y_wave_err[j]!=0):
				chi_squared[j] = ((y_wave[j] - model[j])/y_wave_err[j])

		chi_sq = (np.abs(np.sum(chi_squared[:])))/(len(wave)-len(popt))
		choose_class[i] = float(chi_sq)
		print ('For Extinction Class: ',("%s" %extinction_class_string))
		print ('The chi squared value is: ',("%.6f" %chi_sq))

	print (raj_custom)
	print (choose_class)
	print ('Fitting complete......')


	for i in range(len(raj_custom)):
		if (choose_class[i]==choose_class.min()):
			best_class = str(raj_custom[i])
			best_E_bv_dla_new = float(E_bv_dla_new[i])
			best_flux_redu_new = float(flux_redu_new[i])
			best_E_bv_qso_new = float(E_bv_qso_new[i])
			best_del_beta_fit_new = float(del_beta_fit_new[i])
			best_E_bv_dla_err_new = float(E_bv_dla_new_err[i])
			best_flux_redu_err_new = float(flux_redu_new_err[i])
			best_E_bv_qso_err_new = float(E_bv_qso_new_err[i])
			best_del_beta_fit_err_new = float(del_beta_fit_new_err[i])

	return (best_class, best_E_bv_dla_new, best_flux_redu_new, best_E_bv_qso_new, best_del_beta_fit_new, best_E_bv_dla_err_new, best_flux_redu_err_new, best_E_bv_qso_err_new, best_del_beta_fit_err_new)
	

##########################FITTING_ESSENTIALS##########################











##############################################################################
##########################FITTING_BUTTONS#####################################
##############################################################################


##########################COMPLETE_FIT_BUTTON##########################

fitax = plt.axes([0.91, 0.625, 0.08, 0.02]) 
button1 = Button(fitax, 'Fit', color=axcolor, hovercolor='0.975')
def fitax(event):

	global best_class, best_E_bv_dla_new, best_E_bv_dla_err_new, best_flux_redu_new, best_flux_redu_err_new, best_E_bv_qso_new, best_E_bv_qso_err_new, best_del_beta_fit_new, best_del_beta_fit_err_new, best_A_v_dla_new, best_A_v_dla_err_new, best_A_v_qso_new, best_A_v_qso_err_new, A_dla_bump, A_dla_bump_err, ratio, A_UV_dla, A_UV_qso

	raj_custom = raj
	best_class, best_E_bv_dla_new, best_flux_redu_new, best_E_bv_qso_new, best_del_beta_fit_new, best_E_bv_dla_err_new, best_flux_redu_err_new, best_E_bv_qso_err_new, best_del_beta_fit_err_new = fitting_function(raj_custom)

	x0_dla, gamma_dla, c1_dla, c2_dla, c3_dla, c4_dla, c5_dla, O1_dla, O2_dla, O3_dla, R_v_dla, k_IR_dla = extinction(best_class)
	x0_qso, gamma_qso, c1_qso, c2_qso, c3_qso, c4_qso, c5_qso, O1_qso, O2_qso, O3_qso, R_v_qso, k_IR_qso = extinction('SMC_Bar')

	print ('best_class - %s' %str(best_class))
	print ('best_E_bv_dla_new - %2.2f' %float(best_E_bv_dla_new))
	print ('best_E_bv_dla_err_new - %2.2f' %float(best_E_bv_dla_err_new))
	print ('best_flux_redu_new - %2.2f' %float(best_flux_redu_new))
	print ('best_flux_redu_err_new - %2.2f' %float(best_flux_redu_err_new))
	print ('best_E_bv_qso_new - %2.2f' %float(best_E_bv_qso_new))
	print ('best_E_bv_qso_err_new - %2.2f' %float(best_E_bv_qso_err_new))
	print ('best_del_beta_fit_new - %2.2f' %float(best_del_beta_fit_new))
	print ('best_del_beta_fit_err_new - %2.2f' %float(best_del_beta_fit_err_new))
	
	best_A_v_dla_new = R_v_dla * best_E_bv_dla_new
	best_A_v_dla_err_new = R_v_dla * best_E_bv_dla_err_new
	best_A_v_qso_new = R_v_qso * best_E_bv_qso_new
	best_A_v_qso_err_new = R_v_qso * best_E_bv_qso_err_new
	
	print ('best_A_v_dla_new %2.2f' %float(best_A_v_dla_new))
	print ('best_A_v_dla_err_new %2.2f' %float(best_A_v_dla_err_new))
	print ('best_A_v_qso_new %2.2f' %float(best_A_v_qso_new))
	print ('best_A_v_qso_err_new %2.2f' %float(best_A_v_qso_err_new))
	c3_without_bump = 0.0

	A_dla_bump = (((np.pi)*c3_dla)/(2*gamma_dla*R_v_dla))*best_A_v_dla_new
	A_dla_bump_err = (((np.pi)*c3_dla)/(2*gamma_dla*R_v_dla))*best_A_v_dla_err_new

	print ('A_dla_bump - %2.2f' %float(A_dla_bump))
	print ('A_dla_bump_err - %2.2f' %float(A_dla_bump_err))

	idx_V_band1 = np.searchsorted(wave, 5510)
	idx_UV_band1 = np.searchsorted(wave, 3650)
	idx_V_band2 = np.searchsorted(x_trish, 5510)
	idx_UV_band2 = np.searchsorted(x_trish, 3650)
	ratio = (y_trish[idx_UV_band2] - flux[idx_UV_band1]) / (y_trish[idx_V_band2] - flux[idx_V_band1]) 

	A_UV_dla = ratio * best_A_v_dla_new
	A_UV_qso = ratio * best_A_v_qso_new
	print ('A_UV to A_V ratio - %2.2f'  %float(ratio))
	print ('A_UV_dla - %2.2f' %float(A_UV_dla))
	print ('A_UV_qso - %2.2f' %float(A_UV_qso))

	idx_best_class = np.where(raj == str(best_class))[0][0]
	radio.set_active(idx_best_class)

	sred_func.set_val(best_E_bv_dla_new)	
	sflux_red.set_val(best_flux_redu_new)	
	sred_qso_func.set_val(best_E_bv_qso_new)	
	sdel_beta.set_val(best_del_beta_fit_new)

	model_new_unreddened = (reddening_func2(ynew, wave, best_E_bv_dla_new, R_v_dla, spline_part_new_fit(best_class, c3_dla), best_E_bv_qso_new, R_v_qso, spline_part_new_fit_qso(), best_del_beta_fit_new)/1.0)
	model_new_reddened = (reddening_func2(ynew, wave, best_E_bv_dla_new, R_v_dla, spline_part_new_fit(best_class, c3_dla), best_E_bv_qso_new, R_v_qso, spline_part_new_fit_qso(), best_del_beta_fit_new)/best_flux_redu_new)
	model_new_reddened_without_bump = (reddening_func2(ynew, wave, best_E_bv_dla_new, R_v_dla, spline_part_new_fit(best_class, c3_without_bump), best_E_bv_qso_new, R_v_qso, spline_part_new_fit_qso(), best_del_beta_fit_new)/best_flux_redu_new)
	
	line5.set_ydata(model_new_reddened)

button1.on_clicked(fitax)

##########################COMPLETE_FIT_BUTTON##########################








##########################SINGLE_FIT_BUTTON##########################

fit_sel_ax = plt.axes([0.91, 0.600, 0.08, 0.02]) 
button_sel_1 = Button(fit_sel_ax, 'Fit_Selected', color=axcolor, hovercolor='0.975')
def fit_sel_ax(event):

	global best_class, best_E_bv_dla_new, best_E_bv_dla_err_new, best_flux_redu_new, best_flux_redu_err_new, best_E_bv_qso_new, best_E_bv_qso_err_new, best_del_beta_fit_new, best_del_beta_fit_err_new, best_A_v_dla_new, best_A_v_dla_err_new, best_A_v_qso_new, best_A_v_qso_err_new, A_dla_bump, A_dla_bump_err, ratio, A_UV_dla, A_UV_qso

	raj_custom = np.array([str(radio.value_selected)])
	best_class, best_E_bv_dla_new, best_flux_redu_new, best_E_bv_qso_new, best_del_beta_fit_new, best_E_bv_dla_err_new, best_flux_redu_err_new, best_E_bv_qso_err_new, best_del_beta_fit_err_new = fitting_function(raj_custom)

	x0_dla, gamma_dla, c1_dla, c2_dla, c3_dla, c4_dla, c5_dla, O1_dla, O2_dla, O3_dla, R_v_dla, k_IR_dla = extinction(best_class)
	x0_qso, gamma_qso, c1_qso, c2_qso, c3_qso, c4_qso, c5_qso, O1_qso, O2_qso, O3_qso, R_v_qso, k_IR_qso = extinction('SMC_Bar')

	print ('best_class - %s' %str(best_class))
	print ('best_E_bv_dla_new - %2.2f' %float(best_E_bv_dla_new))
	print ('best_E_bv_dla_err_new - %2.2f' %float(best_E_bv_dla_err_new))
	print ('best_flux_redu_new - %2.2f' %float(best_flux_redu_new))
	print ('best_flux_redu_err_new - %2.2f' %float(best_flux_redu_err_new))
	print ('best_E_bv_qso_new - %2.2f' %float(best_E_bv_qso_new))
	print ('best_E_bv_qso_err_new - %2.2f' %float(best_E_bv_qso_err_new))
	print ('best_del_beta_fit_new - %2.2f' %float(best_del_beta_fit_new))
	print ('best_del_beta_fit_err_new - %2.2f' %float(best_del_beta_fit_err_new))
	
	best_A_v_dla_new = R_v_dla * best_E_bv_dla_new
	best_A_v_dla_err_new = R_v_dla * best_E_bv_dla_err_new
	best_A_v_qso_new = R_v_qso * best_E_bv_qso_new
	best_A_v_qso_err_new = R_v_qso * best_E_bv_qso_err_new
	
	print ('best_A_v_dla_new %2.2f' %float(best_A_v_dla_new))
	print ('best_A_v_dla_err_new %2.2f' %float(best_A_v_dla_err_new))
	print ('best_A_v_qso_new %2.2f' %float(best_A_v_qso_new))
	print ('best_A_v_qso_err_new %2.2f' %float(best_A_v_qso_err_new))
	c3_without_bump = 0.0

	model_new_unreddened = (reddening_func2(ynew, wave, best_E_bv_dla_new, R_v_dla, spline_part_new_fit(best_class, c3_dla), best_E_bv_qso_new, R_v_qso, spline_part_new_fit_qso(), best_del_beta_fit_new)/1.0)
	model_new_reddened = (reddening_func2(ynew, wave, best_E_bv_dla_new, R_v_dla, spline_part_new_fit(best_class, c3_dla), best_E_bv_qso_new, R_v_qso, spline_part_new_fit_qso(), best_del_beta_fit_new)/best_flux_redu_new)
	model_new_reddened_without_bump = (reddening_func2(ynew, wave, best_E_bv_dla_new, R_v_dla, spline_part_new_fit(best_class, c3_without_bump), best_E_bv_qso_new, R_v_qso, spline_part_new_fit_qso(), best_del_beta_fit_new)/best_flux_redu_new)

	A_dla_bump = (((np.pi)*c3_dla)/(2*gamma_dla*R_v_dla))*best_A_v_dla_new
	A_dla_bump_err = (((np.pi)*c3_dla)/(2*gamma_dla*R_v_dla))*best_A_v_dla_err_new

	print ('A_dla_bump - %2.2f' %float(A_dla_bump))
	print ('A_dla_bump_err - %2.2f' %float(A_dla_bump_err))

	idx_V_band1 = np.searchsorted(wave, 5510)
	idx_UV_band1 = np.searchsorted(wave, 3650)
	idx_V_band2 = np.searchsorted(x_trish, 5510)
	idx_UV_band2 = np.searchsorted(x_trish, 3650)
	ratio = (y_trish[idx_UV_band2] - flux[idx_UV_band1]) / (y_trish[idx_V_band2] - flux[idx_V_band1]) 

	A_UV_dla = ratio * best_A_v_dla_new
	A_UV_qso = ratio * best_A_v_qso_new
	print ('A_UV to A_V ratio - %2.2f'  %float(ratio))
	print ('A_UV_dla - %2.2f' %float(A_UV_dla))
	print ('A_UV_qso - %2.2f' %float(A_UV_qso))


	sred_func.set_val(best_E_bv_dla_new)	
	sflux_red.set_val(best_flux_redu_new)	
	sred_qso_func.set_val(best_E_bv_qso_new)	
	sdel_beta.set_val(best_del_beta_fit_new)	
	
	line5.set_ydata(model_new_reddened)

button_sel_1.on_clicked(fit_sel_ax)

##########################SINGLE_FIT_BUTTON##########################


##############################################################################
##########################FITTING_BUTTONS#####################################
##############################################################################












#####################################PLOT_MAIN_WINDOW#####################################

ax.set_xlim((wave.min()-(0.1*wave.mean())), (wave.max()+(0.1*wave.mean())))
ax.set_ylim((flux.min()-(0.1*flux.mean())), (flux.max()+(2*flux.mean())))

#SHOW_WINDOW
plt.show(1)

#####################################PLOT_MAIN_WINDOW#####################################







