# -*- coding: utf-8 -*-
# niobium @ 04.10.2018
# Thanks for stimulating discussions to Vladimir Berzin, Kirill Emelyanov, Nickolay Boldyrev

# calculation of intensity spectra from interferogram (presented in ascii file)
# taking into account apodization and phase correction
# phase correction: calculation of phase_shift(w) from shift of ZPD (zero path difference)
# zero filling for (quasi-)interpolation
# positive peaks (maxima) detection in the assigned ranges
# info about FTIR: krylov_a_s_vtyurin_a_n_gerasimova_yu_v_obrabotka_dannykh_inf.pdf

import math
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
#from scipy import signal
#from sklearn.linear_model import LinearRegression

#IN_FILENAME  = 'interferogram_air.ascii'
IN_FILENAME  = 'interferogram_HCl.ascii'
OUT_FILENAME = 'spectrum_peaks_out.txt'

APODIZE_TYPE = 'BH' # 'tr' - triangle, 'exp' - exponent, 'HG' - Happ-Genzel, 'BH' - Blackmann-Harris

lim_x_min=500 # 1/cm, xlimits for output only
lim_x_max=5000 # 1/cm, xlimits for output only

p_branch_min=2600 #start of range while looking for P branch, 1/cm
p_branch_max=2870 #end of range while looking for P branch, 1/cm
r_branch_min=2900 #start of range while looking for R branch, 1/cm
r_branch_max=3200 #end of range while looking for R branch, 1/cm

# ZPD = zero path difference (nulevaya raznost' khoda)
ZPD_position    = 2200 #2151 2123 # approx index of ZPD in interferogram, try to change me !!!
ZPD_exact_shift = 0. # value from -1.0 to 1.0, if ZPD is really close

usefull_range   = 21500 # length of long side of the interferogram, try to change me !!!
base_range      = 200000 # range for Fourier transform, try to change me !!! 
base_range_half = 100001 # range after Fourier transform = base_range/2 + 1
wavelength      = 632.766 # laser wavelength, nm 

lines_number  = 10   # quantity of required lines in P,R branches
lines_ignored = 2    # quantity of dots|lines to be ignored in fitting
dots_ignored = 2        # quantity of dots|branches to be ignored in fitting
delt = 100.             # peak detect trigger !!!

base_for_filtration = 500 # base for low-pass Fourier filter and peak detect !!!

def max_detect_decr(xx,yy,min_xx,max_xx,how_many):
    maxtab = []
    peaks_num=0
    for i in np.arange(len(xx)):
        j = len(xx)-5-i
        if (xx[j]<max_xx)and(yy[j]>0):
             if (yy[j-1]<yy[j])and(yy[j+1]<yy[j])and(yy[j-3]<yy[j]-delt)and(yy[j+3]<yy[j]-delt):
                 maxtab.append((j, xx[j], yy[j]))
                 peaks_num=peaks_num+1
        elif (xx[j]<min_xx): 
            break
        if(peaks_num==how_many): break
    return maxtab

def max_detect_incr(xx,yy,min_xx,max_xx,how_many):
    maxtab = []
    peaks_num=0
    for j in np.arange(len(xx)):
        if (xx[j]>min_xx)and(yy[j]>0):
             if (yy[j-1]<yy[j])and(yy[j+1]<yy[j])and(yy[j-3]<yy[j]-delt)and(yy[j+3]<yy[j]-delt):
                 maxtab.append((j, xx[j], yy[j]))
                 peaks_num=peaks_num+1
        elif (xx[j]>max_xx): 
            break
        if(peaks_num==how_many): break
    return maxtab

def elongate_apodize(v, fin_range, apod_type):
    elongated = np.zeros(fin_range)
    if (apod_type=='no'):                # no apodization
        for i in np.arange(len(v)):
            elongated[i]=v[i]
    elif (apod_type=='tr')or(apod_type=='tri')or(apod_type=='triangle'): # triangle apodization
        ini_range = len(v)
        for i in np.arange(ini_range):
            i_norm=np.float64(i)/np.float64(fin_range)
            elongated[i]=v[i]*(1.-i_norm)
    elif (apod_type=='HG'):              # Happ-Genzel or Hamming apodization
        ini_range = len(v)
        for i in np.arange(ini_range):
            i_norm=math.pi*i/np.float64(ini_range)
            elongated[i]=v[i]*(0.54+0.46*math.cos(i_norm)) 
    elif (apod_type=='BH'):              # Blackmann-Harris apodization
        ini_range = len(v)
        for i in np.arange(ini_range):
            i_norm=math.pi*i/np.float64(ini_range)
            elongated[i]=v[i]*(0.35875+0.48829*math.cos(i_norm)+0.14128*math.cos(2.*i_norm)+0.01168*math.cos(3.*i_norm))        
    elif (apod_type=='exp'):             # exponential
        ini_range = len(v)
        for i in np.arange(ini_range):
            i_norm=np.float64(i)/np.float64(fin_range)
            elongated[i]=v[i]*(math.exp(-i_norm*5.))
    return elongated

def low_pass_filtration(v, zero_fin_range, fin_range):
    filtered = np.zeros(fin_range)
    for i in np.arange(zero_fin_range):
        i_norm=math.pi*i/np.float64(fin_range)
        filtered[i]=v[i]*(0.35875+0.48829*math.cos(i_norm)+0.14128*math.cos(2.*i_norm)+0.01168*math.cos(3.*i_norm))
#       Blackmann-Harris apodization
    return filtered

def phase_correct_atan(spectrum_in_y):
    spectrum_out_y = []
    leng=len(phase_y)
    for j in np.arange(leng):
        spectrum_out_y.append(np.real(spectrum_in_y[j])*math.cos(phase_y[j])
        +np.imag(spectrum_in_y[j])*math.sin(phase_y[j]))
    return spectrum_out_y

def phase_correct_powspec(spectrum_in_y):
    spectrum_out_y = []
    for j in np.arange(len(phase_y)):
        spectrum_out_y.append(math.sqrt(np.real(spectrum_in_y[j])**2+np.imag(spectrum_in_y[j])**2))
    return spectrum_out_y

def offset_compensate(spectrum_out_y):
    leng=len(spectrum_out_y)
    offset=10
    slope = (spectrum_out_y[offset]-spectrum_out_y[leng-1-offset])/np.float64(leng-2*offset)
    for j in np.arange(leng):
        spectrum_out_y[j]=spectrum_out_y[j]-np.float64(leng-j)*slope-spectrum_out_y[leng-1-offset]

def parabola_shifted(x, x0, a0, k):
    return a0+k*(x-x0)**2

def ZPD_search(v,ZPD_approx):
    j_max = ZPD_approx
    max_val=v[ZPD_approx]
    for j in np.arange(ZPD_approx*2):
        if(v[j]>max_val):
            j_max=j
            max_val=v[j]
    ZPD_search = 0.
    interf_xx = []
    interf_yy = []
    for j in np.arange(7):
        interf_yy.append(v[j_max+j-3])
        interf_xx.append(np.float64(j-3))
    popt,pconv=curve_fit(parabola_shifted, interf_xx, interf_yy, p0=[0.,max_val,0.])
    ZPD_search=popt[0]+np.float64(j_max-ZPD_approx)
    print(interf_yy)
    print(popt)
    print(j_max)
    print(ZPD_search)
    return ZPD_search

# parsing of (ascii) csv file with interferogram
dataframe = pd.read_csv(IN_FILENAME, header=None, comment='#', skiprows=0, sep='\t')

array = dataframe.values

signal_x = array[:,0]
signal_x = signal_x.astype(np.float64)
signal_y = array[:,1]
signal_y = signal_y.astype(np.float64)

ZPD_exact_shift=ZPD_search(signal_y, ZPD_position)
if (int(ZPD_exact_shift)<-5): 
    ZPD_position=int(ZPD_position+ZPD_exact_shift)
    print('Maximum out of range!!! ZPD corrected:'); print(ZPD_position)

# extraction of long-side part of interferogram
index_range = np.arange(usefull_range)+ZPD_position
signal_x_1 = np.take(signal_x,index_range)
signal_y_1 = np.take(signal_y,index_range)

# elongation of interferogram with zero filling and apodization (or not)
signal_x_range = np.arange(base_range)
signal_y_2_ap = elongate_apodize(signal_y_1,base_range,APODIZE_TYPE)
signal_y_2_no = elongate_apodize(signal_y_1,base_range,'no')
signal_y_2_sm = low_pass_filtration(signal_y_1,base_for_filtration,base_range)

# DFT (discrete Fourier transform) of array: N (real) to N/2+1 (complex, Hermitian conjugate at w and -w)
scale_x = 2.0e7/wavelength/np.float64(base_range)
spectrum_x = np.arange(base_range_half)*scale_x
spectrum_y_2_ap = np.fft.rfft(signal_y_2_ap,base_range)
spectrum_y_2_no = np.fft.rfft(signal_y_2_no,base_range)
spectrum_y_2_sm = np.fft.rfft(signal_y_2_sm,base_range)

# phase shift estimation from shift of ZPD (zero path difference)
ZPD_exact_shift=ZPD_search(signal_y, ZPD_position)
    
phase_y=[]
for j in np.arange(len(spectrum_x)):
    ph=(-ZPD_exact_shift)*math.pi*np.float64(j)/np.float64(len(spectrum_x))
#    ph = math.atan2(np.imag(spectrum_y_5[j]),np.real(spectrum_y_5[j]))
    phase_y.append(ph)

# phase correction and compensation of offset and slope
spectrum_y_3_ps=phase_correct_powspec(spectrum_y_2_ap) # power_spectrum
spectrum_y_3_ap=phase_correct_atan(spectrum_y_2_ap) # with apodization
spectrum_y_3_no=phase_correct_atan(spectrum_y_2_no) # no apodization
spectrum_y_3_sm=phase_correct_atan(spectrum_y_2_sm) # smoothed: no high harmonics

offset_compensate(spectrum_y_3_ap)
offset_compensate(spectrum_y_3_no)
offset_compensate(spectrum_y_3_sm)
for_search_y = np.real(spectrum_y_3_sm)-np.real(spectrum_y_3_ap)

# peak detection
max_table = max_detect_incr(spectrum_x,for_search_y,r_branch_min,r_branch_max,lines_number)
r_branch = []
for i in np.arange(len(max_table)):
    ind = int(np.array(max_table)[i,0])
    r_branch.append((spectrum_x[ind],spectrum_y_3_ap[ind]))

max_table = max_detect_decr(spectrum_x,for_search_y,p_branch_min,p_branch_max,lines_number)
p_branch = []
p_branch.append((0,0))
for i in np.arange(len(max_table)):
    ind = int(np.array(max_table)[i,0])
    p_branch.append((spectrum_x[ind],spectrum_y_3_ap[ind]))

# file output and plotting of the results
plt.figure(figsize=(14, 6))
plt.xlabel('step',size=10)
plt.ylabel('intensity, arb. un.',size=10)
plt.xlim(0,ZPD_position)
plt.plot(signal_x_range,signal_y_2_no)
plt.plot(signal_x_range,signal_y_2_ap)
plt.grid()
plt.show()

plt.figure(figsize=(14, 6))
plt.xlabel('wavenumber, 1/cm',size=10)
plt.ylabel('spectral amplitude, arb. un.',size=10)
plt.xlim(lim_x_min,lim_x_max)
#plt.plot(spectrum_x,for_search_y,label='for peak detection')
plt.plot(spectrum_x,spectrum_y_3_ps,label='naive power spectrum',ls='-')
plt.plot(spectrum_x,spectrum_y_3_no,label='no apodization',ls='-')
plt.plot(spectrum_x,spectrum_y_3_ap,label='with apodization',ls='-')
#plt.plot(spectrum_x,spectrum_y_3_sm,label='high freq skipped',ls='--')
if (len(r_branch)>1):
    plt.scatter(np.array(p_branch)[:,0],np.array(p_branch)[:,1],label='P branch (if any)',color='red',marker='o')
    plt.scatter(np.array(r_branch)[:,0],np.array(r_branch)[:,1],label='R branch (if any)',color='green',marker='o')
plt.legend()
plt.grid()
plt.show()


d = open(OUT_FILENAME, 'w')
d.write('line_wavenumber,1/cm intensity,a.u.\n');
if (len(p_branch)>1): 
    for i in np.arange(len(p_branch)):
        out_line = '{:.8}\t'.format(np.array(p_branch)[i,0]) + '{:.6}\t'.format(np.array(p_branch)[i,1]) + '\n'
        d.write(out_line); print(np.array(p_branch)[i,0])                            
d.write('\n');
if (len(r_branch)>1): 
    for i in np.arange(len(r_branch)):
        out_line = '{:.8}\t'.format(np.array(r_branch)[i,0]) + '{:.6}\t'.format(np.array(r_branch)[i,1]) + '\n'
        d.write(out_line); print(np.array(r_branch)[i,0])                            
d.close()
