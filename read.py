import numpy as np
from pylab import *
from scipy import *
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit
import math
from lmfit import Model

c = 299792458 #m/s




#gaussian function defined
def gaussian(x, amp, cen, fwhm):
    return amp * exp(-(x-cen)**2 /(2*(fwhm/2.35482)**2))
gmod = Model(gaussian)


#uploading file
x, y1, y2, y3 = np.loadtxt('will03_aqu3.txt', unpack=True, skiprows=1)
# x-wavelength, y1-flux, y2-continuum, y3-normalized
#np.isnan(x)==False



#choosing only OI line
mask = (x[:] > 62) & (x[:] < 65) 
x = x[mask]
y = y3[mask]
#remove repeating wavelengths
uniq = np.unique(x,return_index=True)
x = x[uniq[1]]
y = y[uniq[1]]
 



#calculating instrumental fwhm, 
#here not automatic, just taken from the file
fwhm = 63.00/3405.0922
lab_wvl = 63.18
print len(x)


#channel size
dist = np.empty(200)
for ii in range(len(x)-1):
	dist[ii] = x[ii+1]-x[ii]
dist = dist[0:ii]
chan_size = np.average(dist)


#velocity resolution and velocity limit
#300 km/s limit for high-velocity wings
#center +- 300km/s - area to measure flux
#>300km/s area to measure rms

print 'Channel size:', chan_size
vel_res = c*1e-3*chan_size/lab_wvl		
print 'Velocity resolution:',vel_res
vel_limit = 300*lab_wvl/(c*1e-3)
print '200 km/s limit:', vel_limit


#setting gauss parameters
gmod.set_param_hint('amp',value=max(y), min=max(y)-1.5, max=max(y)+1.5) #careful???!!!
gmod.set_param_hint('cen',value=63.18)
gmod.set_param_hint('fwhm', value = fwhm, vary=False) #fixed
gmod.make_params()
#parameters with fixed fwhm and very constrained amplitude
#basically we let the fitting to play only with center


#fitting
result = gmod.fit(y,x=x)
#optimal parameters
amp = result.best_values['amp']
cen = result.best_values['cen']
fwhm = result.best_values['fwhm']
sigma = fwhm/2.35482


#calculate field under the gauss function
gaussfield = amp*math.sqrt(2*math.pi*sigma)/1e6 #Jy*m


#residual of data-gauss
residual=y-result.best_fit

#rms
mask_rms = (x[:] > cen+vel_limit) | (x[:] < cen-vel_limit)
std = np.std(y[mask_rms])
std_integrated = std*len(y[mask_rms])/chan_size*1e-6 #Jy*m

#right wing
mask_right = ((x[:] < cen+vel_limit) & (x[:] > cen))
sum_right = sum(residual[mask_right])*len(residual[mask_right])/chan_size*1e-6 #Jy*m

#left wing
mask_left = ((x[:] > cen-vel_limit) & (x[:] < cen))
sum_left = sum(residual[mask_left])*len(residual[mask_left])/chan_size*1e-6 #Jy*m

#integrated channels sum
mask_sum = ((x[:] > cen-vel_limit) & (x[:] < cen+vel_limit))
sum_sum = sum(y[mask_sum])*len(y[mask_sum])/chan_size*1e-6 #Jy*m



#int_std = std*len1/3.0*chan_size*1e-6 
#int_std = int_std*1e-26*c/(1e-6*lab_wvl)**2 #W/m-2
#int_std = int_std*1e-4 #W/cm-2
				
#for smooth gauss on the plot
xx = np.linspace(63.1,63.5,num=1000)



#===============================#
#=========plotti

fig = plt.figure()
ax = fig.add_subplot(111)

bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)

ax.text(0.7,0.5,'3$\sigma$ = '+str(3*std_integrated)[0:8], fontsize=10,transform=ax.transAxes,bbox=bbox_props)
ax.text(0.7,0.55,'right = '+str(sum_right)[0:8], fontsize=10,transform=ax.transAxes,bbox=bbox_props)
ax.text(0.7,0.6,'left = '+str(sum_left)[0:8], fontsize=10,transform=ax.transAxes,bbox=bbox_props)
ax.text(0.7,0.65,'gauss = '+str(gaussfield)[0:8], fontsize=10,transform=ax.transAxes,bbox=bbox_props)
ax.text(0.7,0.7,'sum = '+str(sum_sum)[0:8], fontsize=10,transform=ax.transAxes,bbox=bbox_props)

ax.plot(x[mask_rms],y[mask_rms],'ko')
ax.plot(x[mask_left],residual[mask_left],'bo')
ax.plot(x[mask_right],residual[mask_right],'co')
ax.plot(x,y,drawstyle='steps-mid')
#ax.plot(x,result.init_fit)
ax.plot(xx,gaussian(xx,amp,cen,fwhm))
#ax.plot(x,result.best_fit)
ax.plot(x,residual,drawstyle='steps-mid')

figname1 = 'test1.png'
fig.show()
fig.savefig(figname1)
fig.clf()