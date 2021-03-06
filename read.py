import numpy as np
from pylab import *
from scipy import *
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit
import math
from lmfit import Model
import glob,os

c = 299792458 #m/s


#obs = 'digit'
obs = 'will'
lab_wvl = 63.18

if obs == 'digit':
	output = open('OI_digit.dat','w')
	output2 = open('OI_digit_short.dat','w')
	fwhm_start = 63.00/1601.0513  #digit
	x1 = 62.4
	x2 = 63.8
if obs == 'will':	
	output = open('OI_will.dat','w')
	output2 = open('OI_will_short.dat','w')
	fwhm_start = 63.00/3405.0922	#will
	x1 = 63.1
	x2 = 63.5


print >> output,'#|  name                    |   F_gauss   |   F_total  |    3*sigma   |     L_wing   |    R_wing    |  tag  |'
print >> output,'#|   all                    :    W/cm-2'




#gaussian function defined
def gaussian(x, amp, cen, fwhm, level):
    return level + amp * exp(-(x-cen)**2 /(2*(fwhm/2.35482)**2))
gmod = Model(gaussian)




fitslist=glob.glob('*.txt') #search current catalog for fits files
fitslist.sort()

for f in fitslist:
	root, ext = os.path.splitext(f)
	name=str(os.path.splitext(root)[0])
#uploading file
	x, y1, y2, y3 = np.loadtxt(name+'.txt', unpack=True, skiprows=1)
# x-wavelength, y1-flux, y2-continuum, y3-normalized
#np.isnan(x)==False









#choosing only OI line
	mask = (x[:] > x1) & (x[:] < x2) 
	x = x[mask]
	y = y3[mask]
#remove repeating wavelengths
	uniq = np.unique(x,return_index=True)
	x = x[uniq[1]]
	y = y[uniq[1]]
 



#calculating instrumental fwhm, 
#here not automatic, just taken from the file









#channel size
	dist = np.empty(600)
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
	vel_50 = 60*lab_wvl/(c*1e-3)
	print '200 km/s limit:', vel_limit


#setting gauss parameters
	gmod.set_param_hint('amp',value=max(y), vary=False) #careful???!!!
	gmod.set_param_hint('cen',value=63.18)
	gmod.set_param_hint('level',value=0.0, vary=False) 
	gmod.set_param_hint('fwhm', value = fwhm_start, vary=False) #fixed
	gmod.make_params()
#parameters with fixed fwhm and very constrained amplitude
#basically we let the fitting to play only with center


#fitting
	result = gmod.fit(y,x=x)
#optimal parameters
	amp = result.best_values['amp']
	cen = result.best_values['cen']
	fwhm = result.best_values['fwhm']
	level = result.best_values['level']
	sigma = fwhm/2.35482


#calculate field under the gauss function
	gaussfield = amp*math.sqrt(2*math.pi*sigma**2)*1e-4*c/(cen**2) 


#residual of data-gauss
	residual=y-result.best_fit



#right wing
	mask_right = ((x[:] < cen+vel_limit) & (x[:] > cen+vel_50))
	sum_right = sum(residual[mask_right])*chan_size*1e-4*c/(cen**2) #Jy*m

#left wing
	mask_left = ((x[:] > cen-vel_limit) & (x[:] < cen-vel_50))
	sum_left = sum(residual[mask_left])*chan_size*1e-4*c/(cen**2) #Jy*m

#integrated channels sum
	mask_sum = ((x[:] > cen-vel_limit) & (x[:] < cen+vel_limit))
	sum_sum = sum(y[mask_sum])*chan_size*1e-4*c/(cen**2) #jednostka???

#rms
	mask_rms = (x[:] > cen+vel_limit) | (x[:] < cen-vel_limit)
	std = np.std(y[mask_rms])
	y_rms = y[mask_rms]
	x_rms = x[mask_rms]
	mask_rms2 = (y_rms < 2.5*std)
	std = np.std(y_rms[mask_rms2])
	std_integrated = std*chan_size*1e-4*c/(cen**2)
	#should be multiplied by number of channels


#int_std = std*len1/3.0*chan_size*1e-6 
#int_std = int_std*1e-26*c/(1e-6*lab_wvl)**2 #W/m-2
#int_std = int_std*1e-4 #W/cm-2
				
#for smooth gauss on the plot
	xx = np.linspace(x1,x2,num=1000)

	tag = ''
	if 3*std_integrated<sum_left: 
		tag = 'L'
		clrl = 'r'
		clrr = 'k'
	if 3*std_integrated<sum_right:
		tag = 'R'
		clrr = 'r'
		clrl = 'k'
	if 3*std_integrated<sum_left and 3*std_integrated<sum_right: 
		tag = 'LR'
		clrr = 'r'
		clrl = 'r'
	if 3*std_integrated>sum_left and 3*std_integrated>sum_right: 
		tag = '-'
		clrr = 'k'
		clrl = 'k'

		



#===============================#
#=========plotti

	fig = plt.figure()
	ax = fig.add_subplot(111)


	ax.set_xlabel(r'$\lambda \ (\mu m)$',size="x-large")
	ax.set_ylabel(r'$F_\lambda\ (Jy)$',size="x-large")


	bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)

	ax.text(0.7,0.5,'3$\sigma$ = '+str(3*std_integrated)[0:8], fontsize=10,transform=ax.transAxes,bbox=bbox_props)
	ax.text(0.7,0.55,'right = '+str(sum_right)[0:8], fontsize=10,transform=ax.transAxes,bbox=bbox_props,color=clrr)
	ax.text(0.7,0.6,'left = '+str(sum_left)[0:8], fontsize=10,transform=ax.transAxes,bbox=bbox_props,color=clrl)
	ax.text(0.7,0.65,'gauss = '+str(gaussfield)[0:8], fontsize=10,transform=ax.transAxes,bbox=bbox_props)
	ax.text(0.7,0.7,'sum = '+str(sum_sum)[0:8], fontsize=10,transform=ax.transAxes,bbox=bbox_props)

	ax.plot(x_rms[mask_rms2],y_rms[mask_rms2],'ko')
	ax.plot(x[mask_left],residual[mask_left],'bo')
	ax.plot(x[mask_right],residual[mask_right],'co')
	ax.plot(x,y,drawstyle='steps-mid')
#ax.plot(x,result.init_fit)
	ax.plot(xx,gaussian(xx,amp,cen,fwhm,level))
#ax.plot(x,result.best_fit)
	ax.plot(x,residual,drawstyle='steps-mid')

	figname1 = name+'.eps'
	#fig.show()
	fig.savefig(figname1)
	fig.clf()
	plt.close()


	

	print >> output, "%25s    %10.6f    %10.6f    %10.6f    %10.6f    %10.6f    %6s" %(name,gaussfield, sum_sum, 3*std_integrated, sum_left, sum_right,tag)
	print >> output2, "%25s %6s" %(name, tag)

output.close()
