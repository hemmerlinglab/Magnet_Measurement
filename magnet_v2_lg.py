import numpy as np
import matplotlib.pyplot as plt
import lmfit
from lmfit import Minimizer, Parameters, fit_report

z_in = np.linspace(0,1,201)
z = z_in*25.4e-3

#uo = 1.25663706e-6 # N/A^2
prb = 0.056*25.4e-3
R = 0.75*25.4e-3
L = 0.75*25.4e-3
d = .015*25.4e-3 # difference between 
#a1 = (prb-d)/2
#a2 = (prb+d)/2
z1 = z + L
z2 = z + L

b1f = open('b1l.txt','r')
b2f = open('b2l.txt','r')
B1mT = np.genfromtxt(b1f)
B2mT = np.genfromtxt(b2f)
b1f.close()
b2f.close()

B1 = B1mT*10
B2 = B2mT*10
#print(len(B1),len(B2))
#plt.show()

#B(z) = uo*M/2(z/(sqrt(z^2+R^2))-(z-L)/sqrt((z-L)^2+R^2))
def fcn2min(params,x,data,plot_fit=False):
        M = params['M']
        #y_offset = params['y_offset']
        a = params['a']
        offset = -params['offset']
        if plot_fit == False:
            model = M*((x+a+offset)/(np.sqrt((x+a+offset)**2+R**2))-(x+a-L+offset)/np.sqrt((x+a-L+offset)**2+R**2))
            return model - data
        else:
            x_plot = np.linspace(np.min(x), np.max(x), 500)
            model = M*((x_plot+a+offset)/(np.sqrt((x_plot+a+offset)**2+R**2))-(x_plot+a-L+offset)/np.sqrt((x_plot+a-L+offset)**2+R**2))
            return (x_plot, model)




def fit_B(x,y):
	params = Parameters()
	params.add('M',value=8000,min=0,max=1e6,vary=True)
	params.add('a',value=prb,min=-1,max=1,vary=True)
	params.add('offset',value=0,min=-1,max=1,vary=False)
	#params.add('y_offset',value=0,min=-1e4,max=1e4,vary=True)
	minner = Minimizer(fcn2min,params,fcn_args=(x,y))
	result = minner.minimize()
	con_report = fit_report(result.params)
	(xplot,model) = fcn2min(result.params,x,y,plot_fit=True)
	return (xplot,model,result)




n = len(B1)
#print(n)
(x_fit1,y_fit1,result1) = fit_B(z1[0:n],B1[0:n])
(x_fit2,y_fit2,result2) = fit_B(z2[0:n],B2[0:n])

M1 = result1.params['M'].value
M2 = result2.params['M'].value
M = (M1+M2)/2
#y_off_1 = result1.params['y_offset'].value
#y_off_2 = result2.params['y_offset'].value
#print(y_off_1,y_off_2)
print(M1,M2)
#y_off = (y_off_1+y_off_2)/2
a1 = result1.params['a'].value
a2 = result2.params['a'].value
print('M = {}'.format(int(M)))
#print('y_off = {} Gs'.format(y_off))
print('a1 = {} mm\na2 = {} mm'.format(a1*10**3,a2*10**3))
print('a_diff = {} mm'.format(abs(a2-a1)*10**3))
print('Expecting: {} mm'.format(d*1e3))

plt.subplot(121)
plt.scatter(z1*1e3,B1,color='r')
plt.plot(x_fit1*1e3,y_fit1,linestyle='-')
plt.xlabel('z (mm)')
plt.ylabel('Magnetic Flux (Gs)')
plt.ylim(0,6000)
plt.xlim(np.min(z1*1e3),np.max(z1*1e3))
plt.title('Hole Away')

plt.subplot(122)
plt.scatter(z2*1e3,B2,color='r')
plt.plot(x_fit2*1e3,y_fit2,linestyle='-')
plt.xlabel('z (mm)')
plt.ylim(0,6000)
plt.xlim(np.min(z1*1e3),np.max(z1*1e3))
plt.title('Hole Towards')

#plt.figure()
#plt.plot(x_fit1,y_fit1)

plt.figure()


(xe, ye) = fcn2min(result1.params,np.linspace(0,0.1,100),1,plot_fit=True)

plt.plot((xe - L)*1e2,ye, '-')
plt.plot((z1 - L)*1e2,B1,'o')

plt.axhline(268, linestyle = '--', color = 'r')
plt.axvline(0, linestyle = '--')

plt.yscale('log')

plt.xlabel('Distance from magnet (cm)')
plt.ylabel('Magnetic Field (Gs)')


plt.figure()

dst_between_magnets = L+6.34e-2


xspace = np.linspace(-0.05,L+0.15,100)

(xe, ye) = fcn2min(result1.params,xspace,1,plot_fit=True)

result1.params['offset'].value = dst_between_magnets
(x2e, y2e) = fcn2min(result1.params,xspace,1,plot_fit=True)


plt.plot((x2e - L)*1e2, (ye + y2e), '-')

plt.axhline(356, linestyle = '--', color = 'r')
plt.text(15, 356, 'Yb 200 m/s')

plt.axhline(268, linestyle = '--', color = 'r')
plt.text(15, 268, 'Yb 150 m/s')

plt.axhline(179, linestyle = '--', color = 'r')
plt.text(15, 179, 'Yb 100 m/s')


plt.axvline(0, linestyle = '--', color = 'k')
plt.axvline(-L/1e-2, linestyle = '--', color = 'k')
plt.axvline((dst_between_magnets)/1e-2, linestyle = '--', color = 'k')
plt.axvline((dst_between_magnets-L)/1e-2, linestyle = '--', color = 'k')

plt.yscale('log')

plt.xlabel('Distance from magnet (cm)')
plt.ylabel('Magnetic Field (Gs)')


plt.show()



