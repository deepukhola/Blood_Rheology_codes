#!/usr/bin/env python
# coding: utf-8

# ## BLOOD RHEOLOGY CODES

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import matplotlib


# ### defining functions

# In[2]:


def fx(Hc = np.nan, T = np.nan, Ht = np.nan, x = np.nan):
    ac = 0.070 * np.exp(2.49 * Hc + 1107 * np.exp(-1.69 * Hc) / T)
    y1 = Ht / 0.4 - (Ht/Hc) * (1 + ((1 - Ht / Hc)**2)/((Ht / Hc) * (2 - 2 * Ht / Hc + Ht*(1 - ac*Hc) / Hc)))
    y2 = 1 / (1 - ac*Hc)
    y3 = Ht / 0.4 - (x ** 2) * (1 + (((1 - (x ** 2)) ** 2)/ ((x ** 2) * (2 - 2 * (x ** 2) + (x ** 2)/ y2))))
    af = 0.070 * np.exp(2.49 * 0.4 + 1107 * np.exp(-1.69 * 0.4) / T)
    y4 = (1-0.4*af) / (1 - ((x**4)*ac*Hc))
    return y1, y2, y3, y4

HT = np.linspace(0.335, 0.39, 10)    # HT data list
Temp = [293.15, 310.15]              # temperature values


# ### $H_c$ calculation using bisection method

# In[3]:


def root_range():
    return [0.3, 0.5]


# In[4]:



Hc_range = [0.3, 0.5]      # range of HC
a = Hc_range
Hc_values = []
count = 0
for j in Temp:
    Hc_values.append([])
    for i in HT:
        while True:
            # if guess value is one of the final root value or goes to infinity
            if fx(Hc = Hc_range[0], T = j, Ht = i)[0] == float('inf'):
                Hc_range[0] += sum(Hc_range)/0.000001
            if fx(Hc = Hc_range[0], T = j, Ht = i)[0] == 0:
                root = Hc_range[0]
                break
            if fx(Hc = Hc_range[1], T = j, Ht = i)[0] == float('inf'):
                Hc_range[1] -= sum(Hc_range)/0.000001
            if fx(Hc = Hc_range[1], T = j, Ht = i)[0] == 0:
                root = Hc_range[1]
                break

            #if Hc_range extremes lie on same side of x-axis
            if fx(Hc = Hc_range[0], T = j, Ht = i)[0] / fx(Hc = Hc_range[1], T = j, Ht = i)[0] > 0:
                print("enter Hc_range again")
                Hc_range = [float(input("enter initial value :")), float(input("enter final value : "))]
                a = Hc_range
                continue

            # if mid value is infinity 
            if fx(Hc = sum(Hc_range)/2, T = j, Ht = i)[0] == float('inf'):
                Hc_range[0] += sum(Hc_range)/0.000001
            
            # if evrything is fine
            if fx(Hc = sum(Hc_range)/2, T = j, Ht = i)[0] / fx(Hc = Hc_range[0], T = j, Ht = i)[0] > 0:
                Hc_range[0] = sum(Hc_range)/2
            if fx(Hc = sum(Hc_range)/2, T = j, Ht = i)[0] / fx(Hc = Hc_range[1], T = j, Ht = i)[0] > 0:
                Hc_range[1] = sum(Hc_range)/2
                
            
            # tolerence condition
            if Hc_range[1] - Hc_range[0] < 0.000001:
                root = sum(Hc_range)/2
                break
        Hc_values[count].append(root)
        Hc_range = root_range()                       # redefiening range for fresh iteration

    count += 1

plt.plot(HT, Hc_values[0], "r--", HT, Hc_values[1], "b--")
plt.legend(["T = 20°C","T = 37°C"])
plt.xlabel(r"$H_T$")
plt.ylabel(r"$H_C$")
plt.text(0.38, 0.425, "$H_F$ = 40%")
plt.show()


# ### $\frac{\mu_{c}}{\mu_a}$ Calculation

# In[5]:


# calculation of y2
count = 0
ucua = []
for i in Temp:
    ucua.append([])
    for j,k in zip(Hc_values[count], HT):
        ucua[count].append(fx(Hc = j, T = i, Ht = k, x = np.nan)[1])
    count += 1

plt.plot(HT, ucua[0], "r--", HT, ucua[1], "b--")
plt.legend(["T = 20°C","T = 37°C"])
plt.xlabel(r"$H_T$")
plt.ylabel(r'$\frac{\mu_{c}}{\mu_a}$')
plt.text(0.38, 2.13, "$H_F$ = 40%")
plt.show()


# ### $\sigma$ calculation using bisection method

# In[6]:


def x_root_range():
    return [0.5, 1]


# In[7]:


count = 0
x_range = x_root_range()
x_values = []
for i in Temp:
    x_values.append([]) 
    for j, k in zip(Hc_values[count], HT):
        while True:
            # check if any of the extreme is zero or infinity
            if fx(Hc = j, T = i, Ht = k, x = x_range[0])[2] == float("inf"):
                x_range[0] += (x_range[1]-x_range[0])/100000
                continue
            if fx(Hc = j, T = i, Ht = k, x = x_range[0])[2] == 0:
                root = x_range[0]
                break
            if fx(Hc = j, T = i, Ht = k, x = x_range[1])[2] == float("inf"):
                x_range[1] += (x_range[1]-x_range[0])/100000
                continue
            if fx(Hc = j, T = i, Ht = k, x = x_range[1])[2] == 0:
                root = x_range[1]
                break
                
            # check if both extremes lie on same side of x-axis
            if fx(Hc = j, T = i, Ht = k, x = x_range[0])[2] / fx(Hc = j, T = i, Ht = k, x = x_range[1])[2] > 0:
                print("enter right interval")
                x_range = [float(input("enter first value : ")), float(input("enter last value : "))]
                continue
                
                
            # check if function is infinity at mid value
            if fx(Hc = j, T = i, Ht = k, x = sum(x_range)/2)[2] == float("inf"):
                x_range[0] += (x_range[1]-x_range[0])/100000
                continue
                
                
            # if everything is fine
            if fx(Hc = j, T = i, Ht = k, x = x_range[0])[2] / fx(Hc = j, T = i, Ht = k, x = sum(x_range)/2)[2] > 0:
                x_range[0] = sum(x_range) / 2
            
            if fx(Hc = j, T = i, Ht = k, x = x_range[1])[2] / fx(Hc = j, T = i, Ht = k, x = sum(x_range)/2)[2] > 0:
                x_range[1] = sum(x_range) / 2
                
                
                
            # termination situation
            if x_range[1] - x_range[0] < 0.0000001:
                root = sum(x_range) / 2
                break
        x_values[count].append(root)
        x_range = x_root_range()
    count += 1
            
                
plt.plot(HT, x_values[0], "r--", HT, x_values[1], "b--")  
plt.legend(["T = 20°C","T = 37°C"])
plt.xlabel(r"$H_T$")
plt.ylabel(r'$\sigma=\frac{r*}{R}$')
plt.text(0.34, 0.96, "$H_F$ = 40%")
plt.show()
                


# ### $\frac{\mu_{app}}{\mu_F}$ calculation

# In[8]:


uappuf = []
count = 0
for i in Temp:
    uappuf.append([])
    for j, k, l in zip(x_values[count], Hc_values[0], HT):
        uappuf[count].append(fx(Hc = k, T = i, Ht = l, x = j)[3])
    count += 1

plt.plot(HT, uappuf[0], "r--", HT, uappuf[1], "b--")
plt.legend(["T = 20°C","T = 37°C"])
plt.xlabel(r"$H_T$")
plt.ylabel(r'$\frac{\mu_{app}}{\mu_F}$')
plt.text(0.34, 0.9, "$H_F$ = 40%")
plt.show()


# ## Bifurcation case

# ### $\eta , \lambda$ calculation

# In[9]:


X = np.linspace(145, 1000, 10000)
y1 = []
t = []
for x in X:
    y1.append((1+(1000/x)**4)**(-1))
    t.append(x)

x1 = t
l1 =[]

for i, j in zip(x1,y1):
    z = (-6.96/1000)*np.log(i/1000) + (1+6.98*(1-0.4)/1000) * np.log(((j-0.4/1000)/(1-2*0.4/1000))/(1-(j-0.4/1000)/(1-2*0.4/1000)))
    l1.append(1/(1+1/np.exp(z)))
plt.plot(x1, y1, "r-",x1, l1, "b-" )

plt.loglog()
plt.legend([r"$\eta$",r"$\lambda$"])
plt.xlabel(r"$D_1(\mu m)$")
plt.ylabel(r"$\eta , \lambda$")
plt.text(600, 0.0001, "Bifurcation Case \n D2=DF=1000μm \n $H_F$=40% \n 37$\degree$C")
plt.show()


# ### ${H_T}_i$ calculation

# In[10]:


HD = 0.4
m = []
n = []
for i, j in zip(l1, y1):
    m.append(1-i)
    n.append(1-j)
l = [l1, m]
y = [y1, n]
HT = []
x2 = []
for i in range(0, len(x1)):
    x2.append(1000)
x = [x1,x2]

for i in [0,1]:
    HT.append([])

    for j,k,p in zip(l[i],y[i],x[i]):
        Ht = (HD*j/k) * ((HD*j/k)+(1-(HD*j/k))*(1+(1.7*np.exp(-35*p/100))-0.6*(np.exp(-p/100))))
        HT[i].append(Ht) 


# ### $\sigma_1 , \sigma_2$ calculation

# In[11]:


def sigma(si):

    aci = 0.070*np.exp(2.49*hti/(si**2) + (1107/310)*np.exp(-1.69*hti/(si**2)))
    funs = hti/(0.4*li/yi) - (si**2)*(1 + ((1-(si**2))**2)/((si**2)*(2 - 2*si*si + si*si*(1-aci*hti/(si*si)))))
    return funs

sig = []
for i in [0,1]:
    sig.append([])
    for hti, li, yi in zip(HT[i], l[i], y[i]):
        
        a = fsolve(sigma, 0.9)
        sig[i].append(a)

plt.plot(x[0], sig[0], "r--", x[0], sig[1], "b--")
plt.xscale("log")
plt.legend([r"$\sigma_1$",r"$\sigma_2$"], loc = 1)
plt.xlabel(r"$D_1(\mu m)$")
plt.ylabel(r"$\sigma_1 , \sigma_2$")
plt.text(600, 0.92, "Bifurcation Case \n D2=DF=1000μm \n $H_F$=40% \n 37$\degree$C")
plt.show()


# ### $H_{C1} , H_{C2}$ calculation

# In[12]:


HC = []
for i in [0,1]:
    HC.append([])
    for sigi, hti in zip(sig[1], HT[i]):
        hci = hti/sigi**(0.5)
        HC[i].append(hci)
    
plt.plot(x[0], HC[0], "r--", x[0], HC[1], "b--")
plt.xscale("log")
plt.legend([r"$H_{C1}$",r"$H_{C2}$"], loc = 1)
plt.xlabel(r"$D_1(\mu m)$")
plt.ylabel(r"$H_{C1} , H_{C2}$")
plt.text(600, 0.05, "Bifurcation Case \n D2=DF=1000μm \n $H_F$=40% \n 37$\degree$C")
plt.show()

