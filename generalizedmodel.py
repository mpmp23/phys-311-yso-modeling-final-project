import numpy as np
import math ## needed this for code to work on my version
import matplotlib.pyplot as plt
import re


#****************************** Disk Parameters ****************************************

#Determine the temperature of each star from the peak/lowest wavelength using wien's law
starTempgap=(2.9e-3)/(.22e-6)
print("HD 169142 Star temperature:",starTempgap)
starTempdebris=(2.9e-3)/(.2965e-6)
print("HD 32297 Star temperature:",starTempdebris)
starTempIRX=(2.9e-3)/(.21e-6)
print("HD 31293 Star temperature:",starTempIRX)

#define some values for the Planck function

def PlanckLaw(T, wave):
    c=1.1927*10**-16 #2hc^2
    d=.014394135  #hc/Kb
    if (d/(T*wave))>500:  #Put in this condition to avoid getting overflow of digits in exponential
        return 0
    else:
        B= ((c/wave)**5)*(math.exp(d/(T*wave)) - 1 )**-1
        return B

#**************************** HD 169142 MODEL *****************************************
Rsub169=0.5    #the sublimation radius
#****************************** DUST WALL PARAMETERS ****************************************

Rwall169 = .5 # distance away from star of wall
wallHeight169 = 1 # height of wall (IN WHICH UNITS???)
wallScale169 = .07
Twall169 = 1800*(Rwall169/Rsub169)**(-1/2) # temperature of disk wall (assuming uniformity for simplicity's sake)
Awall169 = 2 * Rwall169 * wallHeight169 # code in lab dropped pi constant, so I will too!
Fwall169 = [] # where we'll store the flux of the wall, intensity times Awall
sumWall169 = 0

#Define the power law, scaling for the disk emission, and number of radial steps from .07 AU to 200. AU, with step size .05 AU
tpow169=-0.5	#defines the power law for the disk temperature profile
Dscale169=0.1   #scales the strength of the disk emission relative to the star
Rstepsize169=0.05   #the width of each annulus in AU
# Rsub169=0.5    #the sublimation radius
Rend169=200  #ending the disk at 200AU
RstepNumber169= int(round( (Rend169-Rsub169)/Rstepsize169) ) #Gives integer number of steps from 0.7-200 AU

#Define numpy arrays for starting radius, area, and temperature of each annulus, then fill the arrays
R169=[] # range(RstepNumber)  # will be inner radius of annulus
Area169=[] #=range(RstepNumber-1) # area of annulus , since total lum  ~ I x area
Temp169=[] #=range(RstepNumber-1) #Average temperature of our annuli

Rdisk169 = 50
for j in range(0,RstepNumber169):
  R169.append(Rdisk169)
  Rdisk169= Rdisk169 + Rstepsize169

for j in range(0,RstepNumber169-1):
  Area169.append((R169[j+1]**2- R169[j]**2))
  Rmid169= (R169[j+1] + R169[j])/2.0
  Temp169.append(1800*(Rmid169/Rsub169)**(tpow169))

#***** FLUX FROM DISK *************

wstart169=10e-9  # starting wavelength in m
wend169=800e-6   #end wavelength in m
wstep169 =100e-9  #step for wavelength in m
wstepNumber169=int(round( (wend169-wstart169)/wstep169) )

wavePlot169=[] # save wavelengths in array for easy plotting
Fdisk169=[]  # will be the total disk luminosity per wavelength, summed up from all the annuli
Fstar169=[]  # the star's luminosity per wavelength

Sumstar169=0
Sumdisk169=0

#Fill the wavelength, disk luminosity, and star luminosity arrays
wave169=wstart169
for i in range(0, wstepNumber169):
    wavePlot169.append(wave169)
    F_annuli169=0
    Fstar169.append(PlanckLaw(starTempgap,wavePlot169[i]))
    Sumstar169=Fstar169[i] + Sumstar169

    Fwall169.append(PlanckLaw(Twall169,wavePlot169[i])*Awall169)
    sumWall169=Fwall169[i] + sumWall169

    # Fhalo169.append(PlanckLaw(Thalo169,wavePlot169[i])*Ahalo169)
    # sumHalo169=Fhalo169[i] + sumHalo169
    
    for j in range(0,RstepNumber169-1):
            F_annuli169=PlanckLaw(Temp169[j],wavePlot169[i])*Area169[j] + F_annuli169

    Fdisk169.append(F_annuli169)


    wave169= wave169 + wstep169
    Sumdisk169=Fdisk169[i]+Sumdisk169

#************************** NORMALIZATION ****************************************

# NORMALIZE ALL THE FLUXES and scale disk flux strength to DScale
Fstar169=np.array(Fstar169)/Sumstar169
Fwall169=wallScale169*np.array(Fwall169)/sumWall169
Fdisk169=Dscale169*np.array(Fdisk169)/Sumdisk169

Ftotal169= Fstar169 + Fdisk169 + Fwall169


#**************************** HD 32297 MODEL *****************************************
Rsub322=0.7    #the sublimation radius
#****************************** DUST WALL PARAMETERS ****************************************

Rwall322 = .7 # distance away from start of wall
wallScale322 = 0.001 # scaling factor for dust wall
wallHeight322 = 1 # height of wall (IN WHICH UNITS???)
Twall322 = 1800*(Rwall322/Rsub322)**(-1/2) # temperature of disk wall (assuming uniformity for simplicity's sake)
Awall322 = 2 * Rwall322 * wallHeight322 # code in lab dropped pi constant, so I will too!
Fwall322 = [] # where we'll store the flux of the wall, intensity times Awall
sumWall322 = 0


#Define the power law, scaling for the disk emission, and number of radial steps from .07 AU to 200. AU, with step size .05 AU
tpow322=-0.75	#defines the power law for the disk temperature profile
Dscale322=0.008   #scales the strength of the disk emission relative to the star
Rstepsize322=0.05   #the width of each annulus in AU
# Rsub322=0.7    #the sublimation radius
Rend322=70  #ending the disk at 70 AU
RstepNumber322= int(round( (Rend322-Rsub322)/Rstepsize322  ) ) #Gives integer number of steps from 0.7-200 AU

#Define numpy arrays for starting radius, area, and temperature of each annulus, then fill the arrays
R322=[] # range(RstepNumber)  # will be inner radius of annulus
Area322=[] #=range(RstepNumber-1) # area of annulus , since total lum  ~ I x area
Temp322=[] #=range(RstepNumber-1) #Average temperature of our annuli

Rdisk322 = 18
for j in range(0,RstepNumber322):
  R322.append(Rdisk322)
  Rdisk322= Rdisk322 + Rstepsize322

for j in range(0,RstepNumber322-1):
  Area322.append((R322[j+1]**2- R322[j]**2))
  Rmid322= (R322[j+1] + R322[j])/2.0
  Temp322.append(1800*(Rmid322/Rsub322)**(tpow322))

#***** FLUX FROM DISK *************

wstart322=100e-9  # starting wavelength in m
wend322=1350e-6   #end wavelength in m
wstep322 =100e-9  #step for wavelength in m
wstepNumber322=int(round( (wend322-wstart322)/wstep322) )

wavePlot322=[] # save wavelengths in array for easy plotting
Fdisk322=[]  # will be the total disk luminosity per wavelength, summed up from all the annuli
Fstar322=[]  # the star's luminosity per wavelength

Sumstar322=0
Sumdisk322=0

#Fill the wavelength, disk luminosity, and star luminosity arrays
wave322=wstart322
for i in range(0, wstepNumber322):
    wavePlot322.append(wave322)
    F_annuli322=0
    Fstar322.append(PlanckLaw(starTempdebris,wavePlot322[i]))
    Sumstar322=Fstar322[i] + Sumstar322

    Fwall322.append(PlanckLaw(Twall322,wavePlot322[i])*Awall322)
    sumWall322=Fwall322[i] + sumWall322

    #Fhalo322.append(PlanckLaw(Thalo322,wavePlot322[i])*Ahalo322)
    #sumHalo322=Fhalo322[i] + sumHalo322
    
    for j in range(0,RstepNumber322-1):
            F_annuli322=PlanckLaw(Temp322[j],wavePlot322[i])*Area322[j] + F_annuli322

    Fdisk322.append(F_annuli322)


    wave322= wave322 + wstep322
    Sumdisk322=Fdisk322[i]+Sumdisk322



#************************** NORMALIZATION ****************************************


# NORMALIZE ALL THE FLUXES and scale disk flux strength to DScale
Fstar322=np.array(Fstar322)/Sumstar322
Fwall322=wallScale322*np.array(Fwall322)/sumWall322
Fdisk322=Dscale322*np.array(Fdisk322)/Sumdisk322

Ftotal322= Fstar322 + Fdisk322 + Fwall322 #+Fhalo322 

#**************************** HD 31293 MODEL *****************************************
Rsub312=0.13    #the sublimation radius

#****************************** DUST WALL PARAMETERS ****************************************

Rwall312 = .13 # distance away from star of wall
wallHeight312 = 1 # height of wall (IN WHICH UNITS???)
Twall312 = 1800*(Rwall312/Rsub312)**(-1/2) # temperature of disk wall (assuming uniformity for simplicity's sake)
Awall312 = 2 * Rwall312 * wallHeight312 # code in lab dropped pi constant, so I will too!
Fwall312 = [] # where we'll store the flux of the wall, intensity times Awall
sumWall312 = 0

#****************************** DUST HALO PARAMETERS ****************************************

Rhalo312 = 50 # distance away from start of halo
Thalo312 = 1800*(Rhalo312/Rsub312)**(-1/2) # blackbody relation from boltzman equation
Hscale312 = 0.05 # fix the contribution of the halo relative to the rest of the SED
Ahalo312 = 4 * Rhalo312**2 # code in lab dropped pi constant, so I will too!
Fhalo312 = [] # where we'll store the flux of the halo, intensity times Awall
sumHalo312 = 0

#***DISK PARAMETERS!***

#Define the power law, scaling for the disk emission, and number of radial steps from .07 AU to 200. AU, with step size .05 AU
tpow312=-0.75	#defines the power law for the disk temperature profile
Dscale312=0.1   #scales the strength of the disk emission relative to the star
Rstepsize312=0.05   #the width of each annulus in AU
# Rsub312=0.5    #the sublimation radius
Rend312=.5  #ending the disk in au
RstepNumber312= int(round( (Rend312-Rsub312)/Rstepsize312  ) ) #Gives integer number of steps from 0.7-200 AU

#Define numpy arrays for starting radius, area, and temperature of each annulus, then fill the arrays
R312=[] # range(RstepNumber)  # will be inner radius of annulus
Area312=[] #=range(RstepNumber-1) # area of annulus , since total lum  ~ I x area
Temp312=[] #=range(RstepNumber-1) #Average temperature of our annuli

Rdisk312 = 0.13
for j in range(0,RstepNumber312):
  R312.append(Rdisk312)
  Rdisk312= Rdisk312 + Rstepsize312

for j in range(0,RstepNumber312-1):
  Area312.append((R312[j+1]**2- R312[j]**2))
  Rmid312= (R312[j+1] + R312[j])/2.0
  Temp312.append(1800*(Rmid312/Rsub312)**(tpow312))

#***** FLUX FROM DISK *************

wstart312=10e-9  # starting wavelength in m
wend312=900e-6   #end wavelength in m
wstep312 =100e-9  #step for wavelength in m
wstepNumber312=int(round( (wend312-wstart312)/wstep312) )

wavePlot312=[] # save wavelengths in array for easy plotting
Fdisk312=[]  # will be the total disk luminosity per wavelength, summed up from all the annuli
Fstar312=[]  # the star's luminosity per wavelength

Sumstar312=0
Sumdisk312=0

#Fill the wavelength, disk luminosity, and star luminosity arrays
wave312=wstart312
for i in range(0, wstepNumber312):
    wavePlot312.append(wave312)
    F_annuli312=0
    Fstar312.append(PlanckLaw(starTempIRX,wavePlot312[i]))
    Sumstar312=Fstar312[i] + Sumstar312

    Fwall312.append(PlanckLaw(Twall312,wavePlot312[i])*Awall312)
    sumWall312=Fwall312[i] + sumWall312

    Fhalo312.append(PlanckLaw(Thalo312,wavePlot312[i])*Ahalo312)
    sumHalo312=Fhalo312[i] + sumHalo312
    
    for j in range(0,RstepNumber312-1):
            F_annuli312=PlanckLaw(Temp312[j],wavePlot312[i])*Area312[j] + F_annuli312

    Fdisk312.append(F_annuli312)


    wave312= wave312 + wstep312
    Sumdisk312=Fdisk312[i]+Sumdisk312



#************************** NORMALIZATION ****************************************


# NORMALIZE ALL THE FLUXES and scale disk flux strength to DScale
Fstar312=np.array(Fstar312)/Sumstar312
# Fwall312=np.array(Fwall312)/sumWall312
Fhalo312=Hscale312*np.array(Fhalo312)/sumHalo312
Fdisk312=Dscale312*np.array(Fdisk312)/Sumdisk312

Ftotal312= Fstar312 + Fhalo312 + Fdisk312 #+Fwall312 


#****************************** DATA FOR OBJECTS ********************************************

##### DATA POINTS FOR HD169142, wavelength in microns, lFl in ergs/s/cm^2
ysowithGap_L=[2.50E-01,0.22009045,0.35311225,0.4312245,0.64901257,1.230097,1.6670855,2.2125103,3.7765126,4.755935,11.851916,24.219994,59.08247,101.90737,455.75897,793.27783]
ysowithGap_F=[6.83E-09,8.19E-09,1.37E-08,2.25E-08,1.25E-08,4.74E-09,3.01E-09,2.05E-09,1.02E-09,7.10E-10,6.92E-10,2.18E-09,1.41E-09,6.71E-10,2.13E-11,1.95E-12]

##### DATA POINTS FOR HD32297, wavelength in microns, lFl in ergs/s/cm^2
debrisDisk_L=[0.2965441,0.37621716,0.8845664,1.2148691,1.5825667,2.5026786,3.5292904,6.2038813,8.90433,9.554837,12.893381,14.846047,18.182312,19.683344,60.82065,79.22889,145.54538,1318.0631]
debrisDisk_F=[2.00712E-11,1.63465E-11,4.54483E-12,2.2844E-12,1.15893E-12,3.34338E-13,1.2877E-13,3.19945E-14,1.67904E-14,1.72601E-14,1.64352E-14,1.77503E-14,3.54268E-14,3.41836E-14,4.18196E-14,3.10704E-14,9.28702E-15,1.16246E-17]

##### DATA POINTS FOR HD31293, wavelength in microns, lFl in ergs/s/cm^2
ysowithIRX_L=[0.21083546,0.3936139,0.4409265,0.54000074,0.5761869,1.1760412,1.6531409,2.177851,3.8105347,4.7430353,12.148064,24.795118,59.03734,97.59703,345.7263,440.92648,740.8296,803.3958,836.6334]
ysowithIRX_F=[6.88E-04,0.001624276,0.001568308,0.001066397,9.43E-04,3.13E-04,3.53E-04,3.72E-04,2.36E-04,1.72E-04,1.30E-04,1.07E-04,1.02E-04,6.57E-05,1.54E-06,4.93E-07,4.02E-08,3.55E-08,2.37E-08]


##need to convert to numpy array in order to be able to multiply by the scaling factor!
ysowithGap_F = np.array(ysowithGap_F)
ysowithGap_L = np.array(ysowithGap_L)
debrisDisk_F = np.array(debrisDisk_F)
debrisDisk_L = np.array(debrisDisk_L)
ysowithIRX_F = np.array(ysowithIRX_F)
ysowithIRX_L = np.array(ysowithIRX_L)

#****************************** CREATING A SCALING FACTOR FOR OBJECTS ********************************************

# C = Ftotal/GMAur_F <- scaling factor we're looking for based on lab - but, each object will have a different wavelength to normalize around!

####### YSO WITH GAP FIRST - HD169142 ########
# we're scaling this to around 1 micron
    #print(wavePlot[10]) <- this gave me 1.01 micrometers! right on target! so Ftotal[10] will give us what we want
    #ysowithGap_L_microns[5] is 1.23, but that's close enough for our purposes. if not, we can re-adjust the model.

Cgap = Ftotal169[10]/ysowithGap_F[5]

ysowithGap_F= Cgap*ysowithGap_F #Scale the fluxes to match model SED fluxes at ~1 micron

####### DEBRIS DISK - HD32297 ########
#modeling this also to 1 micron
#print(debrisDisk_L[3]) <- gave me 1.21 microns! close enough.

Cdebris = Ftotal322[10]/debrisDisk_F[3]

debrisDisk_F= Cdebris*debrisDisk_F #Scale the fluxes to match model SED fluxes at ~1 micron

####### YSO WITH IR EXCESS - HD31293 ########

#IR excess makes me want to model this to 0.5 microns - .21
print(ysowithIRX_L[4]) #<- gave me 2 microns! as good as we're going to get..
print(ysowithIRX_L[7]) #3.8 microns
print(wavePlot312[7]) #<- .41 microns! perfect.

Cirx = Ftotal312[7]/ysowithIRX_F[4] # 1 microns
# 4 2 -.44
# 20 7 - 2 microns
# 10 5 1 micron

ysowithIRX_F= Cirx*ysowithIRX_F #Scale the fluxes to match model SED fluxes at ~.5 microns

#****************************** MODELING EACH OBJECT WITH AN SED ********************************************
#****************************** DERIVING INFORMATION ABOUT EACH SYSTEM ********************************************

# LIR/L* (i.e. the ratio of the disk contribution to luminosity to the star’s luminosity. The value of LIR/L* is often correlated with the age of the object).


# Estimating the mass of the dust (or dust + gas disk). This is definitely easiest to do for HD 32297. To go about this, you’ll want to think about the surface area of emitting dust you need to fit the data, and then you can do something similar to what we did in HW #4 (for Saturn).



# Estimate L_IR/Lstar by integrating lambda F_lambda over d ln(lambda)

def luminosity_ratio(wave_microns, LFstar, *LF_ir_components):
    logwave = np.log(wave_microns)
    Lstar_est = np.trapezoid(LFstar, x=logwave)
    LIR_est = 0
    for comp in LF_ir_components:
        LIR_est += np.trapezoid(comp, x=logwave)
    return LIR_est / Lstar_est

LIR_Lstar_169 = luminosity_ratio(wavePlot169, Fstar169, Fdisk169, Fwall169)
LIR_Lstar_322 = luminosity_ratio(wavePlot322, Fstar322, Fdisk322, Fwall322)
LIR_Lstar_312 = luminosity_ratio(wavePlot312, Fstar312, Fdisk312, Fhalo312)

print("HD 169142 LIR/Lstar:", LIR_Lstar_169)
print("HD 32297 LIR/Lstar:", LIR_Lstar_322)
print("HD 31293 LIR/Lstar:", LIR_Lstar_312)

#dust mass estimate:
sigma = 5.67e-8  # W m^-2 K^-4

def dust_mass_from_LIR(LIR_Lstar, Lstar_W, Tdust_K, grain_radius_m=1e-6, grain_density=3000):
    Ldust = LIR_Lstar * Lstar_W
    A_emit = Ldust / (sigma * Tdust_K**4)
    Mdust = A_emit * grain_radius_m * grain_density / 3
    return Mdust, A_emit

# Example for HD 32297
Lsun = 3.828e26

# You need to estimate/assume Lstar for HD 32297.
# If you do not know it, use your stellar model or an approximate spectral type luminosity.
Lstar_322 = 5 * Lsun      # placeholder example
Tdust_322 = 50            # K, rough if peak is around 60 microns

Mdust_322, Aemit_322 = dust_mass_from_LIR(
    LIR_Lstar_322,
    Lstar_322,
    Tdust_322,
    grain_radius_m=1e-6,
    grain_density=3000
)

print("HD 32297 emitting area:", Aemit_322, "m^2")
print("HD 32297 dust mass:", Mdust_322, "kg")
print("HD 32297 dust mass in Earth masses:", Mdust_322 / 5.97e24)

## using our D/H/W scale values to come up with percentages of contribution!
print("HD 169142 disk contribution percent compared to Fstar: ", Dscale169*100,"%")
print("HD 169142 dust wall contribution percent compared to Fstar: ", wallScale169*100,"%")
print("HD 32297 disk contribution percent compared to Fstar: ", Dscale322*100,"%")
print("HD 32297 dust wall contribution percent compared to Fstar: ", wallScale322*100,"%")
print("HD 31293 disk contribution percent compared to Fstar: ", Dscale312*100,"%")
print("HD 31293 halo contribution percent compared to Fstar: ", Hscale312*100,"%")

# Relative integrated flux strengths!! derived from the scaling factor of each model.

IR_ratio_169 = Dscale169 + wallScale169
IR_ratio_322 = Dscale322 + wallScale322
IR_ratio_312 = Dscale312 + Hscale312

print("HD 169142 IR/star ratio:", IR_ratio_169)
print("HD 32297 IR/star ratio:", IR_ratio_322)
print("HD 31293 IR/star ratio:", IR_ratio_312)

print("HD 169142 / HD 32297 IR ratio:", IR_ratio_169 / IR_ratio_322)
print("HD 31293 / HD 32297 IR ratio:", IR_ratio_312 / IR_ratio_322)
print("HD 169142 / HD 31293 IR ratio:", IR_ratio_169 / IR_ratio_312)








#****************************** MODELING EACH OBJECT WITH AN SED ********************************************

# Convert wavelengths to microns and calculate lambda F_lambda for each object

wavePlot169 = np.array(wavePlot169) * 10**6
LFdisk169 = np.array(Fdisk169) * wavePlot169
LFstar169 = np.array(Fstar169) * wavePlot169
LFwall169 = np.array(Fwall169) * wavePlot169
# LFhalo169 = np.array(Fhalo169) * wavePlot169
LFtotal169 = np.array(Ftotal169) * wavePlot169

wavePlot322 = np.array(wavePlot322) * 10**6
LFdisk322 = np.array(Fdisk322) * wavePlot322
LFstar322 = np.array(Fstar322) * wavePlot322
LFwall322 = np.array(Fwall322) * wavePlot322
LFtotal322 = np.array(Ftotal322) * wavePlot322

wavePlot312 = np.array(wavePlot312) * 10**6
LFdisk312 = np.array(Fdisk312) * wavePlot312
LFstar312 = np.array(Fstar312) * wavePlot312
# LFwall312 = np.array(Fwall312) * wavePlot312
LFhalo312 = np.array(Fhalo312) * wavePlot312
LFtotal312 = np.array(Ftotal312) * wavePlot312


# Creates one figure with 3 separate plots stacked vertically
fig, axes = plt.subplots(3, 1, figsize=(10, 15), sharex=True)


#**************************** PLOT 1: HD 169142 *****************************************

ax = axes[0]

ax.plot(wavePlot169, np.log10(LFdisk169), label="flared disk")
ax.plot(wavePlot169, np.log10(LFstar169), label="star")
ax.plot(wavePlot169, np.log10(LFwall169), label="wall")
# ax.plot(wavePlot169, np.log10(LFhalo169), label="halo")
ax.plot(wavePlot169, np.log10(LFtotal169), label="total model")
ax.plot(ysowithGap_L, np.log10(ysowithGap_F), marker="o", linestyle="none", label="HD169142 data")

ax.set_title("HD 169142: YSO with Gap")
ax.set_xscale("log")
ax.set_xlim([0.1, 2000])
ax.set_ylim([-10, 0])
ax.set_ylabel("log(lambda F_lambda)")
ax.legend(loc="best", prop={"size": 8})
ax.grid(alpha=0.25, which="both")


#**************************** PLOT 2: HD 32297 *****************************************

ax = axes[1]

ax.plot(wavePlot322, np.log10(LFdisk322), label="flat disk")
ax.plot(wavePlot322, np.log10(LFstar322), label="star")
ax.plot(wavePlot322, np.log10(LFwall322), label="wall")
ax.plot(wavePlot322, np.log10(LFtotal322), label="total model")
ax.plot(debrisDisk_L, np.log10(debrisDisk_F), marker="o", linestyle="none", label="HD32297 data")

ax.set_title("HD 32297: Debris Disk")
ax.set_xscale("log")
ax.set_xlim([0.1, 2000])
ax.set_ylim([-10, 0])
ax.set_ylabel("log(lambda F_lambda)")
ax.legend(loc="best", prop={"size": 8})
ax.grid(alpha=0.25, which="both")


#**************************** PLOT 3: HD 31293 *****************************************

ax = axes[2]

ax.plot(wavePlot312, np.log10(LFdisk312), label="flared disk")
ax.plot(wavePlot312, np.log10(LFstar312), label="star")
# ax.plot(wavePlot312, np.log10(LFwall312), label="wall")
ax.plot(wavePlot312, np.log10(LFhalo312), label="halo")
ax.plot(wavePlot312, np.log10(LFtotal312), label="total model")
ax.plot(ysowithIRX_L, np.log10(ysowithIRX_F), marker="o", linestyle="none", label="HD31293 data")

ax.set_title("HD 31293: YSO with IR Excess")
ax.set_xscale("log")
ax.set_xlim([0.1, 2000])
ax.set_ylim([-10, 0])
ax.set_xlabel("wavelength (microns)")
ax.set_ylabel("log(lambda F_lambda)")
ax.legend(loc="best", prop={"size": 8})
ax.grid(alpha=0.25, which="both")


# Make spacing look better
fig.suptitle("SED Models for Three Disk Systems", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()
