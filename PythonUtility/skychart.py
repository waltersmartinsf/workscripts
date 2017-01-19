from astropy import units as u #units package
from astropy.coordinates import SkyCoord #coordinate system package
from astropy.coordinates import EarthLocation, AltAz #obtain a Earth Location Coordinate
from astropy.coordinates import Longitude, Latitude
from astropy.coordinates import Angle #work with angles
from astropy.coordinates import get_sun,get_moon #Sun and Moon Location on Sky
from astropy.time import Time
import pandas as pd #dataframe
import numpy as np
import yaml
from string import split
import os
from matplotlib import pyplot as plt #graphical package



def skychart(table,site,utc_offset,midpointJD_label,ingressJD_label, egressJD_label, RA_deg_label,DEC_deg_label,planet_label,savepath):
    """
    Create a skychart using table dataframe with information of the ingress, egress, RA and DEC of the transit.
    ___
    INPUT:
    
    OUTPUT:
    
    """

    #Start a graphical size for the plots 
    def init_plotting():
        plt.rcParams['figure.figsize'] = (14.0,8.0)
        plt.rcParams['font.size'] = 20
        #plt.rcParams['font.family'] = 'Times New Roman'
        plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
        plt.rcParams['axes.titlesize'] = 0.75*plt.rcParams['font.size']
        plt.rcParams['legend.fontsize'] = 0.65*plt.rcParams['font.size']
        plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
        plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
        plt.rcParams['xtick.major.size'] = 3
        plt.rcParams['xtick.minor.size'] = 3
        plt.rcParams['xtick.major.width'] = 1
        plt.rcParams['xtick.minor.width'] = 1
        plt.rcParams['ytick.major.size'] = 3
        plt.rcParams['ytick.minor.size'] = 3
        plt.rcParams['ytick.major.width'] = 1
        plt.rcParams['ytick.minor.width'] = 1
        plt.rcParams['legend.frameon'] = True
        plt.rcParams['legend.loc'] = 'best'
        plt.rcParams['axes.linewidth'] = 1
    init_plotting()

    #if skychart directory exist, continue; if not, create.
    if not os.path.exists(savepath):
        os.makedirs(savepath)
    if not os.path.exists(savepath+'/skychart'):
        os.makedirs(savepath+'/skychart')
        
    #delta Time array 
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    
    #Remove NAN missing values in ingressJD, RA_deg and DEC_deg
    table = table[pd.isnull(table[midpointJD_label])==False]
    table = table[pd.isnull(table[ingressJD_label])==False]
    table = table[pd.isnull(table[egressJD_label])==False]
    
    #Number of objects
    N = len(table)
    print 'Number of exoplanets = ',N
    
    # observable = []
    for i  in range(N):
        midpoint = Time(Time(table[midpointJD_label].values[i],format='jd',scale='utc'))
        ingress2 = Time(Time(table[ingressJD_label].values[i],format='jd',scale='utc')) - 2*u.hour
        egress2 = Time(Time(table[egressJD_label].values[i],format='jd',scale='utc')) + 2*u.hour
        ingress = Time(Time(table[ingressJD_label].values[i],format='jd',scale='utc'))
        egress = Time(Time(table[egressJD_label].values[i],format='jd',scale='utc'))
        
        calendar = midpoint + np.linspace(-12, 12, 1000)*u.hour
        planet = SkyCoord(ra=table[RA_deg_label].values[i]*u.degree, dec=table[DEC_deg_label].values[i]*u.degree)
        planet_altazs = planet.transform_to(AltAz(obstime=calendar, location=site))
        altazframe = AltAz(obstime=calendar, location=site)
        sunaltazs = get_sun(calendar).transform_to(altazframe)
        moonaltazs = get_moon(calendar).transform_to(altazframe) #get moon coordiantes in sky
        
        info_exoplanet = pd.DataFrame([calendar.value,sunaltazs.alt.value,])
        
        #Time Criterium for Observable Exoplanet:
        #Transit (Ingress -2h ) > MAX(calendar[sunaltazs.alt > 0 and calendar < egress])
        # transit_criterium = calendar.value[(sunaltazs.alt > 0) & (calendar < egress2)]
        # crit1 = ingress2.jd > transit_criterium.max()
        # transit_criterium = calendar.value[(sunaltazs.alt > 0) & (calendar > egress2)]
        # crit2 = egress2.jd < transit_criterium.min()
        # if crit1 == True:
        #     if crit2 == True:
        #         observable.append(True)
        # else:
        #     observable.append(False)
        
        #Create the plot skychart
        plt.figure()
        plt.grid()
        plt.plot(calendar.value, sunaltazs.alt, color='y', label='Sun')
        plt.plot(calendar.value, moonaltazs.alt, 'o',label = 'Moon',color='green')  
        plt.scatter(calendar.value, planet_altazs.alt, c=planet_altazs.az, label=table[planet_label].values[i], lw=0, s=8)
        plt.vlines(midpoint.jd,0,90,color='lightgreen',label='Midpoint Transit')
        plt.vlines(ingress2.jd,0,90,color='red',label='Ingress - 2h')
        plt.vlines(egress2.jd,0,90,color='purple',label='Egress + 2h')
        plt.vlines(ingress.jd,0,90,color='red',label='Ingress',linestyles='--')
        plt.vlines(egress.jd,0,90,color='purple',label='Egress',linestyles='--')
        # plt.title('Ingress ='+str(ingress2.iso)+' Egress ='+str(egress2.iso))
        ##### PRINT INFORMATION ABOUT NIGHT, INGRESS, and EGRESS in UTC and LOCAL TIME
        local_time_ingress = (Time(table[ingressJD_label].values[i],format='jd',scale='utc') + utc_offset).iso
        local_time_egress = (Time(table[ingressJD_label].values[i],format='jd',scale='utc') + utc_offset).iso
        #Choose the correct night for asking time
        # start_time = int((Time(table[ingressJD_label].values[2],format='jd',scale='utc') + utc_offset).iso[11:-10])
        # if (start_time < 24) and (start_time >= 18):
        #     night_date = (Time(table[ingressJD_label].values[2],format='jd',scale='utc') + utc_offset).iso[:10]
        #     # print night_date
        # if (start_time < 8) and (start_time >= 0):
        #     night_date = (Time(table[ingressJD_label].values[2],format='jd',scale='utc') + utc_offset).iso[:10]
        #     night_date = Time(night_date) - 1*u.day
        #     night_date = night_date.iso[:10]
        #     # print night_date
        # plt.title('Site Night = '+night_date+'\n'+'Site Time: '+'Ingress ='+str(local_time_ingress)+' Egress ='+str(local_time_egress)+'\n'+'UTC: '+'Ingress ='+str((Time(table[ingressJD_label].values[i],format='jd',scale='utc')).iso)+' Egress ='+str((Time(table[ingressJD_label].values[i],format='jd',scale='utc')).iso))
        plt.title('Site Time: '+'Ingress ='+str(local_time_ingress)+' Egress ='+str(local_time_egress)+'\n'+'UTC: '+'Ingress ='+str((Time(table[ingressJD_label].values[i],format='jd',scale='utc')).iso)+' Egress ='+str((Time(table[ingressJD_label].values[i],format='jd',scale='utc')).iso))
        #####
        plt.fill_between(calendar.value, 0, 90, sunaltazs.alt < -0*u.deg, color='0.5', zorder=0)  
        plt.fill_between(calendar.value, 0, 90, sunaltazs.alt < -18*u.deg, color='k', zorder=0)  
        plt.colorbar().set_label('Azimuth [deg]')  
        plt.legend(loc='upper left')
        plt.xlim(calendar.value.min(),calendar.value.max()) 
        plt.ylim(0, 90)  
        plt.xlabel('Julian Date')  
        plt.ylabel('Altitude [deg]')
        plt.savefig(savepath+'/skychart/'+str(table[planet_label].values[i])+' '+str(midpoint.jd)+'_.png')
        plt.close()
    print 'Done \n'
    # return observable