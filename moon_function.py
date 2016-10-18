import fitsio
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, get_moon, EarthLocation
from astropy.time import Time
from astropy.io.fits.hdu.hdulist import HDUList

def get_moon_angle(file_name, number_of_files=100):
    df_SKY=pd.read_csv(file_name)
    new_SKY = df_SKY.query('OBJTYPE=="SKY".ljust(16)')
    DEC = np.array(new_SKY['DEC'].tolist())
    FIBERID = np.array(new_SKY['FIBERID'].tolist())
    MJD = np.array(new_SKY['MJD'].tolist())
    RA = np.array(new_SKY['RA'].tolist())
    plates   = list(new_SKY['PLATE'][:number_of_files].values)
    mjds   = list(new_SKY['MJD'][:number_of_files].values)
    fiberids   = list(new_SKY['FIBERID'][:number_of_files].values)
    files = ['spec-%s-%s-%s.fits'%(plate,mjd,str(fiberid).zfill(4)) for plate,mjd,fiberid in zip(plates,mjds,fiberids)]
    #FINDING INTERPLATE SKY TIME
    fin_mean=[] #time for each moon sky observation interplate, used for moon angle
    for r in range(len(files)): #mechanization of ubiquitious exposure count
        sampling=fitsio.read_header(files[r],0) #finding number of exposures
        h_beg=[]
        h_end=[]
        mean_per=[] #time for each moon sky observation intraplate
        for k in range(4, sampling['NEXP']+4): #4 is constant for everything
            h=fitsio.read_header(files[3],k)
            h_beg.append(h['TAI-BEG']) #ONLY VARIES WITH PLATE NUMBER, different for k's
            h_end.append(h['TAI-END']) #ONLY VARIES WITH PLATE NUMBER, different for k's
            mean_times=(h['TAI-BEG']+h['TAI-END'])/2
            mean_per.append(mean_times)
        tot=np.mean(mean_per)
        fin_mean.append(tot)
    #TAI TO MJD
    new_time=np.array(fin_mean)
    time_MJD=new_time/(86400)
    moon_coords1=[]
    moon_coords=[]
    for i in range (len(time_MJD)):
        t=Time(time_MJD[i], format='mjd')
        moon_coords1.append(get_moon(t))
        moon_coords.append(moon_coords1[i].spherical)
    moon_coords=list(moon_coords)
    #print(moon_coords)
    moon_coords2=[]
    for i in range(len(moon_coords)):
        moon_coords2.append(moon_coords[i]._values)
    moon_coords2=np.array(moon_coords2)
    sky_coords=np.vstack((RA,DEC)).T
    #SPLITTING ARRAY
    new_sky_coords=sky_coords[:number_of_files]
    new_sky_RA = [item[0] for item in new_sky_coords]
    new_sky_DEC = [item[1] for item in new_sky_coords]
    moon_RA = [item[0] for item in moon_coords2]
    moon_DEC = [item[1] for item in moon_coords2]
    new_sky_RA=np.array(new_sky_RA)*(np.pi/180)
    new_sky_DEC=np.array(new_sky_DEC)*(np.pi/180)
    moon_RA=np.array(moon_RA)*(np.pi/180)
    moon_DEC=np.array(moon_DEC)*(np.pi/180)
    #MOON-SKY ANGLE CALCULATIONS
    moon_sky_angle=[]
    for i in range(len(new_sky_coords)):
        moon_sky_angle.append(np.arccos((np.sin(new_sky_RA[i])*np.sin(moon_RA[i]))+np.cos(new_sky_RA[i])*np.cos(moon_RA[i])*(np.cos(new_sky_DEC[i]-moon_DEC[i])))*(180/np.pi))
    #list(moon_sky_angle)
    return (np.array(moon_sky_angle))

