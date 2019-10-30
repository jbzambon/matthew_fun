# Program to plot from the new CNAPS OPeNDAP server
#   using Hurricane Matthew as a guinea pig
#
# Joseph B Zambon
# 4 October 2016

#Dependencies
import pandas as pd
from pydap.client import open_url
import httplib2
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime
import numpy as np

# Put your desired start/end dates in here.
#                              YYYY,MM,DD,HH
start_date = datetime.datetime(2016,10, 4, 0)
end_date = datetime.datetime  (2016,10, 7, 0)

# For inline plotting
#get_ipython().magic(u'pylab inline')

# Link OPeNDAP datasets
wrf_cnaps_url = 'http://oceanus.meas.ncsu.edu:8080/thredds/dodsC/fmrc/useast_coawst_wrf/COAWST-WRF_Forecast_Model_Run_Collection_best.ncd'
wrf_dataset = open_url(wrf_cnaps_url)
print('Available WRF data:')
print(wrf_dataset.keys)
roms_cnaps_url = 'http://oceanus.meas.ncsu.edu:8080/thredds/dodsC/fmrc/useast_coawst_roms/COAWST-ROMS_SWAN_Forecast_Model_Run_Collection_best.ncd'
roms_dataset = open_url(roms_cnaps_url)
print('')
print('')
print('Available ROMS/SWAN data:')
print(roms_dataset.keys)

# Let's ingest the latitude and longitude values
wrf_lon=np.array(wrf_dataset['lon'])
wrf_lat=np.array(wrf_dataset['lat'])
roms_lon=np.array(roms_dataset['lon_rho'])
roms_lat=np.array(roms_dataset['lat_rho'])
roms_mask=np.array(roms_dataset['mask_rho'])

#Find WRF time indices
wrf_origin_date = datetime.datetime(2016,9,20,0,0,0)
wrf_time=(np.array(wrf_dataset['time'][:])/24)+datetime.date.toordinal(wrf_origin_date)
wrf_start_index = np.where(wrf_time==datetime.date.toordinal(start_date))
wrf_start_index = wrf_start_index[0][0]
wrf_end_index = np.where(wrf_time==datetime.date.toordinal(end_date))
wrf_end_index=wrf_end_index[0][0]+1

# Print ordinal times to check against ROMS/SWAN output (should match)
print(wrf_time[wrf_start_index:wrf_end_index])

#Find ROMS time indices
roms_origin_date = datetime.datetime(2013,8,30,0,0,0)
roms_time=(np.array(roms_dataset['time'][:]))/24+datetime.date.toordinal(roms_origin_date)

roms_start_index = np.where(roms_time==datetime.date.toordinal(start_date))
roms_start_index = roms_start_index[0][0]
roms_end_index = np.where(roms_time==datetime.date.toordinal(end_date))
roms_end_index=roms_end_index[0][0]+1

# Print ordinal times to check against WRF output (should match)
print(roms_time[roms_start_index:roms_end_index])

#Make some plots

fig1=plt.figure(figsize=(30,20),dpi=300)
#fig1.plt.figsize(30,20)

for t in range(0,roms_end_index-roms_start_index):
    #t=10

    fig1.clf()
    
    #Ingest data
    slp=np.array(wrf_dataset['slp'][wrf_start_index+t,:,:])
    slp=np.squeeze(slp)
    sst=np.array(roms_dataset['temp'][roms_start_index+t,35,:,:])
    sst=np.ma.array(sst,mask=np.isnan(sst))
    sst=np.squeeze(sst)
    u10=np.array(wrf_dataset['u_10m_tr'][wrf_start_index+t,:,:])
    u10=np.squeeze(u10)
    v10=np.array(wrf_dataset['v_10m_tr'][wrf_start_index+t,:,:])
    v10=np.squeeze(v10)
    wnd_mag = (u10**2 + v10**2 ) **0.5
    cp=np.array(wrf_dataset['precip_c'][wrf_start_index+t,:,:])
    cp=np.squeeze(cp)
    gp=np.array(wrf_dataset['precip_g'][wrf_start_index+t,:,:])
    gp=np.squeeze(gp)
    precip = (cp + gp) * 0.0393701  #convert to inches
    precip=np.ma.array(precip,mask=(precip<0.01))
    mdbz=np.array(wrf_dataset['mdbz'][wrf_start_index+t,:,:])
    mdbz=np.squeeze(mdbz)
    mdbz=np.ma.array(mdbz,mask=(mdbz<10))
    wave=np.array(roms_dataset['Hwave'][roms_start_index+t,:,:])
    wave=np.ma.array(wave,mask=np.isnan(wave))
    wave=np.squeeze(wave)

    date_t=datetime.date.fromordinal(wrf_time[wrf_start_index+t].astype(int))
    date_hrs=(wrf_time[wrf_start_index+t]-wrf_time[wrf_start_index+t].astype(int))*24
    date_hrs=date_hrs.astype(int)
    fore_valid = datetime.datetime(date_t.year,date_t.month,date_t.day,date_hrs)

    plt.suptitle('Forecast Valid: ' + fore_valid.strftime("%d %b %Y %H"+"Z UTC"),fontsize=36,family='Helvetica')

    plt.subplot(2,3,1)
    map = Basemap(projection='merc',
        resolution='c',lat_0=((np.max(wrf_lat)-np.min(wrf_lat))/2),
        lon_0=((np.max(wrf_lon)-np.min(wrf_lon))/2),
        llcrnrlon=np.min(wrf_lon),llcrnrlat=np.min(wrf_lat),
        urcrnrlon=np.max(wrf_lon),urcrnrlat=np.max(wrf_lat))
    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()
    map.pcolormesh(wrf_lon,wrf_lat,slp[:,:],cmap='jet',vmin=950,vmax=1030,latlon='true')
    map.colorbar(location='bottom')
    plt.title(('Sea Level Pressure (hPa)'),fontsize=24,family='Helvetica')

    plt.subplot(2,3,2)
    map = Basemap(projection='merc',
        resolution='c',lat_0=((np.max(wrf_lat)-np.min(wrf_lat))/2),
        lon_0=((np.max(wrf_lon)-np.min(wrf_lon))/2),
        llcrnrlon=np.min(wrf_lon),llcrnrlat=np.min(wrf_lat),
        urcrnrlon=np.max(wrf_lon),urcrnrlat=np.max(wrf_lat))
    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()
    map.pcolormesh(roms_lon,roms_lat,sst[:,:],cmap='jet',vmin=5,vmax=32,latlon='true')
    map.colorbar(location='bottom')
    plt.title(('Sea Surface Temperature ($^\circ$C)'),fontsize=24,family='Helvetica')

    plt.subplot(2,3,3)
    map = Basemap(projection='merc',
        resolution='c',lat_0=((np.max(wrf_lat)-np.min(wrf_lat))/2),
        lon_0=((np.max(wrf_lon)-np.min(wrf_lon))/2),
        llcrnrlon=np.min(wrf_lon),llcrnrlat=np.min(wrf_lat),
        urcrnrlon=np.max(wrf_lon),urcrnrlat=np.max(wrf_lat))
    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()
    map.pcolormesh(wrf_lon,wrf_lat,wnd_mag[:,:],cmap='jet',vmin=0,vmax=30,latlon='true')
    map.quiver(wrf_lon[::10,::10],wrf_lat[::10,::10],u10[::10,::10],v10[::10,::10],latlon='true',               width=0.001,scale=1 / 0.001)
    map.colorbar(location='bottom')
    plt.title(('Wind Speed (m/s) + Direction'),fontsize=24,family='Helvetica')

    plt.subplot(2,3,4)
    map = Basemap(projection='merc',
        resolution='c',lat_0=((np.max(wrf_lat)-np.min(wrf_lat))/2),
        lon_0=((np.max(wrf_lon)-np.min(wrf_lon))/2),
        llcrnrlon=np.min(wrf_lon),llcrnrlat=np.min(wrf_lat),
        urcrnrlon=np.max(wrf_lon),urcrnrlat=np.max(wrf_lat))
    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()
    map.pcolormesh(wrf_lon,wrf_lat,precip[:,:],cmap='jet',vmin=0,vmax=1,latlon='true')
    map.colorbar(location='bottom')
    plt.title(('Precipitation (inches)'),fontsize=24,family='Helvetica')

    plt.subplot(2,3,5)
    map = Basemap(projection='merc',
        resolution='c',lat_0=((np.max(wrf_lat)-np.min(wrf_lat))/2),
        lon_0=((np.max(wrf_lon)-np.min(wrf_lon))/2),
        llcrnrlon=np.min(wrf_lon),llcrnrlat=np.min(wrf_lat),
        urcrnrlon=np.max(wrf_lon),urcrnrlat=np.max(wrf_lat))
    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()
    map.pcolormesh(wrf_lon,wrf_lat,mdbz[:,:],cmap='jet',vmin=0,vmax=50,latlon='true')
    map.colorbar(location='bottom')
    plt.title(('Sim. Radar Reflectivity (dBZ)'),fontsize=24,family='Helvetica')

    plt.subplot(2,3,6)
    map = Basemap(projection='merc',
        resolution='c',lat_0=((np.max(wrf_lat)-np.min(wrf_lat))/2),
        lon_0=((np.max(wrf_lon)-np.min(wrf_lon))/2),
        llcrnrlon=np.min(wrf_lon),llcrnrlat=np.min(wrf_lat),
        urcrnrlon=np.max(wrf_lon),urcrnrlat=np.max(wrf_lat))
    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()
    map.pcolormesh(roms_lon,roms_lat,wave[:,:],cmap='jet',vmin=0,vmax=10,latlon='true')
    map.colorbar(location='bottom')
    plt.title(('Sig. Wave Height (m)'),fontsize=24,family='Helvetica')

    # Save image
    plt.savefig('matthew_'+ str(t).zfill(2))

