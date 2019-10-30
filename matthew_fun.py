# Program to plot DINEOF and Original MODIS SST and Chlor-a for Lindsay Dubbs / NCROEP
#    Usage: python -W ignore plot.py YYYY MM DD LON LAT
#    e.g. : python -W ignore plot.py 2017 04 13 -75.70838 34.91828
#
# Joseph B. Zambon
# jbzambon@ncsu.edu
# 29 October 2019
#
# Using conda, create environment, activate, and run code
# conda env create -f ncroep.yml 
# source activate ncroep
# python -W ignore plot.py YYYY MM DD LON LAT

#Dependencies
from pydap.client import open_url
import numpy as np
import datetime
import cmocean
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
import matplotlib.colorbar as cb
import matplotlib.colors as colors

# Define Date Range
date_start = str(sys.argv[1]+'-'+sys.argv[2]+'-'+sys.argv[3])
#date_start = '2017-04-13'
#date_end   = '2017-04-30'

# Define location
#plot_loc = [-75.70838, 35.91828]
plot_loc = [np.float(sys.argv[4]),np.float(sys.argv[5])]

# Define SST colorbar range ºC
sst_range = [16,30]

# OPeNDAP linked datasets
#modis_sst_url = 'http://oceanus.meas.ncsu.edu:8080/thredds/dodsC/secoora/modis/sst.nc'
#modis_chl_url = 'http://oceanus.meas.ncsu.edu:8080/thredds/dodsC/secoora/modis/chla.nc'
#dineof_sst_url = 'http://oceanus.meas.ncsu.edu:8080/thredds/dodsC/secoora/dineof/sst.nc'
#dineof_chl_url = 'http://oceanus.meas.ncsu.edu:8080/thredds/dodsC/secoora/dineof/chla.nc'
modis_sst_url = 'http://oceanus.meas.ncsu.edu:8080/thredds/dodsC/carolinas/coastwatch/sst.nc'
modis_chl_url = 'http://oceanus.meas.ncsu.edu:8080/thredds/dodsC/carolinas/coastwatch/chla.nc'
dineof_sst_url = 'http://oceanus.meas.ncsu.edu:8080/thredds/dodsC/carolinas/coastwatch/dineof_sst.nc'
dineof_chl_url = 'http://oceanus.meas.ncsu.edu:8080/thredds/dodsC/carolinas/coastwatch/dineof_chla.nc'
modis_sst_dataset = open_url(modis_sst_url)
modis_chl_dataset = open_url(modis_chl_url)
dineof_sst_dataset = open_url(dineof_sst_url)
dineof_chl_dataset = open_url(dineof_chl_url)

date_start = datetime.datetime.strptime(date_start,"%Y-%m-%d")
#date_end   = datetime.datetime.strptime(date_end,"%Y-%m-%d")
#num_days = date_end - date_start

# For inline plotting in Jupyter Notebook
#get_ipython().run_line_magic('pylab', 'inline')

parallels = np.arange(0.,90,2.)
meridians = np.arange(180.,360.,2.)

# For inline plotting in Jupyter Notebook
#figsize(22,20)
# For script-based plotting
fig=plt.figure(frameon=False,figsize=(22,20))

#for t in range(0,num_days.days+1):
for t in range(0,1):
    #Assume date ranges are congruent among datasets, just use modis_sst
    time=np.array(modis_sst_dataset['time'])
    t_ind = np.where(time==datetime.datetime.strftime((date_start + datetime.timedelta(t)),"%Y-%m-%d"+"T00:00:00Z"))
    t_ind = t_ind[0]
    if t_ind.size == 0:
        print('No matching time in dataset found.')
        exit(1)
    curr_date = date_start + datetime.timedelta(t)
    # Assume coordinates are congruent among datasets, just use modis_sst
    lon = np.array(modis_sst_dataset['lon'][:])
    lat = np.array(modis_sst_dataset['lat'][:])
    # Locate grid point closest to prescribed lat/lon
    diff_lon = lon - plot_loc[0]
    diff_lat = lat - plot_loc[1]
    min_lon = np.argmin(abs(diff_lon))
    min_lat = np.argmin(abs(diff_lat))
    # OPeNDAP linked datasets
    raw_sst = modis_sst_dataset['sst']
    raw_sst = raw_sst['sst'][int(t_ind),:,:]
    raw_sst = np.squeeze(raw_sst)
    raw_sst = np.ma.filled(raw_sst.astype(float), np.nan)
    raw_sst[raw_sst<-5]=np.nan; raw_sst= np.ma.array(raw_sst,mask=np.isnan(raw_sst))
    if np.isnan(raw_sst[min_lat,min_lon]) == 0:
        raw_sst_loc = str(round(raw_sst[min_lat,min_lon],1))
    else:
        raw_sst_loc = 'NaN'
    # Raw Chlor-a
    # OPeNDAP linked datasets
    raw_chla = modis_chl_dataset['chlor_a']
    raw_chla = raw_chla['chlor_a'][int(t_ind),:,:]
    raw_chla = np.squeeze(raw_chla)
    raw_chla = np.ma.filled(raw_chla.astype(float), np.nan)
    raw_chla[raw_chla<0]=np.nan; raw_chla= np.ma.array(raw_chla,mask=np.isnan(raw_chla))
    if np.isnan(raw_chla[min_lat,min_lon]) == 0:
        raw_chla_loc = str(round(raw_chla[min_lat,min_lon],3))
    else:
        raw_chla_loc = 'NaN'
    # DINEOF SST
    # OPeNDAP linked datasets
    dineof_sst = dineof_sst_dataset['sst']
    dineof_sst = dineof_sst['sst'][int(t_ind),:,:]
    dineof_sst = np.squeeze(dineof_sst)
    dineof_sst = np.ma.filled(dineof_sst.astype(float), np.nan)
    dineof_sst[dineof_sst<-5]=np.nan; dineof_sst= np.ma.array(dineof_sst,mask=np.isnan(dineof_sst))
    if np.isnan(dineof_sst[min_lat,min_lon]) == 0:
        dineof_sst_loc = str(round(dineof_sst[min_lat,min_lon],1))
    else:
        dineof_sst_loc = 'NaN '
    # Raw SST
    # OPeNDAP linked datasets
    dineof_chla = dineof_chl_dataset['chlor_a']
    dineof_chla = dineof_chla['chlor_a'][int(t_ind),:,:]
    dineof_chla = np.squeeze(dineof_chla)
    dineof_chla = np.ma.filled(dineof_chla.astype(float), np.nan)
    dineof_chla[dineof_chla<0]=np.nan; dineof_chla= np.ma.array(dineof_chla,mask=np.isnan(dineof_chla))
    if np.isnan(dineof_chla[min_lat,min_lon]) == 0:
        dineof_chla_loc = str(round(dineof_chla[min_lat,min_lon],3))
    else:
        dineof_chla_loc = 'NaN'
    plt.clf()
    plt.suptitle('1km Observed and Cloud Free: ' + curr_date.strftime("%d %b %Y %H"+"UTC"),fontsize=36,family='Helvetica')
    # Raw SST
    plt.subplot(2,2,1)
    map = Basemap(projection='merc',
      resolution='i',lat_0=((np.max(lat)-np.min(lat))/2),
      lon_0=((np.max(lon)-np.min(lon))/2),
      llcrnrlon=np.min(lon),llcrnrlat=np.min(lat),
      urcrnrlon=np.max(lon),urcrnrlat=np.max(lat))
    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()
    [X,Y] = np.meshgrid(lon,lat)
    map.pcolormesh(X,Y,raw_sst[:,:],cmap=cmocean.cm.thermal,                   vmin=sst_range[0],vmax=sst_range[1],latlon='true')
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=18)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=18)
    x_loc,y_loc = map(plot_loc[0], plot_loc[1])
    map.plot(x_loc, y_loc, 'ro', fillstyle='none',markersize=24)   #location
    tx_loc,ty_loc = map(-79.9, 37.5)
    plt.text(tx_loc,ty_loc,''.join([raw_sst_loc,'$^\circ$C']),fontsize=36,fontweight='bold',color='r')   #location
    plt.title(('Original SST ($^\circ$C)'),fontsize=24,family='Helvetica')
    cbar=map.colorbar(location='right',ticks=np.arange(sst_range[0],sst_range[1]+0.01,2))
    cbar.ax.tick_params(labelsize=20)
    # DINEOF SST
    plt.subplot(2,2,2)
    map = Basemap(projection='merc',
      resolution='i',lat_0=((np.max(lat)-np.min(lat))/2),
      lon_0=((np.max(lon)-np.min(lon))/2),
      llcrnrlon=np.min(lon),llcrnrlat=np.min(lat),
      urcrnrlon=np.max(lon),urcrnrlat=np.max(lat))
    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()
    map.pcolormesh(X,Y,dineof_sst[:,:],cmap=cmocean.cm.thermal,                   vmin=sst_range[0],vmax=sst_range[1],latlon='true')
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=18)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=18)
    x_loc,y_loc = map(plot_loc[0], plot_loc[1])
    map.plot(x_loc, y_loc, 'ro', fillstyle='none',markersize=24)   #location
    tx_loc,ty_loc = map(-79.9, 37.5)
    plt.text(tx_loc,ty_loc,''.join([dineof_sst_loc,'$^\circ$C']),fontsize=36,fontweight='bold',color='r')   #location
    plt.title(('Cloud Free SST ($^\circ$C)'),fontsize=24,family='Helvetica')
    cbar=map.colorbar(location='right',ticks=np.arange(sst_range[0],sst_range[1]+0.01,2))
    cbar.ax.tick_params(labelsize=20)
    # Raw Chlorophyll-a
    plt.subplot(2,2,3)
    map = Basemap(projection='merc',
      resolution='i',lat_0=((np.max(lat)-np.min(lat))/2),
      lon_0=((np.max(lon)-np.min(lon))/2),
      llcrnrlon=np.min(lon),llcrnrlat=np.min(lat),
      urcrnrlon=np.max(lon),urcrnrlat=np.max(lat))
    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()
    map.pcolormesh(X,Y,raw_chla[:,:],norm=LogNorm(vmin=0.01, vmax=100),                    cmap=cmocean.cm.algae,latlon='true')
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=18)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=18)
    x_loc,y_loc = map(plot_loc[0], plot_loc[1])
    map.plot(x_loc, y_loc, 'ro', fillstyle='none',markersize=24)   #location
    tx_loc,ty_loc = map(-79.9, 37.5)
    plt.text(tx_loc,ty_loc,''.join([raw_chla_loc,'mg/m$^3$']),fontsize=36,fontweight='bold',color='r')   #location
    plt.title(('Original Chl-a (mg/m$^3$)'),fontsize=24,family='Helvetica')
    cbar=map.colorbar(location='right',norm=LogNorm(vmin=0.01, vmax=100),                      ticks=[0.01,0.1,1,10,100])
    cbar.ax.tick_params(labelsize=20)
    # DINEOF Chlorophyll-a
    plt.subplot(2,2,4)
    #ax = fig.add_subplot(2,2,4)
    map = Basemap(projection='merc',
      resolution='i',lat_0=((np.max(lat)-np.min(lat))/2),
      lon_0=((np.max(lon)-np.min(lon))/2),
      llcrnrlon=np.min(lon),llcrnrlat=np.min(lat),
      urcrnrlon=np.max(lon),urcrnrlat=np.max(lat))
    map.drawcoastlines()
    map.drawcountries()
    map.drawstates()
    map.pcolormesh(X,Y,dineof_chla[:,:],norm=LogNorm(vmin=0.01, vmax=100),                    cmap=cmocean.cm.algae,latlon='true')
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=18)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=18)
    x_loc,y_loc = map(plot_loc[0], plot_loc[1])
    map.plot(x_loc, y_loc, 'ro', fillstyle='none',markersize=24)   #location
    tx_loc,ty_loc = map(-79.9, 37.5)
    plt.text(tx_loc,ty_loc,''.join([dineof_chla_loc,'mg/m$^3$']),fontsize=36,fontweight='bold',color='r')   #location    plt.title(('Cloud Free Chl-a (mg/m$^3$)'),fontsize=24,family='Helvetica')
    cbar=map.colorbar(location='right',norm=LogNorm(vmin=0.01, vmax=100),                      ticks=[0.01,0.1,1,10,100])
    cbar.ax.tick_params(labelsize=20)

    plt.savefig('ncroep_' + curr_date.strftime("%Y%m%d") + '.png')


