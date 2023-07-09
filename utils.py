#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 11:26:45 2022

@author: mohammad.elaabaribaoune

Ce script a pour but de:

	1. Lire les données du fichier '.nc' notamment le géopotentiel.
	2. Lire les dates pluvieuses depuis un fichier 'txt'.
	3. Séléctinner les données corespondantes aux dates pluvieuses.
	4. Ploter la moyenne de ces donnée sur une carte.

Pour exucuter ce script : 

	 1. Il suppose avoir intallé : numpy, netCDF4, matplotlib, Basemap et datetime.
	 2. IL supose avoir les deux fichiers 'era_xx.nc' et 'Classic_dates' sur le même reprtoire que ce script.
	 
"""

import os

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from sklearn.preprocessing import MinMaxScaler, scale
import numpy as np
import  netCDF4 as nc4
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap as LSCmap
from datetime import datetime
import Nio
#import cartopy.crs as ccrs


def PlotMap(lonin,latin,data,**kwargs):

	mp = Basemap(llcrnrlat=np.min(latin),urcrnrlat=np.max(latin),llcrnrlon=np.min(lonin),urcrnrlon=np.max(lonin))#,resolution='c')
	plt.figure(figsize=(12,12),dpi=400)
	ax = plt.gca()
	# Set axis 
	#lon, lat = np.meshgrid(lon, lat)
	# compute native map projection coordinates of lat/lon grid.
	#x, y = mp(lon, lat)

	# Draw coastlines 
	mp.drawcoastlines(linewidth=0.25)
	#mp.drawcountries(linewidth=0.25)

	# draw lat/lon grid lines every 30 degrees.
#	mp.drawmeridians(np.arange(-180,180,30))
	mp.drawparallels(np.arange(-90,90,30))
	avail_cmap = [cm  for cm in plt.colormaps() if cm not in dir(plt.cm)]
		

	# Data to plot :
	#sc = mp.pcolor(x, y, data, cmap = 'jet', vmin= 270, vmax=300)

	colorbar_bounds = (np.min(data),np.max(data))
	if 'colorbar_bounds' in kwargs:
		colorbar_bounds = kwargs['colorbar_bounds']

	cmap = 'jet'
#	cmap = 'RdBu_r'
	if 'cmap' in kwargs:
		cmap = kwargs['cmap']

	if 'PlotUnderMinWiteColor' in kwargs:
		if kwargs['PlotUnderMinWiteColor'] == True:
			cmap.set_under('w')

	marker = 'o'
	if 'marker' in kwargs:
		marker = kwargs['marker']
	size = 30
	if 'size' in kwargs:
		size = kwargs['size']
	
	scatter = False
	if 'scatter' in kwargs:
		size = kwargs['scatter']
		

	
	parallels = np.arange(0.,81,20.)
	# labels = [left,right,top,bottom]
	mp.drawparallels(parallels,labels=[False,True,True,False],fontsize=6)
	meridians = np.arange(10.,351.,40.)
	mp.drawmeridians(meridians,labels=[True,False,False,True],fontsize=6)
	
	
	
	lon, lat = np.meshgrid(lonin, latin)
	xx, yy = mp(lon, lat)
	# Draw coastlines 
	mp.drawcoastlines(linewidth=1.0)
	#mp.drawcountries(linewidth=0.25)
	# Data to plot :
	#sc = mp.pcolor(xx, yy, data, cmap = 'jet'
		
	cmapMed = LSCmap.from_list('cmapMed', [(0 , 'white'),
										(0.4, 'lime'),
										(0.45, 'green'),
										(0.5, 'blue'),
										(0.52, 'orange'),
										(0.6, 'yellow'),
										(1, 'red')])
	if not 'cmapMed' in  avail_cmap : plt.cm.register_cmap(name='cmapMed', cmap=cmapMed)

	cmapMed_centred_on_white = LSCmap.from_list('cmapMed_centred_on_white', [(0 , 'blue'),
										(0.2, 'green'),
										(0.3, 'lime'),
										(0.48, 'white'),
										(0.5, 'white'),
										(0.52, 'white'),
										(0.7, 'yellow'),
										(0.8, 'orange'),
										(1, 'red')])
	if not 'cmapMed_centred_on_white' in  avail_cmap : plt.cm.register_cmap(name='cmapMed_centred_on_white', cmap=cmapMed_centred_on_white)

	cmapMed_centred_on_white_for_diff = LSCmap.from_list('cmapMed_centred_on_white_for_diff', 
										[(0 , 'blue'),
										(0.2, 'green'),
										(0.3, 'lime'),
										(0.48, 'white'),
										(0.5, 'white'),
										(0.52, 'white'),
										(0.7, 'red'),
										(0.8, 'brown'),
										(1, 'black')])
		  
	if not 'cmapMed_centred_on_white_for_diff' in  avail_cmap : plt.cm.register_cmap(name='cmapMed_centred_on_white_for_diff', cmap=cmapMed_centred_on_white_for_diff)

	cmap=plt.cm.get_cmap('jet') #seismic rainbow nipy_spectral
	
	if 'cmap' in kwargs.keys():
		cmap = kwargs['cmap']
	
	if 'set_under_for_RH' in kwargs.keys():
		set_under = kwargs['set_under_for_RH']
		if set_under:
			cmap = plt.cm.get_cmap(cmap)
			cmap.set_under('white', 0.08)
	levels = 15
	if 'levels' in kwargs.keys():
		levels = kwargs['levels']
	if 'levels_for_diff' in kwargs.keys(): 
		levels = kwargs['levels_for_diff']
	
	sc = mp.contourf(xx, yy, data, cmap =cmap,levels= levels,extend="both")
	#sc = mp.contourf(xx, yy, data, cmap =cmap,levels=15)
#	sc = mp.contourf(xx, yy, data, cmap =cmap,levels=18)
	plt.colorbar(sc,orientation='horizontal',pad=0.04)
	
	sc = mp.contour(xx, yy, data,  colors='k', linewidths=0.3)
	ax.clabel(sc,  inline=True, fmt='%1.0f', fontsize=8, colors='k')#,levels=sc.levels,)
		
	# Add station 
	lat_station = [27.16]
	lon_station = [-13.21]
	lon_station, lat_station = mp(lon_station, lat_station)
	#plt.scatter(lon_station, lat_station, linewidths = 2, marker="^", color='green', edgecolor ="blue") 
	
	# add wind on map 
	if 'plotWind' in kwargs.keys():
		plotWind = kwargs['plotWind']
		if plotWind:
			# plot wind vectors on projection grid.
			# first, shift grid so it goes from -180 to 180 (instead of 0 to 360
			# in longitude).  Otherwise, interpolation is messed up.
			#ugrid,newlons = shiftgrid(180.,kwargs['u'],lon,start=False)
			#vgrid,newlons = shiftgrid(180.,kwargs['v'],lon,start=False)
			# transform vectors to projection grid.
			#uproj,vproj,xx,yy = mp.transform_vector(ugrid,vgrid,newlons,lat,31,31,returnxy=True,masked=True)
				# now plot.
			Q = mp.quiver(xx[::6, ::6],yy[::6, ::6],kwargs['u'][::6, ::6],kwargs['v'][::6, ::6])#,scale=700)#uproj,vproj,scale=700)
			# make quiver key.
			qk = plt.quiverkey(Q, 0.1, 0.1, 20, '20 m/s', labelpos='W')
	

	title = ''
	if 'title' in kwargs.keys():
		title = kwargs['title']

	namefig = 'RR_6stat_NAO+_Classic_neut'
	if 'namefig' in kwargs.keys():
		namefig = kwargs['namefig']
		if namefig[-3:] != 'png':
			namefig=namefig+'.png'


	outdir = os.path.join(os.getcwd(),'Figs')
	if 'outdir' in kwargs.keys():
		outdir = kwargs['outdir']
	if not os.path.exists(outdir) : os.mkdir(outdir) 

	to_save = os.path.join(outdir,namefig)
	print ('  Saving ',to_save)
	plt.title(title, fontsize=11,fontweight='bold')
	plt.savefig(to_save)
	#plt.show()


def plotTwoflds(lonin,latin,data1,data2,**kwargs):

	mp = Basemap(llcrnrlat=np.min(latin),urcrnrlat=np.max(latin),llcrnrlon=np.min(lonin),urcrnrlon=np.max(lonin))#,resolution='c')
	plt.figure(figsize=(12,12),dpi=400)
	ax = plt.gca()
	# Set axis 
	#lon, lat = np.meshgrid(lon, lat)
	# compute native map projection coordinates of lat/lon grid.
	#x, y = mp(lon, lat)

	# Draw coastlines 
	mp.drawcoastlines(linewidth=0.25)
	#mp.drawcountries(linewidth=0.25)

	# draw lat/lon grid lines every 30 degrees.
#	mp.drawmeridians(np.arange(-180,180,30))
	mp.drawparallels(np.arange(-90,90,30))
	avail_cmap = [cm  for cm in plt.colormaps() if cm not in dir(plt.cm)]
		

	# Data to plot :
	#sc = mp.pcolor(x, y, data, cmap = 'jet', vmin= 270, vmax=300)

	colorbar_bounds = (np.min(data1),np.max(data1))
	if 'colorbar_bounds' in kwargs:
		colorbar_bounds = kwargs['colorbar_bounds']

	cmap = 'jet'
#	cmap = 'RdBu_r'
	if 'cmap' in kwargs:
		cmap = kwargs['cmap']

	if 'PlotUnderMinWiteColor' in kwargs:
		if kwargs['PlotUnderMinWiteColor'] == True:
			cmap.set_under('w')

	marker = 'o'
	if 'marker' in kwargs:
		marker = kwargs['marker']
	size = 30
	if 'size' in kwargs:
		size = kwargs['size']
	
	scatter = False
	if 'scatter' in kwargs:
		size = kwargs['scatter']
		

	
	parallels = np.arange(0.,81,20.)
	# labels = [left,right,top,bottom]
	mp.drawparallels(parallels,labels=[False,True,True,False],fontsize=6)
	meridians = np.arange(10.,351.,40.)
	mp.drawmeridians(meridians,labels=[True,False,False,True],fontsize=6)
	
	
	
	lon, lat = np.meshgrid(lonin, latin)
	xx, yy = mp(lon, lat)
	# Draw coastlines 
	mp.drawcoastlines(linewidth=1.0)
	#mp.drawcountries(linewidth=0.25)
	# Data to plot :
	#sc = mp.pcolor(xx, yy, data, cmap = 'jet'
		
	cmapMed = LSCmap.from_list('cmapMed', [(0 , 'white'),
										(0.4, 'lime'),
										(0.45, 'green'),
										(0.5, 'blue'),
										(0.52, 'orange'),
										(0.6, 'yellow'),
										(1, 'red')])
	if not 'cmapMed' in  avail_cmap : plt.cm.register_cmap(name='cmapMed', cmap=cmapMed)

	cmapMed_centred_on_white = LSCmap.from_list('cmapMed_centred_on_white', [(0 , 'blue'),
										(0.2, 'green'),
										(0.3, 'lime'),
										(0.48, 'white'),
										(0.5, 'white'),
										(0.52, 'white'),
										(0.7, 'yellow'),
										(0.8, 'orange'),
										(1, 'red')])
	if not 'cmapMed_centred_on_white' in  avail_cmap : plt.cm.register_cmap(name='cmapMed_centred_on_white', cmap=cmapMed_centred_on_white)

	cmapMed_centred_on_white_for_diff = LSCmap.from_list('cmapMed_centred_on_white_for_diff', 
										[(0 , 'blue'),
										(0.2, 'green'),
										(0.3, 'lime'),
										(0.48, 'white'),
										(0.5, 'white'),
										(0.52, 'white'),
										(0.7, 'red'),
										(0.8, 'brown'),
										(1, 'black')])
		  
	if not 'cmapMed_centred_on_white_for_diff' in  avail_cmap : plt.cm.register_cmap(name='cmapMed_centred_on_white_for_diff', cmap=cmapMed_centred_on_white_for_diff)

	cmap=plt.cm.get_cmap('jet') #seismic rainbow nipy_spectral
	
	if 'cmap' in kwargs.keys():
		cmap = kwargs['cmap']
	
	if 'set_under_for_RH' in kwargs.keys():
		set_under = kwargs['set_under_for_RH']
		if set_under:
			cmap = plt.cm.get_cmap(cmap)
			cmap.set_under('white', 0.08)
	levels = 15
	if 'levels' in kwargs.keys():
		levels = kwargs['levels']
	if 'levels_for_diff' in kwargs.keys(): 
		levels = kwargs['levels_for_diff']
	
	sc = mp.contourf(xx, yy, data1, cmap =cmap,levels= levels,extend="both")
	#sc = mp.contourf(xx, yy, data, cmap =cmap,levels=15)
#	sc = mp.contourf(xx, yy, data, cmap =cmap,levels=18)
	plt.colorbar(sc,orientation='horizontal',pad=0.04)
	
	sc = mp.contour(xx, yy, data2,  colors='k', linewidths=0.3)
	ax.clabel(sc,  inline=True, fmt='%1.0f', fontsize=8, colors='k')#,levels=sc.levels,)
		
	# Add station 
	lat_station = [27.16]
	lon_station = [-13.21]
	lon_station, lat_station = mp(lon_station, lat_station)
	#plt.scatter(lon_station, lat_station, linewidths = 2, marker="^", color='green', edgecolor ="blue") 
	
	# add wind on map 
	if 'plotWind' in kwargs.keys():
		plotWind = kwargs['plotWind']
		if plotWind:
			# plot wind vectors on projection grid.
			# first, shift grid so it goes from -180 to 180 (instead of 0 to 360
			# in longitude).  Otherwise, interpolation is messed up.
			#ugrid,newlons = shiftgrid(180.,kwargs['u'],lon,start=False)
			#vgrid,newlons = shiftgrid(180.,kwargs['v'],lon,start=False)
			# transform vectors to projection grid.
			#uproj,vproj,xx,yy = mp.transform_vector(ugrid,vgrid,newlons,lat,31,31,returnxy=True,masked=True)
				# now plot.
			Q = mp.quiver(xx[::6, ::6],yy[::6, ::6],kwargs['u'][::6, ::6],kwargs['v'][::6, ::6])#,scale=700)#uproj,vproj,scale=700)
			# make quiver key.
			qk = plt.quiverkey(Q, 0.1, 0.1, 20, '20 m/s', labelpos='W')
	

	title = ''
	if 'title' in kwargs.keys():
		title = kwargs['title']

	namefig = 'RR_6stat_NAO+_Classic_neut'
	if 'namefig' in kwargs.keys():
		namefig = kwargs['namefig']
		if namefig[-3:] != 'png':
			namefig=namefig+'.png'


	outdir = os.path.join(os.getcwd(),'Figs')
	if 'outdir' in kwargs.keys():
		outdir = kwargs['outdir']
	if not os.path.exists(outdir) : os.mkdir(outdir) 

	to_save = os.path.join(outdir,namefig)
	print ('  Saving ',to_save)
	plt.title(title, fontsize=11,fontweight='bold')
	plt.savefig(to_save)
	#plt.show()

	
def Plotmap_subplots(lon,lat,ls_data,nber_of_plos_horz,nber_of_plos_vrtcl,**kwargs):
	# Create a map instance using basemap classe. 
	# Remarque: Basemap class constructor have many arguments and all are optional,
	# Dafault values are taken if none of the argumnets was not choosen.
	#map = Basemap(projection='ortho', llcrnrx=-3000000, llcrnry=1000000, urcrnrx=3000000, urcrnry=6000000, \
		#lat_0=10, lon_0=0, resolution='c') 
	fig, axs = plt.subplots(nber_of_plos_horz, nber_of_plos_vrtcl,dpi=400,figsize=(15, 5))
	fig.subplots_adjust(hspace = 1.5, wspace=.1)
	axs = axs.ravel()
	#fig, axes = plt.subplots(figsize=(12,12),dpi=300,nrows=nber_of_plos_horz, ncols=nber_of_plos_vrtcl)
	
	for ax_i, ax in enumerate(axs.flat):
		
		print('Plotting field number : ', ax_i)
		data = ls_data[ax_i]
		#mp = Basemap(ax=ax)#projection='merc')#,llcrnrlat=-90,urcrnrlat=90,\
		mp = Basemap(llcrnrlat=np.min(lat),urcrnrlat=np.max(lat),llcrnrlon=np.min(lon),urcrnrlon=np.max(lon))
		#llcrnrlon=-180,urcrnrlon=180,resolution='c')
		
		# Set axis 
		#lon, lat = np.meshgrid(lon, lat)
		# compute native map projection coordinates of lat/lon grid.
		#x, y = mp(lon, lat)

		# Draw coastlines 
		mp.drawcoastlines(linewidth=0.25)
		#mp.drawcountries(linewidth=0.25)

		# draw lat/lon grid lines every 30 degrees.
		#mp.drawmeridians(np.arange(-180,180,30))
		#mp.drawparallels(np.arange(-90,90,30))
	
		# Data to plot :
		#sc = mp.pcolor(x, y, data, cmap = 'jet', vmin= 270, vmax=300)
	
		colorbar_bounds = (np.min(data),np.max(data))
		if 'colorbar_bounds' in kwargs:
			colorbar_bounds = kwargs['colorbar_bounds']
	
		cmap = 'RdBu_r'
		if 'cmap' in kwargs:
			cmap = kwargs['cmap']
	
		if 'PlotUnderMinWiteColor' in kwargs:
			if kwargs['PlotUnderMinWiteColor'] == True:
				cmap.set_under('w')
	
		marker = 'o'
		if 'marker' in kwargs:
			marker = kwargs['marker']
		size = 30
		if 'size' in kwargs:
			size = kwargs['size']
		
		scatter = False
		if 'scatter' in kwargs:
			size = kwargs['scatter']
			
	
		
		parallels = np.arange(0.,81,20.) #[ int(x) for x in  lat[0:-1:30]]# np.arange(0.,81,20.)
		# labels = [left,right,top,bottom]
		mp.drawparallels(parallels,labels=[False,True,True,False],fontsize=6)
		meridians = np.arange(10.,351.,40.)# [ int(x) for x in lon[0:-1:60]] #np.arange(10.,351.,40.)
		mp.drawmeridians(meridians,labels=[True,False,False,True],fontsize=6)
		
		
		
		lons, lats = np.meshgrid(lon, lat)
		xx, yy = mp(lons, lats)
		# Draw coastlines 
		mp.drawcoastlines(linewidth=0.1)
		#mp.drawcountries(linewidth=0.25)
		avail_cmap = [cm  for cm in plt.colormaps() if cm not in dir(plt.cm)]
		#sc = mp.pcolor(xx, yy, data, cmap = 'jet')
		
		
		
		
		cmapMed = LSCmap.from_list('cmapMed', [(0 , 'white'),
														(0.4, 'lime'),
														(0.45, 'green'),
														(0.5, 'blue'),
														(0.52, 'orange'),
														(0.6, 'yellow'),
														(1, 'red')])
		if not 'cmapMed' in  avail_cmap : plt.cm.register_cmap(name='cmapMed', cmap=cmapMed)
		
		cmapMed_centred_on_white = LSCmap.from_list('cmapMed_centred_on_white', [(0 , 'blue'),
														(0.2, 'green'),
														(0.3, 'lime'),
														(0.48, 'white'),
														(0.5, 'white'),
														(0.52, 'white'),
														(0.7, 'yellow'),
														(0.8, 'orange'),
														(1, 'red')])
		if not 'cmapMed_centred_on_white' in  avail_cmap : plt.cm.register_cmap(name='cmapMed_centred_on_white', cmap=cmapMed_centred_on_white)
	  
		cmapMed_centred_on_white_for_diff = LSCmap.from_list('cmapMed_centred_on_white_for_diff', 
														[(0 , 'blue'),
														(0.2, 'green'),
														(0.3, 'lime'),
														(0.48, 'white'),
														(0.5, 'white'),
														(0.52, 'white'),
														(0.7, 'red'),
														(0.8, 'brown'),
														(1, 'black')])
						  
		if not 'cmapMed_centred_on_white_for_diff' in  avail_cmap : plt.cm.register_cmap(name='cmapMed_centred_on_white_for_diff', cmap=cmapMed_centred_on_white_for_diff)
		cmap = plt.cm.get_cmap('cmapMed_centred_on_white', 32)
		
		#levels=np.arange(-0.1, 0.1, .005) #16; precip(0.03, 0.5, .05) 
		#levels=[-5,-4,-2,-1,-0.5,0.5,1,2,3, 4, 5]
		#levels=np.arange(np.min(data), np.max(data),1)
		#cmap='hsv'#plt.cm.get_cmap('jet', 12) #seismic
		if 'levels' in kwargs.keys():
			levels = kwargs['levels']
		if 'levels_for_diff' in kwargs.keys(): 
			if ('title' in kwargs) and  ('Diff' in kwargs['title'][ax_i]) :
				levels = kwargs['levels_for_diff']
		if 'ls_cmaps' in kwargs.keys():
			cmap = kwargs['ls_cmaps'][ax_i]
			if cmap=='cmapMed':
				cmap = plt.cm.get_cmap('cmapMed')
			if cmap=='cmapMed_centred_on_white':
				cmap = plt.cm.get_cmap('cmapMed_centred_on_white')
			if cmap=='cmapMed_centred_on_white_for_diff':
				cmap = plt.cm.get_cmap('cmapMed_centred_on_white_for_diff')
		  
		if 'set_under_for_RH' in kwargs.keys():
			set_under = kwargs['set_under_for_RH']
			if set_under:
				cmap.set_under('white', 10)
		sc = mp.contourf(xx, yy, data, cmap =cmap,levels= levels,extend="both")
		#sc.cmap.set_under('k')
		#sc.set_clim(190, 315)
		
		fig.colorbar(sc,orientation='horizontal',pad=0.04,ax=ax)
		
		#if 'precip' in kwargs['namefig']: Colorabar nonlinear 
			#fig.colorbar(sc,ticks=[0.5,1,2,3,4,6,8,10,12,15,20], orientation='horizontal',pad=0.04,ax=ax)
		
		sc.levels = np.arange(start=np.min(sc.levels), stop=np.max(sc.levels)-2, step=4)
		sc = mp.contour(xx, yy, data, levels=sc.levels, linewidths=0.3, colors='k')
		
		#sc.levels = np.arange(start=np.min(sc.levels), stop=np.max(sc.levels)-2, step=10)
	
		ax.clabel(sc,  inline=True, fmt='%.1f', fontsize=4, colors='k')#,levels=sc.levels,)
		
		
		#plt.colorbar(sc,cax=ax)
		
		#sc = mp.contourf(xx, yy, data, cmap = 'jet')#,levels=8)
	
		# Colorbar 
	
		
		#norm = mpl.colors.Normalize(vmin=0,vmax=2)
		#sm = plt.cm.ScalarMappable(cmap="jet")
		#sm.set_array([])
	
		#plt.colorbar(sm)
		#mp.colorbar(cmap=cm.jet)
	
		title = ''
		if 'title' in kwargs.keys():
			title = kwargs['title'][ax_i]
		ax.set_title(title)
	
	namefig = 'Fig.png'
	if 'namefig' in kwargs.keys():
		namefig = kwargs['namefig']
		if namefig[-3:] != 'png':
			namefig=namefig+'.png'
	
	outdir = os.getcwd()
	if 'outdir' in kwargs.keys():
		outdir = kwargs['outdir']
			
	to_save = os.path.join(outdir,namefig)
	print ('  Saving ',to_save)
	plt.savefig(to_save)
	#plt.show()

def read_var_from_nc(var, file_name):
	ncfile = Nio.open_file(file_name, format="netcdf")
	data = ncfile.variables[var][:].squeeze()
	lons = ncfile.variables['longitude'][:].squeeze()# Extraire les longitudes 
	lats = ncfile.variables['latitude'][:].squeeze() # Extraire les latitudes 
	temps = ncfile.variables['time'][:].squeeze() # Extraire le temps 
	dtime = nc4.num2date(ncfile.variables['time'],ncfile.variables['time'].units)
	#tt = datetime.strptime(str(dtime[0]),'%Y-%m-%d %H:%M:%S')
	#print(tt.year, tt.month, tt.day, tt.hour, tt.second)
	temps_splitted = []
	[temps_splitted.append(datetime.strptime(str(dtime[elem]),'%Y-%m-%d %H:%M:%S')) for elem in range(len(temps))]
	ncfile.close()
	return(np.squeeze(lons), np.squeeze(lats), np.squeeze(temps_splitted), np.squeeze(data))
	
	
def get_var_from_nc(var, filename, read_by_chunks=False):

	ncfile = nc4.Dataset(filename)
	#print(ncfile.variables.keys())
	lons = ncfile['longitude'][:].data # Extraire les longitudes 
	lats = ncfile['latitude'][:].data # Extraire les latitudes 
	temps = ncfile['time'][:].data # Extraire le temps 
	dtime = nc4.num2date(ncfile.variables['time'],ncfile.variables['time'].units)
	#tt = datetime.strptime(str(dtime[0]),'%Y-%m-%d %H:%M:%S')
	#print(tt.year, tt.month, tt.day, tt.hour, tt.second)
	temps_splitted = []
	[temps_splitted.append(datetime.strptime(str(dtime[elem]),'%Y-%m-%d %H:%M:%S')) for elem in range(len(temps))] 
	
	if read_by_chunks:
		# To chunks: proceed 2000 by 2000
		var_all_chunks  = []
		dim_time = len(temps)
		length_check = 4000
		number_chunks  = dim_time//length_check
		#rest = dim_time%length_check
		print('   +> extracting data by chnks for', number_chunks, ' of length ', length_check)
		for k in range(number_chunks):
			bnd1 = k*length_check
			bnd2 = (k+1)*length_check
			print('   +> extracting data between : ', bnd1, ' and ', bnd2) 

	#        name_var = (ncfile["t"][bnd1:bnd2, :, :].data)-273.15
			var_chunk = ncfile[var][bnd1:bnd2, :, :].data/100
			var_all_chunks.extend(var_chunk)
		data_var = var_all_chunks
			
	else:   
		print (' ')
		data_var = ncfile[var][:]
	
	return(lons, lats, temps_splitted, data_var)

def select_precip_days(field, dates, temps):
	
	data_of_precipitant_days = []
	for date in dates:
		for it in range (len(field)):
			if temps[it].year == int(date[0]) : 
				if temps[it].month == int(date[1]) :
					if temps[it].day == int(date[2]) :
						#print('date : ', date, 'found !')
						data_of_precipitant_days.append(np.squeeze(field[it]))
						break
	return(data_of_precipitant_days)

def read_dates(path_to_dates_file) : 
	classic_dates=[]
	with open(path_to_dates_file) as f:  
		#[print(line.split("	")[-1]) for line in f.readlines()]
		[classic_dates.append(line[:-1].split("	")) for line in f.readlines()]
	return(classic_dates)


def plot_wind(lon,lat,u,v):
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
	ax.quiver(lon, lat, u, v, transform=ccrs.PlateCarree())