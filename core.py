#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 14:57:44 2023

@author: M. El Aabbaribaoune (@_um6p)
"""

import os
import sys
import numpy as np
import utils as ut 


def plot_flds_when_it_rains(plotcfg, inputdir,**kwargs):
	
	varibales = plotcfg['plotWhenItRains']['HorizMeans']['varibales']
	vars_names_in_netcdf = plotcfg['plotWhenItRains']['HorizMeans']['vars_names_in_netcdf']
	file_name_of_precip_days = plotcfg['plotWhenItRains']['file_dates_precip']
	file_name_of_extreme_days = plotcfg['plotWhenItRains']['file_name_of_extreme_days']
	
	bounds = plotcfg['plotWhenItRains']['HorizMeans']['cmap_bounds']
	bounds_diff = plotcfg['plotWhenItRains']['HorizMeans']['cmap_bounds_diff']
	cmaps = plotcfg['plotWhenItRains']['HorizMeans']['cmaps']
	cmaps_diff = plotcfg['plotWhenItRains']['HorizMeans']['cmaps_diff']
	
	plot_extreme = plotcfg['plotWhenItRains']['HorizMeans']['plot_extreme']

	plotWind = plotcfg['plotWhenItRains']['HorizMeans']['plotWindonVars']
	plotWindWithvars = plotcfg['plotWhenItRains']['HorizMeans']['plotWindWithvars']
	levWind = plotcfg['plotWhenItRains']['HorizMeans']['levWind']
	
	plottwovaronsameplot = plotcfg['plotWhenItRains']['HorizMeans']['plottwovaronsameplot']
	two_vars_on_same_plot = plotcfg['plotWhenItRains']['HorizMeans']['two_vars_on_same_plot']
	
	
	set_under_for_RH = False # To not modify
	
	if not os.path.exists(inputdir):
		print(inputdir," does not exists. Exit.")
		sys.exit()
		
	os.chdir(inputdir)
	print("Working inside ",inputdir)
	
	dir_list = os.listdir(inputdir)
	
	for file_name in dir_list : 
		filename = os.path.join(inputdir, file_name)
		for ivar, var in  enumerate(varibales) : 
			if var in filename :
				print(">> Reading : ", var,"from ",filename)
				lons, lats, temps_splitted, data_var = ut.read_var_from_nc(vars_names_in_netcdf[var], filename)#, read_by_chunks=True)
				print("   + Select data according to precipitant days. ")
				dates_precip_days = ut.read_dates(os.path.join(inputdir,file_name_of_precip_days))
				data_of_precip_days = ut.select_precip_days(data_var, dates_precip_days, temps_splitted)
				fld = np.array(data_of_precip_days)
				print("   + Plot mean of ",var," corresponding to the precipitant dates")
				fld_mean = np.mean(fld,axis=0)
				if var == 'MSLP':
					fld_mean = ((fld_mean*0.222441365418949)+98519.2012793173)/100
				if 'RH' in var:
					fld_mean = ((fld_mean*0.00270007307256261)+76.4849890557184)
					set_under_for_RH = True
				if 'T' in var:
					fld_mean = ((fld_mean*0.00133419765949562)+261.775997635057-273.15)# # To °C
				title= 'Average of ' +var +'  over all NAO- (k4) rainy days (RR>= 1 mm)'
				levels=np.arange(bounds[ivar][0], bounds[ivar][1], bounds[ivar][2])
				cmap = cmaps[ivar][0]
				#title='Temperature average at xxx hPa for extreme rainfall (RR>= Q99) at Laayoune'
				outdir = os.path.join(os.getcwd(),'Figs/'+var+'/')
				if not os.path.exists(outdir) : os.mkdir(outdir)
				if plotWind and var in plotWindWithvars :
					lons, lats, temps_splitted, u = ut.read_var_from_nc('u', 'U'+levWind[var]+'.nc')
					lons, lats, temps_splitted, v = ut.read_var_from_nc('v', 'V'+levWind[var]+'.nc')
					u = np.mean(u,axis=0)
					v = np.mean(v,axis=0)
					ut.PlotMap(lons, lats, fld_mean, title=title, namefig=var, outdir=outdir, levels=levels, cmap=cmap,set_under_for_RH= set_under_for_RH, plotWind=plotWind, u=u, v=v )
				else:
					ut.PlotMap(lons, lats, fld_mean, title=title, namefig=var, outdir=outdir, levels=levels, cmap=cmap,set_under_for_RH= set_under_for_RH )
				
				if plot_extreme :
					print("    + Plot ", var, " of extrme days precip with regards to the mean of ",var)
					if 'file_name_of_extreme_days' in kwargs.keys():
						file_name_of_extreme_days = kwargs['file_name_of_extreme_days']
						filename = os.path.join(inputdir, file_name_of_extreme_days)
						dates_precip_days_extreme = ut.read_dates(filename)
						print("  + Select ", var, " of extrme days precip ")
						data_of_precip_days_extrems = ut.select_precip_days(data_var, dates_precip_days_extreme, temps_splitted)
						print("  + Plot mean of ",var," corresponding to the extreme precipitant dates")
						for ixtrm, xtrm in enumerate(data_of_precip_days_extrems):
							fld_1 = fld_mean
							fld_2 = xtrm
							if var == 'MSLP':
								fld_2 = ((fld_2*0.222441365418949)+98519.2012793173)/100
							if 'RH' in var:
								fld_2 = ((fld_2*0.00270007307256261)+76.4849890557184)
								set_under_for_RH = True
							if 'T' in var:
								fld_2 = ((fld_2*0.00133419765949562)+261.775997635057-273.15) # To °C
							fld_2_minus_1 = fld_2-fld_1
							datestr = dates_precip_days_extreme[ixtrm][0]+'-' +dates_precip_days_extreme[ixtrm][1]+'-'+dates_precip_days_extreme[ixtrm][2]
							cmap = 'cmapMed_centred_on_white'
							cmap_of_diff = 'cmapMed_centred_on_white_for_diff'
							outdir = os.path.join(os.getcwd(),'Figs/'+var+'/')
							if not os.path.exists(outdir) : os.mkdir(outdir)
							levels=np.arange(bounds[ivar][0], bounds[ivar][1], bounds[ivar][2])
							cmap = cmaps[ivar][0]
							title = var + ' of the NAO- (k4) day of '+ datestr+ ' with extreme rain (RR>=Q95)'
							ut.PlotMap(lons, lats, fld_2, title=title, namefig=var+'_'+datestr, outdir=outdir, levels=levels, cmap=cmap, set_under_for_RH=set_under_for_RH)
							levels=np.arange(bounds_diff[ivar][0], bounds_diff[ivar][1], bounds_diff[ivar][2])
							cmap = cmaps_diff[ivar][0]
							title = var + ' of the NAO- (k4) day of ' + datestr + ' with extreme rain (RR>=Q95) MINUS the average of '+ var + '\n over all NAO- (k4) rainy days (RR>= 1 mm)'
							ut.PlotMap(lons, lats, fld_2_minus_1, title=title, namefig='Diff_'+var+'_'+datestr, outdir=outdir,levels=levels, cmap=cmap_of_diff)
					


		

