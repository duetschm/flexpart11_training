#!/usr/bin/env python

'''
    plot_FLEX_binary.py

    Using FLEXPART binary output files, creates 2D plots of the data.

    Original Author contact information:

    Delia Arnold
    delia.arnold-arias@zamg.ac.at

    Marie Danielle Mulder
    marie.mulder@zamg.ac.at

    Christian Maurer
    christian.maurer@zamg.ac.at

    Disclaimer:

    This software is presented freely to the public with
    no restrictions on its use.  However, it would be
    appreciated if any use of the software or methods in
    part or in full acknowledges the source. This software 
    has been developed for flexpart workshop 2019 @ ZAMG.
    No liabilty taken by developers.

'''
import ast
import sys
import os
import matplotlib.pyplot as plt
import imageio
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import ticker
import numpy as np

# Importing various functions to read in
# and work with FLEXPART binary output
#   --- based on pflexible (J. BUkhardt)
#   --- possible to use f2py for faster reading
import readHeader_py3_v11 as readH  
import read_grid_py3_v11 as readG

import warnings
warnings.filterwarnings("ignore")


def main(argv = sys.argv):

    # TO EDIT 
    plot_thres = 1.0E-30
    # plot_thres = 0.0
    # DEFAULT values
    
    num_SPECIES = [1]
    factor = 1.e-12
    LEVEL = 0
    RELEASE = 0
    projection = 'cyl'
    depo = False
    plotting_region = []
    style = 'mesh'
    plot_markers = False

    print ('!!!!!!!!!!!!!!!!!DEFAULT values!!!!!!!!!!!!!!')
    print ('Species to be considered: ', num_SPECIES)
    print ('Multiplication factor for output: ', factor)
    print ('Selected output levels ', LEVEL)
    print ('Selected release: ', RELEASE)
    print ('All available date-times are considered')
    print ('Projection: ', projection)
    print ('Deposition: ', depo)
    print ('Plotting region according to OUTGRID')
    print ('Plotting style: ', style)
    print ('Plot markers: ', plot_markers)

    # Minimum 3 arguments (including script name)
    if len(argv) < 3:
        print ('Argvs: <Flexpart output dir> <nest> [<numSPECIES> <LEVEL> <RELEASE> <d_list/alldates> <projection> <depo> <plotting_region (ll_lon,ur_lon,ll_lat,ur_lat)> <style> <plot_markers>]')
        print ('Example: ./output True (or False) 1 0 0 20170925000000,20170925010000 (or alldates) cyl (or lcc) True/False 110,120,20,25 contour True (or False)')
        print ('Marker locations will be taken from file markers.dat')
        print ('MIND: numSPECIES starts with 1, LEVEL and RELEASE with 0 (python style)! d_list and plotting_region elements have to be comma-separated!')      
        sys.exit(0)

    FLEXPARTOutputDir = argv[1]
    isNested = ast.literal_eval(argv[2])
    # Read the general header information common to all files
    H = readH.read_header(FLEXPARTOutputDir, nested = isNested)
    print ('Arguments entered: ',FLEXPARTOutputDir,isNested)

    if len(argv) >= 4:
        num_SPECIES = list(argv[3].split(','))
        LEVEL = int(argv[4])
        RELEASE = int(argv[5])
        if argv[6] != 'alldates':
            d_list = list(argv[6].split(','))
        else:
            # if 'alldates': take list of available date-times
            d_list = H.available_dates
        projection = argv[7]
        depo = ast.literal_eval(argv[8])
        plotting_region = list(argv[9].split(','))
        style = argv[10]
        plot_markers = ast.literal_eval(argv[11])
        
        if plot_markers == True:
            lon_stat = []
            lat_stat = []
            name_stat = []
            if os.path.isfile('./markers.dat'):
                print ('Reading marker locations')
                stations = open('./markers.dat', "r")
                lines = stations.readlines()
                for i in range(len(lines)):
                    print (lines[i].split(' '))
                    lon_stat.append(float(lines[i].split(' ')[1]))
                    lat_stat.append(float(lines[i].split(' ')[2]))
                    name_stat.append(lines[i].split(' ')[0])
                stations.close()
            else:
                print ('ERROR: File markers.dat not found!!!')
                sys.exit(0)

        
        print ('DEFAULT values UPDATED to: ')
        print (num_SPECIES, LEVEL, RELEASE, d_list, projection, depo, plotting_region, style, plot_markers)
        
        
    savepath = FLEXPARTOutputDir

    # Get something needed from header
    outputFiletype = H.unit  # This will determine name of file (conc, time, pptv)
    simDirection = H.direction  # backward or forward

    dates_file = FLEXPARTOutputDir + 'dates'
    instance = 0

    with open(dates_file) as dfile:
        # substitute dfile-list with d_list-list if given
        if len(d_list) > 0:
            dfile = d_list
        for date in dfile:
            current_time = date.strip() 
            print ('Plotting timestamp: ' + str(current_time))

            if os.path.isfile(FLEXPARTOutputDir + 'grid_' + outputFiletype + '_' + str(current_time)+'_001'):
                # Sum over species
                theHorizslice=np.zeros((int(H.numxgrid), int(H.numygrid)))
                theHorizslice_wetdep=np.zeros((int(H.numxgrid), int(H.numygrid)))
                theHorizslice_drydep=np.zeros((int(H.numxgrid), int(H.numygrid)))
                theHorizslice_depo=np.zeros((int(H.numxgrid), int(H.numygrid)))
                for sp in num_SPECIES:
                    if isNested == False:
                        f =  FLEXPARTOutputDir + 'grid_' + outputFiletype + '_' + str(current_time) + '_' + str(sp).zfill(3)
                    else:
                        f =  FLEXPARTOutputDir + 'grid_' + outputFiletype + '_nest_' + str(current_time) + '_' + str(sp).zfill(3)
                    print (f)
                    grid, wetdep, drydep, itime = readG._readgridBF(H, f)
                    theHorizslice+=grid[:,:,LEVEL,RELEASE,0] # last dimension is AGECLASS -> hard coded for one AGECLASS
                    theHorizslice_wetdep+=wetdep[:,:,RELEASE,0]
                    theHorizslice_drydep+=drydep[:,:,RELEASE,0]
                    theHorizslice_depo=theHorizslice_wetdep+theHorizslice_drydep

                theHorizslice=np.transpose(theHorizslice)
                theHorizslice_depo=np.transpose(theHorizslice_depo)
                theHorizslice_wetdep=np.transpose(theHorizslice_wetdep)
                theHorizslice_drydep=np.transpose(theHorizslice_drydep)
                if simDirection == 'forward':
                    # Convert 
                    theHorizslice = theHorizslice*factor
                    theHorizslice_depo=theHorizslice_depo*factor
                    theHorizslice_wetdep=theHorizslice_wetdep*factor
                    theHorizslice_drydep=theHorizslice_drydep*factor

                print ('Max residence time or concentration: ' +str(theHorizslice.max().max()))
                print ('Min residence time or concentration: ' +str(theHorizslice.min().min()))
                if simDirection == 'forward' and depo == True:
                    print ('Max total deposition: ' + str(theHorizslice_depo.max().max()))
                    print ('Min total deposition: ' + str(theHorizslice_depo.min().min()))
                    print ('Max wet deposition: ' + str(theHorizslice_wetdep.max().max()))
                    print ('Min wet deposition: ' + str(theHorizslice_wetdep.min().min()))
                    print ('Max dry deposition: ' + str(theHorizslice_drydep.max().max()))
                    print ('Min dry deposition: ' + str(theHorizslice_drydep.min().min()))

                
                if instance == 0:

                    instance=1

                    # Get number of grid points, longitudes and latitudes
                    numygrid=H.numygrid
                    numxgrid=H.numxgrid
                    lat=H.latitude
                    lon=H.longitude

                    # use confined plotting region if given
                    if len(plotting_region) > 0:
                        lon_ll=float(plotting_region[0])
                        lon_ur=float(plotting_region[1])
                        lat_ll=float(plotting_region[2])
                        lat_ur=float(plotting_region[3])
                    # use default OUTGRID plotting region
                    else: 
                        lat_ll=lat[0]
                        lon_ll=lon[0]
                        lat_ur=lat[numygrid-1]
                        lon_ur=lon[numxgrid-1]

                    # Set up the projection
                    if lon_ll > lon_ur:  # transgression of date boundary
                        lon_ur += 360.0
                        centr_long = 180.0
                    else:
                        centr_long = 0.0

                    if projection == 'cyl':
                        projccrs = ccrs.PlateCarree(central_longitude=centr_long)
                    elif projection == 'lcc':
                        projccrs = ccrs.LambertConformal(central_longitude=lon_ll+(lon_ur-lon_ll)/2.0,\
                        central_latitude=lat_ll+(lat_ur-lat_ll)/2.0)
                    else:
                       raise ValueError(
                       "Projection set to %s, but must be either Lambert 'lcc' or cylindrical 'cyl'",
                       projection)

  
                    if (lon_ll <= -175.0 and lon_ur >= 175.0) and (lat_ll <= -80.0 or lat_ur >= 80.0):
                        delta=30.
                    else:
                        delta=10.
                

                fig = plt.figure(figsize=(12,12))
                domain_map = fig.add_subplot(1, 1, 1, projection=projccrs)
                domain_map.set_extent([lon_ll, lon_ur, lat_ll, lat_ur],crs=ccrs.PlateCarree())
                domain_map.add_feature(cfeature.COASTLINE,zorder=10)
                domain_map.add_feature(cfeature.BORDERS,zorder=10)
                domain_map.add_feature(cfeature.LAKES, alpha=0.5)
                domain_map.add_feature(cfeature.RIVERS)
                domain_map.add_feature(cfeature.LAND)
                parallels = np.arange(lat_ll,lat_ur,delta)
                meridians = np.arange(lon_ll,lon_ur,delta)
                gl = domain_map.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,\
                linewidth=1, color='gray', alpha=0.5, xlocs=meridians, ylocs=parallels,\
                x_inline=False, y_inline=False)

                colors =  [('White')]+[('White')]+[('White')]+[('LavenderBlush')]+[('LavenderBlush')]+[('LavenderBlush')]+[('PowderBlue')]+[('PowderBlue')]+[('LightBlue')]+[('LightBlue')]+[('LightSkyBlue')]+[('LightSkyblue')]+[('SkyBlue')]+[('Skyblue')]+[('DodgerBlue')]+[('DodgerBlue')]+[('RoyalBlue')]+[('RoyalBlue')]+[('MediumBlue')]+[('MediumBlue')]+[(plt.cm.jet(i)) for i in range(20,256)]
                delia_cmap2 = LinearSegmentedColormap.from_list('new_map', colors, N=256)
       
                # Mark release point
                releaseXLOC = (H.xp1[RELEASE] +    \
                H.xp2[RELEASE] ) / 2.0
                releaseYLOC = (H.yp1[RELEASE] +    \
                H.yp2[RELEASE] ) / 2.0

                # additional markers
                if plot_markers == True:
                    domain_map.plot(lon_stat, lat_stat, transform=ccrs.PlateCarree(), \
                            marker = 'o', linestyle = 'None', color = 'darkred', markersize=6)
                    for i in range(len(lon_stat)): # lable markers
                        t = domain_map.text(lon_stat[i]-2.0, lat_stat[i]+1.0,\
                        name_stat[i], clip_on=True, transform=ccrs.PlateCarree())
                        t.clipbox = domain_map.bbox

                
                domain_map.plot(releaseXLOC, releaseYLOC, 'r^', markersize=8, transform=ccrs.PlateCarree())
                if theHorizslice.max().max() >= plot_thres:
                    if style == 'contour':
                        cs = domain_map.contourf(lon, lat, theHorizslice,\
                                cmap=delia_cmap2, locator=ticker.LogLocator(), zorder=1, transform=ccrs.PlateCarree())
                    else:
                        cs = domain_map.pcolormesh(lon, lat, theHorizslice, \
                                cmap=delia_cmap2, vmin=0., vmax=theHorizslice.max().max(), zorder=1, transform=ccrs.PlateCarree())
                    cbar = plt.colorbar(cs, location='bottom', pad=0.1)
                    if simDirection == 'backward':
                        cbar.set_label('Residence Time [s]')
                    else:
                        cbar.set_label('Concentration [kg/m3] or [Bq/m3]')
                plt.title('2D slice at time '+ str(current_time))
                # sum species if several species indices are entered
                if len (num_SPECIES) > 1:
                    num_SPECIES_help = 'SUM'
                else:
                    num_SPECIES_help = num_SPECIES[0]
                filename = savepath + '2D_' + str(current_time) + '_species' + str(num_SPECIES_help) + '_release' + str(RELEASE) + '_level' + str(LEVEL)
                if isNested == False:
                    filename+='.png'
                else:
                    filename+='_nest.png'
                plt.savefig(filename, bbox_inches = 'tight')
                print ('Saving figure as ' + filename)
                plt.close()
                if depo == True:
                    if theHorizslice_depo.max().max() >= plot_thres:
                        # total depo
                        fig = plt.figure(figsize=(12,12))
                        domain_map = fig.add_subplot(1, 1, 1, projection=projccrs)
                        domain_map.set_extent([lon_ll, lon_ur, lat_ll, lat_ur],crs=ccrs.PlateCarree())
                        domain_map.add_feature(cfeature.COASTLINE,zorder=10)
                        domain_map.add_feature(cfeature.BORDERS,zorder=10)
                        domain_map.add_feature(cfeature.LAKES, alpha=0.5)
                        domain_map.add_feature(cfeature.RIVERS)
                        domain_map.add_feature(cfeature.LAND)
                        gl = domain_map.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,\
                        linewidth=1, color='gray', alpha=0.5, xlocs=meridians, ylocs=parallels,\
                        x_inline=False, y_inline=False)
                        domain_map.plot(releaseXLOC, releaseYLOC, 'r^', markersize=8, transform=ccrs.PlateCarree())
                        if style == 'contour':
                            cs = domain_map.contourf(lon, lat, theHorizslice_depo,\
                                    cmap=delia_cmap2, locator=ticker.LogLocator(), zorder=1, transform=ccrs.PlateCarree())
                        else:
                            cs = domain_map.pcolormesh(lon, lat, theHorizslice_depo,\
                                    cmap=delia_cmap2, vmin=0., vmax=theHorizslice_depo.max().max(), zorder=1, transform=ccrs.PlateCarree())
                        cbar = plt.colorbar(cs, location='bottom', pad=0.1)
                        cbar.set_label('Total deposition [kg/m2] or [Bq/m2]')
                        plt.title('2D total deposition slice at time '+ str(current_time))
                        filename = savepath + '2D_depo_' + str(current_time) + '_species' + str(num_SPECIES_help) + '_release' + str(RELEASE) + '_level' + str(LEVEL)
                        if isNested == False:
                            filename+='.png'
                        else:
                            filename+='_nest.png'
                        plt.savefig(filename, bbox_inches = 'tight')
                        print ('Saving figure as ' + filename)
                        plt.close()
                    if theHorizslice_wetdep.max().max() >= plot_thres:
                        # wet depo
                        fig = plt.figure(figsize=(12,12))
                        domain_map = fig.add_subplot(1, 1, 1, projection=projccrs)
                        domain_map.set_extent([lon_ll, lon_ur, lat_ll, lat_ur],crs=ccrs.PlateCarree())
                        domain_map.add_feature(cfeature.COASTLINE,zorder=10)
                        domain_map.add_feature(cfeature.BORDERS,zorder=10)
                        domain_map.add_feature(cfeature.LAKES, alpha=0.5)
                        domain_map.add_feature(cfeature.RIVERS)
                        domain_map.add_feature(cfeature.LAND)
                        gl = domain_map.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,\
                        linewidth=1, color='gray', alpha=0.5, xlocs=meridians, ylocs=parallels,\
                        x_inline=False, y_inline=False)
                        domain_map.plot(releaseXLOC, releaseYLOC, 'r^', markersize=8, transform=ccrs.PlateCarree())
                        if style == 'contour':
                            cs = domain_map.contourf(lon, lat, theHorizslice_wetdep, \
                                    cmap=delia_cmap2, locator=ticker.LogLocator(), zorder=1, transform=ccrs.PlateCarree())
                        else:
                            cs = domain_map.pcolormesh(lon, lat, theHorizslice_wetdep, \
                                    cmap=delia_cmap2, vmin=0., vmax=theHorizslice_wetdep.max().max(), zorder=1, transform=ccrs.PlateCarree())
                        cbar = plt.colorbar(cs, location='bottom', pad=0.1)
                        cbar.set_label('Wet deposition [kg/m2] or [Bq/m2]')
                        plt.title('2D wet deposition slice at time '+ str(current_time))
                        filename = savepath + '2D_wetdep_' + str(current_time) + '_species' + str(num_SPECIES_help) + '_release' + str(RELEASE) + '_level' + str(LEVEL)
                        if isNested == False:
                            filename+='.png'
                        else:
                            filename+='_nest.png'
                        plt.savefig(filename, bbox_inches = 'tight')
                        print ('Saving figure as ' + filename)
                        plt.close()
                    if theHorizslice_drydep.max().max() >= plot_thres:
                        # dry depo
                        fig = plt.figure(figsize=(12,12))
                        domain_map = fig.add_subplot(1, 1, 1, projection=projccrs)
                        domain_map.set_extent([lon_ll, lon_ur, lat_ll, lat_ur],crs=ccrs.PlateCarree())
                        domain_map.add_feature(cfeature.COASTLINE,zorder=10)
                        domain_map.add_feature(cfeature.BORDERS,zorder=10)
                        domain_map.add_feature(cfeature.LAKES, alpha=0.5)
                        domain_map.add_feature(cfeature.RIVERS)
                        domain_map.add_feature(cfeature.LAND)
                        gl = domain_map.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,\
                        linewidth=1, color='gray', alpha=0.5, xlocs=meridians, ylocs=parallels,\
                        x_inline=False, y_inline=False)
                        domain_map.plot(releaseXLOC, releaseYLOC, 'r^', markersize=8, transform=ccrs.PlateCarree())
                        if style == 'contour':
                            cs = domain_map.contourf(lon, lat, theHorizslice_drydep, \
                                    cmap=delia_cmap2, locator=ticker.LogLocator(), zorder=1, transform=ccrs.PlateCarree())
                        else:
                            cs = domain_map.pcolormesh(lon, lat, theHorizslice_drydep,\
                                    cmap=delia_cmap2, vmin=0., vmax=theHorizslice_drydep.max().max(), zorder=1, transform=ccrs.PlateCarree())
                        cbar = plt.colorbar(cs, location='bottom', pad=0.1)
                        cbar.set_label('Dry deposition [kg/m2] or [Bq/m2]')
                        plt.title('2D dry deposition slice at time '+ str(current_time))
                        filename = savepath + '2D_drydep_' + str(current_time) + '_species' + str(num_SPECIES_help) + '_release' + str(RELEASE) + '_level' + str(LEVEL)
                        if isNested == False:
                            filename+='.png'
                        else:
                            filename+='_nest.png'
                        plt.savefig(filename, bbox_inches = 'tight')
                        print ('Saving figure as ' + filename)
                        plt.close()


            else:
                print ('File does not exist: ' + str(FLEXPARTOutputDir) + 'grid_' + outputFiletype + '_(nest_)' + str(current_time) + '_001')
    print ('SUCCESSFULLY FINISHED PLOTTING OF FLEXPART BINARY OUTPUT!')

#----------------------------------------------
if __name__ == "__main__":

    main()

























