#
# create_nhot_climbasedvars.py
# ============================
#
# 2016-08-15 (updated)
# christian.werner@senckenberg.de
#
# - now used bi-weekley PRECsum

import xarray as xr
import math, progressbar
import numpy as np
import pandas as pd
import sys

#cdo griddes ifile > mygrid
## edit mygrid and set xfirst to the new value
#cdo setgrid,mygrid ifile ofile

#lats = np.arange(-89.75, 90, 0.5)
#lons = np.arange(-179.75, 180, 0.5)


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def moving_sum(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

YEARS = range(1990,2001)

if len(sys.argv) > 1:
    YEARS = [int(x) for x in sys.argv[1:]]


def do_this(control): #, parts=['all']):
    if 'all' in parts:
        print("Process all sections")
        return True
    elif control in parts:
        print("Processing:", control)
        return True
    else:
        return False


for year in YEARS:

    parts=['all'] #'freezethaw', 'freezedays']

    print('processing year: %d' % year)

    ds_out = xr.Dataset()
    

    prev = year-1

    ds_prev = xr.open_dataset( "/Users/cwerner/Dropbox/climate_%d.nc" % prev)[['prcp','tmax','tmin']].load()
    ds      = xr.open_dataset( "/Users/cwerner/Dropbox/climate_%d.nc" % year)[['prcp','tmax','tmin']].load()


    lats = ds['latitude'].copy()
    lons = ds['longitude'].copy()
    atime = ds['time'].copy()

    #precip
    prcp_prev = ds_prev['prcp'].squeeze()
    prcp      = ds['prcp'].squeeze()

    tavg_prev = ( ds_prev['tmax'].squeeze() + ds_prev['tmin'].squeeze() ) * 0.5
    tavg      = ( ds['tmax'].squeeze()      + ds['tmin'].squeeze() )      * 0.5

    jdays = range(0, len(prcp))
    
    tavg = tavg.to_masked_array()

    prcpcsum    = np.empty( (len(jdays), 720, 1440) )
    prcpcsum2   = np.empty( (len(jdays), 720, 1440) )
    freezethaw  = np.empty( (len(jdays), 720, 1440) )
    drywet_q04  = np.empty( (len(jdays), 720, 1440) )
    drywet_q08  = np.empty( (len(jdays), 720, 1440) )

    frzdegdays  = np.empty( (len(jdays), 720, 1440) )
    frzdegtemps = np.empty( (len(jdays), 720, 1440) )

    if do_this('precip'):

        #print 'cumulative precip'
        bar = progressbar.ProgressBar(maxval=len(jdays), term_width=80, \
                    widgets=[progressbar.Bar('=', ' %s [' % "cum precip (7d)".rjust(20), ']'), ' ', progressbar.Percentage()]).start()

        # 2-week aggregate now
        cumdays = 7

        for jd in jdays:
            if jd >= cumdays:
                prcpcsum[jd,:,:] = np.sum(prcp[jd-cumdays:jd,:,:], axis=0)
            else:
                # sum 
                x = abs(jd-cumdays)
                x_prev = len(prcp_prev)-x
                sum_prev = np.sum(prcp_prev[x_prev:len(prcp_prev),:,:], axis=0)
                sum_act  = np.sum(prcp[0:jd,:,:], axis=0)
                prcpcsum[jd,:,:] = sum_prev + sum_act

            bar.update(jd)
        
        bar.finish()


        #print 'cumulative precip'
        bar = progressbar.ProgressBar(maxval=len(jdays), term_width=80, \
                    widgets=[progressbar.Bar('=', ' %s [' % "cum precip (14d)".rjust(20), ']'), ' ', progressbar.Percentage()]).start()

        # 2-week aggregate now
        cumdays = 14

        for jd in jdays:
            if jd >= cumdays:
                prcpcsum2[jd,:,:] = np.sum(prcp[jd-cumdays:jd,:,:], axis=0)
            else:
                # sum 
                x = abs(jd-cumdays)
                x_prev = len(prcp_prev)-x
                sum_prev = np.sum(prcp_prev[x_prev:len(prcp_prev),:,:], axis=0)
                sum_act  = np.sum(prcp[0:jd,:,:], axis=0)
                prcpcsum2[jd,:,:] = sum_prev + sum_act

            bar.update(jd)
        
        bar.finish()

    if do_this('freezethaw'):

        bar = progressbar.ProgressBar(maxval=len(jdays), term_width=80, \
                widgets=[progressbar.Bar('=', ' %s [' % "freeze-thaw".rjust(20), ']'), ' ', progressbar.Percentage()]).start()

        # 2-week aggregate now
        cumdays = 14

        for jd in jdays:
    
            if jd >= cumdays:
                subset = tavg[jd-7:jd,:,:]
                bool_subset = subset > -3

                subset = tavg[jd-7:jd,:,:]

            else:
                x = abs(jd-7)
                x_prev = len(prcp_prev)-x

                subset_prev = tavg_prev[x_prev:len(tavg_prev),:,:]
                subset_act  = tavg[0:jd,:,:]
                subset = np.vstack([subset_prev, subset_act])

            bool_subset = subset > -3
                    
            M = np.average(bool_subset, axis=0)
            I = np.average(bool_subset[0:3], axis=0) < np.average(bool_subset[-3:], axis=0)
            F = M * (1 - M) * I

            freezethaw[jd,:,:] = F


            bar.update(jd)

        bar.finish()


    if do_this('freezedays'):


        # freeze days / temps
        # use this to precompute the running vars
        jdays2 = list(range(-60,-1)) + list(jdays)

        bar = progressbar.ProgressBar(maxval=len(jdays2), term_width=80, \
                widgets=[progressbar.Bar('=', ' %s [' % "freeze-thaw".rjust(20), ']'), ' ', progressbar.Percentage()]).start()


        # frz deg days counter
        non_frz_days = np.zeros_like(tavg[0,:,:])
        _frzdegdays  = np.zeros_like(tavg[0,:,:]) 
        _frzdegtemps = np.zeros_like(tavg[0,:,:])

        threshold_temp = 0
        threshold_days = 14



        for cnt, jd in enumerate(jdays2):
            if jd < 0:
                prev_jd = len(tavg_prev) + jd 
                cur_tavg = tavg_prev[prev_jd,:,:]
            else:
                cur_tavg = tavg[jd,:,:]
            
            non_frz_days = np.where(cur_tavg > threshold_temp, non_frz_days + 1, 0)
            _frzdegdays  = np.where(non_frz_days >= threshold_days, 0, _frzdegdays)
            _frzdegtemps = np.where(non_frz_days >= threshold_days, 0, _frzdegtemps)

            # increment
            _frzdegdays  = np.where( (non_frz_days < threshold_days) & (cur_tavg <= 0), _frzdegdays + 1, _frzdegdays)
            _frzdegtemps = np.where( (non_frz_days < threshold_days) & (cur_tavg <= 0), _frzdegtemps + cur_tavg[:,:], _frzdegtemps)

            # if we are in the actual year, write to out array
            if jd >= 0:
                frzdegdays[jd, :, :]  = _frzdegdays
                frzdegtemps[jd, :, :] = _frzdegtemps

            bar.update(cnt)

        bar.finish()


    if do_this('drywet'):

        bar = progressbar.ProgressBar(maxval=len(jdays), term_width=80, \
                widgets=[progressbar.Bar('=', ' %s [' % "dry-wet".rjust(20), ']'), ' ', progressbar.Percentage()]).start()

        for jd in jdays:

            if jd >= 28:
                P   = prcp[jd-28:jd,:,:]

            else:
                x = abs(jd-28)
                x_prev = len(prcp_prev)-x

                subset_prev = prcp_prev[x_prev:len(prcp_prev),:,:]
                subset_act  = prcp[0:jd,:,:]
                P           = np.vstack([subset_prev, subset_act])



            for q in [0.8]:

                Q = np.zeros( P.shape )

                med = np.median( P, axis=0 )

                B = np.where( P <= med, np.where(med == 0, 0, P / med), 1.0)

                Q = ((1.0 - q) / (1.0 - q**28)) * np.sum([q**j * B[27-j,:,:] for j in range(28)], axis=0)

                # this is actually the last weeks (reverse order in list)
                I = np.sum(P[0:14], axis=0) < np.sum(P[-14:], axis=0)

                D = Q * (1-Q)*I

                #if q == 0.4: drywet_q04[jd,:,:] = D
                if q == 0.8: drywet_q08[jd,:,:] = D

            bar.update(jd)

        bar.finish()
    
    print('output')



    COORDS = [('time', atime), ('lat', lats), ('lon', lons)]

    add_attrs = {'complevel':5,
                        'zlib':True,
                        'chunksizes': (10, 40, 20),
                        'shuffle':True }


    # test output
    ##da5a = xr.DataArray(frzdegdays, coords=COORDS)
    ##da5b = xr.DataArray(frzdegtemps, coords=COORDS)
    ##
    ##da5a.encoding.update( add_attrs )
    ##da5b.encoding.update( add_attrs )

    ##ds_out['frzdegdays'] = da5a
    ##ds_out['frzdegtemps'] = da5b

    ##ds_out.to_netcdf( "dummy_%d.nc" % year, format='NETCDF4_CLASSIC' )

    ##exit()

    # regular output

    da0 = xr.DataArray(tavg,       coords=COORDS)
    da1a = xr.DataArray(prcpcsum,   coords=COORDS)
    da1b = xr.DataArray(prcpcsum2,   coords=COORDS)
    da2 = xr.DataArray(freezethaw, coords=COORDS)
    #da3 = xr.DataArray(drywet_q04, coords=COORDS)
    da4 = xr.DataArray(drywet_q08, coords=COORDS)

    da5a = xr.DataArray(frzdegdays, coords=COORDS)
    da5b = xr.DataArray(frzdegtemps, coords=COORDS)



    da0.encoding.update( add_attrs )
    da1a.encoding.update( add_attrs )
    da1b.encoding.update( add_attrs )
    da2.encoding.update( add_attrs )
    #da3.encoding.update( add_attrs )
    da4.encoding.update( add_attrs )

    da5a.encoding.update( add_attrs )
    da5b.encoding.update( add_attrs )

    ds_out['tavg']       = da0
    ds_out['prcp_7d']    = da1a
    ds_out['prcp_14d']   = da1b
    ds_out['freezethaw'] = da2
    #ds_out['drywet_q04'] = da3
    ds_out['drywet_q08'] = da4

    ds_out['frzdegdays'] = da5a
    ds_out['frzdegtemps'] = da5b


    ds_out.to_netcdf( "climdata_new2/climdata_processed_%d.nc" % year, format='NETCDF4_CLASSIC' )


exit()




















# ~/Downloads/Global_N_application
ds = xr.Dataset()


Dlut = {}

for latcnt, lat in enumerate(lats):
    for loncnt, lon in enumerate(lons):
        Dlut[(lat, lon)] = (latcnt, loncnt)


crops = ["barley","maiz","millet","sorghum","soybean","wheat"]

for yearF in range(1985, 2009):
    print('year', yearF)

    year = str(yearF)[2:]

    da     = xr.DataArray(np.zeros( (len(lats), len(lons)) ), coords=[('lat', lats), ('lon', lons)])
    da_cnt = xr.DataArray(np.zeros( (len(lats), len(lons)) ), coords=[('lat', lats), ('lon', lons)])


    da_crops = []

    for crop in crops:

        da_crop = xr.DataArray(np.zeros( (len(lats), len(lons)) ), coords=[('lat', lats), ('lon', lons)])

        
        f = "Global_N_application/%s_fert_allcells.csv" % crop
        df = pd.read_csv( f, error_bad_lines=False, na_values=-9999, sep="," )
        for index, row in df.iterrows():
            val = row['year_%s' % year]
            if math.isnan(val) == False:
                idx, jdx = Dlut[(row['latitude'], row['longitude'])] 
                da_crop[idx,jdx] = (val * 10)
                da[idx,jdx]      += (val * 10)
                da_cnt[idx,jdx]  += 1

        da_crops.append( da_crop )

    # ---
    da = da / da_cnt.where((da_cnt > 0)).to_masked_array(copy=True)

    ds['avg'] = da

    for lcrop, crop in enumerate(crops):
        ds[crops[lcrop]] = da_crops[lcrop]

    ds.to_netcdf( "fert_%d.nc" % yearF, format='NETCDF4_CLASSIC' )




