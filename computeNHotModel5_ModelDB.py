#
# MODIFICATION TO USE THE NEW MODEL
#
# computeNHotModel5 (ModelDB.Rds)
# ===============================
#
#
# 2018-04-23
# christian.werner@senckenberg.de
#
#
# changes:
# 2018-07-23: adapt to final 4 models
#
#
# new variables:
# precip_mm
# 

import math
import os
import sys

import numpy as np
import pandas as pd
import progressbar as pb
import xarray as xr

from optparse import OptionParser

# this is for rpy2 < 2.4
#import pandas.rpy.common as com

# import rpy bridge
import rpy2
import rpy2.robjects as robj
from rpy2.robjects import pandas2ri
pandas2ri.activate()

# progressbar style
pbar_widget = [pb.Bar('=', pb.FormatLabel("JDay: %(value)d ["), ']'),
               ' ', pb.Percentage(),
               ' ', pb.FormatLabel(" RunTime: %(elapsed)s"),
               ' ', pb.ETA()]



def main(options, args, fmodel):

    # pass some commandline args
    # TODO: check, this is actually domain and not month!!!
    #       rename appropriately

    year = int(args[0])
    if len(args) > 1:
        domain = "_" + args[1]
    else:
        domain = ""

    SMODE = False
    if options.coord != None:
        SMODE = True
        x = options.coord.split(',')
        slat = float(x[0])  # site coordinate (lat)
        slon = float(x[1])  # site coordinate (lon)

    def spatialmode():
        if SMODE==False:
            return True
        else:
            return False

    # load data
    ds_soil = xr.open_dataset( "soil/GLOBAL_WISESOIL_D1%s.nc" % domain )

    if spatialmode():
        # regular mode (full netcdf domain)
        if "lay" in ds_soil.dims:
            corg = ds_soil.TOTC.isel(lay=0).to_masked_array()
            pH   = ds_soil.PHAQ.isel(lay=0).to_masked_array()
            clay = ds_soil.CLPC.isel(lay=0).to_masked_array()
        else:
            corg = ds_soil.TOTC.isel(lev=0).to_masked_array()
            pH   = ds_soil.PHAQ.isel(lev=0).to_masked_array()
            clay = ds_soil.CLPC.isel(lev=0).to_masked_array()
    else:
        # site mode, only for specified lat, lon
        # NOTE: to stay consistent with the default mode we retain all dimensions ([] for lat, lon)
        corg = ds_soil.TOTC.isel(lay=0).sel(lat=[slat], lon=[slon], method='nearest').to_masked_array()
        pH   = ds_soil.PHAQ.isel(lay=0).sel(lat=[slat], lon=[slon], method='nearest').to_masked_array()
        clay = ds_soil.CLPC.isel(lay=0).sel(lat=[slat], lon=[slon], method='nearest').to_masked_array()


    lats  = ds_soil['lat'].copy()
    lons  = ds_soil['lon'].copy()

    # convert corg for mask
    corg = np.ma.masked_where(corg < 0, corg, copy=True)

    # copy the mask to other variables
    A_pH   = np.ma.masked_array(pH, corg.mask)
    A_clay = np.ma.masked_array(clay, corg.mask)
    soil_dims = A_pH.shape

    A_pH   = A_pH.ravel()
    A_clay = A_clay.ravel()

    ds_soil.close()

    # encoding for output
    latD = {'units': 'degrees_north', 'long_name': 'latitude'}
    lonD = {'units': 'degrees_east', 'long_name': 'longitude'}


    # loop over latitude bands [0.5deg]
    ds_climate = xr.open_dataset( "climdata_new/climdata_processed_%d%s.nc" % (year, domain) )
    ds_map     = xr.open_dataset( "prcp/MAP%s.nc" % domain )
    ds_prcp    = xr.open_dataset( "prcp/prcp_%d%s.nc" % (year, domain) )

    atime = ds_climate['time'].copy()
    clons = ds_climate['lon'].copy()

    A_output = np.zeros( (len(atime),) + soil_dims )

    jdays = atime.values.astype(np.int64, copy=False).tolist()
    # depending on the level specified on the commandline,
    # run the model for crop 1,2,3,4 or natural (0)

    doagri = np.zeros( soil_dims, 'i' )

    if options.level > 0:
        if year == 1990:
            # firstYear
            prev = 1990
        else:
            prev = year - 1

        if spatialmode():
            ds_prev = xr.open_dataset( "fertilizer/fertilizer_%d%s.nc" % (prev, domain)).sel(lev=options.level)
            ds      = xr.open_dataset( "fertilizer/fertilizer_%d%s.nc" % (year, domain)).sel(lev=options.level)

            alon = xr.open_dataset( "fertilizer/fertilizer_%d%s.nc" % (year, domain))['lon'].copy()
            alat = xr.open_dataset( "fertilizer/fertilizer_%d%s.nc" % (year, domain))['lat'].copy()

        else:
            ds_prev = xr.open_dataset( "fertilizer/fertilizer_%d%s.nc" % (prev, domain)).sel(lev=options.level, lat=[slat], lon=[slon], method='nearest')
            ds      = xr.open_dataset( "fertilizer/fertilizer_%d%s.nc" % (year, domain)).sel(lev=options.level, lat=[slat], lon=[slon], method='nearest')


        fert_prev = np.ma.zeros( A_output.shape, 'i' )
        fert      = np.ma.zeros( A_output.shape, 'i' )

        jdf1 = ds["day_fert1"].to_masked_array().astype('i')
        jdf2 = ds["day_fert2"].to_masked_array().astype('i')

        adf1 = ds["fert1_N"].to_masked_array()
        adf2 = ds["fert2_N"].to_masked_array()

        p_jdf1 = ds_prev["day_fert1"].to_masked_array().astype('i')
        p_jdf2 = ds_prev["day_fert2"].to_masked_array().astype('i')

        p_adf1 = ds_prev["fert1_N"].to_masked_array()
        p_adf2 = ds_prev["fert2_N"].to_masked_array()

        ds_prev.close()
        ds.close()


        bar = pb.ProgressBar(maxval=len(fert[0]),
                             term_width=80,
                             widgets=pbar_widget).start()

        fertcsum = np.zeros( A_output.shape )
        for ix in range(len(fert[0])):
            for jx in range(len(fert[0,0])):
                if np.ma.getmaskarray(jdf1)[ix, jx] == False:
                    jd1 = jdf1[ix,jx]
                    jd2 = jdf2[ix,jx]
                    p_jd1 = p_jdf1[ix,jx]
                    p_jd2 = p_jdf2[ix,jx]

                    fert[ jd1-1, ix, jx ] = adf1[ix,jx]
                    doagri[ix,jx] = 1

                    # same for previous year
                    fert_prev[ p_jd1-1, ix, jx ] = p_adf1[ix,jx]

                    if adf2[ix, jx ]>0:
                        fert[ jd2-1, ix, jx ] = adf2[ix,jx]

                    # same for previous year
                    if p_adf2[ix, jx]>0:
                        fert_prev[ p_jd2-1, ix, jx ] = p_adf2[ix,jx]
            bar.update(ix)
        bar.finish()


        del jdf1, jdf2, adf1, adf2, p_jdf1, p_jdf2, p_adf1, p_adf2

        # compute the aggregate (bi-weekly vars)
        cumdays = 14    # we now have bi-weekly aggregates
        for jd in range(len(jdays)):
            if jd >= cumdays:
                fertcsum[jd,:,:] = np.sum(fert[jd-cumdays:jd,:,:], axis=0)
            else:
                x = abs(jd-cumdays)
                x_prev = len(fert_prev)-x
                f_sum_prev = np.sum(fert_prev[x_prev:len(fert_prev),:,:], axis=0)
                f_sum_act  = np.sum(fert[0:jd,:,:], axis=0)
                fertcsum[jd,:,:] = f_sum_prev + f_sum_act

        del fert, fert_prev, f_sum_prev, f_sum_act

        # output to netcdf file if in regional mode
        #if spatialmode():
        #    # for stesting purposes, dump fertcsum
        #    dout = xr.Dataset()
        #    daof = xr.DataArray( fertcsum, coords=[('time',atime ), ('lat',alat, latD), ('lon',alon, lonD)] )
        #    dout['fertcsum'] = daof
        #    dout.to_netcdf('fertcsum_lev%d_%d%s.nc' % (options.level, year, domain), format='NETCDF4_CLASSIC')



    bar = pb.ProgressBar(maxval=len(jdays),
                         term_width=80,
                         widgets=pbar_widget).start()

    # default arrays
    ONES     = np.ones( A_pH.shape )
    ZEROS    = np.zeros( A_pH.shape )
    CROPMASK = doagri.ravel()

    def expand2d(inarray, size):
        return (np.repeat(np.repeat(inarray, 3, axis=0), 3, axis=1)).ravel()

    if not spatialmode():
        outRows = []

    if spatialmode():
        print "SPATIALMODE"

    for jd_pos, jd in enumerate(jdays):

        ds_slice      = ds_climate.isel(time=jd_pos) #timesubset, lat=latsubset, lon=lonsubset)
        ds_slice_prcp = ds_prcp.isel(time=jd_pos)

        if spatialmode():
            S_avg_temp = expand2d(ds_slice['tavg'].values, 3)
            S_MATEMP = expand2d(ds_slice['freezethaw'].values, 3)
            S_NCP_2week = expand2d(ds_slice['prcp_14d'].values, 3)
            S_MAP = expand2d(ds_slice['drywet_q04'].values, 3)
            S_prcp = expand2d(ds_slice_prcp['prcp'].values, 3)
        else:
            S_avg_temp  = ds_slice['tavg'].sel(lat=[slat], lon=[slon], method='nearest').values
            S_MATEMP    = ds_slice['freezethaw'].sel(lat=[slat], lon=[slon], method='nearest').values
            S_NCP_2week = ds_slice['prcp_14d'].sel(lat=[slat], lon=[slon], method='nearest').values
            S_MAP = ds_slice['drywet_q04'].sel(lat=[slat], lon=[slon], method='nearest').values
            S_prcp = ds_slice_prcp['prcp'].sel(lat=[slat], lon=[slon], method='nearest').values

        if options.level > 0:
            A_NCF = fertcsum[jd_pos].ravel()
        else:
            A_NCF = ZEROS[:]

        # new dataframe
        #print type(S_prcp), S_prcp
        #print type(A_NCF), A_NCF

        df = pd.DataFrame({'mz.one':       ONES * -1,
                           'mz.clay':      A_clay * -1,
                           'mz.pH':        A_pH * -1,
                           'mz.clay.pH':   A_clay * A_pH * -1,
                           'mz.h':         ONES * -1,
                           # IS DELTADOY == 1 ???,-nhota$delta.doy, # ????
                           'mz.h.clay':    ONES * A_clay * -1,
                           'mz.h.pH':      ONES * A_pH * -1,
                           'mz.h.clay.pH': ONES * A_pH * A_clay * -1,
                           #
                           'NCF.2week':       A_NCF,
                           'NCP.2week':       S_NCP_2week,
                           'NMAP':            S_MAP,
                           'Navg_temp':       S_avg_temp,
                           'Nprecip_mm':      S_prcp,
                           'Ncp2week_temp':   S_NCP_2week * S_avg_temp,
                           'Ncp2week_precip': S_NCP_2week * S_prcp,                   # CHECK
                           'Ntemp_MAP':       S_avg_temp * S_MAP,
                           'Nprecip_MAP':     S_prcp * S_MAP,
                           'Ntemp_precip':    S_avg_temp * S_prcp,
                           'Ncp2week_ph':     S_NCP_2week * A_pH,
                           'Nph_MAP':         A_pH * S_MAP,
                           'Nph_temp':        A_pH * S_avg_temp,
                           'Nclay_temp':      A_clay * S_avg_temp,
                           'Nclay_precip':    A_clay * S_prcp,
                           'Ncf2week_temp':   A_NCF * S_avg_temp,
                           'Ncf2week_precip': A_NCF * S_prcp})

        # call the model, select cells
        df_index = np.where(df['mz.pH'].notnull())[0]

        if options.level > 0:
            # also check that we do have agriculture
            df_index = np.intersect1d( df_index, np.where(CROPMASK == 1))
        df = df.iloc[df_index]
        #rdf = com.convert_to_r_dataframe(df)
        rdf = rpy2.robjects.pandas2ri.py2ri(df)
        r_res = robj.r.predict(fmodel, newdata=rdf, type='response')
        res = rpy2.robjects.pandas2ri.ri2py(r_res)

        if not spatialmode():
            outD = [x for x in df.values.ravel().tolist()] #+ [res]
            if jd_pos == 0:
                headerL = [x for x in df.columns.values.tolist()] + ['nhot']
                headerS = '\t'.join( headerL ) + '\n'
                outRows.append(headerS)
            outS = "\t".join( [str(x) for x in outD + res.values.tolist()]  ) + '\n'
            outRows.append( outS )


        A_output_slice = ONES * -9999
        np.put(A_output_slice, df_index, res)

        A_output_slice = A_output_slice.reshape(soil_dims)
        A_output[jd_pos] = A_output_slice

        bar.update(jd_pos)

    bar.finish()

    ds_climate.close()

    #output
    print "producing output"
    if spatialmode():
        ds_out = xr.Dataset()

        A_output = np.ma.masked_where(A_output <= -9999, A_output)

        da = xr.DataArray(A_output,
                          coords=[('time', atime), ('lat', lats, latD), ('lon', lons, lonD)])

        da.encoding.update({'complevel': 5,
                            'zlib': True,
                            'chunksizes': (10, 40, 20),
                            'shuffle':True })  # add compression

        ds_out['nhot']       = da
        ds_out.to_netcdf( "nhot_predict_lev%d_%d%s.nc" % (options.level, year, domain), format='NETCDF4_CLASSIC' )

        ds_out.close()
    else:
        open('siterun2_lev%d_%d_lat%.2f-lon%.2f.txt' % (options.level, year, slat, slon), 'w' ).write( ''.join( outRows ) )



class MyParser( OptionParser ):
    def format_epilog(self, formatter):
        return self.epilog


if __name__ == "__main__":

    parser = MyParser( "usage: %prog year domainslice", epilog=
"""

NHot Project NHot Model (v1)

Use this tool to calculate global nhot peak maps
___________________________________________
2016/01/23, christian.werner@senckenberg.de
""")

    parser.add_option("-s", "--silent", dest="silent",
                      default=False, action='store_true',
                      help="disable print statements")
    parser.add_option("-l", "--level",  dest="level",
                      type="int", default=0,
                      help="level to run (0: no fert, 1-4 fert)")
    parser.add_option("-c", "--coord", dest="coord",
                      default=None,
                      help="site location mode, run given lat,lon coord")
    parser.add_option("-m", "--model", dest="model",
                      default="model1.Rds",
                      help="stats model")


    (options, args) = parser.parse_args()

    print "calculating nhot for level %d (0=natural, 1-4=agri)" % options.level

    if len(args) != 2:
        print parser.usage

    if options.silent:
        f = open(os.devnull, 'w')
        sys.stdout = f
        sys.stderr = f

    # new model iteration provided by David
    fmodel = robj.r.readRDS(options.model)
    print fmodel

    main(options, args, fmodel)
