############################################################## 
# Date: 10/01/16
# Name: Alek Petty
# Description:  Functions/classes used by the BG currents scripts
import numpy as np
from pylab import *
from scipy.io import netcdf
import numpy.ma as ma
from scipy.interpolate import griddata
from glob import glob

def defGrid(m, dxRes=50000):
    nx = int((m.xmax-m.xmin)/dxRes)+1; ny = int((m.ymax-m.ymin)/dxRes)+1
    gridStr=str(int(dxRes/1000))+'km'
    lonsG, latsG, xptsG, yptsG = m.makegrid(nx, ny, returnxy=True)

    return lonsG, latsG, xptsG, yptsG

def getKimuraGrid(m, dataPath):
    lonlat = loadtxt(dataPath+'lonlat_gt60.xy', unpack=True)
    lon = lonlat[0]
    lat = lonlat[1]
    xptsM, yptsM=m(lon, lat)
    return xptsM, yptsM, lon, lat
    
def getGriddedKimura(m, file, lon, xptsG, yptsG, lonsG, latsG, xptsM, yptsM):
    f = loadtxt(file, unpack=True)
    
    # Comes in xy coordinates so need to rotate to UV
    xvel= f[0]
    yvel = f[1]
    alpha = lon*pi/180.
    uvel = yvel*sin(alpha) + xvel*cos(alpha)
    vvel = yvel*cos(alpha) - xvel*sin(alpha) 

    # Re-grid data
    uvelG = griddata((xptsM, yptsM),uvel, (xptsG, yptsG), method='linear')
    vvelG = griddata((xptsM, yptsM),vvel, (xptsG, yptsG), method='linear')

    # Rotate data onto new grid
    xvelG,yvelG = m.rotate_vector(uvelG,vvelG,lonsG,latsG)
    xvelG=ma.masked_invalid(xvelG)
    yvelG=ma.masked_invalid(yvelG)
    return xvelG, yvelG

def getGriddedFowlerFromDaily(m, files, lon, xptsM, yptsM, xptsG, yptsG, lonsG, latsG):
    
    uvelD=ma.masked_all((size(files), lon.shape[0], lon.shape[1]))
    vvelD=ma.masked_all((size(files), lon.shape[0], lon.shape[1]))
    # print 'uvel', uvelD.shape
    x=0
    for file in files:
        fd = open(file, 'rb')
        motionDat = fromfile(file=fd, dtype='<i2')
        motionDat = reshape(motionDat, [361, 361, 3])

        xt = motionDat[:, :, 0]/1000.
        yt = motionDat[:, :, 1]/1000.         
        q = motionDat[:, :, 2]/1000.

        mask = where((q<=0) | (q>1), 0, 1)

        xt = ma.masked_where(mask<0.5, xt)
        yt = ma.masked_where(mask<0.5, yt)

        # Comes in xy coordinates so need to rotate to UV
        #xvel= ma.masked_where(np.isnan(xt), xt)
        #yvel = ma.masked_where(np.isnan(yt), yt)
        #xvel=f[0]
        #yvel=f[1]

        alpha = lon*pi/180.
        uvelT = yt*sin(alpha) + xt*cos(alpha)
        vvelT = yt*cos(alpha) - xt*sin(alpha) 
        uvelD[x]=uvelT
        vvelD[x]=vvelT
        x+=1
        # COULD ROTATE HERE AND DO CURL OF DAILY VARIABLES.
        
        #vvelD=vstack([vvelD, vvelT])
    uvel=ma.mean(uvelD, axis=0)
    vvel=ma.mean(vvelD, axis=0)
    #print uvel
    #if we want to set masked values back to nan for gridding purposes
    uvel[where(ma.getmask(uvel))]=np.nan
    vvel[where(ma.getmask(vvel))]=np.nan
    #print uvel
    # Re-grid data
    # print uvel.flatten().shape, xptsM.flatten().shape, xptsG.shape
    uvelG = griddata((xptsM.flatten(), yptsM.flatten()),uvel.flatten(), (xptsG, yptsG), method='linear')
    vvelG = griddata((xptsM.flatten(), yptsM.flatten()),vvel.flatten(), (xptsG, yptsG), method='linear')

    # Rotate data onto new grid
    xvelG,yvelG = m.rotate_vector(uvelG,vvelG,lonsG,latsG)
    xvelG=ma.masked_invalid(xvelG)
    yvelG=ma.masked_invalid(yvelG)
    return xvelG, yvelG

def getGriddedFowlerCurlFromDaily(m, files, lon, xptsM, yptsM, xptsG, yptsG, lonsG, latsG, dxRes):
    
    xvelG=ma.masked_all((size(files), lonsG.shape[0], lonsG.shape[1]))
    yvelG=ma.masked_all((size(files), lonsG.shape[0], lonsG.shape[1]))
    curlG=ma.masked_all((size(files), lonsG.shape[0], lonsG.shape[1]))
    #print 'uvel', uvelD.shape
    x=0
    for file in files:

        fd = open(file, 'rb')
        motionDat = fromfile(file=fd, dtype='<i2')
        motionDat = reshape(motionDat, [361, 361, 3])

        xt = motionDat[:, :, 0]/1000.
        yt = motionDat[:, :, 1]/1000.         
        q = motionDat[:, :, 2]/1000.

        mask = where((q<=0) | (q>1), 0, 1)
        xt = ma.masked_where(mask<0.5, xt)
        yt = ma.masked_where(mask<0.5, yt)


        alpha = lon*pi/180.
        uvelT = yt*sin(alpha) + xt*cos(alpha)
        vvelT = yt*cos(alpha) - xt*sin(alpha) 
        
    
        # Set masked values back to nan for gridding purposes
        uvelT[where(ma.getmask(uvelT))]=np.nan
        vvelT[where(ma.getmask(vvelT))]=np.nan

        #print uvel
        # Re-grid data
        #print uvel.flatten().shape, xptsM.flatten().shape, xptsG.shape
        uvelG = griddata((xptsM.flatten(), yptsM.flatten()),uvelT.flatten(), (xptsG, yptsG), method='linear')
        vvelG = griddata((xptsM.flatten(), yptsM.flatten()),vvelT.flatten(), (xptsG, yptsG), method='linear')

        # Rotate data onto new grid
        xvelGT,yvelGT = m.rotate_vector(uvelG,vvelG,lonsG,latsG)
        xvelGT=ma.masked_invalid(xvelGT)
        yvelGT=ma.masked_invalid(yvelGT)

        xvelG[x]=xvelGT
        yvelG[x]=yvelGT

        curlG[x] = calcCurlSq2dXYGradient(xvelGT, yvelGT, dxRes)

        # print x, curlG[x]
        x+=1
        # COULD ROTATE HERE AND DO CURL OF DAILY VARIABLES.
    xvelMean=ma.mean(xvelG, axis=0)
    yvelMean=ma.mean(yvelG, axis=0)
    curlMean=ma.mean(curlG, axis=0)
        #vvelD=vstack([vvelD, vvelT])
    
    return xvelMean, yvelMean, curlMean

    #grid data onto 100km grid - matches wind forcing fields
    #Already on XY grid
    #Also put it on a grid that is the right way around (bottom to top)

    #uG, vG = interpUv(ut, vt, q, xptsF, yptsF, xptsG, yptsG)

    #mask = where((q<=0) | (q>1), 0, 1)

    #ut = ma.masked_where(mask<0.5, ut)
    #vt = ma.masked_where(mask<0.5, vt)

    #u_int=interp_data(ut, xpts, ypts, xpts2m, ypts2m)
    #v_int=interp_data(vt, xpts, ypts, xpts2m, ypts2m)

    #mask_int=interp_data_nomask(mask, xpts, ypts, xpts2m, ypts2m)
    #DO 0.5 NOT 1 AS THIS MEANS MIGHT BE INFECTED WITH OTHER DATA ON THE INTERP GRID?
    #u_int = ma.masked_where(mask_int<0.25, u_int)
    #v_int = ma.masked_where(mask_int<0.25, v_int)


    #uvelD=ma.array([]).reshape(0, lon.shape[0])
    #vvelD=ma.array([]).reshape(0, lon.shape[0])

def getGriddedKimuraFromDaily(m, files, lon, xptsG, yptsG, lonsG, latsG, xptsM, yptsM):
    
    uvelD=ma.masked_all((size(files), lon.shape[0]))
    vvelD=ma.masked_all((size(files), lon.shape[0]))
    
    #uvelD=ma.array([]).reshape(0, lon.shape[0])
    #vvelD=ma.array([]).reshape(0, lon.shape[0])
    x=0
    for file in files:
        f = loadtxt(file, unpack=True)
        # Comes in xy coordinates so need to rotate to UV
        xvel= ma.masked_where(np.isnan(f[0]), f[0])
        yvel = ma.masked_where(np.isnan(f[1]), f[1])
        #xvel=f[0]
        #yvel=f[1]

        #put on uv grid to help with regridding (consistent direction)
        alpha = lon*pi/180.
        uvelT = yvel*sin(alpha) + xvel*cos(alpha)
        vvelT = yvel*cos(alpha) - xvel*sin(alpha) 
        uvelD[x]=uvelT
        vvelD[x]=vvelT
        x+=1
        # COULD ROTATE HERE AND DO CURL OF DAILY VARIABLES.

        #vvelD=vstack([vvelD, vvelT])
    uvel=ma.mean(uvelD, axis=0)
    vvel=ma.mean(vvelD, axis=0)

    #if we want to set masked values back to nan for gridding purposes
    uvel[where(ma.getmask(uvel))]=np.nan
    vvel[where(ma.getmask(vvel))]=np.nan
    # Re-grid data
    uvelG = griddata((xptsM, yptsM),uvel, (xptsG, yptsG), method='linear')
    vvelG = griddata((xptsM, yptsM),vvel, (xptsG, yptsG), method='linear')

    # Rotate data onto new grid
    xvelG,yvelG = m.rotate_vector(uvelG,vvelG,lonsG,latsG)
    xvelG=ma.masked_invalid(xvelG)
    yvelG=ma.masked_invalid(yvelG)
    return xvelG, yvelG

def getSavedMonthlyFowlerMeanDrifts(mF, rawdatapath, startYear, endYear):
    # need to use Fowler map to ensure drifts are orientated correctly
    fowlerPath = rawdatapath+'/ICE_DRIFT/FOWLER/'
    lonlatF = loadtxt(fowlerPath+'/north_x_y_lat_lon.txt')
    lonsF = np.reshape(lonlatF[:, 3], (361, 361))
    latsF = np.reshape(lonlatF[:, 2], (361, 361))
    xptsF, yptsF = mF(lonsF, latsF)

    drift_xy = load(fowlerPath+'1980-2014-drift_data_months_xy.txt')
    drift_xy=drift_xy[startYear-1980]
    return drift_xy, latsF, lonsF, xptsF, yptsF

def getFowlerLonLat(mF, rawdatapath, startYear, endYear):
    # need to use Fowler map to ensure drifts are orientated correctly
    fowlerPath = rawdatapath+'/ICE_DRIFT/FOWLER/'
    lonlatF = loadtxt(fowlerPath+'/north_x_y_lat_lon.txt')
    lonsF = np.reshape(lonlatF[:, 3], (361, 361))
    latsF = np.reshape(lonlatF[:, 2], (361, 361))
    xptsF, yptsF = mF(lonsF, latsF)

    return latsF, lonsF, xptsF, yptsF

def getMonthlyFowlerDailyDrifts(rawdatapath, year, month):
    fowlerPath = rawdatapath+'/ICE_DRIFT/FOWLER/DAILY/'
    timeIndex = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
    #if (year>2012):
        #timeIndex[12]=364
    #    files = glob(fowlerPath+str(year)+'/*.vec')
    #else:
    files = glob(fowlerPath+str(year)+'/*.bin')

    return files[timeIndex[month]:timeIndex[month+1]]

def calcCurlSq2dXYGradient(x_vel, y_vel, dx_res):
# calculate the curl (squared) of a field given some input vector grid
    #CALCULATE THE CURL OF AN X/Y VEXTOR FIELD. DX_RES IF THE GRID RESOLUTION OF THIS REGULARLY SPACED GRID.
    #MULTIPLY BY MAGNITUDE TO GET THIS IN SQUARED UNITS ANALAGOUS TO THE WIND STRESS CURL.
    #USE GRADIENT FUNCTION WHICH USES CENTRAL DIFFERENCES IN MIDDLE CELLS AND FIRST DIFFERENCES AT THE BOUNDARIES (GIVES SAME SHAPE AS INPUT)

    mag = sqrt((x_vel**2) + (y_vel**2))

    x_vel_mag = x_vel*mag
    y_vel_mag = y_vel*mag

    #gradient [0] returns row divs (y direction) then [1] gives column divs (x direction)
    dvelydx = np.gradient(y_vel_mag, dx_res)[1]
    dvelxdy = np.gradient(x_vel_mag, dx_res)[0]

    zeta = dvelydx - dvelxdy
    #MASK ARRAY WHERE VALUES ARE NAN
    zeta = ma.array(zeta,mask=np.isnan(zeta))

    return zeta

def interpUv(ut, vt, q, xpts, ypts, xpts2m, ypts2m):

    #TURN UT/VT INTO MASKED ARRAYS
    #IF Q NEGATIVE THEN REMOVE AS IT IS BY COAST
    #IF Q GREATER THAN 0.5 THEN REMOVE AS THIS MEANS NO POINTS NEARBY IN MAKING THE VECTOR
    #IF 0 THEN MASK AS NO DATA BUT THINK THAT"S ALREADY THE CASE

    #negative q is coastal, 0 is no data, greater than 1 means no data closeby
    #create mask where 0 is mask and 1 is data
    mask = where((q<=0) | (q>1), 0, 1)

    ut = ma.masked_where(mask<0.5, ut)
    vt = ma.masked_where(mask<0.5, vt)

    u_int=interp_data(ut, xpts, ypts, xpts2m, ypts2m)
    v_int=interp_data(vt, xpts, ypts, xpts2m, ypts2m)

    mask_int=interp_data_nomask(mask, xpts, ypts, xpts2m, ypts2m)
    #DO 0.5 NOT 1 AS THIS MEANS MIGHT BE INFECTED WITH OTHER DATA ON THE INTERP GRID?
    u_int = ma.masked_where(mask_int<0.25, u_int)
    v_int = ma.masked_where(mask_int<0.25, v_int)


    return u_int, v_int

def plot_var_xy(m, xpts , ypts, var_x, var_y, var_mag, out='./figure', units_lab='units', units_vec=r'm s$^{-1}$',
 minval=1., maxval=1., base_mask=1,res=1, scale_vec=1, vector_val=1, year_string='year', month_string='months', extra='',cbar_type='both', cmap_1=plt.cm.RdBu_r):

        #PLOT SCALAR FIELD WITH OVERLYING VECTORS. 
        #VAR MAG MAY NOT NECCESARRILY BE THE MAGNITUDE OF THE VECTORS (E.G. IN THE CASE OF WIND CURL)

        fig = figure(figsize=(3.5,4))
        ax1 = fig.add_axes([0.0, 0.15, 1.0, 0.85])
        #if etopo==1:
         #       im_etopo = m.pcolormesh(xpts_etopo, ypts_etopo , etopo_var, cmap=plt.cm.Greens_r, vmin=0, vmax=1000, zorder=1)
        if (maxval-minval)<1e-10:
            minval = -round_to_1(np.amax(var_mag))
            maxval = round_to_1(np.amax(var_mag))

        #var_mag=ma.masked_where(var_mag<1e8, var_mag)
        im1 = m.pcolormesh(xpts , ypts, var_mag, cmap=cmap_1,vmin=minval, vmax=maxval,shading='flat', zorder=4)
        # LOWER THE SCALE THE LARGER THE ARROW
        Q = m.quiver(xpts[::res, ::res], ypts[::res, ::res], var_x[::res, ::res], var_y[::res, ::res], units='inches',scale=scale_vec, zorder=5)
        #ASSIGN A LEGEND OF THE VECTOR 
        #m.plot(xpts[191, 100], ypts[191, 100], 'x', zorder=10)
        m.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
        m.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)
        #m.drawmapboundary(fill_color='0.3')
        #m.drawmapboundary(fill_color='0.4' , zorder=1)
        if base_mask==1:
        #m.drawmapboundary(fill_color='0.4' , zorder=1)
            m.fillcontinents(color='0.7',lake_color='grey', zorder=6)
        m.drawcoastlines(linewidth=0.5)

        cax = fig.add_axes([0.1, 0.1, 0.8, 0.04])
        cbar = colorbar(im1,cax=cax, orientation='horizontal', extend=cbar_type, use_gridspec=True)
        cbar.set_label(units_lab)
        xticks = np.linspace(minval, maxval, 5)
        cbar.set_ticks(xticks)
        cbar.formatter.set_powerlimits((-3, 4))
        cbar.formatter.set_scientific(True)
        cbar.update_ticks() 

        xS, yS = m(140, 58)
        qk = quiverkey(Q, xS, yS, vector_val, str(vector_val)+' '+units_vec, fontproperties={'size': 'medium'}, coordinates='data', zorder = 11)   
        
        xS, yS = m(235, 55)
        ax1.text(xS, yS, year_string+'\n'+month_string+'\n'+extra,fontsize=10, zorder = 11)

        subplots_adjust(bottom=0.0, left=0.0, top = 1.0, right=1.0)

        savefig(out+'.png', dpi=300)
        close(fig)
#plot vector map (with vectors in x/y directions)


def interp_data(var, xpts, ypts, xpts2m, ypts2m, int_meth='cubic'):
#interpoalte data onto a regular 2d grid. Used by interp_uv
    data = ma.getdata(var)
    mask = ma.getmask(var)
    index = np.where(mask==False)
    data_masked = data[index]
    xpoints = xpts[index]
    ypoints = ypts[index]
    points = (xpoints, ypoints)
    var_int = griddata(points,data_masked,(xpts2m,ypts2m),method=int_meth)
    var_int = ma.masked_array(var_int,np.isnan(var_int))

    return var_int

def interp_data_nomask(var, xpts, ypts, xpts2m, ypts2m, int_meth='cubic'):
#interpoalte mask of data onto a regular 2d grid. Used by interp_uv
    index = np.where(var>-1)
    xpoints = xpts[index]
    ypoints = ypts[index]
    points = (xpoints, ypoints)
    data_ma = var[index]
    #tck = interpolate.bisplrep(points, data_masked, , s=0)
    #var_int = interpolate.bisplev(xnew[:,0], ynew[0,:], tck)
    var_int = griddata(points,data_ma,(xpts2m,ypts2m),method=int_meth)
    return var_int
