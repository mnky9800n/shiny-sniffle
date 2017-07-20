#!/bin/python

from math import pi, cos, sin, sqrt, asin, acos, log10
import math
import pandas as pd
import numpy as np
from pytz import timezone
from datetime import datetime
import time
import os
import calendar
from decimal import Decimal

def distbaz(la1,lo1,la2,lo2):
    '''
        (1) Computes distance between two spherical coordinates
        (2) Computes back-azimuth two spherical coordinates

        la1 = station latitude
        lo1 = station longitude
        la2 = event latitude
        lo2 = event longitude

        returns distance and back-azimuth
    '''

    d2r = pi/180.0
    r2d = 180.0/pi
    rho = 6372.795
    t1 = (90.0-la1)*d2r
    t2 = (90.0-la2)*d2r

    # avoid math domain error
    if (lo1 == lo2): 
        p1 = (lo1+0.0001)*d2r
    else:
        p1 = lo1*d2r

    p2 = lo2*d2r

    c1, s1 = cos(t1), sin(t1)
    c2, s2 = cos(t2), sin(t2)
    arg1 = (c1*c2)+s2*s1*cos(p1-p2)
    delta = acos(arg1)
    distance = delta*rho
    arg2 = (1/sin(delta))*((c2*s1)-(s2*c1*cos(p2-p1)))
    
    if arg2 < -1 or arg2 > 1:
        arg2 = int(arg2)

        
    baz1 = acos(arg2)*r2d
    arg3 = (1/sin(delta))*(s2*sin(p2-p1))
    if arg3 < -1 or arg3 > 1:
        arg3 = int(arg3)

    baz2 = asin(arg3)*r2d
    if baz2 < 0:
        baz = 360 - baz1    
    else:
        baz = baz1
    
#    baz1 = acos(arg2)*r2d
#    arg3 = (1/sin(delta))*(s2*sin(p2-p1))
#    baz2 = asin(arg3)*r2d
#    if baz2 < 0:
#        baz = 360 - baz1    
#    else:
#        baz = baz1

    return (distance, baz)

def ndistbaz(stlo,stla,evlo,evla):
    '''
        Computes distance and back-azimuth for multiple events.
		Adds this information to a dataframe
    '''

    # INPUT
    # stlo, stla, evlo, evla are panda Series of equal length
    # stlo = station longitude
    # stla = station latitude
    # evlo = event longitude
    # evla = event latitude

    # OUTPUT
    # panda Dataframe of distance and back-azimuth from station to event

    # begin output
    out = pd.DataFrame({'stlo':stlo,'stla':stla,'evlo':evlo,'evla':evla})

    # initialize
    d2r = math.pi/180.0
    r2d = 180.0/math.pi
    rho = 6372.795

    # avoid math domain error
    loDifference = stlo - evlo
    indices = loDifference[loDifference==0].index.tolist()
    if indices:
        stlo.ix[indices] = stlo.ix[indices] + 0.000001

    laDifference = stla - evla
    indices = laDifference[laDifference==0].index.tolist()
    if indices:
        stla.ix[indices] = stla.ix[indices] + 0.000001

    # compute sines and cosines
    theta1 = d2r*(90*np.ones(len(stla)) - stla) # orientation matters
    theta2 = d2r*(90*np.ones(len(evla)) - evla)
    phi1 = stlo*d2r
    phi2 = evlo*d2r

    cos1, sin1 = np.cos(theta1), np.sin(theta1)
    cos2, sin2 = np.cos(theta2), np.sin(theta2)

    # compute part one of the equation (Stein and Wysession, Introduction to Seismology)
    arg1 = (cos1*cos2) + sin2*sin1*np.cos(phi1-phi2)
    delta = np.arccos(arg1)

    # compute the distance between the station and event
    distance = delta*rho*np.ones(len(delta))

    # compute part two of the equation
    arg2 = (1/np.sin(delta))*((cos2*sin1)-(sin2*cos1*np.cos(phi2-phi1)))

    # avoid math domain error
    indices = arg2[(arg2 < -1) | (arg2 > 1)].index.tolist()
    if indices:
        arg2[arg2 < -1] = -1.0
        arg2[arg2 > 1 ] = 1.0

    # compute back-azimuth
    baz = np.arccos(arg2)*r2d

    # compute part three of the equation
    arg3 = (1/np.sin(delta))*(sin2*np.sin(phi2-phi1))

    # avoid math domain error
    indices = arg3[(arg3 < -1) | (arg3 > 1)].index.tolist()
    if indices:
        arg3[arg3 < -1] = -1.0
        arg3[arg3 > 1 ] = 1.0

    # check the sign of the back-azimuth
    baz2 = np.arcsin(arg3)*r2d
    indices = baz2[(baz2<0)].index.tolist()
    if indices:
        baz[baz2<0] = 360 - baz[baz2<0]

    # dataframe out
    out['distance'] = pd.Series(distance)
    out['baz'] = pd.Series(baz)

    # if the station location is the same as an event set it to zero
    row_ids = out[(out['stlo'] == out['evlo']) & (out['stla'] == out['evla'])].index.tolist()
    if row_ids:
        out.distance.ix[row_ids] = 0.0
        out.baz.ix[row_ids] = 0.0

    return out

def near_stress(f,g):
    
    # compute PGV for events within 800 km
    # m = magnitude
    # r = hypocentral distance in km
    # returns dynamic stress in kPa    
    #
    # basic equation stress = shear rigidity * PGV / phase velocity
    # where shear rigidity is 35 GPa and phase velocity is 3.5 km/s
    # log10(PGV) = c1 + c2*M - c3*log10(sqrt(r^2+c4^2))
    # where c1 = -2.29, c2 = 0.85, c3 = 1.29, c4 = 0
    
    P = math.pow(10,(-2.29+(0.85*f)-(1.29*log10(g))))
#   P = math.pow(10,-2.29)*math.pow(10,0.85*f)*math.pow(10,-1.29*math.log10(g))
    sigma = P*100
    return sigma 

def far_stress(p,q):

    # compute PGV for events beyond 800 km
    # p = magnitude
    # q = epicentral distance in degrees
    # returns dynamic stress in kPa
    #
    # basic equation stress = shear rigidity * PGV / phase velocity
    # where shear rigidity is 35 GPa and phase velocity is 3.5 km/s
    # PGV is computed as 2 * pi * A /20  where A is the amplitude of the 20-s period surface wave
    # computed empirical as below.
    # see van der Elst and Brodksy (2010)

    A = math.pow(10,(p-(1.66*log10(q))-6)) # convert to meters
#    A = math.pow(10,p)*math.pow(10,-1.66*math.log10(q))*math.pow(10,-6)
    sigma = A*pi*10 # simplified
    return sigma	

def decimal_split(df): 
    # splits one column data frame (i.e. a series) of decimal values
    # into a two column data frame
    # df is series
    #
    # EXAMPLE USAGE:
    # S = pd.DataFrame({'sc' : df3['sc'][:]})
    # integers, decimals = calc.decimal_split(S)
    #
    df.columns = ['a']
    dfSplit = df.a.apply(math.modf)
    num = len(dfSplit)
    integers = []
    decimals = []
    for p in np.arange(0,num,1):
        decimal, integer = dfSplit[p]
        decimals.append(decimal)
        integers.append(integer)
    return (integers, decimals)

def roundup(x,i):
    '''
    x = number to be rounded
    i = increment to be round to
    
    EXAMPLE:
    >>> x = 130
    >>> i = 100
    >>> out = roundup(x,i)
    >>> print(out)
    200
    
    '''
    return x if x % i == 0 else x + i - x % i

def fmd(df,binSize):
    ''' 
    computes Mc value of a catalog 

    Mc = magnitude with maximum value
    of the first derivative of the freq-mag 
    distribution (may or may not be bin with
    the highest frequency of events).

    input is a data frame with at least magnitudes

    '''
    
    # remove bad magnitude events
    df2 = df.copy()
    df2 = df2[df2.mg<=10]

    df2 = df2.replace([np.inf,-np.inf],np.nan).dropna()

    # compute the frequencies for each magnitude
    # get mid-point of bins
    minMag = df2.mg.min()
    maxMag = df2.mg.max()

    bins = np.arange(minMag,maxMag+2*binSize,binSize)
    data = df2.mg

    y, x = np.histogram(data,bins)

    midBin = []
    for x in np.arange(1,len(bins),1):
        mid = 0.5*(bins[x]+bins[x-1])
        midBin.append(mid)

    # identify Mc value from greatest derivative
    indices = np.arange(0,len(midBin),1)
    fmd = pd.DataFrame(columns=['m','N'],index=indices)
    fmd.m = pd.Series(midBin,index=indices)
    fmd.N = pd.Series(y,index=indices)

    fmdDIFF = fmd.diff()
    Mc = fmd.m[fmdDIFF.N.argmax()]

    # drop rows with a N-value of 0
    fmd = fmd[fmd.N!=0]
  
    # compute the cumulative sum of events
    fmd = fmd.sort(['m'],ascending=[0])
    fmd = fmd.reset_index(drop=True)
    fmd['cumSum'] = fmd.N.cumsum()   

    return fmd, Mc

def convert_time(df,q,r):
    # 
    # converts UTC times for formatted catalogs
    # df = data frame with at least 'ep' column
    # q = time in modified catalog (e.g. 'UTC')
    # r = time to convert to (e.g. 'Japan')
    #
    # returns data frame 'out' with new date/time/epoch
    #
    # for a list of timezones, try the following:
    #
    # import pytz
    # for tz in pytz.all_timezones:
    #     print tz

    current_time = datetime.now(timezone(q))
    new_time = current_time.astimezone(timezone(r))    

    fmt = '%Y %m %d %H %M %S %Z%z'
    tFormat = new_time.strftime(fmt)

    tSplit = tFormat.split('+')
    tDiff = int(float(tSplit[1])/100.0)*3600

    new_epochs = [df.ep[x]+tDiff for x in np.arange(0,len(df),1)]

    pattern = '%Y %m %d %H %M %S'
    new_date_time = [ time.strftime(pattern,time.gmtime(new_epochs[x])) \
                      for x in np.arange(0,len(df),1)] 

    date_split = [ new_date_time[x].split(' ') \
                   for x in np.arange(0,len(df),1) ]

    out = df.copy()
    out.yr = [date_split[x][0] for x in np.arange(0,len(df),1)]
    out.mo = [date_split[x][1] for x in np.arange(0,len(df),1)]
    out.dy = [date_split[x][2] for x in np.arange(0,len(df),1)]
    out.hr = [date_split[x][3] for x in np.arange(0,len(df),1)]
    out.ep = [new_epochs[x] for x in np.arange(0,len(df),1)] 
    
    return out

def norm_min_max(A):
    # normalizes a data set A
    # by the absolute maximum

    if type(A) in ['pandas.core.series.Series']:
        minA = A.min()
        maxA = A.max()

        if minA > maxA:
            dA = [z/minA for z in A]
        else:
            dA = [z/maxA for z in A]
    else:
        minA = min(A)*-1
        maxA = max(A)

        if minA > maxA:
            dA = A/minA
        else:
            dA = A/maxA

    return dA

def MAD(s):
    # computes the median absolute deviation of the dataset (s)

    med = np.median(s)
    dev = s - med*np.ones_like(s)
    aDev = np.abs(dev)
    MAD = np.median(aDev)

    return (MAD)

def projection(origin_lon,origin_lat,bearing,df):

    # this function transforms longitude and latitude into
    # cartesian coordinates and then projects them along
    # the bearing provided relative to an origin point.

    # input:
    # df = a dataframe with at least 'la','lo',and 'dp' columns
    # origin_lon = origin longitude about which to rotate the coordinates 
    # origin_lat = origin latitude about which to rotate the coordinates
    # bearing = direction from north to rotate the coordinates

    # rename the variables
    oLo = origin_lon
    oLa = origin_lat
    bearingRot = 90 - bearing # rotate so that 0 deg is east (not north)

    # create dataframe for computation
    cols = ['oLo','oLa','oX','oY','eLo','eLa','eX','eY','eXshift','eYshift',
            'eXrot','eYrot','xOffset','yOffset','xx','yy']
    dfOut = pd.DataFrame(columns=cols,index=df.index)

        # add column of origin lon/lat
    dfOut.oLo = pd.Series(oLo*np.ones(len(df)),index=df.index)
    dfOut.oLa = pd.Series(oLa*np.ones(len(df)),index=df.index)
        # add column of event lons/lats
    dfOut.eLo = pd.Series([df.lo[i] for i in df.index],index=df.index)
    dfOut.eLa = pd.Series([df.la[i] for i in df.index],index=df.index)

    # convert to cartesian coordinates
    rho = 6372.795
    deg2rad = math.pi/180.0

        # add column of event cartesian coordinates    
    x = [rho*math.cos(dfOut.eLo[i]*deg2rad)*math.sin((dfOut.eLa[i])*deg2rad) for i in df.index]
    y = [rho*math.sin(dfOut.eLo[i]*deg2rad)*math.sin((dfOut.eLa[i])*deg2rad) for i in df.index]
    dfOut.eX = pd.Series(x,index=df.index)
    dfOut.eY = pd.Series(y,index=df.index)
        # add column of origin cartesian coordinates
    x = [rho*math.cos(dfOut.oLo[i]*deg2rad)*math.sin((dfOut.oLa[i])*deg2rad) for i in df.index]
    y = [rho*math.sin(dfOut.oLo[i]*deg2rad)*math.sin((dfOut.oLa[i])*deg2rad) for i in df.index]
    dfOut.oX = pd.Series(x,index=df.index)
    dfOut.oY = pd.Series(y,index=df.index)

    # translate/shift the coordinate system so that the origin point is (0,0)
        # x translation
    xShift = [dfOut.eX[i]-dfOut.oX[i] for i in df.index]
    dfOut.eXshift = pd.Series(xShift,index=df.index)
        # y translation
    yShift = [dfOut.eY[i]-dfOut.oY[i] for i in df.index]
    dfOut.eYshift = pd.Series(yShift,index=df.index)

    # rotate the coordinates according to the provided bearing
        # x rotation
    eXrot = [dfOut.eXshift[i]*math.cos(bearingRot*deg2rad)-dfOut.eYshift[i]*math.sin(bearingRot*deg2rad) 
             for i in df.index]
    dfOut.eXrot = pd.Series(eXrot,index=df.index)
        # y rotation
    eYrot = [dfOut.eXshift[i]*math.sin(bearingRot*deg2rad)+dfOut.eYshift[i]*math.cos(bearingRot*deg2rad) 
             for i in df.index]
    dfOut.eYrot = pd.Series(eYrot,index=df.index)

    # determine the offset (distance of rotated coordinate from where the coordinate should be)
        # x offset
    xOffset = [dfOut.eXrot[i]-dfOut.eXshift[i] for i in df.index]
    dfOut.xOffset = pd.Series(xOffset,index=df.index)
        # y offset
    yOffset = [dfOut.eYrot[i]-dfOut.eYshift[i] for i in df.index]
    dfOut.yOffset = pd.Series(yOffset,index=df.index)

    #  compute the location in the new reference frame (origin_lon,original_lat,bearing from north)
        # x final location
    xFinal = [dfOut.eXrot[i]-dfOut.xOffset[i] for i in df.index]
    dfOut.xx = pd.Series(xFinal,index=df.index)
        # y final location
    yFinal = [dfOut.eYrot[i]-dfOut.yOffset[i] for i in df.index]
    dfOut.yy = pd.Series(yFinal,index=df.index)

    # add the translated/rotated cartesian coordinates to the input dataframe
    df['xx'] = pd.Series([i for i in dfOut.xx],index=df.index)
    df['yy'] = pd.Series([i for i in dfOut.yy],index=df.index)

    return df

