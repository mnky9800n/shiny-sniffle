#!/bin/python

import pandas as pd
import numpy as np
import math
import myCalc as calc

def geo_project(df,originLon,originLat,bearing,dist,maxDist=0):

    # projects longitude and latitude of events in df 
    # into a cartesian coordinate system based on 
    # origin longitude (originLon), origin latitude 
    # (originLat), and bearing from the origin point
    # at a given distance (dist)
    #
    # input:
    # df = dataframe of events to be projected, must have
    #      at least 'la' and 'lo' columns
    # bearing = degrees from north
    # originLon = origin longitude
    # originLat = origin latitude
    # dist = distance from origin for transect end point calculation
    
    def transect_endPts(lo,la,dist,bearing):
    
    # epicentral distance in km from lo, la
        # bearing in degrees from north
        # returns: lons and lats of two points of a transect segment

        # check bearing
        if bearing > 180:
            bearing = bearing-180

        # assess the extent of the possible longitudes and latitudes
        deg2rad = math.pi/180.0
        bearingRot = (90 - bearing)*deg2rad
        dx = 1.5*dist*math.cos(bearingRot)
        dy = 1.5*dist*math.sin(bearingRot)

        km2deg = 1.0/111.19
        dLon = math.fabs(dx*km2deg)
        dLat = math.fabs(dy*km2deg)

        lonMin = lo-dLon
        lonMax = lo+dLon
        latMin = la-dLat
        latMax = la+dLat

        loTest = np.arange(lonMin,lonMax,0.01)    
        laTest = np.arange(latMin,latMax,0.01)

        # construct a dataframe for computing the distance and back-azimuth
        # to each possible point
        m = len(loTest)
        n = len(laTest)

        cols = ['lon1','lat1','lon2','lat2','dist','baz','bazDiffabs']
        df = pd.DataFrame(columns=cols,index=np.arange(0,m*n,1))    
        df.lon1 = pd.Series(lo*np.ones(m*n),index=df.index)
        df.lat1 = pd.Series(la*np.ones(m*n),index=df.index)

        lon2 = []
        lat2 = []
        for i in np.arange(0,m,1):
            for j in np.arange(0,n,1):
                lon2.append(loTest[i])
                lat2.append(laTest[j])

        df.lon2 = pd.Series(lon2,index=df.index)
        df.lat2 = pd.Series(lat2,index=df.index)

        # compute the distances and back-azimuths
        out = calc.ndistbaz(df.lon1,df.lat1,df.lon2,df.lat2)

        df.dist = pd.Series(out.distance.values,index=df.index)
        df.baz = pd.Series(out.baz.values,index=df.index)

        # find the first point which corresponds to the provided bearing
        bazDiff = [math.fabs(bearing - i) for i in df.baz]
        df.bazDiffabs = pd.Series(bazDiff,index=df.index)

        df2 = df.copy()
        df2 = df2[df2.bazDiffabs==df2.bazDiffabs.min()]
        df2 = df2.reset_index(drop=True)

        point1 = [df2.lon2[0],df2.lat2[0]]

        # find the second point which corresponds to the provided bearing +/- 180 degrees
        if bearing > 180:
            bearingOpposite = bearing-180
        elif bearing < 180:
            bearingOpposite = bearing + 180
        elif bearing == 180:
            bearingOpposite = 0

        bazDiff = [math.fabs(bearingOpposite - i) for i in df.baz]
        df.bazDiffabs = pd.Series(bazDiff,index=df.index)

        df2 = df.copy()
        df2 = df2[df2.bazDiffabs==df2.bazDiffabs.min()]
        df2 = df2.reset_index(drop=True)

        point2 = [df2.lon2[0],df2.lat2[0]]

        return (point1,point2)
        
    def reduce_dataframe(df,point1,point2):

        # data frame must have header 'lo' and 'la'

        p = point1
        q = point2
        s = df.copy()

        if p[0] < q[0]:
            pIsWest = 1
        else:
            pIsWest = 0

        if p[1] < q[1]:
            pIsSouth = 1
        else:
            pIsSouth = 0

        # check that the point is within the bounds of the transect
        if pIsWest == 1 and pIsSouth == 1:
            s = s[(s.lo>p[0]) & (s.lo<q[0]) & (s.la>p[1]) & (s.la<q[1])]
        elif pIsWest == 0 and pIsSouth == 1:
            s = s[(s.lo>q[0]) & (s.lo<p[0]) & (s.la>p[1]) & (s.la<q[1])]
        elif pIsWest == 1 and pIsSouth == 0:
            s = s[(s.lo>p[0]) & (s.lo<q[0]) & (s.la>q[1]) & (s.la<p[1])]
        elif pIsWest == 0 and pIsSouth == 0:
            s = s[(s.lo>q[0]) & (s.lo<p[0]) & (s.la>q[1]) & (s.la<p[1])]

        s = s.reset_index(drop=True)

        return s
    
    def singlePointDistanceRelative2Transect(p,q,r,s):

        # p, q, r, s must be lists with
        # longitude is first entry and latitude is second entry of each element
        #
        # p and q cannot not have the same latitude or same longitude 
        # (no horizontal/vertical transects)
        #
        # p = start point of transect
        # q = end point of transect
        # r = origin of transect
        # s = point of interest

        # p, q, r - must all line on the same great circle path 
        # i.e., q and r must have same back-azimuth relative to p

        # check the orientation of the transect
        # compute the backazimuth
        d, b = calc.distbaz(p[0],p[1],q[0],q[1])
        if b > 180:
            a = p
            b = q
            p = b
            q = a

        if p[0] < q[0]:
            pIsWest = 1
        else:
            pIsWest = 0

        if p[1] < q[1]:
            pIsSouth = 1
        else:
            pIsSouth = 0

        # check that the point is within the bounds of the transect
        if pIsWest == 1 and pIsSouth == 1:
            if (s[0]>=p[0]) and (s[0]<=q[0]) and (s[1]>=p[1]) and (s[1]<=q[1]):
                inBounds = 1
            else:
                inBounds = 0
        elif pIsWest == 0 and pIsSouth == 1:
            if (s[0]>=q[0]) and (s[0]<=p[0]) and (s[1]>=p[1]) and (s[1]<=q[1]):
                inBounds = 1
            else:
                inBounds = 0
        elif pIsWest == 1 and pIsSouth == 0:
            if (s[0]>=p[0]) and (s[0]<=q[0]) and (s[1]>=q[1]) and (s[1]<=p[1]):
                inBounds = 1
            else:
                inBounds = 0
        elif pIsWest == 0 and pIsSouth == 0:
            if (s[0]>=q[0]) and (s[0]<=p[0]) and (s[1]>=q[1]) and (s[1]<=p[1]):
                inBounds = 1
            else:
                inBounds = 0

        if inBounds == 1:

            if p == s:
                xDistanceComponent = 0
                yDistanceComponent = 0
            else:
                # compute distance in kilometers -> from P to Q (transect)
                # compute the back-azimuth angle -> from P to Q (transect)
                transectDistance, transectBackAzimuth = calc.distbaz(p[1],p[0],q[1],q[0])

                # compute distance in kilometers -> from P to A (transect start to event)
                # compute the back-azimuth angle -> from P to A (transect start to event)
                eventDistance, eventBackAzimuth = calc.distbaz(p[1],p[0],s[1],s[0])

                # compute the angle the event makes with the transect
                # convert to radians
                deg2rad = math.pi/180.0
                angleWithTransect = (transectBackAzimuth - eventBackAzimuth)*deg2rad

                # compute the distance along the transect to the point
                # take the absolute value (disregard sign for now)
                xDistanceComponentMagnitude = math.fabs(eventDistance*math.cos(math.fabs(angleWithTransect)))

                # compute the distance away from the transect to the point
                # take the absolute value (disregard sign for now)
                yDistanceComponentMagnitude = math.fabs(eventDistance*math.sin(math.fabs(angleWithTransect)))

                # determine the sign of the yDistanceComponent
                if yDistanceComponentMagnitude == 0:
                    yDistanceComponent = 0
                elif angleWithTransect < 0:
                    yDistanceComponent = yDistanceComponentMagnitude*(-1)
                elif angleWithTransect > 0:
                    yDistanceComponent = yDistanceComponentMagnitude

                # shift the xDistanceComponent relative to the provided origin, r
                if p == r:
                    xDistanceComponent = xDistanceComponentMagnitude
                else:
                    # compute the distance in kilometers -> p to r
                    originDistance, originBackAzimuth = calc.distbaz(p[1],p[0],r[1],r[0])

                    if int(originBackAzimuth) != int(transectBackAzimuth):
                        print('Error: The transect does not run through origin!')
                        return ([NaN,NaN])
                    else:
                        xDistanceComponent = xDistanceComponentMagnitude - originDistance
                        return ([xDistanceComponent,yDistanceComponent])
        else:
            print('Error: The point of interest is not in the bounds of the transect!')
            return ([NaN,NaN])

    def multiPointDistanceRelative2Transect(p,q,r,s):

        # p, q, r must be lists of length 1 with
        # longitude is first entry and latitude is second entry of each element
        #
        # p and q cannot not have the same latitude or same longitude 
        # (no horizontal/vertical transects)
        #
        # p = start point of transect
        # q = end point of transect
        # r = origin of transect
        # s = dataframe of events with at least 'lo' and 'la' columns

        # p, q, r - must all line on the same great circle path 
        # i.e., q and r must have same back-azimuth relative to p

        # returns a list with the xx, yy distances

        # check the orientation of the transect (make it <= 180)
        # compute the backazimuth
        d, b = calc.distbaz(p[0],p[1],q[0],q[1])
        if b > 180:
            a = p
            b = q
            p = b
            q = a

        s2 = reduce_dataframe(s,p,q)

        XY = []
        for i,j in zip(s2.lo,s2.la):

            point = [i,j]
            xy = singlePointDistanceRelative2Transect(p,q,r,point)

            XY.append(xy)

        return(XY)
        
    def reduce_dataframe_radially(df,r,dist):

        # reduces the dataframe to keep only events within
        # the provided radial distance from point r
        # r is a list with longitude, latitude
        # df is a dataframe with at least lo, la
        # dist is distance in kilometer from point r

        cols = ['lat1','lat2','lon1','lon2','dist']
        df2 = pd.DataFrame(columns=cols,index=df.index)
        df2.lat1 = pd.Series(r[1]*np.ones(len(df2)),index=df2.index)
        df2.lon1 = pd.Series(r[0]*np.ones(len(df2)),index=df2.index)
        df2.lat2 = pd.Series(df.la,index=df2.index)
        df2.lon2 = pd.Series(df.lo,index=df2.index)

        out = calc.ndistbaz(df2.lon1,df2.lat1,df2.lon2,df2.lat2)

        df2.dist = pd.Series(out.distance.values,index=df2.index)

        df2 = df2[df2.dist<=dist]

        df3 = df.copy()
        df3 = df3.iloc[df2.index]

        return df3

    # find the transect endpoints
    lo = originLon
    la = originLat
    pt1, pt2 = transect_endPts(lo,la,dist,bearing)

    # compute the xx and yy distances to the origin along the bearing/transect
    XY = multiPointDistanceRelative2Transect(pt1,pt2,[lo,la],df)    
    
    # reduce the original data frame to match the xx, yy output
    dfSlice = reduce_dataframe(df,pt1,pt2)
    
    # add the xx, yy values to the dataframe
    xx = [XY[i][0] for i in np.arange(0,len(XY),1)]
    dfSlice['xx'] = pd.Series(xx,index=dfSlice.index)

    yy = [XY[i][1] for i in np.arange(0,len(XY),1)]
    dfSlice['yy'] = pd.Series(yy,index=dfSlice.index)
    
    # reduce events to only those within a specified radially distance
    if maxDist != 0:
        dfOut = reduce_dataframe_radially(dfSlice,[lo,la],maxDist)
        dfOut = dfOut.reset_index(drop=True)
    else:
        dfOut = dfSlice.copy()

    return(dfOut)

def transect_endPts(lo,la,dist,bearing):
    
    # epicentral distance in km from lo, la
    # bearing in degrees from north
    # returns: lons and lats of two points of a transect segment

    # check bearing
    if bearing > 180:
        bearing = bearing-180

    # assess the extent of the possible longitudes and latitudes
    deg2rad = math.pi/180.0
    bearingRot = (90 - bearing)*deg2rad
    dx = 1.5*dist*math.cos(bearingRot)
    dy = 1.5*dist*math.sin(bearingRot)

    km2deg = 1.0/111.19
    dLon = math.fabs(dx*km2deg)
    dLat = math.fabs(dy*km2deg)

    lonMin = lo-dLon
    lonMax = lo+dLon
    latMin = la-dLat
    latMax = la+dLat

    loTest = np.arange(lonMin,lonMax,0.01)    
    laTest = np.arange(latMin,latMax,0.01)

    # construct a dataframe for computing the distance and back-azimuth
    # to each possible point
    m = len(loTest)
    n = len(laTest)

    cols = ['lon1','lat1','lon2','lat2','dist','baz','bazDiffabs']
    df = pd.DataFrame(columns=cols,index=np.arange(0,m*n,1))    
    df.lon1 = pd.Series(lo*np.ones(m*n),index=df.index)
    df.lat1 = pd.Series(la*np.ones(m*n),index=df.index)

    lon2 = []
    lat2 = []
    for i in np.arange(0,m,1):
        for j in np.arange(0,n,1):
            lon2.append(loTest[i])
            lat2.append(laTest[j])

    df.lon2 = pd.Series(lon2,index=df.index)
    df.lat2 = pd.Series(lat2,index=df.index)

    # compute the distances and back-azimuths
    out = calc.ndistbaz(df.lon1,df.lat1,df.lon2,df.lat2)

    df.dist = pd.Series(out.distance.values,index=df.index)
    df.baz = pd.Series(out.baz.values,index=df.index)

    # find the first point which corresponds to the provided bearing
    bazDiff = [math.fabs(bearing - i) for i in df.baz]
    df.bazDiffabs = pd.Series(bazDiff,index=df.index)

    df2 = df.copy()
    df2 = df2[df2.bazDiffabs==df2.bazDiffabs.min()]
    df2 = df2.reset_index(drop=True)

    point1 = [df2.lon2[0],df2.lat2[0]]

    # find the second point which corresponds to the provided bearing +/- 180 degrees
    if bearing > 180:
        bearingOpposite = bearing-180
    elif bearing < 180:
        bearingOpposite = bearing + 180
    elif bearing == 180:
        bearingOpposite = 0

    bazDiff = [math.fabs(bearingOpposite - i) for i in df.baz]
    df.bazDiffabs = pd.Series(bazDiff,index=df.index)

    df2 = df.copy()
    df2 = df2[df2.bazDiffabs==df2.bazDiffabs.min()]
    df2 = df2.reset_index(drop=True)

    point2 = [df2.lon2[0],df2.lat2[0]]

    return (point1,point2)
    
