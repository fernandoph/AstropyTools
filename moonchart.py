# -*- coding: utf-8 -*-
import sys
import astropy.coordinates
import csv
from astropy.coordinates import AltAz
from astropy.time import Time, TimezoneInfo
from astropy.time import TimeDelta
import astropy.units as u
from astropy.coordinates import EarthLocation
import astroplan
from astropy.coordinates import solar_system_ephemeris
solar_system_ephemeris.set('jpl') 
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import datetime

# Posicion del observador (IAR)
lat = -34.8601
lon = -58.1385

h   = 20
IAR = EarthLocation(lat=lat, lon=lon, height=h*u.m)

observer = astroplan.Observer(location = IAR, timezone='America/Buenos_Aires')
localtz  =  TimezoneInfo(utc_offset=-3 * u.hour)

# Limites de la antena
EASTLIMIT  = -29
WESTLIMIT  = +29
NORTHLIMIT = -10
SOUTHLIMIT = -90


# Algunos TimeDeltas
# 1 Minuto
oneMin  = TimeDelta(val = 60, format='sec')
# 5 Minutos
fiveMins = TimeDelta(val = 300, format='sec')
# 10 Minutos
tenMins  = TimeDelta(val = 600, format='sec')
# 1 Hora
oneHour = TimeDelta(val = 3600, format='sec')
# 30 secs
treintaSegs = TimeDelta(val = 30, format='sec')
# -3 hrs
ART = TimeDelta(val = 10800, format='sec')

###############################################################################
###############################################################################
###############################################################################

def plotMoon(csvFile):
    """
    Plotea muchas cosas del archivo csv generado por makeTableNoContrains
    """
    fd = open(csvFile, 'r')

    vals = csv.reader(fd)

    l = []

    for row in vals:
        l.append(row)

    fd.close()

    x    = []

    for row in l:
        x.append(datetime.datetime.strptime(row[0], '%Y-%m-%d %H:%M:%S'))

    date = x[0]

    year  = date.strftime("%Y")
    month = date.strftime("%B")
        
    ydec = []

    for row in l:
        ydec.append(float(row[5]))

    yHA  = []

    for row in l:
        yHA.append(float(row[2]))

    yElev = []

    for row in l:
        yElev.append(float(row[7]))


    ydeclim  = [NORTHLIMIT] * len(x)
    ydecz    = [lat] * len(x)
    yElevlim = [60] * len(x)
    
    plt.style.use('dark_background')
    plt.plot (x, yHA,     label='Hour Angle',          color='darkslateblue' )
    plt.plot (x, yElev,   label='Elev',        color='green')
    plt.plot (x, ydeclim, label='Lim Dec -10', color='red'       )
    plt.plot (x, ydecz,   label='Zenith',      color='salmon'    )
    plt.plot (x, ydec,    label='Dec',         color='white'    )
    plt.plot (x, yElevlim, label='East Limit', color='crimson')

    plt.title("Moon Chart for IAR %s %s" % (month, year))
    plt.xlabel('Date/Time (UTC)')
    plt.ylabel('Declination/Elevation (deg)')
    plt.grid(color='gray')

    plt.legend()

    fname = "moonchart-%s%s.pdf" % (month, year)
    plt.savefig(fname, format='pdf', quality=95)
    plt.show()

    
def meanElev(csvFile):
    """
    Calcula las distribuciones de elevacion y declinacion de la luna en el IAR
    para un archivo csv generado con makeTableNoConstrains
    """
    fd = open(csvFile, 'r')

    vals = csv.reader(fd)

    l = []

    for row in vals:
        l.append(row)

    fd.close()

    dataelev = []
    datadec  = []

    for row in l:
        if -10 > float(row[5]) > -90:
            print row
            if (-30 < float(row[2]) < -27):# or (27 < float(row[2]) < 30):

#            if (float(row[2]) > -30 and float(row[2]) < -27) or (float(row[2]) > 27 and float(row[2]) < 30):
                dataelev.append(float(row[7]))
                datadec.append(float(row[5]))


    mu, std = norm.fit(dataelev)

    print mu, std

    # Plotear distribucion
    plt.style.use('dark_background')
    hx, hy, _ = plt.hist(dataelev,color="lightblue", label = 'Elevation Distribution')
    hx, hy, _ = plt.hist(datadec, color='red',       label = 'Declination Distribution')

    plt.ylim(0.0,max(hx)+0.05)
    
    plt.title("Moon Elevation for -29 and +29 Hour Angle for IAR")
    plt.xlabel('Declination/Elevation (deg)')
    plt.ylabel('Vals')

    plt.grid(color='gray')
    plt.grid()
    plt.show()

    return 
            
def getVisibleDaysOfMonth(year = None, month = None):
    """
    Devuelve una lista con todos los días en los que la luna está visible
    para la antena del IAR.
    Los parámentros son el año y el mes
    """
    year  = year
    month = month
    day   = 1

    t0 = Time("%d-%02d-%02d 00:00:00" % (year, month, day), format = 'iso',
              scale = 'utc')

    titer = t0

    visibledays = []
    
    while titer.datetime.month == month:
    
        moonpos = astropy.coordinates.get_moon(titer, location=IAR)
        moonHA  = observer.target_hour_angle(target = moonpos, time = titer)
        moonHA  = moonHA.deg
    
        if 180 < moonHA < 360:
	    moonHA = moonHA - 360.
        
        if -10 > moonpos.dec.deg > -90:
            if 29 > moonHA > -29:
                if visibledays.count(titer.datetime.day) == 0:
                    #print "Ventana visibilidad: %d/%02d/%02d" % (titer.datetime.year, titer.datetime.month, titer.datetime.day)
                    visibledays.insert(0, titer.datetime.day)

        titer += oneHour


    return visibledays

def makeTableElevConstrains(year, month):
    day = 1
    
    t0 = Time("%d-%02d-%02d 00:00:00" % (year, month, day), format = 'iso',
              scale = 'utc')

    titer = t0

    while titer.datetime.month == month:
        moonpos = astropy.coordinates.get_moon(titer, location=IAR)
        moonHA  = observer.target_hour_angle(target = moonpos, time = titer)
        moonHA  = moonHA.deg
        moonAltAz = moonpos.transform_to(AltAz(obstime = titer, location = IAR))
        if 180 < moonHA < 360:
	    moonHA = moonHA - 360.

        if (-29 <= moonHA <= -27) or (27 <= moonHA <= 29):
            print "%s,%02.02f,%0.03f,%0.02f,%0.02f" % (titer.datetime,
                                                       moonHA,
                                                       moonpos.dec.deg,
                                                       moonAltAz.az.deg,
                                                       moonAltAz.alt.deg)

        titer += fiveMins

def makeTableNoConstrains(year, month):
    day = 1
    
    t0 = Time("%d-%02d-%02d 00:00:00" % (year, month, day), format = 'iso',
              scale = 'utc')

    titer = t0

    while titer.datetime.month == month:
        moonpos = astropy.coordinates.get_moon(titer, location=IAR)
        moonHA  = observer.target_hour_angle(target = moonpos, time = titer)
        moonHA  = moonHA.deg
        moonAltAz = moonpos.transform_to(AltAz(obstime = titer, location = IAR))
        if 180 < moonHA < 360:
	    moonHA = moonHA - 360.


        print "%s,%f,%02.02f,%02d:%02d:%02d,%0.03f,%0.03f,%0.02f,%0.02f" % (titer.datetime,
                                                                            titer.unix,
                                                                            moonHA,
                                                                            moonpos.ra.hms[0],
                                                                            moonpos.ra.hms[1],
                                                                            moonpos.ra.hms[2],
                                                                            moonpos.ra.deg,
                                                                            moonpos.dec.deg,
                                                                            moonAltAz.az.deg,
                                                                            moonAltAz.alt.deg)

        titer += fiveMins
    

def makeTable(year, month, visibledays):
    """
    Genera una tabla para el año, mes y días visibles con la posición de la luna
    con intervalos de tiempo de 5mins
    """
    while len(visibledays) > 0:
        day = visibledays.pop()
        t0 = Time("%d-%02d-%02d 00:00:00" % (year, month, day), format = 'iso',
              scale = 'utc')
        
        titer = t0
        while titer.datetime.day == day:
            moonpos = astropy.coordinates.get_moon(titer, location=IAR)
            moonHA  = observer.target_hour_angle(target = moonpos, time = titer)
            moonHA  = moonHA.deg
            moonAltAz = moonpos.transform_to(AltAz(obstime = titer, location = IAR))
            if 180 < moonHA < 360:
	        moonHA = moonHA - 360.
        
            if -10 > moonpos.dec.deg > -90:
                if 29 > moonHA > -29:
                    # Datetime tiempoUnix HALuna RAh RAm RAs RAdegs Dec
                    print "%s,%f,%02.02f,%02d:%02d:%02d,%0.03f,%0.03f,%0.02f,%0.02f" % (titer.datetime,
                                                                                        titer.unix,
                                                                                        moonHA,
                                                                                        moonpos.ra.hms[0],
                                                                                        moonpos.ra.hms[1],
                                                                                        moonpos.ra.hms[2],
                                                                                        moonpos.ra.deg,
                                                                                        moonpos.dec.deg,
                                                                                        moonAltAz.az.deg,
                                                                                        moonAltAz.alt.deg)

            titer += fiveMins

def printMoonRises(year, month, visibledays):
    """
    Devuelve las horas UTC de los días en los cuales la luna
    está visible para la antena del IAR.
    """
    pass


###############################################################################
###############################################################################
###############################################################################

def main():
    year  = 2019
    month = 12

    #isibledays = getVisibleDaysOfMonth(year, month)

    #print visibledays
    makeTableNoConstrains(year, month)
    #makeTableElevConstrains(year, month)
    #akeTable(year, month, visibledays)
    
    return True

if __name__ == "__main__":
    sys.exit(main())
