# -*- coding: utf-8 -*-
import sys
import astropy.coordinates
from astropy.coordinates import AltAz
from astropy.time import Time, TimezoneInfo
from astropy.time import TimeDelta
import astropy.units as u
from astropy.coordinates import EarthLocation
import astroplan
from astropy.coordinates import solar_system_ephemeris
solar_system_ephemeris.set('jpl') 

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
# 1 Hora
oneHour = TimeDelta(val = 3600, format='sec')
# 30 secs
treintaSegs = TimeDelta(val = 30, format='sec')
# -3 hrs
ART = TimeDelta(val = 10800, format='sec')

###############################################################################
###############################################################################
###############################################################################

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
    month = 11

    visibledays = getVisibleDaysOfMonth(year, month)

    #print visibledays
    
    makeTable(year, month, visibledays)
    
    return True

if __name__ == "__main__":
    sys.exit(main())
