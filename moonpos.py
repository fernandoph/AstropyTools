import sys
import astropy.coordinates
from astropy.time import Time, TimezoneInfo
from astropy.time import TimeDelta
import astropy.units as u
from astropy.coordinates import EarthLocation
import astroplan
from time import sleep

# Observatory Location (default IAR)
lat = -34.8601
lon = -58.1385
# Google:
#lat = -34.8661
#lon = -58.1385
h   = 20
IAR = EarthLocation(lat=lat, lon=lon, height=h*u.m)

def main():
    observer = astroplan.Observer(location = IAR, timezone='America/Buenos_Aires')
    localtz =  TimezoneInfo(utc_offset=-3 * u.hour)
    
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


    year  = 2019
    month = 10
    day   = 1

    t0 = Time("%d-%02d-%02d 00:00:00" % (year, month, day), format = 'iso',
              scale = 'utc')

    titer = t0

    while True:
        moonpos = astropy.coordinates.get_moon(Time.now(), location = IAR)
        moonHA  = observer.target_hour_angle(target = moonpos, time=Time.now())
        moonHA  = moonHA.deg

        print '================================================================'
        print 'RA: ', moonpos.ra.hms[0], moonpos.ra.hms[1], moonpos.ra.hms[2]
        print 'DEC: ', moonpos.dec.dms[0], moonpos.dec.dms[1], moonpos.dec.dms[2]
        print 'HA: ', moonHA

        sleep(1)
    
    

    


if __name__ == "__main__":
    sys.exit(main())
