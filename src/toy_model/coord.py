from astropy.coordinates import SkyCoord
import datetime
import astropy.time
import astropy.units as u
from astropy.coordinates import get_sun

if __name__ == '__main__':
    c = SkyCoord(0, 0, unit="deg", frame="galactic")
    # print(c.galactic)
    # print(c.icrs)
    # print(c.barycentricmeanecliptic)
    cg = c.barycentricmeanecliptic.lon.deg
    print(cg)
    tz = astropy.time.TimezoneInfo(9 * u.hour)
    toki = datetime.datetime(2028, 1, 11, 0, 0, 0, tzinfo=tz)
    toki0 = astropy.time.Time(toki, scale="tcb")
    toki = toki0
    # print(toki)
    dt = astropy.time.TimeDelta(5801 * u.s)
    observable = 0
    orbit = 0
    while toki < toki0 + astropy.time.TimeDelta(1 * u.yr):
        taiyou = get_sun(toki)
        taiyou_rel = (taiyou.geocentricmeanecliptic.lon.deg - cg) % 180
        orbit = orbit + 1
        if 45 < taiyou_rel < 135:
            observable = observable + 1
            print(toki)
        toki = toki + dt
    print("orbit: " + str(orbit) + "  observable:" + str(observable))
