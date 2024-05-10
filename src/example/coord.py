from astropy.coordinates import SkyCoord
import datetime
import astropy.time
import astropy.units as u
from astropy.coordinates import get_sun


def gc_in_ecliptic():
    c = SkyCoord(0, 0, unit="deg", frame="galactic")
    # print(c.galactic)
    # print(c.icrs)
    # print(c.barycentricmeanecliptic)
    cg = c.barycentricmeanecliptic
    # print(str(cg.lon.deg) + "," + str(cg.lat.deg))
    return cg


def print_observation_region_in_icrs():
    c = SkyCoord(0, 0, unit="deg", frame="galactic")
    c0 = SkyCoord(-1.4, 0.6, unit="deg", frame="galactic")
    c1 = SkyCoord(-1.4, -0.6, unit="deg", frame="galactic")
    c2 = SkyCoord(0.7, -0.6, unit="deg", frame="galactic")
    c3 = SkyCoord(0.7, 0.6, unit="deg", frame="galactic")
    ci = c.icrs
    print("galactic center: " + str(ci.ra.deg) + "," + str(ci.dec.deg))
    ci = c0.icrs
    print("galactic center: " + str(ci.ra.deg) + "," + str(ci.dec.deg))
    ci = c1.icrs
    print("galactic center: " + str(ci.ra.deg) + "," + str(ci.dec.deg))
    ci = c2.icrs
    print("galactic center: " + str(ci.ra.deg) + "," + str(ci.dec.deg))
    ci = c3.icrs
    print("galactic center: " + str(ci.ra.deg) + "," + str(ci.dec.deg))


def print_observable_orbit(year: int):
    tz = astropy.time.TimezoneInfo(9 * u.hour)
    toki = datetime.datetime(year, 1, 11, 0, 0, 0, tzinfo=tz)
    toki0 = astropy.time.Time(toki, scale="tcb")
    toki = toki0
    dt = astropy.time.TimeDelta(5801 * u.s)
    observable = 0
    orbit = 0
    cg = gc_in_ecliptic()
    while toki < toki0 + astropy.time.TimeDelta(1 * u.yr):
        taiyou = get_sun(toki)
        taiyou_rel = (taiyou.geocentricmeanecliptic.lon.deg - cg.lon.deg) % 180
        orbit = orbit + 1
        if 45 < taiyou_rel < 135:
            observable = observable + 1
            print(toki)
        toki = toki + dt
    print("orbit: " + str(orbit) + "  observable:" + str(observable))


if __name__ == '__main__':
    c = SkyCoord(0, 90, unit="deg", frame="galactic")
    cg = c.barycentricmeanecliptic
    print(str(cg.lon.deg) + "," + str(cg.lat.deg))
    c = SkyCoord(266.84, 24.46, unit="deg", frame="barycentricmeanecliptic")
    cg = c.galactic
    print(str(cg.l.deg) + "," + str(cg.b.deg))