from astropy.time import Time, TimeDelta
from astropy.coordinates import solar_system_ephemeris, EarthLocation, SkyCoord
from astropy.coordinates import get_body
import astropy.units as u


def regions() -> tuple:
    c0: SkyCoord = SkyCoord(-1.4, 0.6, unit="deg", frame="galactic")
    c1: SkyCoord = SkyCoord(-1.4, -0.6, unit="deg", frame="galactic")
    c2: SkyCoord = SkyCoord(0.7, -0.6, unit="deg", frame="galactic")
    c3: SkyCoord = SkyCoord(0.7, 0.6, unit="deg", frame="galactic")
    return c0, c1, c2, c3


def planetTrajectory(times: list, planet: str) -> None:
    loc: EarthLocation = EarthLocation.of_site('greenwich')
    with solar_system_ephemeris.set('builtin'):
        for tt in times:
            print(get_body(planet, tt, loc).galactic)


if __name__ == '__main__':
    t0 = Time("2028-01-01T01:00:00.00", scale="tcb")
    dt = TimeDelta(30 * u.day)
    t: list = [t0]
    for i in range(10):
        t0 = t0 + dt
        t.append(t0)
    print(t)
    planetTrajectory(t, 'jupiter')
