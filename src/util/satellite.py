import math
from astropy.coordinates import get_sun, SkyCoord
from datetime import datetime, timedelta
import astropy.constants
import astropy.time
import numpy as np


def sky_coord_to_vector(_coord: SkyCoord):
    longitude = _coord.geocentricmeanecliptic.lon.rad
    latitude = _coord.geocentricmeanecliptic.lat.rad
    x_vec = math.cos(longitude) * math.cos(latitude)
    y_vec = math.sin(longitude) * math.cos(latitude)
    z_vec = math.sin(latitude)
    return np.array([x_vec, y_vec, z_vec])


def vector_to_sky_coord(x):
    v_norm = x / np.linalg.norm(x)
    lat = math.asin(v_norm[2])
    lon = math.atan2(x[1], x[0])
    return SkyCoord(lon, lat, unit='rad', frame='geocentricmeanecliptic')


class Satellite:
    def __init__(self, _obs_time: datetime, _phase: float):
        self._base2 = None
        self._base1 = None
        self.satellite_velocity = None
        self._obs_date = _obs_time
        self._phase = _phase
        self.gc_vector = sky_coord_to_vector(SkyCoord(0, 0, unit="deg", frame="galactic"))
        self._revolution_vector = None
        self._orbit_angle = 0
        self.altitude = 600000
        self._fix_date()
        self._revolution_velocity()
        self.orbital_base()

    def _fix_date(self):
        """
        If the specified observation date is out of the date we observe the galactic center, the date is fixed to
         the date when the next observation starts.
        """
        dt = timedelta(days=1)
        inner_product = 1
        while 45 > math.acos(inner_product) * 180 / math.pi or math.acos(inner_product) * 180 / math.pi > 135:
            t0 = astropy.time.Time(self._obs_date)
            sun_vector = sky_coord_to_vector(get_sun(t0))
            inner_product = self.gc_vector.dot(sun_vector)
            if 45 <= math.acos(inner_product) * 180 / math.pi <= 135:
                break
            self._obs_date = self._obs_date + dt
        print(self._obs_date)

    def _revolution_velocity(self):
        """
        Revolution velocity is approximately defined as that of the circular orbit.
        :return: direction is the direction (sky coordinate object). revolution_speed is the speed defined above.
        """
        revolution_speed = math.sqrt(astropy.constants.G.value * astropy.constants.M_sun.value
                                     / astropy.constants.au.value)
        v_longitude = get_sun(astropy.time.Time(self._obs_date)).geocentricmeanecliptic.lon.rad - math.pi / 2
        direction = SkyCoord(v_longitude, 0, frame='geocentricmeanecliptic', unit='rad')
        self._revolution_vector = revolution_speed * sky_coord_to_vector(direction)
        self._orbit_angle = direction.lon.deg \
            - SkyCoord(0, 0, unit="deg", frame="galactic").geocentricmeanecliptic.lon.deg

    def orbital_base(self):
        self._revolution_velocity()
        sun_vector = sky_coord_to_vector(get_sun(astropy.time.Time(self._obs_date)))
        self._base1 = np.cross(self.gc_vector, sun_vector)
        self._base1 = self._base1 / np.linalg.norm(self._base1)
        self._base2 = np.cross(sun_vector, self._base1)
        if self._base2.dot(self.gc_vector) < 0:
            self._base2 = - self._base2

    def calc_satellite_velocity(self, _obs_time=None):
        # inclination = 1.7067125274157757
        # obliquity_of_the_ecliptic = 0.4090926006005829
        if _obs_time is not None:
            time0 = self.obs_date
            phase0 = self.phase
            if _obs_time > time0:
                dt = (_obs_time - time0).seconds
            else:
                dt = - (time0 - _obs_time).seconds
            self._phase = phase0 + dt / 5801 * 2 * math.pi
            self._obs_date = _obs_time
        self.orbital_base()
        altitude = 600000
        orbital_speed = math.sqrt(astropy.constants.G.value * astropy.constants.M_earth.value
                                  / (astropy.constants.R_earth.value + altitude))
        orbital_velocity = (math.sin(self.phase) * self._base1 + math.cos(self.phase) * self._base2) * orbital_speed
        self.satellite_velocity = orbital_velocity + self._revolution_vector
        return self.satellite_velocity

    @property
    def obs_date(self):
        return self._obs_date

    @property
    def phase(self):
        return self._phase

    @property
    def revolution_vector(self):
        return self._revolution_vector

    @property
    def orbit_angle(self):
        return self._orbit_angle

    @property
    def base1(self):
        return self._base1

    @property
    def base2(self):
        return self._base2

    def title_string(self):
        return str(self.obs_date.month) + "/" + str(self.obs_date.day) + "(" + "{:.1f}".format(self.orbit_angle) + \
            "°)" + " phase=" + "{:.0f}".format(np.rad2deg(self.phase)) + "°" + ":snapshot"

    def orbital_speed(self):
        return math.sqrt(astropy.constants.G.value * astropy.constants.M_earth.value
                         / (astropy.constants.R_earth.value + self.altitude))

    def orbital_velocity(self):
        return (math.sin(self._phase) * self._base1 + math.cos(self._phase) * self._base2) * self.orbital_speed()
