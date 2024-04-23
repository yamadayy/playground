import math
from datetime import datetime, timedelta
import astropy.time
import astropy.constants
import cartopy.crs as ccrs
import numpy as np
from astropy.coordinates import get_sun, SkyCoord
import astropy.units as u
from matplotlib import pyplot as plt
from util.vector_field_2d import VectorField2d


def calc_distortion_in_degree(_v, _vector_field):
    v_unity = _v / np.linalg.norm(_v)
    c = vector_to_sky_coord(v_unity)
    _x0 = c.galactic.l.rad
    _y0 = c.galactic.b.rad
    print(str(np.rad2deg(_x0)) + "," + str(np.rad2deg(_y0)), end=",")
    _v0 = np.linalg.norm(_v) / astropy.constants.c.value
    _x, _y = _vector_field.get_grid()
    len_cos = np.sin(np.deg2rad(_y))*np.sin(_y0) + np.cos(np.deg2rad(_y)) * np.cos(_y0) * np.cos(np.deg2rad(_x) - _x0)
    len_c = np.arccos(len_cos)
    ang_b_sin = np.sin(_x0 - np.deg2rad(_x)) * np.cos(np.deg2rad(_y)) / np.sin(len_c)
    ang_b_cos = (np.sin(_y0) - np.sin(np.deg2rad(_y)) * len_cos) / np.sin(len_c) / np.cos(np.deg2rad(_y))
    _shift = _v0 * np.sin(len_c) * (1 + 2.25 * _v0 * _v0) - 0.25 * _v0 * _v0 * np.sin(2 * len_c) \
        + _v0 * _v0 * _v0 * np.sin(3 * len_c) / 12.0
    _uu = _shift * ang_b_sin * 180 / math.pi
    _vv = _shift * ang_b_cos * 180 / math.pi
    return _uu, _vv


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


class RelativisticAberration:
    def __init__(self, _sat: Satellite):
        """
        Constructor of the class RelativisticAberration
        """
        self.satellite = _sat
        # self.revolution_vector = None
        # self.base1 = None
        # self.base2 = None
        # self.orbit_angle = 0

    def _figure_title(self):
        fig = plt.figure(tight_layout=True)
        fig.suptitle(str(self.satellite.obs_date.month) + "/" + str(self.satellite.obs_date.day) + "("
                     + "{:.1f}".format(self.satellite.orbit_angle)
                     + "°)" + " phase=" + "{:.0f}".format(np.rad2deg(self.satellite.phase)) + "°" + ":snapshot")
        return fig

    def figure3(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
        fig2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
        self.satellite.calc_satellite_velocity()
        v = self.satellite.calc_satellite_velocity()
        u0, v0 = calc_distortion_in_degree(v, vector_field)
        vector_field.draw(u0, v0, fig1, "aberration field", 15)
        vector_field.draw(u0 - np.average(u0), v0 - np.average(v0), fig2, "difference from mean", 0.1)
        plt.draw()
        plt.savefig('fig3.png')
        plt.show()

    def figure9_left_down(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(2, 2, 1, projection=ccrs.PlateCarree())
        fig2 = fig.add_subplot(2, 2, 2, projection=ccrs.PlateCarree())
        fig3 = fig.add_subplot(2, 2, 3, projection=ccrs.PlateCarree())
        fig4 = fig.add_subplot(2, 2, 4, projection=ccrs.PlateCarree())
        self.satellite.calc_satellite_velocity()
        v = self.satellite.satellite_velocity
        u0, v0 = calc_distortion_in_degree(v, vector_field)
        vector_field.draw(u0, v0, fig1, "aberration field", 15)
        vector_field.draw(u0 - np.average(u0), v0 - np.average(v0), fig2, "difference from mean", 0.1)
        altitude = 600000
        orbital_speed = math.sqrt(astropy.constants.G.value * astropy.constants.M_earth.value
                                  / (astropy.constants.R_earth.value + altitude))
        orbital_velocity = (math.sin(self.satellite.phase) * self.satellite.base1 + math.cos(self.satellite.phase)
                            * self.satellite.base2) * orbital_speed
        c = vector_to_sky_coord(orbital_velocity)
        dl = 1.2 / 3600  # unit degree, /3600 means the unit in arc seconds
        ov_with_dl = sky_coord_to_vector(SkyCoord(l=c.galactic.l.deg + dl, b=c.galactic.b.deg, unit='deg',
                                                  frame='galactic')) * orbital_speed
        v_with_dl = ov_with_dl + self.satellite.revolution_vector
        u1, v1 = calc_distortion_in_degree(v_with_dl, vector_field)
        vector_field.draw(u1 - u0, v1 - v0, fig3, "δl=" + str(dl * 3600) + "arc sec", 2e-5)  # error of direction
        dv = 0.022  # m/s
        ov_with_dv = orbital_velocity * (orbital_speed + dv) / orbital_speed
        v_with_dv = ov_with_dv + self.satellite.revolution_vector
        u1, v1 = calc_distortion_in_degree(v_with_dv, vector_field)
        vector_field.draw(u1 - u0, v1 - v0, fig4, "d|v|=" + str(dv) + "m/s", 2e-5)  # error of speed
        plt.draw()
        plt.savefig('fig9ld.png')
        plt.show()

    def figure7(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
        fig2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
        self.satellite.calc_satellite_velocity()
        v = self.satellite.satellite_velocity
        u0, v0 = calc_distortion_in_degree(v, vector_field)
        altitude = 600000
        orbital_speed = math.sqrt(astropy.constants.G.value * astropy.constants.M_earth.value
                                  / (astropy.constants.R_earth.value + altitude))
        orbital_velocity = (math.sin(self.satellite.phase) * self.satellite.base1
                            + math.cos(self.satellite.phase) * self.satellite.base2) * orbital_speed
        c = vector_to_sky_coord(orbital_velocity)
        dl = 160 / 3600  # unit degree, /3600 means the unit in arc seconds
        ov_with_dl = sky_coord_to_vector(SkyCoord(l=c.galactic.l.deg + dl, b=c.galactic.b.deg, unit='deg',
                                                  frame='galactic')) * orbital_speed
        v_with_dl = ov_with_dl + self.satellite.revolution_vector
        u1, v1 = calc_distortion_in_degree(v_with_dl, vector_field)
        vector_field.draw(u1 - u0 - np.average(u1) + np.average(u0), v1 - v0 - np.average(v1) + np.average(v0), fig1,
                          "δl=" + str(dl * 3600) + "arc sec", 2e-5)  # error of direction
        dv = 15  # m/s
        ov_with_dv = orbital_velocity * (orbital_speed + dv) / orbital_speed
        v_with_dv = ov_with_dv + self.satellite.revolution_vector
        u1, v1 = calc_distortion_in_degree(v_with_dv, vector_field)
        vector_field.draw(u1 - u0 - np.average(u1) + np.average(u0), v1 - v0 - np.average(v1) + np.average(v0), fig2,
                          "d|v|=" + str(dv) + "m/s", 2e-5)  # error of speed
        plt.draw()
        plt.savefig('fig7.png')
        plt.show()

    def figure5(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(2, 2, 1, projection=ccrs.PlateCarree())
        fig2 = fig.add_subplot(2, 2, 2, projection=ccrs.PlateCarree())
        # fig3 = fig.add_subplot(3, 2, 3, projection=ccrs.PlateCarree())
        # fig4 = fig.add_subplot(3, 2, 4, projection=ccrs.PlateCarree())
        fig5 = fig.add_subplot(2, 2, 3, projection=ccrs.PlateCarree())
        fig6 = fig.add_subplot(2, 2, 4, projection=ccrs.PlateCarree())
        t0 = self.satellite.obs_date
        self.satellite.calc_satellite_velocity()
        um, vm = calc_distortion_in_degree(self.satellite.satellite_velocity, vector_field)
        vector_field.draw(um, vm, fig1, "u(tm): tm=(t1+t2)/2", 15.0)  # aberration field 1. itself as top left
        _phase = self.satellite.phase
        dt = timedelta(seconds=6.25)
        self.satellite.calc_satellite_velocity(t0 - dt)
        us, vs = calc_distortion_in_degree(self.satellite.satellite_velocity, vector_field)
        self.satellite.calc_satellite_velocity(t0 + dt)
        ue, ve = calc_distortion_in_degree(self.satellite.satellite_velocity, vector_field)
        vector_field.draw(ue - us, ve - vs, fig2, "u(t2) - u(t1)", 0.05)
        # vector_field.draw(ue - np.average(um), ve - np.average(vm), fig3, "v(t2) - <v(tm)>", 0.1)
        # vector_field.draw(ue - um, ve - vm, fig4, "v(t2) - v(tm)", 0.005)
        vector_field.draw((ue + us) * 0.5 - um, (ve + vs) * 0.5 - vm, fig5, "<u> - u(tm). 12.5sec", 0.00015)
        dt = timedelta(seconds=2.499)
        self.satellite.calc_satellite_velocity(t0 - dt)
        us, vs = calc_distortion_in_degree(self.satellite.satellite_velocity, vector_field)
        self.satellite.calc_satellite_velocity(t0 + dt)
        ue, ve = calc_distortion_in_degree(self.satellite.satellite_velocity, vector_field)
        vector_field.draw((ue + us) * 0.5 - um, (ve + vs) * 0.5 - vm, fig6, "<u> - u(tm). 5sec", 0.000020)
        plt.draw()
        plt.savefig('fig5.png')
        plt.show()

    def figure10(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(2, 2, 1, projection=ccrs.PlateCarree())
        fig2 = fig.add_subplot(2, 2, 2, projection=ccrs.PlateCarree())
        fig5 = fig.add_subplot(2, 2, 3, projection=ccrs.PlateCarree())
        fig6 = fig.add_subplot(2, 2, 4, projection=ccrs.PlateCarree())
        t0 = self.satellite.obs_date
        self.satellite.calc_satellite_velocity()
        um, vm = calc_distortion_in_degree(self.satellite.satellite_velocity, vector_field)
        vector_field.draw(um - np.average(um), vm - np.average(vm), fig1, "u(tm): tm=(t1+t2)/2", 0.1)
        _phase = self.satellite.phase
        dt = timedelta(seconds=6.25)
        self.satellite.calc_satellite_velocity(t0 - dt)
        us, vs = calc_distortion_in_degree(self.satellite.satellite_velocity, vector_field)
        self.satellite.calc_satellite_velocity(t0 + dt)
        ue, ve = calc_distortion_in_degree(self.satellite.satellite_velocity, vector_field)
        vector_field.draw(ue - us - np.average(ue - us), ve - vs - np.average(ve - vs), fig2, "u(t2) - u(t1)", 0.00035)
        # vector_field.draw(ue - np.average(um), ve - np.average(vm), fig3, "v(t2) - <v(tm)>", 0.1)
        # vector_field.draw(ue - um, ve - vm, fig4, "v(t2) - v(tm)", 0.005)
        vector_field.draw((ue + us) * 0.5 - um - np.average((ue + us) * 0.5 - um),
            (ve + vs) * 0.5 - vm - np.average((ve + vs) * 0.5 - vm), fig5, "<u> - u(tm). 12.5sec", 2e-7)
        dt = timedelta(seconds=60)
        self.satellite.calc_satellite_velocity(t0 - dt)
        us, vs = calc_distortion_in_degree(self.satellite.satellite_velocity, vector_field)
        self.satellite.calc_satellite_velocity(t0 + dt)
        ue, ve = calc_distortion_in_degree(self.satellite.satellite_velocity, vector_field)
        vector_field.draw((ue + us) * 0.5 - um - np.average((ue + us) * 0.5 - um),
            (ve + vs) * 0.5 - vm - np.average((ve + vs) * 0.5 - vm), fig6, "<u> - u(tm). 120sec", 2e-5)
        plt.draw()
        plt.savefig('fig10.png')
        plt.show()


def global_map(x0, y0, _v0):
    _fig = plt.figure()
    _fig1 = _fig.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
    _fig1.set_global()
    _fig1.gridlines(draw_labels=True)
    _x = 10 * np.arange(36 + 1) - 180.0
    _y = 90. - 10 * np.arange(18 + 1)
    _X, _Y = np.meshgrid(_x, _y)
    _len_cos = np.sin(np.deg2rad(_Y)) * np.sin(y0) + np.cos(np.deg2rad(_Y)) * np.cos(y0) * np.cos(np.deg2rad(_X) - x0)
    _len_c = np.arccos(_len_cos)
    _ang_b_sin = np.sin(x0 - np.deg2rad(_X)) * np.cos(np.deg2rad(_Y)) / np.sin(_len_c)
    _ang_b_cos =  (np.sin(y0) - np.sin(np.deg2rad(_Y)) * _len_cos) / np.sin(_len_c) / np.cos(np.deg2rad(_Y))
    _shift = _v0 * np.sin(_len_c) * (1 + 2.25 * _v0 * _v0) - 0.25 * _v0 * _v0 * np.sin(2 * _len_c) \
        + _v0 * _v0 * _v0 * np.sin(3 * _len_c) / 12.0
    _u = _shift * _ang_b_sin
    _v = _shift * _ang_b_cos
    _fig1.quiver(_X, _Y, _u, _v, color='b', transform=ccrs.PlateCarree())
    plt.draw()
    plt.savefig("global.png")
    plt.show()


global_map(0, 0, 0.0001)

a = VectorField2d(0.3, 0.3, 0.1)
tz = astropy.time.TimezoneInfo(0 * u.hour)
s = Satellite(datetime(2024, 2, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(120))
r = RelativisticAberration(s)
r.figure3(a)
r.figure9_left_down(a)
r.figure5(a)
r.figure10(a)
r = RelativisticAberration(Satellite(datetime(2024, 3, 17, 12, 0, 0, tzinfo=tz), np.deg2rad(110)))
r.figure7(a)
