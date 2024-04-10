import math
from datetime import datetime, timedelta
import astropy.time
import astropy.constants
import astropy.units as u
import numpy as np
from astropy.coordinates import get_sun, SkyCoord
from matplotlib import pyplot as plt


class VectorField2d:
    """
    2 dimensional vector field class. Constructor gives the mesh grid of x and y coordinate.
    """

    def __init__(self, _lx: float, _ly: float, _grid_width: float):
        """

        :param _lx: maximum x coordinate value, region is -_lx to _lx
        :param _ly: maximum y coordinate value, region is -_ly to _ly
        :param _grid_width: grid width
        """
        self.LX = _lx
        self.LY = _ly
        self.grid_width = _grid_width
        self.X, self.Y = np.meshgrid(np.arange(-self.LX, self.LX * 1.01, self.grid_width),
                                     np.arange(-self.LY, self.LY * 1.01, self.grid_width))

    def get_grid(self):
        """

        :return: array of mesh grid X and Y
        """
        return self.X, self.Y

    def draw_limits(self, _plt):
        """
        draw_limits function set the limits of x and y coordinate.
        :param _plt: plotter instance
        """
        _plt.set_xticks([-0.3, -0.15, 0, 0.15, 0.3])
        _plt.set_yticks([-0.3, -0.15, 0, 0.15, 0.3])
        _plt.set_xticklabels([-0.3, -0.15, 0, 0.15, 0.3], rotation=90)
        _plt.set_xlim([-self.LX*1.1, self.LX*1.1])
        _plt.set_ylim([-self.LY*1.1, self.LY*1.1])

    def draw(self, _u, _v, _plt, _title_string, _arrow_scale):
        """
        Draw vector field specified in the array _u and _v. scale is automatically defined and printed to the console.
        :param _u: x component of the vector field stored in the same type array as X
        :param _v: y component of the vector field stored in the same type array as Y
        :param _plt: plotter instance
        :param _title_string: graph title
        :param _arrow_scale: scale of arrow drawn on the graph
        """
        # print("max = " + str(np.max(_u)) + "deg," + str(np.max(_u) * 3600) + "as")
        print(str(np.max(_u) * 3600), end=',')
        _scale = max(np.max(_u), -np.min(_u), np.max(_v), -np.min(_v)) / self.grid_width
        _arrow = _arrow_scale / _scale / 3600
        # self.draw_limits(_plt)
        if _arrow_scale < 1e-3:
            arrow_string = str(_arrow_scale * 1e6) + "uas"
        elif _arrow_scale < 1:
            arrow_string = str(_arrow_scale * 1e3) + "mas"
        else:
            arrow_string = str(_arrow_scale) + "as"
        _plt.set_aspect('equal')
        self.draw_limits(_plt)
        _plt.set_xlabel("l")
        _plt.set_ylabel("b")
        _plt.arrow(x=0, y=-self.LY, dx=_arrow, dy=0, head_width=0.25 * _arrow, length_includes_head=True)
        _plt.text(_arrow, -self.LY, arrow_string)
        _plt.quiver(self.X, self.Y, _u, _v, color='red', angles='xy', scale_units='xy', scale=_scale)
        _plt.set_title(_title_string)
        _plt.grid()

    def draw_diff(self, _u, _v, _plt, _title_string, _arrow_scale):
        """
        Draw difference of vector field itself and mean.
        :param _u: x component of the vector field stored in the same type array as X
        :param _v: y component of the vector field stored in the same type array as Y
        :param _plt: plotter instance
        :param _title_string: graph title
        :param _arrow_scale: scale of arrow drawn on the graph
        """
        _uu = _u - np.average(_u)
        _vv = _v - np.average(_v)
        # print("max = " + str(np.max(_uu)) + "deg," + str(np.max(_uu) * 3600) + "as")
        print(str(np.max(_uu) * 3600), end=',')
        _scale = max(np.max(_uu), -np.min(_uu), np.max(_vv), -np.min(_vv)) / self.grid_width
        _arrow = _arrow_scale / _scale / 3600
        if _arrow_scale < 1e-3:
            arrow_string = str(_arrow_scale * 1e6) + "uas"
        elif _arrow_scale < 1:
            arrow_string = str(_arrow_scale * 1e3) + "mas"
        else:
            arrow_string = str(_arrow_scale) + "as"
        _plt.set_aspect('equal')
        self.draw_limits(_plt)
        _plt.set_xlabel("l")
        _plt.set_ylabel("b")
        _plt.arrow(x=0, y=-self.LY * 1.05, dx=_arrow, dy=0, head_width=0.25 * _arrow, length_includes_head=True)
        _plt.text(_arrow, -self.LY * 1.05, arrow_string)
        _plt.quiver(self.X, self.Y, _uu, _vv, color='red', angles='xy', scale_units='xy', scale=_scale)
        _plt.set_title(_title_string)
        _plt.grid()


def calc_distortion_in_degree(_v, _vector_field):
    """

    :param _v: velocity vector in 3D geocentric mean ecliptic coordinate and in m/s
    :param _vector_field: instance of the VectorField2D class
    :return: aberration in degree
    """
    v_unity = _v / np.linalg.norm(_v)
    c = vector_to_sky_coord(v_unity)
    _x0 = c.galactic.l.deg
    _y0 = c.galactic.b.deg
    # print(str(_x0) + "," + str(_y0))
    _v0 = np.linalg.norm(_v) / astropy.constants.c.value

    _x, _y = _vector_field.get_grid()
    r0 = np.sqrt((_x - _x0) ** 2 + (_y - _y0) ** 2)
    r1 = r0 * math.pi / 180.0  # r0 in radian
    _shift = _v0 * np.sin(r1) * (1 + 2.25 * _v0 * _v0) - 0.25 * _v0 * _v0 * np.sin(2 * r1)\
        + _v0 * _v0 * _v0 * np.sin(3 * r1) / 12.0
    _u = - _shift * (_x - _x0) / r1
    _v = - _shift * (_y - _y0) / r1
    return _u, _v


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


class RelativisticAberration:
    def __init__(self, _obs_time: datetime, _phase: float):
        """
        Constructor of the class RelativisticAberration
        :param _obs_time: observation time as datetime object.
        """
        self.satellite_velocity = None
        # tz = astropy.time.TimezoneInfo(0 * u.hour)
        self.obs_date = _obs_time
        self.phase = _phase
        self.gc_vector = sky_coord_to_vector(SkyCoord(0, 0, unit="deg", frame="galactic"))
        self.revolution_vector = None
        self.base1 = None
        self.base2 = None
        self.orbit_angle = 0
        self._fix_date()
        self._revolution_velocity()
        self._orbital_base()

    def _fix_date(self):
        """
        If the specified observation date is out of the date we observe the galactic center, the date is fixed to
         the date when the next observation starts.
        """
        dt = timedelta(days=1)
        inner_product = 1
        while 45 > math.acos(inner_product) * 180 / math.pi or math.acos(inner_product) * 180 / math.pi > 135:
            t0 = astropy.time.Time(self.obs_date)
            sun_vector = sky_coord_to_vector(get_sun(t0))
            inner_product = self.gc_vector.dot(sun_vector)
            if 45 <= math.acos(inner_product) * 180 / math.pi <= 135:
                break
            self.obs_date = self.obs_date + dt
        print(self.obs_date)

    def _revolution_velocity(self):
        """
        Revolution velocity is approximately defined as that of the circular orbit.
        :return: direction is the direction (sky coordinate object). revolution_speed is the speed defined above.
        """
        revolution_speed = math.sqrt(astropy.constants.G.value * astropy.constants.M_sun.value
                                     / astropy.constants.au.value)
        v_longitude = get_sun(astropy.time.Time(self.obs_date)).geocentricmeanecliptic.lon.rad - math.pi / 2
        direction = SkyCoord(v_longitude, 0, frame='geocentricmeanecliptic', unit='rad')
        self.revolution_vector = revolution_speed * sky_coord_to_vector(direction)
        self.orbit_angle = direction.lon.deg\
                           - SkyCoord(0, 0, unit="deg", frame="galactic").geocentricmeanecliptic.lon.deg

    def _orbital_base(self):
        sun_vector = sky_coord_to_vector(get_sun(astropy.time.Time(self.obs_date)))
        self.base1 = np.cross(self.gc_vector, sun_vector)
        self.base1 = self.base1 / np.linalg.norm(self.base1)
        # print("base1 = " + str(base1))
        self.base2 = np.cross(sun_vector, self.base1)
        if self.base2.dot(self.gc_vector) < 0:
            self.base2 = - self.base2

    def _calc_satellite_velocity(self, _obs_time=None):
        # inclination = 1.7067125274157757
        # obliquity_of_the_ecliptic = 0.4090926006005829
        if _obs_time is not None:
            time0 = self.obs_date
            phase0 = self.phase
            if _obs_time > time0:
                dt = (_obs_time - time0).seconds
            else:
                dt = - (time0 - _obs_time).seconds
            self.phase = phase0 + dt / 5801 * 2 * math.pi
            self.obs_date = _obs_time
        self._revolution_velocity()
        self._orbital_base()
        altitude = 600000
        orbital_speed = math.sqrt(astropy.constants.G.value * astropy.constants.M_earth.value
                                  / (astropy.constants.R_earth.value + altitude))
        orbital_velocity = (math.sin(self.phase) * self.base1 + math.cos(self.phase) * self.base2) * orbital_speed
        self.satellite_velocity = orbital_velocity + self.revolution_vector

    def _figure_title(self):
        fig = plt.figure(tight_layout=True)
        fig.suptitle(str(self.obs_date.month) + "/" + str(self.obs_date.day) + "(" + "{:.1f}".format(self.orbit_angle)
                     + "°)" + " phase=" + "{:.0f}".format(np.rad2deg(self.phase)) + "°" + ":snapshot")
        return fig

    def figure2(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(2, 2, 1)
        fig2 = fig.add_subplot(2, 2, 2)
        self._calc_satellite_velocity()
        v = self.satellite_velocity
        u0, v0 = calc_distortion_in_degree(v, vector_field)
        vector_field.draw(u0, v0, fig1, "aberration field", 1)
        vector_field.draw_diff(u0, v0, fig2, "difference from mean", 0.1)
        # fig.tight_layout()
        plt.draw()
        plt.savefig('fig2.png')
        plt.show()

    def figure4(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(2, 2, 1)
        fig2 = fig.add_subplot(2, 2, 2)
        v = self.satellite_velocity
        u0, v0 = calc_distortion_in_degree(v, vector_field)
        altitude = 600000
        orbital_speed = math.sqrt(astropy.constants.G.value * astropy.constants.M_earth.value
                                  / (astropy.constants.R_earth.value + altitude))
        orbital_velocity = (math.sin(self.phase) * self.base1 + math.cos(self.phase) * self.base2) * orbital_speed
        # self.satellite_velocity = orbital_velocity + self.revolution_vector
        c = vector_to_sky_coord(orbital_velocity)
        dl = 250 / 3600  # unit degree, /3600 means the unit in arc seconds
        ov_with_dl = sky_coord_to_vector(SkyCoord(l=c.galactic.l.deg + dl, b=c.galactic.b.deg, unit='deg',
                                                  frame='galactic')) * orbital_speed
        v_with_dl = ov_with_dl + self.revolution_vector
        u1, v1 = calc_distortion_in_degree(v_with_dl, vector_field)
        # vector_field.draw(u1 - u0, v1 - v1, fig1, "δl=" +str(dl * 3600) + "arc sec", 2.5e-5)  # error of direction
        vector_field.draw_diff(u1 - u0, v1 - v1, fig1, "δl=" +str(dl * 3600) + "arc sec", 2.5e-5)  # error of direction
        dv = 80  # m/s
        ov_with_dv = orbital_velocity * (orbital_speed + dv) / orbital_speed
        v_with_dv = ov_with_dv + self.revolution_vector
        u1, v1 = calc_distortion_in_degree(v_with_dv, vector_field)
        # vector_field.draw(u1 - u0, v1 - v1, fig2, "d|v|=" + str(dv) + "m/s", 2.5e-5)  # error of speed
        vector_field.draw_diff(u1 - u0, v1 - v1, fig2, "d|v|=" + str(dv) + "m/s", 2.5e-5)  # error of speed
        plt.draw()
        plt.savefig('fig4.png')
        plt.show()

    def figure5(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(2, 2, 1)
        fig2 = fig.add_subplot(2, 2, 2)
        # fig3 = fig.add_subplot(3, 2, 3)
        # fig4 = fig.add_subplot(3, 2, 4)
        fig5 = fig.add_subplot(2, 2, 3)
        fig6 = fig.add_subplot(2, 2, 4)
        t0 = self.obs_date
        self._calc_satellite_velocity()
        um, vm = calc_distortion_in_degree(self.satellite_velocity, vector_field)
        vector_field.draw(um, vm, fig1, "v(tm): tm=(t1+t2)/2", 7.0)  # aberration field 1. itself as top left
        _phase = self.phase
        dt = timedelta(seconds=6.25)
        self._calc_satellite_velocity(t0 - dt)
        us, vs = calc_distortion_in_degree(self.satellite_velocity, vector_field)
        self._calc_satellite_velocity(t0 + dt)
        ue, ve = calc_distortion_in_degree(self.satellite_velocity, vector_field)
        vector_field.draw(ue - us, ve - vs, fig2, "v(t2) - v(t1)", 0.02)
        # vector_field.draw(ue - np.average(um), ve - np.average(vm), fig3, "v(t2) - <v(tm)>", 0.1)
        # vector_field.draw(ue - um, ve - vm, fig4, "v(t2) - v(tm)", 0.005)
        vector_field.draw((ue + us) * 0.5 - um, (ve + vs) * 0.5 - vm, fig5, "<v> - v(tm). 12.5sec", 0.00010)
        dt = timedelta(seconds=3)
        self._calc_satellite_velocity(t0 - dt)
        us, vs = calc_distortion_in_degree(self.satellite_velocity, vector_field)
        self._calc_satellite_velocity(t0 + dt)
        ue, ve = calc_distortion_in_degree(self.satellite_velocity, vector_field)
        vector_field.draw((ue + us) * 0.5 - um, (ve + vs) * 0.5 - vm, fig6, "<v> - v(tm). 6sec", 0.000025)
        plt.draw()
        plt.savefig('fig5.png')
        plt.show()


a = VectorField2d(0.3, 0.3, 0.1)
tz = astropy.time.TimezoneInfo(0 * u.hour)
r = RelativisticAberration(datetime(2024, 2, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(120))
r.figure2(a)
r.figure4(a)
r = RelativisticAberration(datetime(2024, 3, 24, 12, 0, 0, tzinfo=tz), np.deg2rad(110))
r.figure5(a)
"""
print("fig2")
r.figure2(a)
print()
print("fig4")
r.figure4(a)
"""
# r = RelativisticAberration(datetime(2024, 2, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(120))
# r.figure2(a)
