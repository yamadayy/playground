import math
from datetime import datetime, timedelta
import astropy.time
import astropy.constants
import cartopy.crs as ccrs
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib import pyplot as plt
from util.vector_field_2d import VectorField2d
from util.satellite import vector_to_sky_coord, sky_coord_to_vector, Satellite


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
        fig.suptitle(self.satellite.title_string())
        return fig

    def figure3(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
        fig2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
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
        v = self.satellite.satellite_velocity
        u0, v0 = calc_distortion_in_degree(v, vector_field)
        vector_field.draw(u0, v0, fig1, "aberration field", 15)
        vector_field.draw(u0 - np.average(u0), v0 - np.average(v0), fig2, "difference from mean", 0.1)
        c = vector_to_sky_coord(self.satellite.orbital_velocity())
        dl = 1.2 / 3600  # unit degree, /3600 means the unit in arc seconds
        ov_with_dl = sky_coord_to_vector(SkyCoord(l=c.galactic.l.deg + dl, b=c.galactic.b.deg, unit='deg',
                                                  frame='galactic')) * self.satellite.orbital_speed()
        v_with_dl = ov_with_dl + self.satellite.revolution_vector
        u1, v1 = calc_distortion_in_degree(v_with_dl, vector_field)
        vector_field.draw(u1 - u0, v1 - v0, fig3, "δl=" + str(dl * 3600) + "arc sec", 2e-5)  # error of direction
        dv = 0.022  # m/s
        ov_with_dv = self.satellite.orbital_velocity() * (self.satellite.orbital_speed() + dv) / self.satellite.orbital_speed()
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
        v = self.satellite.calc_satellite_velocity()
        u0, v0 = calc_distortion_in_degree(v, vector_field)
        c = vector_to_sky_coord(self.satellite.orbital_velocity())
        dl = 160 / 3600  # unit degree, /3600 means the unit in arc seconds
        ov_with_dl = sky_coord_to_vector(SkyCoord(l=c.galactic.l.deg + dl, b=c.galactic.b.deg, unit='deg',
                                                  frame='galactic')) * self.satellite.orbital_speed()
        v_with_dl = ov_with_dl + self.satellite.revolution_vector
        u1, v1 = calc_distortion_in_degree(v_with_dl, vector_field)
        vector_field.draw(u1 - u0 - np.average(u1) + np.average(u0), v1 - v0 - np.average(v1) + np.average(v0), fig1,
                          "δl=" + str(dl * 3600) + "arc sec", 2e-5)  # error of direction
        dv = 15  # m/s
        ov_with_dv = self.satellite.orbital_velocity() * (self.satellite.orbital_speed() + dv) / self.satellite.orbital_speed()
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
    _ang_b_cos = (np.sin(y0) - np.sin(np.deg2rad(_Y)) * _len_cos) / np.sin(_len_c) / np.cos(np.deg2rad(_Y))
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
