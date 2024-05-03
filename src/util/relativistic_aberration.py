from util.satellite import vector_to_sky_coord, sky_coord_to_vector, Satellite
from astropy.coordinates import SkyCoord
import numpy as np
import astropy.constants
import math
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from datetime import timedelta


def subtract_trend(_u, _x):
    u_0 = _u[0, 0]
    u_1 = _u[-1, -1]
    x0 = _x[0, 0]
    x1 = _x[-1, -1]
    _a = (u_1 - u_0) / (x1 - x0)
    _b = (u_0 * x1 - u_1 * x0) / (x1 - x0)
    return _u - _a * _x - _b


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

    def calc_distortion_in_degree(self, _vector_field, _dl=0.0, _dv=0.0, _obs_time=None):
        self.satellite.calc_satellite_velocity(_obs_time=_obs_time)
        _v = self.satellite.satellite_velocity
        if _dl != 0:
            c = vector_to_sky_coord(self.satellite.orbital_velocity())
            ov_with_dl = sky_coord_to_vector(SkyCoord(l=c.galactic.l.deg + _dl, b=c.galactic.b.deg, unit='deg',
                                                      frame='galactic')) * self.satellite.orbital_speed()
            _v = ov_with_dl + self.satellite.revolution_vector
        if _dv != 0:
            ov_with_dv = self.satellite.orbital_velocity() * (
                        self.satellite.orbital_speed() + _dv) / self.satellite.orbital_speed()
            _v = ov_with_dv + self.satellite.revolution_vector
        v_unity = _v / np.linalg.norm(_v)
        c = vector_to_sky_coord(v_unity)
        _x0 = c.galactic.l.rad
        _y0 = c.galactic.b.rad
        print(str(np.rad2deg(_x0)) + "," + str(np.rad2deg(_y0)), end=",")
        _v0 = np.linalg.norm(_v) / astropy.constants.c.value
        _x, _y = _vector_field.get_grid()
        len_cos = np.sin(np.deg2rad(_y)) * np.sin(_y0) + np.cos(np.deg2rad(_y)) * np.cos(_y0) * np.cos(
            np.deg2rad(_x) - _x0)
        len_c = np.arccos(len_cos)
        ang_b_sin = np.sin(_x0 - np.deg2rad(_x)) * np.cos(np.deg2rad(_y)) / np.sin(len_c)
        ang_b_cos = (np.sin(_y0) - np.sin(np.deg2rad(_y)) * len_cos) / np.sin(len_c) / np.cos(np.deg2rad(_y))
        _shift = _v0 * np.sin(len_c) * (1 + 2.25 * _v0 * _v0) - 0.25 * _v0 * _v0 * np.sin(2 * len_c) \
            + _v0 * _v0 * _v0 * np.sin(3 * len_c) / 12.0
        _uu = _shift * ang_b_sin * 180 / math.pi
        _vv = _shift * ang_b_cos * 180 / math.pi
        return _uu, _vv

    def _figure_title(self):
        fig = plt.figure(tight_layout=True)
        fig.suptitle(self.satellite.title_string())
        return fig

    def figure3(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
        fig2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
        u0, v0 = self.calc_distortion_in_degree(vector_field)
        vector_field.draw(u0, v0, fig1, "aberration field", 15)
        u1 = u0 - np.average(u0)
        v1 = v0 - np.average(v0)
        vector_field.draw(u1, v1, fig2, "difference from mean", 0.1)
        plt.draw()
        plt.savefig('fig3.png')
        plt.show()

    def figure_3a(self, vector_field, file_name):
        fig = self._figure_title()
        fig1 = fig.add_subplot(1, 3, 1, projection=ccrs.PlateCarree())
        fig2 = fig.add_subplot(1, 3, 2, projection=ccrs.PlateCarree())
        fig3 = fig.add_subplot(1, 3, 3, projection=ccrs.PlateCarree())
        u0, v0 = self.calc_distortion_in_degree(vector_field)
        vector_field.draw(u0, v0, fig1, "aberration field", 15)
        u1 = u0 - np.average(u0)
        v1 = v0 - np.average(v0)
        vector_field.draw(u1, v1, fig2, "difference from mean", 0.1)
        X, Y = vector_field.get_grid()
        u2 = subtract_trend(u1, X)
        v2 = subtract_trend(v1, Y)
        vector_field.draw(u2, v2, fig3, "subtract trend", 0.0003)
        plt.draw()
        plt.savefig(file_name)
        plt.show()

    def figure9_left_down(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(2, 2, 1, projection=ccrs.PlateCarree())
        fig2 = fig.add_subplot(2, 2, 2, projection=ccrs.PlateCarree())
        fig3 = fig.add_subplot(2, 2, 3, projection=ccrs.PlateCarree())
        fig4 = fig.add_subplot(2, 2, 4, projection=ccrs.PlateCarree())
        u0, v0 = self.calc_distortion_in_degree(vector_field)
        vector_field.draw(u0, v0, fig1, "aberration field", 15)
        vector_field.draw(u0 - np.average(u0), v0 - np.average(v0), fig2, "difference from mean", 0.1)
        dl = 1.2 / 3600  # unit degree, /3600 means the unit in arc seconds
        u1, v1 = self.calc_distortion_in_degree(vector_field, _dl=dl)
        vector_field.draw(u1 - u0, v1 - v0, fig3, "δl=" + str(dl * 3600) + "arc sec", 2e-5)  # error of direction
        dv = 0.022  # m/s
        u1, v1 = self.calc_distortion_in_degree(vector_field, _dv=dv)
        vector_field.draw(u1 - u0, v1 - v0, fig4, "d|v|=" + str(dv) + "m/s", 2e-5)  # error of speed
        plt.draw()
        plt.savefig('fig9ld.png')
        plt.show()

    def figure7(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
        fig2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
        u0, v0 = self.calc_distortion_in_degree(vector_field)
        dl = 160 / 3600  # unit degree, /3600 means the unit in arc seconds
        u1, v1 = self.calc_distortion_in_degree(vector_field, _dl=dl)
        vector_field.draw(u1 - u0 - np.average(u1) + np.average(u0), v1 - v0 - np.average(v1) + np.average(v0), fig1,
                          "δl=" + str(dl * 3600) + "arc sec", 2e-5)  # error of direction
        dv = 15  # m/s
        u1, v1 = self.calc_distortion_in_degree(vector_field, _dv=dv)
        vector_field.draw(u1 - u0 - np.average(u1) + np.average(u0), v1 - v0 - np.average(v1) + np.average(v0), fig2,
                          "d|v|=" + str(dv) + "m/s", 2e-5)  # error of speed
        plt.draw()
        plt.savefig('fig7.png')
        plt.show()

    def figure5(self, vector_field):
        fig = self._figure_title()
        fig1 = fig.add_subplot(2, 2, 1, projection=ccrs.PlateCarree())
        fig2 = fig.add_subplot(2, 2, 2, projection=ccrs.PlateCarree())
        fig5 = fig.add_subplot(2, 2, 3, projection=ccrs.PlateCarree())
        fig6 = fig.add_subplot(2, 2, 4, projection=ccrs.PlateCarree())
        t0 = self.satellite.obs_date
        um, vm = self.calc_distortion_in_degree(vector_field)
        vector_field.draw(um, vm, fig1, "u(tm): tm=(t1+t2)/2", 15.0)  # aberration field 1. itself as top left
        _phase = self.satellite.phase
        dt = timedelta(seconds=6.25)
        us, vs = self.calc_distortion_in_degree(vector_field, _obs_time=t0 - dt)
        ue, ve = self.calc_distortion_in_degree(vector_field, _obs_time=t0 + dt)
        vector_field.draw(ue - us, ve - vs, fig2, "u(t2) - u(t1)", 0.05)
        # vector_field.draw(ue - np.average(um), ve - np.average(vm), fig3, "v(t2) - <v(tm)>", 0.1)
        # vector_field.draw(ue - um, ve - vm, fig4, "v(t2) - v(tm)", 0.005)
        vector_field.draw((ue + us) * 0.5 - um, (ve + vs) * 0.5 - vm, fig5, "<u> - u(tm). 12.5sec", 0.00015)
        dt = timedelta(seconds=2.499)
        us, vs = self.calc_distortion_in_degree(vector_field, _obs_time=t0 - dt)
        ue, ve = self.calc_distortion_in_degree(vector_field, _obs_time=t0 + dt)
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
        um, vm = self.calc_distortion_in_degree(vector_field)
        vector_field.draw(um - np.average(um), vm - np.average(vm), fig1, "u(tm): tm=(t1+t2)/2", 0.1)
        _phase = self.satellite.phase
        dt = timedelta(seconds=6.25)
        us, vs = self.calc_distortion_in_degree(vector_field, _obs_time=t0 - dt)
        ue, ve = self.calc_distortion_in_degree(vector_field, _obs_time=t0 + dt)
        vector_field.draw(ue - us - np.average(ue - us), ve - vs - np.average(ve - vs), fig2, "u(t2) - u(t1)", 0.00035)
        vector_field.draw((ue + us) * 0.5 - um - np.average((ue + us) * 0.5 - um),
                          (ve + vs) * 0.5 - vm - np.average((ve + vs) * 0.5 - vm), fig5, "<u> - u(tm). 12.5sec", 2e-7)
        dt = timedelta(seconds=60)
        us, vs = self.calc_distortion_in_degree(vector_field, _obs_time=t0 - dt)
        ue, ve = self.calc_distortion_in_degree(vector_field, _obs_time=t0 + dt)
        vector_field.draw((ue + us) * 0.5 - um - np.average((ue + us) * 0.5 - um),
                          (ve + vs) * 0.5 - vm - np.average((ve + vs) * 0.5 - vm), fig6, "<u> - u(tm). 120sec", 2e-5)
        plt.draw()
        plt.savefig('fig10.png')
        plt.show()
