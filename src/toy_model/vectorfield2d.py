from datetime import datetime
import astropy.time
import astropy.constants
import cartopy.crs as ccrs
import numpy as np
import astropy.units as u
from matplotlib import pyplot as plt
from util.vector_field_2d import VectorField2d
from util.satellite import Satellite
from util.relativistic_aberration import RelativisticAberration


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
r = RelativisticAberration(Satellite(datetime(2024, 2, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(120)))
r.figure3(a)
r.figure9_left_down(a)
r.figure5(a)
r.figure10(a)
r = RelativisticAberration(Satellite(datetime(2024, 3, 17, 12, 0, 0, tzinfo=tz), np.deg2rad(110)))
r.figure7(a)
