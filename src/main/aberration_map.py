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
    # _fig1.gridlines(draw_labels=False)
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


def aberration_statistics(_dx):
    a = VectorField2d(_dx * 30, _dx * 30, _dx)
    tz = astropy.time.TimezoneInfo(0 * u.hour)
    r = RelativisticAberration(Satellite(datetime(2024, 8, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(120)))
    fig = r._figure_title()
    fig1 = fig.add_subplot(2, 2, 1, projection=ccrs.PlateCarree())
    fig2 = fig.add_subplot(2, 2, 2, projection='3d')
    fig3 = fig.add_subplot(2, 2, 3)
    fig4 = fig.add_subplot(2, 2, 4)
    X, Y = a.get_grid()
    u0, v0 = r.calc_distortion_in_degree(a)
    a.draw(u0 - np.average(u0), v0 - np.average(v0), fig1, "difference from mean", 0.1)
    mag = np.sqrt((u0 - np.average(u0)) ** 2 + (v0 - np.average(v0)) ** 2) * 3600
    fig2.plot_surface(X, Y, mag)
    fig2.set_title("relative aberration")
    fig2.set_xlabel("lon(deg)")
    fig2.set_ylabel("lat(deg)")
    fig2.set_zlabel("|aberration|(arcsec)")
    R = np.sqrt(X * X + Y * Y)
    fig3.scatter(np.ravel(R), np.ravel(mag), s=1)
    fig3.set_title("distance vs relative aberration")
    fig3.set_xlabel("distance")
    fig3.set_ylabel("|aberration|")
    angx = np.rad2deg(np.arctan2(X, Y))
    angm = np.rad2deg(np.arctan2(- u0 + np.average(u0), - v0 + np.average(v0)))
    fig4.scatter(angx, angm, s=1)
    fig4.set_title("angle of coord vs aberration")
    fig4.set_xlim(-180, 180)
    fig4.set_ylim(-180, 180)
    fig4.set_xticks([-180, -90, 0, 90, 180])
    fig4.set_yticks([-180, -90, 0, 90, 180])
    fig4.set_xlabel("position")
    fig4.set_ylabel("aberration")
    plt.draw()
    plt.savefig('fig3.png')
    plt.show()


def normal_figures():
    a = VectorField2d(0.3, 0.3, 0.1)
    tz = astropy.time.TimezoneInfo(0 * u.hour)
    r = RelativisticAberration(Satellite(datetime(2024, 2, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(120)))
    r.figure3(a)
    r.figure9_left_down(a)
    r.figure5(a)
    r.figure10(a)
    r = RelativisticAberration(Satellite(datetime(2024, 3, 17, 12, 0, 0, tzinfo=tz), np.deg2rad(110)))
    r.figure7(a)


def sub_trend():
    a = VectorField2d(0.3, 0.3, 0.1)
    tz = astropy.time.TimezoneInfo(0 * u.hour)
    # r = RelativisticAberration(Satellite(datetime(2024, 2, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(120)))
    # r.figure_3a(a, "fig3a1.png")
    #r = RelativisticAberration(Satellite(datetime(2024, 3, 17, 12, 0, 0, tzinfo=tz), np.deg2rad(120)))
    #r.figure_3a(a, "fig3a2.png")
    r = RelativisticAberration(Satellite(datetime(2024, 5, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(0)))
    r.figure_3b(a, "fig3a3-000.png")
    r = RelativisticAberration(Satellite(datetime(2024, 5, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(30)))
    r.figure_3b(a, "fig3a3-030.png")
    r = RelativisticAberration(Satellite(datetime(2024, 5, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(60)))
    r.figure_3b(a, "fig3a3-060.png")
    r = RelativisticAberration(Satellite(datetime(2024, 5, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(90)))
    r.figure_3b(a, "fig3a3-090.png")
    r = RelativisticAberration(Satellite(datetime(2024, 5, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(120)))
    r.figure_3b(a, "fig3a3-120.png")
    r = RelativisticAberration(Satellite(datetime(2024, 5, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(150)))
    r.figure_3b(a, "fig3a3-150.png")
    r = RelativisticAberration(Satellite(datetime(2024, 5, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(180)))
    r.figure_3b(a, "fig3a3-180.png")
    #r = RelativisticAberration(Satellite(datetime(2024, 8, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(120)))
    #r.figure_3a(a, "fig3a4.png")
    #r = RelativisticAberration(Satellite(datetime(2024, 9, 20, 12, 0, 0, tzinfo=tz), np.deg2rad(120)))
    #r.figure_3a(a, "fig3a5.png")
    #r = RelativisticAberration(Satellite(datetime(2024, 11, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(120)))
    #r.figure_3a(a, "fig3a6.png")


def sub_trend_c():
    a = VectorField2d(0.3, 0.3, 0.1)
    tz = astropy.time.TimezoneInfo(0 * u.hour)
    r = RelativisticAberration(Satellite(datetime(2024, 2, 20, 12, 0, 0, tzinfo=tz), np.deg2rad(120)))
    r.figure_3b(a, "fig3b3.png")


def change_during_exposure():
    a = VectorField2d(0.3, 0.3, 0.1)
    tz = astropy.time.TimezoneInfo(0 * u.hour)
    r = RelativisticAberration(Satellite(datetime(2024, 2, 1, 12, 0, 0, tzinfo=tz), np.deg2rad(120)))
    r.figure11(a)


def trend(_u, _x):
    u_0 = _u[0, 0]
    u_1 = _u[-1, -1]
    x0 = _x[0, 0]
    x1 = _x[-1, -1]
    _a = (u_1 - u_0) / (x1 - x0)
    _b = (u_0 * x1 - u_1 * x0) / (x1 - x0)
    return _a * _x + _b

# global_map(0, 0, 0.0001)
# aberration_statistics(0.01)
normal_figures()
# sub_trend()
# change_during_exposure()
