from matplotlib import pyplot as plt
import numpy as np
import astropy.time
import astropy.units as u
import cartopy.crs as ccrs
from astropy.coordinates import SkyCoord
from datetime import datetime, timedelta
from util.satellite import Satellite, vector_to_sky_coord


def observation_region():
    coord_galactic = np.array([[-1.4, -0.6], [-1.4, 0.6], [0.7, 0.6], [0.7, -0.6], [-1.4, -0.6]])
    coord_ecliptic = np.empty(coord_galactic.shape)
    for _i in range(coord_galactic.shape[0]):
        c0 = SkyCoord(coord_galactic[_i][0], coord_galactic[_i][1], unit="deg", frame="galactic")
        coord_ecliptic[_i] = [c0.barycentricmeanecliptic.lon.deg, c0.barycentricmeanecliptic.lat.deg]
    return coord_ecliptic


def galactic_coordinate(ax):
    coord_galactic_b = np.array([[0.0, 0.0], [0.0, 0.8]])
    coord_galactic_l = np.array([[0.0, 0.0], [0.8, 0.0]])
    coord_ecliptic_b = np.empty(coord_galactic_b.shape)
    coord_ecliptic_l = np.empty(coord_galactic_l.shape)
    for _i in range(coord_galactic_l.shape[0]):
        c0 = SkyCoord(coord_galactic_l[_i][0], coord_galactic_l[_i][1], unit="deg", frame="galactic")
        coord_ecliptic_l[_i] = [c0.barycentricmeanecliptic.lon.deg, c0.barycentricmeanecliptic.lat.deg]
    for _i in range(coord_galactic_b.shape[0]):
        c0 = SkyCoord(coord_galactic_b[_i][0], coord_galactic_b[_i][1], unit="deg", frame="galactic")
        coord_ecliptic_b[_i] = [c0.barycentricmeanecliptic.lon.deg, c0.barycentricmeanecliptic.lat.deg]
    x = (coord_ecliptic_l[0][0], coord_ecliptic_l[1][0])
    y = (coord_ecliptic_l[0][1], coord_ecliptic_l[1][1])
    ax.plot(x, y, color="red", transform=ccrs.PlateCarree())
    x = (coord_ecliptic_b[0][0], coord_ecliptic_b[1][0])
    y = (coord_ecliptic_b[0][1], coord_ecliptic_b[1][1])
    ax.plot(x, y, color="red", transform=ccrs.PlateCarree())


def base_map():
    _fig = plt.figure()
    _fig1 = _fig.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
    _fig1.set_global()
    _fig1.gridlines(draw_labels=True)
    # _x = 10 * np.arange(36 + 1) - 180.0
    # _y = 90. - 10 * np.arange(18 + 1)
    # _X, _Y = np.meshgrid(_x, _y)
    _fig1.invert_xaxis()
    return _fig1


def velocity_plot(_ax, _t_start):
    dt = timedelta(minutes=5)
    for j in range(len(_t_start)):
        sat = Satellite(_t_start[j], 0.0)
        N = 10
        lon = np.empty(N)
        lat = np.empty(N)
        for i in range(N):
            t = _t_start[j] + dt * i
            v = sat.calc_satellite_velocity(_obs_time=t)
            c = vector_to_sky_coord(v)
            lon[i] = c.barycentricmeanecliptic.lon.deg
            lat[i] = c.barycentricmeanecliptic.lat.deg
        _ax.scatter(lon, lat, transform=ccrs.PlateCarree())


fig1 = base_map()
ce = observation_region()
fig1.plot(ce[:, 0], ce[:, 1], color="green", transform=ccrs.PlateCarree())
galactic_coordinate(fig1)
t_start = [datetime(2024, 2, 1, 12, 0, 0, tzinfo=astropy.time.TimezoneInfo(0 * u.hour)),
           datetime(2024, 3, 20, 12, 0, 0, tzinfo=astropy.time.TimezoneInfo(0 * u.hour)),
           datetime(2024, 5, 1, 12, 0, 0, tzinfo=astropy.time.TimezoneInfo(0 * u.hour)),
           datetime(2024, 8, 1, 12, 0, 0, tzinfo=astropy.time.TimezoneInfo(0 * u.hour)),
           datetime(2024, 9, 20, 12, 0, 0, tzinfo=astropy.time.TimezoneInfo(0 * u.hour)),
           datetime(2024, 11, 1, 12, 0, 0, tzinfo=astropy.time.TimezoneInfo(0 * u.hour))]
velocity_plot(fig1, t_start)
plt.draw()
plt.savefig("velmap.png")
plt.show()
