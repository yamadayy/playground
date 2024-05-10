import numpy as np
from astropy.coordinates import SkyCoord
from matplotlib import pyplot as plt
import cartopy.crs as ccrs


def convert_galactic_to_ecliptic(_array_galactic):
    _array_ecliptic = np.empty(_array_galactic.shape)
    for _i in range(_array_galactic.shape[0]):
        c0 = SkyCoord(_array_galactic[_i][0], _array_galactic[_i][1], unit="deg", frame="galactic")
        _array_ecliptic[_i] = [c0.barycentricmeanecliptic.lon.deg, c0.barycentricmeanecliptic.lat.deg]
    return _array_ecliptic


def observation_region():
    coord_galactic = np.array([[-1.4, -0.6], [-1.4, 0.6], [0.7, 0.6], [0.7, -0.6], [-1.4, -0.6]])
    return convert_galactic_to_ecliptic(coord_galactic)


def detectors_location():
    _focal_length = 4.3704
    _detector_size = 0.01952
    _detector_gap = 0.0033
    _li = np.rad2deg(_detector_gap / 2 / _focal_length) # location inner
    _lo = np.rad2deg((_detector_gap / 2 + _detector_size) / _focal_length) # location outer
    de1 = convert_galactic_to_ecliptic(np.array([[_li, _li], [_lo, _li], [_lo, _lo], [_li, _lo], [_li, _li]]))
    de2 = convert_galactic_to_ecliptic(np.array([[-_li, _li], [-_lo, _li], [-_lo, _lo], [-_li, _lo], [-_li, _li]]))
    de3 = convert_galactic_to_ecliptic(np.array([[-_li, -_li], [-_lo, -_li], [-_lo, -_lo], [-_li, -_lo], [-_li, -_li]]))
    de4 = convert_galactic_to_ecliptic(np.array([[_li, -_li], [_lo, -_li], [_lo, -_lo], [_li, -_lo], [_li, -_li]]))
    return de1, de2, de3, de4


def global_map():
    ce = observation_region()
    fig = plt.figure()
    fig1 = fig.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
    fig1.set_global()
    fig1.invert_xaxis()
    fig1.gridlines(draw_labels=True)
    fig1.plot(ce[:, 0], ce[:, 1], color="green", transform=ccrs.PlateCarree())
    plt.draw()
    # plt.savefig("global_map.jpg")
    plt.show()


def local_map():
    ce = observation_region()
    fig = plt.figure()
    proj = ccrs.PlateCarree()
    fig1 = fig.add_subplot(1, 1, 1, projection=proj)
    fig1.set_extent([265, 268, -7.5, -4.5], proj)
    fig1.invert_xaxis()
    fig1.gridlines(draw_labels=True)
    fig1.plot(ce[:, 0], ce[:, 1], color="green", transform=proj)
    plt.draw()
    # plt.savefig("local_map.jpg")
    plt.show()


def detector():
    ce = observation_region()
    de1, de2, de3, de4 = detectors_location()
    fig = plt.figure()
    proj = ccrs.PlateCarree()
    fig1 = fig.add_subplot(1, 1, 1, projection=proj)
    fig1.set_extent([265, 268, -7.5, -4.5], proj)
    fig1.gridlines(draw_labels=True)
    fig1.plot(ce[:, 0], ce[:, 1], color="green", transform=proj)
    fig1.plot(de1[:, 0], de1[:, 1], color="green", transform=proj)
    fig1.plot(de2[:, 0], de2[:, 1], color="green", transform=proj)
    fig1.plot(de3[:, 0], de3[:, 1], color="green", transform=proj)
    fig1.plot(de4[:, 0], de4[:, 1], color="green", transform=proj)
    plt.draw()
    # plt.savefig("detector_map.jpg")
    plt.show()


if __name__ == '__main__':
    global_map()
    local_map()
    detector()
