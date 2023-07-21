import math

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time

from jasmine_toolkit.datamodel.catalogue_entry import CatalogueEntry


def test_ecliptic_to_eplictic():
    coord = SkyCoord(0.0 * u.deg, 0.0 * u.deg, distance=1.0 * u.kpc, pm_lat=1.0 * u.mas / u.yr,
                     pm_lon_coslat=1.0 * u.mas / u.yr, frame='barycentricmeanecliptic')
    c = CatalogueEntry(1, coord, 12.5, Time('2022-01-01T00:00:00'))
    assert math.isclose(c.lat.to('deg').value, 0.0, abs_tol=1.0e-3)
    assert math.isclose(c.lon.to('deg').value, 0.0, abs_tol=1.0e-3)
    assert math.isclose(c.pm_lat.value, 1.0, abs_tol=1.0e-3)
    assert math.isclose(c.pm_lon_coslat.value, 1.0, abs_tol=1.0e-3)
    assert math.isclose(c.parallax.value, 1.0, abs_tol=1.0e-3)


def test_icrs_to_ecliptic():
    coord = SkyCoord(90.0 * u.deg, 0.0 * u.deg, pm_dec=1.0 * u.mas / u.yr, pm_ra_cosdec=1.0 * u.mas / u.yr,
                     distance=1.0 * u.kpc, frame='icrs')
    c = CatalogueEntry(1, coord, 12.5, Time('2022-01-01T00:00:00'))
    assert math.isclose(c.lat.to('deg').value, -23.439, abs_tol=1.0e-3)
    assert math.isclose(c.lon.to('deg').value, 90.0, abs_tol=1.0e-3)
    assert math.isclose(c.pm_lat.value, 1.0, abs_tol=1.0e-3)
    assert math.isclose(c.pm_lon_coslat.value, 1.0, abs_tol=1.0e-3)
    assert math.isclose(c.parallax.value, 1.0, abs_tol=1.0e-3)


def test_galactic_to_ecliptic():
    coord = SkyCoord(0.0 * u.deg, 0.0 * u.deg, pm_b=1.0 * u.mas / u.yr, pm_l_cosb=1.0 * u.mas / u.yr,
                     distance=1.0 * u.kpc, frame='galactic')
    c = CatalogueEntry(1, coord, 12.5, Time('2022-01-01T00:00:00'))
    assert math.isclose(c.lat.to('deg').value, -5.5363, abs_tol=1.0e-3)
    assert math.isclose(c.lon.to('deg').value, 266.8395, abs_tol=1.0e-3)
    assert math.isclose(c.pm_lat.value, 1.366, abs_tol=1.0e-3)
    assert math.isclose(c.pm_lon_coslat.value, -0.3668, abs_tol=1.0e-3)
    assert math.isclose(c.parallax.value, 1.0, abs_tol=1.0e-3)


def test_ecliptic_to_eplictic_without_pm():
    coord = SkyCoord(0.0 * u.deg, 0.0 * u.deg, frame='barycentricmeanecliptic')
    c = CatalogueEntry(1, coord, 12.5, Time('2022-01-01T00:00:00'))
    assert math.isclose(c.lat.to('deg').value, 0.0, abs_tol=1.0e-3)
    assert math.isclose(c.lon.to('deg').value, 0.0, abs_tol=1.0e-3)
    assert math.isclose(c.pm_lat.value, 0.0, abs_tol=1.0e-3)
    assert math.isclose(c.pm_lon_coslat.value, 0.0, abs_tol=1.0e-3)
    assert math.isclose(c.parallax.value, 0.0, abs_tol=1.0e-3)


def test_icrs_to_ecliptic_without_pm():
    coord = SkyCoord(90.0 * u.deg, 0.0 * u.deg, frame='icrs')
    c = CatalogueEntry(1, coord, 12.5, Time('2022-01-01T00:00:00'))
    assert math.isclose(c.lat.to('deg').value, -23.439, abs_tol=1.0e-3)
    assert math.isclose(c.lon.to('deg').value, 90.0, abs_tol=1.0e-3)
    assert math.isclose(c.pm_lat.value, 0.0, abs_tol=1.0e-3)
    assert math.isclose(c.pm_lon_coslat.value, 0.0, abs_tol=1.0e-3)
    assert math.isclose(c.parallax.value, 0.0, abs_tol=1.0e-3)
