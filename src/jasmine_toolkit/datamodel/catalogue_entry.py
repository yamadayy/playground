from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u


class CatalogueEntry:
    def __init__(self, stellar_id: int, coordinate: SkyCoord, magnitude: float = 0.0,
                 datetime: Time = Time('2022-01-01T00:00:00')):
        self.__stellar_id = stellar_id
        self.__coordinate = coordinate
        self.__fix_pm()
        self.__magnitude = magnitude
        self.__datetime = datetime

    def __fix_pm(self):
        try:
            self.__coordinate.barycentricmeanecliptic.pm_lat
        except TypeError:
            tmp = SkyCoord(self.__coordinate.barycentricmeanecliptic.lon, self.__coordinate.barycentricmeanecliptic.lat,
                           self.__coordinate.distance, pm_lat=0.0 * u.mas / u.yr, pm_lon_coslat=0.0 * u.mas / u.yr,
                           frame='barycentricmeanecliptic')
            self.__coordinate = tmp
        else:
            self.__coordinate = self.__coordinate.barycentricmeanecliptic

    @property
    def stellar_id(self):
        return self.__stellar_id

    @property
    def lon(self):
        return self.__coordinate.lon

    @property
    def lat(self):
        return self.__coordinate.lat

    @property
    def magnitude(self):
        return self.__magnitude

    @property
    def datetime(self):
        return self.__datetime

    @property
    def pm_lon_coslat(self):
        return self.__coordinate.pm_lon_coslat

    @property
    def pm_lat(self):
        return self.__coordinate.pm_lat

    @property
    def parallax(self):
        if self.__coordinate.distance.unit == "":
            return 0.0 * u.mas
        else:
            return 1 / self.__coordinate.distance.to('kpc').value * u.mas
