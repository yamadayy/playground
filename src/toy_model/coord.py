from astropy.coordinates import SkyCoord

if __name__ == '__main__':
    c = SkyCoord(0, 0, unit="deg", frame="galactic")
    print(c.galactic)
    print(c.icrs)
    print(c.barycentricmeanecliptic)
