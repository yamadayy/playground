from astropy.coordinates import SkyCoord

c = SkyCoord(0, 0, unit="deg", frame="galactic")
print(c.galactic)
print(c.icrs)
print(c.barycentricmeanecliptic)
