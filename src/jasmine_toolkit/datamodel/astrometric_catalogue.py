import numpy as np

from jasmine_toolkit.datamodel.catalogue_entry import CatalogueEntry


class AstrometricCatalogue:
    def __init__(self):
        self.__catalogue_entries = np.empty(0)

    def add_entry(self, entry: CatalogueEntry):
        pass

    def get_entry(self, stellar_id: int):
        pass

    def get_catalogue(self):
        return self.__catalogue_entries


if __name__ == '__main__':
    c = AstrometricCatalogue()
    print(c.get_catalogue())
