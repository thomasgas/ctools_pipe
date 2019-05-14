from astropy.table import Table
import sys
import yaml


class IRFPicker:
    def __init__(self, irf_container):
        """
        Load informations. Dictionary created from yaml file: see example.
        :param irf_container: dictionary to choose IRF name
        """
        self.irf_info = irf_container['IRF']
        self.zenith = str(self.irf_info['zenith'])
        self.prod_info = self.irf_info['prod']
        self.prod_version = self.prod_info['version']
        self.prod_number = self.prod_info['number']
        self.irf_site = self.irf_info['site']
        self.irf_time = self.irf_info['time']
        self.irf_pointing = self.irf_info['pointing']
        self.irf_subarray = self.irf_info['subarray']
        self.TS = str(self.irf_info['TS'])

    def irf_pick(self):
        """
        Select between the supported IRFs.
        - prod3b-v1 and prod3b-v2 are the names of those provided in the public cta-performance page.
        - prod3b can be downloaded from the redmine page. There are many more IRFs there.
        :return: the irf fits table
        """
        if self.prod_number == "3b" and self.prod_version == 2:
            irf_table = Table.read("irf_fits/prod3b-v2.fits")
            irf = irf_table[
                (irf_table['site'] == self.irf_site) &
                (irf_table['zenith'] == self.zenith) &
                (irf_table['pointing'] == self.irf_pointing) &
                (irf_table['time'] == self.irf_time)
            ]['name']

        elif self.prod_number == "3b" and self.prod_version == 1:
            irf_table = Table.read("irf_fits/prod3b-v1.fits")
            irf = irf_table[
                (irf_table['site'] == self.irf_site) &
                (irf_table['pointing'] == self.irf_pointing) &
                (irf_table['zenith'] == self.zenith) &
                (irf_table['time'] == self.irf_time)
            ]['name']

        elif self.prod_number == "3b" and self.prod_version == 0:
            irf_table = Table.read("irf_fits/prod3b.fits")
            irf = irf_table[
                (irf_table['site'] == self.irf_site) &
                (irf_table['pointing'] == self.irf_pointing) &
                (irf_table['zenith'] == self.zenith) &
                (irf_table['time'] == self.irf_time) &
                (irf_table['subarray'] == self.irf_subarray) &
                (irf_table['TS'] == self.TS)
            ]['name']
        else:
            irf = None
            print("Cannot open any prod fits file. Try with other numbers.\n")

        try:
            return irf[0]
        except IndexError:
            print("No IRF available with those inputs.")
            print("------------------------")
            print(yaml.dump(self.irf_info))
            print("------------------------")
            print("Check IRF availability")
            sys.exit()





