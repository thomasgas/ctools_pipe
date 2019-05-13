import gammalib
import ctools
import cscripts
import sys
import os
import yaml

from astropy.io.misc import yaml as yl_as
import astropy.units as u


def simulate_background(input_yaml):
    config_in = yaml.safe_load(open(input_yaml))
    print(yaml.dump(config_in))

    irf_info = config_in['IRF']

    # zenith = irf_info['zenith']
    # prod_info = irf_info['prod']
    # prod_version = prod_info['version']
    # prod_number = prod_info['number']
    # irf_site = irf_info['site']
    # irf_time = u.Quantity(irf_info['time'])
    # irf_pointing = irf_info['pointing']
    # irf_subarray = irf_info['subarray']

    # parse the IRF input to get a list of IRFs compliant
    # with the one present in the ctools folders

    # sim = ctools.ctobssim()
    # sim['inmodel'] =
    # sim['caldb'] = f'prod{prod_number}-v{prod_version}'
    # sim['irf'] =
    # sim['ra'] =
    # sim['dec'] =
    # sim['rad'] =
    # sim['tmin'] =
    # sim['tmax'] =
    # sim['emin'] =
    # sim['emax'] =
    # sim['outevents'] =
    # sim['debug'] =
    # sim['seed'] =
    # sim.execute()


