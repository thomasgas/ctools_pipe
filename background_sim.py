import gammalib
import ctools
import cscripts
import sys
import os
import yaml
from irf_handler import IRFPicker

from astropy.io.misc import yaml as yl_as
import astropy.units as u


def simulate_background(input_yaml):
    config_in = yaml.safe_load(open(input_yaml))
    # print(yaml.dump(config_in))

    irf = IRFPicker(config_in)
    name_irf = irf.irf_pick()
    print(name_irf)

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


if __name__ == '__main__':
    simulate_background(sys.argv[1])
