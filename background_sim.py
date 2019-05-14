import gammalib
import ctools
import cscripts
import sys
import os
import yaml
from irf_handler import IRFPicker
import numpy as np

from astropy.io.misc import yaml as yl_as
import astropy.units as u


def simulate_background(input_yaml, count, seed):
    config_in = yaml.safe_load(open(input_yaml))

    # find proper IRF name
    irf = IRFPicker(config_in)
    name_irf = irf.irf_pick()

    if irf.prod_version == "3b" and irf.prod_number == 0:
        caldb = "prod3b"
    else:
        caldb = f'prod{irf.prod_number}-v{irf.prod_version}'

    # output file
    out_path = config_in['exe']['path']
    if not out_path.endswith("/"):
        out_path += "/"

    # simulation details
    sim_details = config_in['sim']

    seed = np.random.randint(1000)

    # do the simulation
    sim = ctools.ctobssim()
    sim['inmodel'] = "models/bkg_only_model.xml"
    sim['caldb'] = caldb
    sim['irf'] = name_irf
    sim['ra'] = 0
    sim['dec'] = 0
    sim['rad'] = sim_details['radius']
    sim['tmin'] = u.Quantity(sim_details['time']['t_min']).to_value(u.s)
    sim['tmax'] = u.Quantity(sim_details['time']['t_max']).to_value(u.s)
    sim['emin'] = u.Quantity(sim_details['energy']['e_min']).to_value(u.TeV)
    sim['emax'] = u.Quantity(sim_details['energy']['e_max']).to_value(u.TeV)
    sim['outevents'] = f"{out_path}/background_z-{irf.zenith}_site-{irf.irf_site}_{str(count).zfill(2)}_{seed}.fits"
    sim['seed'] = seed
    sim.execute()


if __name__ == '__main__':
    simulate_background(sys.argv[1], sys.argv[2])
