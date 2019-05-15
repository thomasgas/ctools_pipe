import gammalib
import ctools
import cscripts
import sys
import os
import yaml
from irf_handler import IRFPicker
import numpy as np
import astropy.units as u
from environs import Env

def simulate_background(input_yaml, count):
    config_in = yaml.safe_load(open(input_yaml))

    # find proper IRF name
    irf = IRFPicker(config_in)
    name_irf = irf.irf_pick()

    if irf.prod_version == "3b" and irf.prod_number == 0:
        caldb = "prod3b"
    else:
        caldb = f'prod{irf.prod_number}-v{irf.prod_version}'

    # output file
    env = Env()
    out_path = config_in['exe']['path']

    # try just to interpret the env variable (removing the dollar)
    # if it's $VAR/folder then split 
    if out_path.startswith("$"):
        out_path = out_path[1:] 
        try:
            out_path = env(out_path)
        except:
            env_var = out_path.split('/', 1)[0]
            folders = out_path.split('/', 1)[1]
            out_path = env(env_var) + '/' + folders

    if out_path.endswith("/"):   # just used to remove backslash if applied
        out_path = out_path[:-1]

    print("out_path")
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    
    # simulation details
    sim_details = config_in['sim']

    seed = np.random.randint(10000)

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
