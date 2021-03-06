import ctools
import sys
import yaml
from irf_handler import IRFPicker
import astropy.units as u

from utils import create_path


def simulate_background(input_yaml, jobs_yaml, count):
    config_in = yaml.safe_load(open(input_yaml))
    jobs_config = yaml.safe_load(open(jobs_yaml))

    try:
        ctools_pipe_path = create_path(jobs_config['exe']['software_path'])
    except KeyError:
        ctools_pipe_path = "."

    # find proper IRF name
    irf = IRFPicker(config_in, ctools_pipe_path)
    name_irf = irf.irf_pick()

    if irf.prod_number == "3b" and irf.prod_version == 0:
        caldb = "prod3b"
    else:
        caldb = f'prod{irf.prod_number}-v{irf.prod_version}'

    out_path = create_path(f"{jobs_config['exe']['path']}/back_sim/{irf.prod_number}_{irf.prod_version}_{name_irf}")

    # simulation details
    sim_details = config_in['sim']

    seed = int(count) * 10

    # do the simulation
    sim = ctools.ctobssim()
    sim['inmodel'] = f"{ctools_pipe_path}/models/bkg_only_model.xml"
    sim['caldb'] = caldb
    sim['irf'] = name_irf
    sim['ra'] = 0
    sim['dec'] = 0
    sim['rad'] = sim_details['radius']
    sim['tmin'] = u.Quantity(sim_details['time']['t_min']).to_value(u.s)
    sim['tmax'] = u.Quantity(sim_details['time']['t_max']).to_value(u.s)
    sim['emin'] = u.Quantity(sim_details['energy']['e_min']).to_value(u.TeV)
    sim['emax'] = u.Quantity(sim_details['energy']['e_max']).to_value(u.TeV)
    sim['outevents'] = f"{out_path}/background_z-{irf.zenith}_site-{irf.irf_site}_{str(count).zfill(2)}_seed{seed}.fits"
    sim['seed'] = seed
    sim.execute()


if __name__ == '__main__':
    simulate_background(sys.argv[1], sys.argv[2], sys.argv[3])
