import sys
import ctools
import gammalib
import yaml
from utils import create_path
from irf_handler import IRFPicker


def sim_select_like(sim_in, jobs_in, model_xml, counter):
    """

    :param sim_in:
    :param jobs_in:
    :param model_xml:
    :param counter:
    :return:
    """
    print("----------------------------")
    print(model_xml, counter)

    config_in = yaml.safe_load(open(jobs_in))
    try:
        ctools_pipe_path = create_path(config_in['exe']['software_path'])
    except KeyError:
        ctools_pipe_path = "."

    print(ctools_pipe_path)

    # # find proper IRF name
    # irf = IRFPicker(config_in, ctools_pipe_path)
    # name_irf = irf.irf_pick()
    #
    # if irf.prod_version == "3b" and irf.prod_number == 0:
    #     caldb = "prod3b"
    # else:
    #     caldb = f'prod{irf.prod_number}-v{irf.prod_version}'
    #
    # obs_bkg = gammalib.GObservations(background_file)
    #
    # seed = int(count)*10
    # sim_details = config_in['sim']
    #
    # # do the simulation
    # sim = ctools.ctobssim()
    # sim['inmodel'] = f"{ctools_pipe_path}/models/bkg_only_model.xml"
    # sim['caldb'] = caldb
    # sim['irf'] = name_irf
    # sim['ra'] = 0
    # sim['dec'] = 0
    # sim['rad'] = sim_details['radius']
    # sim['tmin'] = u.Quantity(sim_details['time']['t_min']).to_value(u.s)
    # sim['tmax'] = u.Quantity(sim_details['time']['t_max']).to_value(u.s)
    # sim['emin'] = u.Quantity(sim_details['energy']['e_min']).to_value(u.TeV)
    # sim['emax'] = u.Quantity(sim_details['energy']['e_max']).to_value(u.TeV)
    # sim['outevents'] = f"{out_path}/background_z-{irf.zenith}_site-{irf.irf_site}_{str(count).zfill(2)}_seed{seed}.fits"
    # sim['seed'] = seed
    # sim.execute()


if __name__ == '__main__':

    sim_select_like(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
