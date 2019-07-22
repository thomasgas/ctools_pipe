import sys
import gammalib
import ctools
import cscripts
import yaml
from utils import create_path
from irf_handler import IRFPicker
import astropy.units as u


def grb_simulation(sim_in, config_in, model_xml, background_fits, counter):
    """

    :param sim_in:
    :param config_in:
    :param model_xml:
    :param background_fits:
    :param counter:
    :return:
    """

    print(model_xml.split('/')[-1], counter)

    ctools_pipe_path = create_path(config_in['exe']['software_path'])

    ctobss_params = sim_in['ctobssim']

    # find proper IRF name
    irf = IRFPicker(sim_in, ctools_pipe_path)
    name_irf = irf.irf_pick()

    if irf.prod_number == "3b" and irf.prod_version == 0:
        caldb = "prod3b"
    else:
        caldb = f'prod{irf.prod_number}-v{irf.prod_version}'

    background_id = f"{str(int(counter) + 1).zfill(6)}"

    obs_back = gammalib.GCTAObservation(background_fits)

    seed = int(counter)*10

    # source simulation
    sim = ctools.ctobssim()
    sim['inmodel'] = model_xml
    sim['caldb'] = caldb
    sim['irf'] = name_irf
    sim['ra'] = 0.5
    sim['dec'] = 0.5
    sim['rad'] = ctobss_params['radius']
    sim['tmin'] = u.Quantity(ctobss_params['time']['t_min']).to_value(u.s)
    sim['tmax'] = u.Quantity(ctobss_params['time']['t_max']).to_value(u.s)
    sim['emin'] = u.Quantity(ctobss_params['energy']['e_min']).to_value(u.TeV)
    sim['emax'] = u.Quantity(ctobss_params['energy']['e_max']).to_value(u.TeV)
    sim['seed'] = seed
    sim.run()

    obs = sim.obs().copy()

    # append all background events to GRB ones ==> there's just one observation and not two
    for event in obs_back.events():
        obs[0].events().append(event)

    # delete all 70+ models from the obs def file...not needed any more
    obs.models(gammalib.GModels())

    select = ctools.ctselect(obs)
    select['rad'] = 5
    select['tmin'] = 0
    select['tmax'] = 100
    select['emin'] = 0.05
    select['emax'] = 1.0
    select.run()

    onoff_sim = cscripts.csphagen(select.obs())
    onoff_sim['inmodel'] = 'NONE'
    onoff_sim['ebinalg'] = 'LOG'
    onoff_sim['emin'] = 0.05
    onoff_sim['emax'] = 1.0
    onoff_sim['enumbins'] = 30
    onoff_sim['coordsys'] = 'CEL'
    onoff_sim['ra'] = 0.0
    onoff_sim['dec'] = 0.0
    onoff_sim['rad'] = 0.2
    onoff_sim['bkgmethod'] = 'REFLECTED'
    onoff_sim['use_model_bkg'] = False
    onoff_sim['stack'] = False
    onoff_sim.run()

    # CSPHAGEN creates a "Dummy" source: here I create a mock model of the GRB to be used in ctlike
    onoff_sim.obs().models()[0].name('GRB')
    expplaw = gammalib.GModelSpectralExpPlaw()
    expplaw['Prefactor'].value(onoff_sim.obs().models()[0]['Prefactor'].value())
    expplaw['Index'].value(onoff_sim.obs().models()[0]['Index'].value())
    expplaw['PivotEnergy'].value(onoff_sim.obs().models()[0]['PivotEnergy'].value())
    expplaw['CutoffEnergy'].value(1.e6)
    onoff_sim.obs().models()[0].spectral(expplaw)

    like = ctools.ctlike(onoff_sim.obs())
    like.run()

    # TODO: extract likelihood informations


if __name__ == '__main__':

    sim_yaml_file = yaml.safe_load(open(sys.argv[1]))
    jobs_yaml_file = yaml.safe_load(open(sys.argv[2]))

    if sim_yaml_file['source']['type'] == "GRB":
        grb_simulation(sim_yaml_file, jobs_yaml_file, sys.argv[3], sys.argv[4], sys.argv[5])
    elif sim_yaml_file['source']['type'] == "GW":
        gw_simulation(sim_yaml_file, jobs_yaml_file, sys.argv[3], sys.argv[4], sys.argv[5])

