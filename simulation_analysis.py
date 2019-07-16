import sys
import ctools
import gammalib
import yaml
from utils import create_path
from irf_handler import IRFPicker
import astropy.units as u


def sim_select_like(sim_yaml, jobs_in, model_xml, background_fits, counter):
    """

    :param sim_yaml:
    :param jobs_in:
    :param model_xml:
    :param background_fits:
    :param counter:
    :return:
    """
    #print("----------------------------")
    print(model_xml.split('/')[-1], counter)

    config_in = yaml.safe_load(open(jobs_in))
    ctools_pipe_path = create_path(config_in['exe']['software_path'])

    sim_in = yaml.safe_load(open(sim_yaml))
    ctobss_params = sim_in['ctobssim']

    # find proper IRF name
    irf = IRFPicker(sim_in, ctools_pipe_path)
    name_irf = irf.irf_pick()

    if irf.prod_version == "3b" and irf.prod_number == 0:
        caldb = "prod3b"
    else:
        caldb = f'prod{irf.prod_number}-v{irf.prod_version}'

    # loading background (this is a way of doing it without saving any file)
    # output.save("name.xml") to save the file
    obs_def = gammalib.GObservations()

    background_id = f"{str(int(counter) + 1).zfill(6)}"

    output = gammalib.GXml()
    level0 = gammalib.GXmlElement('observation_list title="observation library"')
    level1 = gammalib.GXmlElement(f'observation name="name_source" id="{background_id}" instrument="CTA"')
    level2 = gammalib.GXmlElement(f'parameter name="EventList"  file="{background_fits}"')
    level1.append(level2)
    level0.append(level1)
    output.append(level0)

    obs_def.read(output)

    seed = int(counter)*10

    # do the simulation
    sim = ctools.ctobssim()
    sim['inmodel'] = model_xml
    sim['caldb'] = caldb
    sim['irf'] = name_irf
    sim['ra'] = 0
    sim['dec'] = 0
    sim['rad'] = ctobss_params['radius']
    sim['tmin'] = u.Quantity(ctobss_params['time']['t_min']).to_value(u.s)
    sim['tmax'] = u.Quantity(ctobss_params['time']['t_max']).to_value(u.s)
    sim['emin'] = u.Quantity(ctobss_params['energy']['e_min']).to_value(u.TeV)
    sim['emax'] = u.Quantity(ctobss_params['energy']['e_max']).to_value(u.TeV)
    sim['seed'] = seed
    sim.run()

    obs_def.append(sim.obs()[0])

    # for obs in obs_def:
    #     print(obs.id(), obs.nobserved())

    select = ctools.ctselect(obs_def)
    select['rad'] = 3
    select['tmin'] = 0
    select['tmax'] = 50
    select['emin'] = 0.05
    select['emax'] = 1.0
    select.run()

    # print(select.obs().models())

if __name__ == '__main__':

    sim_select_like(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
