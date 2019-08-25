import sys
import gammalib
import ctools
import cscripts
import yaml
from utils import create_path
from irf_handler import IRFPicker
import astropy.units as u
from gammapy.stats import significance_on_off
import numpy as np


def grb_simulation(sim_in, config_in, model_xml, counter, background_fits):
    """

    :param sim_in:
    :param config_in:
    :param model_xml:
    :param background_fits:
    :param counter:
    :return:
    """

    grb_name = model_xml.split('/')[-1].split('_')[1][:-4]

    print(grb_name, counter)

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

    sim_t_min = u.Quantity(ctobss_params['time']['t_min']).to_value(u.s)
    sim_t_max = u.Quantity(ctobss_params['time']['t_max']).to_value(u.s)
    sim_e_min = u.Quantity(ctobss_params['energy']['e_min']).to_value(u.TeV)
    sim_e_max = u.Quantity(ctobss_params['energy']['e_max']).to_value(u.TeV)
    sim_rad = ctobss_params['radius']

    # source simulation
    sim = ctools.ctobssim()
    sim['inmodel'] = model_xml
    sim['caldb'] = caldb
    sim['irf'] = name_irf
    sim['ra'] = 0.0
    sim['dec'] = 0.0
    sim['rad'] = sim_rad
    sim['tmin'] = sim_t_min
    sim['tmax'] = sim_t_max
    sim['emin'] = sim_e_min
    sim['emax'] = sim_e_max
    sim['seed'] = seed
    sim.run()

    obs = sim.obs().copy()

    # append all background events to GRB ones ==> there's just one observation and not two
    for event in obs_back.events():
        obs[0].events().append(event)

    # delete all 70+ models from the obs def file...not needed any more
    obs.models(gammalib.GModels())

    obs_pl = obs.copy()

    # CTSELECT
    select_time = sim_in['ctselect']['time_cut']
    slices = select_time['t_slices']
    if slices == 0:
        times = [sim_t_min, sim_t_max]
        times_start = times[:-1]
        times_end = times[1:]
    elif slices > 0:
        time_mode = select_time['mode']
        t_min = select_time['t_min']
        t_max = select_time['t_max']
        if time_mode == "log":
            times = np.logspace(np.log10(t_min), np.log10(t_max), slices + 1, endpoint=True)
        elif time_mode == "lin":
            times = np.linspace(t_min, t_max, slices + 1, endpoint=True)
        else:
            print(f"{time_mode} not valid. Use 'log' or 'lin' ")
            sys.exit()

        if select_time['obs_mode'] == "iter":
            times_start = times[:-1]
            times_end = times[1:]
        elif select_time['obs_mode'] == "cumul":
            times_start = np.repeat(times[0], slices)         # this is to use the same array structure for the loop
            times_end = times[1:]
        elif select_time['obs_mode'] == "all":
            begins, ends = np.meshgrid(times[:-1], times[1:])
            mask_times = begins < ends
            times_start = begins[mask_times].ravel()
            times_end = ends[mask_times].ravel()
        else:
            print(f"obs_mode: {select_time['obs_mode']} not supported")
            sys.exit()

    else:
        print(f"value {slices} not supported...check yaml file")
        sys.exit()

    for t_in, t_end in zip(times_start, times_end):

        # perform time selection for this specific time bin
        select = ctools.ctselect(obs)
        select['rad'] = sim_rad
        select['tmin'] = t_in
        select['tmax'] = t_end
        select['emin'] = sim_e_min
        select['emax'] = sim_e_max
        select.run()

        sel_obs_copy = select.obs().copy()

        # TODO: add option for ON/OFF or normal simulation
        onoff_sim = cscripts.csphagen(select.obs())
        onoff_sim['inmodel'] = 'NONE'
        onoff_sim['ebinalg'] = 'LOG'
        onoff_sim['emin'] = sim_e_min
        onoff_sim['emax'] = sim_e_max
        onoff_sim['enumbins'] = 20
        onoff_sim['coordsys'] = 'CEL'
        onoff_sim['ra'] = 0.0
        onoff_sim['dec'] = 0.5
        onoff_sim['rad'] = 0.2
        onoff_sim['bkgmethod'] = 'REFLECTED'
        onoff_sim['use_model_bkg'] = False
        onoff_sim['stack'] = False
        onoff_sim.run()

        on_counts = onoff_sim.obs()[0].on_spec().counts()
        off_counts = onoff_sim.obs()[0].off_spec().counts()
        # on_scale = onoff_sim.obs()[0].on_spec().backscal_spectrum()
        # off_scale = onoff_sim.obs()[0].off_spec().backscal_spectrum()

        alpha = onoff_sim.obs()[0].on_spec().backscal(0)

        sigma_onoff = significance_on_off(
            n_on=on_counts,
            n_off=off_counts,
            alpha=alpha,
            method='lima'
        )

        print(f"sigma ON/OFF: {sigma_onoff:.2f}")

        # Low Energy PL fitting
        # TODO: TIME INTERVAL TO BE THE ONE IN THE TIME LOOP in the first CTSELECT
        obs_pl = obs.copy()
        select_pl = ctools.ctselect(obs_pl)
        select_pl['rad'] = 3
        select_pl['tmin'] = 0
        select_pl['tmax'] = 1000
        select_pl['emin'] = 0.030
        select_pl['emax'] = 0.080
        select_pl.run()


        # CREATE POWER LAW MODEL
        models_pl_temp = gammalib.GModels()

        # create test source
        src_dir = gammalib.GSkyDir()
        src_dir.radec_deg(0, 0.5)
        spatial = gammalib.GModelSpatialPointSource(src_dir)

        # create and append source spectral model
        spectral = gammalib.GModelSpectralPlaw()
        spectral['Prefactor'].value(5.5e-16)
        spectral['Prefactor'].scale(1e-16)
        spectral['Index'].value(-2.6)
        spectral['Index'].scale(-1.0)
        spectral['PivotEnergy'].value(50000)
        spectral['PivotEnergy'].scale(1e3)
        model_src = gammalib.GModelSky(spatial, spectral)
        model_src.name('PL_fit_temp')
        model_src.tscalc(True)
        models_pl_temp.append(model_src)

        # create and append background
        back_model = gammalib.GCTAModelIrfBackground()
        back_model.instruments('CTA')
        back_model.name('Background')
        spectral_back = gammalib.GModelSpectralPlaw()
        spectral_back['Prefactor'].value(1.0)
        spectral_back['Prefactor'].scale(1.0)
        spectral_back['Index'].value(0)
        spectral_back['PivotEnergy'].value(300000)
        spectral_back['PivotEnergy'].scale(1e6)
        back_model.spectral(spectral_back)
        models_pl_temp.append(back_model)

        # save models
        xmlmodel_PL_temp = 'test_model_PL.xml'
        models_pl_temp.save(xmlmodel_PL_temp)

        # -------------------------------------------------------------
        # --------------------PL FIT---------------------------------
        # -------------------------------------------------------------

        like = ctools.ctlike(select_pl.obs())
        like['inmodel'] = xmlmodel_PL_temp
        like['refit'] = True
        like.run()

        pref_out_pl = like.obs().models()[0]['Prefactor'].value()
        index_out_pl = like.obs().models()[0]['Index'].value()
        pivot_out_pl = like.obs().models()[0]['PivotEnergy'].value()

        # CSPHAGEN creates a "Dummy" source: here I create a model of the GRB to be used in ctlike
        onoff_sim.obs().models()[0].name(grb_name)
        expplaw = gammalib.GModelSpectralExpPlaw()
        expplaw['Prefactor'].value(pref_out_pl)
        expplaw['Index'].value(index_out_pl)
        expplaw['PivotEnergy'].value(pivot_out_pl)
        expplaw['CutoffEnergy'].value(80e3)
        onoff_sim.obs().models()[0].tscalc(True)
        onoff_sim.obs().models()[0].spectral(expplaw)

        like = ctools.ctlike(onoff_sim.obs())
        #like['refit'] = True
        like.run()

        ts_like = like.obs().models()[0].ts()
        E_cut_off = like.obs().models()[0]['CutoffEnergy'].value()
        E_cut_off_error = like.obs().models()[0]['CutoffEnergy'].error()

        print(f"sqrt(TS): {np.sqrt(ts_like):.2f}")
        print(f"E_cut_off: {E_cut_off:.2f} +- {E_cut_off_error:.2f}")

        # -------------------------------------------------------------
        # --------------------MODEL CTLIKE-----------------------------
        # -------------------------------------------------------------

        models_fit_ctlike = gammalib.GModels()

        # create test source
        src_dir = gammalib.GSkyDir()
        src_dir.radec_deg(0, 0.5)
        spatial = gammalib.GModelSpatialPointSource(src_dir)

        # create and append source spectral model
        spectral = gammalib.GModelSpectralExpPlaw()
        spectral['Prefactor'].value(pref_out_pl)
        spectral['Index'].value(index_out_pl)
        spectral['PivotEnergy'].value(pivot_out_pl)
        spectral['CutoffEnergy'].value(80e3)

        model_src = gammalib.GModelSky(spatial, spectral)
        model_src.name('Source_fit')
        model_src.tscalc(True)
        models_fit_ctlike.append(model_src)

        # create and append background
        back_model = gammalib.GCTAModelIrfBackground()
        back_model.instruments('CTA')
        back_model.name('Background')
        spectral_back = gammalib.GModelSpectralPlaw()
        spectral_back['Prefactor'].value(1.0)
        spectral_back['Prefactor'].scale(1.0)
        spectral_back['Index'].value(0)
        spectral_back['PivotEnergy'].value(300000)
        spectral_back['PivotEnergy'].scale(1e6)
        back_model.spectral(spectral_back)
        models_fit_ctlike.append(back_model)

        # save models
        input_ctlike_xml = "model_GRB_fit_ctlike_in.xml"
        models_fit_ctlike.save(input_ctlike_xml)

        # -------------------------------------------------------------
        # --------------------CTLIKE STD-------------------------------
        # -------------------------------------------------------------

        like_std = ctools.ctlike(sel_obs_copy)
        like_std['inmodel'] = input_ctlike_xml
        # like_std['refit'] = True
        like_std.run()

        ts_like_std = like_std.obs().models()[0].ts()

        E_cut_off_ctlike = like_std.obs().models()[0]['CutoffEnergy'].value()
        E_cut_off_ctlike_error = like_std.obs().models()[0]['CutoffEnergy'].error()

        print(f"sqrt(TS) STD: {np.sqrt(ts_like_std):.2f}")
        print(f"E_cut_off: {E_cut_off_ctlike:.2f} +- {E_cut_off_ctlike_error:.2f}")


def gw_simulation(sim_in, config_in, model_xml, counter, background_fits):
    """

    :param sim_in:
    :param config_in:
    :param model_xml:
    :param background_fits:
    :param counter:
    :return:
    """

    print("Working on it...")
    sys.exit()


if __name__ == '__main__':

    sim_yaml_file = yaml.safe_load(open(sys.argv[1]))
    jobs_yaml_file = yaml.safe_load(open(sys.argv[2]))

    if sim_yaml_file['source']['type'] == "GRB":
        grb_simulation(sim_yaml_file, jobs_yaml_file, sys.argv[3], sys.argv[4], sys.argv[5])
    elif sim_yaml_file['source']['type'] == "GW":
        gw_simulation(sim_yaml_file, jobs_yaml_file, sys.argv[3], sys.argv[4], sys.argv[5])
