import sys
import os
import gammalib
import ctools
import cscripts
import yaml
from utils import create_path
from irf_handler import IRFPicker
import astropy.units as u
from astropy.io import fits
from gammapy.stats import significance_on_off
import numpy as np
import glob
from utils import Observability
from astropy.time import Time
import pandas as pd


def grb_simulation(sim_in, config_in, model_xml, fits_header_0, counter):
    """
    Function to handle the GRB simulation.
    :param sim_in: the yaml file for the simulation (unpacked as a dict of dicts)
    :param config_in: the yaml file for the job handling (unpacked as a dict of dicts)
    :param model_xml: the XML model name for the source under analysis
    :param fits_header_0: header for the fits file of the GRB model to use. Used in the visibility calculation
    :param counter: integer number. counts the id of the source realization
    :return: significance obtained with the activated detection methods
    """

    src_name = model_xml.split('/')[-1].split('_')[1][:-4]
    print(src_name, counter)

    ctools_pipe_path = create_path(config_in['exe']['software_path'])
    ctobss_params = sim_in['ctobssim']

    seed = int(counter)*10

    # PARAMETERS FROM THE CTOBSSIM
    sim_t_min = u.Quantity(ctobss_params['time']['t_min']).to_value(u.s)
    sim_t_max = u.Quantity(ctobss_params['time']['t_max']).to_value(u.s)
    sim_e_min = u.Quantity(ctobss_params['energy']['e_min']).to_value(u.TeV)
    sim_e_max = u.Quantity(ctobss_params['energy']['e_max']).to_value(u.TeV)
    sim_rad = ctobss_params['radius']

    output_path = create_path(sim_in['output'] + '/' + src_name)

    save_simulation = ctobss_params['save_simulation']

    with open(f"{output_path}/GRB-{src_name}_seed-{seed}.txt", "w") as f:
        f.write(f"GRB,seed,time_start,time_end,sigma_lima,sqrt_TS_onoff,sqrt_TS_std\n")
        # VISIBILITY PART
        # choose between AUTO mode (use visibility) and MANUAL mode (manually insert IRF)
        simulation_mode = sim_in['IRF']['mode']

        if simulation_mode == "auto":
            print("using visibility to get IRFs")

            # GRB information from the fits header
            ra = fits_header_0['RA']
            dec = fits_header_0['DEC']
            t0 = Time(fits_header_0['GRBJD'])

            irf_dict = sim_in['IRF']
            site = irf_dict['site']
            obs_condition = Observability(site=site)
            obs_condition.set_irf(irf_dict)

            t_zero_mode = ctobss_params['time']['t_zero'].lower()

            if t_zero_mode == "VIS":
                # check if the source is visible one day after the onset of the source
                print("time starts when source becomes visible")
                obs_condition.Proposal_obTime = 86400
                condition_check = obs_condition.check(RA=ra, DEC=dec, t_start=t0)

            elif t_zero_mode == "ONSET":
                print("time starts from the onset of the GRB")
                condition_check = obs_condition.check(RA=ra, DEC=dec, t_start=t0, t_min=sim_t_min, t_max=sim_t_max)

            else:
                print(f"Choose some proper mode between 'VIS' and 'ONSET'. {t_zero_mode} is not a valid one.")
                sys.exit()

            # NO IRF in AUTO mode ==> No simulation! == EXIT!
            if len(condition_check) == 0:
                f.write(f"{src_name},{seed}, -1, -1, -1, -1, -1\n")
                sys.exit()

        elif simulation_mode == "manual":
            print("manual picking IRF")

            # find proper IRF name
            irf = IRFPicker(sim_in, ctools_pipe_path)
            name_irf = irf.irf_pick()

            backgrounds_path = create_path(ctobss_params['bckgrnd_path'])
            fits_background_list = glob.glob(
                f"{backgrounds_path}/{irf.prod_number}_{irf.prod_version}_{name_irf}/background*.fits")

            if len(fits_background_list) == 0:
                print(f"No background for IRF {name_irf}")
                sys.exit()

            background_fits = fits_background_list[int(counter) - 1]
            obs_back = gammalib.GCTAObservation(background_fits)

        else:
            print(f"wrong input for IRF - mode. Input is {simulation_mode}. Use 'auto' or 'manual' instead")
            sys.exit()

        if irf.prod_number == "3b" and irf.prod_version == 0:
            caldb = "prod3b"
        else:
            caldb = f'prod{irf.prod_number}-v{irf.prod_version}'

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

        obs = sim.obs()

        # # move the source photons from closer to (RA,DEC)=(0,0), where the background is located
        # for event in obs[0].events():
        #     # ra_evt = event.dir().dir().ra()
        #     dec_evt = event.dir().dir().dec()
        #     ra_evt_deg = event.dir().dir().ra_deg()
        #     dec_evt_deg = event.dir().dir().dec_deg()
        #
        #     ra_corrected = (ra_evt_deg - ra_pointing)*np.cos(dec_evt)
        #     dec_corrected = dec_evt_deg - dec_pointing
        #     event.dir().dir().radec_deg(ra_corrected, dec_corrected)

        # append all background events to GRB ones ==> there's just one observation and not two
        for event in obs_back.events():
            obs[0].events().append(event)

        # ctselect to save data on disk
        if save_simulation:
            select_time = ctools.ctselect(obs)
            select_time['rad'] = sim_rad
            select_time['tmin'] = sim_t_min
            select_time['tmax'] = sim_t_max
            select_time['emin'] = sim_e_min
            select_time['emax'] = sim_e_max
            event_list_path = create_path(f"{ctobss_params['output_path']}/{src_name}/seed-{seed:03}/")
            sim['outobs'] = f"{event_list_path}/event_list_source-{src_name}_seed-{seed:03}.fits"
            select_time.execute()
            sys.exit()

        # delete all 70+ models from the obs def file...not needed any more
        obs.models(gammalib.GModels())

        # CTSELECT
        select_time = sim_in['ctselect']['time_cut']
        slices = int(select_time['t_slices'])

        if slices == 0:
            times = [sim_t_min, sim_t_max]
            times_start = times[:-1]
            times_end = times[1:]
        elif slices > 0:
            time_mode = select_time['mode']
            if time_mode == "log":
                times = np.logspace(np.log10(sim_t_min), np.log10(sim_t_max), slices + 1, endpoint=True)
            elif time_mode == "lin":
                times = np.linspace(sim_t_min, sim_t_max, slices + 1, endpoint=True)
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

        # ------------------------------------
        # ----- TIME LOOP STARTS HERE --------
        # ------------------------------------

        ctlike_mode = sim_in['detection']
        mode_1 = ctlike_mode['counts']
        mode_2 = ctlike_mode['ctlike-onoff']
        mode_3 = ctlike_mode['ctlike-std']

        for t_in, t_end in zip(times_start, times_end):
            sigma_onoff = 0
            sqrt_ts_like_onoff = 0
            sqrt_ts_like_std = 0
            print("-----------------------------")
            print(f"t_in: {t_in:.2f}, t_end: {t_end:.2f}")

            # different ctlikes (onoff or std) need different files.
            # will be appended here and used later on for the final likelihood
            dict_obs_select_time = {}

            # perform time selection for this specific time bin
            select_time = ctools.ctselect(obs)
            select_time['rad'] = sim_rad
            select_time['tmin'] = t_in
            select_time['tmax'] = t_end
            select_time['emin'] = sim_e_min
            select_time['emax'] = sim_e_max
            select_time.run()

            if mode_1:
                fits_temp_title = f"skymap_{seed}_{t_in:.2f}_{t_end:.2f}.fits"
                pars_counts = ctlike_mode['pars_counts']
                scale = float(pars_counts['scale'])
                npix = 2*int(sim_rad/scale)
                skymap = ctools.ctskymap(select_time.obs().copy())
                skymap['emin'] = sim_e_min
                skymap['emax'] = sim_e_max
                skymap['nxpix'] = npix
                skymap['nypix'] = npix
                skymap['binsz'] = scale
                skymap['proj'] = 'TAN'
                skymap['coordsys'] = 'CEL'
                skymap['xref'] = 0
                skymap['yref'] = 0
                skymap['bkgsubtract'] = 'RING'
                skymap['roiradius'] = pars_counts['roiradius']
                skymap['inradius'] = pars_counts['inradius']
                skymap['outradius'] = pars_counts['outradius']
                skymap['iterations'] = pars_counts['iterations']
                skymap['threshold'] = pars_counts['threshold']
                skymap['outmap'] = fits_temp_title
                skymap.execute()

                input_fits = fits.open(fits_temp_title)
                datain = input_fits[2].data
                datain[np.isnan(datain)] = 0.0
                datain[np.isinf(datain)] = 0.0

                sigma_onoff = np.max(datain)
                os.remove(fits_temp_title)

            if mode_3:
                dict_obs_select_time['std'] = select_time.obs().copy()

            if mode_2:
                onoff_time_sel = cscripts.csphagen(select_time.obs().copy())
                onoff_time_sel['inmodel'] = 'NONE'
                onoff_time_sel['ebinalg'] = 'LOG'
                onoff_time_sel['emin'] = sim_e_min
                onoff_time_sel['emax'] = sim_e_max
                onoff_time_sel['enumbins'] = 30
                onoff_time_sel['coordsys'] = 'CEL'
                onoff_time_sel['ra'] = 0.0
                onoff_time_sel['dec'] = 0.5
                onoff_time_sel['rad'] = 0.2
                onoff_time_sel['bkgmethod'] = 'REFLECTED'
                onoff_time_sel['use_model_bkg'] = False
                onoff_time_sel['stack'] = False
                onoff_time_sel.run()

                dict_obs_select_time['onoff'] = onoff_time_sel.obs().copy()

                del onoff_time_sel

                # print(f"sigma ON/OFF: {sigma_onoff:.2f}")

            if mode_2 or mode_3:

                # Low Energy PL fitting
                # to be saved in this dict
                dict_pl_ctlike_out = {}

                e_min_pl_ctlike = 0.030
                e_max_pl_ctlike = 0.080

                # simple ctobssim copy and select for ctlike-std
                select_pl_ctlike = ctools.ctselect(select_time.obs().copy())
                select_pl_ctlike['rad'] = 3
                select_pl_ctlike['tmin'] = t_in
                select_pl_ctlike['tmax'] = t_end
                select_pl_ctlike['emin'] = e_min_pl_ctlike
                select_pl_ctlike['emax'] = e_max_pl_ctlike
                select_pl_ctlike.run()

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

                spectral_back = gammalib.GModelSpectralPlaw()
                spectral_back['Prefactor'].value(1.0)
                spectral_back['Prefactor'].scale(1.0)
                spectral_back['Index'].value(0)
                spectral_back['PivotEnergy'].value(300000)
                spectral_back['PivotEnergy'].scale(1e6)

                if mode_2:
                    back_model = gammalib.GCTAModelIrfBackground()
                    back_model.instruments('CTAOnOff')
                    back_model.name('Background')
                    back_model.spectral(spectral_back.copy())

                    onoff_pl_ctlike_lima = cscripts.csphagen(select_pl_ctlike.obs().copy())
                    onoff_pl_ctlike_lima['inmodel'] = 'NONE'
                    onoff_pl_ctlike_lima['ebinalg'] = 'LOG'
                    onoff_pl_ctlike_lima['emin'] = e_min_pl_ctlike
                    onoff_pl_ctlike_lima['emax'] = e_max_pl_ctlike
                    onoff_pl_ctlike_lima['enumbins'] = 30
                    onoff_pl_ctlike_lima['coordsys'] = 'CEL'
                    onoff_pl_ctlike_lima['ra'] = 0.0
                    onoff_pl_ctlike_lima['dec'] = 0.5
                    onoff_pl_ctlike_lima['rad'] = 0.2
                    onoff_pl_ctlike_lima['bkgmethod'] = 'REFLECTED'
                    onoff_pl_ctlike_lima['use_model_bkg'] = False
                    onoff_pl_ctlike_lima['stack'] = False
                    onoff_pl_ctlike_lima.run()

                    onoff_pl_ctlike_lima.obs().models(gammalib.GModels())
                    onoff_pl_ctlike_lima.obs().models().append(model_src.copy())
                    onoff_pl_ctlike_lima.obs().models().append(back_model.copy())

                    like_pl = ctools.ctlike(onoff_pl_ctlike_lima.obs())
                    like_pl['refit'] = True
                    like_pl.run()
                    dict_pl_ctlike_out['onoff'] = like_pl.obs().copy()
                    del onoff_pl_ctlike_lima
                    del like_pl

                if mode_3:
                    models_ctlike_std = gammalib.GModels()
                    models_ctlike_std.append(model_src.copy())
                    back_model = gammalib.GCTAModelIrfBackground()
                    back_model.instruments('CTA')
                    back_model.name('Background')
                    back_model.spectral(spectral_back.copy())
                    models_ctlike_std.append(back_model)

                    # save models
                    xmlmodel_PL_ctlike_std = 'test_model_PL_ctlike_std.xml'
                    models_ctlike_std.save(xmlmodel_PL_ctlike_std)
                    del models_ctlike_std

                    like_pl = ctools.ctlike(select_pl_ctlike.obs().copy())
                    like_pl['inmodel'] = xmlmodel_PL_ctlike_std
                    like_pl['refit'] = True
                    like_pl.run()
                    dict_pl_ctlike_out['std'] = like_pl.obs().copy()
                    del like_pl

                del spatial
                del spectral
                del model_src
                del select_pl_ctlike

                # EXTENDED CTLIKE
                for key in dict_obs_select_time.keys():
                    likelihood_pl_out = dict_pl_ctlike_out[key]
                    selected_data = dict_obs_select_time[key]

                    pref_out_pl = likelihood_pl_out.models()[0]['Prefactor'].value()
                    index_out_pl = likelihood_pl_out.models()[0]['Index'].value()
                    pivot_out_pl = likelihood_pl_out.models()[0]['PivotEnergy'].value()

                    expplaw = gammalib.GModelSpectralExpPlaw()
                    expplaw['Prefactor'].value(pref_out_pl)
                    expplaw['Index'].value(index_out_pl)
                    expplaw['PivotEnergy'].value(pivot_out_pl)
                    expplaw['CutoffEnergy'].value(80e3)

                    if key == "onoff":
                        selected_data.models()[0].name(src_name)
                        selected_data.models()[0].tscalc(True)
                        selected_data.models()[0].spectral(expplaw.copy())

                        like = ctools.ctlike(selected_data)
                        like['refit'] = True
                        like.run()
                        ts = like.obs().models()[0].ts()
                        if ts > 0:
                            sqrt_ts_like_onoff = np.sqrt(like.obs().models()[0].ts())
                        else:
                            sqrt_ts_like_onoff = 0

                        del like

                    if key == "std":
                        models_fit_ctlike = gammalib.GModels()

                        # create test source
                        src_dir = gammalib.GSkyDir()
                        src_dir.radec_deg(0, 0.5)
                        spatial = gammalib.GModelSpatialPointSource(src_dir)

                        # append spatial and spectral models
                        model_src = gammalib.GModelSky(spatial, expplaw.copy())
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
                        del models_fit_ctlike

                        like = ctools.ctlike(selected_data)
                        like['inmodel'] = input_ctlike_xml
                        like['refit'] = True
                        like.run()
                        ts = like.obs().models()[0].ts()
                        if ts > 0:
                            sqrt_ts_like_std = np.sqrt(like.obs().models()[0].ts())
                        else:
                            sqrt_ts_like_std = 0

                        del like

                    # E_cut_off = like.obs().models()[0]['CutoffEnergy'].value()
                    # E_cut_off_error = like.obs().models()[0]['CutoffEnergy'].error()

                    # print(f"sqrt(TS) {key}: {np.sqrt(ts_like):.2f}")
                    # print(f"E_cut_off {key}: {E_cut_off:.2f} +- {E_cut_off_error:.2f}")
                del dict_pl_ctlike_out

            f.write(f"{src_name},{seed},{t_in:.2f},{t_end:.2f},{sigma_onoff:.2f},{sqrt_ts_like_onoff:.2f},{sqrt_ts_like_std:.2f}\n")
            del dict_obs_select_time
            del select_time


def gw_simulation(sim_in, config_in, model_xml, fits_model, counter):
    """

    :param sim_in:
    :param config_in:
    :param model_xml:
    :param fits_model:
    :param counter:
    :return:
    """

    src_name = fits_model.split("/")[-1][:-5]

    fits_header_0 = fits.open(fits_model)[0].header
    ra_src = fits_header_0['RA']
    dec_src = fits_header_0['DEC']

    src_yaml = sim_in['source']

    point_path = create_path(src_yaml['pointings_path'])
    opt_point_path = f"{point_path}/optimized_pointings"

    ctools_pipe_path = create_path(config_in['exe']['software_path'])
    ctobss_params = sim_in['ctobssim']

    seed = int(counter)*10

    # # PARAMETERS FROM THE CTOBSSIM
    sim_e_min = u.Quantity(ctobss_params['energy']['e_min']).to_value(u.TeV)
    sim_e_max = u.Quantity(ctobss_params['energy']['e_max']).to_value(u.TeV)

    sim_rad = ctobss_params['radius']
    output_path = create_path(sim_in['output'] + '/' + src_name)

    irf_dict = sim_in['IRF']
    site = irf_dict['site']

    detection = sim_in['detection']
    significance_map = detection['counts']
    srcdetect_ctlike = detection['detect_like']

    save_simulation = ctobss_params['save_simulation']

    with open(f"{output_path}/GW-{src_name}_seed-{seed:03}_site-{site}.txt", "w") as f:
        f.write(f"GW_name,RA_src,DEC_src,seed,pointing_id,ra_point,dec_point,radius,time_start,time_end,sigma\n")
        run_id, merger_id = src_name.split('_')
        pointing_data = pd.read_csv(
            f"{opt_point_path}/{run_id}_Merger{merger_id}_SuggestedPointings_GWOptimisation_test.txt", header=0,
            sep=" ")

        RA_data = pointing_data['RA(deg)']
        DEC_data = pointing_data['DEC(deg)']
        times = pointing_data['Observation Time UTC']

        # LOOP OVER POINTINGS
        for index in range(0, len(pointing_data) - 1):
            RA_point = RA_data[index]
            DEC_point = DEC_data[index]
            t_in_point = Time(times[index])
            t_end_point = Time(times[index + 1])

            obs_condition = Observability(site=site)
            obs_condition.set_irf(irf_dict)
            obs_condition.Proposal_obTime = 10
            obs_condition.TimeOffset = 0
            obs_condition.Steps_observability = 10
            condition_check = obs_condition.check(RA=RA_point, DEC=DEC_point, t_start=t_in_point)

            if len(condition_check) == 0:
                f.write(
                    f"{src_name},{ra_src},{dec_src},{seed},{index},{RA_point},{DEC_point},{sim_rad},{t_in_point},{t_end_point}, -1\n")
                continue

            name_irf = condition_check['IRF_name'][0]
            irf = condition_check['IRF'][0]
            # model loading

            if irf.prod_number == "3b" and irf.prod_version == 0:
                caldb = "prod3b"
            else:
                caldb = f'prod{irf.prod_number}-v{irf.prod_version}'

            # simulation
            sim = ctools.ctobssim()
            sim['inmodel'] = model_xml
            sim['caldb'] = caldb
            sim['irf'] = name_irf
            sim['ra'] = RA_point
            sim['dec'] = DEC_point
            sim['rad'] = sim_rad
            sim['tmin'] = t_in_point.value
            sim['tmax'] = t_end_point.value
            sim['emin'] = sim_e_min
            sim['emax'] = sim_e_max
            sim['seed'] = seed

            if save_simulation:
                event_list_path = create_path(f"{ctobss_params['output_path']}/{src_name}/seed-{seed:03}/")
                sim['outevents'] = f"{event_list_path}/event_list_source-{src_name}_seed-{seed:03}_pointingID-{index}.fits"
                sim.execute()
                continue
            else:
                sim.run()

            obs = sim.obs()
            # ctskymap

            skymap = ctools.ctskymap(obs)
            skymap['proj'] = 'CAR'
            skymap['coordsys'] = 'CEL'
            skymap['xref'] = RA_point
            skymap['yref'] = DEC_point
            skymap['binsz'] = 0.02
            skymap['nxpix'] = 200
            skymap['nypix'] = 200
            skymap['emin'] = sim_e_min
            skymap['emax'] = sim_e_max
            skymap['bkgsubtract'] = 'NONE'
            skymap.run()

            # cssrcdetect
            srcdetect = cscripts.cssrcdetect(skymap.skymap().copy())
            srcdetect['srcmodel'] = 'POINT'
            srcdetect['bkgmodel'] = 'NONE'  # we will determine the background model in a later step
            srcdetect['threshold'] = 5
            srcdetect['corr_rad'] = 0.1
            srcdetect.run()

            models = srcdetect.models()
            if len(models) > 0:
                hotspot = models['Src001']
                ra_hotspot = hotspot['RA'].value()
                dec_hotspot = hotspot['DEC'].value()
            else:
                ra_hotspot = 0
                dec_hotspot = 0

            f.write(
                f"{src_name},{ra_src},{dec_src},{seed},{index},{RA_point},{DEC_point},{sim_rad},{t_in_point},{t_end_point},0\n")

    #     # VISIBILITY PART
    #     # choose between AUTO mode (use visibility) and MANUAL mode (manually insert IRF)
    #     # simulation_mode = sim_in['IRF']['mode']
    #     #
    #     # if simulation_mode == "auto":
    #     #     print("using visibility to get IRFs")
    #     #
    #     #     # GRB information from the fits header
    #     #     ra = fits_header_0['RA']
    #     #     dec = fits_header_0['DEC']
    #     #     t0 = Time(fits_header_0['GRBJD'])
    #     #
    #     #     irf_dict = sim_in['IRF']
    #     #     site = irf_dict['site']
    #     #     obs_condition = Observability(site=site)
    #     #     obs_condition.set_irf(irf_dict)
    #     #
    #     #     t_zero_mode = ctobss_params['time']['t_zero'].lower()
    #     #     print("time starts when source becomes visible")
    #     #     obs_condition.Proposal_obTime = 86400
    #     #     condition_check = obs_condition.check(RA=ra, DEC=dec, t_start=t0)
    #     #
    #     # else:
    #     #     print(f"wrong input for IRF - mode. Input is {simulation_mode}. Only 'AUTO' is supported in GW analysis")
    #     #     sys.exit()
    #


if __name__ == '__main__':

    sim_yaml_file = yaml.safe_load(open(sys.argv[1]))
    jobs_yaml_file = yaml.safe_load(open(sys.argv[2]))
    xml_model_input = sys.argv[3]
    fits_model_input = sys.argv[4]
    realization_id = sys.argv[5]

    if sim_yaml_file['source']['type'] == "GRB":
        header_GRB = fits.open(fits_model_input)[0].header
        grb_simulation(sim_yaml_file, jobs_yaml_file, xml_model_input, header_GRB, realization_id)
    elif sim_yaml_file['source']['type'] == "GW":
        gw_simulation(sim_yaml_file, jobs_yaml_file, xml_model_input, fits_model_input, realization_id)
