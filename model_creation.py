from astropy.io import fits
import numpy as np
import astropy.units as u
from astropy.units import UnitsWarning
import os
import sys
import yaml
import glob
import subprocess
import warnings
from utils import create_path


def create_models(input_yaml, jobs_yaml):
    """
    Create fits and txt files for the generation of time variable spectra
    :param input_yaml: YAML configuration file
    :param jobs_yaml: YAML configuration file for job details (software path & co)
    :return: None (data saved to disk)
    """
    config_in = yaml.safe_load(open(input_yaml))
    jobs_in = yaml.safe_load(open(jobs_yaml))
    ctools_pipe_path = create_path(jobs_in['exe']['software_path'])

    models = config_in['models']
    input_data = create_path(models['input_data']['path'])
    output_data = create_path(jobs_in['exe']['path'])

    model_folder = f"{output_data}"
    try:
        os.mkdir(model_folder)
    except FileExistsError:
        print(f'{model_folder} already created. will be use for output source models')

    list_files = glob.glob(f"{input_data}/*fits")

    # create example fits
    time_vals = np.zeros(4)
    norm_vals = np.array([0, 1, 1, 0])

    time_fits = fits.Column(name='TIME', array=time_vals, format='1D')
    norm_fits = fits.Column(name='NORM', array=norm_vals, format='1D')
    output = fits.BinTableHDU.from_columns([time_fits, norm_fits])

    hdr = fits.Header()
    output.header['EXTNAME'] = 'Time Profile'
    output.header.comments['EXTNAME'] = "Name of this extension."
    output.header['MJDREFI'] = 51544
    output.header.comments['MJDREFI'] = "[days] Integer part of time reference MJD."
    output.header['MJDREFF'] = 0.5
    output.header.comments['MJDREFF'] = "[days] Fractional part of time reference MJD."
    output.header['TIMEUNIT'] = 's'
    output.header.comments['TIMEUNIT'] = "Time units."
    output.header['TIMESYS'] = 'TT'
    output.header.comments['TIMESYS'] = "Time system."
    output.header['TIMEREF'] = 'LOCAL'
    output.header.comments['TIMEREF'] = "Time reference."
    output.header['TUNIT1'] = 's'

    primary_hdu = fits.PrimaryHDU(header=hdr)
    hdul = fits.HDUList([primary_hdu, output])

    for file_name in list_files:
        grb_name = file_name.split('/')[-1][:-5]
        grb_test = fits.open(f"{file_name}")
        try:
            os.mkdir(f'{model_folder}/{grb_name}')
        except FileExistsError:
            print(f'{grb_name} already created')

        try:
            os.mkdir(f'{model_folder}/{grb_name}/spectra')
            os.mkdir(f'{model_folder}/{grb_name}/lightcv')
        except FileExistsError:
            print(f'spectra and lightcv already created')

        header_prim = grb_test['PRIMARY'].header

        header_en = grb_test['ENERGIES'].header

        # energies (with proper units)
        energy_unit_fits = u.Unit(header_en['TUNIT1'])
        energy_unit_ctools = u.MeV
        energies = u.Quantity(
            grb_test['ENERGIES'].data['Energies'],
            energy_unit_fits
        ).to_value(energy_unit_ctools)

        # times (with proper units)
        header_time = grb_test['TIMES'].header
        time_unit_fits = u.Unit(header_time['TUNIT1'])
        time_unit_ctools = u.s
        times = u.Quantity(
            grb_test['TIMES'].data['Times'],
            time_unit_fits
        ).to_value(time_unit_ctools)

        # spectra: try to use EBL absorbed spectra (GRB-like)...otherwise use unabsorbed (more GW-like)
        # Units will be added during the loop otherwise it's too slow
        try:
            header_spec = grb_test['EBL-ABS. SPECTRA'].header
            spectra = grb_test['EBL-ABS. SPECTRA'].data
        except KeyError:
            header_spec = grb_test['SPECTRA'].header
            spectra = grb_test['SPECTRA'].data

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UnitsWarning)
            flux_unit_fits = u.Unit(header_spec['UNITS'])

        flux_unit_ctools = u.Unit('ph / (MeV * s *cm^2)')

        model_title = f"{model_folder}/{grb_name}/model_{grb_name}.txt"
        TS = 1
        RA = 0
        DEC = 0
        with open(model_title, 'w') as file:
            for counter, (time_in, time_out) in enumerate(zip(times[:-1], times[1:])):
                fits_name = f"{model_folder}/{grb_name}/lightcv/lc_{str(counter).zfill(3)}_tin-{time_in:.3f}_tend-{time_out:.3f}.fits"
                spec_name = f"{model_folder}/{grb_name}/spectra/spec_{str(counter).zfill(3)}_tin-{time_in:.3f}_tend-{time_out:.3f}.txt"

                save_spectra = u.Quantity(spectra.field(counter), flux_unit_fits).to_value(flux_unit_ctools)
                save_spectra = np.clip(save_spectra, a_min=1e-200, a_max=1.0)
                np.savetxt(spec_name, np.c_[energies, save_spectra], newline="\n")

                deltat = time_out - time_in

                # empirical way to create a square signal according to the delta,
                # which varies between one bin and the other
                delta_scale = 1000.
                fits_time = [time_in - deltat / (2 * delta_scale),
                             time_in,
                             time_out,
                             time_out + deltat / (2 * delta_scale)]

                hdul['Time Profile'].data['TIME'] = np.array(fits_time)
                hdul.writeto(fits_name, overwrite=True)

                source = f"{grb_name}_{counter} Point  {TS}  {RA} {DEC} 0 0 0 FUNC 1.0  {spec_name} 1.0 {fits_name} \n"
                file.write(source)

        # create XML model
        subprocess.Popen([
            "python",
            f"{ctools_pipe_path}/scriptModel_variable.py",
            model_title]
        )


if __name__ == '__main__':
    create_models(sys.argv[1], sys.argv[2])
