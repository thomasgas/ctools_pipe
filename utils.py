from environs import Env
import os

import numpy as np
from astropy import units as u
from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, get_sun, get_moon)

from astropy.time import Time

from irf_handler import IRFPicker



# # Some defifinitions used in func observability -----------
# # Array N coordinates
# CTAN_latitude = 28.762164
# CTAN_longitude = -17.892005
#
# # Array S coordinates
# CTAS_latitude = -24.683428
# CTAS_longitude = -70.316344
#
# Proposal_obTime = 144000.0  # (s) checking for 4h of observation!
# Steps_observability = 600.0  # (s) the observability is checked with steps of 300s
# Alert_zenith_max = 66.0  # (deg)
# Sun_zenith_min = 105.0  # (deg)  (105-108 ?)
# Moon_separation = 30.0  # (deg)
#
# # TimeOffset = 0.0  # (in days) 50-years time shift
# TimeOffset = 18262.0  # (in days) 50-years time shift
#
# irfprod = '3b'  # IRF production
# irfversion = 2  # IRF version
# sub = 'FULL'  # options for the array are LST/MST/SST/MSTSST/FULL/TS
# thr = 0  # threshold array (y/n)
# timeref = '0.5h'  # Observation Time


# ----------------------------------------------------------


def create_path(path_to_expand):
    env = Env()
    out_path = path_to_expand

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

    if out_path.endswith("/"):  # just used to remove backslash if applied
        out_path = out_path[:-1]

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    return out_path


class Observability:
    def __init__(self, site, longitude=0, latitude=0, height=0):
        """
        Build the site either from the name "North" and "South" of the CTA arrays or another array.
        Some parameters are set as default. Change them before using the CHECK function (if needed).
        :param site: site name
        :param longitude: longitude in degrees
        :param latitude: latitude in degrees
        :param height: height above sea level in meters
        """
        self.site = site
        if site == "North":
            self.longitude = -17.892005
            self.latitude = 28.762164
            self.height = 2200 * u.m
        elif site == "South":
            self.longitude = -70.316344
            self.latitude = -24.683428
            self.height = 2150 * u.m
        else:
            self.latitude = latitude
            self.longitude = longitude
            self.height = height * u.m

        self.earth_location = EarthLocation(lon=self.longitude, lat=self.latitude, height=self.height)

        self.Proposal_obTime = 144000.0  # (s) checking for 4h of observation!
        self.Steps_observability = 600.0  # (s) the observability is checked with steps of 300s
        self.Alert_zenith_max = 66.0  # (deg)
        self.Sun_zenith_min = 105.0  # (deg)  (105-108 ?)
        self.Moon_separation = 30.0  # (deg)

        self.TimeOffset = 18262.0  # (in days) 50-years time shift
        self.prod = None
        self.prod_version = None
        self.window = None
        self.subarray = None
        self.threshold = None

    def set_irf(self, input_dict):
        self.prod = input_dict['prod']['number']
        self.prod_version = input_dict['prod']['version']
        self.window = input_dict['time']
        self.subarray = input_dict['subarray']
        self.threshold = input_dict['TS']

    def check(self, RA, DEC, t_start):

        alert_coord = SkyCoord(RA, DEC, frame='fk5', unit='deg')

        alert_time = Time(t_start - self.TimeOffset, format='jd', scale='utc')

        # We check the observability for 4h
        allocated_obtime = np.arange(0, self.Proposal_obTime + self.Steps_observability, self.Steps_observability) * u.s

        observing_time = alert_time + allocated_obtime

        obs_frame = AltAz(location=self.earth_location, obstime=observing_time)

        # We calculate  (az, alt) in deg [in this order!]
        local = alert_coord.transform_to(obs_frame)

        alert_zenith = 90. - local.alt.deg
        alert_azimuth = local.az.deg

        Sun_coos = get_sun(observing_time).transform_to(obs_frame)
        Sun_Zenith = 90.0 - Sun_coos.alt.deg

        moon_coos = get_moon(observing_time).transform_to(obs_frame)

        # calculating separation btw GRB and moon
        safe_moon = local.separation(moon_coos).deg

        mask = np.logical_and(np.logical_and(np.logical_and(alert_zenith >= 0, alert_zenith <= self.Alert_zenith_max),
                                             Sun_Zenith > self.Sun_zenith_min), safe_moon > self.Moon_separation)

        observable_mask = np.argwhere(mask == True)

        observable = observing_time[mask]

        IRF_dict = {}

        for obs_id in observable_mask:
            time_in_window = observing_time[obs_id]

            if 0.0 < alert_zenith[obs_id] <= 30.0:
                zenithIRF = 20
            elif 30.0 < alert_zenith[obs_id] <= 50.0:
                zenithIRF = 40
            elif 50.0 < alert_zenith[obs_id] <= 66.0:
                zenithIRF = 60
            else:
                zenithIRF = None

            if (0.0 < alert_azimuth[obs_id] <= 45.0) or (
                    315.0 < alert_azimuth[obs_id] <= 360.0):
                pointingIRF = 'N'
            elif 135.0 < alert_azimuth[obs_id] <= 225.0:
                pointingIRF = 'S'
            elif (45.0 < alert_azimuth[obs_id] <= 135.0) or (
                    225.0 < alert_azimuth[obs_id] <= 315.0):
                pointingIRF = 'average'
            else:
                pointingIRF = None

            IRFstruct = {'IRF': {'prod': {'number': self.prod, 'version': self.prod_version},
                                 'zenith': zenithIRF,
                                 'site': 'North',
                                 'time': self.window,
                                 'pointing': pointingIRF,
                                 'subarray': self.subarray,
                                 'TS': self.threshold}}

            irf = IRFPicker(IRFstruct)
            IRF_name = irf.irf_pick()
            IRF_dict[f"{time_in_window[0].value:.4f}"] = irf

        return IRF_dict
