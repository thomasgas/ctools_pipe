import gammalib
import ctools
import cscripts
import sys
import os
import yaml

from astropy.io.misc import yaml as yl_as
# from astropy

#
# sim = ctools.ctobssim()
# sim['inmodel'] = list_obssim[0]
# sim['caldb'] = list_obssim[1]
# sim['irf'] = list_obssim[2]
# sim['ra'] = list_obssim[3]
# sim['dec'] = list_obssim[4]
# sim['rad'] = list_obssim[5]
# sim['tmin'] = list_obssim[6]
# sim['tmax'] = list_obssim[7]
# sim['emin'] = list_obssim[8]
# sim['emax'] = list_obssim[9]
# sim['outevents'] = nome_fits
# sim['debug'] = debug_sim
# sim['seed'] = seed_in
# sim.execute()

test = yaml.safe_load(open("conf.yaml"))
print(yaml.dump(test))
