import os
import pandas as pd
from astropy.table import Table


try:
    caldb_path = os.environ['CALDB']
    print(caldb_path)
except KeyError:
    print("activate ctools environment!!!")

path_cta = caldb_path+"/data/cta/"
irf_folder_list = os.listdir(path_cta)

try:
    prod3b = os.listdir(path_cta + "prod3b/bcf/")
except:
    print("prod3b not present...not a big deal as long as you have one set of IRF")

try:
    prod3b_v1 = os.listdir(path_cta + "prod3b-v1/bcf/")
except:
    print("prod3b_v1 are not present...not a big deal as long as you have one set of IRF")

try:
    prod3b_v2 = os.listdir(path_cta + "prod3b-v2/bcf/")
except:
    print("prod3b_v2 are not present...not a big deal as long as you have one set of IRF")

# prod3b
try:
    prod3b = os.listdir(path_cta + "prod3b/bcf/")
    d = pd.DataFrame(index=prod3b, columns=["site", "zenith", "pointing", "TS", "subarray", "time", "name"])
    d['name'] = prod3b

    for irf in prod3b:
        d.loc[irf]['site'] = irf.split("_")[0]
        d.loc[irf]['zenith'] = irf.split("_")[1][1:]
        d.loc[irf]['pointing'] = irf.split("_")[2]
        d.loc[irf]['time'] = irf.split("_")[-1]
        if irf.split("_")[3] == "TS":
            d.loc[irf]['TS'] = "1"
        else:
            d.loc[irf]['TS'] = "0"

        if irf.split("_")[-2] == "MSTSST":
            d.loc[irf]['subarray'] = "MSTSST"
        elif irf.split("_")[-2] == "MST":
            d.loc[irf]['subarray'] = "MST"
        elif irf.split("_")[-2] == "SST":
            d.loc[irf]['subarray'] = "SST"
        elif irf.split("_")[-2] == "LST":
            d.loc[irf]['subarray'] = "LST"
        elif irf.split("_")[-2] == "TS":
            d.loc[irf]['subarray'] = "TS"
        else:
            d.loc[irf]['subarray'] = "FULL"

    t = Table.from_pandas(d)
    t.write('prod3b.fits')
except:
    print("prod3b not present...not a big deal as long as you have one set of IRF")

#########################
# prod3b_v1
#########################

try:
    prod3b_v1 = os.listdir(path_cta + "prod3b-v1/bcf/")
    d = pd.DataFrame(index=prod3b_v1, columns=["site", "zenith", "pointing", "time", "name"])
    d['name'] = prod3b_v1

    for irf in prod3b_v1:
        split_list = irf.split("_")
        d.loc[irf]['site'] = split_list[0]
        d.loc[irf]['zenith'] = split_list[1][1:]
        d.loc[irf]['pointing'] = split_list[2]
        d.loc[irf]['time'] = split_list[-1]

    t = Table.from_pandas(d)
    t.write('irf_fits/prod3b-v1.fits')
except:
    print("prod3b_v1 are not present...not a big deal as long as you have one set of IRF")

#########################
# prod3b_v2
#########################

try:
    d = pd.DataFrame(index=prod3b_v2, columns=["site", "zenith", "pointing", "time", "name"])
    d['name'] = prod3b_v2

    for irf in prod3b_v2:
        split_list = irf.split("_")
        d.loc[irf]['site'] = split_list[0]
        d.loc[irf]['zenith'] = split_list[1][1:]
        d.loc[irf]['pointing'] = split_list[2]
        d.loc[irf]['time'] = split_list[-1]

        if split_list[2] == "N":
            d.loc[irf]['pointing'] = "N"
        elif split_list[2] == "S":
            d.loc[irf]['pointing'] = "S"
        else:
            d.loc[irf]['pointing'] = "average"

    t = Table.from_pandas(d)
    t.write('irf_fits/prod3b-v2.fits', overwrite=True)
except:
    print("prod3b_v2 are not present...not a big deal as long as you have one set of IRF")
