import argparse
import os
import healpy as hp
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from utils import read_flattened_skymap, plot_skymap, flatten_skymap
from ingest import get_event_coordinates
from schedule import load_schedule_file, get_fields_coordinates, get_survey_field
from lightcurves import load_lightcurve, calculate_detection_fractions

if __name__ == '__main__':
    skymap_path = '/Users/viraj/winter/pycbc/realistic/realistic_191_VK.fits'
    schedule_file = '/Users/viraj/winter/gwemopt_sims/output_parallel/pycbc/7days/realistic/191.0/schedule_WINTER.dat'
    flatten = True
    flatten_files_dir = '/Users/viraj/emgw_trigger_criteria/data/'
    flatten_file_path = os.path.join(flatten_files_dir,skymap_path.split('/')[-1].split('.fits')[0]+'_flat.fits')
    plot = True
    fields_file = '/Users/viraj/winter/gwemopt_sims/input/WINTER_fields.txt'
    filt = 'J'
    lcfile = f'/Users/viraj/gwem_test/gwemlightcurves/output/bns_Bulla_parameter_grid_{filt}band.dat'
    num_trans = 10000

    if flatten:
        flatten_skymap(skymap_path=skymap_path, flatten_file_path=flatten_file_path, update_object_name=True)
    (prob, distmu, distsigma, distnorm), hdr = read_flattened_skymap(flatten_file_path)

    hdr_keys = np.array([x[0] for x in hdr])
    hdr_vals = np.array([x[1] for x in hdr])
    mjdobs = float(hdr_vals[hdr_keys == 'MJD-OBS'][0])

    if plot:
        plot_skymap(prob)

    ras, decs, distances_mpc = get_event_coordinates(prob, distnorm, distmu, distsigma, num_trans=num_trans)
    crds = SkyCoord(ra=ras, dec=decs, unit=(u.deg, u.deg))
    distances_dm = 5 * np.log10(distances_mpc * 1e6 / 10)

    schedule = load_schedule_file(schedule_file=schedule_file, mjdobs=mjdobs)
    survey_fields_coords = get_fields_coordinates(fields_file=fields_file)
    survey_fields = get_survey_field(ras, decs, survey_fields_coords)

    print('Loading lightcurves')
    lc = load_lightcurve(lcfile=lcfile)
    print('Calculating detection fractions')
    print('Survey fields', survey_fields)
    kn_detection_fractions = calculate_detection_fractions(survey_fields, distances_dm, schedule, lc)
    print(np.mean(kn_detection_fractions))