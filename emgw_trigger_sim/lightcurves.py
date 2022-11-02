from schedule import get_obs_phase
import numpy as np
import pickle


def load_lightcurve(lcfile):
    with open(lcfile, 'rb') as f:
        lc = pickle.load(f)
    peakmags = []
    for row in lc:
        peakmags.append(np.min(row['mag'][0]))
    lc['peak_mag'] = peakmags
    return lc


def calculate_detection_fractions(survey_fields, distances_dm, schedule, lc):
    kn_detection_fractions = []
    for ind in range(len(distances_dm)):
        kn_dm = distances_dm[ind]
        kn_field = survey_fields[ind]
        kn_t_obs, kn_obs_lim_mag = get_obs_phase(kn_field, schedule)

        if np.all(kn_t_obs < 0):
            kn_detection_fractions.append(0)
        else:
            kn_obs_inds = np.array([np.argmin(np.abs(lc['t'][0] - x)) for x in kn_t_obs])
            kn_obs_mags = np.array([lc['mag'][:, 0, x] + kn_dm for x in kn_obs_inds]).T
            kn_det_mask = kn_obs_mags<kn_obs_lim_mag
            kn_det_frac = np.sum(np.any(kn_det_mask, axis=1))/kn_obs_mags.shape[0]
            # kn_det_frac = np.sum([np.any(kn_obs_mags[x] < kn_obs_lim_mag) for x in range(len(kn_obs_mags))]) / kn_obs_mags.shape[0]
            kn_detection_fractions.append(kn_det_frac)
            # print(ind, len(kn_obs_inds))

    kn_detection_fractions = np.array(kn_detection_fractions)

    return kn_detection_fractions