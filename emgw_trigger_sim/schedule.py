from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u


def load_schedule_file(schedule_file, mjdobs=0):
    schedule = ascii.read(schedule_file,
                          names=['fieldid', 'ra', 'dec', 'mjd', 'limiting_mag', 'col6', 'col7', 'col8', 'band',
                                 'col10'])
    schedule['phase'] = schedule['mjd'] - mjdobs

    return schedule


def get_fields_coordinates(fields_file):
    survey_fields = ascii.read(fields_file)
    survey_fields_crds = SkyCoord(ra=survey_fields['col2'], dec=survey_fields['col3'], unit=(u.deg, u.deg))

    return survey_fields_crds


def get_survey_field(ra, dec, field_crds):
    crds = SkyCoord(ra=ra,dec=dec,unit=(u.deg,u.deg))
    idx,d2d,d3d = crds.match_to_catalog_sky(field_crds)
    return idx+1


def get_obs_phase(kn_field,schedule):
    obs_sched = schedule[(schedule['fieldid']==kn_field)]
    if len(obs_sched)==0:
        return -99, -99
    else:
        return obs_sched['phase'], obs_sched['limiting_mag']
