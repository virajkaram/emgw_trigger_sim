import healpy as hp
import numpy as np
from ligo.skymap.postprocess import find_greedy_credible_levels
from scipy.stats import norm
import matplotlib.pyplot as plt
import pickle
from astropy.io import fits
from astropy.table import Table
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
import os
from astropy.io import ascii


def read_skymap_fits(filename):
    h = fits.open(filename)
    data = Table(h[1].data)
    prob = np.array(data['PROB'])
    distmu = np.array(data['DISTMU'])
    distsigma = np.array(data['DISTSIGMA'])
    distnorm = np.array(data['DISTNORM'])
    return prob, distmu, distsigma, distnorm


def get_credible_levels(prob):
    npix = len(prob)
    credible_levels = find_greedy_credible_levels(prob)
    nside = hp.npix2nside(npix)
    return credible_levels


def dist_function(r, distnorm, distmu, distsigma):
    dp_dr = r ** 2 * distnorm * norm(distmu, distsigma).pdf(r)
    return dp_dr


def flatten_skymap(skymap_path, flatten_file_path, update_object_name=True):
    cmd = f'ligo-skymap-flatten {skymap_path} {flatten_file_path}'
    print(cmd)
    os.system(cmd)

    if update_object_name:
        h = fits.open(flatten_file_path, 'update')
        h[1].header['OBJECT'] = 'S' + str(h[1].header['OBJECT'])
        h.writeto(flatten_file_path, overwrite=True)


def read_flattened_skymap(flatten_file_path):
    (prob, distmu, distsigma, distnorm), hdr = hp.fitsfunc.read_map(flatten_file_path, verbose=False, partial=False,
                                                                    field=range(4), h=True)

    return (prob, distmu, distsigma, distnorm), hdr


def plot_skymap(prob):
    fig = plt.figure()
    hp.mollview(prob, fig=fig, cmap='Oranges')
    plt.savefig(r'data/mollweide.png', bbox_inches='tight')