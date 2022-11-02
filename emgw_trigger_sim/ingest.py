import healpy as hp
import numpy as np
from ligo.skymap.postprocess import find_greedy_credible_levels
from scipy.stats import norm
from utils import dist_function
import matplotlib.pyplot as plt
import pickle
from astropy.io import fits
from astropy.table import Table
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
import os
from astropy.io import ascii


def return_distsamples(Rlow_mpc, Rhigh_mpc, distnorm, distmu, distsigma):
    # least bright is in bronze sample at -11
    # brightest literature is -16

    distsamples = []
    Rs = np.linspace(Rlow_mpc, Rhigh_mpc, 100)
    dp_drs = dist_function(Rs, distnorm, distmu, distsigma)
    #     print(dp_drs)
    scale = 1e4 / np.median(dp_drs)
    for ind, r in enumerate(Rs):
        distsamples.extend([r for j in range(int(scale * dp_drs[ind]))])

    distsamples = np.array(distsamples)
    np.random.shuffle(distsamples)

    return distsamples


def get_2d_pixel_samples(prob):
    credible_levels = find_greedy_credible_levels(prob)
    npix = len(prob)
    nside = hp.npix2nside(npix)

    prob_vals = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
    Ntot = 1000

    all_pixs = []

    for ind in range(len(prob_vals) - 1):
        low_lim = prob_vals[ind]
        up_lim = prob_vals[ind + 1]

        pixnums = np.where((credible_levels < up_lim) & (credible_levels < up_lim))[0]
        num_samples = round(Ntot * (up_lim - low_lim))
        #         all_pixs.append(list(pixnums)*num_samples)
        randinds = np.random.randint(0, len(pixnums), num_samples)
        all_pixs.append(pixnums[randinds])

    weighted_pixels = np.concatenate(all_pixs)

    return weighted_pixels


def get_event_coordinates(prob, distnorm, distmu, distsigma, r_low_mpc=0, r_high_mpc=4000, num_trans=1):
    npix = len(prob)
    nside = hp.npix2nside(npix)

    pixel_samples = get_2d_pixel_samples(prob)
    randindex = np.random.randint(0, len(pixel_samples), num_trans)

    event_pixel = pixel_samples[randindex]

    theta, phi = hp.pix2ang(nside, event_pixel)
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)
    print(ra,dec)
    rs = np.linspace(r_low_mpc, r_high_mpc, 100)
    dp_dr = np.array(
        [r ** 2 * distnorm[event_pixel] * norm(distmu[event_pixel], distsigma[event_pixel]).pdf(r) for r in rs]).T

    scale = 1e3 / np.median(dp_dr)

    sums = np.array([np.sum(dp_dr[ind]) for ind in range(len(dp_dr))])
    print(sums)
    dp_dr = dp_dr[sums > 0]
    ra = ra[sums > 0]
    dec = dec[sums > 0]
    distances_mpc = np.array([np.random.choice(rs, p=dp_dr[ind] / np.sum(dp_dr[ind])) for ind in range(len(dp_dr))])

    return ra, dec, distances_mpc


def get_event_coordinates_loop(prob, r_low_mpc=0, r_high_mpc=4000, num_trans=1):
    npix = len(prob)
    nside = hp.npix2nside(npix)

    randindex = []
    for i in range(num_trans):
        pixel_samples = get_2d_pixel_samples(prob)
        randindex.append(np.random.randint(0, len(pixel_samples), 1))

    event_pixel = pixel_samples[randindex]

    theta, phi = hp.pix2ang(nside, event_pixel)
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)

    distances_mpc = []
    for pix in event_pixel:
        distance_samples = return_distsamples(r_low_mpc, r_high_mpc, distnorm[pix], distmu[pix], distsigma[pix])
        randindex = np.random.randint(len(distance_samples))
        distance_mpc = distance_samples[randindex]

        distances_mpc.append(distance_mpc)
    distances_mpc = np.array(distances_mpc)
    return ra, dec, distances_mpc