"""
Vetting code for non-localized events with transients
"""

import time
import io
import logging

import numpy as np
import pandas as pd

from astropy.utils.introspection import minversion

from astropy import units as u
from astropy.time import Time

import sqlalchemy as sa
from sqlalchemy.orm import Session

from django.conf import settings

from trove_mpc import Transient
from tom_targets.models import Target, TargetExtra
from tom_nonlocalizedevents.models import (
    # EventCandidate,
    EventLocalization,
    # SkymapTile,
    NonLocalizedEvent,
)
from tom_nonlocalizedevents.healpix_utils import (
    sa_engine,
    SaSkymapTile,
    # uniq_to_bigintrange,
    # update_all_credible_region_percents_for_candidates
)
from tom_dataproducts.models import ReducedDatum

from candidate_vetting.public_catalogs.static_catalogs import (
    # DesiSpec,
    Cosmicflows4,
    GladePlus,
    Gwgc,
    Hecate1,
    Hecate2,
    LsDr9North,
    LsDr10South,
    Ps1Galaxy,
    Sdss12Photoz,
    AsassnVariableStar,
    Gaiadr3Variable,
    ZtfVarStar,
    Ps1PointSource,
    Milliquas,
    NedLvs,
    # TwoMass,
    DesiDr1,
    ExtendedVirgoClusterCatalog,
)
from candidate_vetting.public_catalogs.dynamic_catalogs import UserGalaxy

if minversion(np, "2.0.0"):
    np_trapz_fn = np.trapezoid
else:
    np_trapz_fn = np.trapz  # np.trapz is deprecated in numpy >2.0.0

cosmo = settings.COSMO
logger = logging.getLogger(__name__)

HOST_DF_COLMAP = {
    "trove_uniq": "troveID",
    "name": "ID",
    "pcc": "PCC",
    "offset": "Offset",
    "ra": "RA",
    "dec": "Dec",
    "lumdist": "Dist",
    "lumdist_err": "DistErr",
    "z": "z",
    "z_err": "zErr",
    "z_type": "z_type",
    "default_mag": "Mags",
    "catalog": "Source",
    "submitter": "Submitter",
}
HOST_DF_COLMAP_INVERSE = {v: k for k, v in HOST_DF_COLMAP.items()}

# Host, point source, and AGN association radii
HOST_ASSOC_RADIUS = 5 * 60  # 5 arcmin = 300 arcsec, as used in Franz+25 and Vieira+26
PS_ASSOC_RADIUS = 2  # 2 arcsec, as used in Franz+25 and Vieira+26
AGN_ASSOC_RADIUS = 2  # 2 arcsec, as used in Franz+25 and Vieira+26

# After we order the dataframe by the Pcc score, remove any host matches with a greater
# Pcc score than this
PCC_THRESHOLD = 0.15  # this is the value used in Rastinejad+2022

# rank order of the galaxy catalogs for getting the "default" distance to this transient
# this is kinda arbitrary, but generally I consider
# 1) is this a redshift catalog or a galaxy distance catalog? An actual galaxy distance
#    catalog is preferred over a general redshift catalog
# 2) Does this catalog have spec-z's or photo-z's? A spec-z catalog is preferred.
GALAXY_CATALOGS = [
    UserGalaxy,
    ExtendedVirgoClusterCatalog,
    GladePlus,
    Gwgc,
    Hecate2,
    DesiDr1,
    NedLvs,
    Cosmicflows4,
    LsDr9North,
    LsDr10South,
    Ps1Galaxy,
    Sdss12Photoz,
]


def localization_sequence_from_name(nonlocalized_event_name):

    nle = NonLocalizedEvent.objects.get(event_id=nonlocalized_event_name)

    seqs = nle.sequences.all()

    latest_seq = seqs[0]
    for seq in seqs:
        curr_latest_time = Time(latest_seq.details["time"])
        test_latest_time = Time(seq.details["time"])
        if test_latest_time > curr_latest_time:
            latest_seq = seq

    return seq


def save_score_to_targetextra(target, key, score):
    """
    Saves the scores that don't change to a TargetExtra object rather than a ScoreFactor
    This is for:
    1. point source score
    2. MPC score
    Since they are independent of the NLE that we are vetting the target against
    """

    # first delete the host galaxy key for this target if it already exists
    te = TargetExtra.objects.filter(target_id=target.id, key=key)
    if te.exists():
        te.delete()

    # then save the new score
    TargetExtra.objects.update_or_create(target=target, key=key, value=score)


def _save_host_galaxy_df(df, target):

    # first delete the host galaxy key for this target if it already exists
    if TargetExtra.objects.filter(target_id=target.id, key="Host Galaxies").exists():
        TargetExtra.objects.filter(target_id=target.id, key="Host Galaxies").delete()

    newdf = df[
        [
            "trove_uniq",
            "name",
            "pcc",
            "offset",
            "ra",
            "dec",
            "lumdist",
            "z",
            "z_type",
            "default_mag",
            "catalog",
            "submitter",
        ]
    ].copy()
    newdf["z_err"] = [
        [neg, pos]
        if neg != pos  # errors are asymmetric
        else neg  # errors are not assymetric
        for neg, pos in zip(df.z_neg_err, df.z_pos_err)
    ]
    newdf["lumdist_err"] = [
        [neg, pos]
        if neg != pos  # errors are asymmetric
        else neg  # errors are not assymetric
        for neg, pos in zip(df.lumdist_neg_err, df.lumdist_pos_err)
    ]
    newdf = newdf.rename(columns=HOST_DF_COLMAP)
    TargetExtra.objects.update_or_create(target=target, key="Host Galaxies", value=newdf.to_json(orient="records"))


def _save_associated_agn_df(df, target):

    # first delete the associated AGN key for this target if it already exists
    if TargetExtra.objects.filter(target_id=target.id, key="Associated AGN").exists():
        TargetExtra.objects.filter(target_id=target.id, key="Associated AGN").delete()

    col_map = {
        "name": "ID",
        # "pcc":"PCC",
        # "offset":"Offset",
        "ra": "RA",
        "dec": "Dec",
        "lumdist": "Dist",
        "lumdist_err": "DistErr",
        "z": "z",
        "z_err": "zErr",
        # "default_mag":"Mags",
        "catalog": "Source",
    }
    newdf = df[
        [
            "name",
            # "pcc",
            # "offset",
            "ra",
            "dec",
            "lumdist",
            "z",
            # "default_mag",
            "catalog",
        ]
    ]
    newdf["z_err"] = [
        [neg, pos]
        if neg != pos  # errors are asymmetric
        else neg  # errors are not assymetric
        for neg, pos in zip(df.z_neg_err, df.z_pos_err)
    ]
    newdf["lumdist_err"] = [
        [neg, pos]
        if neg != pos  # errors are asymmetric
        else neg  # errors are not assymetric
        for neg, pos in zip(df.lumdist_neg_err, df.lumdist_pos_err)
    ]
    newdf = newdf.rename(columns=col_map)
    TargetExtra.objects.update_or_create(target=target, key="Associated AGN", value=newdf.to_json(orient="records"))


def pcc(r: list[float], m: list[float]):
    """
    Probability of chance coincidence calculation (originally from
    Bloom et al. 2002 and re-calibrated in Berger2010)

    PARAMETERS
    ----------
    r : transient-galaxy offsets, array of floats
        arcseconds
    m : magnitudes of galaxies, array of floats

    RETURNS
    -------
    Pcc values : array of floats [0,1]
    """
    sigma = (1 / (0.33 * np.log(10))) * 10 ** (0.33 * (m - 24) - 2.44)
    prob = 1 - np.exp(-(np.pi * (r**2) * sigma))

    return prob


def host_association(
    target_id: int,
    radius: float = HOST_ASSOC_RADIUS,
    pcc_threshold: float = PCC_THRESHOLD,
    _verbose: bool = False,
):
    """
    Find all of the potential hosts associated with this target
    """

    target = Target.objects.filter(id=target_id)[0]
    ra, dec = target.ra, target.dec

    start = time.time()
    res = []
    for catalog in GALAXY_CATALOGS:
        cat = catalog()
        catname = str(cat)
        if _verbose:
            logger.info(f"Querying {cat}...")
        query_set = cat.pcc_filter(ra, dec, radius=radius, pcc_max=pcc_threshold)
        if _verbose:
            logger.info(f"Found {query_set.count()} matches in {catname}")

        # if no queries are returned we can skip this catalog
        if query_set.count() == 0:
            continue

        # convert to a dataframe and standardize the column names
        cols = list(cat.ogcols) + ["ang_dist", "pcc"]
        rows = query_set.values_list(*cols)
        df = pd.DataFrame.from_records(rows, columns=cols)
        df = cat.to_standardized_catalog(df)

        # some extra cleaning before continuing
        df = df.dropna(subset=["default_mag", "ra", "dec"])  # drop rows without the information we need
        df["trove_uniq"] = df["trove_uniq"].astype(int)  # set to an int

        # copy the ang_dist column to a column called "offset" for
        # backwards compatability
        # and convert to arcsec from degrees
        df["offset"] = 3600 * df.ang_dist

        # now save the cleaned dataset
        df["catalog"] = catname
        res.append(df)

    if not res:  # if no host matches
        cols = list(HOST_DF_COLMAP.keys()) + ["z_neg_err", "z_pos_err", "lumdist_neg_err", "lumdist_pos_err"]
        rows = []
        res.append(pd.DataFrame(rows, columns=cols))  # append empty df with appropriate columns

    # concact results of individual catalogs into one dataframe
    df = pd.concat(res).reset_index(drop=True)

    # TODO: We will need to put some deduplication code for the galaxy dataframe
    #       here at some point. For now it seems to work without it though!

    # sort inversely by pcc
    ret_df = df.sort_values("pcc", ascending=True)

    end = time.time()
    print(f"Galaxy table queries finished in {end - start}s")

    # save the host galaxy dataframe to the TargetExtra "Host Galaxies" keyword
    _save_host_galaxy_df(ret_df, target)
    return ret_df


def point_source_association(target_id: int, radius: float = 2):

    target = Target.objects.get(id=target_id)
    ra, dec = target.ra, target.dec

    point_source_catalogs = [
        ("source_id", AsassnVariableStar),
        ("source_id", Gaiadr3Variable),
        ("objid", Ps1PointSource),
        # ZtfVarStar,
        # this is the 2MASS point source catalog
        # I'm leaving it commented out because we need to test it a bit more before
        # using it!
        # TwoMass
    ]

    matches = {}
    for name_column, catalog in point_source_catalogs:
        cat = catalog()
        query_set = cat.query(ra, dec, radius)

        # if no matches returned, good! We can check another PS catalog
        if query_set.count() == 0:
            continue

        matches[cat.catalog_model.__name__] = (
            [ps_match.__dict__[name_column] for ps_match in query_set],
            [ps_match.ang_dist for ps_match in query_set],
        )

    return matches


def agn_association_2d(target_id: int, radius: float = AGN_ASSOC_RADIUS):
    """
    This searches the AGN catalogs for a match for this target
    """

    target = Target.objects.get(id=target_id)
    ra, dec = target.ra, target.dec

    agn_catalogs = [Milliquas]  # there is currently only one, but this should help to "future proof" the code

    agn_matches = None
    res = []
    start = time.time()
    for catalog in agn_catalogs:
        cat = catalog()
        query_set = cat.query(ra, dec, radius)

        # no match found here! let's check another catalog!
        if query_set.count() == 0:
            continue

        if agn_matches is None:
            agn_matches = query_set
        else:
            agn_matches |= query_set  # this will perform a SQL UNION on the query sets

        # convert to a dataframe and standardize the column names
        df = pd.DataFrame(list(agn_matches.values()))
        df = cat.to_standardized_catalog(df)

        # some extra cleaning before continuing
        df = df.dropna(subset=["default_mag", "ra", "dec"])  # drop rows without the information we need

        # now save the cleaned dataset
        df["catalog"] = cat.__class__.__name__
        res.append(df)

    if len(res) > 0:  # when no matches, nothing to concatenate
        df = pd.concat(res).reset_index(drop=True)
    else:  # return an empty dataframe
        return pd.DataFrame({})

    # put any more cleaning up / filtering here; none for now
    ret_df = df.copy()

    end = time.time()
    logger.info(f"AGN catalog queries finished in {end - start}s")

    # save the host galaxy dataframe to the TargetExtra "Associated AGN" keyword
    _save_associated_agn_df(ret_df, target)

    return ret_df


def run_mpc(target_id: int) -> None:

    target = Target.objects.get(id=target_id)

    # get photometry, throwing out limiting mags, phot with no error, and phot with SNR < 5
    phot = ReducedDatum.objects.filter(
        target_id=target_id,
        data_type="photometry",
        value__magnitude__isnull=False,
        value__error__isnull=False,
        value__error__lte=2.5 / np.log(10) / 5,
    )
    # if more than (5-sigma) 1 detection, likely not a MPC object
    if phot.exists() and len(phot) > 1:
        logger.warn("This candidate has more than 1 >5-sigma detection, " + "skipping MPC!")
        mpc_match = None
    # if only 1 detection, run MPC match
    elif phot.exists() and len(phot) == 1:
        latest_det = phot.latest()
        date = Time(latest_det.timestamp).mjd
        t = Transient(target.ra, target.dec)
        mpc_match = t.minor_planet_match(date)
    # no detections --> can't do MPC match
    else:
        logger.warn("This candidate has no photometry, skipping MPC!")
        mpc_match = None

    if mpc_match is not None:
        # update the score factor information
        save_score_to_targetextra(target, "mpc_match_name", mpc_match.match_name)
        save_score_to_targetextra(target, "mpc_match_sep", mpc_match.distance)
        save_score_to_targetextra(target, "mpc_match_date", latest_det.timestamp)
    else:
        save_score_to_targetextra(target, "mpc_match_name", None)
