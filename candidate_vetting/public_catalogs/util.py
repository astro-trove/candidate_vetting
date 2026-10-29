"""
Some useful variables that will be used throughout this entire directory
"""

import numpy as np

from django.db.models import (
    Func,
    BooleanField,
    FloatField,
    ExpressionWrapper,
)
from django.db.models import F, Value
from django.db.models.functions import Pi, Exp, Power

#######################################
#Dependencies for _cluster_galaxy_check
import logging
import pandas as pd
import numpy as np                      
import astropy.units as u
from astropy.coordinates import angular_separation
from tom_targets.models import Target, TargetExtra
from django.db.models import F, Case, When, Value, FloatField  

logger = logging.getLogger(__name__)
RH_TO_KRON = 2.5
SHAPE_R_FLOOR_KPC = 5.0
ARCSEC_PER_RADIAN = 206265.0
EVCC_DF_COLMAP = {
    "name": "ID", "offset_arcsec": "Offset", "offset_kron": "OffsetKron",
    "rad": "Kron", "r50": "R50", "rmag": "rMag", "morphology": "Type", "cz_kms": "cz",
    }
FORNAX_DORADO_DF_COLMAP = {
    "name": "ID", "offset_arcsec": "Offset", "offset_kron": "OffsetKron",
    "eff_shape_r": "Rh", "r_kron": "rKron", "dered_mag_r": "rMag",
    "type": "Type", "radvel": "radvel", "dist_mpc": "Dist_Mpc",
    }
#######################################


RADIUS_ARCSEC = 2.0


## threshold for PS1 point sources
# 0.83, as used is Rastinejad+22 and based on Tachibana & Miller 2018
# morphological classification
# yields true positive rate of ~0.995, false positive rate of ~0.005
PS1_TB18_POINT_SOURCE_THRESHOLD = 0.83


## Beck+21 boundary for STRM classification of PS1 stars, galaxies, and QSOs
PS1_B21_DECISION_BOUNDARY = 0.7


_QUERY_METHOD_DOCSTRING = f"""
Parameters
----------
ra : float
    The RA in degrees
dec : float
    The declination in degrees
radius : float, default={RADIUS_ARCSEC}
    The search radius in arcseconds

Returns
-------
QueryResult
    The result of the query, either a Galaxy or PointSource class
"""


class ConeSearch(Func):
    function = "q3c_join"
    arity = 5
    output_field = BooleanField()


class AngDist(Func):
    function = "q3c_dist"
    arity = 4
    output_field = FloatField()


class Greatest(Func):
    function = "GREATEST"
    arity = 2
    output_field = FloatField()


def cone_search_q3c(queryset, ra, dec, radius=RADIUS_ARCSEC, ra_colname="ra", dec_colname="dec"):
    f"""Do a cone search with q3c on the provided QuerySet

    Parameters
    ----------
    queryset : django.QuerySet
        A query set to filter with the cone search
    ra : float
        The RA in degrees
    dec : float
        The dec in degrees
    radius : float, default={RADIUS_ARCSEC}
        The search radius in arcseconds
    ra_colname: str, default='ra'
        The RA column name in the SQL table
    dec_colname: str, default="dec"
        The dec column name in the SQL table

    Returns
    -------
    django.QuerySet
        The filtered query set
    """

    return queryset.filter(ConeSearch(ra, dec, ra_colname, dec_colname, radius / 3600)).annotate(
        ang_dist=AngDist(ra_colname, dec_colname, ra, dec)
    )


def pcc_q3c(queryset, ra, dec, pcc_max, mag_colname, mag_ab_offset=0, ra_colname="ra", dec_colname="dec"):
    """Do a cut on Pcc with q3c on the provided QuerySet. Pcc is from Bloom+2002 and
    we use the re-calibration from Berger2010.

    Parameters
    ----------
    queryset : django.QuerySet
        A query set to filter with the cone search
    ra : float
        The RA in degrees
    dec : float
        The dec in degrees
    pcc_max : float
        The maximum Pcc to consider
    mag_colname : str
        The magnitude column name to use for Pcc
    mag_ab_offset : float
        The magnitude system offset from AB. For example, if the catalog provides J-band
        magnitudes in the Vega system, you would pass in 0.91. 
    ra_colname: str, default='ra'
        The RA column name in the SQL table
    dec_colname: str, default="dec"
        The dec column name in the SQL table

    Returns
    -------
    django.QuerySet
        The filtered query set
    """
    pcc_prefactor = 1 / (0.33 * np.log(10))
    return queryset.annotate(
        ang_dist=AngDist(ra_colname, dec_colname, ra, dec),
        pcc_sigma=ExpressionWrapper(
            pcc_prefactor
            * Power(
                10,
                Greatest(0.33 * ((F(mag_colname) + mag_ab_offset) - 24) - 2.44, Value(-10)),
            ),
            output_field=FloatField(),
        ),
        pcc=ExpressionWrapper(
            1
            - Exp(
                Greatest(
                    -1 * Pi() * Power(3600 * F("ang_dist"), 2) * F("pcc_sigma"),
                    -700,  # anything less than e^-700 will be zero anyways!
                ),
            ),
            output_field=FloatField(),
        ),
    ).filter(pcc__lte=pcc_max)


def citation(doi="", ads_bibcode="", version=None, data_url=""):
    """
    A citation decorator for a catalog class to allow us to keep track
    """

    def decorator(cls):

        doi_list = doi
        if not isinstance(doi, list):
            doi_list = [doi]

        ads_bibcode_list = ads_bibcode
        if not isinstance(ads_bibcode, list):
            ads_bibcode_list = [ads_bibcode]

        cls.dois = [d for d in doi_list]
        cls.bibcodes = [b for b in ads_bibcode_list]

        base_ads_url = "https://ui.adsabs.harvard.edu/abs/"
        cls.ads_urls = [base_ads_url + b for b in cls.bibcodes]

        cls.catalog_version = version
        cls.data_url = data_url

        return cls

    return decorator


def _cluster_galaxy_check(
    target_id: int,
    *,
    model,
    key: str,
    cluster: str,
    dist_mpc: float,
    center_ra: float,
    center_dec: float,
    bounding_radius: float,
    max_kron_arcsec: float,
    kron_scale: float = 2.0,
) -> str:
    """
    Shared geometric association against EVCC or LS DR10 cluster catalogs. A target is
    'within' a galaxy when its angular offset < kron_scale * the galaxy's Kron radius (rad) for EVCC 
    and angular offset < kron_scale * r_kron, where r_kron = RH_TO_KRON * shape_r (half-light radius), 
    and shape_r <= 0 (PSF-fit galaxies) is floored to SHAPE_R_FLOOR_KPC (for LS DR10). Writes matches 
    to the given TargetExtra `key` as JSON records and returns that JSON string ("[]" when none). 
    Used by evcc_galaxy_check, fornax_galaxy_check, and dorado_galaxy_check.
    Some of the fixed constants and the column names of EVCC, Fornax+Dorado are below
    """
    target = Target.objects.get(id=target_id)
    ra, dec = target.ra, target.dec

    # cheap footprint guard: no DB query for the vast majority of sky.
    cone_radius  = kron_scale * max_kron_arcsec
    guard_radius = bounding_radius + cone_radius / 3600.0

    sep_deg = (
        angular_separation(
            ra * u.deg,
            dec * u.deg,
            center_ra * u.deg,
            center_dec * u.deg,
        )
        .to(u.deg)
        .value
    )
    if sep_deg > guard_radius:
        TargetExtra.objects.update_or_create(
            target=target,
            key=key,
            defaults={"value": "[]"},
        )
        return "[]"

    if cluster == 'EVCC':
        qs = cone_search_q3c(model.objects.exclude(rad=None), ra, dec,
                             radius=cone_radius, ra_colname="ra", dec_colname="dec"
             ).filter(ang_dist__lt=kron_scale * F("rad") / 3600.0)
        if not qs.exists():
            TargetExtra.objects.update_or_create(target=target, key=key, defaults={"value": "[]"})
            return "[]"
        df = pd.DataFrame(qs.values("evcc", "vcc", "ngc", "ang_dist", "rad", "r50", "rmag",
                "pmorph", "srvel"))

        df["offset_arcsec"] = (df["ang_dist"] * 3600.0).round(2)
        df["offset_kron"]   = (df["ang_dist"] * 3600.0 / df["rad"]).round(2)
        df["name"] = df.apply(
        lambda r: f"NGC {r['ngc']}" if r["ngc"]
            else f"VCC {r['vcc']}" if r["vcc"]
            else f"EVCC {r['evcc']}",
            axis=1,
        )
        df = df.rename(columns={"pmorph": "morphology", "srvel": "cz_kms"})
        df = df.sort_values("offset_arcsec")

        for _, r in df.iterrows():
            logger.info(
                f"{target.name} within {kron_scale}x Kron of {r['name']}: "
                f"offset={r['offset_arcsec']:.1f}\" ({r['offset_kron']:.2f} Kron radii)"
            )

        newdf = df[list(EVCC_DF_COLMAP.keys())].rename(columns=EVCC_DF_COLMAP)
        result = newdf.to_json(orient="records")   # pandas serializes NaN -> null
        TargetExtra.objects.filter(target_id=target.id, key=key).delete()
        TargetExtra.objects.create(target=target, key=key, value=result)
        return result

    else:
        # effective half-light radius: measured shape_r, or the SHAPE_R_FLOOR_KPC
        # physical floor (converted to arcsec at this cluster's distance) for the
        # PSF-fit galaxies (shape_r <= 0). Applied in-DB so the containment filter
        # and the cone both use the floored value.
        rh_floor_arcsec = SHAPE_R_FLOOR_KPC / (dist_mpc * 1e3) * ARCSEC_PER_RADIAN
        eff_shape_r = Case(When(shape_r__gt=0, then=F("shape_r")),
                           default=Value(rh_floor_arcsec), output_field=FloatField())
        qs = cone_search_q3c(model.objects.all(), ra, dec,
                             radius=cone_radius, ra_colname="ra", dec_colname="dec"
             ).annotate(eff_shape_r=eff_shape_r
             ).filter(ang_dist__lt=kron_scale * RH_TO_KRON * F("eff_shape_r") / 3600.0)
        if not qs.exists():
            TargetExtra.objects.update_or_create(target=target, key=key, defaults={"value": "[]"})
            return "[]"

        df = pd.DataFrame(qs.values("name", "ang_dist", "eff_shape_r", "type", "dered_mag_r", "radvel"))

        df["r_kron"] = RH_TO_KRON * df["eff_shape_r"]
        df["offset_arcsec"] = df["ang_dist"] * 3600.0
        df["offset_kron"] = (df["offset_arcsec"] / df["r_kron"]).round(2)

        df["r_kron"] = df["r_kron"].round(2)
        df["offset_arcsec"] = df["offset_arcsec"].round(2)

        df["radvel"] = df["radvel"].replace(0, np.nan)  # 0 = missing
        df["dist_mpc"] = dist_mpc
        df = df.sort_values("offset_arcsec")

        for _, r in df.iterrows():
            logger.info(
                f"{target.name} within {kron_scale}x r_kron of {cluster} {r['name']}: "
                f'offset={r["offset_arcsec"]:.1f}" ({r["offset_kron"]:.2f} r_kron)'
            )

        newdf = df[list(FORNAX_DORADO_DF_COLMAP.keys())].rename(columns=FORNAX_DORADO_DF_COLMAP)
        result = newdf.to_json(orient="records")  # pandas serializes NaN -> null
        TargetExtra.objects.filter(target_id=target.id, key=key).delete()
        TargetExtra.objects.create(target=target, key=key, value=result)
        return result