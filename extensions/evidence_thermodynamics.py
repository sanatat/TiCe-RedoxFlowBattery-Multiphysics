"""Source-traceable solubility evidence and *medium-limited* Ce sulfate speciation.

This is NOT a validated precipitation-equilibrium model for Ti/Ce in MSA/H2SO4.
It does not calculate Ksp, precipitation rates, mineral identity, or a
saturation index. It distinguishes empirical Ce solubility observations
from conditional sulfate complex formation at fixed ionic medium (2 M NaClO4).

Sources (retrieved 2026-11-11):
 * Patent US12119527B2, Table 1: Ce oxidation-state-specific solubility in
   MSA/H2SO4, temperature and equilibration protocol not specified.
   https://patents.google.com/patent/US12119527B2/en
 * Xie et al. Energy & Fuels (2011), DOI 10.1021/ef200354b: 1 M Ce
   mixed-acid stability; valence-specific interpretation not established by abstract.
 * Rare Earth Elements Equilibria in Aqueous Media, Table 11.2:
   conditional Ce/SO4 log beta, 25 C, background 2 M NaClO4 only.
   Background ion medium CANNOT be transferred directly to 2 M MSA.
 * Ti(IV) hydrolysis/speciation discrepancies: DOI 10.1039/D1RA04284J;
   temperature dependences DOI 10.1016/j.jct.2018.11.001.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from math import isfinite
from typing import Optional

PATENT_URI = "https://patents.google.com/patent/US12119527B2/en"
ACS_URI = "https://doi.org/10.1021/ef200354b"
THERMO_REVIEW_URI = "https://pubs.rsc.org/en/content/articlehtml/2021/ra/d1ra04284j"
SULFATE_URI = "https://ftp.idu.ac.id/wp-content/uploads/ebook/tdg/BUKU%20REE/Rare%20Earth%20Elements%20in%20Human%20and%20Environmental%20Health%20At%20the%20Crossroads%20Between%20Toxicity%20and%20Safety%20by%20Pagano%2C%20Giovanni%20%28z-lib.org%29.pdf"


@dataclass(frozen=True)
class SolubilityObservation:
    msa_M: float
    h2so4_M: float
    ce3_limit_M: float
    ce4_limit_M: Optional[float]
    ce4_qualitative: Optional[str] = None
    method: str = "patent_table_1"
    temperature_C: Optional[float] = None
    source: str = PATENT_URI


# Exact acid compositions from patent Table 1. NO interpolated/extrapolated values.
PATENT_SOLUBILITY = (
    SolubilityObservation(4.0, 0.00, 0.9, 0.9),
    SolubilityObservation(3.0, 0.25, 0.9, 0.5),
    SolubilityObservation(2.0, 0.50, 0.8, 0.2),
    SolubilityObservation(1.0, 0.75, 0.7, None, "reported_insoluble_no_numeric_limit"),
    SolubilityObservation(0.0, 1.00, 0.5, 0.5),
)


def lookup_exact_mixture(msa_M: float, h2so4_M: float, tolerance: float = 1e-8):
    """Return evidence ONLY for measured acid combinations, otherwise None."""
    if not all(isfinite(x) and x >= 0 for x in (msa_M, h2so4_M)):
        raise ValueError("Acid concentrations must be finite and nonnegative")
    if not isfinite(tolerance) or tolerance < 0 or tolerance > 1e-3:
        raise ValueError("Tolerance must be nonnegative and at most 1e-3 M")
    for point in PATENT_SOLUBILITY:
        if abs(msa_M-point.msa_M) <= tolerance and abs(h2so4_M-point.h2so4_M) <= tolerance:
            return point
    return None


def _limit_assessment(concentration_M: float, limit_M: Optional[float], qualifier=None):
    if limit_M is None:
        return {"status":"qualitative_only_no_quantitative_limit", "reported_limit_M":None,
                "margin_M":None, "ratio_to_reported_limit":None,
                "qualification":qualifier}
    return {"status":("above_reported_limit" if concentration_M > limit_M+1e-10
                       else "at_reported_limit" if abs(concentration_M-limit_M)<=1e-10
                       else "below_reported_limit"),
            "reported_limit_M":limit_M,
            "margin_M":limit_M-concentration_M,
            "ratio_to_reported_limit":concentration_M/limit_M,
            "qualification":None}


def screen_ce_evidence(msa_M: float, h2so4_M: float, ce_total_M: float,
                       ce4_fraction: float, temperature_C: Optional[float] = None) -> dict:
    """Check Ce(III/IV) against OBSERVED limits for exact acid mixtures.

    ce4_fraction is the Ce oxidation-state fraction (not necessarily full-cell SOC).
    It is NOT a precipitation-onset prediction, Ksp, or chemical equilibrium solver.
    Patent measurements give no verified temperature, salt preparation or equilibration.
    """
    if not isfinite(ce_total_M) or ce_total_M < 0:
        raise ValueError("Ce loading must be finite and nonnegative")
    if not isfinite(ce4_fraction) or not (0<=ce4_fraction<=1):
        raise ValueError("Ce(IV) fraction must be within [0,1]")
    if temperature_C is not None and not isfinite(temperature_C):
        raise ValueError("Temperature must be finite")
    ce3 = ce_total_M*(1-ce4_fraction)
    ce4 = ce_total_M*ce4_fraction
    obs=lookup_exact_mixture(msa_M,h2so4_M)
    if obs is None:
        return {"evidence_status":"not_measured_in_published_table", "msa_M":msa_M,
                "h2so4_M":h2so4_M,"ce_total_M":ce_total_M,"ce4_fraction":ce4_fraction,
                "ce3_M":ce3,"ce4_M":ce4,"temperature_C":temperature_C,
                "ce3":None,"ce4":None,"is_thermodynamic_saturation_index":False,
                "assessment":"No defensible solubility prediction at this composition. No interpolation."}
    a3=_limit_assessment(ce3, obs.ce3_limit_M)
    a4=_limit_assessment(ce4, obs.ce4_limit_M,obs.ce4_qualitative)
    return {"evidence_status":"patent_observational_limits_unverified_temperature",
            "msa_M":msa_M,"h2so4_M":h2so4_M,"ce_total_M":ce_total_M,
            "ce4_fraction":ce4_fraction,"ce3_M":ce3,"ce4_M":ce4,
            "temperature_C":temperature_C,"ce3":a3,"ce4":a4,
            "source":obs.source,"is_thermodynamic_saturation_index":False,
            "assessment":"Evidence-based flag at an exactly reported acid mixture; NOT a prediction of a solid phase, precipitation onset, or long-term safety. Temperature and ionic activities not established."}


# LOG BETA are *conditional* formation constants in 2 M NaClO4 at 25 C.
# They are not portable to MSA/H2SO4 and must not be used for battery precipitation.
CE_SULFATE_LOG_BETA_2M_NACLO4_25C = {
    "CeIII_SO4": 1.24,
    "CeIV_SO4": 4.62,
    "CeIV_2SO4": 8.00,
    "CeIV_3SO4": 10.38,
}


def conditional_ce_sulfate_fractions(free_so4_M:float, medium:str="2M_NaClO4_25C") -> dict:
    """Illustrate literature-supported *conditional complex fractions*, NOT Ksp.

    Free sulfate input is supplied; this function does not solve HSO4 acidity
    or the SO4 material balance. Valid medium is explicitly the one measured.
    """
    if medium!="2M_NaClO4_25C":
        raise ValueError("Literature log beta cannot be transferred to other acid/ionic media")
    if not isfinite(free_so4_M) or free_so4_M<0:
        raise ValueError("Free sulfate must be finite/nonnegative")
    b=CE_SULFATE_LOG_BETA_2M_NACLO4_25C
    w3=10**b['CeIII_SO4']*free_so4_M
    w41=10**b['CeIV_SO4']*free_so4_M
    w42=10**b['CeIV_2SO4']*free_so4_M**2
    w43=10**b['CeIV_3SO4']*free_so4_M**3
    s3=1+w3;s4=1+w41+w42+w43
    return {"free_SO4_M":free_so4_M,"medium":medium,"temperature_C":25.0,
            "CeIII_free_fraction":1/s3,"CeIII_SO4_fraction":w3/s3,
            "CeIV_free_fraction":1/s4,"CeIV_SO4_fraction":w41/s4,
            "CeIV_2SO4_fraction":w42/s4,"CeIV_3SO4_fraction":w43/s4,
            "is_precipitation_prediction":False}


def equilibrium_data_gap_audit() -> dict:
    """Publish explicit criteria for a genuine thermodynamic precipitation model."""
    return {"completed":[
        "Oxidation-state-resolved measured Ce solubility evidence for five exact nominal acid mixtures (patent)",
        "Ce(III)/Ce(IV) sulfate conditional complexation at 25 C in 2 M NaClO4 (NOT mixed acid)",
        "Qualitative acid/temperature/SOC stability constraints with published references"],
        "missing":[
        "Reliable Ce(III/IV) and Ti(III/IV) methanesulfonate/sulfate/hydrolysis formation constants in the target mixed-acid ionic medium",
        "Verified identity and Ksp/chemical potentials of all precipitating solid phases across relevant acidity and oxidation state",
        "Validated high-ionic-strength activity model (Pitzer/SIT) for mixed MSA/H2SO4/Ti/Ce, temperature dependent",
        "Temperature-dependent solubility measurements for both oxidation states at controlled SOC",
        "Precipitation kinetics, induction times, and metastable versus equilibrium boundaries",
        "Raw experimental acid/free sulfate/pH/temperature data; co-dissolved Ti effects"],
        "scientific_verdict":"Partial evidence-based screening only; quantitative Ti/Ce mixed-acid saturation indices remain unvalidated"}
