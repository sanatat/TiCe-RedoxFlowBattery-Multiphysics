"""Numerical / provenance controls; NOT experimental predictive validation."""
from math import isclose
import pytest
from extensions.evidence_thermodynamics import (
    lookup_exact_mixture,screen_ce_evidence,conditional_ce_sulfate_fractions,
    CE_SULFATE_LOG_BETA_2M_NACLO4_25C,equilibrium_data_gap_audit)


def test_exact_record_no_interpolation():
    p=lookup_exact_mixture(2,.5)
    assert p.ce3_limit_M==.8 and p.ce4_limit_M==.2
    assert lookup_exact_mixture(2.5,.4) is None


def test_ce_oxidation_state_inventory():
    for frac in (0,.25,.5,.9,1):
        r=screen_ce_evidence(2,.5,.5,frac)
        assert r['ce3_M']+r['ce4_M']==pytest.approx(.5)
        assert r['ce4_M']==pytest.approx(frac*.5)
        assert r['is_thermodynamic_saturation_index'] is False


def test_ce_iv_risk_at_exact_mixture():
    r=screen_ce_evidence(2,.5,.5,.5)
    assert r['ce3']['status']=='below_reported_limit'
    assert r['ce4']['status']=='above_reported_limit'
    assert r['ce4']['ratio_to_reported_limit']==pytest.approx(1.25)


def test_no_extrapolation_is_unknown():
    r=screen_ce_evidence(2.0,.1,.5,.5)
    assert r['evidence_status']=='not_measured_in_published_table'
    assert r['ce4'] is None


def test_all_msa_below_reported_ce_limits():
    r=screen_ce_evidence(4,0,.5,.5)
    assert r['ce4']['status']=='below_reported_limit'
    assert r['ce3']['status']=='below_reported_limit'


def test_qualitative_insoluble_not_fake_numeric_limit():
    r=screen_ce_evidence(1,.75,.5,.8)
    assert r['ce4']['status']=='qualitative_only_no_quantitative_limit'
    assert r['ce4']['reported_limit_M'] is None


def test_negative_and_invalid_inputs_rejected():
    with pytest.raises(ValueError):screen_ce_evidence(2,.5,-.3,.5)
    with pytest.raises(ValueError):screen_ce_evidence(2,.5,.5,1.2)
    with pytest.raises(ValueError):lookup_exact_mixture(-1,2)
    with pytest.raises(ValueError):conditional_ce_sulfate_fractions(-1)


def test_conditional_ce_mass_conservation():
    for sulfate in (0,1e-7,1e-5,.001,.1):
        x=conditional_ce_sulfate_fractions(sulfate)
        sum3=x['CeIII_free_fraction']+x['CeIII_SO4_fraction']
        sum4=sum(x[k] for k in ('CeIV_free_fraction','CeIV_SO4_fraction',
                               'CeIV_2SO4_fraction','CeIV_3SO4_fraction'))
        assert sum3==pytest.approx(1,abs=1e-12)
        assert sum4==pytest.approx(1,abs=1e-12)


def test_conditional_speciation_correct_order():
    lo=conditional_ce_sulfate_fractions(1e-7)
    hi=conditional_ce_sulfate_fractions(1e-3)
    assert lo['CeIV_free_fraction']>hi['CeIV_free_fraction']
    assert hi['CeIV_2SO4_fraction']>lo['CeIV_2SO4_fraction']
    assert CE_SULFATE_LOG_BETA_2M_NACLO4_25C['CeIV_2SO4']==8.00


def test_forbid_medium_transfer():
    with pytest.raises(ValueError):conditional_ce_sulfate_fractions(.01,'2M_MSA')


def test_audit_does_not_claim_precipitation_equilibrium():
    a=equilibrium_data_gap_audit()
    assert len(a['missing'])>=5
    assert 'unvalidated' in a['scientific_verdict']
