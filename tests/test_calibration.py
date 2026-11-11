"""New tests: computational consistency, not proof of a real experimental digital twin."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from calibrate_transport import (discharge_voltage,load_digitized,fit_reduced,
  transport_sensitivity,effective_transport_envelope,felt_structural_sensitivity,FARADAY,
  MEMBRANE_ASR_OHM_CM2)

@pytest.fixture(scope='module')
def dat():
    d=load_digitized();j=d.j_mA_cm2.to_numpy();v=d.V_discharge_approx_V.to_numpy()
    mask=~d.held_out_for_assessment.to_numpy()
    p=fit_reduced(j,v,mask)
    return d,j,v,mask,p

@pytest.fixture(scope='module')
def tx():return transport_sensitivity()

def test_digitization_clearly_tagged(dat):
    d,*_=dat
    assert d.source.str.contains('visual estimate').all()
    assert (d.approx_digitization_uncertainty_V>.02).all()

def test_calibration_holdout_is_not_train(dat):
    d,j,v,train,p=dat
    assert np.sum(~train)==4 and len(j)==14 and np.all(np.diff(j)>0)
    assert np.sqrt(np.mean((discharge_voltage(j[~train],p)-v[~train])**2))<.06

def test_reduced_model_voltage_monotone(dat):
    _,j,v,train,p=dat
    a=discharge_voltage(np.linspace(0,275,100),p)
    assert np.all(np.diff(a)<0)
    assert np.isfinite(a).all()

def test_no_physical_claim_made_on_surface_area(dat):
    _,j,v,train,p=dat
    # A lumped 3-param fit does not report specific electroactive area as a parameter.
    assert len(p)==3
    assert p[2]>j.max()

def test_fvm_conservation_for_all_conditions(tx):
    assert len(tx)==18
    assert (tx.mass_balance_relative_error.abs()<1e-9).all()
    assert (tx.k_channel_m_s>0).all()

def test_fvm_flow_sensitivity(tx):
    for spec in ('Ti3','Ce4'):
        s=tx[(tx.species==spec)&(tx.D_scale==1)].sort_values('Q_mL_min')
        assert np.all(np.diff(s.k_channel_m_s)>0)

def test_fvm_diffusion_sensitivity(tx):
    s=tx[(tx.species=='Ti3')&(tx.Q_mL_min==20)].sort_values('D_scale')
    assert np.all(np.diff(s.k_channel_m_s)>0)

def test_transport_envelope_soc_and_dimensions(dat,tx):
    _,j,v,train,p=dat
    env=effective_transport_envelope(p,tx)
    t=env[(env.species=='Ti3')&(env.D_scale==1)].sort_values('SOC_assumed')
    assert np.all(np.diff(t.k_required_m_s)<0)
    expected=p[2]*10/(FARADAY*250.)
    assert np.isclose(t[t.SOC_assumed==.5].iloc[0].k_required_m_s,expected)

def test_felt_area_thiele_trends():
    t=felt_structural_sensitivity()
    a=t[(t.D_apparent_multiplier==1)&(t.intrinsic_k_surface_m_s_assumed==1e-6)]
    assert np.all(np.diff(a.sort_values('a_s_m2_m3_assumed').k_felt_proxy_m_s)>0)
    assert (a.eta_felt>0).all() and (a.eta_felt<=1).all()

def test_asr_real_value():
    assert np.isclose(MEMBRANE_ASR_OHM_CM2,0.0942)
    # At j=90 mA/cm2 membrane drop is <10 mV, NOT the entire discharge voltage loss.
    assert 90*.001*MEMBRANE_ASR_OHM_CM2<.01
