from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from calibrate_transport import transport_sensitivity,discharge_voltage

def test_fvm_ratio_plus_flow_increases_conditional_limiting_current():
    d=transport_sensitivity()
    s=d[(d.species=='Ti3')&(d.D_scale==1)].sort_values('Q_mL_min')
    p=[1.7,3.9,280.]
    vals=[]
    for _,row in s.iterrows():
        m=row.k_channel_m_s/s.iloc[1].k_channel_m_s
        v=discharge_voltage(np.array([150.]),[p[0],p[1],p[2]*m])[0]
        vals.append(v)
    assert np.all(np.diff(vals)>0)

def test_fixed_j_flux_per_pass_decreases_when_flow_increases():
    F=96485.33212
    for flow in [10,20,40]:
        utilization=.9/F/(flow*1e-6/60*500.)
        assert utilization>0
    assert .9/F/(40e-6/60*500.) < .9/F/(10e-6/60*500.)
