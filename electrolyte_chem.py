"""Source-derived electrolyte coordination energy tables and ideal proton sensitivity.

These electronic binding energies are NOT conditional Gibbs formation constants;
they cannot predict Ti/Ce speciation, dissolution/precipitation, or mixed acids.
"""
import numpy as np
R=8.314462618;F=96485.33212
# SI pp 30-31 Table S2, atomic Hartree units; ligand count x.
S2={
 'Ce3':{'ligands':[0,1,2,3,4,5,6,7,8],
        'binding_Ha':[-.15502,-.32603,-.48887,-.61594,-.69434,-.77739,-1.44978,-1.50243,-1.56718]},
 'Ce4':{'ligands':[0,1,2,3,4,5,6,7,8,9,10],
        'binding_Ha':[-.14266,-.32178,-.47149,-.59553,-.66155,-.74229,-.81383,-.84463,-.89600,-.95216,-.98824]},
 'Ti3':{'ligands':[0,1,2,3,4,5,6],
        'binding_Ha':[-.10948,-.28154,-.43616,-.55578,-.62292,-.69152,-.74768]},
 'Ti4_binuclear':{'ligands':[0,2,4,6,8],
        'binding_Ha':[-.32412,-.62123,-.73081,-.70035,-.71578]}}


def formal_cell_voltage_ideal_proton(soc=.5,T=298.15,relative_H_activity=1.0):
    """Formal Nernst for Ce(IV)/Ce(III) and TiO2+ +2H+ +e- <-> Ti3+ +H2O.

    Assumes Ce/Ti redox active concentrations same SOC and references proton
    activity at 1.0 relative, *not* pH measured, and ignores complexation.
    """
    if not 0<soc<1 or relative_H_activity<=0:raise ValueError('positive H, 0<SOC<1')
    base=1.51+2*R*T/F*np.log(soc/(1-soc))
    return base-2*R*T/F*np.log(relative_H_activity)
