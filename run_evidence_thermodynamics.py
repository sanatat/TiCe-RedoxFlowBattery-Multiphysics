"""Generate traceable mixed-acid solubility evidence, Ce-SO4 benchmark, and audit.
Run from any cwd: python run_evidence_thermodynamics.py
"""
from __future__ import annotations
from pathlib import Path
import csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from extensions.evidence_thermodynamics import (
    PATENT_SOLUBILITY, CE_SULFATE_LOG_BETA_2M_NACLO4_25C,
    screen_ce_evidence, conditional_ce_sulfate_fractions, equilibrium_data_gap_audit,
)
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'thermo_evidence'/'outputs'
FIG=ROOT/'thermo_evidence'/'figures'

def dump_csv(filename,rows):
    with (OUT/filename).open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def run():
    OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    evidence_rows=[]
    for p in PATENT_SOLUBILITY:
        evidence_rows.append(dict(msa_M=p.msa_M,h2so4_M=p.h2so4_M,
                                  ce3_limit_M=p.ce3_limit_M,
                                  ce4_limit_M='' if p.ce4_limit_M is None else p.ce4_limit_M,
                                  ce4_qualitative=p.ce4_qualitative or '',
                                  source=p.source,temperature_C='not_specified',
                                  level='patent-reported_obs_not_equilibrium'))
    dump_csv('01_observed_ce_solubility_patent_table.csv',evidence_rows)
    cases=[]
    for p in PATENT_SOLUBILITY:
        for soc in np.linspace(0,1,21):
            r=screen_ce_evidence(p.msa_M,p.h2so4_M,0.5,float(soc),25.0)
            cases.append(dict(msa_M=p.msa_M,h2so4_M=p.h2so4_M,ce_total_M=0.5,
                              CeIV_fraction=float(soc),ce3_M=r['ce3_M'],ce4_M=r['ce4_M'],
                              CeIII_flag=r['ce3']['status'],CeIV_flag=r['ce4']['status'],
                              evidence_basis='patent oxid-state specific; temperature unverified'))
    dump_csv('02_ce_oxidation_fraction_evidence_flags.csv',cases)
    # A categorical state map avoids presenting 'precipitation pressure' or
    # falsely precise saturation index from high ionic strength solutions.
    matrix=np.full((len(PATENT_SOLUBILITY),21),np.nan)
    for i,p in enumerate(PATENT_SOLUBILITY):
        for k,soc in enumerate(np.linspace(0,1,21)):
            r=screen_ce_evidence(p.msa_M,p.h2so4_M,0.5,float(soc))
            if r['ce4']['reported_limit_M'] is None:
                matrix[i,k]=np.nan
            else:
                matrix[i,k]=int(any(r[s]['status']=='above_reported_limit' for s in ('ce3','ce4')))
    fig,ax=plt.subplots(figsize=(9,4.5)); im=ax.imshow(matrix,cmap='RdYlGn_r',vmin=0,vmax=1,aspect='auto',interpolation='nearest')
    ax.set(xticks=[0,5,10,15,20],xticklabels=['0','.25','.50','.75','1.0'],
           yticks=list(range(len(PATENT_SOLUBILITY))),
           yticklabels=[f'{p.msa_M:g} M MSA + {p.h2so4_M:g} M H2SO4' for p in PATENT_SOLUBILITY],
           xlabel='Fraction of dissolved Ce in +IV state (assumed)',
           title='At 0.50 M total Ce: exceeding reported Ce(III/IV) limits?\nPatent evidence only; temperature unspecified; not a precipitate prediction')
    fig.text(.5,.015,'Green: below reported limits  |  Red: above at least one  |  White: Ce(IV) limit not numerical',
             ha='center',fontsize=9)
    fig.subplots_adjust(left=.28,right=.98,top=.80,bottom=.18)
    fig.savefig(FIG/'01_Ce_solubility_SOC_evidence_matrix.png',dpi=170);plt.close(fig)
    # Fixed-ionic-medium sulfate conditional beta only; deliberately separate.
    sulfate=np.logspace(-7,-1,160)
    spec=[conditional_ce_sulfate_fractions(float(s)) for s in sulfate]
    dump_csv('03_ce_sulfate_conditional_2M_NaClO4_only.csv',spec)
    fig,ax=plt.subplots(figsize=(9,5))
    for key,short in [('CeIV_free_fraction','Ce(IV) unbound'),('CeIV_SO4_fraction','CeSO4(2+)'),
                      ('CeIV_2SO4_fraction','Ce(SO4)2'),('CeIV_3SO4_fraction','Ce(SO4)3(2-)'),
                      ('CeIII_SO4_fraction','Ce(III) sulfate bound')]:
        ax.semilogx(sulfate,[s[key] for s in spec],label=short,lw=2)
    ax.set(xlabel='Prescribed free SO4(2-) concentration (M)',ylabel='Conditional Ce species fraction',
           ylim=(-.02,1.04),title='Published Ce sulfate log beta, 25 C, 2 M NaClO4 ONLY\nNOT directly usable for MSA/H2SO4 battery electrolyte')
    ax.grid(True,alpha=.25);ax.legend(fontsize=9)
    fig.tight_layout();fig.savefig(FIG/'02_Ce_sulfate_conditional_2M_NaClO4.png',dpi=170);plt.close(fig)
    audit=equilibrium_data_gap_audit()
    summary={'status':'evidence-based Ce screening, NOT verified Ti/Ce mixed-acid precipitation thermodynamics',
             'source_notes':{
               'patent':'https://patents.google.com/patent/US12119527B2/en',
               'published_Ce_mixed_acid':'https://doi.org/10.1021/ef200354b',
               'Ti_speciation':'https://doi.org/10.1039/D1RA04284J',
               'Ti_thermohydrolysis':'https://doi.org/10.1016/j.jct.2018.11.001',
               'Ce_background':'Rare Earth Elements Equilibria in Aqueous Media, Table 11.2 (2 M NaClO4 at 25C)'},
             'ce_sulfate_conditional_log_beta':CE_SULFATE_LOG_BETA_2M_NACLO4_25C,
             'worked_example':screen_ce_evidence(2,0.5,0.5,0.5),
             'audit':audit}
    (OUT/'audit.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Wrote 3 evidence CSV files, 2 plots, and audit.json in',OUT.parent)
    return summary

if __name__=='__main__':run()
