# Ti–Ce mixed-acid precipitation: evidence-based extension (not an equilibrium Ksp model)

## What can and cannot be claimed

A **genuine thermodynamic saturation index** requires (i) a precise solid phase and its equilibrium constant, (ii) full aqueous complexation/hydrolysis constants for Ti(III/IV), Ce(III/IV), MSA and sulfate at the actual high ionic strengths, and (iii) an activity model fitted to this specific mixed-acid electrolyte (typically SIT/Pitzer at these ionic strengths). **The article + SI do not provide these inputs**. The repository therefore does NOT report a numerical Ti/Ce saturation index, mineral precipitation mass, temperature-dependent Ksp, or an experimentally optimized acid recipe.

**What is now implemented:**

1. **Published, oxidation-state-resolved Ce(III)/Ce(IV) solubility limits** at five *exact nominal acid compositions* from the measured-observation table in US Patent 12119527 B2 (Table 1). The oxidation-state fraction is a supplied input; the software compares Ce(III) and Ce(IV) concentrations separately against reported limits. It does NOT interpolate to other acid compositions. **The patent table does not establish measurement temperature or replicate error bars**, and is not a substitute for peer-reviewed equilibrium measurements.
2. **Conditional Ce/sulfate complex speciation from published formation constants** at **25 °C in 2 M NaClO4 medium only** (Table 11.2 of *Rare Earth Elements Equilibria in Aqueous Media*). This supports mass-action and stoichiometric constraints, **not precipitation**. The model forbids transferring these coefficients to MSA/H2SO4.
3. An explicit **data-gap audit** and literature links for an eventual defensible high-ionic-strength chemical-equilibrium model.

### Why the oxidized state matters

For example, the patent-reported acid mix 2 M MSA + 0.5 M H2SO4 has Ce(III) reported solubility 0.8 M, but Ce(IV) only 0.2 M. If we assume 0.5 M total dissolved Ce at 50% Ce(IV), the dissolved Ce(IV) component is 0.25 M, **above the reported Ce(IV) solubility**. This is a literature evidence flag, *not proof of a particular precipitate forming*.

Importantly, **another published study** (Xie et al. 2011, *Energy & Fuels*) reports that a 1 M Ce solution in 2 M MSA + 0.5 M H2SO4 is sufficiently stable for >1 month up to 313 K. Its abstract does not explicitly establish identical Ce(IV) fraction, mixing protocol, or precipitation-equilibration conditions. **These results must not be treated as equivalent or force-fitted to one Ksp**. Literature on 0.6 M Ce mixed-acid solutions also reports precipitation risks at elevated temperature and high Ce(IV) fraction. These differences make oxidation state, experimental history and time necessary validation dimensions.

### Titanium remains a genuine unresolved chemistry gap

Ti(IV) hydrolyzes strongly and forms sulfate and oxo/polymeric species. Published Ti(IV) hydrolysis and formation constants differ substantially among authors, even in simpler aqueous acids. The available redox-flow-battery studies demonstrate that TiOSO4 can be soluble in strongly acidic sulfate media, but this does *not* yield a verified Ti(solid)/Ti(aq) saturation calculation for a co-dissolved Ti/Ce, MSA/H2SO4 mixture.

**Do not construct an imaginary PHREEQC database** from DFT binding energies in the Fu et al. SI. Those gas-/cluster-level electronic energies are not the required aqueous standard Gibbs energies / activities and cannot be re-labelled log Ksp.

## Reproduce

```bash
python run_evidence_thermodynamics.py
python -m pytest -q tests/test_evidence_thermodynamics.py
```

Outputs in `thermo_evidence/`:

- `outputs/01_observed_ce_solubility_patent_table.csv`: directly transcribed published numbers with provenance;
- `outputs/02_ce_oxidation_fraction_evidence_flags.csv`: derived, explicitly conditional classification at 0.5 M total Ce;
- `outputs/03_ce_sulfate_conditional_2M_NaClO4_only.csv`: conditional complex fractions in the **published ionic medium only**;
- `figures/01_Ce_solubility_SOC_evidence_matrix.png`: evidence flags; **not measured precipitation map**;
- `figures/02_Ce_sulfate_conditional_2M_NaClO4.png`: valid-medium conditional mass action; **not MSA speciation**;
- `outputs/audit.json`: source URLs, an example and explicit missing parameters.

## References and evidence grade

| Evidence | Source | Grade / applicability |
|---|---|---|
| Ce(III), Ce(IV) mixed-acid solubility Table 1 | [US12119527B2](https://patents.google.com/patent/US12119527B2/en) | Patent-reported solubility values, **not** a verified equilibrium Ksp dataset; temperature unspecified |
| Mixed acid 1 M Ce stability >1 month, up to 313 K | [Xie et al., *Energy & Fuels* (2011), 10.1021/ef200354b](https://doi.org/10.1021/ef200354b) | Peer-reviewed experimental observation but oxidation-state specification from abstract insufficient |
| Cerium mixed-acid stability sensitive to temperature and SOC | [Effect of Mixed Acid Media on the Positive Side of Hybrid Zn–Ce RFB (2014)](https://www.sciencedirect.com/science/article/abs/pii/S0013468614013589) | Peer-reviewed; experimental conditions not identical to Ti–Ce |
| Ce sulfate conditional beta in 2 M NaClO4, 25 C | *Rare Earth Elements Equilibria in Aqueous Media*, Table 11.2 (2016 book chapter) | **Chemically meaningful only within this background medium**; no transfer to MSA/H2SO4 |
| Ti(IV) hydrolysis and sulfate complexation uncertainties | [RSC Advances (2021), 10.1039/D1RA04284J](https://doi.org/10.1039/D1RA04284J) | Peer-reviewed discussion of incompatible published constants |
| Temperature-dependent Ti hydrolysis | [J. Chem. Thermodyn. (2019), 10.1016/j.jct.2018.11.001](https://doi.org/10.1016/j.jct.2018.11.001) | Peer-reviewed study in sulfate/chloride; cannot directly transfer to mixed MSA |
| Fu et al. battery electrolyte, original SI | [ChemSusChem (2026), 10.1002/cssc.202502732](https://doi.org/10.1002/cssc.202502732) | Exact 2 M MSA battery experimental parameters, not high-acid equilibrium database |

## Path to full validation

To make a **thermodynamically defensible quantitative precipitation model** rather than a screening tool, secure oxidation-state-resolved measurements at several acid compositions and temperatures; measure acid/sulfate activities, soluble complexes, and solid identities (XRD/Raman/FTIR) where precipitation actually occurs; choose an ionic-activity model and independently fit & validate all equilibrium constants; then incorporate the validated reactions and phases into PHREEQC or an equivalent Gibbs minimizer. An experiment-derived solubility surface also needs separate uncertainty bounds and long-term kinetic/stability checks.