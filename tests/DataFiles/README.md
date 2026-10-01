# Example EBSD data

| File | Format | Phases (crystal system) | Size | ebsdlab |
|---|---|---|---|---|
| `EBSD.ang` | EDAX/TSL `.ang`, hex grid | Copper (cubic, FCC) | 2.2 MB | loads |
| `EBSD.osc` | EDAX/TSL binary `.osc` | Copper (cubic, FCC), same scan as `EBSD.ang` | 1.0 MB | loads; pass `symmetry='cubic'` |
| `EBSD.oim` | EDAX OIM project | refers to `EBSD.ang` | 1.9 MB | not read |
| `Catillopecten.crc` + `.cpr` | Oxford binary `.crc` | Calcite (trigonal), Aragonite (orthorhombic) | 0.4 MB | loads |
| `Eclogite.crc` + `.cpr` | Oxford binary `.crc`, 6208 points | Garnet (cubic), Omphacite, Hornblende, Clinozoisite (monoclinic), Albite (triclinic), Quartz (trigonal) | 0.2 MB | loads |
| `Ti_ZrN.ctf` | Oxford text `.ctf`, 100 x 100 | Ti alpha (hexagonal), Ti beta (cubic, BCC), ZrN (cubic, FCC) | 0.6 MB | loads |
| `W_TKD.ctf` | Oxford text `.ctf`, TKD, 600 x 399 | Tungsten (cubic, BCC), 58 % not indexed | 13 MB | loads |


## Sources

### `EBSD.ang`, `EBSD.osc`, `EBSD.oim`
- Source: measured by Steffen Brinckmann
- License: MIT, as ebsdlab

### `W_TKD.ctf`
- Original name: `2000nm_005s_2TKD.ctf`
- Source: transmission Kikuchi diffraction map of a notched tungsten crystal
- Authors: Jin Wang, IMD-1, FZ Jülich
- License: MIT, as ebsdlab

### `Catillopecten.crc`, `Catillopecten.cpr`
- Original names: `A checa Catillopecten Site 24 Map Data 92.crc` / `.cpr`
- Source: Zenodo, "Raw crystallographic (EBSD) data for the article: Crystallographic control of the fabrication of an
  extremely sophisticated shell surface microornament: tuning into the aerials of the glass scallop Catillopecten"
- URL: https://zenodo.org/records/6375610, DOI 10.5281/zenodo.6375610
- Authors: A. G. Checa, C. Salas, F. M. Varela-Feria, A. B. Rodríguez-Navarro, C. Grenier, G. M. Kamenev, E. M. Harper
- License: CC-BY-4.0

### `Eclogite.crc`, `Eclogite.cpr`
- Original names: `stagescanS63.crc` / `.cpr`
- Source: Zenodo, "Punta Telcio Zermatt-Saas Eclogite Mineral Data - Chemistry and Crystallography"
- URL: https://zenodo.org/records/7837199, DOI 10.5281/zenodo.7837199
- Authors: D. D. McNamara, J. Wheeler, M. Pearce, D. Prior
- License: CC-BY-4.0

### `Ti_ZrN.ctf`
- Original name: `FIgure_10_EBSD.ctf`
- Source: Zenodo, "Dataset for paper entitled, 'The potential for grain refinement of Wire-Arc Additive Manufactured
  (WAAM) Ti-6Al-4V by ZrN and TiN inoculation'"
- URL: https://zenodo.org/records/5708619, DOI 10.5281/zenodo.5708619
- Authors: J. Kennedy, A. Davis, A. Caballero
- License: CC-BY-4.0

## MTEX examples, not included
MTEX ships further examples in `data/EBSD/`:
https://github.com/mtex-toolbox/mtex/tree/c836b404a6729ef339857e216ff4adda143d38fb/data/EBSD

They are not copied here because MTEX is licensed GPL-2.0, while ebsdlab and the files above are MIT or CC-BY-4.0.
`tests/test_mtex.py` downloads some of them into `tests/mtex_cache/` (not committed) and compares plots with
`tests/baseline/test_mtex_*.png`

| File | Content | Test |
|---|---|---|
| `DC06_2uniax.ang` | Iron alpha (cubic, BCC), square grid | IPF baseline |
| `olivineopticalmap.ang` | Olivine, Enstatite (orthorhombic), Dolomite (trigonal), Chalcopyrite (tetragonal) | phase-map baseline |
| `copper.osc` | Copper (cubic, FCC) | IPF baseline |
| `EDXLMDTi64.crc` + `.cpr` | Ti alpha (hexagonal), Ti beta (cubic, BCC) | IPF baseline |
| `twins.ctf` | Magnesium (hexagonal), twinned | IPF baseline |
