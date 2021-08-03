# ChemTree

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs)

Calculate similarity between chemical compounds using Molecular Assembly theory.

## Usage / Examples

#### Generate molecular assembly fragments

```python
python generate_mas.py
```

Generates molecular assembly fragments for a variety of compounds. Currently supports KEGG compounds, including BRITE classifications, and opioids.

Assumes KEGG compounds and BRITE classifications are in `KeggData\`, and drug compounds are in `DrugData\`

#### Train molecular assembly fragments using Word-to-Vector

In progress

#### Find similarity between compounds

```python
python find_similarity.py
```

Finds similarity between compounds. Currently finds the Tanimoto similarity between calculated RDKit molecular fingerprints. Visualizations can 
be found in the `tanimoto_tree.ipynb` jupyter notebook.
