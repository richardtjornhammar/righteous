# Righteous Pathway Factor Analysis
Decomposes a set of expressions into a group expression using
factor analysis. The expression regulation can be studied via 
an ANOVA that relates it to the observables in the journal 
file. The final p values are then FDR corrected and the resulting
q values are produced.

The journal and analyte expression file must be ordered
the same way with respect to the samples that are
positioned on the columns.

Visit the active code via :
https://github.com/richardtjornhammar/righteous
https://pypi.org/project/righteous-fa/
https://zenodo.org/record/3373405

Download the data used in the examples from:
https://zenodo.org/record/3407557

That can be recreated by using the scripts and instructions in:
https://zenodo.org/record/3406460

Install with :
pip install righteous-fa

Examples can be run by issuing:

>>> from righteous.examples import *
>>> run_MS()

