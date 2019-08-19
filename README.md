# APhA
Code to automatically assemble phantoms from DICOM CT scans and RT STRUCT files.

If you find this code helpful, please consider citing our work:
[1] *Machine learning for automatic construction of pseudo-realisticpediatric abdominal phantoms*. Marco Virgolin, Ziyuan Wang, Tanja Alderliesten and Peter A.N. Bosman.

## Details
This code is an anonymized version of the one actually used to generate phantoms for [1]. 
In the lookup table, information on predictions of machine learning algorithms is collected to assemble the phantoms. 
For a practical use, lookup tables must be replaced by actual predictions of pre-trained models, that can be generated on the fly from patient features.
