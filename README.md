# Automatic Phantom Assembling
Code to automatically assemble phantoms from DICOM CT scans and RT STRUCT files using machine learning models.

If you find this code helpful, please consider citing our work:
 
```
@article{virgolin2020machine,
  title={Machine learning for the prediction of pseudorealistic pediatric abdominal phantoms for radiation dose reconstruction},
  author={Virgolin, Marco and Wang, Ziyuan and Alderliesten, Tanja and Bosman, Peter AN},
  journal={Journal of Medical Imaging},
  volume={7},
  number={4},
  pages={046501},
  year={2020},
  publisher={SPIE}
}
```

## Details
This code is an anonymized version of the one actually used to generate phantoms for [1]. 
In the lookup table, information on predictions of machine learning algorithms is pre-collected to assemble the phantoms. 
For a practical use, lookup tables must be replaced by actual predictions of pre-trained models, so models need to be integrated and called-upon when needed, to generate predictions on the fly given the features of the input patient.
