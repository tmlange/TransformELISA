# TransformELISA
Transformation of raw optical density (OD) values from an enzyme-linked immunosorbent assay (ELISA) to estimate the corresponding logarithmised protein concentration in the samples. 

TransformELISA.R:
R script to perform the transformation. The script expects a column sample that contains the string bc for each buffer control, the string sample for each sample in the trial and the string serialdilution and the corresponding number of the sample of the serial dilution for each sample of the serial dilution.
Furthermore, the user needs to set the technical limit of the ELISA machine, thus, the maximum OD value that can be measured. This will be used as the top asymptote of the non-linear regression model. 

TrialData.csv:
Example of a data set which fulfills all the requirements of the R script TransformELISA.R. 
The data set was collected in a trial with eleven ELISA plates and a serial dillution containing twelve samples.

# Requirements
R 4.1.2 or higher and the R package minpack.lm

# Licenses
The software is licensed under the MIT license. The data set is licensed under the Creative Commons Legal Code (CC0 1.0 Universal). For further details see the LICENSES file. 

# Citation
If you use this software, please cite it as below:

DOI: 
https://doi.org/10.1186/s12985-022-01804-3

APA:
Lange, T. M., Rotärmel, M., Müller, D., Mahone, G. S., Kopisch-Obuch, F., Keunecke, H., & Schmitt, A. O. (2022). Non-linear transformation of enzyme-linked immunosorbent assay (ELISA) measurements allows usage of linear models for data analysis. Virology Journal, 19(1), 1-11. 

BibTeX:
@article{,
  title={Non-linear transformation of enzyme-linked immunosorbent assay ({ELISA}) measurements allows usage of linear models for data analysis},
  author={Lange, Thomas M. and Rot{\"a}rmel, Maria and M{\"u}ller, Dominik and Mahone, Gregory S. and Kopisch-Obuch, Friedrich and Keunecke, Harald and Schmitt, Armin O.},
  journal={Virology Journal},
  volume={19},
  number={1},
  pages={1--11},
  year={2022},
  publisher={BioMed Central}
}
