## Title

Graded sensorimotor retraining compared to sham control for people with chronic low back pain: The causal effect of following the protocol 

Joanna Diong<sup>1,2</sup>, Joy Shi<sup>3,4</sup>, Matthew K Bagg<sup>5,6,7</sup>, Benedict M Wand<sup>6</sup>, Harrison J Hansford<sup>5,7</sup>, Aidan G Cashin<sup>5,7</sup>, James H McAuley<sup>5,7</sup>

1. Sydney School of Public Health, Faculty of Medicine and Health, The University of Sydney, Sydney, Australia
2. Charles Perkins Centre, The University of Sydney, Sydney, Australia
3. Mongan Institute Health Policy Research Centre, Massachusetts General Hospital, Boston, USA
4. CAUSALab, Department of Epidemiology, Harvard T.H. Chan School of Public Health, Harvard University, Boston, USA
5. Neuroscience Research Australia (NeuRA), Sydney, Australia
6. School of Health Sciences, Faculty of Medicine, Nursing and Midwifery and Health Sciences, The University of Notre Dame Australia, Fremantle, Australia
7. School of Health Sciences, Faculty of Medicine and Health, University of New South Wales, Sydney, Australia

## Suggested citation

Diong J, Shi J, Bagg M, Wand BM, Hansford HJ, Cashin AG, McAuley JH (2026) Graded sensorimotor retraining compared to sham control for people with chronic low back pain: The causal effect of following the protocol. 

## Corresponding author

__Manuscript__\
Dr Joanna Diong\
Sydney School of Public Health\
Faculty of Medicine and Health\
The University of Sydney\
joanna.diong@sydney.edu.au\
+61 2 8627 5491

__Code__\
Dr Joy Shi\
Mongan Institute Health Policy Research Centre\
Massachusetts General Hospital\
joyshi@hsph.harvard.edu

## Data

De-identified data are available on reasonable request.

## R code 

R code is made open for transparency, and researchers who may wish to read or apply similar procedures. 

The R script will call the raw data file stored in a folder labelled _data_ in the current directory, to perform statistical analysis. 

R code to perform causal per-protocol analysis was written by Joy Shi (R v4.1): 

* _RESOLVE_Bootstrap_Function_20250815.R_: R script to import dataset and perform statistical analyses

To run the R script, create a folder labelled _data_ in the current directory. Place the dataset within that folder.

In the script file line 38, set the file path for project to the current directory.

Note, paths are currently set for Linux or Mac operating systems using forward slashes `/`. 
They will need to be updated for Windows using backward slashes `\`, 
or for different project locations on different machines. 

Run the R script in e.g. R Studio to replicate the analysis.

## Latex code

LaTeX code to generate Fig 1 was written by Joanna Diong:

* _fig1/fig1.tex_, which requires the dependency _myStyle.sty_

Run _fig1.tex_ in e.g. Texmaker to generate a PDF of the figure and legend.


