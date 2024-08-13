# Gender and anxiety reveal distinct computational sources of underconfidence

 
  #################### #################### #################### 
####################

All code and data for reproducing analyses and figures for:

Sucharit Katyal, and Stephen M. Fleming. 2024. “Gender and anxiety reveal distinct computational sources of underconfidence.” PsyArXiv. May 24. 
https://psyarxiv.com/qcg92/


 
  #################### #################### #################### 
####################

Analyses and figure generation was performed in R v4.4.1 in RStudio v2024.04.2 
(see: _code/r_). For reproducing analyses and figures, step through the following code:

    analysis and figure 2:  code/confTime_analysis.R
    
    figure 3A–3F:           code/confTime_simulateModel.R

    figure 3G and 3I:       code/confTime_fitBehaviour.R
    
    figure 3H–3J:           code/plot_metaDDM.R

NOTE: Before performing analyses, please unzip the zipped .csv files _data/Katyal_et_al_2023/mbsExp2_noperfexcl.csv.zip_

  #################### ####################

The computational model is available in the _code/meta-ddm-model_ directory

To run the model call the files _code/meta-ddm-model/fit_metaDDM_Exp*.R_ separately for each of Exps 1–4. 

For Exps 1 and 2, the model was run separately for the perception and memory tasks, which can be chosen by setting the _task2fit_ parameter equal to 1 and 2 respectively.

Use _code/meta-ddm-model/fit_metaDDM_recover.R_ to perform model simulation and parameter recovery.

Model parameters are plotted using the _code/plot_metaDDM.R_ file.

Model fitting functions for Exp 1 & 2 fit continuous confidence ratings while Exp 3 & 4 fit discrete confidence ratings.

  #################### #################### #################### 
####################

To run the models on your own data, modify one of the fit_metaDDM_Exp* functions by pointing it to your dataset containing trial-by-trial data on subject id, choice, accuracy, choice response times, confidence, and confidence response times.

For help with running the model on your data, contact: ska [@] psy [dot] ku [dot] dk


  #################### #################### #################### 
  
  License

This code is being released with a permissive open-source license. Please feel free to use or adapt the code as long as you follow the terms of the license enumerated below. If you use the model in a publication, we ask that you cite the following paper:



Copyright (c) 2024, Sucharit Katyal

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
