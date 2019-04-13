# EPSDC

Ensemble Prediction of Synergistic Drug Combinations
The code of EPSDC was implemented in R with Rstudio. The supported vector machine (SVM) and Navie Bayes (NB) are implemented by functions from library "e1071". Meanwhile, We employed function trapz from library "caTools" to calculate the area under the ROC curve. All experiments are carried out with Dell Precision T3600 (3.2 GHz Intel Xeon, 16G memory).
The main objective of the task is to accuratly and efficiently screen synergistic drug combinaitons by integrating multi-sources knoledge. The related data in the experiment were provided inliterature "Ensemble Prediction of Synergistic Drug Combinations Incorporating Biological, Chemical, Pharmacological and Network Knowledge".
The framework EPSDC has some highlight as follows:
1.	The biological, chemical, pharmacological and network data could provide more abundant knowledge that enables EPSDC to predict potential drug combinations with multi perspectives;
2.	Network-based base predictor using transductive learning could prioritize potential drug combinations without negative drug combination sample, which is hard to be obtained. Meanwhile, integrating network-based predictor could correct the error caused by negative sample selection in feature-based predictor;
3.	Our method outperform the state-of-the-art method under five-fold cross validation scheme.


Note: supplementary material2.xlsx and supplementary material3.xlsx are the corrected dataset.
