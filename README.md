# Risk Calculator Type 3/4 TRUS/DRE(/Phi)

Risk calculator Type 3/4 from the JavaScript WebGUI rewritten to Python 3.7.6.

## Description

The original JavaScript code was found on [Prostatecancer-Riskcalculator.com](http://www.prostatecancer-riskcalculator.com/seven-prostate-cancer-risk-calculators) and is based on logistic regression ([more information on testing/validation of the original tool can be found here](http://swopresearch.nl/wp-content/uploads/2018/04/Tekst_proefschrift_dr._Van_Vugt_-_Impact_studie.pdf
)).
This tool gives an indication for finding (clinical significant) prostate cancer if a biopsy would be performed, based on several features.

By rewriting the JavaScript code to Python, processing of large amounts of sample data is 
automized and therefore less time-consuming.


## Data


Supplied artificial sample data (CSV) for demonstration purposes. This data contains the following features and target:

* Age
* PSA (µg/L)
* Prostate volume (mL)
* PIRADS score (1-5)
* PCa (0/1)

## Settings

The tool includes two extra features which are not included in our sample data:

* DRE_outcome (abnormal/normal) defined as `1` or `0` respectively.
* Previous biopsy (positive/negative) defined as functions with trailing `_pos` or `_neg` names respectively.

It is up to the user to decide what options to choose in case of missing data, or change the function arguments manually if these features are known for each sample.

## Results

The tool returns two types of probabilities of finding (clinical significant) prostate cancer):

* Risk of finding prostate cancer
* Risk of finding high-grade/significant prostate cancer

These results are appended in a `sample-by-2` array named `risk`, with probabilities between `0` and `1` (rounded to 3 decimal places).

## Note

Risk calculators can give an indication for prostate cancer by determining probabilities for each sample, but these tools are not 100% accurate and should not be solely used for diagnosing prostate cancer.
