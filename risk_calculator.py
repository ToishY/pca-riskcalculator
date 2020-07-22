# -*- coding: utf-8 -*-
"""
Created on Tue Jun 2 20:33:55 2020
Modified on Wed Jul 22 15:49:28 2020

@author: ToishY

"Risk calculator" / "Prostaatwijzer" Type 3/4 TRUS/DRE(/Phi)

This script was rewritten in Python with the purpose for automated testing for 
large amount of data of (possible) prostate cancer patients.

The risk calculator is based on logistic regression and returns two values:
Chance of detecting cancer & the chance of detecting clinical significant cancer

"""

import math
import numpy as np
    
def read_data(file, deli=',',header=1, y_col=-1):
    """ Read data and return without possible header """
    
    raw = np.genfromtxt(file, delimiter=deli)
    return raw[header:,:y_col]

### Prostaatwijzer defs ###

def volume_to_volume_classes(volume):
    """ Categorize volume into 3 classes """
    
    if (volume < 30.0):
        volume_class = 25 
    elif (volume >= 30.0 and volume < 50):
        volume_class = 40;
    else:
        volume_class = 60;
  
    return volume_class                   

def pcrc_formula_prev_biopsy_pos(psa, dre_outcome, volume, age, pirads):
    """ Get chance prediction if previous biopsy was positive previously """
    
    # Volume class
    volume_class = volume
    wpirads = 0.0
    
    """
    The JS webtool offsets PIRADS score with -1; So in the tool, an
    actual value of PIRADS 1, is passed as 0 in the argument; 5 is passed as 4.
    This was changed back to remove confusion, so argument pi-rads 1, is 1.
    Note: the results for PIRADS score 1 and 2 do not differ due to the pirads
    weight being 0.0. 
    """
    
    if(int(pirads) == 3):
        wpirads = -0.38169
    elif(int(pirads) == 4):
        wpirads = 0.61290
    elif(int(pirads) == 5):
        wpirads = 2.01363
    
    # Based on logistic regression
    psaLog2 = math.log(psa) / math.log(2)
    volumeLog2 = math.log(volume_class) / math.log(2)

    # Weights/bias
    aux =  ((-1.826 + 1.024 * (psaLog2 - 2.0)  - 1.50 * (volumeLog2 - 5.4) + 0.992 * dre_outcome) * 1.139 ) + 1.325

    val_prior = -3.52647 + (-1.165) + aux * 0.78810 + age * 0.04832 + wpirads;
    
    """
    The PHI score here is not the same as the PHI scored described in other 
    literature, because the PHI = (p2PSA/fPSA) x (square root of PSA) [1]
    The "prostaatwijzer" tool does not take p2PSA or free PSA as arguments,
    so I don't think this `PposRC45DRE_PHI` is really meant to be PHI.
    Seems like plain logistic regression: (1/(1+e^-(b0+b1*x1)))
    
    [1] DOI: 10.1590/S1677-5538.IBJU.2016.0256
    """
    
    PposRC45DRE_PHI = 1/(1+math.exp(-(val_prior)))

    return PposRC45DRE_PHI

def pcrc_formula_prev_biopsy_pos_serious(psa, dre_outcome, volume, age, pirads):
    """ Get chance prediction if previous biopsy was positive previously serious """
    
    volume_class = volume
    wpirads = 0.0
    
    if(int(pirads) == 3):
        wpirads = 0.44510
    elif(int(pirads) == 4):
        wpirads = 1.47778
    elif(int(pirads) == 5):
        wpirads = 2.61551
        
    calc_age = age * 0.04055
    calc_psaLog2 = ( ( math.log(psa) / math.log(2) ) - 2.0) * 1.177
    calc_dre = 1.813 * dre_outcome
    calc_volume = (math.log(volume_class) / math.log(2) - 5.4) * -1.526
    intercept = -3.457
		
    aux = -4.05287 + (-1.61) + ( ( (intercept + calc_psaLog2 + calc_dre + calc_volume) * 0.718 + 1.069 ) * 0.69846) +  calc_age + wpirads
    PposRC45DRE_PHI_ser = math.exp(aux)/(1+math.exp(aux))
    
    return PposRC45DRE_PHI_ser

def pcrc_formula_prev_biopsy_neg(psa, dre_outcome, volume, age, pirads):
    """ Get chance prediction if previous biopsy was negative previously """
    
    priorbiopsy = 1.0
    # Volume class
    volume_class = volume
    wpirads = 0.0
    
    if(int(pirads) == 3):
        wpirads = 0.08258
    elif(int(pirads) == 4):
        wpirads = 1.40371
    elif(int(pirads) == 5):
        wpirads = 1.91967
    
    psaLog2 = math.log(psa) / math.log(2)
    volumeLog2 = math.log(volume_class) / math.log(2)

    val_prior = -2.62636 + (-0.745) + ((((-1.470 - 0.677 * priorbiopsy + (0.576 - 0.423 * priorbiopsy) * (psaLog2 -2)-1.043 * (volumeLog2 - 5.5) + 0.68 * dre_outcome))*1.4857+2.3426)) * 0.78749+age*0.02467+wpirads;

    PposRC45DRE_PHI = 1/(1+math.exp(-(val_prior)))
	
    return PposRC45DRE_PHI

def pcrc_formula_prev_biopsy_neg_serious(psa, dre_outcome, volume, age, pirads):
    """ Get chance prediction if previous biopsy was negative previously serious"""
   
    priorbiopsy = 1.0
    # Volume class
    volume_class = volume
    wpirads = 0.0
    
    if(int(pirads) == 3):
        wpirads = 0.94909
    elif(int(pirads) == 4):
        wpirads = 2.67536
    elif(int(pirads) == 5):
        wpirads = 3.05243
        
    psaLog2 = math.log(psa) / math.log(2)
    volumeLog2 = math.log(volume_class) / math.log(2)

    val_prior_ser = -5.19208+ (-0.745) +((((-3.489 - 1.136 *  priorbiopsy + (1.075 - 0.434 * priorbiopsy) * (psaLog2 -2) - 1.501 * (volumeLog2 - 5.5) + 1.311 * dre_outcome)) * 0.9420+ 2.1452))* 0.74027+ age*0.04286+wpirads

    PposRC45DRE_PHI_ser= 1/(1+math.exp(-(val_prior_ser)))

    return PposRC45DRE_PHI_ser

# %% Main
if __name__ == "__main__":
    """ Main """
    
    # Read in CSV `read_data`
    data = read_data("sample_data.csv")
    
    # Init normal/serious cancer risk; sample-by-2
    risk = np.empty([data.shape[0],2])
    
    # Loop over data
    for idx, item in enumerate(data):
        # Assuming previous biopsy was negative (unknown beforehand) for all entries
        # Assuming the DRE outcome was normal (unknown beforehand) for all entries
        detectable_cancer_risk = round( pcrc_formula_prev_biopsy_neg(item[1],0,volume_to_volume_classes(item[2]),item[0],item[3]), 3)
        significant_cancer_risk = round( pcrc_formula_prev_biopsy_neg_serious(item[1],0,volume_to_volume_classes(item[2]),item[0],item[3]), 3)
        
        # Append
        risk[idx] = np.array([detectable_cancer_risk,significant_cancer_risk])