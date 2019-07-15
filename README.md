# rotavirusMalawi

This website contains the MATLAB code and aggregate data on rotavirus gastroenteritis (RVGE) cases from Queen Elizabeth Central Hospital in Blantyre, Malawi needed to reproduce the results presented in Pitzer VE et al, Science Transl Med, 2019.

DATA

The file 'malawidata.mat' contains the following variables:
blantyre_pop - Monthly population estimates for Blantyre (columns: (1) year, (2) month, (3) Blantyre rural pop, (4) Blantyre city pop, (5) Blantyre district pop (total), (6) proportion of pop in Blantyre city
blantyre_pop2 - Estimates of Blantyre city population by year (1998-2014, rows) and 5-year age group (0-4y to 80+y, columns)
datepop - Date of Malawi population census data (columns: year, month, day, hour, min, sec)
malawiagedist - Proportion of the Malawi population in 5-year age groups (0-4y to 80+y, columns)
malawi_cbr - Crude birth rate in Malawi (live births per 1,000 total population) for years in 'datepop'
malawi_cdr - Crude death rate in Malawi (deaths per 1,000 total population) for years in 'datepop'
malpop - Population of Southern Region of Malawi by age group (<1y, 1-4y, 5-9y,...,75y+) according to 2008 census
malVS3date - Date (year, month, day) of post-vaccination surveillance data
malweek08 - Date (in MATLAB number format) of pre-vaccination surveillance data from 2008-2009
malweeknum - Date (in MATLAB number format) of pre-vaccination surveillance data from 1997-2007
negtest08_age - Number of rotavirus-negative AGE cases by week (rows) and age group (columns: 0-11m, 12-23m, 2-4y) for 2008-09
negtest08_movavg - 105-week moving average of the number of rotavirus-negative AGE cases by week for 2008-09
negtest_age - Number of rotavirus-negative AGE cases by week (rows) and age group (columns: 0-11m, 12-23m, 2-4y) for 1997-2007
negtest_movavg - 105-week moving average of the number of rotavirus-negative AGE cases by week for 1997-2007
negtestage08_movavg - 105-week moving average of the number of rotavirus-negative AGE cases by week and age group for 2008-09
negtestage_movavg - 105-week moving average of the number of rotavirus-negative AGE cases by week and age group for 1997-2007
negtestVS3 - Number of rotavirus-negative AGE cases by week (rows) and age group (columns: 0-11m, 12-23m, 2-4y) for 2012-17
negtestVS3_movavg - 105-week moving average of the number of rotavirus-negative AGE cases by week and age group for 2012-17
pfit_mal - maximum a posterior estimates of parameters (rows) for model fit to pre-vaccination data assuming different relative infeciousness of subsequent infections (0.1-0.5, columns)
repeffV3 - relative reporting effort by age group (columns: <1m,...,23m,1y,...,4y,5-9y,...,75+y) for the post-vaccination period, as determined by moving average of the rotavirus-negative AGE cases compared to pre-vaccination period
rotamal08T2 - Number of RVGE cases by week (rows) and age group (columns: <1m,...,23m,1y,...,4y,5-9y,...,75+y) for 2008-09
rotamalT2 - Number of RVGE cases by week (rows) and age group (columns: <1m,...,23m,1y,...,4y,5-9y,...,75+y) for 1997-2007
rotamalVS3 - Number of RVGE cases by week (rows) and age group (columns: <1m,...,23m,1y,...,4y,5-9y,...,75+y) for 2012-2018
rotamalVS3_unvacc - Number of RVGE cases among unvaccinated individuals by week and age group for 2012-2018
rotamalVS3_unvacc - Number of RVGE cases among vaccinated individuals by week and age group for 2012-2018
vcov_mavg - 27-week moving average of the 1- and 2-dose rotavirus vaccine coverage estimates for 2012-2018

MODEL FITTING 

MAP: Contains the code needed for maximum a posteriori estimation
malawimodelfit - Start here. M-file calls on other code to identify best-fit parameters (by minimizing log-likelihood).

rasisM1 - Function file containing differential equations for pre-vaccination model (called by 'rotafitM1.m', etc.)
rasisV1 - Function file containing differential equations for post-vaccination models 1-2 
rasisV1w - Function file containing differential equations for post-vaccination models 3-4
rotafitM1 - Function file that calculates the negative log-likelihood for model fit to pre-vaccination data
rotafitM1nb - Function file that calculates the negative LL for model fit to pre-vaccination data (assuming negative binomial likelihood distribution)
rotafitMV1 - Function file that calculates the negative log-likelihood for Model 1 fit to post-vaccination data
rotafitMV1nr - Function file that calculates the negative log-likelihood for Model 2 fit to post-vaccination data
rotafitMV1w - Function file that calculates the negative log-likelihood for Model 3 fit to post-vaccination data
rotafitMV1wnr - Function file that calculates the negative log-likelihood for Model 4 fit to post-vaccination data

MCMC: Contains the code needed for MCMC estimation of 95% credible intervals for model parameters

rasisM1 - Function file containing differential equations for pre-vaccination model 
rasisV1 - Function file containing differential equations for post-vaccination models 1-2 
rasisV1w - Function file containing differential equations for post-vaccination models 3-4
rotamal_LL - Function file that calculates log-likelihood of pre-vaccination model fit to data
rotamalawi_mcmc1 - M-file used for fitting model to pre-vaccination data (NOTE: takes up to 1 week to run)
rotamalawi_mcmcV_bc - M-file used for fitting Model 1 to post-vaccination data (NOTE: takes up to 1 week to run)
rotamalawi_mcmcV_nr - M-file used for fitting Model 2 to post-vaccination data (NOTE: takes up to 1 week to run)
rotamalawi_mcmcVw - M-file used for fitting Model 3 to post-vaccination data (NOTE: takes up to 1 week to run)
rotamalawi_mcmcVw_nr - M-file used for fitting Model 4 to post-vaccination data (NOTE: takes up to 1 week to run)
rotamalV_LL - Function file that calculates log-likelihood of post-vaccination model fit to data

'sensitivity analyses' folder

rasisM1 - Function file containing differential equations for pre-vaccination model 
rotamal_LLnb - Function file that calculates log-likelihood of pre-vaccination model fit to data assuming neg. binomial distribution
rotamalawi_mcmc1-5 - M-file used for fitting model via MCMC to pre-vaccination data assuming different values of the relative infectiousness of subsequent infections from 0.1 to 0.5 (NOTE: takes up to 1 week to run)

MODEL SIMULATION: Contains the code used to simulate the models based on best-fit parameters (estimated using code above) and calculating model output for vaccine effectiveness and impact

rasisM1 - Function file containing differential equations for pre-vaccination model (called by rotagemodelMV, etc.)
rasisV1 - Function file containing differential equations for post-vaccination models 1-2 
rasisV1_direct - Function file containing differential equations for post-vaccination models 1-2 ASSUMING NO REDUCTION IN THE TRANSMISSION RATE following vaccine introduction (i.e. vaccine confers direct protection only)
rasisV1im - Function file containing differential equations for post-vaccination models 1-2 ASSUMING AN IMPROVEMENT IN THE PROPORTION WHO RESPOND TO EACH VACCINE DOSE BEGINNING in 2018
rasisV1redR0 - Function file containing differential equations for post-vaccination models 1-2 ASSUMING A REDUCTION IN R0 BEGINNING in 2018
rasisV1w - Function file containing differential equations for post-vaccination models 3-4
rasisV1_direct - Function file containing differential equations for post-vaccination models 3-4 ASSUMING NO REDUCTION IN THE TRANSMISSION RATE following vaccine introduction (i.e. vaccine confers direct protection only)
rasisV1im - Function file containing differential equations for post-vaccination models 3-4 ASSUMING AN IMPROVEMENT IN THE PROPORTION WHO RESPOND TO EACH VACCINE DOSE BEGINNING in 2018
rasisV1redR0 - Function file containing differential equations for post-vaccination models 3-4 ASSUMING A REDUCTION IN R0 
rotagemodelMV - M-file used to simulate predicted impact of vaccination under Model 1
rotagemodelMV_nr - M-file used to simulate predicted impact of vaccination under Model 2
rotagemodelMV_samp - M-file used to simulate predicted impact of vaccination under Model 1 while sampling from posterior distribution of model parameters (1,000 samples; NOTE: will take a few days to run)
rotagemodelMVdirect - M-file used to simulate predicted impact of vaccination under Model 1 ASSUMING NO REDUCTION IN THE TRANSMISSION RATE following vaccine introduction (i.e. vaccine confers direct protection only)
rotagemodelMVdirect_nr - M-file used to simulate predicted impact of vaccination under Model 2 ASSUMING NO REDUCTION IN THE TRANSMISSION RATE following vaccine introduction (i.e. vaccine confers direct protection only)
rotagemodelMVim - M-file used to simulate predicted impact of vaccination under Model 1 ASSUMING AN IMPROVEMENT IN THE PROPORTION WHO RESPOND TO EACH VACCINE DOSE BEGINNING in 2018
rotagemodelMVim_nr - M-file used to simulate predicted impact of vaccination under Model 2 ASSUMING AN IMPROVEMENT IN THE PROPORTION WHO RESPOND TO EACH VACCINE DOSE BEGINNING in 2018
rotagemodelMVnr_samp - M-file used to simulate predicted impact of vaccination under Model 2 while sampling from posterior distribution of model parameters (1,000 samples; NOTE: will take a few days to run)
rotagemodelMVredR0 - M-file used to simulate predicted impact of vaccination under Model 1 ASSUMING A REDUCTION IN R0 BEGINNING IN 2018
rotagemodelMVredR0_nr - M-file used to simulate predicted impact of vaccination under Model 2 ASSUMING A REDUCTION IN R0 BEGINNING IN 2018
rotagemodelMVredR0w - M-file used to simulate predicted impact of vaccination under Model 3 ASSUMING A REDUCTION IN R0 BEGINNING IN 2018
rotagemodelMVredR0w_nr - M-file used to simulate predicted impact of vaccination under Model 4 ASSUMING A REDUCTION IN R0 BEGINNING IN 2018
rotagemodelMVw - M-file used to simulate predicted impact of vaccination under Model 3 
rotagemodelMVw_direct - M-file used to simulate predicted impact of vaccination under Model 3 ASSUMING NO REDUCTION IN THE TRANSMISSION RATE following vaccine introduction (i.e. vaccine confers direct protection only)
rotagemodelMVw_direct_nr - M-file used to simulate predicted impact of vaccination under Model 4 ASSUMING NO REDUCTION IN THE TRANSMISSION RATE following vaccine introduction (i.e. vaccine confers direct protection only)
rotagemodelMVw_nr - M-file used to simulate predicted impact of vaccination under Model 4 
rotagemodelMVw_samp - M-file used to simulate predicted impact of vaccination under Model 3 while sampling from posterior distribution of model parameters (1,000 samples; NOTE: will take a few days to run)
rotagemodelMVwim - M-file used to simulate predicted impact of vaccination under Model 3 ASSUMING AN IMPROVEMENT IN THE PROPORTION WHO RESPOND TO EACH VACCINE DOSE BEGINNING in 2018
rotagemodelMVwim_nr - M-file used to simulate predicted impact of vaccination under Model 4 ASSUMING AN IMPROVEMENT IN THE PROPORTION WHO RESPOND TO EACH VACCINE DOSE BEGINNING in 2018
rotagemodelMVwnr_samp - M-file used to simulate predicted impact of vaccination under Model 4 while sampling from posterior distribution of model parameters (1,000 samples; NOTE: will take a few days to run)
