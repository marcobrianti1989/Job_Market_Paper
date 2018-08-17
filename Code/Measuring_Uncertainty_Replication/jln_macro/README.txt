Measuring Uncertainty: - Kyle Jurado, Sydney Ludvigson and Serena Ng
September 2014.

E

**** A ******
Files to create macro-uncertainty index using data 1960:1-2011:12

Steps:

1. Prepare jlnrawdata.xlsx.
     Data for 132 Macro Series in Sheet 1
     Data for 157 Financial Series in Sheet 2
     Series 6:15 in Sheet 2 are already included or are highly correlated with variables in Sheet1.
     These variables  are omitted in Step 3, hence 147 series.
2. Transform data and create matlab data: generate_jlndata.m  -> creates jlndata.mat
    The windows and mac version of xlsread treat the data column differently.
   The default is set to windows, line 14.
3. Generate_ferrors.m
      Estimate factors from 132+147 series (See step1 why financial series 6:15 are dropped)
      Compute forecast errors for each of the factors using AR(4) model
      Compute forecast errors for each of the 132 macro series using factor-augmented regression
      output: vyt.txt and vft.txt
4. Estimate stochastic volatility in the forecast errors using 'stochvol' (R package):
     generate_svfdraws.R -> svfmeans.txt
     generate_svydraws.R -> svymeans.txt
5. Use svfmeans.txt and svymenas.txt as input to produce 132 uncertainty estimates
    generate_ut.m -> ut.mat
6. Aggregate ut into macro uncertainty
    generate_aggu.m -> aggu.mat (also includes the time series of uncertainty for each of the 132 series, and for h=1..12)
    


Firm level uncertainty is similar. Files are in the replication_firm folder.    
generate_arut.m  gives ut without factors as predictors
generate_nput.m  gives ut without predictors
    
*** B ****    
Figures:

Figure 1: plot_aggu.m
Figure 2: plot_utpredictor.m
Figure 3: plot_keyseries.m
Figure 4: plot_usp500.m
Figure 5: plot_tsvxo_csa.m
Figure 6: generate_cee05var_csa.m
Figure 7: generate_var_csa_levels.m
Figure 8: plot_firmu.m
Figure 9: generate_cee05var_proxies.m


**** Data ***

The macro data are collected from Global Insight, with call numbers listed in the
supplementary file.
The financial data are mostly taken from Fama-French website, with 5 series taken from CRSP.
The Cochrane-Piassezi factor is based on our calculations.

We cannot make Global-Insight and CRSP available. To implements Step 1, you need
to collect your own data.
We only post jlndata.mat, the transformed data that comes out of Step 2.
