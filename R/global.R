fitdata_varnames_noref <- 
    c("id", "date", "base",
    "base_trend+meanShift", "base_mean_adj", 
    "base_anomaly",             # 6, base - mean annual cycle
    "fit_of_base_anomaly",      # 7
    "annualC+meanShifts",       # 8
    # "deseason_base_trend", "deseason_base_adj_mean",
    # "base_minus_annualC", "fit_of_base_minus_annualC",
    "QM_adjusted",              # 9
    "fit_of_deseason_base")     # 10

# fitdata_varnames_ref <- 
#     c("id", "date", "base", 
#     "meanhatD", "adjB", "meanhat+EB1", "adj", "base",
#     "meanhat", "meanhat+EB", "B", "meanhat0")
fitdata_varnames_ref <- 
    c("id", "date", "base", 
    # "meanhatD", "adjB", 
    "base_trend+meanShift.ref", "base_mean_adj.ref", # shift estimated by `base-ref`
    "base_trend+meanShift.base", "base_mean_adj.base", # shift estimated by deseasonal base
    "base_anomaly", # 7, 8
    "fit_of_base_anomaly", 
    "annualC+meanShifts", #"meanhat+EB",  # 10
    "QM_adjusted",  # B, 11
    "fit_of_deseason_base") # 12


fitprcp_varnames_noref <- 
    c("id", "date", "base",
    "trend+meanShift_org", 
    
    "QM_adjusted" , #5
    "mean_adj_orj", #6

    "trend+meanShift_QM_adjusted", #7
    "BoxCox_transf_org", #8

    "trend+meanShift_BoxCox", #9
    "mean_adj BoxCoxTransf_org" , #10

    "trend+meanShift_mean_adj BoxCox", #11
    "QM_adjusted - base",              # 12, (5-3)
    "mean_adj - base") # 13     # 10
