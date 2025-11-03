bagel_slope_varknown.cpp         # C++: slope model, known variance (seasonality supported)
bagel_slope_varunknown.cpp       # C++: slope model, unknown variance (seasonality supported, unfinished)
KL_approx_varknown.cpp           # C++: KL-based pruning helpers for slope model with known var
find_min_err_index_list_matrix.cpp # C++: exact pruning helper for the mean_change_known_model
bagel_slope_varknown.R           # R reference: slope, known variance
mean_known_paper_matrix_version.R# R reference: mean change only

Still need to write:
- KL_approx_varunknown.cpp
- mean_known.cpp
- mean_varunknown (by changing the design matrix in slope varunknown h=(1,t) to h=(1,1))