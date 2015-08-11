% +ROINETS
%
% Files
%   apply_function_to_structure              - applies function to every field of a structure
%   BayesGLasso_Columnwise                   - Efficient Bayesian Graphical Lasso MCMC sampler 
%   call_fsl_wrapper                         - wrapper on call_fsl function to check for errors. 
%   check_inputs                             - Checks properties of Settings structure
%   cholinv                                  - matrix inverse via cholesky decomposition
%   closest_orthogonal_matrix                - Computes closest orthogonal matrix
%   col_sum                                  - Sum for each column.
%   cols                                     - The number of columns.
%   convert_correlations_to_normal_variables - converts correlations to z-scores
%   convert_precision_to_pcorr               - converts to partial correlation
%   demean                                   - Remove mean value
%   do_group_level_statistics                - add group level inference 
%   do_pairwise_calculation                  - pairwise source-leakage corrected network matrix
%   dp_glasso                                - graphical lasso 
%   envelope_data                            - applies Hilbert envelope to data, without normalisation
%   estimate_AR_coeffs                       - fits an AR model to voxel data and estimates coefficients
%   example                                  - example analysis for a single subject
%   example_external_use                     - EXAMPLE  example analysis for a single subject, if you're external to OHBA
%   example_many_subj                        - EXAMPLE  example analysis for a single subject
%   false_discovery_rate                     - converts standard z-scores to q-scores
%   fast_svds                                - RAM/time-efficient version of SVDS, singular value decomposition
%   find_empirical_dev                       - Build up correlations from empirical null distribtion
%   find_empirical_H0_distribution_width     - FIND_EMPIRICAL_H0_DISTRIBUTON_WIDTH
%   Fisher_r_to_z                            - Converts correlations to z-scores
%   get_node_tcs                             - extracts ROI time-courses
%   glasso_cv                                - K-fold cross-validation for shrinkage parameter in glasso
%   glasso_frequentist                       - graphical lasso for regularized precision matrix estimation
%   householder_orthogonalise                - orthogonalisation using householder method
%   isposdef                                 - Test for positive definite matrix.
%   logdet                                   - Computation of logarithm of determinant of a matrix
%   make_directory                           - Makes directory if not already existing - wrapper on mkdir
%   nii_parcel_quicksave                     - Saves data in parcels as nifti
%   p_to_z_two_tailed                        - Z_TO_P_TWO_TAILED convert p-value to standard z-value in two-tailed test
%   randgamma                                - GC_RANDGAMMA random sample from gamma distribution with shape and scale 
%   reformat_results                         - move session correlation mats to frequency band mats
%   regression                               - Solves multivariate regression y = X*b + e using fast mex binaries
%   remove_source_leakage                    - correct ROI time-courses for source leakage
%   retrieve_analysis_from_cluster           - collects and formats results after
%   row_sum                                  - Sum for each row.
%   rows                                     - The number of rows.
%   run_correlation_analysis                 - runs various correlations on node data
%   run_individual_network_analysis          - runs a single session network analysis
%   scale_cols                               - Scale each column of a matrix.
%   scale_rows                               - Scale each row of a matrix.
%   scomponents                              - Compute the strongly connected components of a graph
%   setdiff_pos_int                          - Set difference of two sets of positive integers (much faster than built-in setdiff)
%   submit_analysis_to_cluster               - submit a set of sessions for network analysis
%   symmetric_orthogonalise                  - closest orthogonal matrix
%   test_me                                  - test that pipeline runs
%   z_to_p_two_tailed                        - convert standard z-value to p-value in two-tailed test
