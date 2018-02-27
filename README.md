# bc-fisheries
A repo to test the Fisheries goal model for OHI BC


Scripts are run in this order:

1. 1_data_prep_fis_ram
2. 2_data_prep_fis_ram_spatial
3. 3_calc_fis_scores
4. 4_plot_fis_scores


Supporting scripts that are sourced in the above scripts include:

`kobe_fxns.R` - functions for creating Kobe plots  
`fis_function.R` - Fisheries goal model

The base model we are using considers the following:

1. There is no score weighting by catch in any region
2. We are not applying a penalty for unassessed catch (stocks without B/Bmsy)
3. We are keeping the under/overfishing penalties for assessed stocks