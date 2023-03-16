///base model
cd ../../

///run model using cluster100 for each drug
foreach clus in cluster100 cluster12 cluster25 cluster50 {
foreach drug in AMI INH RIF EMB KAN LEV MXF ETH BDQ CFZ LZD RFB DLM {
	
cls
clear

///import data
import delimited Input_data/Initial_model/`drug'_no_interaction_equal_end_short_final.csv

ds

local var_list = r(varlist)

///remove extra variables
local exclude = "uplog2mic log2mic platedesign lineage siteid uniqueid cluster12 cluster25 cluster50 cluster100"

///create final variable list
local var_list : list var_list - exclude

///run mixed effect linear regression with lineage and siteid as fixed effect variables and clusterID as random effect variable. Lineage 4 and site 2 were used as baseline.
meintreg log2mic uplog2mic `var_list' b4.lineage b2.siteid || `clus': , iterate(100)

///save all output
translate @Results Output_files/no_int/STATA/equal_end_range/`clus'/`drug'_no_interaction_equal_end.txt, replace

///write output results to file with statistics
regsave, tstat pval ci
outsheet var coef stderr tstat pval ci_lower ci_upper using Output_files/no_int/STATA/equal_end_range/`clus'/`drug'_no_interaction_equal_end.csv , comma replace

}
}

///rerun model now accounting for interactions

foreach drug in INH EMB CFZ RFB ETH MXF LEV KAN RIF {
foreach clus in cluster100 {
	
cls
clear

//import data with mutation interaction variables
import delimited Input_data/Interaction_model/`drug'_interaction.csv

ds

local var_list = r(varlist)

local exclude = "uplog2mic log2mic platedesign lineage siteid uniqueid cluster12 cluster25 cluster50 cluster100"

local var_list : list var_list - exclude

//run model
meintreg log2mic uplog2mic `var_list' b4.lineage b2.siteid || `clus': , iterate(100)

///save data
translate @Results Output_files/int/STATA/`clus'/`drug'_interaction.txt, replace

regsave, tstat pval ci
outsheet var coef stderr tstat pval ci_lower ci_upper using Output_files/int/STATA/`clus'/`drug'_interaction.csv , comma replace

}
}


///  linear hypothesis testing

///cluster25 cluster50
///AMI INH RIF EMB  KAN LEV MXF ETH BDQ CFZ LZD RFB DLM 
///AMI INH RIF EMB KAN LEV MXF ETH BDQ CFZ LZD RFB DLM

foreach clus in cluster100 {
foreach drug in INH {
	
cls
clear

///import data for linear hypothesis testing
import delimited Input_files/Initial_model/`drug'_no_interaction_equal_end_short_final.csv

ds

local var_list = r(varlist)

local exclude = "uplog2mic log2mic platedesign lineage siteid uniqueid cluster12 cluster25 cluster50 cluster100"

local var_list : list var_list - exclude

///run mixed effect linear regression with lineage and siteid as fixed effect variables and clusterID as random effect variable. Lineage 4 and site 2 were used as baseline.
meintreg log2mic uplog2mic `var_list' b4.lineage b2.siteid || `clus': , iterate(100)

}
}
///run linear hypothesis test
test (inha_i21v=fabg1_c15t) (inha_i21v=fabg1_l203l), mtest(sidak)
