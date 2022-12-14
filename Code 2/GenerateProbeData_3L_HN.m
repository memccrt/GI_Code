% Genereate good probes for 3 laser-high noise
clear
load SingleProbeVariables.mat %all probes
%load SingleProbeVariables_short.mat %only 1 & 2-way probes
% Single Probes
[OverallMCC_probes_HN_3l,ProbesRemoved_probes_HN_3l,ProbesRemoved_names_probes_HN_3l,Opt_Thresh_probes_HN_3l,Good_probes_probes_HN_3l,R_Cond_final_probes_HN_3l] = RemoveProbes_onebyone(probes,R_3l,high_noise,thresh,names)
% Barcodes 
[OverallMCC_barcodes_HN_3l,ProbesRemoved_barcodes_HN_3l,Opt_Thresh_barcodes_HN_3l,Good_probes_barcodes_HN_3l,R_Cond_final_barcodes_HN_3l,ProbesRemoved_names_barcodes_HN_3l] = RemoveBarcodes_onebyone(Good_probes_probes_HN_3l,R_Cond_final_probes_HN_3l,high_noise,thresh,names);
% Two Barcodes
[OverallMCC_twobarcodes_HN_3l,ProbesRemoved_twobarcodes_HN_3l,Opt_Thresh_twobarcodes_HN_3l,Good_probes_twobarcodes_HN_3l,R_Cond_final_twobarcodes_HN_3l,ProbesRemoved_names_twobarcodes_HN_3l] = RemoveTwoBarcodes_onebyone(Good_probes_barcodes_HN_3l,R_Cond_final_barcodes_HN_3l,high_noise,thresh,names);
save('ProbeData_3l_HN')