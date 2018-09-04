% function [uptakes] = uptakes(a,b,c1,f1,p1)

%UPTAKES is a function which calculates a vector of uptake rates for different diets by multiplying basal uptake
%rates by random percentages of macronutrients i.e. %carbohydrate (sugars and fibers), %fat and %protein. 
%This vector of uptake rates is then used to inform each novel diet formulation (with adjusted lower bounds) prior to conducting FBA.  

rng(0,'twister');

a = 5; %lower and upper limits for percentages
b = 95;

c = (b-a).*(randi(100,100,1) + a)/100;
c = round(c);
f = (b-a).*(randi(100,100,1) + a)/100;
f = round(f);
p = (b-a).*(randi(100,100,1) + a)/100;
p = round(p);
% 
% 
% c = round(randi(100,1)); %generate random values for carbohydrates, fats and proteins
% f = round(randi(100,1));
% p = round(randi(100,1));
% c1 = c(100);
% f1 = f(100);
% p1 = p(100);
% A = [c1, f1, p1];
% S = sum(A); % calculate sum of three random values
% c1 = [c1/S];  % divide each value by the total sum to normalise (obtain random percentages stored in an array)
% c1 = round(c1*100);
% c1 = c1/100;
% f1 = [f1/S]; 
% f1 = (round(f1*100));
% f1 = f1/100;
% p1 = [p1/S];
% p1 = (round(p1*100));
% p1 = p1/100;
% sumok = c1 + f1 + p1 % check sum = 1.00, i.e. percentages add up to 100 exactly
% %separate array elements back into carbohydrate, fat and protein content
% 
% % load('basal.mat'); %load basal uptake rates
% 
% carbuptake = c1* basal(1:21)/-0.1; % multiply basal rates for relevant macronutrients by their respective random percentages
% fatuptake = f1* basal(22:33)*2/-0.1; % fats and protein rates are multiplied by 2 again to scale up to 100%
% protuptake = p1* basal(34:53)*2/-0.1; % divide by -0.1 to rescale uptake values(since metabolites are taken up, values should be negative)
% 
% alluptakes4 = [carbuptake; fatuptake; protuptake]; % recombine into single vector of uptakes


% For each iteration, incorporate new uptake rates to simulate a new diet

% We start by creating indices for all macronutrient exchanges
% ix_fru = find(ismember(modelJoint.rxnNames,'EX_fru[u]')==1);
% ix_glcD = find(ismember(modelJoint.rxnNames,'EX_glc-D[u]')==1);
% ix_gal = find(ismember(modelJoint.rxnNames,'EX_gal[u]')==1);
% ix_man = find(ismember(modelJoint.rxnNames,'EX_man[u]')==1);
% ix_lcts = find(ismember(modelJoint.rxnNames,'EX_lcts[u]')==1);
% ix_malt = find(ismember(modelJoint.rxnNames,'EX_malt[u]')==1);
% ix_sucr = find(ismember(modelJoint.rxnNames,'EX_sucr[u]')==1);
% ix_melib = find(ismember(modelJoint.rxnNames,'EX_melib[u]')==1);
% ix_strch1 = find(ismember(modelJoint.rxnNames,'EX_strch1[u]')==1);
% ix_inulin = find(ismember(modelJoint.rxnNames,'EX_inulin[u]')==1);
% ix_levan1000 = find(ismember(modelJoint.rxnNames,'EX_levan1000[u]')==1);
% ix_starch1200 = find(ismember(modelJoint.rxnNames,'EX_starch1200[u]')==1);
% ix_arabinogal = find(ismember(modelJoint.rxnNames,'EX_arabinogal[u]')==1);
% ix_pect = find(ismember(modelJoint.rxnNames,'EX_pect[u]')==1);
% ix_pullulan1200 = find(ismember(modelJoint.rxnNames,'EX_pullulan1200[u]')==1)
% ix_amylose300 = find(ismember(modelJoint.rxnNames,'EX_amylose300[u]')==1);
% ix_lmn30 = find(ismember(modelJoint.rxnNames,'EX_lmn30[u]')==1);
% ix_raffin = find(ismember(modelJoint.rxnNames,'EX_raffin[u]')==1);
% ix_stys = find(ismember(modelJoint.rxnNames,'EX_stys[u]')==1);
% ix_oligofru4 = find(ismember(modelJoint.rxnNames,'EX_oligofru4[u]')==1);
% ix_dextran40 = find(ismember(modelJoint.rxnNames,'EX_dextran40[u]')==1);
% ix_arachd = find(ismember(modelJoint.rxnNames,'EX_arachd[u]')==1);
% ix_chsterol = find(ismember(modelJoint.rxnNames,'EX_chsterol[u]')==1);
% ix_glyc = find(ismember(modelJoint.rxnNames,'EX_glyc[u]')==1);
% ix_hdca = find(ismember(modelJoint.rxnNames,'EX_hdca[u]')==1);
% ix_hdcea = find(ismember(modelJoint.rxnNames,'EX_hdcea[u]')==1);
% ix_lnlc = find(ismember(modelJoint.rxnNames,'EX_lnlc[u]')==1);
% ix_lnlnca = find(ismember(modelJoint.rxnNames,'EX_lnlnca[u]')==1);
% ix_lnlncg = find(ismember(modelJoint.rxnNames,'EX_lnlncg[u]')==1);
% ix_ocdca = find(ismember(modelJoint.rxnNames,'EX_ocdca[u]')==1);
% ix_ocdcea = find(ismember(modelJoint.rxnNames,'EX_ocdcea[u]')==1);
% ix_octa = find(ismember(modelJoint.rxnNames,'EX_octa[u]')==1);
% ix_ttdca = find(ismember(modelJoint.rxnNames,'EX_ttdca[u]')==1);
% ix_ala = find(ismember(modelJoint.rxnNames,'EX_ala-L[u]')==1);
% ix_ser = find(ismember(modelJoint.rxnNames,'EX_ser-L[u]')==1);
% ix_cys = find(ismember(modelJoint.rxnNames,'EX_cys-L[u]')==1);
% ix_arg = find(ismember(modelJoint.rxnNames,'EX_arg-L[u]')==1);
% ix_ile = find(ismember(modelJoint.rxnNames,'EX_ile-L[u]')==1);
% ix_leu = find(ismember(modelJoint.rxnNames,'EX_leu-L[u]')==1);
% ix_lys = find(ismember(modelJoint.rxnNames,'EX_lys-L[u]')==1);
% ix_his = find(ismember(modelJoint.rxnNames,'EX_his-L[u]')==1);
% ix_asn = find(ismember(modelJoint.rxnNames,'EX_asn-L[u]')==1);
% ix_asp = find(ismember(modelJoint.rxnNames,'EX_asp-L[u]')==1);
% ix_thr = find(ismember(modelJoint.rxnNames,'EX_thr-L[u]')==1);
% ix_glu = find(ismember(modelJoint.rxnNames,'EX_glu-L[u]')==1);
% ix_met = find(ismember(modelJoint.rxnNames,'EX_met-L[u]')==1);
% ix_gln = find(ismember(modelJoint.rxnNames,'EX_gln-L[u]')==1);
% ix_pro = find(ismember(modelJoint.rxnNames,'EX_pro-L[u]')==1);
% ix_val = find(ismember(modelJoint.rxnNames,'EX_val-L[u]')==1);
% ix_phe = find(ismember(modelJoint.rxnNames,'EX_phe-L[u]')==1);
% ix_tyr = find(ismember(modelJoint.rxnNames,'EX_tyr-L[u]')==1);
% ix_gly = find(ismember(modelJoint.rxnNames,'EX_gly[u]')==1);
% ix_trp = find(ismember(modelJoint.rxnNames,'EX_trp-L[u]')==1);

%% Set new uptake rates in the model using these indices
modelJoint.lb(ix_fru) = alluptakes4(1)
modelJoint.lb(ix_glcD) = alluptakes4(2)
modelJoint.lb(ix_gal) = alluptakes4(3)
modelJoint.lb(ix_man) = alluptakes4(4)
modelJoint.lb(ix_lcts) = alluptakes4(5)
modelJoint.lb(ix_malt) = alluptakes4(6)
modelJoint.lb(ix_sucr) = alluptakes4(7)
modelJoint.lb(ix_melib) = alluptakes4(8)
modelJoint.lb(ix_strch1) = alluptakes4(9)
modelJoint.lb(ix_inulin) = alluptakes4(10)
modelJoint.lb(ix_levan1000) = alluptakes4(11)
modelJoint.lb(ix_starch1200) = alluptakes4(12)
modelJoint.lb(ix_arabinogal) = alluptakes4(13)
modelJoint.lb(ix_pect) = alluptakes4(14)
modelJoint.lb(ix_pullulan1200) = alluptakes4(15)
modelJoint.lb(ix_amylose300) = alluptakes4(16)
modelJoint.lb(ix_lmn30) = alluptakes4(17)
modelJoint.lb(ix_raffin) = alluptakes4(18)
modelJoint.lb(ix_stys) = alluptakes4(19)
modelJoint.lb(ix_oligofru4) = alluptakes4(20)
modelJoint.lb(ix_dextran40) = alluptakes4(21)
modelJoint.lb(ix_arachd) = alluptakes4(22)
modelJoint.lb(ix_chsterol) = alluptakes4(23)
modelJoint.lb(ix_glyc) = alluptakes4(24)
modelJoint.lb(ix_hdca) = alluptakes4(25)
modelJoint.lb(ix_hdcea) = alluptakes4(26)
modelJoint.lb(ix_lnlc) = alluptakes4(27)
modelJoint.lb(ix_lnlnca) = alluptakes4(28)
modelJoint.lb(ix_lnlncg) = alluptakes4(29)
modelJoint.lb(ix_ocdca) = alluptakes4(30)
modelJoint.lb(ix_ocdcea) = alluptakes4(31)
modelJoint.lb(ix_octa) = alluptakes4(32)
modelJoint.lb(ix_ttdca) = alluptakes4(33)
modelJoint.lb(ix_ala) = alluptakes4(34)
modelJoint.lb(ix_ser) = alluptakes4(47)
modelJoint.lb(ix_cys) = alluptakes4(53)
modelJoint.lb(ix_arg) = alluptakes4(35)
modelJoint.lb(ix_ile) = alluptakes4(38)
modelJoint.lb(ix_leu) = alluptakes4(39)
modelJoint.lb(ix_lys) = alluptakes4(40)
modelJoint.lb(ix_his) = alluptakes4(45)
modelJoint.lb(ix_asn) = alluptakes4(43)
modelJoint.lb(ix_asp) = alluptakes4(44)
modelJoint.lb(ix_thr) = alluptakes4(48)
modelJoint.lb(ix_glu) = alluptakes4(36)
modelJoint.lb(ix_met) = alluptakes4(41)
modelJoint.lb(ix_gln) = alluptakes4(51)
modelJoint.lb(ix_pro) = alluptakes4(42)
modelJoint.lb(ix_val) = alluptakes4(52)
modelJoint.lb(ix_phe) = alluptakes4(46)
modelJoint.lb(ix_tyr) = alluptakes4(50)
modelJoint.lb(ix_gly) = alluptakes4(37)
modelJoint.lb(ix_trp) = alluptakes4(49)