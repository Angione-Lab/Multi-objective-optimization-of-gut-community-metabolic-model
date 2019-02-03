% % SORT FLUX DISTRIBUTIONS ACCORDING TO SUBSYSTEMS, SPECIES, OR BOTH
% sort fluxes according to subsystems

% list_subsystems = unique(modelJoint.subSystems); % list all subsystems
% ixs_subsystems = cell(numel(list_subsystems),1); % create cell variables
% cardinality_subsystems = zeros(numel(list_subsystems),1);
% sum_flux_subsys= zeros(numel(list_subsystems),1);
% 
% for i=1:numel(list_subsystems)
%     ixs_subsystems{i} = find(strcmp(list_subsystems{i},modelJoint.subSystems));
%     sum_flux_subsys(i) = sum(abs(firstround(ixs_subsystems{i}))); % select each round of fluxes individually e.g. firstround, secondround, thirdround, etc...
%     cardinality_subsystems(i) = numel(ixs_subsystems{i});
% end
% 
% %sort fluxes according to species

% list_species = unique(modelJoint.species); % list all species
% ixs_species = cell(numel(list_species),1); % create cell variables
% cardinality_species = zeros(numel(list_species),1); 
% sum_flux_spec= zeros(numel(list_species),1);

% for i=1:numel(list_species)
%     ixs_species{i} = find(strcmp(list_species{i},modelJoint.species));
%     sum_flux_spec(i) = sum(abs(firstround(ixs_species{i}))); % select each round of fluxes individually e.g. firstround, secondround, thirdround, etc...
%     cardinality_species(i) = numel(ixs_species{i});
% end

%sort flux distribution according to subsystems within species
list_subsystems_by_species = unique(modelJoint.subSystemsbyspecies); % list all species-wise subsystems
ixs_subsystems_by_species = cell(numel(list_subsystems_by_species),1); % create cell variables
cardinality_subsystems_by_species = zeros(numel(list_subsystems_by_species),1);
sum_flux_subsys_by_species= zeros(numel(list_subsystems_by_species),1);
current_vec = firstround(1:28199,1); % select diets 1-100 in turn for each round of fluxes i.e. firstround, secondround, thirdround, fourthround
for i=1:numel(list_subsystems_by_species)
    ixs_subsystems_by_species{i} = find(strcmp(list_subsystems_by_species{i},modelJoint.subSystemsbyspecies));
    sum_flux_subsys_by_species(i) = sum(abs(current_vec(ixs_subsystems_by_species{i})));
    cardinality_subsystems_by_species(i) = numel(ixs_subsystems_by_species{i});
end

averaged_r1 = (sum_flux_subsys_by_species)./(cardinality_subsystems_by_species); % divide sum of fluxes for each subsystem by no. of reactions in that subsystem
% averaged_r2 = (sum_flux_subsys_by_species)./(cardinality_subsystems_by_species); % divide sum of fluxes for each subsystem by no. of reactions in that subsystem
% averaged_r3 = (sum_flux_subsys_by_species)./(cardinality_subsystems_by_species); % divide sum of fluxes for each subsystem by no. of reactions in that subsystem
% averaged_r4 = (sum_flux_subsys_by_species)./(cardinality_subsystems_by_species); % divide sum of fluxes for each subsystem by no. of reactions in that subsystem
                                                                                 
first_corr = corr(averaged_r1(1:573,1:100)',gamma_1); % calculate Pearson correlation for first round of fluxes
second_corr = corr(averaged_r2(1:573,1:100)',gamma_2); %                         ...,,... second round of fluxes
third_corr = corr(averaged_r3(1:573,1:100)',gamma_3); %                          ...,,,... third round of fluxes
fourth_corr = corr(averaged_r4(1:573,1:100)',gamma_4); %                         ...,,,... fourth round of fluxes