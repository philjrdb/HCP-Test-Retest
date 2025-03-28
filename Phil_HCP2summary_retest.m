%% HCP summary v3
% Run after HCPextract

test_names = {'Test1' 'Retest'}; % order of tests 
n_blocks = 6;
unp_blocks = 1:2;
pun_blocks = 3:5;
pun2_blocks = 6;

bin_var = HCP_parameters.BinLabels;

map_sz = [64 40];
loc_lim = []; % color limits for click heatmap (normalized to # subjs+blocks)

%% Prep variables
clear HCPsum
close all

for t = 1:length(test_names)
  fprintf([' - ' test_names{t} ' -\n']);
  valid_idx = find(~cellfun(@isempty,{HCP_aggr(:).(['raw_' test_names{t}])}));
  groups = unique({HCP_aggr(valid_idx).(['Group_' test_names{t}])});
  n_groups = length(groups);
  
  clear blockstruct
  blockstruct = vertcat(HCP_aggr(valid_idx).(['Block_' test_names{t}]));
  
  n_subj = length(valid_idx);
  exclude_idx = vertcat(HCP_aggr(valid_idx).(['Exclude_' test_names{t}]));

  ValInf_vars = fieldnames(HCP_aggr(valid_idx(1)).(['Block_' test_names{t}]));
  ValInf_vars = ValInf_vars(endsWith(ValInf_vars,{'Val' 'Inf' 'Pred' 'Reval'}));
  
  HCPsum.(test_names{t}).nPunShldAv = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).nUnpShldAv = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).pcPunShldAv = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).pcUnpShldAv = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).PunShldRT = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).UnpShldRT = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).pcPunShldAv_PunB = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).pcUnpShldAv_PunB = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).PunShldRT_PunB = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).UnpShldRT_PunB = NaN(n_subj,n_blocks);
  
  HCPsum.(test_names{t}).nPunShldTaken = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).nUnpShldTaken = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).pcPunShldTaken = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).pcUnpShldTaken = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).pcPunShldTaken_PunB = NaN(n_subj,1);
  HCPsum.(test_names{t}).pcUnpShldTaken_PunB = NaN(n_subj,1);
  HCPsum.(test_names{t}).Pref_pcShldTaken_PunB = NaN(n_subj,1);
  HCPsum.(test_names{t}).PunS_suppr = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).UnpS_suppr = NaN(n_subj,n_blocks);
  HCPsum.(test_names{t}).PunS_suppr_PunB = NaN(n_subj,1);
  HCPsum.(test_names{t}).UnpS_suppr_PunB = NaN(n_subj,1);
  HCPsum.(test_names{t}).PrefS_suppr_PunB = NaN(n_subj,1);
  
  HCPsum.(test_names{t}).EndPoints = vertcat(blockstruct(:).EndPoints);
  
  %% Fill HCPsum
  n_var = length(ValInf_vars);
  for v = 1:n_var
    HCPsum.(test_names{t}).(ValInf_vars{v}) = vertcat(blockstruct.(ValInf_vars{v}));
    
    % Pun Blocks
    HCPsum.(test_names{t}).([ValInf_vars{v} '_PunB']) = mean(HCPsum.(test_names{t}).(ValInf_vars{v})(:,pun_blocks),2);
    HCPsum.(test_names{t}).([ValInf_vars{v} '_PunB2']) = mean(HCPsum.(test_names{t}).(ValInf_vars{v})(:,pun2_blocks),2);
    HCPsum.(test_names{t}).([ValInf_vars{v} '_PunBPunB2ratio']) = HCPsum.(test_names{t}).([ValInf_vars{v} '_PunB2'])./...
      (HCPsum.(test_names{t}).([ValInf_vars{v} '_PunB'])+HCPsum.(test_names{t}).([ValInf_vars{v} '_PunB2']));
    
    % Unp Blocks + PunB:UnpB ratio
    if any(HCPsum.(test_names{t}).(ValInf_vars{v})(:,unp_blocks),'all')
      HCPsum.(test_names{t}).([ValInf_vars{v} '_UnpB']) = mean(HCPsum.(test_names{t}).(ValInf_vars{v})(:,unp_blocks),2);
      HCPsum.(test_names{t}).([ValInf_vars{v} '_PunUnpBratio']) = HCPsum.(test_names{t}).([ValInf_vars{v} '_PunB'])./...
        (HCPsum.(test_names{t}).([ValInf_vars{v} '_PunB'])+HCPsum.(test_names{t}).([ValInf_vars{v} '_UnpB']));
    end
  end
  
  BinTime = cell2mat(vertcat(blockstruct.BinTime));
  zero_idx = BinTime == 0;
  BinPunRate = cell2mat(vertcat(blockstruct.BinPunRate));
  BinPunRate(zero_idx) = 0;
  BinPunClicks = BinPunRate.*BinTime;
  BinUnpRate = cell2mat(vertcat(blockstruct.BinUnpRate));
  BinUnpRate(zero_idx) = 0;
  BinUnpClicks = BinUnpRate.*BinTime;
  for r = 1:n_subj
    blockstruct(r).BinShieldRate{1} = zeros(13,1);
    blockstruct(r).BinShieldRate{2} = zeros(13,1);
  end
  BinShieldRate = cell2mat(vertcat(blockstruct.BinShieldRate));
  BinShieldRate(zero_idx) = 0;
  BinShieldClicks = BinShieldRate.*BinTime;
  
  n_var = length(bin_var);
  for v = 1:n_var
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time']) = BinTime(v:13:end,:);
    
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunClicks']) = BinPunClicks(v:13:end,:);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate']) = BinPunRate(v:13:end,:);
    
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpClicks']) = BinUnpClicks(v:13:end,:);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate']) = BinUnpRate(v:13:end,:);
    
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ShieldClicks']) = BinShieldClicks(v:13:end,:);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ShieldRate']) = BinShieldRate(v:13:end,:);
    
    % Pun Blocks
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time_PunB']) = sum(BinTime(v:13:end,pun_blocks),2);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time_PunB2']) = sum(BinTime(v:13:end,pun2_blocks),2);
    
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunClicks_PunB']) = sum(BinPunClicks(v:13:end,pun_blocks),2);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunClicks_PunB'])...
      ./HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time_PunB']);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunClicks_PunB2']) = sum(BinPunClicks(v:13:end,pun2_blocks),2);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB2']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunClicks_PunB2'])...
      ./HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time_PunB2']);
    
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpClicks_PunB']) = sum(BinUnpClicks(v:13:end,pun_blocks),2);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpClicks_PunB'])...
      ./HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time_PunB']);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpClicks_PunB2']) = sum(BinUnpClicks(v:13:end,pun2_blocks),2);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB2']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpClicks_PunB2'])...
      ./HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time_PunB2']);
    
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ShieldClicks_PunB']) = sum(BinShieldClicks(v:13:end,pun_blocks),2);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ShieldRate_PunB']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ShieldClicks_PunB'])...
      ./HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time_PunB']);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ShieldClicks_PunB2']) = sum(BinShieldClicks(v:13:end,pun2_blocks),2);
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ShieldRate_PunB2']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ShieldClicks_PunB2'])...
      ./HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time_PunB2']);
  
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunBPunB2ratio']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB2'])...
      ./(HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB'])+HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB2']));
  
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunBPunB2ratio']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB2'])...
      ./(HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])+HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB2']));
  
    if v == 1 % ITI for UnpB + PunB:UnpB ratio
      % Unp Blocks
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time_UnpB']) = sum(BinTime(v:13:end,unp_blocks),2);
  
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunClicks_UnpB']) = sum(BinPunClicks(v:13:end,unp_blocks),2);
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_UnpB']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunClicks_UnpB'])...
        ./HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time_UnpB']);
  
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpClicks_UnpB']) = sum(BinUnpClicks(v:13:end,unp_blocks),2);
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_UnpB']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpClicks_UnpB'])...
        ./HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time_UnpB']);
  
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ShieldClicks_UnpB']) = sum(BinShieldClicks(v:13:end,unp_blocks),2);
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ShieldRate_UnpB']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ShieldClicks_UnpB'])...
        ./HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_Time_UnpB']);
  
      % Pun:Unp Block ratio
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunUnpBratio']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB'])...
        ./(HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB'])+HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_UnpB']));
  
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunUnpBratio']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])...
        ./(HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])+HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_UnpB']));
    else % Bin:ITI ratio   
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ITIratio_Pun']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate'])./...
        (HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate']) + HCPsum.(test_names{t}).ITI_PunRate);
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ITIratio_Unp']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate'])./...
        (HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate']) + HCPsum.(test_names{t}).ITI_UnpRate);
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ITIratio_Both']) = ...
        (HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate'])+HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate']))./...
        (HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate'])+HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate'])...
        + HCPsum.(test_names{t}).ITI_PunRate+HCPsum.(test_names{t}).ITI_UnpRate);
  
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ITIratio_Pun_PunB']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB'])./...
        (HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB']) + HCPsum.(test_names{t}).ITI_PunRate_PunB);
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ITIratio_Unp_PunB']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])./...
        (HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB']) + HCPsum.(test_names{t}).ITI_UnpRate_PunB);
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ITIratio_Both_PunB']) = ...
        (HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB'])+HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB']))./...
        (HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])+HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB'])...
        + HCPsum.(test_names{t}).ITI_PunRate_PunB + HCPsum.(test_names{t}).ITI_UnpRate_PunB);
      
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ITIratio_Pun_PunB2']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB2'])./...
        (HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB2']) + HCPsum.(test_names{t}).ITI_PunRate_PunB2);
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ITIratio_Unp_PunB2']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB2'])./...
        (HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB2']) + HCPsum.(test_names{t}).ITI_UnpRate_PunB2);
      HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_ITIratio_Both_PunB2']) = ...
        (HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB2'])+HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB2']))./...
        (HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB2'])+HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB2'])...
        + HCPsum.(test_names{t}).ITI_PunRate_PunB2 + HCPsum.(test_names{t}).ITI_UnpRate_PunB2);  
    end
    
    % UnpITI ratio
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunUnpITIratio']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB'])...
      ./(HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_PunRate_PunB'])+HCPsum.(test_names{t}).ITI_PunRate_UnpB);
  
    HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunUnpITIratio']) = HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])...
      ./(HCPsum.(test_names{t}).([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])+HCPsum.(test_names{t}).ITI_UnpRate_UnpB);
  end
  
  %% Unshielded rates + suppr
  HCPsum.(test_names{t}).PunS_Unshld_Time = (BinTime(2:13:end,:)+BinTime(3:13:end,:)+BinTime(5:13:end,:));
  HCPsum.(test_names{t}).PunS_Unshld_PunRate = ...
    (BinPunClicks(2:13:end,:)+BinPunClicks(3:13:end,:)+BinPunClicks(5:13:end,:))./HCPsum.(test_names{t}).PunS_Unshld_Time;
  HCPsum.(test_names{t}).PunS_Unshld_UnpRate = ...
    (BinUnpClicks(2:13:end,:)+BinUnpClicks(3:13:end,:)+BinUnpClicks(5:13:end,:))./HCPsum.(test_names{t}).PunS_Unshld_Time;
  HCPsum.(test_names{t}).PunS_Unshld_ShieldRate = ...
    (BinShieldClicks(2:13:end,:)+BinShieldClicks(3:13:end,:)+BinShieldClicks(5:13:end,:))./HCPsum.(test_names{t}).PunS_Unshld_Time;
  
  HCPsum.(test_names{t}).PunS_Unshld_PunRate_PunB = ...
    sum((BinPunClicks(2:13:end,pun_blocks)+BinPunClicks(3:13:end,pun_blocks)+BinPunClicks(5:13:end,pun_blocks)),2)...
    ./sum(HCPsum.(test_names{t}).PunS_Unshld_Time(:,pun_blocks),2);
  HCPsum.(test_names{t}).PunS_Unshld_UnpRate_PunB = ...
    sum((BinUnpClicks(2:13:end,pun_blocks)+BinUnpClicks(3:13:end,pun_blocks)+BinUnpClicks(5:13:end,pun_blocks)),2)...
    ./sum(HCPsum.(test_names{t}).PunS_Unshld_Time(:,pun_blocks),2);
  HCPsum.(test_names{t}).PunS_Unshld_PunRate_PunB2 = ...
    sum((BinPunClicks(2:13:end,pun2_blocks)+BinPunClicks(3:13:end,pun2_blocks)+BinPunClicks(5:13:end,pun2_blocks)),2)...
    ./sum(HCPsum.(test_names{t}).PunS_Unshld_Time(:,pun2_blocks),2);
  HCPsum.(test_names{t}).PunS_Unshld_UnpRate_PunB2 = ...
    sum((BinUnpClicks(2:13:end,pun2_blocks)+BinUnpClicks(3:13:end,pun2_blocks)+BinUnpClicks(5:13:end,pun2_blocks)),2)...
    ./sum(HCPsum.(test_names{t}).PunS_Unshld_Time(:,pun2_blocks),2);
  
  HCPsum.(test_names{t}).UnpS_Unshld_Time = (BinTime(8:13:end,:)+BinTime(9:13:end,:)+BinTime(11:13:end,:));
  HCPsum.(test_names{t}).UnpS_Unshld_PunRate = ...
    (BinPunClicks(8:13:end,:)+BinPunClicks(9:13:end,:)+BinPunClicks(11:13:end,:))./HCPsum.(test_names{t}).UnpS_Unshld_Time;
  HCPsum.(test_names{t}).UnpS_Unshld_UnpRate = ...
    (BinUnpClicks(8:13:end,:)+BinUnpClicks(9:13:end,:)+BinUnpClicks(11:13:end,:))./HCPsum.(test_names{t}).UnpS_Unshld_Time;
  HCPsum.(test_names{t}).UnpS_Unshld_ShieldRate = ...
    (BinShieldClicks(8:13:end,:)+BinShieldClicks(9:13:end,:)+BinShieldClicks(11:13:end,:))./HCPsum.(test_names{t}).UnpS_Unshld_Time;
  
  HCPsum.(test_names{t}).UnpS_Unshld_PunRate_PunB = ...
    sum((BinPunClicks(8:13:end,pun_blocks)+BinPunClicks(9:13:end,pun_blocks)+BinPunClicks(11:13:end,pun_blocks)),2)...
    ./sum(HCPsum.(test_names{t}).UnpS_Unshld_Time(:,pun_blocks),2);
  HCPsum.(test_names{t}).UnpS_Unshld_UnpRate_PunB = ...
    sum((BinUnpClicks(8:13:end,pun_blocks)+BinUnpClicks(9:13:end,pun_blocks)+BinUnpClicks(11:13:end,pun_blocks)),2)...
    ./sum(HCPsum.(test_names{t}).UnpS_Unshld_Time(:,pun_blocks),2);
  HCPsum.(test_names{t}).UnpS_Unshld_PunRate_PunB2 = ...
    sum((BinPunClicks(8:13:end,pun2_blocks)+BinPunClicks(9:13:end,pun2_blocks)+BinPunClicks(11:13:end,pun2_blocks)),2)...
    ./sum(HCPsum.(test_names{t}).UnpS_Unshld_Time(:,pun2_blocks),2);
  HCPsum.(test_names{t}).UnpS_Unshld_UnpRate_PunB2 = ...
    sum((BinUnpClicks(8:13:end,pun2_blocks)+BinUnpClicks(9:13:end,pun2_blocks)+BinUnpClicks(11:13:end,pun2_blocks)),2)...
    ./sum(HCPsum.(test_names{t}).UnpS_Unshld_Time(:,pun2_blocks),2);
  
  HCPsum.(test_names{t}).PunS_suppr = (HCPsum.(test_names{t}).PunS_Unshld_PunRate+HCPsum.(test_names{t}).PunS_Unshld_UnpRate)./...
    ((HCPsum.(test_names{t}).PunS_Unshld_PunRate+HCPsum.(test_names{t}).PunS_Unshld_UnpRate)+(HCPsum.(test_names{t}).ITI_PunRate+HCPsum.(test_names{t}).ITI_UnpRate));
  HCPsum.(test_names{t}).UnpS_suppr = (HCPsum.(test_names{t}).UnpS_Unshld_PunRate+HCPsum.(test_names{t}).UnpS_Unshld_UnpRate)./...
    ((HCPsum.(test_names{t}).UnpS_Unshld_PunRate+HCPsum.(test_names{t}).UnpS_Unshld_UnpRate)+(HCPsum.(test_names{t}).ITI_PunRate+HCPsum.(test_names{t}).ITI_UnpRate));
  HCPsum.(test_names{t}).PunS_suppr_PunB = (HCPsum.(test_names{t}).PunS_Unshld_PunRate_PunB+HCPsum.(test_names{t}).PunS_Unshld_UnpRate_PunB)./...
    (HCPsum.(test_names{t}).PunS_Unshld_PunRate_PunB+HCPsum.(test_names{t}).PunS_Unshld_UnpRate_PunB+...
    HCPsum.(test_names{t}).ITI_PunRate_PunB+HCPsum.(test_names{t}).ITI_UnpRate_PunB);
  HCPsum.(test_names{t}).UnpS_suppr_PunB = (HCPsum.(test_names{t}).UnpS_Unshld_PunRate_PunB+HCPsum.(test_names{t}).UnpS_Unshld_UnpRate_PunB)./...
    (HCPsum.(test_names{t}).UnpS_Unshld_PunRate_PunB+HCPsum.(test_names{t}).UnpS_Unshld_UnpRate_PunB+...
    HCPsum.(test_names{t}).ITI_PunRate_PunB+HCPsum.(test_names{t}).ITI_UnpRate_PunB);
  HCPsum.(test_names{t}).PrefS_suppr_PunB = HCPsum.(test_names{t}).PunS_suppr_PunB./(HCPsum.(test_names{t}).PunS_suppr_PunB+HCPsum.(test_names{t}).UnpS_suppr_PunB);
  
  %% ITI ratios
  HCPsum.(test_names{t}).Pref = HCPsum.(test_names{t}).ITI_PunRate./(HCPsum.(test_names{t}).ITI_PunRate+HCPsum.(test_names{t}).ITI_UnpRate);
  HCPsum.(test_names{t}).Pref_UnpB = mean(HCPsum.(test_names{t}).Pref(:,unp_blocks),2);
  HCPsum.(test_names{t}).Pref_PunB = mean(HCPsum.(test_names{t}).Pref(:,pun_blocks),2);
  HCPsum.(test_names{t}).Pref_PunB2 = mean(HCPsum.(test_names{t}).Pref(:,pun2_blocks),2);
  HCPsum.(test_names{t}).Pref_PunUnpBratio = HCPsum.(test_names{t}).Pref_PunB./(HCPsum.(test_names{t}).Pref_UnpB+HCPsum.(test_names{t}).Pref_PunB);
  HCPsum.(test_names{t}).Pref_PunBPunB2ratio = HCPsum.(test_names{t}).Pref_PunB2./(HCPsum.(test_names{t}).Pref_PunB+HCPsum.(test_names{t}).Pref_PunB2);
  
  HCPsum.(test_names{t}).PrefVal = HCPsum.(test_names{t}).PunPlanetVal./(HCPsum.(test_names{t}).PunPlanetVal+HCPsum.(test_names{t}).UnpPlanetVal);
  HCPsum.(test_names{t}).PrefVal_UnpB = mean(HCPsum.(test_names{t}).PrefVal(:,unp_blocks),2);
  HCPsum.(test_names{t}).PrefVal_PunB = mean(HCPsum.(test_names{t}).PrefVal(:,pun_blocks),2);
  HCPsum.(test_names{t}).PrefVal_PunB2 = mean(HCPsum.(test_names{t}).PrefVal(:,pun2_blocks),2);
  HCPsum.(test_names{t}).PrefVal_PunUnpBratio = HCPsum.(test_names{t}).PrefVal_PunB./(HCPsum.(test_names{t}).PrefVal_UnpB+HCPsum.(test_names{t}).PrefVal_PunB);
  HCPsum.(test_names{t}).PrefVal_PunBPunB2ratio = HCPsum.(test_names{t}).PrefVal_PunB2./(HCPsum.(test_names{t}).PrefVal_PunB2+HCPsum.(test_names{t}).PrefVal_PunB);
  
  HCPsum.(test_names{t}).PrefRewInf = HCPsum.(test_names{t}).PunPlanet_RewInf./(HCPsum.(test_names{t}).PunPlanet_RewInf+HCPsum.(test_names{t}).UnpPlanet_RewInf);
  HCPsum.(test_names{t}).PrefRewInf_UnpB = mean(HCPsum.(test_names{t}).PrefRewInf(:,unp_blocks),2);
  HCPsum.(test_names{t}).PrefRewInf_PunB = mean(HCPsum.(test_names{t}).PrefRewInf(:,pun_blocks),2);
  HCPsum.(test_names{t}).PrefRewInf_PunB2 = mean(HCPsum.(test_names{t}).PrefRewInf(:,pun2_blocks),2);
  HCPsum.(test_names{t}).PrefRewInf_PunUnpBratio = HCPsum.(test_names{t}).PrefRewInf_PunB./(HCPsum.(test_names{t}).PrefRewInf_UnpB+HCPsum.(test_names{t}).PrefRewInf_PunB);
  HCPsum.(test_names{t}).PrefRewInf_PunBPunB2ratio = HCPsum.(test_names{t}).PrefRewInf_PunB2./(HCPsum.(test_names{t}).PrefRewInf_PunB2+HCPsum.(test_names{t}).PrefRewInf_PunB);
  
  %HCPsum.(test_names{t}).PrefRewInfconf = HCPsum.(test_names{t}).PunPlanet_RewInfconf./(HCPsum.(test_names{t}).PunPlanet_RewInfconf+HCPsum.(test_names{t}).UnpPlanet_RewInfconf);
  %HCPsum.(test_names{t}).PrefRewInfconf_UnpB = mean(HCPsum.(test_names{t}).PrefRewInfconf(:,unp_blocks),2);
  %HCPsum.(test_names{t}).PrefRewInfconf_PunB = mean(HCPsum.(test_names{t}).PrefRewInfconf(:,pun_blocks),2);
  %HCPsum.(test_names{t}).PrefRewInfconf_PunB2 = mean(HCPsum.(test_names{t}).PrefRewInfconf(:,pun2_blocks),2);
  %HCPsum.(test_names{t}).PrefRewInfconf_PunUnpBratio = HCPsum.(test_names{t}).PrefRewInfconf_PunB./(HCPsum.(test_names{t}).PrefRewInfconf_UnpB+HCPsum.(test_names{t}).PrefRewInfconf_PunB);
  %HCPsum.(test_names{t}).PrefRewInfconf_PunBPunB2ratio = HCPsum.(test_names{t}).PrefRewInfconf_PunB2./(HCPsum.(test_names{t}).PrefRewInfconf_PunB2+HCPsum.(test_names{t}).PrefRewInfconf_PunB);
  
  HCPsum.(test_names{t}).PrefPAtkInf = HCPsum.(test_names{t}).PunPlanet_AttackInf./(HCPsum.(test_names{t}).PunPlanet_AttackInf+HCPsum.(test_names{t}).UnpPlanet_AttackInf);
  HCPsum.(test_names{t}).PrefPAtkInf_PunB = mean(HCPsum.(test_names{t}).PrefPAtkInf(:,pun_blocks),2);
  HCPsum.(test_names{t}).PrefPAtkInf_PunB2 = mean(HCPsum.(test_names{t}).PrefPAtkInf(:,pun2_blocks),2);
   
  HCPsum.(test_names{t}).PunTsuppr = HCPsum.(test_names{t}).ITI_PunRate./(HCPsum.(test_names{t}).ITI_PunRate+HCPsum.(test_names{t}).ITI_PunRate(:,2));
  HCPsum.(test_names{t}).UnpTsuppr = HCPsum.(test_names{t}).ITI_UnpRate./(HCPsum.(test_names{t}).ITI_UnpRate+HCPsum.(test_names{t}).ITI_UnpRate(:,2));
  HCPsum.(test_names{t}).PunValChange = HCPsum.(test_names{t}).PunPlanetVal./(HCPsum.(test_names{t}).PunPlanetVal+HCPsum.(test_names{t}).PunPlanetVal(:,2));
  HCPsum.(test_names{t}).UnpValChange = HCPsum.(test_names{t}).UnpPlanetVal./(HCPsum.(test_names{t}).UnpPlanetVal+HCPsum.(test_names{t}).UnpPlanetVal(:,2));
  
  %% Ships & Shields (n, percent, rt)
  HCPsum.(test_names{t}).nPunShips = cellfun(@(x) size(x,1), vertcat(blockstruct.PunShields));
  PunShields = cellfun(@(x) sum(~isnan(x),1), vertcat(blockstruct.PunShields),'UniformOutput',false);
  idx = cellfun(@(x) ~isempty(x), PunShields);
  HCPsum.(test_names{t}).nPunShldAv(~idx) = 0;
  HCPsum.(test_names{t}).nPunShldAv(idx) = cellfun(@(x) x(1), PunShields(idx));
  HCPsum.(test_names{t}).nPunShldTaken(idx) = cellfun(@(x) x(2), PunShields(idx));
  HCPsum.(test_names{t}).nPunShldTaken(~idx) = 0;
  HCPsum.(test_names{t}).pcPunShldAv = HCPsum.(test_names{t}).nPunShldAv./HCPsum.(test_names{t}).nPunShips;
  HCPsum.(test_names{t}).pcPunShldTaken = HCPsum.(test_names{t}).nPunShldTaken./HCPsum.(test_names{t}).nPunShldAv;
  PunShields = cellfun(@(x) mean(x,1,'omitnan'),vertcat(blockstruct.PunShields),'UniformOutput',false);
  HCPsum.(test_names{t}).PunShldRT(idx) = cellfun(@(x) x(2), PunShields(idx));
  
  HCPsum.(test_names{t}).pcPunShldAv_PunB = sum(HCPsum.(test_names{t}).nPunShldAv(:,pun_blocks),2)./sum(HCPsum.(test_names{t}).nPunShips(:,pun_blocks),2);
  HCPsum.(test_names{t}).pcPunShldTaken_PunB = sum(HCPsum.(test_names{t}).nPunShldTaken(:,pun_blocks),2)./sum(HCPsum.(test_names{t}).nPunShldAv(:,pun_blocks),2);
  HCPsum.(test_names{t}).PunShldRT_PunB = sum((HCPsum.(test_names{t}).PunShldRT(:,pun_blocks).*HCPsum.(test_names{t}).nPunShldTaken(:,pun_blocks)),2,'omitnan')...
    ./sum(HCPsum.(test_names{t}).nPunShldTaken(:,pun_blocks),2);
  HCPsum.(test_names{t}).pcPunShldAv_PunB2 = sum(HCPsum.(test_names{t}).nPunShldAv(:,pun2_blocks),2)./sum(HCPsum.(test_names{t}).nPunShips(:,pun2_blocks),2);
  HCPsum.(test_names{t}).pcPunShldTaken_PunB2 = sum(HCPsum.(test_names{t}).nPunShldTaken(:,pun2_blocks),2)./sum(HCPsum.(test_names{t}).nPunShldAv(:,pun2_blocks),2);
  HCPsum.(test_names{t}).PunShldRT_PunB2 = sum((HCPsum.(test_names{t}).PunShldRT(:,pun2_blocks).*HCPsum.(test_names{t}).nPunShldTaken(:,pun2_blocks)),2,'omitnan')...
    ./sum(HCPsum.(test_names{t}).nPunShldTaken(:,pun2_blocks),2);
  
  HCPsum.(test_names{t}).nUnpShips = cellfun(@(x) size(x,1), vertcat(blockstruct.UnpShields));
  UnpShields = cellfun(@(x) sum(~isnan(x),1), vertcat(blockstruct.UnpShields),'UniformOutput',false);
  idx = cellfun(@(x) ~isempty(x), UnpShields);
  HCPsum.(test_names{t}).nUnpShldAv(idx) = cellfun(@(x) x(1), UnpShields(idx));
  HCPsum.(test_names{t}).nUnpShldAv(~idx) = 0;
  HCPsum.(test_names{t}).nUnpShldTaken(idx) = cellfun(@(x) x(2), UnpShields(idx));
  HCPsum.(test_names{t}).nUnpShldTaken(~idx) = 0;
  HCPsum.(test_names{t}).pcUnpShldAv = HCPsum.(test_names{t}).nUnpShldAv./HCPsum.(test_names{t}).nUnpShips;
  HCPsum.(test_names{t}).pcUnpShldTaken = HCPsum.(test_names{t}).nUnpShldTaken./HCPsum.(test_names{t}).nUnpShldAv;
  UnpShields = cellfun(@(x) mean(x,1,'omitnan'),vertcat(blockstruct.UnpShields),'UniformOutput',false);
  HCPsum.(test_names{t}).UnpShldRT(idx) = cellfun(@(x) x(2), UnpShields(idx));
  
  HCPsum.(test_names{t}).pcUnpShldAv_PunB = sum(HCPsum.(test_names{t}).nUnpShldAv(:,pun_blocks),2)./sum(HCPsum.(test_names{t}).nUnpShips(:,pun_blocks),2);
  HCPsum.(test_names{t}).pcUnpShldTaken_PunB = sum(HCPsum.(test_names{t}).nUnpShldTaken(:,pun_blocks),2)./sum(HCPsum.(test_names{t}).nUnpShldAv(:,pun_blocks),2);
  HCPsum.(test_names{t}).UnpShldRT_PunB = mean((HCPsum.(test_names{t}).UnpShldRT(:,pun_blocks).*HCPsum.(test_names{t}).nUnpShldTaken(:,pun_blocks)),2,'omitnan')...
    ./sum(HCPsum.(test_names{t}).nUnpShldTaken(:,pun_blocks),2);
  HCPsum.(test_names{t}).pcUnpShldAv_PunB2 = sum(HCPsum.(test_names{t}).nUnpShldAv(:,pun2_blocks),2)./sum(HCPsum.(test_names{t}).nUnpShips(:,pun2_blocks),2);
  HCPsum.(test_names{t}).pcUnpShldTaken_PunB2 = sum(HCPsum.(test_names{t}).nUnpShldTaken(:,pun2_blocks),2)./sum(HCPsum.(test_names{t}).nUnpShldAv(:,pun2_blocks),2);
  HCPsum.(test_names{t}).UnpShldRT_PunB2 = mean((HCPsum.(test_names{t}).UnpShldRT(:,pun2_blocks).*HCPsum.(test_names{t}).nUnpShldTaken(:,pun2_blocks)),2,'omitnan')...
    ./sum(HCPsum.(test_names{t}).nUnpShldTaken(:,pun2_blocks),2);
  
  %% Planet-Attack chain inference estimate
  HCPsum.(test_names{t}).PunP_AtkInfEst = ...
    (HCPsum.(test_names{t}).PunPlanet_PunShipInf.*HCPsum.(test_names{t}).PunShip_AttackInf/100)+...
    (HCPsum.(test_names{t}).PunPlanet_UnpShipInf.*HCPsum.(test_names{t}).UnpShip_AttackInf/100);
  HCPsum.(test_names{t}).UnpP_AtkInfEst = ...
    (HCPsum.(test_names{t}).UnpPlanet_PunShipInf.*HCPsum.(test_names{t}).PunShip_AttackInf/100)+...
    (HCPsum.(test_names{t}).UnpPlanet_UnpShipInf.*HCPsum.(test_names{t}).UnpShip_AttackInf/100);
  
  HCPsum.(test_names{t}).PunP_PunS_AtkInfEst = HCPsum.(test_names{t}).PunPlanet_PunShipInf.*HCPsum.(test_names{t}).PunShip_AttackInf/100;
  HCPsum.(test_names{t}).PunP_UnpS_AtkInfEst = HCPsum.(test_names{t}).PunPlanet_UnpShipInf.*HCPsum.(test_names{t}).UnpShip_AttackInf/100;
  HCPsum.(test_names{t}).UnpP_UnpS_AtkInfEst = HCPsum.(test_names{t}).UnpPlanet_UnpShipInf.*HCPsum.(test_names{t}).UnpShip_AttackInf/100;
  HCPsum.(test_names{t}).UnpP_PunS_AtkInfEst = HCPsum.(test_names{t}).UnpPlanet_PunShipInf.*HCPsum.(test_names{t}).PunShip_AttackInf/100;
  
  HCPsum.(test_names{t}).PunP_AtkInfEst_cap = HCPsum.(test_names{t}).PunP_AtkInfEst;
  HCPsum.(test_names{t}).PunP_AtkInfEst_cap(HCPsum.(test_names{t}).PunP_AtkInfEst_cap>100) = 100;
  HCPsum.(test_names{t}).UnpP_AtkInfEst_cap = HCPsum.(test_names{t}).UnpP_AtkInfEst;
  HCPsum.(test_names{t}).UnpP_AtkInfEst_cap(HCPsum.(test_names{t}).UnpP_AtkInfEst>100) = 100;
  
  %% Figures
  for n = 1:n_groups
    %% Click locations
    clicks = vertcat(blockstruct.NormClickLoc);
    
    incl_idx = ismember(vertcat({HCP_aggr(valid_idx).(['Group_' test_names{t}])})',groups(n)) & ~exclude_idx;
    
    left_idx = ismember({HCP_aggr(valid_idx).(['PunPlanet_' test_names{t}])},'left')';
  
    avgUnpBleftP_clicks = vertcat(clicks{left_idx&incl_idx,unp_blocks});
    avgUnpBleftP_clickmap = histcounts2(avgUnpBleftP_clicks(:,1),avgUnpBleftP_clicks(:,2)...
      ,-1:(2/map_sz(1)):1,-1:(2/map_sz(2)):1)/(length(unp_blocks)*sum(left_idx&incl_idx));
    avgUnpBrightP_clicks = vertcat(clicks{~left_idx&incl_idx,unp_blocks});
    avgUnpBrightP_clickmap = histcounts2(avgUnpBrightP_clicks(:,1),avgUnpBrightP_clicks(:,2)...
      ,-1:(2/map_sz(1)):1,-1:(2/map_sz(2)):1)/(length(unp_blocks)*sum(~left_idx&incl_idx));
    avgPunBleftP_clicks = vertcat(clicks{left_idx&incl_idx,pun_blocks});
    avgPunBleftP_clickmap = histcounts2(avgPunBleftP_clicks(:,1),avgPunBleftP_clicks(:,2)...
      ,-1:(2/map_sz(1)):1,-1:(2/map_sz(2)):1)/(length(pun_blocks)*sum(left_idx&incl_idx));
    avgPunBrightP_clicks = vertcat(clicks{~left_idx&incl_idx,pun_blocks});
    avgPunBrightP_clickmap = histcounts2(avgPunBrightP_clicks(:,1),avgPunBrightP_clicks(:,2)...
      ,-1:(2/map_sz(1)):1,-1:(2/map_sz(2)):1)/(length(pun_blocks)*sum(~left_idx&incl_idx));
    avgPunB2leftP_clicks = vertcat(clicks{left_idx&incl_idx,pun2_blocks});
    avgPunB2leftP_clickmap = histcounts2(avgPunB2leftP_clicks(:,1),avgPunB2leftP_clicks(:,2)...
      ,-1:(2/map_sz(1)):1,-1:(2/map_sz(2)):1)/(length(pun2_blocks)*sum(left_idx&incl_idx));
    avgPunB2rightP_clicks = vertcat(clicks{~left_idx&incl_idx,pun2_blocks});
    avgPunB2rightP_clickmap = histcounts2(avgPunB2rightP_clicks(:,1),avgPunB2rightP_clicks(:,2)...
      ,-1:(2/map_sz(1)):1,-1:(2/map_sz(2)):1)/(length(pun2_blocks)*sum(~left_idx&incl_idx));
  
    if isempty(loc_lim)
      minofmax = min([max(avgUnpBleftP_clickmap,[],'all') max(avgUnpBrightP_clickmap,[],'all') ...
        max(avgPunBleftP_clickmap,[],'all') max(avgPunBrightP_clickmap,[],'all')]);
      loc_lim = [0 minofmax];
    end
  
    figure; 
    subplot(2,2,1); hold on
    imagesc(avgUnpBleftP_clickmap',loc_lim)
    %hist3(vertcat(clicks{left_idx&~exclude_idx,unp_blocks}),'Ctrs',{-1:(2/80):1 -1:(2/50):1},'CdataMode','auto')
    title(['Group' groups{n} ' Pre-pun average click locations (PunPlanet=left)']);
  
    subplot(2,2,2); hold on
    imagesc(avgPunBleftP_clickmap',loc_lim)
    %hist3(vertcat(clicks{left_idx&~exclude_idx,pun_blocks}),'Ctrs',{-1:(2/80):1 -1:(2/50):1},'CdataMode','auto')
    title(['Group' groups{n} ' Pun phase average click locations (PunPlanet=left)']);
  
    subplot(2,2,3); hold on
    imagesc(avgUnpBrightP_clickmap',loc_lim)
    %hist3(vertcat(clicks{~left_idx&~exclude_idx,unp_blocks}),'Ctrs',{-1:(2/80):1 -1:(2/50):1},'CdataMode','auto')
    title(['Group' groups{n} ' Pre-pun average click locations (PunPlanet=right)']);
  
    subplot(2,2,4); hold on
    imagesc(avgPunBrightP_clickmap',loc_lim)
    %hist3(vertcat(clicks{~left_idx&~exclude_idx,pun_blocks}),'Ctrs',{-1:(2/80):1 -1:(2/50):1},'CdataMode','auto')
    title(['Group' groups{n} ' Pun phase average click locations (PunPlanet=right)']);
  
    for p = 1:4
      subplot(2,2,p); hold on
      axis tight off 
      colormap('hot');
      colorbar;
    end
  
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[work_folder '\Figs\' test_names{t} ' Group' groups{n} ' Click heatmap.png']);
    close all
  
    % Location change
    LvRdiff_leftP_UnpB = avgUnpBleftP_clickmap(1:map_sz(1)/2,:)-flipud(avgUnpBleftP_clickmap(map_sz(1)/2+1:end,:));
    LvRdiff_leftP_PunB = avgPunBleftP_clickmap(1:map_sz(1)/2,:)-flipud(avgPunBleftP_clickmap(map_sz(1)/2+1:end,:));
    LvRdiff_rightP_UnpB = avgUnpBrightP_clickmap(1:map_sz(1)/2,:)-flipud(avgUnpBrightP_clickmap(map_sz(1)/2+1:end,:));
    LvRdiff_rightP_PunB = avgPunBrightP_clickmap(1:map_sz(1)/2,:)-flipud(avgPunBrightP_clickmap(map_sz(1)/2+1:end,:));
    minofmax = max([abs(min([max(LvRdiff_rightP_UnpB,[],'all') max(LvRdiff_rightP_PunB,[],'all')])) ...
      abs(max([min(LvRdiff_rightP_UnpB,[],'all') min(LvRdiff_rightP_PunB,[],'all')]))]);
    loc_lim = [-minofmax minofmax];
  
    figure; 
    subplot(2,2,1); hold on
    imagesc(vertcat(LvRdiff_leftP_UnpB,LvRdiff_leftP_PunB)',loc_lim)
    title(['Group' groups{n} ' Left-Right [UnpBlock / PunBlock] (PunPlanet=left)'],'interpreter','none');
    plot([map_sz(1)/2+0.5 map_sz(1)/2+0.5],[0.5 map_sz(2)+0.5],'k--');
    drawnow;
  
    subplot(2,2,3); hold on
    imagesc(vertcat(LvRdiff_rightP_UnpB,LvRdiff_rightP_PunB)',loc_lim)
    title(['Group' groups{n} ' Left-Right [UnpBlock / PunBlock] (PunPlanet=right)'],'interpreter','none');
    plot([map_sz(1)/2+0.5 map_sz(1)/2+0.5],[0.5 map_sz(2)+0.5],'k--');
    drawnow;
  
    UnpPunDiff_leftP = avgPunBleftP_clickmap-avgUnpBleftP_clickmap;
    UnpPunDiff_rightP = avgPunBrightP_clickmap-avgUnpBrightP_clickmap;
    minofmax = max([abs(min([max(UnpPunDiff_leftP,[],'all') max(UnpPunDiff_rightP,[],'all')])) ...
      abs(max([min(UnpPunDiff_leftP,[],'all') min(UnpPunDiff_rightP,[],'all')]))]);
    loc_lim = [-minofmax minofmax];
  
    subplot(2,2,2); hold on
    imagesc(UnpPunDiff_leftP',loc_lim)
    title(['Group' groups{n} ' Pun-UnpBlock click location (PunPlanet=left)'],'interpreter','none');
  
    subplot(2,2,4); hold on
    imagesc(UnpPunDiff_rightP',loc_lim)
    title(['Group' groups{n} ' Pun-UnpBlock click location (PunPlanet=right)'],'interpreter','none');
    for p = 1:4
      subplot(2,2,p); hold on
      axis tight off 
      colormap(whitejet(100));
      colorbar;
    end
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[work_folder '\Figs\' test_names{t} ' Group' groups{n} ' click change heatmap.png']);
    close all
  
  
    %% Planet clicks (per sess)
    % ITI rates
    figure; 
    subplot(2,3,1); hold on
    errorbar(mean(HCPsum.(test_names{t}).ITI_PunRate(incl_idx,:)*60,1),sem(HCPsum.(test_names{t}).ITI_PunRate(incl_idx,:)*60),...
      'k','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).ITI_UnpRate(incl_idx,:)*60,1),sem(HCPsum.(test_names{t}).ITI_UnpRate(incl_idx,:)*60),...
      'k','CapSize',2);
    p1 = plot(mean(HCPsum.(test_names{t}).ITI_PunRate(incl_idx,:)*60,1),'Color',col_rep(2),'LineWidth',2);
    p2 = plot(mean(HCPsum.(test_names{t}).ITI_UnpRate(incl_idx,:)*60,1),'Color',col_rep(3),'LineWidth',2);
    xlabel('Block');
    xlim([0.5 n_blocks+0.5]);
    xticks(1:n_blocks);
    ylims = ylim;
    ylims(1) = 0;
    ylim(ylims);
    ylabel('Clicks/min','interpreter','none');
    title(['Group' groups{n} ' ITI planet click rates'])
    legend([p1 p2],{['PunPlanet (n=' num2str(sum(incl_idx)) ')'] 'UnpPlanet'},...
      'Location','best','interpreter','none');
    legend('boxoff');
  
    % T suppr
    subplot(2,3,2); hold on
    title('ITI ratio (relative to Blk2)');
    incl_idx2 = incl_idx & ~isnan(HCPsum.(test_names{t}).PunTsuppr(:,2)) & ~isnan(HCPsum.(test_names{t}).UnpTsuppr(:,2));
    errorbar(mean(HCPsum.(test_names{t}).PunTsuppr(incl_idx2,:),1),sem(HCPsum.(test_names{t}).PunTsuppr(incl_idx2,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).UnpTsuppr(incl_idx2,:),1),sem(HCPsum.(test_names{t}).UnpTsuppr(incl_idx2,:)),...
      'k.','CapSize',2);
    plot(HCPsum.(test_names{t}).PunTsuppr(incl_idx2,:)','Color',col_rep(2)+((1-col_rep(2)).*.5),'LineWidth',0.5);
    plot(HCPsum.(test_names{t}).UnpTsuppr(incl_idx2,:)','Color',col_rep(3)+((1-col_rep(3)).*.5),'LineWidth',0.5);
    p1 = plot(mean(HCPsum.(test_names{t}).PunTsuppr(incl_idx2,:),1),'Color',col_rep(2),'LineWidth',2);
    p2 = plot(mean(HCPsum.(test_names{t}).UnpTsuppr(incl_idx2,:),1),'Color',col_rep(3),'LineWidth',2);
    xlabel('Block');
    xlim([0.5 n_blocks+0.5]);
    xticks(1:n_blocks);
    plot(xlim,[0.5 0.5],'k:');
    ylim([0 1]);
    legend([p1 p2],{['Pun T_suppr (n=' num2str(sum(incl_idx2)) ')'] 'Unp T_suppr'},...
      'Location','best','interpreter','none');
    legend('boxoff');
  
    % Ratios
    subplot(2,3,3); hold on
    title('Suppression Ratios');
    plot(HCPsum.(test_names{t}).Pref(incl_idx,:)','Color',col_rep(14),'LineWidth',0.5);
    errorbar(mean(HCPsum.(test_names{t}).Pref(incl_idx,:),1),sem(HCPsum.(test_names{t}).Pref(incl_idx,:)),...
      'k.','CapSize',2);
    p1 = plot(mean(HCPsum.(test_names{t}).Pref(incl_idx,:),1),'k','LineWidth',1.5);
    incl_idx2 = incl_idx & ~any(isnan(HCPsum.(test_names{t}).PunS_suppr(:,[pun_blocks, pun2_blocks])),2);
    errorbar(mean(HCPsum.(test_names{t}).PunS_suppr(incl_idx2,:),1),sem(HCPsum.(test_names{t}).PunS_suppr(incl_idx2,:)),...
      'k.','CapSize',2);
    p2 = plot(mean(HCPsum.(test_names{t}).PunS_suppr(incl_idx2,:),1),'d--','Color',col_rep(2),'LineWidth',2,...
      'MarkerFaceColor',col_rep(2));
    incl_idx2 = incl_idx & ~any(isnan(HCPsum.(test_names{t}).UnpS_suppr(:,[pun_blocks, pun2_blocks])),2);
    errorbar(mean(HCPsum.(test_names{t}).UnpS_suppr(incl_idx2,:),1),sem(HCPsum.(test_names{t}).UnpS_suppr(incl_idx2,:)),...
      'k.','CapSize',2);
    p3 = plot(mean(HCPsum.(test_names{t}).UnpS_suppr(incl_idx2,:),1),'o--','Color',col_rep(3),'LineWidth',2,...
      'MarkerFaceColor',col_rep(3));
    incl_idx2 = incl_idx & ~any(isnan(HCPsum.(test_names{t}).pcPunShldTaken(:,[pun_blocks, pun2_blocks])),2);
    errorbar(mean(HCPsum.(test_names{t}).pcPunShldTaken(incl_idx2,:),1,'omitnan'),sem(HCPsum.(test_names{t}).pcPunShldTaken(incl_idx2,:)),...
      'k.','CapSize',2);
    p4 = plot(mean(HCPsum.(test_names{t}).pcPunShldTaken(incl_idx2,:),1,'omitnan'),'d-','Color',col_rep(1),'LineWidth',2,...
      'MarkerFaceColor',col_rep(2));
    incl_idx2 = incl_idx & ~any(isnan(HCPsum.(test_names{t}).pcUnpShldTaken(:,[pun_blocks, pun2_blocks])),2);
    errorbar(mean(HCPsum.(test_names{t}).pcUnpShldTaken(incl_idx2,:),1,'omitnan'),sem(HCPsum.(test_names{t}).pcUnpShldTaken(incl_idx2,:)),...
      'k.','CapSize',2);
    p5 = plot(mean(HCPsum.(test_names{t}).pcUnpShldTaken(incl_idx2,:),1,'omitnan'),'o-','Color',col_rep(1),'LineWidth',2,...
      'MarkerFaceColor',col_rep(3));
    xlabel('Block');
    xticks(1:n_blocks);
    xlim([0.5 n_blocks+0.5]);
    plot(xlim,[0.5 0.5],'k:');
    ylim([0 1]);
    ylabel('Ratio');
    legend([p1 p2 p3 p4 p5],{'ITI PunP v UnpP' 'PunS Suppr' 'UnpS Suppr' 'PunS Shld %' 'UnpS Shld %'},...
      'Location','best','interpreter','none');
    legend('boxoff')
  
    % Values
    subplot(2,3,4); hold on
    title('Values');
    errorbar(mean(HCPsum.(test_names{t}).PunPlanetVal(incl_idx,:),1),sem(HCPsum.(test_names{t}).PunPlanetVal(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).UnpPlanetVal(incl_idx,:),1),sem(HCPsum.(test_names{t}).UnpPlanetVal(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).PunShipVal(incl_idx,:),1),sem(HCPsum.(test_names{t}).PunShipVal(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).UnpShipVal(incl_idx,:),1),sem(HCPsum.(test_names{t}).UnpShipVal(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).RewVal(incl_idx,:),1),sem(HCPsum.(test_names{t}).RewVal(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).AttackVal(incl_idx,:),1),sem(HCPsum.(test_names{t}).AttackVal(incl_idx,:)),...
      'k.','CapSize',2);
    p1 = plot(mean(HCPsum.(test_names{t}).PunPlanetVal(incl_idx,:),1),'Color',col_rep(2),'LineWidth',2);
    p2 = plot(mean(HCPsum.(test_names{t}).UnpPlanetVal(incl_idx,:),1),'Color',col_rep(3),'LineWidth',2);
    p3 = plot(mean(HCPsum.(test_names{t}).PunShipVal(incl_idx,:),1),'d--','Color',col_rep(2),'LineWidth',2,...
      'MarkerFaceColor',col_rep(2));
    p4 = plot(mean(HCPsum.(test_names{t}).UnpShipVal(incl_idx,:),1),'o--','Color',col_rep(3),'LineWidth',2,...
      'MarkerFaceColor',col_rep(3));
    p5 = plot(mean(HCPsum.(test_names{t}).RewVal(incl_idx,:),1),'k+:','LineWidth',1.5,'MarkerSize',10);
    p6 = plot(mean(HCPsum.(test_names{t}).AttackVal(incl_idx,:),1),'x-','Color',col_rep(2),'LineWidth',1.5,...
      'MarkerFaceColor',col_rep(2),'MarkerSize',10);
    xlabel('Block');
    xlim([0.5 n_blocks+0.5]);
    xticks(1:n_blocks);
    ylim([0 100]);
    ylabel('Value rating');
    legend([p1 p2 p3 p4 p5 p6],{'PunPlanet' 'UnpPlanet' 'PunShip' 'UnpShip' 'Reward' 'Attack'},...
      'Location','best','interpreter','none');
    legend('boxoff');
  
    % Inferences
    subplot(2,3,5); hold on
    title('Inferences');
    errorbar(mean(HCPsum.(test_names{t}).PunPlanet_RewInf(incl_idx,:),1),sem(HCPsum.(test_names{t}).PunPlanet_RewInf(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).UnpPlanet_RewInf(incl_idx,:),1),sem(HCPsum.(test_names{t}).UnpPlanet_RewInf(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).PunPlanet_PunShipInf(incl_idx,:),1),sem(HCPsum.(test_names{t}).PunPlanet_PunShipInf(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).UnpPlanet_UnpShipInf(incl_idx,:),1),sem(HCPsum.(test_names{t}).UnpPlanet_UnpShipInf(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).PunPlanet_UnpShipInf(incl_idx,:),1),sem(HCPsum.(test_names{t}).PunPlanet_UnpShipInf(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).UnpPlanet_PunShipInf(incl_idx,:),1),sem(HCPsum.(test_names{t}).UnpPlanet_PunShipInf(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).PunShip_AttackInf(incl_idx,:),1),sem(HCPsum.(test_names{t}).PunShip_AttackInf(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).UnpShip_AttackInf(incl_idx,:),1),sem(HCPsum.(test_names{t}).UnpShip_AttackInf(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).PunPlanet_AttackInf(incl_idx,:),1),sem(HCPsum.(test_names{t}).PunPlanet_AttackInf(incl_idx,:)),...
      'k.','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).UnpPlanet_AttackInf(incl_idx,:),1),sem(HCPsum.(test_names{t}).UnpPlanet_AttackInf(incl_idx,:)),...
      'k.','CapSize',2);
    p1 = plot(mean(HCPsum.(test_names{t}).PunPlanet_RewInf(incl_idx,:),1),'+:','Color',col_rep(2),...
      'LineWidth',1.5,'MarkerSize',10);
    p2 = plot(mean(HCPsum.(test_names{t}).UnpPlanet_RewInf(incl_idx,:),1),'+:','Color',col_rep(3),...
      'LineWidth',1.5,'MarkerSize',10);
    p3 = plot(mean(HCPsum.(test_names{t}).PunPlanet_PunShipInf(incl_idx,:),1),'d--','Color',col_rep(2),...
      'LineWidth',2,'MarkerFaceColor',col_rep(2));
    p9 = plot(mean(HCPsum.(test_names{t}).PunPlanet_UnpShipInf(incl_idx,:),1),'o--','Color',col_rep(2),...
      'LineWidth',2,'MarkerFaceColor',col_rep(3));
    p4 = plot(mean(HCPsum.(test_names{t}).UnpPlanet_UnpShipInf(incl_idx,:),1),'o--','Color',col_rep(3),...
      'LineWidth',2,'MarkerFaceColor',col_rep(3));
    p10 = plot(mean(HCPsum.(test_names{t}).UnpPlanet_PunShipInf(incl_idx,:),1),'d--','Color',col_rep(3),...
      'LineWidth',2,'MarkerFaceColor',col_rep(2));
    p5 = plot(mean(HCPsum.(test_names{t}).PunShip_AttackInf(incl_idx,:),1),'x--','Color',col_rep(2),...
      'LineWidth',1.5,'MarkerSize',10);
    p6 = plot(mean(HCPsum.(test_names{t}).UnpShip_AttackInf(incl_idx,:),1),'x--','Color',col_rep(3),...
      'LineWidth',1.5,'MarkerSize',10);
    p7 = plot(mean(HCPsum.(test_names{t}).PunPlanet_AttackInf(incl_idx,:),1),'x:','Color',col_rep(2),...
      'LineWidth',1.5,'MarkerSize',10);
    p8 = plot(mean(HCPsum.(test_names{t}).UnpPlanet_AttackInf(incl_idx,:),1),'x:','Color',col_rep(3),...
      'LineWidth',1.5,'MarkerSize',10);
    xlabel('Block');
    xlim([0.5 n_blocks+0.5]);
    xticks(1:n_blocks);
    ylim([0 100]);
    ylabel('% likelihood');
    legend([p1 p2 p3 p9 p4 p10 p5 p6 p7 p8],...
      {'PunP-Rew' 'UnpP-Rew' 'PunP-PunS' 'PunP-UnpS' 'UnpP-UnpS' 'UnpP-PunS' ...
      'PunS-Attack' 'UnpS-Attack' 'PunP-Attack' 'UnpP-Attack'},...
      'Location','best','interpreter','none');
    legend('boxoff');
  
    % R1:R2 prediction/revaluation
    subplot(2,3,6); hold on
    title('R1:R2 prediction/revaluation','interpreter','none');
    errorbar(mean(HCPsum.(test_names{t}).PrefPred(incl_idx,:),1),sem(HCPsum.(test_names{t}).PrefPred(incl_idx,:)),...
      'k','CapSize',2);
    errorbar(mean(HCPsum.(test_names{t}).PrefReval(incl_idx,:),1),sem(HCPsum.(test_names{t}).PrefReval(incl_idx,:)),...
      'k','CapSize',2);
    p1 = plot(mean(HCPsum.(test_names{t}).PrefPred(incl_idx,:),1),'k','LineWidth',2);
    p2 = plot(mean(HCPsum.(test_names{t}).PrefReval(incl_idx,:),1),'Color',col_rep(1),'LineWidth',2);
    xlabel('Block');
    xlim([0.5 n_blocks+0.5]);
    xticks(1:n_blocks);
    ylim([0 100]);
    plot(xlim,[50 50],'k:');
    ylabel('Pref %','interpreter','none');
    legend([p1 p2],{'Pref prediction' 'Pref reval'},...
      'Location','best','interpreter','none');
    legend('boxoff');
  
    set(gcf,'Position',get(0,'Screensize'));
    saveas(figure(1),[work_folder '\Figs\' test_names{t} ' Group' groups{n} ' Aggregate summary.png']);
    close all
  end % n (group) loop
end % t (test) loop

clearvars -except HCP* work_folder *keep;
fprintf('Done.\n');