%% HCP Matlab-SPSS unravel
% Puts HCPsum and questionnaire variables 

test_names = {'Test1' 'Retest'}; % order of tests 

n_blocks = 6;
pun_trials = 3:6;

RT_lowcutoff = .8; % valid RT lower cut-off in secs (any RT/#qs < low cutoff = invalid RT)
RT_highcutoff = 30; % valid RT upper cut-off in secs (any RT/#qs > high cutoff = invalid RT)

ques_names = {'cfi_alt' 'cfi_ctrl' 'htq_com' 'htq_reg' 'htq_nov' 'audit_haz' 'audit_dep' 'audit_harm' ...
  'cfi_overall' 'htq_overall' 'audit_overall'}; % questionnaire names

%% Prep variables
clear SPSS_labels SPSS_mat

ag_leng = length(HCP_aggr); 

f_names = fieldnames(HCP_aggr);
q_vars = union(f_names(contains(f_names,'ques')),ques_names);
q_leng = length(q_vars);

vars = fieldnames(HCPsum.(test_names{1}));
v_leng = length(vars);
n_tests = length(test_names);
n_pre_labels = 1+n_tests*12;

SPSS_labels{n_pre_labels+1+(v_leng*n_blocks)*n_tests+q_leng} = [];
SPSS_mat = NaN(ag_leng,n_pre_labels+1+(v_leng*n_blocks)*n_tests+q_leng);

%% Demographics/parameters
SPSS_labels{1} = 'Idx';
SPSS_mat(:,1) = [1:ag_leng]';

for t=1:n_tests
  SPSS_labels{t+1} = [test_names{t} '_Idx'];
  SPSS_mat(:,t+1) = ~cellfun(@isempty,{HCP_aggr(:).(['raw_' test_names{t}])});
  valid_idx = find(SPSS_mat(:,t+1));
  
  SPSS_labels{t+n_tests+1} = ['json_' test_names{t}];
  u_list = unique({HCP_aggr(valid_idx).(['json_' test_names{t}])},'stable');
  for u = 1:length(u_list)
    SPSS_mat(valid_idx(ismember({HCP_aggr(valid_idx).(['json_' test_names{t}])},u_list{u})),t+n_tests+1) = u;
  end
  
  SPSS_labels{t+n_tests*2+1} = ['Gender_' test_names{t}]; 
  SPSS_mat(valid_idx(ismember({HCP_aggr(valid_idx).(['Gender_' test_names{t}])},'male')),t+n_tests*2+1) = 1;
  SPSS_mat(valid_idx(ismember({HCP_aggr(valid_idx).(['Gender_' test_names{t}])},'female')),t+n_tests*2+1) = 2;
  SPSS_mat(valid_idx(ismember({HCP_aggr(valid_idx).(['Gender_' test_names{t}])},'other')),t+n_tests*2+1) = 3;
  
  SPSS_labels{t+n_tests*3+1} = ['Age_' test_names{t}];
  SPSS_mat(valid_idx,t+n_tests*3+1) = vertcat(HCP_aggr(valid_idx).(['Age_' test_names{t}]));
  
  SPSS_labels{t+n_tests*4+1} = ['Language_' test_names{t}]; % 1 = includes english, 2 = other
  idx = cellfun(@(x) any(x),(regexpi({HCP_aggr(valid_idx).(['Language_' test_names{t}])},'english')));
  SPSS_mat(valid_idx(idx),t+n_tests*4+1) = 1;
  SPSS_mat(valid_idx(~idx),t+n_tests*4+1) = 2;
  
  SPSS_labels{t+n_tests*5+1} = ['PunPlanet_' test_names{t}]; % 1 = left, 2 = right
  SPSS_mat(valid_idx(ismember({HCP_aggr(valid_idx).(['PunPlanet_' test_names{t}])},'left')),t+n_tests*5+1) = 1;
  SPSS_mat(valid_idx(ismember({HCP_aggr(valid_idx).(['PunPlanet_' test_names{t}])},'right')),t+n_tests*5+1) = 2;
  
  SPSS_labels{t+n_tests*6+1} = ['PunShip_' test_names{t}]; % 1 = TypeI, 2 = TypeII
  SPSS_mat(valid_idx(ismember({HCP_aggr(valid_idx).(['PunShip_' test_names{t}])},'TypeI')),t+n_tests*6+1) = 1;
  SPSS_mat(valid_idx(ismember({HCP_aggr(valid_idx).(['PunShip_' test_names{t}])},'TypeII')),t+n_tests*6+1) = 2;
  
  SPSS_labels{t+n_tests*7+1} = ['Catch_' test_names{t}]; % 1 = failed catch questions
  SPSS_mat(valid_idx,t+n_tests*7+1) = any(horzcat(HCP_aggr(valid_idx).(['Catch_' test_names{t}])),1);
  
  SPSS_labels{t+n_tests*8+1} = ['MinValInfRT_' test_names{t}];
  SPSS_mat(valid_idx,t+n_tests*8+1) = cellfun(@(x) min(x,[],'all','omitnan'),...
    {HCP_aggr(valid_idx).(['ValInf_AvgRT_' test_names{t}])});
  
  SPSS_labels{t+n_tests*9+1} = ['MaxValInfRT_' test_names{t}];
  SPSS_mat(valid_idx,t+n_tests*9+1) = cellfun(@(x) max(x,[],'all','omitnan'),...
    {HCP_aggr(valid_idx).(['ValInf_AvgRT_' test_names{t}])});
  
  SPSS_labels{t+n_tests*10+1} = ['Exclude_' test_names{t}];
  SPSS_mat(valid_idx,t+n_tests*10+1) = any(horzcat(SPSS_mat(valid_idx,t+n_tests*7+1),...
    SPSS_mat(valid_idx,t+n_tests*8+1)<RT_lowcutoff,SPSS_mat(valid_idx,t+n_tests*9+1)>RT_highcutoff),2);
  
  try
    groups = unique({HCP_aggr(valid_idx).(['Group_' test_names{t}])});
    n_groups = length(groups);
  
    SPSS_labels{t+n_tests*11+1} = ['Group' strjoin(groups,'_') '_' test_names{t}];
    for n = 1:n_groups
      SPSS_mat(valid_idx(ismember(vertcat({HCP_aggr(valid_idx).(['Group_' test_names{t}])})',...
        groups(n))),t+n_tests*11+1) = n;
    end
  catch
    warning('No groups');
  end
end % test loop

cur_col = 1+n_pre_labels; % # previous SPSS columns
  
%% HCPsum vars
for t = 1:n_tests
  fprintf([' - ' test_names{t} ' -\n']);
  valid_idx = find(SPSS_mat(:,t+1));
  n_val = length(valid_idx);

  vars = fieldnames(HCPsum.(test_names{t}));
  v_leng = length(vars);

  for v = 1:v_leng
    if size(HCPsum.(test_names{t}).(vars{v}),2) == n_blocks 
      SPSS_mat(valid_idx,cur_col+1:cur_col+n_blocks) = HCPsum.(test_names{t}).(vars{v});
  
      for p = 1:n_blocks
        SPSS_labels{cur_col+p} = [vars{v} '_B' num2str(p) '_' test_names{t}];
      end
      
      cur_col = cur_col+n_blocks;
    elseif size(HCPsum.(test_names{t}).(vars{v}),2) == 1 ...
        && size(HCPsum.(test_names{t}).(vars{v}),1) == n_val  
      SPSS_mat(valid_idx,cur_col+1) = HCPsum.(test_names{t}).(vars{v});
      
      SPSS_labels{cur_col+1} = [vars{v} '_' test_names{t}];
      cur_col = cur_col+1;
    end
  end
  
  %% Questionnaire vars
  % Questionnaire items
  q_vars = f_names(contains(f_names,'ques') & contains(f_names,test_names{t}));
  q_leng = length(q_vars);

  for v = 1:q_leng
    var_size = size(HCP_aggr(valid_idx(1)).(q_vars{v}));
  
    if var_size(1)>var_size(2)
      SPSS_mat(valid_idx,cur_col+1:cur_col+var_size(1)) = ...
        horzcat(HCP_aggr(valid_idx).(q_vars{v}))';
  
	    SPSS_labels{cur_col+var_size(1)} = []; 
      for s = 1:var_size(1)
        SPSS_labels{cur_col+s} = [q_vars{v} '_' int2str(s)];
      end
  
      cur_col = cur_col+var_size(1);
  
    elseif var_size(2)>var_size(1)
      SPSS_mat(valid_idx,cur_col+1:cur_col+var_size(2)) = ...
        vertcat(HCP_aggr(valid_idx).(q_vars{v}));
  
      SPSS_labels{cur_col+var_size(2)} = [];
      for s = 1:var_size(2)
        SPSS_labels{cur_col+s} = [q_vars{v} '_' int2str(s)];
      end
  
      cur_col = cur_col+var_size(2);
      
    elseif var_size(2) == 1
      SPSS_labels{cur_col+1} = [q_vars{v} '_' test_names{t}];
      SPSS_mat(valid_idx,cur_col+1) = vertcat(HCP_aggr(valid_idx).(q_vars{v}));
 
      cur_col = cur_col+1;
    end
  end % v (questionnaire) loop

  % Questionnaire subscales
  q_vars = f_names(contains(f_names,ques_names) & contains(f_names,test_names{t}));
  q_leng = length(q_vars);

  for v = 1:q_leng
    SPSS_labels{cur_col+1} = [q_vars{v}];
    SPSS_mat(valid_idx,cur_col+1) = vertcat(HCP_aggr(valid_idx).(q_vars{v}));

    cur_col = cur_col+1;
  end % v (questionnaire) loop
end % test loop

%% Make empty/NaN cells -999
any_idx = any(SPSS_mat,1);
SPSS_labels = SPSS_labels(any_idx);
SPSS_mat = SPSS_mat(:,any_idx);
SPSS_mat(isnan(SPSS_mat)) = -999;
SPSS_mat(isempty(SPSS_mat)) = -999;

% Make excluded data = -999
SPSS_mat_inclOnly = SPSS_mat;
for t = 1:n_tests
  excl_idx = SPSS_mat(:,ismember(SPSS_labels,['Exclude_' test_names{t}])) == 1;
  test_idx = contains(SPSS_labels,test_names{t});
  SPSS_mat_inclOnly(excl_idx,test_idx) = -999;
end

%% Write into .csv
T = array2table(SPSS_mat,'VariableNames',SPSS_labels);
writetable(T,[work_folder '\SPSSdata.csv']);

T_inclOnly = array2table(SPSS_mat_inclOnly,'VariableNames',SPSS_labels);
writetable(T_inclOnly,[work_folder '\SPSSdata_inclOnly.csv']);

clearvars -except HCP* T* work_folder *keep SPSS*;

fprintf('Saving... ');
save([work_folder '\HCP.mat']);

fprintf(['\nSaved (' char(datetime) '). ']);
diary OFF