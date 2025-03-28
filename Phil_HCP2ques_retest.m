%% HumanCondPun questionnaires

test_names = {'Test1' 'Retest'}; % order of tests 

% questionnaire subscale/overall names (must start with overall scale name as listed in HCP_aggr)
ques_names = {'cfi_alt' 'cfi_ctrl' 'htq_com' 'htq_reg' 'htq_nov' 'audit_haz' 'audit_dep' 'audit_harm' ...
  'cfi_overall' 'htq_overall' 'audit_overall'}; 

% Scoring matrix (items included in each subscale)
ques{1} = [1  0 1  0 1 1  0 1  0 1  0 1 1 1 0 1  0 1 1 1];
ques{2} = [0 -1 0 -1 0 0 -1 0 -1 0 -1 0 0 0 1 0 -1 0 0 0];
ques{3} = [1 1 1 1 0 0 0 0 0 0 0];
ques{4} = [0 0 0 0 1 1 1 1 0 0 0];
ques{5} = [0 0 0 0 0 0 0 0 -1 -1 -1];
ques{6} = [1 1 1 0 0 0 0 0 0 0];
ques{7} = [0 0 0 1 1 1 0 0 0 0];
ques{8} = [0 0 0 0 0 0 1 1 2 2];
ques{9} = [1 -1 1 -1 1 1 -1 1 -1 1 -1 1 1 1 1 1 -1 1 1 1];
ques{10} = [1 1 1 1 1 1 1 1 -1 -1 -1];
ques{11} = [1 1 1 1 1 1 1 1 2 2];

% Options per subscale (n-1 added to negatively scored items to properly reverse score)
rescore = [7 7 7 7 7 5 5 5 7 7 5];

% Adds this to all items (per ques) to restore item base value (assumes raw data scored as 0)
base_score = [1 1 0 0 0 0 0 0 1 0 0];

a_leng = length(HCP_aggr);

%%
for q = 1:length(ques_names)  
  n_idx = regexp(ques_names{q},'_');
  ques_field = ['ques_' ques_names{q}(1:n_idx-1)];
  q_idx = find(ques{q} ~= 0);
  q_scoring = ques{q}(q_idx);

  q_rescore = zeros(size(ques{q}));
  if any(ques{q} < 0)
    r_idx = find(ques{q} < 0);
    q_rescore(r_idx) = rescore(q)-1;
  end

  for t = 1:length(test_names)
  for r = 1:a_leng
    if ~isempty(HCP_aggr(r).([ques_field '_' test_names{t}]))
      HCP_aggr(r).([ques_names{q} '_' test_names{t}]) = ...
        sum(HCP_aggr(r).([ques_field '_' test_names{t}])(q_idx).*q_scoring+q_rescore(q_idx)+base_score(q));
    end
  end
  end
end

clearvars -except HCP* work_folder *keep;