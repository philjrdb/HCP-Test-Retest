n_spss = length(subjID_idx); %subjID pulled from SPSS Idx variable

followup1 = NaN(n_spss,1);
followup2 = NaN(n_spss,1);
followup3 = NaN(n_spss,1);

f1_name = 'slider-response_followup_q1';

%% Loop through SPSS rows
for s = 1:n_spss
  fprintf(['SPSS row ' int2str(s) ' (subjID/HCP row ' int2str(subjID_idx(s)) ')']);
  if ~isempty(HCP_aggr(subjID_idx(s)).raw_Retest)
    is_idx = find(cellfun(@(x) ~isempty(x),HCP_aggr(subjID_idx(s)).raw_Retest(:,1)));

    f1_idx = is_idx(ismember(HCP_aggr(subjID_idx(s)).raw_Retest(is_idx,1),f1_name)); 
    followup1(s) = HCP_aggr(subjID_idx(s)).raw_Retest{f1_idx,3}.val;
    followup2(s) = HCP_aggr(subjID_idx(s)).raw_Retest{f1_idx+1,3}.val;
    followup3(s) = HCP_aggr(subjID_idx(s)).raw_Retest{f1_idx+2,3}.val;
    fprintf(['\n  Mem: ' int2str(followup1(s)) ...
      ' / Change: ' int2str(followup2(s)) ' / Improve: ' int2str(followup2(s)) '\n']);
  else
    fprintf(' - NO RETEST DATA\n');
    followup1(s) = -999;
    followup2(s) = -999;
    followup3(s) = -999;
  end
end