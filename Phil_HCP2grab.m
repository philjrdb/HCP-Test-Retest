%% Phil_HCPgrab
% Grabs rows meeting criterion from HCPsum

%%
n_blocks = 6;
pun_trials = 3:6;

RT_lowcutoff = .8; % valid RT lower cut-off in secs (any RT/#qs < low cutoff = invalid RT)
RT_highcutoff = 30; % valid RT upper cut-off in secs (any RT/#qs > high cutoff = invalid RT)

exclude_idx = any(horzcat(SPSS_mat(:,9),SPSS_mat(:,10)<RT_lowcutoff,SPSS_mat(:,11)>RT_highcutoff),2);

grab = HCPsum;
flds = fieldnames(grab);

for f = 1:length(flds)
  grab.(flds{f}) = grab.(flds{f})(~exclude_idx,:);
end

ques_names = {'cfi_alt' 'cfi_ctrl' 'htq_com' 'htq_reg' 'htq_nov' 'audit_haz' 'audit_dep' 'audit_harm' ...
  'cfi_overall' 'htq_overall' 'audit_overall'}; % questionnaire names
for f = 1:length(ques_names)
  data = vertcat(HCP_aggr(:).(ques_names{f}));
  grab.(ques_names{f}) = data(~exclude_idx,:);
end

clearvars -except grab T HCP* work_folder *keep SPSS*;
save([work_folder '\HCP.mat']);

%%
% PunShldRT{66,5} = [];
% selected_PunShldRT{4} = [];
% UnpShldRT{66,5} = [];
% selected_UnpShldRT{4} = [];
% 
% data = vertcat(blockstruct.PunShields);
% idx = cellfun(@(x) ~isempty(x), data);
% PunShldRT(idx) = cellfun(@(x) x(:,2), data(idx),'UniformOutput',false);
% PunShldRT = PunShldRT(~exclude_idx,:);
% for b = 1:3
%   selected_PunShldRT{b} = vertcat(PunShldRT{:,pun_trials(b)});
%   selected_PunShldRT{b} = selected_PunShldRT{b}(~isnan(selected_PunShldRT{b}));
% end
% selected_PunShldRT{4} = vertcat(PunShldRT{:});
% selected_PunShldRT{4} = selected_PunShldRT{4}(~isnan(selected_PunShldRT{4}));
% 
% data = vertcat(blockstruct.UnpShields);
% idx = cellfun(@(x) ~isempty(x), data);
% UnpShldRT(idx) = cellfun(@(x) x(:,2), data(idx),'UniformOutput',false);
% UnpShldRT = UnpShldRT(~exclude_idx,:);
% for b = 1:3
%   selected_UnpShldRT{b} = vertcat(UnpShldRT{:,pun_trials(b)});
%   selected_UnpShldRT{b} = selected_UnpShldRT{b}(~isnan(selected_UnpShldRT{b}));
% end
% selected_UnpShldRT{4} = vertcat(UnpShldRT{:});
% selected_UnpShldRT{4} = selected_UnpShldRT{4}(~isnan(selected_UnpShldRT{4}));
% 
% % data = vertcat(blockstruct.UnpShields);
% % idx = cellfun(@(x) ~isempty(x), data);
% % UnpShldRT = cellfun(@(x) x(:,2), data(idx),'UniformOutput',false);
% % UnpShldRT = vertcat(UnpShldRT{:});
% % UnpShldRT = UnpShldRT(~isnan(UnpShldRT));
% 
% %%
% p = histogram(selected_PunShldRT{4},linspace(0,3,30),'FaceColor',col_rep(2),'FaceAlpha',0.5); hold on
% u = histogram(selected_UnpShldRT{4},linspace(0,3,30),'FaceColor',col_rep(1),'FaceAlpha',0.5);
% legend([p u],{'PunShieldRT' 'UnpShieldRT'});
% xlim([0 3]);
% 
% figure;
% subplot(3,1,1); hold on
% p1 = histogram(selected_PunShldRT{1},linspace(0,3,30),'FaceColor',col_rep(2),'FaceAlpha',0.2);
% u1 = histogram(selected_UnpShldRT{1},linspace(0,3,30),'FaceColor',col_rep(1),'FaceAlpha',0.2);
% legend([p1 u1],...
%   {'PunShldRT B3' 'UnpShieldRT B3'});
% xlim([0 3]);
% 
% subplot(3,1,2); hold on
% p2 = histogram(selected_PunShldRT{2},linspace(0,3,30),'FaceColor',col_rep(2),'FaceAlpha',0.3); 
% u2 = histogram(selected_UnpShldRT{2},linspace(0,3,30),'FaceColor',col_rep(1),'FaceAlpha',0.3);
% legend([p2 u2],...
%   {'PunShldRT B5' 'UnpShieldRT B5'});
% xlim([0 3]);
% 
% subplot(3,1,3); hold on
% p3 = histogram(selected_PunShldRT{3},linspace(0,3,30),'FaceColor',col_rep(2),'FaceAlpha',0.4); 
% u3 = histogram(selected_UnpShldRT{3},linspace(0,3,30),'FaceColor',col_rep(1),'FaceAlpha',0.4);
% legend([p3 u3],...
%   {'PunShldRT B5' 'UnpShieldRT B5'});
% xlim([0 3]);
% 
% %% RT figure
% figure; hold on
% histogram(vertcat(HCP_aggr(:).ValInf_AvgRT),...
%   [0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 6 7 8 9 10 15 20 25 30 60 1205], 'Normalization', 'pdf');
