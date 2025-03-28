%% HumanCondPun extraction
% Extracts json files within "Raw" folder (work_folder\Raw\[test_names]\jatos###.txt) into HCP_aggr
% variable in workspace. Preload Raw folder with subfolder name and json files

pun_trials = 3:6;

end_trial_wait = 1000;
block_duration = 180*1000+end_trial_wait;

HCP_parameters.signal_time = 2000;
HCP_parameters.ship_duration = 6000;
HCP_parameters.ship_attack_time = 6000;
HCP_parameters.feedback_duration = 2500;
HCP_parameters.shield_charging_time = 3000;
HCP_parameters.shield_cost = -50; % point loss (i.e. -#)

HCP_parameters.BinLabels = {'ITI'; 'PunS_PreShield'; 'PunS_ShieldAv'; ...
        'PunS_Shielded'; 'PunS_NoShield'; 'PunS_ShieldedOutcome'; 'PunS_UnshieldedOutcome'; ...
        'UnpS_PreShield'; 'UnpS_ShieldAv'; 'UnpS_Shielded'; 'UnpS_NoShield'; ...
        'UnpS_ShieldedOutcome'; 'UnpS_UnshieldedOutcome'}; 
      % bins for HCP_windowbins.m 

test_names = {'Test1' 'Retest'}; % order of tests 
phase_names = {'phase1' 'Phase2'}; % {'phase1' 'Phase2'} 2nd entry = pun block
ques_names = {'ques_cfi' 'ques_htq' 'ques_audit'}; % questionnaire names
ques_checks = {[7;0] [] [7;3]}; % [index;correct answer] for catch questions
RT_check_grid = [3 3 6 6 6 6; 1 1 4 4 4 4; 1 1 4 4 4 4; 0 0 1 1 1 1; 0 0 1 1 1 1; ...
  1 1 1 1 1 1; 1 1 1 1 1 1]; % # post-block val/inf questions
RT_lowcutoff = .8; % valid RT lower cut-off in secs (any RT/#qs < low cutoff = invalid RT)
RT_highcutoff = 30; % valid RT upper cut-off in secs (any RT/#qs > high cutoff = invalid RT)

work_folder = ...
  ['G:\Current\Staff folders\Phil\Post-sub\Experiments\Collaborations\Human Conditioned Punishment\Prolific1 (Test-Retest)'];

ev_ylabel = {'PunShip' 'UnpShip' 'Shield' 'PunClick' 'UnpClick' 'ShieldClick' 'OtherClick' '+100' 'Loss'};
ev_ytick = [12 11 10 8 7 6 5 2 1];
ev_binning_tolerance = 100; %ms

%% House-keeping
clc
close all 
fclose('all');
start_tic_keep = tic;

date = char(datetime('today'));
n_tests = length(test_names);
n_bins = length(HCP_parameters.BinLabels);
[~,ytick_idx] = sort(ev_ytick);

% Make folders, save scripts/functions to Scripts folder
mkdir(work_folder, 'Figs');
mkdir([work_folder '\Figs'],'Click location');
mkdir([work_folder '\Figs'],'Summary');
mkdir([work_folder '\Figs'],'Timecourse');
mkdir([work_folder '\Scripts'],date);
copy_file = {'Phil_HCP2extract_retest.m' 'bin_rate2.m' 'HCP2_sumfig2.m' ...
  'HCP_windowbins.m' 'weighted_rate.m' 'trigger_idx.m' 'col_rep.m'};
for c = 1:length(copy_file)
   copyfile(copy_file{c},[work_folder '\Scripts\' date])
end

% Block structure fieldnames
cell_vars = {'NormClickLoc' 'PunClicks' 'UnpClicks' 'ShieldClicks' 'OtherClicks' 'PunShips' ...
  'UnpShips' 'PunShields' 'UnpShields' 'Attack' 'ShieldCost' 'PunRew' 'UnpRew' ...
  'BinPunRate' 'BinUnpRate' 'BinShieldRate' 'BinTime'};
d_vars = {'ExpPun_ShipRate' 'ExpUnp_ShipRate' ...
  'EndPoints' 'End' 'RewVal' 'PunPlanetVal' 'UnpPlanetVal' 'PunShipVal' 'UnpShipVal' 'AttackVal' ...
  'PunPlanet_RewInf' 'UnpPlanet_RewInf' ...
  'PunPlanet_PunShipInf' 'PunPlanet_UnpShipInf' ...
  'PunPlanet_AttackInf' 'UnpPlanet_PunShipInf' ...
  'UnpPlanet_UnpShipInf' 'UnpPlanet_AttackInf' ...
  'PunShip_AttackInf' 'UnpShip_AttackInf' ...
  'ValCheckRT' 'InfCheck1RT' 'InfCheck2RT' 'InfCheckSh1RT' 'InfCheckSh2RT'};

diaryfile = [work_folder '\Scripts\' date '\HCP extract (' date ').txt'];
diary(diaryfile);
 
%% Searches work_folder for new data - load 
cur_row = 0;

for t = 1:n_tests
  fprintf([' - Extract ' test_names{t} ' -\n']);
  %Set block list to blocks not within HCP_meta
  datasets = dir([work_folder '\Raw\' test_names{t}]);
  datasets = datasets(endsWith({datasets.name},'.txt'));
  dataset_new = {datasets.name}';
  dataset_new = cellfun(@(x) x(1:end-4), dataset_new, 'UniformOutput', false);
  
  if ~isempty(dataset_new)
    fprintf('Add: \n');
    disp(dataset_new)
  
    %Preallocate extra cells (doesn't affect prev aggr cells)
    n_datasets = length(dataset_new);
  
    for d = 1:n_datasets
      data_name = dataset_new{d};
      fprintf(['Extracting ' data_name '. ']);
  
      json_text = importdata([work_folder '\Raw\' test_names{t} '\' data_name '.txt']);
  
      n_j = length(json_text);
      fprintf([num2str(n_j) ' set(s) found\n']);
  
      if n_j > 0
        if t == 1
          HCP_aggr(cur_row+n_j).(['raw_' test_names{t}]){1,3} = [];
        else
          test1_IDs = {HCP_aggr(:).(['SubjID_' test_names{1}])};
        end

        for s = 1:n_j
          try
            data = jsondecode(json_text{s});
          catch
            if ~matches(json_text{s}(end),']')
              data = jsondecode([json_text{s} ']']);
              warning(['json data (#' int2str(s) ') missing "]" on end. "]" added to data.'])
            else
              error('jsondecode error');
            end
          end

          %% Find HCP_aggr row
          if t == 1
            cur_row = cur_row+1;
          else
            retest_SubjID = data{1}.subject_id;

            if ~isempty(retest_SubjID) 
              match_idx = find(ismember(test1_IDs,retest_SubjID));

              if ~isempty(match_idx)
                if length(match_idx) > 1
                  warning('TAKE NOTE!!!');
                  warning(['Multiple matching subjIDs for ' retest_SubjID ' in ' test_names{1} ...
                    ' (HCP_aggr rows:' strjoin(cellfun(@(x) int2str(x),num2cell(match_idx), ...
                    'UniformOutput',false),',') '). ' data_name ',set' int2str(s) ' data added to both rows']);
                else
                  fprintf(['subjIDs ' retest_SubjID ' added to HCP_aggr row' int2str(match_idx) '\n']);
                end
                cur_row = match_idx;

              else % no matching subjIDs
                cur_row = length(HCP_aggr)+1;
                warning(['No ' test_names{1} ' subjID (' retest_SubjID '), json data (#' int2str(s) ...
                  ')!! ' retest_SubjID ' ' test_names{s} ' data added to HCP_aggr row' int2str(cur_row)]);
                HCP_aggr(cur_row).(['SubjID_' test_names{1}]) = ...
                  [test_names{t} '_NoMatch_' int2str(cur_row)];
              end
            else % no subjID for re-test
              cur_row = length(HCP_aggr)+1;
              warning(['No subjID for json data (#' int2str(s) ...
                ')!! Data added to HCP_aggr row' int2str(cur_row)]);
              HCP_aggr(cur_row).(['SubjID_' test_names{1}]) = ...
                ['MISSING_' test_names{t} '_subjID_' int2str(cur_row)];
            end % re-test subjID check
          end

          data_leng = length(data);
          for c = 1:length(cur_row)
            HCP_aggr(cur_row(c)).(['raw_' test_names{t}]){data_leng,3} = [];
            HCP_aggr(cur_row(c)).(['raw_' test_names{t}])(:,3) = data;
            
            HCP_aggr(cur_row(c)).(['json_' test_names{t}]) = data_name;
            HCP_aggr(cur_row(c)).(['Set_' test_names{t}]) = s;
            HCP_aggr(cur_row(c)).(['Sampling_' test_names{t}]) = HCP_aggr(cur_row(c)).(['raw_' test_names{t}]){1,3}.sample;
            if ~isempty(HCP_aggr(cur_row(c)).(['raw_' test_names{t}]){1,3}.subject_id)
              HCP_aggr(cur_row(c)).(['SubjID_' test_names{t}]) = ...
                HCP_aggr(cur_row(c)).(['raw_' test_names{t}]){1,3}.subject_id;
            else
              HCP_aggr(cur_row(c)).(['SubjID_' test_names{t}]) = ...
                ['MISSING_' test_names{t} '_subjID_' int2str(cur_row(c))];
            end
    
            for p = 1:data_leng
               if isfield(data{p},'phase')
                  HCP_aggr(cur_row(c)).(['raw_' test_names{t}]){p,1} = data{p}.phase;
               end
               if isfield(data{p},'block_number')
                  HCP_aggr(cur_row(c)).(['raw_' test_names{t}]){p,2} = data{p}.block_number;
               end 
            end
          end % c loop
          
          fprintf(['Extracted #' int2str(s) '. '])
          toc(start_tic_keep);
        end % s (json set) loop
      else
        warning('No json files data found - skipped.');
      end
    end % d (json file) loop
  else
    fprintf('None new.\n');
  end % check new datasets
  
  fprintf('\n');
end % test loop

%% Process data into HCP_aggr structure
fprintf('Processing data...\n');
for t = 1:n_tests
  fprintf(['\n -- ' test_names{t} '--\n']);
for r = 1:length(HCP_aggr)
  fprintf(['\n - aggr' int2str(r) '-\n']);
     
  if ~isempty(HCP_aggr(r).(['raw_' test_names{t}]))
  HCP_aggr(r).(['Group_' test_names{t}]) = HCP_aggr(r).(['raw_' test_names{t}]){1,3}.group;
  
  % Demographics
  is_idx = find(cellfun(@(x) ~isempty(x),HCP_aggr(r).(['raw_' test_names{t}])(:,1)));
  
  p = is_idx(ismember(HCP_aggr(r).(['raw_' test_names{t}])(is_idx,1),'demographics'));
  [~,g_eIdx] = regexp(HCP_aggr(r).(['raw_' test_names{t}]){p,3}.responses,'"gender":"');
  [a_sIdx,a_eIdx] = regexp(HCP_aggr(r).(['raw_' test_names{t}]){p,3}.responses,'","age":"');
  [l_sIdx,l_eIdx] = regexp(HCP_aggr(r).(['raw_' test_names{t}]){p,3}.responses,'","language":"');
  
  HCP_aggr(r).(['Gender_' test_names{t}]) = ...
    HCP_aggr(r).(['raw_' test_names{t}]){p,3}.responses(g_eIdx+1:a_sIdx-1);
  HCP_aggr(r).(['Age_' test_names{t}]) = ...
    str2num(HCP_aggr(r).(['raw_' test_names{t}]){p,3}.responses(a_eIdx+1:l_sIdx-1));
  if isempty(HCP_aggr(r).(['Age_' test_names{t}]))
    HCP_aggr(r).(['Age_' test_names{t}]) = NaN;
  end
  HCP_aggr(r).(['Language_' test_names{t}]) = ...
    HCP_aggr(r).(['raw_' test_names{t}]){p,3}.responses(l_eIdx+1:end-2);
  

  fprintf([HCP_aggr(r).(['Gender_' test_names{t}]) ',' ...
    num2str(HCP_aggr(r).(['Age_' test_names{t}])) ',' ...
    HCP_aggr(r).(['Language_' test_names{t}]) '\nGroup: ' HCP_aggr(r).(['Group_' test_names{t}]) '\n']);
  
  %% Prep variables
  phase_idx = is_idx(ismember(HCP_aggr(r).(['raw_' test_names{t}])(is_idx,1),phase_names));
  n_phases = length(phase_idx);
  
  punP = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(1),3}.pun_planet_side;
  if punP == '1'
     unpP = '0';
     punP_n = 1;
     unpP_n = 0;
     HCP_aggr(r).(['PunPlanet_' test_names{t}]) = 'right';
     fprintf('PunPlanet = right |');
  else
     unpP = '1';
     punP_n = 0;
     unpP_n = 1;
     HCP_aggr(r).(['PunPlanet_' test_names{t}]) = 'left';
     fprintf('PunPlanet = left | ');
  end
  [~,idx] = regexp(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(end),3}.pun_ship,'ship');
  punS = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(end),3}.pun_ship(idx+1);
  if punS == '1'
    unpS = '2';
     HCP_aggr(r).(['PunShip_' test_names{t}]) = 'TypeI';
    fprintf('PunShip = Type I\n');
  elseif punS == '2'
    unpS = '1';
     HCP_aggr(r).(['PunShip_' test_names{t}]) = 'TypeII';
    fprintf('PunShip = Type II\n');
  else
    error('Ship types?');
  end

  click_loc{n_phases} = [];
  pun_clicks{n_phases} = [];
  unp_clicks{n_phases} = [];
  shield_clicks{n_phases} = [];
  other_clicks{n_phases} = [];
  pun_ships{n_phases} = [];
  unp_ships{n_phases} = [];
  
  % Prep saved variable structure
  for v = 1:length(cell_vars)
    HCP_aggr(r).(['Block_' test_names{t}]).(cell_vars{v}){n_phases} = [];
  end
  for v = 1:length(d_vars)
    HCP_aggr(r).(['Block_' test_names{t}]).(d_vars{v}) = NaN(1,n_phases);
  end
  
  click_fig = figure;
  normclick_fig = figure;
  t_fig = figure; hold on
  ylim([min(ev_ytick)*.5 max(ev_ytick)+min(ev_ytick)*.5]);
  yticks(ev_ytick(ytick_idx));
  yticklabels(ev_ylabel(ytick_idx));
  xlabel('Seconds');
  prev_times = 0;

  %% Loop through phases
  for p = 1:n_phases
    fprintf(['Block' num2str(p)]);
    
    click_loc{p} = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_clicks.loc;
    if ~isempty(click_loc{p})
      view_size = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.viewport_size./2;
      if HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.viewport_size(1) > ...
          HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.viewport_size(2)
        HCP_aggr(r).(['Block_' test_names{t}]).NormClickLoc{p} = click_loc{p};
          HCP_aggr(r).(['Block_' test_names{t}]).NormClickLoc{p}(:,1) = ...
            (HCP_aggr(r).(['Block_' test_names{t}]).NormClickLoc{p}(:,1)-view_size(1))./view_size(1);
          HCP_aggr(r).(['Block_' test_names{t}]).NormClickLoc{p}(:,2) = ...
            -(HCP_aggr(r).(['Block_' test_names{t}]).NormClickLoc{p}(:,2)-view_size(2))./(view_size(2));
      else
        HCP_aggr(r).(['Block_' test_names{t}]).NormClickLoc{p} = fliplr(click_loc{p});
        HCP_aggr(r).(['Block_' test_names{t}]).NormClickLoc{p}(:,1) = ...
          (HCP_aggr(r).(['Block_' test_names{t}]).NormClickLoc{p}(:,1)-view_size(2))./view_size(2);
        HCP_aggr(r).(['Block_' test_names{t}]).NormClickLoc{p}(:,2) = ...
          -(HCP_aggr(r).(['Block_' test_names{t}]).NormClickLoc{p}(:,2)-view_size(1))./(view_size(1));
      end
    end % skip if no clicks
    
    clicks = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_clicks.element;
    if ~isempty(clicks)
      clicks(cellfun(@(x) isempty(x),HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_clicks.element)) = {'none'};
    end
    
    %% Clicks (timestamps [per type])
    pun_clicks{p} = ismember(clicks,['planet-' punP]);
    HCP_aggr(r).(['Block_' test_names{t}]).PunClicks{p} = ...
      HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_clicks.timestamp(pun_clicks{p})/1000;
    unp_clicks{p} = ismember(clicks,['planet-' unpP]);
    HCP_aggr(r).(['Block_' test_names{t}]).UnpClicks{p} = ...
      HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_clicks.timestamp(unp_clicks{p})/1000;    
    shield_clicks{p} = ismember(clicks,'ship-shield-button');
    HCP_aggr(r).(['Block_' test_names{t}]).ShieldClicks{p} = ...
      HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_clicks.timestamp(shield_clicks{p})/1000;
    other_clicks{p} = ismember(clicks,'none');
    HCP_aggr(r).(['Block_' test_names{t}]).OtherClicks{p} = ...
      HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_clicks.timestamp(other_clicks{p})/1000;
    
    if ismember(p,pun_trials)
      %% Ships [pun v unp] (col1 = start, col2 = outcome (disappear = col2 + feedback_duration))
      pun_ships{p} = ismember(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.type,punP);
      HCP_aggr(r).(['Block_' test_names{t}]).PunShips{p} = zeros(sum(pun_ships{p}),2);
      HCP_aggr(r).(['Block_' test_names{t}]).PunShips{p}(:,1) = ...
        HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.time_appear(pun_ships{p})/1000;
      HCP_aggr(r).(['Block_' test_names{t}]).PunShips{p}(:,2) = ...
        (HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.time_outcome(pun_ships{p}))/1000;
      
      unp_ships{p} = ismember(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.type,unpP);
      HCP_aggr(r).(['Block_' test_names{t}]).UnpShips{p} = zeros(sum(unp_ships{p}),2);
      HCP_aggr(r).(['Block_' test_names{t}]).UnpShips{p}(:,1) = ...
        HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.time_appear(unp_ships{p})/1000;
      HCP_aggr(r).(['Block_' test_names{t}]).UnpShips{p}(:,2) = ...
        (HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.time_outcome(unp_ships{p}))/1000;

      %% Shields (col1 = start, col2 = rt activated, col3 = outcome (disappear = col3 + feedback_duration))
      idx = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.shield_available(pun_ships{p});
      HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p} = NaN(sum(pun_ships{p}),3);
      HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p}(idx,1) = ...
        HCP_aggr(r).(['Block_' test_names{t}]).PunShips{p}(idx,1)+HCP_parameters.shield_charging_time/1000;
      HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p}(idx,2) = ...
        (HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.rt_shield_activated(idx))/1000;
      HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p}(idx,3) = ...
        HCP_aggr(r).(['Block_' test_names{t}]).PunShips{p}(idx,2);

      idx = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.shield_available(unp_ships{p});
      HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p} = NaN(sum(unp_ships{p}),3);
      HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p}(idx,1) = ...
        HCP_aggr(r).(['Block_' test_names{t}]).UnpShips{p}(idx,1)+HCP_parameters.shield_charging_time/1000;
      HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p}(idx,2) = ...
        (HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.rt_shield_activated(idx))/1000;
      HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p}(idx,3) = ...
        HCP_aggr(r).(['Block_' test_names{t}]).UnpShips{p}(idx,2);
    end
    
    %% Point Outcomes (col1 = amount, col2 = time [per cause])
    HCP_aggr(r).(['Block_' test_names{t}]).EndPoints(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.points_total;
    
    loss_idx = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_outcomes.outcome < 0;  
    if any(loss_idx)
      % Attack
      idx = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.outcome < 0;
      HCP_aggr(r).(['Block_' test_names{t}]).Attack{p} = zeros(sum(idx),2);
      HCP_aggr(r).(['Block_' test_names{t}]).Attack{p}(:,1) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.outcome(idx);
      HCP_aggr(r).(['Block_' test_names{t}]).Attack{p}(:,2) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.ships.time_outcome(idx);
      
      % Shield cost
      if HCP_parameters.shield_cost < 0
        points = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_outcomes.outcome(loss_idx);
        times = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_outcomes.time_outcome(loss_idx)/1000;
        idx = points == HCP_parameters.shield_cost;
        HCP_aggr(r).(['Block_' test_names{t}]).ShieldCost{p} = zeros(sum(idx),2);
        HCP_aggr(r).(['Block_' test_names{t}]).ShieldCost{p}(:,1) = points(idx);
        HCP_aggr(r).(['Block_' test_names{t}]).ShieldCost{p}(:,2) = times(idx);
      else
        HCP_aggr(r).(['Block_' test_names{t}]).ShieldCost{p} = [];
      end     
    else
      HCP_aggr(r).(['Block_' test_names{t}]).ShieldCost{p} = [];
      HCP_aggr(r).(['Block_' test_names{t}]).Attack{p} = [];
    end
    
    points = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.planets.outcome;
    times = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.planets.time_outcome/1000;
    planet_idx = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.planets.select;

    punP_idx = planet_idx==punP_n;
    unpP_idx = planet_idx==unpP_n;

    HCP_aggr(r).(['Block_' test_names{t}]).PunRew{p} = zeros(sum(punP_idx),2);
    HCP_aggr(r).(['Block_' test_names{t}]).PunRew{p}(:,1) = points(punP_idx);
    HCP_aggr(r).(['Block_' test_names{t}]).PunRew{p}(:,2) = times(punP_idx);
    HCP_aggr(r).(['Block_' test_names{t}]).UnpRew{p} = zeros(sum(unpP_idx),2);
    HCP_aggr(r).(['Block_' test_names{t}]).UnpRew{p}(:,1) = points(unpP_idx);
    HCP_aggr(r).(['Block_' test_names{t}]).UnpRew{p}(:,2) = times(unpP_idx);
    
    HCP_aggr(r).(['Block_' test_names{t}]).End(p) = max([block_duration ...
      max(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_outcomes.time_outcome)+HCP_parameters.feedback_duration ...
      max(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_clicks.timestamp)])/1000;
    
    %% Value & inference checks
    HCP_aggr(r).(['Block_' test_names{t}]).ValCheckRT(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.rt/1000;
    HCP_aggr(r).(['Block_' test_names{t}]).InfCheck1RT(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rt/1000;
    HCP_aggr(r).(['Block_' test_names{t}]).InfCheck2RT(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rt/1000;
    
    HCP_aggr(r).(['Block_' test_names{t}]).RewVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_1);

    if ismember(p,pun_trials) % if punished phase
      HCP_aggr(r).(['Block_' test_names{t}]).InfCheckSh1RT(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+4,3}.rt/1000;
      HCP_aggr(r).(['Block_' test_names{t}]).InfCheckSh2RT(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+5,3}.rt/1000;
      
      HCP_aggr(r).(['Block_' test_names{t}]).PredCheckRT(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+6,3}.rt/1000;
      HCP_aggr(r).(['Block_' test_names{t}]).RevalCheckRT(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+7,3}.rt/1000;
     
      % Planet values (Pun phase)
      if punP == '0' %change val 2,3 -> val 3,4
        HCP_aggr(r).(['Block_' test_names{t}]).PunPlanetVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_3);
        HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanetVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_4);
        
        HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_RewInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_1);
        %HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_RewInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_1);
        HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_RewInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_1);   
        %HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_RewInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_1);  
      elseif punP == '1' %change val 3,2 -> val 4,3
        HCP_aggr(r).(['Block_' test_names{t}]).PunPlanetVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_4);
        HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanetVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_3);
        
        HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_RewInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_1);
        % HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_RewInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_1);
        HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_RewInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_1);   
        % HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_RewInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_1);  
      else
        error('Huh planet?')
      end

      % Attack value (Pun phase) [03/02: changed val 6 to val 2]
      HCP_aggr(r).(['Block_' test_names{t}]).AttackVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_2);
      
      % Ship value (Pun phase)
      if punS == '1'%[03/02: changed val 4 to val 5, val 5 to val 6]
        HCP_aggr(r).(['Block_' test_names{t}]).PunShipVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_5);
        HCP_aggr(r).(['Block_' test_names{t}]).UnpShipVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_6);
        
        HCP_aggr(r).(['Block_' test_names{t}]).PunShip_AttackInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+4,3}.rate_1);
        % HCP_aggr(r).(['Block_' test_names{t}]).PunShip_AttackInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+4,3}.conf_1);
        HCP_aggr(r).(['Block_' test_names{t}]).UnpShip_AttackInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+5,3}.rate_1);
        % HCP_aggr(r).(['Block_' test_names{t}]).UnpShip_AttackInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+5,3}.conf_1);
        
        if punP == '0'
          % R1:R2 response prediction/revaluation
          HCP_aggr(r).(['Block_' test_names{t}]).PrefPred(p) = -HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+6,3}.val+100;    
          HCP_aggr(r).(['Block_' test_names{t}]).PrefReval(p) = -HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+7,3}.val+100;    

          HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_PunShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_3);
          %HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_2);
          HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_4);
          %HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_3);
          HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_AttackInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_2);
          %HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_4);
          
          HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_PunShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_3);
          %HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_2);
          HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_4);
          %HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_3);
          HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_AttackInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_2);
          %HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_4);
        elseif punP == '1'
          % R1:R2 response prediction/revaluation
          HCP_aggr(r).(['Block_' test_names{t}]).PrefPred(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+6,3}.val;   
          HCP_aggr(r).(['Block_' test_names{t}]).PrefReval(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+7,3}.val;   

          HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_PunShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_3);
          % HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_2);
          HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_4);
          % HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_3);
          HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_AttackInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_2);
          % HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_4);
          
          HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_PunShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_3);
          % HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_2);
          HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_4);
          % HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_3);
          HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_AttackInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_2);
          % HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_4);
        else
          error('Huh planet?')
        end
        
      elseif punS == '2'%change val 5,4 -> val 6,5
        HCP_aggr(r).(['Block_' test_names{t}]).PunShipVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_6);
        HCP_aggr(r).(['Block_' test_names{t}]).UnpShipVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_5);
        
        HCP_aggr(r).(['Block_' test_names{t}]).PunShip_AttackInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+5,3}.rate_1);
        %HCP_aggr(r).(['Block_' test_names{t}]).PunShip_AttackInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+5,3}.conf_1);
        HCP_aggr(r).(['Block_' test_names{t}]).UnpShip_AttackInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+4,3}.rate_1);
        %HCP_aggr(r).(['Block_' test_names{t}]).UnpShip_AttackInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+4,3}.conf_1);
      
        if punP == '0'
          HCP_aggr(r).(['Block_' test_names{t}]).PrefPred(p) = -HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+6,3}.val+100;    
          HCP_aggr(r).(['Block_' test_names{t}]).PrefReval(p) = -HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+7,3}.val+100;    

          HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_PunShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_4);
         %HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_3);
          HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_3);
         %HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_2);
          HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_AttackInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_2);
         %HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_4);
          
          HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_PunShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_4);
         %HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_3);
          HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_3);
         %HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_2);
          HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_AttackInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_2);
         %HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_4);
         
        elseif punP == '1'        
          HCP_aggr(r).(['Block_' test_names{t}]).PrefPred(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+6,3}.val;   
          HCP_aggr(r).(['Block_' test_names{t}]).PrefReval(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+7,3}.val;   

          HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_PunShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_4);
         %HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_3);
          HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_3);
         %HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_2);
          HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_AttackInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_2);
         %HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_4);
          
          HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_PunShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_4);
         %HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_3);
          HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_3);
         %HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_2);
          HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_AttackInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_2);
         %HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_4);
        else
          error('Huh planet?')
        end
      else
        error('Huh ship?')
      end
      
      %% Click rate binning
      % Ship/shield bins (all else = ITI) 
      % {1-6} PunShip: PreShield, ShieldAvail, Shielded, NoShield, ShieldedOutcome, UnshieldedOutcome     
      % {7-12} UnpShip: PreShield, ShieldAvail, Shielded, NoShield, ShieldedOutcome, UnshieldedOutcome 
      windows = HCP_windowbins(HCP_aggr(r).(['Block_' test_names{t}]).PunShips{p},...
        HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p},...
        HCP_aggr(r).(['Block_' test_names{t}]).UnpShips{p},...
        HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p},...
        HCP_parameters.shield_charging_time/1000,HCP_parameters.feedback_duration/1000);
           
      % Set BinTime defaults
      ITI_time = HCP_aggr(r).(['Block_' test_names{t}]).End(p);
      times = zeros(length(windows),1);
      
      % Bin PunPlanet clicks
      if ~isempty(HCP_aggr(r).(['Block_' test_names{t}]).PunClicks{p})
        [ITI_rate, ITI_time, rates, times] = bin_rate2(HCP_aggr(r).(['Block_' test_names{t}]).PunClicks{p},...
          HCP_aggr(r).(['Block_' test_names{t}]).End(p), windows);
        HCP_aggr(r).(['Block_' test_names{t}]).BinPunRate{p} = [ITI_rate; vertcat(rates{:})];
      else
        HCP_aggr(r).(['Block_' test_names{t}]).BinPunRate{p} = zeros(n_bins,1);
      end
      % Bin UnpPlanet clicks
      if ~isempty(HCP_aggr(r).(['Block_' test_names{t}]).UnpClicks{p})
        [ITI_rate, ITI_time, rates, times] = bin_rate2(HCP_aggr(r).(['Block_' test_names{t}]).UnpClicks{p},...
          HCP_aggr(r).(['Block_' test_names{t}]).End(p), windows);
        HCP_aggr(r).(['Block_' test_names{t}]).BinUnpRate{p} = [ITI_rate; vertcat(rates{:})];
      else
        HCP_aggr(r).(['Block_' test_names{t}]).BinUnpRate{p} = zeros(n_bins,1);
      end      
      % Bin Shield button clicks
      if ~isempty(HCP_aggr(r).(['Block_' test_names{t}]).ShieldClicks{p})
        [ITI_rate, ITI_time, rates, times] = bin_rate2(HCP_aggr(r).(['Block_' test_names{t}]).ShieldClicks{p},...
          HCP_aggr(r).(['Block_' test_names{t}]).End(p), windows);
        HCP_aggr(r).(['Block_' test_names{t}]).BinShieldRate{p} = [ITI_rate; vertcat(rates{:})];
      else
        HCP_aggr(r).(['Block_' test_names{t}]).BinShieldRate{p} = zeros(n_bins,1);
      end      
      
      HCP_aggr(r).(['Block_' test_names{t}]).BinTime{p} = [ITI_time; times];
      
      %% Experienced Response-Ship contingency (Ships per ITI press)
      HCP_aggr(r).(['Block_' test_names{t}]).ExpPun_ShipRate(p) = size(HCP_aggr(r).(['Block_' test_names{t}]).PunShips{p},1)/...
        (HCP_aggr(r).(['Block_' test_names{t}]).BinPunRate{p}(1)*HCP_aggr(r).(['Block_' test_names{t}]).BinTime{p}(1));
      HCP_aggr(r).(['Block_' test_names{t}]).ExpUnp_ShipRate(p) = size(HCP_aggr(r).(['Block_' test_names{t}]).UnpShips{p},1)/...
        (HCP_aggr(r).(['Block_' test_names{t}]).BinUnpRate{p}(1)*HCP_aggr(r).(['Block_' test_names{t}]).BinTime{p}(1));
%       fprintf([' - experienced Pun/UnpShip probability: ' num2str(round(HCP_aggr(r).(['Block_' test_names{t}]).ExpPun_ShipRate(p),3)) ...
%         '/' num2str(round(HCP_aggr(r).(['Block_' test_names{t}]).ExpUnp_ShipRate(p),3))]);
      
    else % not pun trial

      % Planet values (Non-pun phase)
      if punP == '0'
        HCP_aggr(r).(['Block_' test_names{t}]).PunPlanetVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_2);
        HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanetVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_3);
        
        HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_RewInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_1);
        %HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_RewInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_1);
        HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_RewInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_1);   
        %HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_RewInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_1);  
      elseif punP == '1' 
        HCP_aggr(r).(['Block_' test_names{t}]).PunPlanetVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_3);
        HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanetVal(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+1,3}.val_2);
        
        HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_RewInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.rate_1);
        % HCP_aggr(r).(['Block_' test_names{t}]).PunPlanet_RewInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+3,3}.conf_1);
        HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_RewInf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.rate_1);   
        % HCP_aggr(r).(['Block_' test_names{t}]).UnpPlanet_RewInfconf(p) = str2double(HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+2,3}.conf_1);  
      else
        error('Huh planet?')
      end

      % R1:R2 response prediction/revaluation
      HCP_aggr(r).(['Block_' test_names{t}]).PredCheckRT(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+4,3}.rt/1000;
      HCP_aggr(r).(['Block_' test_names{t}]).RevalCheckRT(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+5,3}.rt/1000;

      if punP == '0'
        HCP_aggr(r).(['Block_' test_names{t}]).PrefPred(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+4,3}.val;
        HCP_aggr(r).(['Block_' test_names{t}]).PrefReval(p) = HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+5,3}.val;
      elseif punP == '1'
        HCP_aggr(r).(['Block_' test_names{t}]).PrefPred(p) = -HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+4,3}.val+100;
        HCP_aggr(r).(['Block_' test_names{t}]).PrefReval(p) = -HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p)+5,3}.val+100;
      end

      % PunRates
      HCP_aggr(r).(['Block_' test_names{t}]).BinPunRate{p} = NaN(n_bins,1);
      HCP_aggr(r).(['Block_' test_names{t}]).BinUnpRate{p} = NaN(n_bins,1);
      
      if ~isempty(HCP_aggr(r).(['Block_' test_names{t}]).PunClicks{p})
        HCP_aggr(r).(['Block_' test_names{t}]).BinPunRate{p}(1) = ...
          size(HCP_aggr(r).(['Block_' test_names{t}]).PunClicks{p},1)/HCP_aggr(r).(['Block_' test_names{t}]).End(p);
      else
        HCP_aggr(r).(['Block_' test_names{t}]).BinPunRate{p}(1) = 0;
      end
      % UnpRates (only ITI)
      if ~isempty(HCP_aggr(r).(['Block_' test_names{t}]).UnpClicks{p})
        HCP_aggr(r).(['Block_' test_names{t}]).BinUnpRate{p}(1) = ...
          size(HCP_aggr(r).(['Block_' test_names{t}]).UnpClicks{p},1)/HCP_aggr(r).(['Block_' test_names{t}]).End(p);
      else
        HCP_aggr(r).(['Block_' test_names{t}]).BinUnpRate{p}(1) = 0;
      end
      
      HCP_aggr(r).(['Block_' test_names{t}]).BinTime{p} = [HCP_aggr(r).(['Block_' test_names{t}]).End(p); zeros(12,1)];
    end
    
    %% Plotting
    % Click location figure
    if ~isempty(click_loc{p})
      figure(click_fig); hold on
      scatter(click_loc{p}(:,1),-click_loc{p}(:,2),[],col_rep(p)); 

      figure(normclick_fig); hold on
      scatter(HCP_aggr(r).(['Block_' test_names{t}]).NormClickLoc{p}(:,1),HCP_aggr(r).(['Block_' test_names{t}]).NormClickLoc{p}(:,2),[],col_rep(p)); 
    end
    
    %% Plot intra-phase data
    figure(t_fig);
    subplot(2,1,1); hold on
    
    %Clicks 
    plot((HCP_aggr(r).(['Block_' test_names{t}]).PunClicks{p}+prev_times),ones(sum(pun_clicks{p}),1)*ev_ytick(4),...
      '.','Color',col_rep(2),'LineStyle','none');
    plot((HCP_aggr(r).(['Block_' test_names{t}]).UnpClicks{p}+prev_times),ones(sum(unp_clicks{p}),1)*ev_ytick(5),...
      '.','Color',col_rep(3),'LineStyle','none');
    plot((HCP_aggr(r).(['Block_' test_names{t}]).ShieldClicks{p}+prev_times),ones(sum(shield_clicks{p}),1)*ev_ytick(6),...
      '.','Color',col_rep(1),'LineStyle','none');
    plot((HCP_aggr(r).(['Block_' test_names{t}]).OtherClicks{p}+prev_times),ones(sum(other_clicks{p}),1)*ev_ytick(7),...
      'k.','LineStyle','none');
    
    %Rewards
    idx = HCP_aggr(r).(['Block_' test_names{t}]).PunRew{p}(:,1)==100;
    plot((HCP_aggr(r).(['Block_' test_names{t}]).PunRew{p}(idx,2)+prev_times),ones(sum(idx),1)*ev_ytick(8)-.1,...
      '+','Color',col_rep(2),'LineStyle','none');
    idx = HCP_aggr(r).(['Block_' test_names{t}]).UnpRew{p}(:,1)==100;
    plot((HCP_aggr(r).(['Block_' test_names{t}]).UnpRew{p}(idx,2)+prev_times),ones(sum(idx),1)*ev_ytick(8)+.1,...
      '+','Color',col_rep(3),'LineStyle','none');
    if ~isempty(HCP_aggr(r).(['Block_' test_names{t}]).Attack{p})
      plot((HCP_aggr(r).(['Block_' test_names{t}]).Attack{p}(:,2)+prev_times),...
        ones(size(HCP_aggr(r).(['Block_' test_names{t}]).Attack{p},1),1)*ev_ytick(9),...
        'x','Color',col_rep(2),'LineStyle','none');
    end
    if ~isempty(HCP_aggr(r).(['Block_' test_names{t}]).ShieldCost{p})
    plot((HCP_aggr(r).(['Block_' test_names{t}]).ShieldCost{p}(:,2)+prev_times),...
      ones(size(HCP_aggr(r).(['Block_' test_names{t}]).ShieldCost{p},1),1)*ev_ytick(9),...
      'k+','LineStyle','none');
    end
    
    if ismember(p,pun_trials)
    %Ships & shields
    for s = 1:size(HCP_aggr(r).(['Block_' test_names{t}]).PunShips{p},1)
      plot([(HCP_aggr(r).(['Block_' test_names{t}]).PunShips{p}(s,1)+prev_times) ...
        (HCP_aggr(r).(['Block_' test_names{t}]).PunShips{p}(s,2)+prev_times)], ...
      	[ev_ytick(1) ev_ytick(1)],'Color',col_rep(2),'LineWidth',4);
      plot([(HCP_aggr(r).(['Block_' test_names{t}]).PunShips{p}(s,2)+prev_times) ...
        (HCP_aggr(r).(['Block_' test_names{t}]).PunShips{p}(s,2)+HCP_parameters.feedback_duration/1000 ...
        +prev_times)], [ev_ytick(1) ev_ytick(1)],'Color',col_rep(14),'LineWidth',4);        
    end 
    for s = 1:size(HCP_aggr(r).(['Block_' test_names{t}]).UnpShips{p},1)
      plot([(HCP_aggr(r).(['Block_' test_names{t}]).UnpShips{p}(s,1)+prev_times) ...
        (HCP_aggr(r).(['Block_' test_names{t}]).UnpShips{p}(s,2)+prev_times)], ...
      	[ev_ytick(2) ev_ytick(2)],'Color',col_rep(3),'LineWidth',4);
      plot([(HCP_aggr(r).(['Block_' test_names{t}]).UnpShips{p}(s,2)+prev_times) ...
        (HCP_aggr(r).(['Block_' test_names{t}]).UnpShips{p}(s,2)+HCP_parameters.feedback_duration/1000 ...
        +prev_times)], [ev_ytick(2) ev_ytick(2)],'Color',col_rep(14),'LineWidth',4);        
    end 
    for s = 1:size(HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p},1)
      if ~isnan(HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p}(s,1))
        plot([(HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p}(s,1)+prev_times) ...
          (HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p}(s,3)+prev_times)], ...
          [ev_ytick(3) ev_ytick(3)],'Color',col_rep(1),'LineWidth',4);
        if ~isnan(HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p}(s,2))
          plot([(HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p}(s,1)+...
            HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p}(s,2)+prev_times) ...
            (HCP_aggr(r).(['Block_' test_names{t}]).PunShields{p}(s,3)+prev_times)], ...
            [ev_ytick(3) ev_ytick(3)],'Color',col_rep(3),'LineWidth',4);
        end
      end
    end  
    for s = 1:size(HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p},1)
      if ~isnan(HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p}(s,1))
        plot([(HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p}(s,1)+prev_times) ...
          (HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p}(s,3)+prev_times)], ...
          [ev_ytick(3) ev_ytick(3)],'Color',col_rep(1),'LineWidth',4);
        if ~isnan(HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p}(s,2))
          plot([(HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p}(s,1)+HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p}(s,2)+prev_times) ...
            (HCP_aggr(r).(['Block_' test_names{t}]).UnpShields{p}(s,3)+prev_times)], ...
            [ev_ytick(3) ev_ytick(3)],'Color',col_rep(3),'LineWidth',4);
        end
      end
    end 
    end
    
    % Points
    subplot(2,1,2); hold on
    plot((HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_outcomes.time_outcome/1000+prev_times),...
      HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(p),3}.all_outcomes.total,'d-','Color',col_rep(1));   

    prev_times = prev_times + HCP_aggr(r).(['Block_' test_names{t}]).End(p) + 1;
    
    fprintf('\n');
  end % phase loop
  
  %% Report experience instrumental contingency
  fprintf(['Experienced ITI Pun/Unp -> Ship probability: ' ...
    num2str(round(sum(cellfun(@(x) size(x,1),HCP_aggr(r).(['Block_' test_names{t}]).PunShips(pun_trials)))/sum(cellfun(@(x) x(1),...
    HCP_aggr(r).(['Block_' test_names{t}]).BinPunRate(pun_trials)).*cellfun(@(x) x(1),HCP_aggr(r).(['Block_' test_names{t}]).BinTime(pun_trials))),3)) '/'...
    num2str(round(sum(cellfun(@(x) size(x,1),HCP_aggr(r).(['Block_' test_names{t}]).UnpShips(pun_trials)))/sum(cellfun(@(x) x(1),...
    HCP_aggr(r).(['Block_' test_names{t}]).BinUnpRate(pun_trials)).*cellfun(@(x) x(1),HCP_aggr(r).(['Block_' test_names{t}]).BinTime(pun_trials))),3)) '\n']);
  
  fprintf(['Experienced Pun/Unp (all "trades") -> Ship probability: ' ...
    num2str(round(sum(cellfun(@(x) size(x,1),HCP_aggr(r).(['Block_' test_names{t}]).PunShips(pun_trials)))/sum(cellfun(@(x) size(x,1),...
    HCP_aggr(r).(['Block_' test_names{t}]).PunRew(pun_trials))),3)) '/'...
    num2str(round(sum(cellfun(@(x) size(x,1),HCP_aggr(r).(['Block_' test_names{t}]).UnpShips(pun_trials)))/sum(cellfun(@(x) size(x,1),...
    HCP_aggr(r).(['Block_' test_names{t}]).UnpRew(pun_trials))),3)) '\n']);
  
  %% Average RT for Val/Inf questions
  HCP_aggr(r).(['ValInf_AvgRT_' test_names{t}]) = vertcat(HCP_aggr(r).(['Block_' test_names{t}]).ValCheckRT,...
    HCP_aggr(r).(['Block_' test_names{t}]).InfCheck1RT, HCP_aggr(r).(['Block_' test_names{t}]).InfCheck2RT,...
    HCP_aggr(r).(['Block_' test_names{t}]).InfCheckSh1RT, HCP_aggr(r).(['Block_' test_names{t}]).InfCheckSh2RT, ...
    HCP_aggr(r).(['Block_' test_names{t}]).PredCheckRT, HCP_aggr(r).(['Block_' test_names{t}]).RevalCheckRT)./RT_check_grid;
  fprintf(['Val/Inf average RT per question (per screen): mean = ' ...
    num2str(mean(HCP_aggr(r).(['ValInf_AvgRT_' test_names{t}]),'all','omitnan')) ...
    's (range = ' num2str(min(HCP_aggr(r).(['ValInf_AvgRT_' test_names{t}]),[],'all','omitnan')) ...
    ' - ' num2str(max(HCP_aggr(r).(['ValInf_AvgRT_' test_names{t}]),[],'all','omitnan')) 's)\n']);
  
  %% Questionnaires
  HCP_aggr(r).(['Catch_' test_names{t}]) = NaN(length(ques_names),1);

  for q = 1:length(ques_names)
    q_idx = is_idx(ismember(HCP_aggr(r).(['raw_' test_names{t}])(is_idx,1),ques_names{q}));

    idx = regexp(HCP_aggr(r).(['raw_' test_names{t}]){q_idx,3}.responses,':');
    HCP_aggr(r).([ques_names{q} '_' test_names{t}]) = HCP_aggr(r).(['raw_' test_names{t}]){q_idx,3}.responses(idx+1) - '0';

    if ~isempty(ques_checks{q}) % Has catch question
      if HCP_aggr(r).([ques_names{q} '_' test_names{t}])(ques_checks{q}(1)) ==  ques_checks{q}(2)
        HCP_aggr(r).(['Catch_' test_names{t}])(q) = false;
      else
        HCP_aggr(r).(['Catch_' test_names{t}])(q) = true;
      end
      
      % Remove catch answer from questionnaire vector
      idx = 1:length(HCP_aggr(r).([ques_names{q} '_' test_names{t}]));
      idx = setdiff(idx,ques_checks{q}(1));
      HCP_aggr(r).([ques_names{q} '_' test_names{t}]) = HCP_aggr(r).([ques_names{q} '_' test_names{t}])(idx);
    end
  end

  %% Failed checks line
  if ~any(HCP_aggr(r).(['Catch_' test_names{t}])) && ...
      ~any(HCP_aggr(r).(['ValInf_AvgRT_' test_names{t}])<RT_lowcutoff,'all') ...
      && ~any(HCP_aggr(r).(['ValInf_AvgRT_' test_names{t}])>RT_highcutoff,'all')
    cq_line = '';
    HCP_aggr(r).(['Exclude_' test_names{t}]) = false;
  elseif any(HCP_aggr(r).(['Catch_' test_names{t}])) ...
      && any(any(HCP_aggr(r).(['ValInf_AvgRT_' test_names{t}])<RT_lowcutoff,'all')|...
      any(HCP_aggr(r).(['ValInf_AvgRT_' test_names{t}])>RT_highcutoff,'all'))
    fprintf('FAILED BOTH CHECKS\n');
    cq_line = ' (Catch_RT_Fail)';
    HCP_aggr(r).(['Exclude_' test_names{t}]) = true;
  elseif any(HCP_aggr(r).(['Catch_' test_names{t}]))
    fprintf('FAILED CATCH QUESTIONS\n');
    cq_line = ' (Catch_Fail)';
    HCP_aggr(r).(['Exclude_' test_names{t}]) = true;
  elseif any(HCP_aggr(r).(['ValInf_AvgRT_' test_names{t}])<RT_lowcutoff,'all') || ...
      any(HCP_aggr(r).(['ValInf_AvgRT_' test_names{t}])>RT_highcutoff,'all')
    fprintf('INVALID AVG Val/Inf RT VALUES\n');  
    cq_line = ' (RT_Fail)';
    HCP_aggr(r).(['Exclude_' test_names{t}]) = true;
  else
    error('Check check fail?');
  end
  
  %% 
  figure(t_fig);
  subplot(2,1,1); hold on
  prev_times = 0;
  for p = 1:n_phases
    prev_times = prev_times + HCP_aggr(r).(['Block_' test_names{t}]).End(p) + 1;
    plot([prev_times-1 prev_times-1],ylim,'k--');
    plot([prev_times prev_times],ylim,'k--');
  end
  xlim([0 prev_times]);
  xlabel('Seconds');
  ylim([min(ev_ytick)*.5 max(ev_ytick)+min(ev_ytick)*.5]);
  yticks(ev_ytick(ytick_idx));
  yticklabels(ev_ylabel(ytick_idx));
  
  subplot(2,1,2); hold on
  prev_times = 0;
  for p = 1:n_phases
    ana_struct = ['Block' num2str(p)]; 
    prev_times = prev_times + HCP_aggr(r).(['Block_' test_names{t}]).End(p) + 1;
    plot([prev_times-1 prev_times-1],ylim,'k--');
    plot([prev_times prev_times],ylim,'k--');
  end
  xlim([0 prev_times]);
  xlabel('Seconds');
  ylabel('Total Points');
    
  set(t_fig,'Position',get(0,'Screensize'));
  saveas(t_fig,[work_folder '\Figs\Timecourse\aggr' int2str(r) test_names{t} cq_line ...
  	' timecourse (' HCP_aggr(r).(['json_' test_names{t}]) ...
    ',set' num2str(HCP_aggr(r).(['Set_' test_names{t}])) ').png']);
  close(t_fig);

  %% Click figures
  figure(click_fig); hold on   
  ylim([-HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(1),3}.viewport_size(2) 0]);
  xlim([0 HCP_aggr(r).(['raw_' test_names{t}]){phase_idx(1),3}.viewport_size(1)]);
  legend({'Block1' 'Block2' 'Block3' 'Block4' 'Block5'});
  title(['Click locations (PunP = ' HCP_aggr(r).(['PunPlanet_' test_names{t}]) ')']);
  set(click_fig,'Position',get(0,'Screensize'));
  saveas(click_fig,[work_folder '\Figs\Click location\aggr' int2str(r) test_names{t} cq_line ...
  	' click location (' HCP_aggr(r).(['json_' test_names{t}]) ...
    ',set' num2str(HCP_aggr(r).(['Set_' test_names{t}])) ').png']);
  close(click_fig);
  
  figure(normclick_fig); hold on   
  ylim([-1 1]);
  xlim([-1 1]);
  title(['Normalised click locations (PunP = ' HCP_aggr(r).(['PunPlanet_' test_names{t}]) ')']);
  legend({'Block1' 'Block2' 'Block3' 'Block4' 'Block5'});
  set(normclick_fig,'Position',get(0,'Screensize'));
  saveas(normclick_fig,[work_folder '\Figs\Click location\aggr' int2str(r) test_names{t} cq_line ...
  	' norm click location (' HCP_aggr(r).(['json_' test_names{t}]) ...
    ',set' num2str(HCP_aggr(r).(['Set_' test_names{t}])) ').png']);
  close(normclick_fig);

  %% Summary plot    
  sum_fig = HCP2_sumfig2(HCP_aggr(r).(['Block_' test_names{t}]),pun_trials);
  
  % Save
  set(sum_fig,'Position',get(0,'Screensize'));
  saveas(figure(1), [work_folder '\Figs\Summary\aggr' int2str(r) test_names{t} cq_line ...
  	' summary (' HCP_aggr(r).(['json_' test_names{t}]) ...
    ',set' num2str(HCP_aggr(r).(['Set_' test_names{t}])) ').png']);
  close all
  else
    fprintf(['HCP_aggr row' int2str(r) ',' test_names{t} ' EMPTY - skipped\n']);
  end
end % r (HCP_aggr row) loop
end % t (test) loop

%% Finish up
N = length(HCP_aggr);
excl = sum(vertcat(HCP_aggr(:).(['Exclude_' test_names{t}])));
fprintf(['\nHCP extracted. N = ' num2str(N) ' (-' ...
  num2str(excl) ' with failed checks). Valid N = ' num2str(N-excl) '. ']);

toc(start_tic_keep);
fprintf('Saving... ');

clearvars -except HCP* work_folder *keep;
save([work_folder '\HCP.mat']);

fprintf(['\nSaved (' char(datetime) '). ']);
toc(start_tic_keep);
clearvars start_tic_keep;

diary OFF