function [dist_rsq] = EyeDistanceSplit(MovieEye,Neuro_MovList,MovFR,ClipFR,Clip_blocks) % ClipEye,

if isempty(MovieEye) % && isempty(ClipEye)
    MovieEye = Spice07_MovEyeTrans;
    % ClipEye = Spice07_ClipEyeTrans;
    % EyeData = Spice07_MovEyeLate;
    % EyeData = Tor01_MovEyeTrans;
    % EyeData = Tor01_MovEyeLate;
end
% if ~exist('MovEye_Blocks','var') || isempty(MovEye_Blocks)
%     MovEye_Blocks = Spice07_Movie05Eyes.blocks;
% end
% if ~exist('ClpEye_Blocks','var') || isempty(ClpEye_Blocks)
%     ClpEye_Blocks = Spice07_05Clip800Eyes.blocks;
% end

% if ~exist('MovieData','var') || isempty(MovieData)
%     MovieData = AF_Movie05Data;
% end

if ~exist('Neuro_MovList','var') || isempty(Neuro_MovList)
    Neuro_MovList = AF_Movie05Data(94).mov_spikes(:,1);
end
% if ~exist('ClipData','var') || isempty(ClipData)
%     ClipData = AF_05Clip800Data;
% end

% Neuro_MovBlocks = MovieData.blocks;
% Neuro_ClpBlocks = ClipData.blocks;

neuro_clptrials = sum(isnan(ClipFR));

%%
neuro_movtrials = zeros(size(Neuro_MovList));
for num_t = 1:length(Neuro_MovList)
    neuro_movtrials(num_t) = ~isempty(Neuro_MovList{num_t});
end

close_mov_trials = zeros(300,floor(sum(double(neuro_movtrials))/2));
close_clp_trials = zeros(size(close_mov_trials));

distal_mov_trials = zeros(size(close_mov_trials));
distal_clp_trials = zeros(size(close_mov_trials));

for ijk = 1:300
%     cur_clpFRs = ClipFR(ijk,~isnan(ClipFR(ijk,:)));
    cur_clpFRs = ClipFR(ijk,:);
    
    cur_dist_mat = MovieEye(ijk).distmat(logical(neuro_movtrials),:);
    if sum(~cellfun(@isempty,Clip_blocks)) < 4
        close_mov_trials(ijk,:) = NaN;
        distal_mov_trials(ijk,:) = NaN;
        distal_clp_trials(ijk,:) = NaN;
        close_clp_trials(ijk,:) = NaN;
        continue
    end
    
    if sum(~cellfun(@isempty,Clip_blocks)) ~= length(cur_clpFRs) % ~= size(cur_dist_mat,2)
        try
            cur_dist_mat = cur_dist_mat(:,~cellfun(@isempty,Clip_blocks));
        catch
            disp(ijk)
        end
    end
    cur_clpFRs(isnan(cur_clpFRs)) = [];
    
    if length(cur_clpFRs) ~= size(cur_dist_mat,2)
        nan_columns = sum(isnan(cur_dist_mat));
        cur_dist_mat(:,nan_columns>3) = [];
    end
    
    mean_distances = MovieEye(ijk).mdist(logical(neuro_movtrials));
    
    md_dist = median(mean_distances);
    [~,mov_low_trials] = find(mean_distances <= md_dist); 
    [~,mov_high_trials] = find(mean_distances > md_dist); 
%     [~,indx_max] = max(cur_dist_mat');
%     [~,indx_min] = min(cur_dist_mat');
%     close_mov_trials(:,ijk) = mov_low_trials(1:size(close_mov_trials,1));
%     distal_mov_trials(:,ijk) = mov_high_trials(1:size(distal_mov_trials,1));

    close_mov_trials(ijk,:) = MovFR(ijk,mov_low_trials(1:size(close_mov_trials,2)));
    distal_mov_trials(ijk,:) = MovFR(ijk,mov_high_trials(1:size(distal_mov_trials,2)));

    pick_low_clip = zeros(1,size(close_mov_trials,2));
    for iii = 1:size(close_mov_trials,2)
        tmp_dists = cur_dist_mat(mov_low_trials(iii),:);
        [~,tmp_min] = min(tmp_dists);
        if iii > 1 && ~ismember(tmp_min,pick_low_clip)
            pick_low_clip(iii) = tmp_min;
        elseif iii == 1
            pick_low_clip(iii) = tmp_min;
        else
%             disp(ijk)
%             disp(iii)
%             pause
            tmp_dists(pick_low_clip(1:iii-1)) = NaN;
            [~,tmp_min] = nanmin(tmp_dists);
            pick_low_clip(iii) = tmp_min;
        end
        clear tmp*
    end
%     try
    close_clp_trials(ijk,:) = cur_clpFRs(pick_low_clip);
%     catch
%         disp('why!')
%     end
    
    pick_hgh_clip = zeros(1,size(distal_mov_trials,2));
    for iii = 1:size(distal_mov_trials,2)
        tmp_dists = cur_dist_mat(mov_high_trials(iii),:);
        [~,tmp_max] = max(tmp_dists);
        if iii > 1 && ~ismember(tmp_max,pick_hgh_clip) && ~ismember(tmp_max,pick_low_clip)
            pick_hgh_clip(iii) = tmp_max;
        elseif iii == 1
            pick_hgh_clip(iii) = tmp_max;
        else
            tmp_dists(pick_hgh_clip(1:iii-1)) = NaN;
            [~,tmp_max] = nanmax(tmp_dists);
            if ismember(tmp_max,pick_low_clip)
                tmp_dists(tmp_max) = NaN;
                [~,tmp_max] = nanmax(tmp_dists);
            end
            pick_hgh_clip(iii) = tmp_max;
        end
        clear tmp*
    end
%     try
    distal_clp_trials(ijk,:) = cur_clpFRs(pick_hgh_clip);
%     catch
%         disp('why!')
%     end
    clear md_dist indx_* pick_* cur_dist* mov_low_trials mov_high_trials cur_clpFRs mean_distances
end

mn_close_mov = mean(close_mov_trials,2);
mn_close_clp = mean(close_clp_trials,2);

mn_distal_mov = mean(distal_mov_trials,2);
mn_distal_clp = mean(distal_clp_trials,2);

dist_rsq = zeros(1,3)*NaN;

[~,~,~,~,tmp_reg_stats] = regress(mn_close_clp,[ones(size(mn_close_mov)) mn_close_mov]);
dist_rsq(1,1) = tmp_reg_stats(1);

[~,~,~,~,tmp_reg_stats] = regress(mn_distal_mov,[ones(size(mn_close_mov)) mn_close_mov]);
dist_rsq(1,2) = tmp_reg_stats(1);

[~,~,~,~,tmp_reg_stats] = regress(mn_distal_clp,[ones(size(mn_distal_mov)) mn_distal_mov]);
dist_rsq(1,3) = tmp_reg_stats(1);

figure(300); clf; hold on;
plot(1:300,mn_close_mov,'g');
plot(1:300,mn_distal_mov,'b');
plot(1:300,mn_close_clp,'r');
plot(1:300,mn_distal_clp,'c');

pause(.5)


clear EyeData ijk iii