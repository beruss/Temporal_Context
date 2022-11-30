%% Load AF clip800
AF_05Clip800Data = struct([]);
AF_Movie05Data = struct([]);
% AF_05FrameData = struct([]);
tic;

disp('AF')
% curdir = cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/torset01_data');
curdir = cd('/Volumes/BER_MovieTiling/MovieTiling_results/torset01_data');
load('Tor01_05Clip800Data');
load('Tor01_Movie05Data');
% load('Tor01_05FrameData');
AF_05Clip800Data = [AF_05Clip800Data Tor01_05Clip800Data];
AF_Movie05Data = [AF_Movie05Data Tor01_Movie05Data];
% AF_05FrameData = [AF_05FrameData Tor01_05FrameData];
clear Tor01*
cd(curdir)
toc/60

% cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/rhomset01_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/rhomset01_data');
load('Rhom01_05Clip800Data');
load('Rhom01_Movie05Data');
% load('Rhom01_05FrameData');
AF_05Clip800Data = [AF_05Clip800Data Rhom01_05Clip800Data];
AF_Movie05Data = [AF_Movie05Data Rhom01_Movie05Data];
% AF_05FrameData = [AF_05FrameData Rhom01_05FrameData];
clear Rhom01* 
cd(curdir)
toc/60

% cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/spiceset06_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/spiceset06_data');
load('Spice06_05Clip800Data');
load('Spice06_Movie05Data');
% load('Spice06_05FrameData');
AF_05Clip800Data = [AF_05Clip800Data Spice06_05Clip800Data];
AF_Movie05Data = [AF_Movie05Data Spice06_Movie05Data];
% AF_05FrameData = [AF_05FrameData Spice06_05FrameData];
clear Spice06* 
cd(curdir)
toc/60

% cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/spiceset07_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/spiceset07_data');
load('Spice07_05Clip800Data');
load('Spice07_Movie05Data');
% load('Spice07_05FrameData');
AF_05Clip800Data = [AF_05Clip800Data Spice07_05Clip800Data];
AF_Movie05Data = [AF_Movie05Data Spice07_Movie05Data];
% AF_05FrameData = [AF_05FrameData Spice07_05FrameData];
clear Spice07* 
cd(curdir)
toc/60

rmcells = [3 24 26 38 49 52 57];
AF_05Clip800Data(rmcells) = [];
AF_Movie05Data(rmcells) = [];
% AF_05FrameData(rmcells) = [];

%%
% SavedPath = '/Volumes/BER_MovieTiling/MovieTiling/';
%%
CombineBlocks_Spice07_220614
[Spice07_Movie05Eyes] = MovieEyeExtraction1('SpiceSet07',AllBehav,AllNeuro,5);
[Spice07_05Clip800Eyes] = ClipEyeExtraction1('SpiceSet07',AllBehav,AllNeuro,5,800);
[Spice07_05FrameEyes] = FrameEyeExtraction1('SpiceSet07',AllBehav,AllNeuro,2);

clear tanklist AllBehav AllNeuro 
%%
CombineBlocks_Spice06_220624
[Spice06_Movie05Eyes] = MovieEyeExtraction1('SpiceSet6',AllBehav,AllNeuro,5);
[Spice06_05Clip800Eyes] = ClipEyeExtraction1('SpiceSet6',AllBehav,AllNeuro,5,800);
[Spice06_05FrameEyes] = FrameEyeExtraction1('SpiceSet6',AllBehav,AllNeuro,2);

clear tanklist AllBehav AllNeuro FirstNeuro firsttank ii SavedPath
%% 
CombineBlocks_Tor01_220624
[Tor01_Movie05Eyes] = MovieEyeExtraction1('TorSet01',AllBehav,AllNeuro,5);
[Tor01_05Clip800Eyes] = ClipEyeExtraction1('TorSet01',AllBehav,AllNeuro,5,800);
[Tor01_05FrameEyes] = FrameEyeExtraction1('TorSet01',AllBehav,AllNeuro,2);

clear tanklist AllBehav AllNeuro FirstNeuro firsttank ii SavedPath
%%
CombineBlocks_Rhom01_220624
[Rhom01_Movie05Eyes] = MovieEyeExtraction1('RhomSet01',AllBehav,AllNeuro,5);
[Rhom01_05Clip800Eyes] = ClipEyeExtraction1('RhomSet01',AllBehav,AllNeuro,5,800);
[Rhom01_05FrameEyes] = FrameEyeExtraction1('RhomSet01',AllBehav,AllNeuro,2);

clear tanklist AllBehav AllNeuro FirstNeuro firsttank ii SavedPath

% %% combine data together with neural data
% for ijk = 1:length(AF_Movie05Data)
%     
%     cur_monk = AF_Movie05Data(ijk).monk;
%     cur_mov_blocks = AF_Movie05Data(ijk).blocks;
%     
%     cur_clp_blocks = AF_05Clip800Data(ijk).blocks;
%     cur_frm_blocks = AF_05FrameData(ijk).blocks;
% 
%     switch cur_monk
%         case 'TorSet01'
%             MovieEyeData = Tor01_Movie05Eyes;
%             ClipEyeData  = Tor01_05Clip800Eyes;
%             FrameEyeData = Tor01_05FrameEyes;
%         case 'RhomSet01'
%             MovieEyeData = Rhom01_Movie05Eyes;
%             ClipEyeData  = Rhom01_05Clip800Eyes;
%             FrameEyeData = Rhom01_05FrameEyes;
%         case 'SpiceSet6'
%             MovieEyeData = Spice06_Movie05Eyes;
%             ClipEyeData  = Spice06_05Clip800Eyes;
%             FrameEyeData = Spice06_05FrameEyes;
%         case 'SpiceSet07'
%             MovieEyeData = Spice07_Movie05Eyes;
%             ClipEyeData  = Spice07_05Clip800Eyes;
%             FrameEyeData = Spice07_05FrameEyes;
%     end
% 
%     eye_mov_blocks = MovieEyeData.blocks;
%     eye_mov_data = MovieEyeData.mov_eyes;
%    
%     AF_Movie05Data(ijk).mov_eyes = eye_mov_data;
%     
%     eye_clp_blocks = ClipEyeData.blocks;
%     eye_clp_data = ClipEyeData.clip_eyes;
% 
%     AF_05Clip800Data(ijk).clip_eyes = eye_clp_data;
%     
%     eye_frm_blocks = FrameEyeData.blocks;
%     eye_frm_data = FrameEyeData.clip_eyes;
%     
%     AF_05FrameData(ijk).clip_eyes = eye_frm_data;
%     
%     clear cur* MovieEyeData ClipEyeData FrameEyeData eye_*
% end

%% Lets look at the eye position

cur_mov_eyes = Spice07_Movie05Eyes.mov_eyes;

figure(400); clf; hold on;
for ijk = 1:size(cur_mov_eyes,1)
    cur_mov = cur_mov_eyes{ijk,3}/1000;
    cur_mov(cur_mov< -3500) = NaN;
    cur_mov(cur_mov> 3500) = NaN;
    tmp_mat = zeros(3500*2,3500*2);
    tmp_mov = zeros(size(cur_mov))*NaN;
    for iii = 1:length(cur_mov(1,:))
        if ~isnan(cur_mov(1,iii)) && ~isnan(cur_mov(2,iii))
            tmp_mov(:,iii) = cur_mov(:,iii);
            tmp_mat(round(cur_mov(1,iii))+3501,round(cur_mov(2,iii))+3501) = tmp_mat(round(cur_mov(1,iii))+3501,round(cur_mov(2,iii))+3501)+1;
        end
    end
    bad_times= isnan(tmp_mov(1,:));
    tmp_mov(:,isnan(tmp_mov(1,:))) = [];

    figure(400);
    plot(cur_mov(1,:),cur_mov(2,:));
    set(gca,"XLim",[-3500 3500],"YLim",[-3500 3500])
    cur_frame = getframe(gca);
    cur_frame = double(rgb2gray(frame2im(cur_frame)));
%     all_frames(:,:,ijk) = cur_frame;
    figure(600); clf;
    imagesc(cur_frame)

    figure(500); clf;
    imagesc(tmp_mat)

    disp(ijk)
    pause   
    clear cur_mov 
    clear cur_frame
end
clear cur_mov_eyes

figure(333); dscatter(tmp_mov(1,:)',tmp_mov(2,:)');
% figure(44); imshow(all_frames(:,:,1));

%%
cur_clp_eyes = Spice07_05Clip800Eyes.clip_eyes;

sess = unique([cur_clp_eyes{:,5}]);

figure(500); clf; hold on;
for ijk = 1:size(cur_clp_eyes,1)
    cur_sess = cur_clp_eyes{ijk,5};
    
    cur_clp = cur_clp_eyes{ijk,4};
    cur_clp(cur_clp< -3500000) = NaN;
    cur_clp(cur_clp> 3500000) = NaN;
    plot(cur_clp(1,:),cur_clp(2,:));
    
    if ijk == size(cur_clp_eyes,1)
        disp(ijk)
        pause   
    elseif cur_clp_eyes{ijk,5} ~= cur_clp_eyes{ijk+1,5}
        disp(ijk)
        pause   
    end
        
end

%%
prev_set = 'tempmonk_set';
for iii = 1:length(AF_Movie05Data)
    
    cur_set = AF_Movie05Data(iii).monk;
    
    if ~strncmp(cur_set,prev_set,10) 
        figure(4000+iii); clf; hold on;
%         figure(6000+iii); clf; hold on;
%         cur_mov_eyes = Spice07_Movie05Eyes.mov_eyes;
        cur_mov_eyes = AF_Movie05Data(iii).mov_eyes;
        disp(cur_set)
        
        for ijk = 1:size(cur_mov_eyes,1)
            cur_mov = cur_mov_eyes{ijk,3};
            eye_std = std(cur_mov(:));
            eye_mn = mean(cur_mov(:));
            bigeyes = eye_mn+(eye_std*1);
            
            cur_mov(cur_mov< bigeyes*-1) = NaN;
            cur_mov(cur_mov> bigeyes) = NaN;
            
            figure(4000+iii); 
            plot(cur_mov(1,:),cur_mov(2,:));
            
%             figure(6000+iii); 
%             plot(double(cur_mov(1,:)),double(cur_mov(2,:)));
            
            if nanmax(cur_mov(:)) > 10
%                 bigeyes
%                 cur_mov = cur_mov/nanmax(abs(cur_mov(:)));
                cur_mov = cur_mov/bigeyes;
                figure(5000+iii); hold on;
                plot(cur_mov(1,:),cur_mov(2,:));
                set(gca,'Ylim',[-1.5 1.5], 'Xlim',[-1.5 1.5])
                figure(4000+iii);
            end
%             cur_mov(cur_mov< -3500000) = NaN;
%             cur_mov(cur_mov> 3500000) = NaN;

        end
        pause
    end
     
    prev_set = AF_Movie05Data(iii).monk;
end

%% Extract Correct Trials from the Neural data (done and saved)
for ijk = 1:length(AF_05Clip800Data)
    disp([AF_05Clip800Data(ijk).monk '_' num2str(AF_05Clip800Data(ijk).chan(1)) '_' num2str(AF_05Clip800Data(ijk).chan(2))])
    SaveFRData_220711_temp(AF_05Clip800Data(ijk),AF_Movie05Data(ijk),-300,200,1,0,'800Trans_AT');%800Trans200
    SaveFRData_220711_temp(AF_05Clip800Data(ijk),AF_Movie05Data(ijk),-100,600,1,0,'800Late_AT');
    disp('Saved') 
end

%% Extract Eye data segments (done and saved) (redone 220714)

[Spice07_ClipEyeTrans,Spice07_MovEyeTrans] = ExtractEyeSeg_220629(Spice07_05Clip800Eyes,Spice07_Movie05Eyes,[],-300,200,3);
[Spice07_ClipEyeLate,Spice07_MovEyeLate] = ExtractEyeSeg_220629(Spice07_05Clip800Eyes,Spice07_Movie05Eyes,[],-100,600,4);

[Tor01_ClipEyeTrans,Tor01_MovEyeTrans] = ExtractEyeSeg_220629(Tor01_05Clip800Eyes,Tor01_Movie05Eyes,[],-300,200,1);
[Tor01_ClipEyeLate,Tor01_MovEyeLate] = ExtractEyeSeg_220629(Tor01_05Clip800Eyes,Tor01_Movie05Eyes,[],-100,600,2);

[Spice06_ClipEyeTrans,Spice06_MovEyeTrans] = ExtractEyeSeg_220629(Spice06_05Clip800Eyes,Spice06_Movie05Eyes,[],-300,200,33);
[Spice06_ClipEyeLate,Spice06_MovEyeLate] = ExtractEyeSeg_220629(Spice06_05Clip800Eyes,Spice06_Movie05Eyes,[],-100,600,44);

[Rhom01_ClipEyeTrans,Rhom01_MovEyeTrans] = ExtractEyeSeg_220629(Rhom01_05Clip800Eyes,Rhom01_Movie05Eyes,[],-300,200,11);
[Rhom01_ClipEyeLate,Rhom01_MovEyeLate] = ExtractEyeSeg_220629(Rhom01_05Clip800Eyes,Rhom01_Movie05Eyes,[],-100,600,22);

%% Rebuild Trial Cor Data for AF extended response (excluding transition) 
% tmpdir = '/Users/brian/Documents/Analyses/MovieTiling/data/TrialCors/';
loaddir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/Workspace';
AF_saves = [loaddir filesep 'AF_names.mat'];
load(AF_saves)
clear loaddir AF_saves

tmpdir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/TrialCors/';
AF_clpFRsLate = cell(length(AF_names),1);
AF_movFRsLate = cell(length(AF_names),1);

AF_clpFRsTrans = cell(length(AF_names),1);                      
AF_movFRsTrans = cell(length(AF_names),1); 

for ii = 1:length(AF_names)
    disp(ii)
    
%     cellname = [AF_05Clip800Data(ii).monk '_' num2str(AF_05Clip800Data(ii).chan(1)) '_' ...
%         num2str(AF_05Clip800Data(ii).chan(2)) '_cordata_notran' ];
    cellname = [AF_names{ii} '_cordata_800Late_AT' ];
    cellfile = [tmpdir cellname];
    tmpfile = matfile(cellfile);
    AF_clpFRsLate{ii} = tmpfile.crnt_clipFR_trials;
    AF_movFRsLate{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellname cellfile
    
    cellsuffix = '_cordata_800Trans_AT';
    cellfile = [tmpdir AF_names{ii} cellsuffix];
    tmpfile = matfile(cellfile);
    AF_clpFRsTrans{ii} = tmpfile.crnt_clipFR_trials;
    AF_movFRsTrans{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellsuffix cellfile
    
end
disp('done')

%% save the midlevel results (done)
datadir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/Workspace';
workdir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/Analyses/MovieTiling';

cd(datadir)
save('Tor01_EyeData','Tor01_*');
save('Rhom01_EyeData','Rhom01_*');
save('Spice06_EyeData','Spice06_*');
save('Spice07_EyeData','Spice07_*');
save('AF_Eyedata','AF_movFRs*','AF_clpFRs*','AF_Dist_*');
cd(workdir)
%% Load midlevel data (if crash or cleared) 
datadir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/Workspace';
workdir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/Analyses/MovieTiling';

cd(datadir)
load('Tor01_EyeData.mat');
load('Rhom01_EyeData.mat');
load('Spice06_EyeData.mat');
load('Spice07_EyeData.mat');
load('AF_Eyedata.mat')
cd(workdir)

%% lets try splitting the data some regions trials of high correspondence versus low 

% AF_Movie05Data(94).mov_spikes(:,1)
% EyeDistanceSplit(MovieEye,ClipEye,Neuro_MovList,MovFR,ClipFR)
% EyeDistanceSplit(Spice07_MovEyeTrans,AF_Movie05Data(94).mov_spikes(:,1),...
%                     AF_movFRsTrans{94},AF_clpFRsTrans{94});

AF_Dist_rsqs_Trans = zeros(length(AF_movFRsTrans),3)*NaN;
for ijk = 1:length(AF_movFRsTrans)
%     sprintf('%s/n',AF_Movie05Data(ijk).monk)
    switch AF_Movie05Data(ijk).monk
        case 'TorSet01'
%             sprintf('%s/n','Tor01')
            tmp_MovEyeData = Tor01_MovEyeTrans;
        case 'RhomSet01'
%             sprintf('%s/n','Rhom01')
            tmp_MovEyeData = Rhom01_MovEyeTrans;
        case 'SpiceSet6'
%             sprintf('%s/n','Spice06')
            tmp_MovEyeData = Spice06_MovEyeTrans;
        case 'SpiceSet07'
%             sprintf('%s/n','Spice07')
            tmp_MovEyeData = Spice07_MovEyeTrans;
    end
    
%     Dist_rsqs(ijk,:) = EyeDistanceSplit(tmp_MovEyeData,AF_Movie05Data(ijk).mov_spikes(:,1),...
%                                         AF_movFRsTrans{ijk},AF_clpFRsTrans{ijk},AF_05Clip800Data(ijk).blocks);
    AF_Dist_rsqs_Trans(ijk,:) = EyeDistanceSplit(tmp_MovEyeData,AF_Movie05Data(ijk).mov_spikes(:,1),...
                                        AF_movFRsTrans{ijk},AF_clpFRsTrans{ijk},AF_05Clip800Data(ijk).blocks);
    clear tmp_MovEyeData
end

bigRs = AF_Dist_rsqs_Trans(:,2) >.2;
figure(4010); clf; hold on; 
for ijk = 1:108
    if AF_Dist_rsqs_Trans(ijk,2) > .2 
        if ijk == 1 || ijk == 94 
            plot([1 2 3],[AF_Dist_rsqs_Trans(ijk,2) AF_Dist_rsqs_Trans(ijk,1) AF_Dist_rsqs_Trans(ijk,3)],'xk-')
        else
            plot([1 2 3],[AF_Dist_rsqs_Trans(ijk,2) AF_Dist_rsqs_Trans(ijk,1) AF_Dist_rsqs_Trans(ijk,3)],'o-')
        end
    end
end
set(gca,'xlim',[.5 3.5],'ylim',[0 1])
% saveas(4010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AF_Rsq_Trans_Paired_Thres2.pdf','pdf');

[~,p,~,t_stats] = ttest(AF_Dist_rsqs_Trans(bigRs,1),AF_Dist_rsqs_Trans(bigRs,3))

thist_CMc_Tr2 = histcounts(AF_Dist_rsqs_Trans(bigRs,1),[0:.025:1],'normalization','probability')*100;
thist_MMb_Tr2 = histcounts(AF_Dist_rsqs_Trans(bigRs,2),[0:.025:1],'normalization','probability')*100;
thist_CMd_Tr2 = histcounts(AF_Dist_rsqs_Trans(bigRs,3),[0:.025:1],'normalization','probability')*100;

figure(4020); clf; hold on;
bar([0.0125:.025:0.9875],thist_CMc_Tr2,'k')
set(gca,'xlim',[-.05 1.05],'ylim',[0 40])
% saveas(4020,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AF_Rsq_CMc_Trans.pdf','pdf');

figure(4024); clf; hold on;
bar([0.0125:.025:0.9875],thist_MMb_Tr2,'k')
set(gca,'xlim',[-.05 1.05],'ylim',[0 40])
% saveas(4024,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AF_Rsq_MMb_Trans.pdf','pdf');

figure(4028); clf; hold on;
bar([0.0125:.025:0.9875],thist_CMd_Tr2,'k')
set(gca,'xlim',[-.05 1.05],'ylim',[0 40])
% saveas(4028,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AF_Rsq_CMd_Trans.pdf','pdf');
clear thist* p t_stats 

%% now lets do the Late response
AF_Dist_rsqs_Late = zeros(length(AF_movFRsLate),3)*NaN;
for ijk = 1:length(AF_movFRsTrans)
%     sprintf('%s/n',AF_Movie05Data(ijk).monk)
    switch AF_Movie05Data(ijk).monk
        case 'TorSet01'
%             sprintf('%s/n','Tor01')
            tmp_MovEyeData = Tor01_MovEyeLate;
        case 'RhomSet01'
%             sprintf('%s/n','Rhom01')
            tmp_MovEyeData = Rhom01_MovEyeLate;
        case 'SpiceSet6'
%             sprintf('%s/n','Spice06')
            tmp_MovEyeData = Spice06_MovEyeLate;
        case 'SpiceSet07'
%             sprintf('%s/n','Spice07')
            tmp_MovEyeData = Spice07_MovEyeLate;
    end
    
%     Dist_rsqs(ijk,:) = EyeDistanceSplit(tmp_MovEyeData,AF_Movie05Data(ijk).mov_spikes(:,1),...
%                                         AF_movFRsTrans{ijk},AF_clpFRsTrans{ijk},AF_05Clip800Data(ijk).blocks);
    AF_Dist_rsqs_Late(ijk,:) = EyeDistanceSplit(tmp_MovEyeData,AF_Movie05Data(ijk).mov_spikes(:,1),...
                                        AF_movFRsLate{ijk},AF_clpFRsLate{ijk},AF_05Clip800Data(ijk).blocks);
    clear tmp_MovEyeData
end

bigRs = AF_Dist_rsqs_Late(:,2) >.2;
figure(5010); clf; hold on; 
for ijk = 1:108
    if AF_Dist_rsqs_Late(ijk,2) > .2 
        if ijk == 1 || ijk == 94 
            plot([1 2 3],[AF_Dist_rsqs_Late(ijk,2) AF_Dist_rsqs_Late(ijk,1) AF_Dist_rsqs_Late(ijk,3)],'xk-')
        else
            plot([1 2 3],[AF_Dist_rsqs_Late(ijk,2) AF_Dist_rsqs_Late(ijk,1) AF_Dist_rsqs_Late(ijk,3)],'o-')
        end
    end
end
set(gca,'xlim',[.5 3.5],'ylim',[0 1])
% saveas(5010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AF_Rsq_Late_Paired_Thres2.pdf','pdf');

[~,p,~,t_stats] = ttest(AF_Dist_rsqs_Late(bigRs,1),AF_Dist_rsqs_Late(bigRs,3))

thist_CMc_Tr2 = histcounts(AF_Dist_rsqs_Late(bigRs,1),[0:.025:1],'normalization','probability')*100;
thist_MMb_Tr2 = histcounts(AF_Dist_rsqs_Late(bigRs,2),[0:.025:1],'normalization','probability')*100;
thist_CMd_Tr2 = histcounts(AF_Dist_rsqs_Late(bigRs,3),[0:.025:1],'normalization','probability')*100;

figure(5020); clf; hold on;
bar([0.0125:.025:0.9875],thist_CMc_Tr2,'k')
set(gca,'xlim',[-.05 1.05])
% saveas(5020,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AF_Rsq_CMc_Late.pdf','pdf');

figure(5024); clf; hold on;
bar([0.0125:.025:0.9875],thist_MMb_Tr2,'k')
set(gca,'xlim',[-.05 1.05])
% saveas(5024,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AF_Rsq_MMb_Late.pdf','pdf');

figure(5028); clf; hold on;
bar([0.0125:.025:0.9875],thist_CMd_Tr2,'k')
set(gca,'xlim',[-.05 1.05])
% saveas(5028,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AF_Rsq_CMd_Late.pdf','pdf');
clear thist* p t_stats 
%%
% Spice07_Tran_trialssets={close_mov_trials; close_clp_trials; distal_mov_trials; distal_clp_trials};
% clear close_* distal_*

% Spice07_Late_trialssets={close_mov_trials; close_clp_trials; distal_mov_trials; distal_clp_trials};
% clear close_* distal_*

% Tor01_Tran_trialssets={close_mov_trials; close_clp_trials; distal_mov_trials; distal_clp_trials};
% clear close_* distal_*
% 
% Tor01_Late_trialssets={close_mov_trials; close_clp_trials; distal_mov_trials; distal_clp_trials};
% clear close_* distal_*

%% AM!!!

%%
%% Create Saves without External Harddrive
% loaddir = '/Users/brian/Documents/Analyses/MovieTiling/data/Workspace';
loaddir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/Workspace';

AM_saves = [loaddir filesep 'AM_names.mat'];
if ~exist(AM_saves,'file')
    AM_names = cell(length(AM_Movie05Data),1);
    for ii = 1:length(AM_names)
        AM_names{ii} = [AM_05Clip800Data(ii).monk '_' num2str(AM_05Clip800Data(ii).chan(1)) '_' num2str(AM_05Clip800Data(ii).chan(2))];

    end
    save(AM_saves,'AM_names');
    
    [AM_movbase] = MovieBaselineFR(AM_Movie05Data);
    save([loaddir filesep 'AM_movbase.mat'],'AM_movbase')
    
    [~,AM_fdprime] = PlotCatHeatmap(AM_CatData,0);
    save([loaddir filesep 'AM_fdprime.mat'],'AM_fdprime')
    close all
else
    AM_names = {};
    load(AM_saves)
    load([loaddir filesep 'AM_movbase.mat'])
    load([loaddir filesep 'AM_fdprime.mat'])
end
%%
%%
% load AM data
disp('AM')
AM_05Clip800Data = struct([]);
AM_Movie05Data = struct([]);

tic;
% curdir = cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/matset03_data');
curdir = cd('/Volumes/BER_MovieTiling/MovieTiling_results/matset03_data');
load('Mat03_05Clip800Data');
load('Mat03_Movie05Data');
AM_05Clip800Data = [AM_05Clip800Data Mat03_05Clip800Data];
AM_Movie05Data = [AM_Movie05Data Mat03_Movie05Data];
clear Mat03* 
cd(curdir)

% curdir = cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/dangoset01_data');
curdir = cd('/Volumes/BER_MovieTiling/MovieTiling_results/dangoset01_data');
load('Dango01_05Clip800Data');
load('Dango01_Movie05Data');
AM_05Clip800Data = [AM_05Clip800Data Dango01_05Clip800Data];
AM_Movie05Data = [AM_Movie05Data Dango01_Movie05Data];
clear Dango01* 
cd(curdir)
toc/60

rmcells = [29 38 72 73 75 76 77 78  80 83 86 87];
AM_05Clip800Data(rmcells) = [];
AM_Movie05Data(rmcells) = [];

rmcells2 = [36 38 47 67 73 75];
AM_05Clip800Data(rmcells2) = [];
AM_Movie05Data(rmcells2) = [];

clear curdir ans rmcells rmcells2

%%
CombineBlocks_Matcha03_220715
[Mat03_Movie05Eyes] = MovieEyeExtraction1('Matcha03',AllBehav,AllNeuro,5);
[Mat03_05Clip800Eyes] = ClipEyeExtraction1('Matcha03',AllBehav,AllNeuro,5,800);
[Mat03_05FrameEyes] = FrameEyeExtraction1('Matcha03',AllBehav,AllNeuro,2);

clear tanklist AllBehav AllNeuro FirstNeuro firsttank ii SavedPath 

CombineBlocks_Dango01_071522
[Dango01_Movie05Eyes] = MovieEyeExtraction1('Dango01',AllBehav,AllNeuro,5);
[Dango01_05Clip800Eyes] = ClipEyeExtraction1('Dango01',AllBehav,AllNeuro,5,800);
[Dango01_05FrameEyes] = FrameEyeExtraction1('Dango01',AllBehav,AllNeuro,2);

clear tanklist AllBehav AllNeuro FirstNeuro firsttank ii SavedPath

%% Look at AM eye data
figure(4000); clf; hold on;
%         figure(6000+iii); clf; hold on;
%         cur_mov_eyes = Spice07_Movie05Eyes.mov_eyes;
% cur_mov_eyes = Dango01_Movie05Eyes.mov_eyes;
cur_mov_eyes = Mat03_05Clip800Eyes.clip_eyes;
eyeindex = 4; % 3=mov 4=clips
for ijk = 1:size(cur_mov_eyes,1)
    cur_mov = cur_mov_eyes{ijk,eyeindex};
    eye_std = std(cur_mov(:));
    eye_mn = mean(cur_mov(:));
    bigeyes = eye_mn+(eye_std*1);

    cur_mov(cur_mov< bigeyes*-1) = NaN;
    cur_mov(cur_mov> bigeyes) = NaN;

    figure(4000); 
    plot(cur_mov(1,:),cur_mov(2,:));
%     pause
%     if nanmax(cur_mov(:)) > 10
%         cur_mov = cur_mov/bigeyes;
%         figure(5000); hold on;
%         plot(cur_mov(1,:),cur_mov(2,:));
%         set(gca,'Ylim',[-1.5 1.5], 'Xlim',[-1.5 1.5])
%         figure(4000);
%     end
end
clear eye_mn bigeyes cur_mov eye_std eyeindex ijk loaddir cur_mov_eyes curdir ans 
%% Extract Correct Trials from the Neural data (done and saved)
for ijk = 1:length(AM_05Clip800Data)
    disp([AM_05Clip800Data(ijk).monk '_' num2str(AM_05Clip800Data(ijk).chan(1)) '_' num2str(AM_05Clip800Data(ijk).chan(2))])
    SaveFRData_220711_temp(AM_05Clip800Data(ijk),AM_Movie05Data(ijk),-300,200,1,0,'800Trans_AT');%800Trans200
    SaveFRData_220711_temp(AM_05Clip800Data(ijk),AM_Movie05Data(ijk),-100,600,1,0,'800Late_AT');
    disp('Saved') 
end

%% Extract Eye data segments (done and saved) (redone 220714)

[Mat03_ClipEyeTrans,Mat03_MovEyeTrans] = ExtractEyeSeg_220629(Mat03_05Clip800Eyes,Mat03_Movie05Eyes,[],-300,200,3);
[Mat03_ClipEyeLate,Mat03_MovEyeLate] = ExtractEyeSeg_220629(Mat03_05Clip800Eyes,Mat03_Movie05Eyes,[],-100,600,4);

[Dango01_ClipEyeTrans,Dango01_MovEyeTrans] = ExtractEyeSeg_220629(Dango01_05Clip800Eyes,Dango01_Movie05Eyes,[],-300,200,1);
[Dango01_ClipEyeLate,Dango01_MovEyeLate] = ExtractEyeSeg_220629(Dango01_05Clip800Eyes,Dango01_Movie05Eyes,[],-100,600,2);

%% save the midlevel results (done)
datadir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/Workspace';
workdir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/Analyses/MovieTiling';

cd(datadir)
save('Dango01_EyeData','Dango01_*');
save('Mat03_EyeData','Mat03_*');
save('AM_Eyedata','AM_movFRs*','AM_clpFRs*','AM_Dist_*');
cd(workdir)
%% Load midlevel data (if crash or cleared) 
datadir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/Workspace';
workdir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/Analyses/MovieTiling';

cd(datadir)
load('Dango01_EyeData.mat');
load('Mat03_EyeData.mat');
load('AM_EyeData.mat');
cd(workdir)

%% Rebuild Trial Cor Data for AF extended response (excluding transition) 
% tmpdir = '/Users/brian/Documents/Analyses/MovieTiling/data/TrialCors/';
tmpdir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/TrialCors/';
AM_clpFRsLate = cell(length(AM_names),1);
AM_movFRsLate = cell(length(AM_names),1);

AM_clpFRsTrans = cell(length(AM_names),1);                      
AM_movFRsTrans = cell(length(AM_names),1); 

for ii = 1:length(AM_names)
    disp(ii)
    
%     cellname = [AM_05Clip800Data(ii).monk '_' num2str(AM_05Clip800Data(ii).chan(1)) '_' ...
%         num2str(AM_05Clip800Data(ii).chan(2)) '_cordata_notran' ];
    cellname = [AM_names{ii} '_cordata_800Late_AT' ];
    cellfile = [tmpdir cellname];
    tmpfile = matfile(cellfile);
    AM_clpFRsLate{ii} = tmpfile.crnt_clipFR_trials;
    AM_movFRsLate{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellname cellfile
    
    cellsuffix = '_cordata_800Trans_AT';
    cellfile = [tmpdir AM_names{ii} cellsuffix];
    tmpfile = matfile(cellfile);
    AM_clpFRsTrans{ii} = tmpfile.crnt_clipFR_trials;
    AM_movFRsTrans{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellsuffix cellfile
    
end
disp('done')

%% lets try splitting the data some regions trials of high correspondence versus low 

% AF_Movie05Data(94).mov_spikes(:,1)
% EyeDistanceSplit(MovieEye,ClipEye,Neuro_MovList,MovFR,ClipFR)
% EyeDistanceSplit(Spice07_MovEyeTrans,AF_Movie05Data(94).mov_spikes(:,1),...
%                     AF_movFRsTrans{94},AF_clpFRsTrans{94});

AM_Dist_rsqs_Trans = zeros(length(AM_movFRsTrans),3)*NaN;
for ijk = 1:length(AM_movFRsTrans)
%     sprintf('%s/n',AF_Movie05Data(ijk).monk)
    switch AM_Movie05Data(ijk).monk
        case 'Matcha03'
%             sprintf('%s/n','Tor01')
            tmp_MovEyeData = Mat03_MovEyeTrans;
        case 'Dango01'
%             sprintf('%s/n','Rhom01')
            tmp_MovEyeData = Dango01_MovEyeTrans;
    end
    
%     Dist_rsqs(ijk,:) = EyeDistanceSplit(tmp_MovEyeData,AM_Movie05Data(ijk).mov_spikes(:,1),...
%                                         AM_movFRsTrans{ijk},AM_clpFRsTrans{ijk},AM_05Clip800Data(ijk).blocks);
    AM_Dist_rsqs_Trans(ijk,:) = EyeDistanceSplit(tmp_MovEyeData,AM_Movie05Data(ijk).mov_spikes(:,1),...
                                        AM_movFRsTrans{ijk},AM_clpFRsTrans{ijk},AM_05Clip800Data(ijk).blocks);
    clear tmp_MovEyeData
end

%%
bigRs = AM_Dist_rsqs_Trans(:,2) >.2;
figure(5010); clf; hold on; 
for ijk = 1:72
    if AM_Dist_rsqs_Trans(ijk,2) > .2 
        if ijk == 1 || ijk == 94 
            plot([1 2 3],[AM_Dist_rsqs_Trans(ijk,2) AM_Dist_rsqs_Trans(ijk,1) AM_Dist_rsqs_Trans(ijk,3)],'xk-')
        else
            plot([1 2 3],[AM_Dist_rsqs_Trans(ijk,2) AM_Dist_rsqs_Trans(ijk,1) AM_Dist_rsqs_Trans(ijk,3)],'o-')
        end
    end
end
set(gca,'xlim',[.5 3.5],'ylim',[0 1])
% saveas(5010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AM_Rsq_Trans_Paired_Thres2.pdf','pdf');

[~,p,~,t_stats] = ttest(AM_Dist_rsqs_Trans(bigRs,1),AM_Dist_rsqs_Trans(bigRs,3))

thist_CMc_Tr2 = histcounts(AM_Dist_rsqs_Trans(bigRs,1),[0:.025:1],'normalization','probability')*100;
thist_MMb_Tr2 = histcounts(AM_Dist_rsqs_Trans(bigRs,2),[0:.025:1],'normalization','probability')*100;
thist_CMd_Tr2 = histcounts(AM_Dist_rsqs_Trans(bigRs,3),[0:.025:1],'normalization','probability')*100;

figure(5020); clf; hold on;
bar([0.0125:.025:0.9875],thist_CMc_Tr2,'k')
set(gca,'xlim',[-.05 1.05])
% saveas(5020,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AM_Rsq_CMc_Trans.pdf','pdf');

figure(5024); clf; hold on;
bar([0.0125:.025:0.9875],thist_MMb_Tr2,'k')
set(gca,'xlim',[-.05 1.05])
% saveas(5024,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AM_Rsq_MMb_Trans.pdf','pdf');

figure(5028); clf; hold on;
bar([0.0125:.025:0.9875],thist_CMd_Tr2,'k')
set(gca,'xlim',[-.05 1.05])
% saveas(5028,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AM_Rsq_CMd_Trans.pdf','pdf');
clear thist* p t_stats 

%%
AM_Dist_rsqs_Late = zeros(length(AM_movFRsLate),3)*NaN;
for ijk = 1:length(AM_movFRsTrans)
%     sprintf('%s/n',AM_Movie05Data(ijk).monk)
    switch AM_Movie05Data(ijk).monk
       case 'Matcha03'
%             sprintf('%s/n','Tor01')
            tmp_MovEyeData = Mat03_MovEyeTrans;
        case 'Dango01'
%             sprintf('%s/n','Rhom01')
            tmp_MovEyeData = Dango01_MovEyeTrans;
    end
    
%     Dist_rsqs(ijk,:) = EyeDistanceSplit(tmp_MovEyeData,AM_Movie05Data(ijk).mov_spikes(:,1),...
%                                         AM_movFRsTrans{ijk},AM_clpFRsTrans{ijk},AM_05Clip800Data(ijk).blocks);
    AM_Dist_rsqs_Late(ijk,:) = EyeDistanceSplit(tmp_MovEyeData,AM_Movie05Data(ijk).mov_spikes(:,1),...
                                        AM_movFRsLate{ijk},AM_clpFRsLate{ijk},AM_05Clip800Data(ijk).blocks);
    clear tmp_MovEyeData
end

bigRs = AM_Dist_rsqs_Late(:,2) >.2;
figure(5010); clf; hold on; 
for ijk = 1:72
    if AM_Dist_rsqs_Late(ijk,2) > .2 
        if ijk == 1 || ijk == 94 
            plot([1 2 3],[AM_Dist_rsqs_Late(ijk,2) AM_Dist_rsqs_Late(ijk,1) AM_Dist_rsqs_Late(ijk,3)],'xk-')
        else
            plot([1 2 3],[AM_Dist_rsqs_Late(ijk,2) AM_Dist_rsqs_Late(ijk,1) AM_Dist_rsqs_Late(ijk,3)],'o-')
        end
    end
end
set(gca,'xlim',[.5 3.5],'ylim',[0 1])
% saveas(5010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AM_Rsq_Late_Paired_Thres2.pdf','pdf');

[~,p,~,t_stats] = ttest(AM_Dist_rsqs_Late(bigRs,1),AM_Dist_rsqs_Late(bigRs,3))

thist_CMc_Tr2 = histcounts(AM_Dist_rsqs_Late(bigRs,1),[0:.025:1],'normalization','probability')*100;
thist_MMb_Tr2 = histcounts(AM_Dist_rsqs_Late(bigRs,2),[0:.025:1],'normalization','probability')*100;
thist_CMd_Tr2 = histcounts(AM_Dist_rsqs_Late(bigRs,3),[0:.025:1],'normalization','probability')*100;

figure(5020); clf; hold on;
bar([0.0125:.025:0.9875],thist_CMc_Tr2,'k')
set(gca,'xlim',[-.05 1.05])
% saveas(5020,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AM_Rsq_CMc_Late.pdf','pdf');

figure(5024); clf; hold on;
bar([0.0125:.025:0.9875],thist_MMb_Tr2,'k')
set(gca,'xlim',[-.05 1.05])
% saveas(5024,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AM_Rsq_MMb_Late.pdf','pdf');

figure(5028); clf; hold on;
bar([0.0125:.025:0.9875],thist_CMd_Tr2,'k')
set(gca,'xlim',[-.05 1.05])
% saveas(5028,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/AM_Rsq_CMd_Late.pdf','pdf');
clear thist* p t_stats 

%% Create Combined Figures
Both_Late_Rsq = [AF_Dist_rsqs_Late; AM_Dist_rsqs_Late];

bigRs = Both_Late_Rsq(:,2) >.2;
figure(5010); clf; hold on; 
for ijk = 1:72
    if Both_Late_Rsq(ijk,2) > .2 
        if ijk == 1 || ijk == 94 
            plot([1 2 3],[Both_Late_Rsq(ijk,2) Both_Late_Rsq(ijk,1) Both_Late_Rsq(ijk,3)],'xk-')
        else
            plot([1 2 3],[Both_Late_Rsq(ijk,2) Both_Late_Rsq(ijk,1) Both_Late_Rsq(ijk,3)],'o-')
        end
    end
end
set(gca,'xlim',[.5 3.5],'ylim',[0 1])
% saveas(5010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Both_Rsq_Late_Paired_Thres2.pdf','pdf');

anova_context = [ones(size(Both_Late_Rsq(bigRs,1),1),1)*1; ones(size(Both_Late_Rsq(bigRs,1),1),1)*2; ones(size(Both_Late_Rsq(bigRs,1),1),1)*3]; 
anova_data = [Both_Late_Rsq(bigRs,2); Both_Late_Rsq(bigRs,1); Both_Late_Rsq(bigRs,3)];
[p,atab] = anova1(anova_data,anova_context); 
% [p,atab] = anovan(anova_data,anova_context,'varnames',{'Context'},'model','full','sstype', 3)

[~,p,~,t_stats] = ttest(Both_Late_Rsq(bigRs,1),Both_Late_Rsq(bigRs,3))

thist_CMc_Tr2 = histcounts(Both_Late_Rsq(bigRs,1),[0:.025:1],'normalization','probability')*100;
thist_MMb_Tr2 = histcounts(Both_Late_Rsq(bigRs,2),[0:.025:1],'normalization','probability')*100;
thist_CMd_Tr2 = histcounts(Both_Late_Rsq(bigRs,3),[0:.025:1],'normalization','probability')*100;

figure(5020); clf; hold on;
bar([0.0125:.025:0.9875],thist_CMc_Tr2,'k')
set(gca,'xlim',[-.05 1.05],'ylim',[0 28])
% saveas(5020,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Both_Rsq_CMc_Late.pdf','pdf');

figure(5024); clf; hold on;
bar([0.0125:.025:0.9875],thist_MMb_Tr2,'k')
set(gca,'xlim',[-.05 1.05],'ylim',[0 28])
% saveas(5024,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Both_Rsq_MMb_Late.pdf','pdf');

figure(5028); clf; hold on;
bar([0.0125:.025:0.9875],thist_CMd_Tr2,'k')
set(gca,'xlim',[-.05 1.05],'ylim',[0 28])
% saveas(5028,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Both_Rsq_CMd_Late.pdf','pdf');
clear thist* p t_stats 

%% And Transition 
Both_Trans_Rsq = [AF_Dist_rsqs_Trans; AM_Dist_rsqs_Trans];

bigRs = Both_Trans_Rsq(:,2) >.2;
figure(5010); clf; hold on; 
for ijk = 1:72
    if Both_Trans_Rsq(ijk,2) > .2 
        if ijk == 1 || ijk == 94 
            plot([1 2 3],[Both_Trans_Rsq(ijk,2) Both_Trans_Rsq(ijk,1) Both_Trans_Rsq(ijk,3)],'xk-')
        else
            plot([1 2 3],[Both_Trans_Rsq(ijk,2) Both_Trans_Rsq(ijk,1) Both_Trans_Rsq(ijk,3)],'o-')
        end
    end
end
set(gca,'xlim',[.5 3.5],'ylim',[0 1])
% saveas(5010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Both_Rsq_Trans_Paired_Thres2.pdf','pdf');

anova_context = [ones(size(Both_Trans_Rsq(bigRs,1),1),1)*1; ones(size(Both_Trans_Rsq(bigRs,1),1),1)*2; ones(size(Both_Trans_Rsq(bigRs,1),1),1)*3]; 
anova_data = [Both_Trans_Rsq(bigRs,2); Both_Trans_Rsq(bigRs,1); Both_Trans_Rsq(bigRs,3)];
[p,atab] = anova1(anova_data,anova_context);

[~,p,~,t_stats] = ttest(Both_Trans_Rsq(bigRs,1),Both_Trans_Rsq(bigRs,3))

thist_CMc_Tr2 = histcounts(Both_Trans_Rsq(bigRs,1),[0:.025:1],'normalization','probability')*100;
thist_MMb_Tr2 = histcounts(Both_Trans_Rsq(bigRs,2),[0:.025:1],'normalization','probability')*100;
thist_CMd_Tr2 = histcounts(Both_Trans_Rsq(bigRs,3),[0:.025:1],'normalization','probability')*100;

figure(5020); clf; hold on;
bar([0.0125:.025:0.9875],thist_CMc_Tr2,'k')
set(gca,'xlim',[-.05 1.05],'ylim',[0 28])
% saveas(5020,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Both_Rsq_CMc_Trans.pdf','pdf');

figure(5024); clf; hold on;
bar([0.0125:.025:0.9875],thist_MMb_Tr2,'k')
set(gca,'xlim',[-.05 1.05],'ylim',[0 28])
% saveas(5024,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Both_Rsq_MMb_Trans.pdf','pdf');

figure(5028); clf; hold on;
bar([0.0125:.025:0.9875],thist_CMd_Tr2,'k')
set(gca,'xlim',[-.05 1.05],'ylim',[0 28])
% saveas(5028,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Both_Rsq_CMd_Trans.pdf','pdf');
clear thist* p t_stats 

%% save the Combined EyePosition results (done)
datadir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/Workspace';
workdir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/Analyses/MovieTiling';

cd(datadir)
save('Both_EyeAnalysis','Both_*');
save('Eye_Rsqs_All','AF_Dist_rsqs*','AM_Dist_rsqs*')
cd(workdir)
%% Load midlevel data (if crash or cleared) 
datadir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/Workspace';
workdir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/Analyses/MovieTiling';

cd(datadir)
load('Both_EyeAnalysis.mat');
load('Eye_Rsqs_All.mat')
cd(workdir)

%% get the amount of faces in the movie 5
rgr_dir = '/Volumes/GoogleDrive/My Drive/NIH/Regressors';
workdir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/Analyses/MovieTiling';

cd(rgr_dir)
load Mov1to6_10fps_rgr.mat
cd(work_dir)

figure(12); clf; hold on;
m5 = RAWRGR.rgrs{5};
face1 = m5{7};
face2 = m5{8};
face3 = face1 + face2;
plot([1:1200]/4,face1,'-g')
plot([1:1200]/4,face2,'-y')
plot([1:1200]/4,face3,'-r')

sum(face3 == 0)/ 1200

%% Plot density of eye positions

%% Lets look at the eye position

cur_mov_eyes = Spice07_Movie05Eyes.mov_eyes;
all_frames = zeros(840,1120,size(cur_mov_eyes,1));

[~,numSamples] = cellfun(@size,cur_mov_eyes(:,3));
all_trials = zeros(2,sum(numSamples))*NaN;
% figure(400); clf; hold on;
sampleIDX = 0;
for ijk = 1:size(cur_mov_eyes,1)
    cur_mov = cur_mov_eyes{ijk,3}/10000;
    cur_mov(cur_mov< -350) = NaN;
    cur_mov(cur_mov> 350) = NaN;
    mov_center = nanmean(cur_mov,2);

    cur_mov(1,:) = cur_mov(1,:) - mov_center(1);
    cur_mov(2,:) = cur_mov(2,:) - mov_center(2);

    cur_mov_scaled = (cur_mov/180)*7.5;
%     tmp_mat = zeros(350*2,350*2);
    for iii = 1:length(cur_mov(1,:))
        sampleIDX = sampleIDX+1;
        if ~isnan(cur_mov(1,iii)) && ~isnan(cur_mov(2,iii))
            all_trials(:,sampleIDX) = cur_mov_scaled(:,iii);
%             tmp_mat(round(cur_mov(1,iii))+351,round(cur_mov(2,iii))+351) = tmp_mat(round(cur_mov(1,iii))+351,round(cur_mov(2,iii))+351)+1;
        end
    end


    disp(ijk)
%     pause   
    clear cur_mov tmp_mat
end
clear cur_mov_eyes numSamples

all_trials(:,isnan(all_trials(1,:))) = [];

figure(333); clf; dscatter(all_trials(1,:)',all_trials(2,:)');
% figure(44); imshow(all_frames(:,:,1));
rectangle('Position',[-7.5 -5.6  15 11.25],'EdgeColor',[1 1 1])

set(gca,'Xlim',[-20 20],'YLim',[-15 15])
set(gca,'Color','k')
colormap('hot')
saveas(333,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Spice_EyeDensity_FreeView20_hot.tiff','tiffn');
%%
figure(444); clf; hold on;
x_std = std(all_trials(1,:));
plot(1:length(all_trials),all_trials(1,:),'-G')
plot(1:length(all_trials),ones(1,length(all_trials))*x_std*2.5,'-R')
plot(1:length(all_trials),ones(1,length(all_trials))*x_std*-2.5,'-R')

% plot(1:length(all_trials),all_trials(2,:),'-B')
set(gca,'Xlim',[-5 1000000])



