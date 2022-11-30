%% NEW SCRATCH Pad to clean up analyses 
% cd /Users/brian/Documents/Analyses/MovieTiling/Analyses/MovieTiling/
cd '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/Analyses/MovieTiling'

%% AF  --  Create Saves without External Harddrive
% loaddir = '/Users/brian/Documents/Analyses/MovieTiling/data/Workspace';
loaddir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/Workspace';

AF_saves = [loaddir filesep 'AF_names.mat'];
if ~exist(AF_saves,'file')
    AF_names = cell(length(AF_Movie05Data),1);
    for ii = 1:length(AF_names)
        AF_names{ii} = [AF_05Clip800Data(ii).monk '_' num2str(AF_05Clip800Data(ii).chan(1)) '_' num2str(AF_05Clip800Data(ii).chan(2))];

    end
    save(AF_saves,'AF_names');
    
    [AF_movbase] = MovieBaselineFR(AF_Movie05Data);
    save([loaddir filesep 'AF_movbase.mat'],'AF_movbase')
    
    [~,AF_fdprime] = PlotCatHeatmap(CatData,0);
    save([loaddir filesep 'AF_fdprime.mat'],'AF_fdprime')
else
    AF_names = {};
    load(AF_saves)
    load([loaddir filesep 'AF_movbase.mat'])
    load([loaddir filesep 'AF_fdprime.mat'])
end

clear AF_saves loaddir   

%% AM -- Create Saves without External Harddrive
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
    
    [~,AM_fdprime] = PlotCatHeatmap(CatData,0);
    save([loaddir filesep 'AM_fdprime.mat'],'AM_fdprime')
    close all
else
    AM_names = {};
    load(AM_saves)
    load([loaddir filesep 'AM_movbase.mat'])
    load([loaddir filesep 'AM_fdprime.mat'])
end

%% AF --  Rebuild Trial Cor Data for AF extended response (excluding transition) 
% tmpdir = '/Users/brian/Documents/Analyses/MovieTiling/data/TrialCors/';
tmpdir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/TrialCors/';
AF_clpFRs800 = cell(length(AF_names),1);
AF_movFRs800 = cell(length(AF_names),1);
AF_clpFRs250 = cell(length(AF_names),1);
AF_movFRs250 = cell(length(AF_names),1);
AF_clpFRs100 = cell(length(AF_names),1);
AF_movFRs100 = cell(length(AF_names),1);
AF_clpFRsFRM = cell(length(AF_names),1);
AF_movFRsFRM = cell(length(AF_names),1);

% AF_clpFRs800t000 = cell(length(AF_names),1);
% AF_movFRs800t000 = cell(length(AF_names),1);

AF_clpFRs800t200 = cell(length(AF_names),1);                      
AF_movFRs800t200 = cell(length(AF_names),1); 

for ii = 1:length(AF_names)
    disp(ii)
    
%     cellname = [AF_05Clip800Data(ii).monk '_' num2str(AF_05Clip800Data(ii).chan(1)) '_' ...
%         num2str(AF_05Clip800Data(ii).chan(2)) '_cordata_notran' ];
    cellname = [AF_names{ii} '_cordata_notran' ];
    cellfile = [tmpdir cellname];
    tmpfile = matfile(cellfile);
    AF_clpFRs800{ii} = tmpfile.crnt_clipFR_trials;
    AF_movFRs800{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellname cellfile
    
    cellname = [AF_names{ii} '_cordata_250notran' ];
    cellfile = [tmpdir cellname];
    tmpfile = matfile(cellfile);
    AF_clpFRs250{ii} = tmpfile.crnt_clipFR_trials;
    AF_movFRs250{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellname cellfile tmpout
    
    cellname = [AF_names{ii} '_cordata_100notran' ];
    cellfile = [tmpdir cellname];
    tmpfile = matfile(cellfile);
    AF_clpFRs100{ii} = tmpfile.crnt_clipFR_trials;
    AF_movFRs100{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellname cellfile tmpout
    
    cellname = [AF_names{ii} '_cordata_FRMnotran' ];
    cellfile = [tmpdir cellname];
    tmpfile = matfile(cellfile);
    AF_clpFRsFRM{ii} = tmpfile.crnt_clipFR_trials;
    AF_movFRsFRM{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellname cellfile tmpout
    
    
    cellsuffix = '_cordata_800Trans200';
    cellfile = [tmpdir AF_names{ii} cellsuffix];
    tmpfile = matfile(cellfile);
    AF_clpFRs800t200{ii} = tmpfile.crnt_clipFR_trials;
    AF_movFRs800t200{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellsuffix cellfile
    
%     cellsuffix = '_cordata_800wtran';
%     cellfile = [tmpdir AF_names{ii} cellsuffix];
%     tmpfile = matfile(cellfile);
%     AF_clpFRs800t000{ii} = tmpfile.crnt_clipFR_trials;
%     AF_movFRs800t000{ii} = tmpfile.crnt_movFR_trials;
%     clear tmpfile cellname cellfile tmpout
end
disp('done')



%% AM --  Rebuild Trial Cor Data for AM extended response (excluding transition) 
tmpdir = '/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/TrialCors/';

AM_clpFRs800 = cell(length(AM_names),1);
AM_movFRs800 = cell(length(AM_names),1);
AM_clpFRs250 = cell(length(AM_names),1);
AM_movFRs250 = cell(length(AM_names),1);
AM_clpFRs100 = cell(length(AM_names),1);
AM_movFRs100 = cell(length(AM_names),1);
AM_clpFRsFRM = cell(length(AM_names),1);
AM_movFRsFRM = cell(length(AM_names),1);
AM_clpFRs800t200 = cell(length(AM_names),1);                      
AM_movFRs800t200 = cell(length(AM_names),1); 

for ii = 1:length(AM_names)
    disp(ii)

    cellname = [AM_names{ii} '_cordata_notran' ];
    cellfile = [tmpdir cellname];
    tmpfile = matfile(cellfile);
    AM_clpFRs800{ii} = tmpfile.crnt_clipFR_trials;
    AM_movFRs800{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellname cellfile
    
    cellname = [AM_names{ii} '_cordata_250notran' ];
    cellfile = [tmpdir cellname];
    tmpfile = matfile(cellfile);
    AM_clpFRs250{ii} = tmpfile.crnt_clipFR_trials;
    AM_movFRs250{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellname cellfile tmpout
    
    cellname = [AM_names{ii} '_cordata_100notran' ];
    cellfile = [tmpdir cellname];
    tmpfile = matfile(cellfile);
    AM_clpFRs100{ii} = tmpfile.crnt_clipFR_trials;
    AM_movFRs100{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellname cellfile tmpout
    
    cellname = [AM_names{ii} '_cordata_FRMnotran' ];
    cellfile = [tmpdir cellname];
    tmpfile = matfile(cellfile);
    AM_clpFRsFRM{ii} = tmpfile.crnt_clipFR_trials;
    AM_movFRsFRM{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellname cellfile tmpout
    
%     cellsuffix = '_cordata_800wtran';
%     cellfile = [tmpdir AM_names{ii} cellsuffix];
%     tmpfile = matfile(cellfile);
%     AM_clpFRs800t000{ii} = tmpfile.crnt_clipFR_trials;
%     AM_movFRs800t000{ii} = tmpfile.crnt_movFR_trials;
%     clear tmpfile cellname cellfile
    
    cellsuffix = '_cordata_800Trans200';
    cellfile = [tmpdir AM_names{ii} cellsuffix];
    tmpfile = matfile(cellfile);
    AM_clpFRs800t200{ii} = tmpfile.crnt_clipFR_trials;
    AM_movFRs800t200{ii} = tmpfile.crnt_movFR_trials;
    clear tmpfile cellsuffix cellfile
end
disp('done')

%% Combine! 
clpFRs800 = [AF_clpFRs800; AM_clpFRs800];
clpFRs100 = [AF_clpFRs100; AM_clpFRs100];
clpFRs250 = [AF_clpFRs250; AM_clpFRs250];
clpFRsFRM = [AF_clpFRsFRM; AM_clpFRsFRM];
% clpFRs800t000 = [AF_clpFRs800t000; AM_clpFRs800t000'];
clpFRs800t200 = [AF_clpFRs800t200; AM_clpFRs800t200];


movFRs800 = [AF_movFRs800; AM_movFRs800];
movFRs100 = [AF_movFRs100; AM_movFRs100];
movFRs250 = [AF_movFRs250; AM_movFRs250];
movFRsFRM = [AF_movFRsFRM; AM_movFRsFRM];
% movFRs800t000 = [AF_movFRs800t000; AM_movFRs800t000'];
movFRs800t200 = [AF_movFRs800t200; AM_movFRs800t200];

Face_dprime = [AF_fdprime; AM_fdprime];
MovieBase = [AF_movbase; AM_movbase];

All_names = [AF_names; AM_names];

clear AF* AM* 

%% Get correlation for figure 2 figure
c94_movFRM = nanmean(movFRsFRM{94},2);
c94_clpFRM = nanmean(clpFRsFRM{94},2);
[c94_frm_R,c94_frm_p] = corrcoef(c94_movFRM,c94_clpFRM);

c94_mov800 = nanmean(movFRs800{94},2);
c94_clp800 = nanmean(clpFRs800{94},2);
[c94_800_R,c94_800_p] = corrcoef(c94_mov800,c94_clp800);

c94_mov800t200 = nanmean(movFRs800t200{94},2);
c94_clp800t200 = nanmean(clpFRs800t200{94},2);
[c94_t200_R,c94_t200_p] = corrcoef(c94_mov800t200,c94_clp800t200)
%% Get Rsqaured
% lims = [-.25 1.1 35 25 0];

if ~exist('MovieBase','var')
    load([loaddir filesep 'AF_movbase.mat'])
    load([loaddir filesep 'AM_movbase.mat'])
    MovieBase = [AF_movbase; AM_movbase];
    clear AF_movbase AM_movbase
end
lims = [-.25 1.1 15 12 0];

[~,Half_Rsq800,Full_Rsq800,ClipSet] = SlopeStats_half(movFRs800,clpFRs800,MovieBase,[],'800',22,lims);
% saveas(7522,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_800_hist_nothres.pdf','pdf');
% saveas(7622,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_half_Rsq_800_hist_nothres.pdf','pdf');

lims = [-.25 1.1 25 30 0];
[~,Half_Rsq800t200,Full_Rsq800t200] = SlopeStats_half(movFRs800t200,clpFRs800t200,MovieBase,[],'800t200',2,lims);
% saveas(7502,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_800t200_hist_nothres.pdf','pdf');

lims = [-.25 1.1 25 60 0];
[~,Half_RsqFRM,Full_RsqFRM] = SlopeStats_half(movFRsFRM,clpFRsFRM,MovieBase,[],'FRAME',4,lims);
% saveas(7504,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_FRAME_hist_nothres.pdf','pdf');

lims = [-.25 1.1 25 30 0];
[~, Half_Rsq250, Full_Rsq250] = SlopeStats_half(movFRs250,clpFRs250,MovieBase,[],'250',25,lims);
% saveas(7525,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_250_hist_nothres.pdf','pdf');

lims = [-.25 1.1 25 30 0];
[~, Half_Rsq100, Full_Rsq100] = SlopeStats_half(movFRs100,clpFRs100,MovieBase,[],'100',10,lims);
% saveas(7510,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_100_hist_nothres.pdf','pdf');

close all

%% Now Make combined figures 

bigRs = Half_Rsq800t200(:,2) >.2;

thist_hlf_r2 = histcounts(Half_Rsq800t200(bigRs,1),[0:.025:1],'normalization','probability')*100;
thist_hlf_mov_r2 = histcounts(Half_Rsq800t200(bigRs,2),[0:.025:1],'normalization','probability')*100;
thist_hlf_clp_r2 = histcounts(Half_Rsq800t200(bigRs,3),[0:.025:1],'normalization','probability')*100;

thist_frm_r2 = histcounts(Half_RsqFRM(bigRs,1),[0:.025:1],'normalization','probability')*100;
thist_frm_mov_r2 = histcounts(Half_RsqFRM(bigRs,2),[0:.025:1],'normalization','probability')*100;

figure(2010); clf; hold on; 
for ijk = 1:180
    if Half_Rsq800t200(ijk,2) >.2
        plot([1 2],[Half_Rsq800t200(ijk,2) Half_Rsq800t200(ijk,1)],'o-')
    end
end
disp(sum(Half_Rsq800t200(:,2) >.2))
disp(sum(Half_Rsq800t200(:,3) >.2))
set(gca,'xlim',[.5 2.5],'ylim',[0 1])
% saveas(2010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_800t200_Paired_Thres2.pdf','pdf');

figure(7010); hold on;
bar([0.0125:.025:0.9875],thist_hlf_mov_r2,'k')
bar([0.0125:.025:0.9875],thist_hlf_r2,'r')
% bar([0.0125:.025:0.9875],thist_frm_r2,'c')
%     set(gca,'xlim',[-10 10],'ylim',[0 20])
set(gca,'xlim',[-.05 1.05],'ylim',[0 45])
title(['Reg Rsq Threshold .2'])
% saveas(7010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_800t200_halves_Thres2.pdf','pdf');

figure(7015); hold on;
bar([0.0125:.025:0.9875],thist_frm_mov_r2,'k')
bar([0.0125:.025:0.9875],thist_frm_r2,'r')
% bar([0.0125:.025:0.9875],thist_frm_r2,'c')
%     set(gca,'xlim',[-10 10],'ylim',[0 20])
set(gca,'xlim',[-.05 1.05],'ylim',[0 50])
title(['Reg FRM Rsq Threshold .2'])
% saveas(7015,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_FRM_halves_hst_Thres2.pdf','pdf');

figure(2020); clf; hold on;
for ijk = 1:180
    plot([1 2],[Half_Rsq800t200(ijk,2) Half_Rsq800t200(ijk,1)],'o-')
end
set(gca,'xlim',[.5 2.5],'ylim',[0 1])
% saveas(2020,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_800t200_Paired_noThres.pdf','pdf');


figure(3010); clf; hold on; 
for ijk = 1:180
    if Half_Rsq800t200(ijk,2) >.2
        if ijk == 1 || ijk == 94 
            plot([1 2],[Half_RsqFRM(ijk,2) Half_RsqFRM(ijk,1)],'xk-')
        else
            plot([1 2],[Half_RsqFRM(ijk,2) Half_RsqFRM(ijk,1)],'o-')
        end
    end
end
set(gca,'xlim',[.5 2.5],'ylim',[0 1])
% saveas(3010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_Frame_Paired_Thres2.pdf','pdf');

figure(3020); clf; hold on;
for ijk = 1:180
    plot([1 2],[Half_RsqFRM(ijk,2) Half_RsqFRM(ijk,1)],'o-')
end
set(gca,'xlim',[.5 2.5],'ylim',[0 1])
% saveas(3020,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_Frame_Paired_noThres.pdf','pdf');

% figure(5020); clf; hold on;
% bar(-.975:.05:.975,thist_indx_r2,'k')
% set(gca,'xlim',[-.75 .75],'ylim',[0 30])
% title(['Index Threshold .2'])
% saveas(5020,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Indx_800Early_Thres2.pdf','pdf');
% [h,p] = ttest(AF_Indx800t200(bigRs,1))


[~,p,~,t_stats] = ttest(Half_Rsq800t200(bigRs,2),Half_Rsq800t200(bigRs,1))

[~,p,~,t_stats] = ttest(Half_RsqFRM(bigRs,2),Half_RsqFRM(bigRs,1))

mean(Full_Rsq800t200(:,1))
std(Full_Rsq800t200(:,1))

mean(Full_RsqFRM(:,1))
std(Full_RsqFRM(:,1))

%% Late response

bigRs = Half_Rsq800(:,2) >.2;

thist_hlf_r2 = histcounts(Half_Rsq800(bigRs,1),[0:.025:1],'normalization','probability')*100;
thist_hlf_mov_r2 = histcounts(Half_Rsq800(bigRs,2),[0:.025:1],'normalization','probability')*100;
thist_hlf_clp_r2 = histcounts(Half_Rsq800(bigRs,3),[0:.025:1],'normalization','probability')*100;

figure(4010); clf; hold on; 
for ijk = 1:108
    if Half_Rsq800(ijk,2) >.2
        if ijk == 1 || ijk == 94 
            plot([1 2],[Half_Rsq800(ijk,2) Half_Rsq800(ijk,1)],'xk-')
        else
            plot([1 2],[Half_Rsq800(ijk,2) Half_Rsq800(ijk,1)],'o-')
        end
    end
end
disp(sum(Half_Rsq800(:,2) >.2))
disp(sum(Half_Rsq800(:,3) >.2))
set(gca,'xlim',[.5 2.5],'ylim',[0 1])
saveas(4010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_800Late_Paired_Thres2.pdf','pdf');

figure(8010); hold on;
bar([0.0125:.025:0.9875],thist_hlf_mov_r2,'k')
bar([0.0125:.025:0.9875],thist_hlf_r2,'r')
bar([0.0125:.025:0.9875],thist_hlf_clp_r2,'c')
%     set(gca,'xlim',[-10 10],'ylim',[0 20])
set(gca,'xlim',[-.05 1.05],'ylim',[0 15])
title(['Reg Rsq Threshold .2'])
saveas(8010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_800Late_halves_Thres2.pdf','pdf');

% figure(5010); hold on;
% bar(-.975:.05:.975,thist_indx_r2,'k')
% set(gca,'xlim',[-.6 .6],'ylim',[0 30])
% title(['Index Threshold .2'])
% saveas(5010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Indx_800Late_Thres2.pdf','pdf');
% [h,p] = ttest(AF_Indx800(bigRs,1))

[~,p,~,t_stats] = ttest(Half_Rsq800(bigRs,2),Half_Rsq800(bigRs,1))

mean(Full_Rsq800(:,1))
std(Full_Rsq800(:,1))

%% Anova test
anova_time = [ones(size(Half_Rsq800,1)*2,1); ones(size(Half_Rsq800,1)*2,1)*2]; % 1 = early response; 2 = late response
anova_context = [ones(size(Half_Rsq800,1),1); ones(size(Half_Rsq800,1),1)*2; ones(size(Half_Rsq800,1),1); ones(size(Half_Rsq800,1),1)*2]; % 1 = in-contex; 2 = between 
anova_data = [Half_Rsq800t200(:,2); Half_Rsq800t200(:,1);  Half_Rsq800(:,2); Half_Rsq800(:,1)];

[p,atab] = anovan(anova_data,[anova_context anova_time],'varnames',{'Context','Time'},'model','full','sstype', 3)
    
anova_time = [ones(size(Half_Rsq800,1)*2,1); ones(size(Half_Rsq800,1)*2,1)*2]; % 1 = early response; 2 = Frame
anova_context = [ones(size(Half_Rsq800,1),1); ones(size(Half_Rsq800,1),1)*2; ones(size(Half_Rsq800,1),1); ones(size(Half_Rsq800,1),1)*2]; % 1 = in-contex; 2 = between 
anova_data = [Half_Rsq800t200(:,2); Half_Rsq800t200(:,1);  Half_RsqFRM(:,2); Half_RsqFRM(:,1)];

[p,atab] = anovan(anova_data,[anova_context anova_time],'varnames',{'Context','Time'},'model','full','sstype', 3)

%% 250 and 100 
bigRs = Half_Rsq800t200(:,2) >.2;

thist_hlf_r2 = histcounts(Half_Rsq250(bigRs,1),[0:.025:1],'normalization','probability')*100;
thist_hlf_mov_r2 = histcounts(Half_Rsq250(bigRs,2),[0:.025:1],'normalization','probability')*100;
thist_hlf_clp_r2 = histcounts(Half_Rsq250(bigRs,3),[0:.025:1],'normalization','probability')*100;
thist_250_r2 = histcounts(Full_Rsq250(:,1),[0:.025:1],'normalization','probability')*100;

figure(2250); clf; hold on; 
for ijk = 1:108
    if Half_Rsq800t200(ijk,2) >.2
        if ijk == 1 || ijk == 94 
            plot([1 2],[Half_Rsq250(ijk,2) Half_Rsq250(ijk,1)],'xk-')
        else
            plot([1 2],[Half_Rsq250(ijk,2) Half_Rsq250(ijk,1)],'or-')
        end
    end
end
set(gca,'xlim',[.5 2.5],'ylim',[0 1])
saveas(2250,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_250_Paired_Thres2.pdf','pdf');

figure(3250); clf; hold on;
bar([0.0125:.025:0.9875],thist_250_r2,'k')
set(gca,'xlim',[-.05 1.05])
saveas(3250,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_250all.pdf','pdf');

thist_hlf_r2 = histcounts(Half_Rsq100(bigRs,1),[0:.025:1],'normalization','probability')*100;
thist_hlf_mov_r2 = histcounts(Half_Rsq100(bigRs,2),[0:.025:1],'normalization','probability')*100;
thist_hlf_clp_r2 = histcounts(Half_Rsq100(bigRs,3),[0:.025:1],'normalization','probability')*100;
thist_100_r2 = histcounts(Full_Rsq100(:,1),[0:.025:1],'normalization','probability')*100;

figure(2100); clf; hold on; 
for ijk = 1:108
    if Half_Rsq800t200(ijk,2) >.2
        if ijk == 1 || ijk == 94 
            plot([1 2],[Half_Rsq100(ijk,2) Half_Rsq100(ijk,1)],'xk-')
        else
            plot([1 2],[Half_Rsq100(ijk,2) Half_Rsq100(ijk,1)],'or-')
        end
    end
end
set(gca,'xlim',[.5 2.5],'ylim',[0 1])
saveas(2100,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_100_Paired_Thres2.pdf','pdf');

figure(3100); clf; hold on;
bar([0.0125:.025:0.9875],thist_100_r2,'k')
set(gca,'xlim',[-.05 1.05])
saveas(3100,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_100all.pdf','pdf');

%%
Comb_bigR = Half_Rsq800(:,2) >.2;

WindowSlopeIndex_mac(All_names,Comb_bigR,800,1);
% saveas(4801,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_Wind_800_T2.pdf');
saveas(4801,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Comb_Rsq_Wind_800_T2_wClip.pdf');


WindowSlopeIndex_mac(All_names,Comb_bigR,250,25);
% saveas(4825,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Rsq_Wind_250_T2.pdf');

%% Now for Faces! 

if ~exist('Comb_fdprime','var')
%     [~,AF_fdprime] = PlotCatHeatmap(AF_CatData,0);
    load([loaddir filesep 'AF_fdprime.mat'])
    load([loaddir filesep 'AM_fdprime.mat'])
    Comb_fdprime = [AF_fdprime; AM_fdprime];
    clear AF_fdprime AM_fdprime
end
fbottom = Comb_fdprime < .65;
ftop  = Comb_fdprime >= .65;

figure(341); clf;
plot(Comb_fdprime,Half_Rsq800t200(:,1),'.k','MarkerSize',15)
lsline;
set(gca,'xlim',[-6 8],'ylim',[-.05 1],'xtick',-6:2:8,'ytick',0:.2:1)
% saveas(341,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_FaceEarly_Rsq.pdf');


figure(441); clf;
plot(Comb_fdprime,Half_Rsq800(:,1),'.k','MarkerSize',15)
lsline;
set(gca,'xlim',[-6 8],'ylim',[-.05 1],'xtick',-6:2:8,'ytick',0:.2:1)
% saveas(441,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_FaceLate_Rsq.pdf');

figure(678); clf;
thist_dprime = histcounts(Comb_fdprime,[-6:.25:8],'normalization','probability')*100;
bar([-5.875:.25:7.875],thist_dprime,'k')
set(gca,'xlim',[-6 8])
% saveas(678,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_Dprime_hist.pdf');


[AF_fr,AF_frp] = corrcoef(Comb_fdprime,Full_Rsq800(:,1));
disp(['800Late Face Corr=' num2str(AF_fr(1,2)) '; p=' num2str(AF_frp(1,2))])

[AF200_fr,AF200_frp] = corrcoef(Comb_fdprime,Half_Rsq800t200(:,1));
disp(['800Trans Face Corr=' num2str(AF200_fr(1,2)) '; p=' num2str(AF200_frp(1,2))])
close all
%% Scatter plots for figure 2
Plot_AllFR_scatter(movFRs800(65),clpFRs800(65));
saveas(9000,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/AF65_spice07782_late_scatter.pdf');
close(9000)
Plot_AllFR_scatter(movFRs800t200(65),clpFRs800t200(65));
saveas(9000,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/AF65_spice07782_early_scatter.pdf');
close(9000)
Plot_AllFR_scatter(movFRsFRM(65),clpFRsFRM(65));
saveas(9000,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/AF65_spice07782_frame_scatter.pdf');
close(9000)

%% Plot Figure 2 example

Plot_AcrossCond_FR(movFRs800{94},clpFRs800{94},[1:300],948,[],18)
Plot_AcrossCond_FR(movFRs800t200{94},clpFRs800t200{94},[1:300],942,[],18)
Plot_AcrossCond_FR(movFRsFRM{94},clpFRsFRM{94},[1:300],940,[],18)
saveas(7948,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/AF94_SDF_LateResponse.pdf');
saveas(7942,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/AF94_SDF_Transition.pdf');
saveas(7940,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/AF94_SDF_Frame.pdf');

PlotCatData(CatData(94),94,1,0,5);
saveas(3294,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/AF94_CatBars.pdf');


Plot_AllFR_scatter(movFRs800(94),clpFRs800(94));

Plot_AllFR_scatter(movFRs800t200(94),clpFRs800t200(94));
saveas(9000,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/AF94_Scatter_Transition.pdf');
tmpclp = nanmean(clpFRs800t200{94},2);
tmpmov = nanmean(movFRs800t200{94},2);
tmpdif = tmpmov - tmpclp;
tmp_clphist=histcounts(tmpdif,[-21:.5:21],'normalization','probability')*100;
figure(7000);clf; hold on;
bar([-20.75:.5:20.75],tmp_clphist,'k')
saveas(7000,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/AF94_Bar_Transition.pdf');



Plot_AllFR_scatter(movFRsFRM(94),clpFRsFRM(94));
% PlotCatData(CatCellData,fignum,plotmain,plotextra,ylim)
% PlotCatData(AF_CatData(65),65,1,0,5);
% saveas(3265,'/Users/brian/Documents/Analyses/MovieTiling/Figures/AF_cell65_CatBars.pdf');

% Plot_AcrossCond_FR(movFRs800t200{65},clpFRs800t200{65},[1:300],200,[],10)
% saveas(7200,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs14/AF_Cell65_SDF_Transition.pdf','pdf');
% 
% Plot_AcrossCond_FR(movFRs800{65},clpFRs800{65},[1:300],800,[],10)
% saveas(7800,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs14/AF_Cell65_SDF_LateResponse.pdf','pdf');
% 
% Plot_AcrossCond_FR(AF_movFRsFRM{65},AF_clpFRsFRM{65},[1:300],0,[],10)
% saveas(7200,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs14/AF_Cell65_SDF_Frame.pdf','pdf');

%% Plot Spikes for Clips 
%% 11/22/17
% binS    = '-350';
% binL    = '800';
% ii = 94;
% PlotClipData211227(AF_05Clip800Data(ii),AF_Movie05Data(ii),ii,...
%                    '1','1','301','0',binS,binL,'0',ii,'full');
%     PlotClipData_slurm(ClipCellData,MovieCellData,fignum,...
%         plotmain,plotextra,whichclip,plotmovsdf,bin_strt,bin_lngth,splthlv,cellnum,savesuffix)

binS    = -350;
binL    = 800;
TPS = 12:18;
PlotIndClipTC_211227(AF_05Clip800Data(ii),AF_Movie05Data(ii),AF_05FrameData(ii),ii,...
    TPS,ii,20,binS,binL,'set1')

TPS = 60:68;
PlotIndClipTC_211227(AF_05Clip800Data(ii),AF_Movie05Data(ii),AF_05FrameData(ii),ii,...
    TPS,ii,20,binS,binL,'set2')

TPS = 105:113;
PlotIndClipTC_211227(AF_05Clip800Data(ii),AF_Movie05Data(ii),AF_05FrameData(ii),ii,...
    TPS,ii,20,binS,binL,'set3')

TPS = 150:155;
PlotIndClipTC_211227(AF_05Clip800Data(ii),AF_Movie05Data(ii),AF_05FrameData(ii),ii,...
    TPS,ii,20,binS,binL,'set4')

TPS = 150:155;
PlotIndClipTC_211227(AF_05Clip800Data(ii),AF_Movie05Data(ii),AF_05FrameData(ii),ii,...
    TPS,ii,20,binS,binL,'set5')

TPS = 178:185;
PlotIndClipTC_211227(AF_05Clip800Data(ii),AF_Movie05Data(ii),AF_05FrameData(ii),ii,...
    TPS,ii,20,binS,binL,'set6')

TPS = [260 261 262 275 276 277];
PlotIndClipTC_211227(AF_05Clip800Data(ii),AF_Movie05Data(ii),AF_05FrameData(ii),ii,...
    TPS,ii,20,binS,binL,'set7')
%     PlotIndClipTC_slurm(Clipfile,Moviefile,Framefile,cellnum,
%         TPS,fignum,sdf_ylim,bin_strt,bin_lngth,filesuf)
%%
TPS     = '240:300';
cellnum = '1'; %'94'; 
plotmov = '1';
ylim    = '60';
binS    = '-350';
binL    = '800';
PlotFullClipTC4_slurm(AF_05Clip800Data,AF_Movie05Data,cellnum,TPS,cellnum,plotmov,ylim,binS,binL,'min5');
%% Category Data
CatData = struct([]);
tic;

disp('AF')
% curdir = cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/torset01_data');
curdir = cd('/Volumes/BER_MovieTiling/MovieTiling_results/torset01_data');
load('Tor01_CatData');
CatData = [CatData Tor01_CatData];
clear Tor01*
cd(curdir)
toc/60

% cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/rhomset01_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/rhomset01_data');
load('Rhom01_CatData');
CatData = [CatData Rhom01_CatData];
clear Rhom01* 
cd(curdir)
toc/60

% cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/spiceset06_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/spiceset06_data');
load('Spice06_CatData');
CatData = [CatData Spice06_CatData];
clear Spice06* 
cd(curdir)
toc/60

% cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/spiceset07_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/spiceset07_data');
load('Spice07_CatData');
CatData = [CatData Spice07_CatData];
clear Spice07* 
cd(curdir)
toc/60

rmcells = [3 24 26 38 49 52 57];
CatData(rmcells) = [];

clear rmcells

AM_CatData = struct([]);
% load AM data
disp('AM')
tic;
% curdir = cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/matset03_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/matset03_data');
load('Mat03_CatData');
AM_CatData = [AM_CatData Mat03_CatData];
clear Mat03* 
cd(curdir)

% curdir = cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/dangoset01_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/dangoset01_data');
load('Dango01_CatData');
AM_CatData = [AM_CatData Dango01_CatData];
clear Dango01* 
cd(curdir)
toc/60

rmcells = [29 38 72 73 75 76 77 78  80 83 86 87];
AM_CatData(rmcells) = [];

rmcells2 = [36 38 47 67 73 75];
AM_CatData(rmcells2) = [];

CatData = [CatData AM_CatData];

% save('/Users/brian/Library/Mobile Documents/com~apple~CloudDocs/Documents/Analyses/MovieTiling/data/Workspace/CatData.mat','CatData','-v7.3');

clear curdir ans rmcells rmcells2 AM_CatData
%%
PlotCatHeatmap(CatData,0);
saveas(890,'/Users/brian/Google Drive/Writing/MovieTiling Paper/RoughFigs18/Comb_FaceHeatMap.pdf');
ftop  = Comb_fdprime >= .65;
sum(ftop)

%% Latencies
Cat_latencies = zeros(size(All_names));
for ijk = 1:length(All_names)
    [~,Cat_latencies(ijk)] = CatRespLatency_220726(CatData,ijk,1,1,30,0);
    pause(.5)
end

nanmean(Cat_latencies)
nanstd(Cat_latencies)

thist_latency = histcounts(Cat_latencies(~isnan(Cat_latencies)),[0:10:250],'normalization','probability')*100;
figure(2020); clf;
bar([5:10:245],thist_latency,'k')
set(gca,'xlim',[0 260],'ylim',[0 25])
saveas(2020,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Comb_Cat_Latencies.pdf','pdf');


%% Category Data
All_FrameData = struct([]);
tic;

disp('AF')
% curdir = cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/torset01_data');
curdir = cd('/Volumes/BER_MovieTiling/MovieTiling_results/torset01_data');
load('Tor01_05FrameData');
All_FrameData = [All_FrameData Tor01_05FrameData];
clear Tor01*
cd(curdir)
toc/60

% cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/rhomset01_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/rhomset01_data');
load('Rhom01_05FrameData');
All_FrameData = [All_FrameData Rhom01_05FrameData];
clear Rhom01* 
cd(curdir)
toc/60

% cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/spiceset06_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/spiceset06_data');
load('Spice06_05FrameData');
All_FrameData = [All_FrameData Spice06_05FrameData];
clear Spice06* 
cd(curdir)
toc/60

% cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/spiceset07_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/spiceset07_data');
load('Spice07_05FrameData');
All_FrameData = [All_FrameData Spice07_05FrameData];
clear Spice07* 
cd(curdir)
toc/60

rmcells = [3 24 26 38 49 52 57];
All_FrameData(rmcells) = [];

clear rmcells

AM_05FrameData = struct([]);
% load AM data
disp('AM')
tic;
% curdir = cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/matset03_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/matset03_data');
load('Mat03_05FrameData');
AM_05FrameData = [AM_05FrameData Mat03_05FrameData];
clear Mat03* 
cd(curdir)

% curdir = cd('/data/NIF/procdata/russbe/Physiology/MovieTiling_results/dangoset01_data');
cd('/Volumes/BER_MovieTiling/MovieTiling_results/dangoset01_data');
load('Dango01_05FrameData');
AM_05FrameData = [AM_05FrameData Dango01_05FrameData];
clear Dango01* 
cd(curdir)
toc/60

rmcells = [29 38 72 73 75 76 77 78  80 83 86 87];
AM_05FrameData(rmcells) = [];

rmcells2 = [36 38 47 67 73 75];
AM_05FrameData(rmcells2) = [];

All_FrameData = [All_FrameData AM_05FrameData];


clear curdir ans rmcells rmcells2 AM_05FrameData
%% Also for the clips
[All_FrameLatencies] = ClipLatency_mac210628(All_FrameData,0);
mn_lat = nanmean(abs(All_FrameLatencies),2);
std_lat = nanstd(abs(All_FrameLatencies'))';
lats = [mn_lat std_lat];
mean(lats)

thist_latency = histcounts(lats(:,1),[0:10:180],'normalization','probability')*100;
figure(1010); clf;
bar([5:10:175],thist_latency,'k')
set(gca,'xlim',[0 200],'ylim',[0 40])
saveas(1010,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Comb_Frame_Latencies.pdf','pdf');

%% Get Slopes (again) 
if ~exist('MovieBase','var')
    load([loaddir filesep 'AF_movbase.mat'])
    load([loaddir filesep 'AM_movbase.mat'])
    MovieBase = [AF_movbase; AM_movbase];
    clear AF_movbase AM_movbase
end
lims = [-.25 1.1 15 12 0];

[Slope_800_late] = SlopeStats_half(movFRs800,clpFRs800,MovieBase,[],'800',22,lims);
thist_slope = histcounts(Slope_800_late,[-.2:.1:2],'normalization','probability')*100;
figure(3030); clf;
bar([-.15:.1:1.95],thist_slope,'k')
set(gca,'xlim',[-.5 2],'ylim',[0 25])
saveas(3030,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Comb_Slope_800Late.pdf','pdf');

lims = [-.25 1.1 25 30 0];
[Slope_800_tran] = SlopeStats_half(movFRs800t200,clpFRs800t200,MovieBase,[],'Trans',2,lims);
thist_slope = histcounts(Slope_800_tran,[-.2:.1:2],'normalization','probability')*100;
figure(4040); clf;
bar([-.15:.1:1.95],thist_slope,'k')
set(gca,'xlim',[-.5 2],'ylim',[0 25])
saveas(3030,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Comb_Slope_800Trans.pdf','pdf');
clear thist*

%% Check for consistency of response based on faciness 
[MovLate_face_cor,MovLate_face_p] = corrcoef(Face_dprime,Half_Rsq800(:,2));
disp(['800Late Movie Corr=' num2str(MovLate_face_cor(1,2)) '; p=' num2str(MovLate_face_p(1,2))])

[MovTrans_face_cor,MovTrans_face_p] = corrcoef(Face_dprime,Half_Rsq800t200(:,2));
disp(['800Trans Movie Corr=' num2str(MovTrans_face_cor(1,2)) '; p=' num2str(MovTrans_face_p(1,2))])

[ClipLate_face_cor,ClipLate_face_p] = corrcoef(Face_dprime,Half_Rsq800(:,3));
disp(['800Late Clip Corr=' num2str(ClipLate_face_cor(1,2)) '; p=' num2str(ClipLate_face_p(1,2))])

[ClipTrans_face_cor,ClipTrans_face_p] = corrcoef(Face_dprime,Half_Rsq800t200(:,3));
disp(['800Trans Clip Corr=' num2str(ClipTrans_face_cor(1,2)) '; p=' num2str(ClipTrans_face_p(1,2))])

figure(341); clf; hold on;
plot(Face_dprime,Half_Rsq800t200(:,2),'.r','MarkerSize',15)
lsline;
plot(Face_dprime,Half_Rsq800t200(:,3),'.b','MarkerSize',15)
lsline;
set(gca,'xlim',[-6 8],'ylim',[-.05 1],'xtick',-6:2:8,'ytick',0:.2:1)
saveas(341,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Comb_Face_Trans_MovClipSep.pdf');



figure(441); clf; hold on;
plot(Face_dprime,Half_Rsq800(:,2),'.r','MarkerSize',15)
lsline;
plot(Face_dprime,Half_Rsq800(:,3),'.b','MarkerSize',15)
lsline;
set(gca,'xlim',[-6 8],'ylim',[-.05 1],'xtick',-6:2:8,'ytick',0:.2:1)
saveas(441,'/Users/brian/Google Drive/Writing/MovieTiling Paper/Neuron Sub/Revision/Comb_Face_Late_MovClipSep.pdf');


