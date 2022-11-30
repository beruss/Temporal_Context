function [resp_lat,abslat] = CatRespLatency_220726(Catfile,cellnum,plotmain,plotextra,ylim,savedata)

% CellData = Comb_Clips100Data(1);

CatCellData  = Catfile(cellnum);

if ~exist('plotmain','var') || isempty(plotmain)
    plotmain = 0;
end

if ~exist('plotextra','var') || isempty(plotextra)
    plotextra = 0;
end

if ~exist('ylim','var') || isempty(ylim)
    ylim = 30;
end

if ~exist('savedata','var') || isempty(savedata)
    savedata = 0;
end

%%
bins = -50:1:250; %msec bins for creating sdf 

neuro_trials = CatCellData.cat_spikes;
types =CatCellData.cat_stim;

if plotmain && plotextra
    f2 = figure(900+cellnum); clf; hold on;
end

cat_means = zeros(1,6)*NaN;
cat_stde = zeros(1,6)*NaN;
    
resp_lat = cell(3,6); % response latencies (latency, direction, significant response)

for iii = 1:6
        %     cat_set = find(TDTdata2(ii).stimtimes(1,:)==iii);
    cat_set = find([neuro_trials{:,1}]==iii);

    spike_bins = zeros(length(cat_set),length(neuro_trials{cat_set(1),4}));
    base_bins = zeros(length(cat_set),76); % baseline currently -50:25
    resp_bins = zeros(length(cat_set),225); % baseline currently 26:250
    for jjj = 1:length(cat_set)
        spike_bins(jjj,:) = neuro_trials{cat_set(jjj),4};
        base_bins(jjj,:) = spike_bins(jjj,1:76);
        resp_bins(jjj,:) = spike_bins(jjj,77:301);
    end

    mn_spike = mean(spike_bins,1)*1000;
    cat_sdf = makeSDF(mn_spike,3);

    mn_base = mean(base_bins,2)*1000;
    mn_resp = mean(resp_bins,2)*1000;
        
    if plotmain == 1 && plotextra == 1
        subplot(2,3,iii); cla; hold on;

        to_plot_list = randperm(length(cat_set));
        to_plot_trials = cat_set(to_plot_list);
        for jjj = 1:50
            x = neuro_trials{to_plot_trials(jjj),3};
            y = ones(size(x))*jjj;
%                     h1 = plot(x,y,'o'); set(h1,'MarkerEdgeColor','none','MarkerFaceColor','k');
            for jk = 1:length(x)
                line([x(jk) x(jk)],[y(jk) y(jk)+1],'Color',[.2 .2 .2],'LineWidth',1.5);
            end
            set(gca,'XLim',[-50 250],'YLim',[0 55])
            clear x y h1
        end
        plot(bins,cat_sdf,'r-','LineWidth',2.5);
        xlabel(types{iii},'FontSize',10);
        line([0 0],[0 5],'Color',[0 0 0],'LineWidth',2);
        line([0 0],[0 10],'Color',[0 0 0],'LineWidth',2);
        line([50 50],[0 2],'Color',[0 0 0],'LineWidth',2);
        line([100 100],[0 2],'Color',[0 0 0],'LineWidth',2);
        line([150 150],[0 2],'Color',[0 0 0],'LineWidth',2);
        line([200 200],[0 2],'Color',[0 0 0],'LineWidth',2);
    end
    
    twinds = [1:25:251];
    catP = zeros(length(twinds)+1,2)*NaN;
    catP(1:length(twinds),1) = twinds';
    for hij = 1:length(twinds)
        si = twinds(hij);
        se = si+50;
        [~,catP(hij,2)] = ttest2(mn_spike(1:51),mn_spike(si:se));
    end
    [~,catP(hij+1,2)] = ttest2(mn_base,mn_resp);
    
    Pstdbase = mean(mn_spike(1:50)) + std(mn_spike(1:50))*2; % calculate 3 STD above mean 
    Prespstd = [mn_spike(51:end) >= Pstdbase]'; % logical of positive responses
    Nstdbase = mean(mn_spike(1:50)) - std(mn_spike(1:50))*2; % calculate 3 STD below mean 
    Nrespstd = [mn_spike(51:end) <= Nstdbase]'; % logical of negative responses
    for ijk = 15:length(Prespstd)-20
        if Prespstd(ijk)==1 && sum(Prespstd(ijk:ijk+20)) >= 20*.80 % check the response for 20msecs 
            resp_lat{1,iii} = ijk;
            resp_lat{2,iii} = 1;
            resp_lat{3,iii} = catP;
            
            line([ijk ijk],[0 12],'Color',[0 .3 0],'LineWidth',.75);
            break
        elseif sum(Nrespstd(ijk:ijk+20)) >= 20*.80
            resp_lat{1,iii} = ijk;
            resp_lat{2,iii} = -1;
            resp_lat{3,iii} = catP;
            
            line([ijk ijk],[0 12],'Color',[.3 0 0],'LineWidth',.75);
            break
        end
    end

%     resp_lat;
%     
%     switch iii
%         case {1, 2}
%             face_data = [face_data; mn_resp];
%         case {3, 4, 5, 6}
%             other_data = [other_data; mn_resp];
%     end
%     
%     anova_data = [anova_data; mn_base; mn_resp];
%     anova_group1 = [anova_group1; (ones(length(mn_base),1)*iii); (ones(length(mn_resp),1)*iii)];
%     anova_group2 = [anova_group2; (ones(length(mn_base),1)*1); (ones(length(mn_resp),1)*2)];

    
    cat_means(iii) = nanmean(mn_resp);
    cat_stde(iii) = nanstd(mn_resp)/sqrt(length(mn_resp));

    clear cat_sdf mn_spike spike_bins cat_set mn_base mn_resp base_bins resp_bins
end


abslat = nanmean([resp_lat{1,:}]);

% [~,faceP,~,faceStats] = ttest2(face_data,other_data);
% faceT = faceStats.tstat;

if plotmain == 1 && plotextra == 1
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
        'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

    celltitle = ['\bf Cell ' num2str(cellnum)];
    text(0.5, 1,celltitle,'HorizontalAlignment','center','VerticalAlignment',...
        'top','FontWeight','bold','FontSize',14)
end
