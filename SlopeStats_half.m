function [outslope,inout_rsq_halves,inout_rsq,ClipSet,inout_pstat,inout_index,inout_clipIDX] = ...
    SlopeStats_half(movFRs,clpFRs,movBase,ClipSet,clen,FigNum,lims)

if nargin < 6
    lims = [-2 2 40 0];
end

if length(lims) < 3
    lims(3:5) = [20 20 0];
elseif length(lims) < 5
    lims(5) = 0;
end

if isempty(ClipSet)
    ClipSet = zeros(size(movFRs,1),300)*NaN;
end
%% Stats800ms
inout_slope = zeros(size(movFRs,1),2)*NaN;
inout_index = zeros(size(movFRs,1),4)*NaN;
inout_clipIDX = zeros(size(movFRs,1),300)*NaN;
% inout_index = zeros(size(movFRs,1),7)*NaN;
inout_test = zeros(size(movFRs,1),1)*NaN;
inout_Stat = zeros(size(movFRs,1),1)*NaN;
inout_diff = zeros(size(movFRs,1),1)*NaN;
inout_rsq  = zeros(size(movFRs,1),1)*NaN;
inout_pstat  = zeros(size(movFRs,1),2)*NaN;

inout_rsq_halves = zeros(size(movFRs,1),6)*NaN;
inout_slope_halves = zeros(size(movFRs,1),6)*NaN;

% inout_dist = zeros(size(movFRs,1),7)*NaN;
% indx_ttest = zeros(size(movFRs,1),4)*NaN;
store_means = zeros(size(movFRs,1),2)*NaN;
% numhigh = zeros(size(movFRs,1),1)*NaN;

% anovagroups = zeros(600*size(movFRs,1),3)*NaN;
% anovadata = zeros(600*size(movFRs,1),1)*NaN;
% 
% movFR_mn = zeros(size(movFRs,1),1)*NaN;
% clpFR_mn = zeros(size(movFRs,1),1)*NaN;    
% 
% unity = [0,0,0; 10000,10000,0];
% dist_a = unity(1,:) - unity(2,:);
for ijk = 1:length(movFRs)
    tmpCLP = clpFRs{ijk};
    tmpMOV = movFRs{ijk};
%     numresp = [(sum(tmpMOV' ~= 0)/size(tmpMOV,2)); (sum(tmpCLP' ~= 0)/size(tmpCLP,2))];
    
    crnt_clipFR_mn = mean(tmpCLP,2);
    crnt_movFR_mn  = mean(tmpMOV,2);
    
%     crnt_movFR_std = std(crnt_movFR_mn);
%     crnt_clpFR_std = std(crnt_clipFR_mn);
    mn_crnt_movFR = mean(crnt_movFR_mn);
    mn_crnt_clpFR = mean(crnt_clipFR_mn);
    store_means(ijk,:) = [mn_crnt_movFR mn_crnt_clpFR];
    
%     resp_mov = crnt_movFR_mn > (mn_crnt_movFR-(crnt_movFR_std*.5));
%     resp_clp = crnt_clipFR_mn > (mn_crnt_clpFR-(crnt_clpFR_std*.5));
%     resp_bth = resp_mov == 1 & resp_clp == 1;
%     resp_ethr = resp_mov == 1 | resp_clp == 1;
%     
%     if isnan(ClipSet(ijk,1))
%         if isempty(movBase)
%             midd_mov = crnt_movFR_mn > (mn_crnt_movFR+(crnt_movFR_std*.5));
%             high_mov = crnt_movFR_mn > (mn_crnt_movFR+(crnt_movFR_std*2));
%         else 
%             midd_mov = crnt_movFR_mn > (movBase(ijk,1)+(movBase(ijk,2)*.5));
%             high_mov = crnt_movFR_mn > (movBase(ijk,1)+(movBase(ijk,2)*2));
%         end
%         ClipSet(ijk,:) = high_mov;
%     elseif size(ClipSet,1) > 1
%         if isempty(movBase)
%             midd_mov = crnt_movFR_mn > (mn_crnt_movFR+(crnt_movFR_std*.5));
%             high_mov = logical(ClipSet(ijk,:));
%         else 
%             midd_mov = crnt_movFR_mn > (movBase(ijk,1)+(movBase(ijk,2)*.5));
%             high_mov = logical(ClipSet(ijk,:));
%         end
%     else 
%         if isempty(movBase)
%             midd_mov = crnt_movFR_mn > (mn_crnt_movFR+(crnt_movFR_std*.5));
%             high_mov = ClipSet(1,:);
%         else 
%             midd_mov = crnt_movFR_mn > (movBase(ijk,1)+(movBase(ijk,2)*.5));
%             high_mov = ClipSet(1,:);
%         end
%     end
%     numhigh(ijk) = sum(high_mov);
    
    [inout_slope(ijk,:),~,~,~,tmp_reg_stats] = regress(crnt_clipFR_mn,[ones(size(crnt_movFR_mn)) crnt_movFR_mn]);
    inout_rsq(ijk,:) = tmp_reg_stats(1);
    
    perm_number = 100;
    slope_perm = zeros(perm_number,1)*nan;
    rsq_perm = zeros(perm_number,1)*nan;
    pv_perm_clp2mov = zeros(perm_number,1)*nan;
    
    slope_perm_clp = zeros(perm_number,1)*nan;
    rsq_perm_clp = zeros(perm_number,1)*nan;
    
    slope_perm_mov = zeros(perm_number,1)*nan;
    rsq_perm_mov = zeros(perm_number,1)*nan;
    pv_perm_mov = zeros(perm_number,1)*nan;
    
    for iii = 1:perm_number
        hlf_clp_all = randperm(size(tmpCLP,2));
        hlf_clp1 = hlf_clp_all(1:ceil(length(hlf_clp_all)/2));
        hlf_clp2 = hlf_clp_all(ceil(length(hlf_clp_all)/2)+1:end);
        clp1_mn = mean(tmpCLP(:,hlf_clp1),2);
        clp2_mn = mean(tmpCLP(:,hlf_clp2),2);
        
        hlf_mov_all = randperm(size(tmpMOV,2));
        hlf_mov1 = hlf_mov_all(1:ceil(length(hlf_mov_all)/2));
        hlf_mov2 = hlf_mov_all(ceil(length(hlf_mov_all)/2)+1:end);
        mov1_mn = mean(tmpMOV(:,hlf_mov1),2);
        mov2_mn = mean(tmpMOV(:,hlf_mov2),2);
        [tmp_slope,~,~,~,tmp_reg_stats] = regress(clp1_mn,[ones(size(mov1_mn)) mov1_mn]);
        slope_perm(iii,:) = tmp_slope(2);
        rsq_perm(iii,:) = tmp_reg_stats(1);
        pv_perm_clp2mov(iii,:) = tmp_reg_stats(3);
        
        [tmp_slope,~,~,~,tmp_reg_stats] = regress(clp2_mn,[ones(size(clp1_mn)) clp1_mn]);
        slope_perm_clp(iii,:) = tmp_slope(2);
        rsq_perm_clp(iii,:) = tmp_reg_stats(1);
        
        [tmp_slope,~,~,~,tmp_reg_stats] = regress(mov2_mn,[ones(size(mov1_mn)) mov1_mn]);
        slope_perm_mov(iii,:) = tmp_slope(2);
        rsq_perm_mov(iii,:) = tmp_reg_stats(1);
        pv_perm_mov(iii,:) = tmp_reg_stats(3);
    end
    
    inout_slope_halves(ijk,:) = [mean(slope_perm) mean(slope_perm_mov) mean(slope_perm_clp) std(slope_perm) std(slope_perm_mov) std(slope_perm_clp)];
    inout_rsq_halves(ijk,:) = [mean(rsq_perm) mean(rsq_perm_mov) mean(rsq_perm_clp) std(rsq_perm) std(rsq_perm_mov) std(rsq_perm_clp)];
    inout_pstat(ijk,:) = [mean(pv_perm_clp2mov) std(pv_perm_mov)];
    
    tmpindex = ((crnt_movFR_mn - crnt_clipFR_mn) ./ (crnt_movFR_mn + crnt_clipFR_mn));
    [~,tmp_indx_test,~,tmpstat] = ttest(tmpindex,0);
    inout_index(ijk,:) = [nanmean(tmpindex) nanstd(tmpindex) tmpstat.tstat tmp_indx_test];
    inout_clipIDX(ijk,:) = tmpindex;
%     tmpindex(isnan(tmpindex)) = 0;
%     
%     inout_index(ijk,2) = mean(tmpindex(resp_ethr));
%     inout_index(ijk,3) = mean(tmpindex(resp_bth));
%     inout_index(ijk,4) = mean(tmpindex(resp_mov));
%     inout_index(ijk,5) = mean(tmpindex(resp_clp));
%     inout_index(ijk,6) = mean(tmpindex(midd_mov));
% 	inout_index(ijk,7) = mean(tmpindex(high_mov));
%     [~,tmpp1,~,tmpstat1] = ttest(tmpindex(midd_mov));
%     [~,tmpp2,~,tmpstat2] = ttest(tmpindex(high_mov));
%     indx_ttest(ijk,:) = [tmpp1, tmpstat1.tstat, tmpp2, tmpstat2.tstat];
    
%     mov_clp_pnts = [crnt_movFR_mn, crnt_clipFR_mn, zeros(size(crnt_movFR_mn))];
%     cur_dists = zeros(size(crnt_movFR_mn))*NaN;
%     for jjj = 1:length(crnt_movFR_mn)
%         dist_b = mov_clp_pnts(jjj,:) - unity(2,:);
%         cur_dists(jjj) = norm(cross(dist_a,dist_b)) / norm(dist_a);
%         if crnt_movFR_mn(jjj) > crnt_clipFR_mn(jjj)
%             cur_dists(jjj) = cur_dists(jjj)*-1;
%         end
%     end
%     inout_dist(ijk,1) = mean(cur_dists);
%     inout_dist(ijk,2) = mean(cur_dists(resp_ethr));
%     inout_dist(ijk,3) = mean(cur_dists(resp_bth));
%     inout_dist(ijk,4) = mean(cur_dists(resp_mov));
%     inout_dist(ijk,5) = mean(cur_dists(resp_clp));
%     inout_dist(ijk,6) = mean(cur_dists(midd_mov));
%     inout_dist(ijk,7) = mean(cur_dists(high_mov));
%     
%     inout_diff(ijk) = mean(crnt_movFR_mn - crnt_clipFR_mn);
    [~,inout_test(ijk),~,tmpstat] = ttest(crnt_movFR_mn,crnt_clipFR_mn);
    inout_Stat(ijk) = tmpstat.tstat;
    
%     anovagroups(1+(600*(ijk-1)):600*ijk,1) = ones(600,1)*ijk;
%     anovagroups(1+(600*(ijk-1)):600*ijk,2) = [ones(300,1); ones(300,1)*2];
%     anovagroups(1+(600*(ijk-1)):600*ijk,3) = [1:300 1:300]';
%     
%     anovadata(1+(600*(ijk-1)):600*ijk,1) = [crnt_movFR_mn; crnt_clipFR_mn];

    clear tmpCLP tmpMOV 
end

% disp(mean(numhigh));

alpha = .05/length(movFRs);

if FigNum > 0
%     figure(4100 + FigNum); 
%     thist = histc(inout_Stat,[-30:1:30]);
%     bar([-29.5:1:29.5],thist(1:end-1))
%     title(['Stats ' clen])
% 
%     figure(4200 + FigNum); 
%     thist = histc(inout_index,[-1:.05:1]);
%     bar([-.975:.05:.975],thist(1:end-1))
%     title(['Index ' clen])
% 

    disp(['slope range =' num2str(min(inout_slope(:,2))) ' ' num2str(max(inout_slope(:,2)))])
    figure(4300 + FigNum); 
    if lims(5)< 1
        thist = histcounts(inout_slope(:,2),[-2:.05:2],'normalization','probability')*100;
    else
        thist = (histcounts(inout_slope(:,2),[-2:.05:2],'normalization','count')/lims(5))*100;
    end
    bar([-1.975:.05:1.975],thist)
    set(gca,'xlim',lims(1:2),'ylim',[0 lims(3)])
    title(['Slope ' clen])
    
    figure(5300 + FigNum); 
    if lims(5)< 1
        thist = histcounts(inout_index(:,1),[-1:.05:1],'normalization','probability')*100;
    else
        thist = (histcounts(inout_index(:,1),[-1:.05:1],'normalization','count')/lims(5))*100;
    end
    bar([-.975:.05:.975],thist)
    set(gca,'xlim',[-.7 .7],'ylim',[0 30])
    title(['Index ' clen])

%     disp(['indx range =' num2str(min(inout_index(:,7))) ' ' num2str(max(inout_index(:,7)))])
%     figure(5500 + FigNum); 
%     if lims(5)< 1
%         thist = histcounts(inout_index(:,7),[-1:.05:1],'normalization','probability')*100;
%     else
%         thist = (histcounts(inout_index(:,7),[-1:.05:1],'normalization','count')/lims(5))*100;
%     end
%     bar([-.975:.05:.975],thist)
%     set(gca,'xlim',[-.5 .7],'ylim',[0 25])
%     set(gca,'xlim',[-.5 1],'ylim',[0 25])
%     title(['High Mov Index ' clen])
    
%     figure(6300 + FigNum); 
%     if lims(5)< 1
%         thist = histcounts(inout_dist(:,1),[-10:.25:10],'normalization','probability')*100;
%     else
%         thist = (histcounts(inout_dist(:,1),[-10:.25:10],'normalization','count')/lims(5))*100;
%     end
%     bar([-9.875:.25:9.875],thist)
%     set(gca,'xlim',[-10 10],'ylim',[0 30])
%     title(['Distance ' clen])
%     
%     figure(6400 + FigNum); 
%     if lims(5)< 1
%         thist = histcounts(inout_dist(:,3),[-10:.25:10],'normalization','probability')*100;
%     else
%         thist = (histcounts(inout_dist(:,3),[-10:.25:10],'normalization','count')/lims(5))*100;
%     end
%     bar([-9.875:.25:9.875],thist)
%     set(gca,'xlim',[-10 10],'ylim',[0 30])
%     title(['Distance ' clen])
%     disp(['dist range =' num2str(min(inout_dist(:,7))) ' ' num2str(max(inout_dist(:,7)))])
%     figure(6600 + FigNum); 
%     if lims(5)< 1
%         thist = histcounts(inout_dist(:,7),[-20:.25:20],'normalization','probability')*100;
%     else
%         thist = (histcounts(inout_dist(:,7),[-20:.25:20],'normalization','count')/lims(5))*100;
%     end
%     bar([-19.875:.25:19.875],thist)
%     set(gca,'xlim',[-10 10],'ylim',[0 20])
%     set(gca,'xlim',[-20 5],'ylim',[0 15])
%     title(['High Mov Distance ' clen])
    
    
    if lims(5)< 1
        thist = histcounts(inout_rsq(:,1),[0:.025:1],'normalization','probability')*100;
        thist_hlf = histcounts(inout_rsq_halves(:,1),[0:.025:1],'normalization','probability')*100;
        thist_hlf_mov = histcounts(inout_rsq_halves(:,2),[0:.025:1],'normalization','probability')*100;
        thist_hlf_clp = histcounts(inout_rsq_halves(:,3),[0:.025:1],'normalization','probability')*100;
        
        thist_indx = histcounts(inout_index(:,1),[-.25:.05:1],'normalization','probability')*100;
    else
        thist = (histcounts(inout_rsq(:,1),[0:.025:1],'normalization','count')/lims(5))*100;
        thist_hlf = (histcounts(inout_rsq_halves(:,1),[0:.025:1],'normalization','count')/lims(5))*100;
        thist_hlf_mov = (histcounts(inout_rsq_halves(:,2),[0:.025:1],'normalization','count')/lims(5))*100;
        thist_hlf_clp = (histcounts(inout_rsq_halves(:,3),[0:.025:1],'normalization','count')/lims(5))*100;
        
        thist_indx = (histcounts(inout_index(:,1),[-.25:.05:1],'normalization','count')/lims(5))*100;
    end
    figure(7500 + FigNum); 
    bar([0.0125:.025:0.9875],thist)
%     set(gca,'xlim',[-10 10],'ylim',[0 20])
    set(gca,'xlim',[-.05 1.05],'ylim',[0 lims(4)])
    title(['Reg Rsq ' clen])
    
    figure(7600 + FigNum);hold on;
    bar([0.0125:.025:0.9875],thist_hlf,'k')
    bar([0.0125:.025:0.9875],thist_hlf_mov,'b')
    bar([0.0125:.025:0.9875],thist_hlf_clp,'y')
%     set(gca,'xlim',[-10 10],'ylim',[0 20])
    set(gca,'xlim',[-.05 1.05],'ylim',[0 lims(4)])
    title(['Reg Rsq Halves ' clen])
    
    mn_r2_halfs = mean(inout_rsq_halves);
    std_r2_halfs = std(inout_rsq_halves)/sqrt(length(inout_rsq_halves(:,1)));
    figure(7700+FigNum); hold on;
    bar([.5 1.5 2.5],mn_r2_halfs(1:3),'k')
    errorbar([.5 1.5 2.5],mn_r2_halfs(1:3),std_r2_halfs(1:3),'.k')
%     plotSpread(inout_rsq_halves(:,1:3),'xValues',[.5 1.5 2.5],'spreadWidth',.75)
    set(gca,'xlim',[-0 3],'ylim',[0 1])
    title(['Reg Rsq Halves ' clen])
    
    
%     figure(7700 + FigNum); 
%     bar([-.225:.05:0.975],thist_indx)
% %     set(gca,'xlim',[-10 10],'ylim',[0 20])
%     if max(thist_indx) < 25
%         set(gca,'xlim',[-.25 1.05],'ylim',[0 25])
%     else
%         set(gca,'xlim',[-.25 1.05],'ylim',[0 40])
%     end
%     title(['Index All ' clen])
end

disp([ clen ' :' num2str(sum(inout_test<alpha)) ' (p = ' ...
    num2str(binopdf(sum(inout_test<alpha),length(movFRs),.5)) ')']);

% sig_diff = inout_diff(inout_test<alpha);
% 
% disp(['Within larger = ' num2str(sum(sig_diff>0)) ';'])

disp([ clen ' :' num2str(sum(inout_pstat(:,1)>alpha)) ' (p = ' ...
    num2str(binopdf(sum(inout_pstat(:,1)>alpha),length(movFRs),.5)) ')']);

% disp(['Index .5:' num2str(sum(indx_ttest(:,1)<alpha)) ' (p = ' ...
%     num2str(binopdf(sum(indx_ttest(:,1)<alpha),length(movFRs),.5)) ')']);
% disp(['Index 2:' num2str(sum(indx_ttest(:,3)<alpha)) ' (p = ' ...
%     num2str(binopdf(sum(indx_ttest(:,3)<alpha),length(movFRs),.5)) ')']);

outslope = inout_slope(:,2);

clear thist crnt* anovag* anovad* 