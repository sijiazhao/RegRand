clear;close all;
% addpath('./Functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
message_sigtime = {};

pct = 0;
spx = 7;
spy = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
message_sigtime = [message_sigtime; '************* Expt1 REG10-RAND20 ***************'];

yrange = [-.2 .8];
tw_stats = [0 2];
trange = [-1 2];

darkgreen = [8 129 54];
lightgreen = [192 212 149];
lightblue = [147 173 215];
darkblue = [45 80 158];
lightred = [242 150 154];
darkred = [231 57 30];
colourmap = [darkgreen;lightgreen;darkblue;lightblue;darkred;lightred]/255;
colourmap2 = [darkgreen;darkblue;darkred]/255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subplot A: PDR for Expt1A (Expt6)
expIDX = '1A';
pct = pct+1;
subplot(spx,spy,pct);

plotPDR_Expt1;
title('Expt1A');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subplot B: PDR for Expt1A (Expt6)
expIDX = '1B';
pct = pct+1;
subplot(spx,spy,pct);

plotPDR_Expt1;
title('Expt1B');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
message_sigtime = [message_sigtime; '************* Expt2 REG10-RAND10d ***************'];
yrange = [-.2 .5];
tw_stats = [0 3];
trange = [-1 3];

lightblue = [147 173 215];%reg10
tulipblue = [0,128,255];%reg10-rand10
mint = [0,255,255]; %reg10-rand10d
darkblue = [45 80 158]; %reg10-rand20

lightred = [242 150 154]; %rand20
darkred = [231 57 30]; %rand20-reg10
lightpink = [255,153,255]; %rand10
orange = [240 127 25]; %rand10-reg10d
hotpink = [255 0 255]; %rand10-reg10

colourmap = [lightblue;tulipblue;mint;darkblue;...
    lightpink;hotpink;orange;...
    lightred;darkred]/255;
colourmap2 = [0 128 255;0 255 255; 0 0 255;...  % REG10
    255 0 255; 255 128 0;...                % RAND10
    255 0 0;...                             % RAND20
    0 204 204; 0 102 204; 0 0 153; ...      % REG10-RAND10 vs other two REG-RAND conditions; REG10-RAND10d vs REG10-RAND20
    255 204 153;...                         % RAND10-REG10 vs RAND10-REG10d
    ]/255; %cluster hoizontal bar

condsName = {
    'REG10';'REG10-RAND10';'REG10-RAND10d';'REG10-RAND20';...
    'RAND10';'RAND10-REG10';'RAND10-REG10d';...
    'RAND20';'RAND20-REG10'};
statpairs=[
    2,1; 3,1; 4,1;...                       % REG10: 1 2 3
    6,5;7,5;...                             % RAND10: 4 5
    9,8;...                                 % RAND20: 6
    2,3; 2,4; 3,4; ...                      % REG10-RAND10 vs other two REG-RAND conditions; REG10-RAND10d vs REG10-RAND20
    6,7];                                   % RAND10-REG10 vs RAND10-REG10d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subplot A: PDR for Expt2 (Expt11)

pct = pct+1+1;
subplot(spx,spy,pct);

expIDX = '2';
plotlistCondition = [1,4,8,9];
plotlistStats = [3,7];
plotPDR_Expt2;
title('Expt2');

pct = pct+1;
subplot(spx,spy,pct);
expIDX = '2';
plotlistCondition = [1,2,3,4];
plotlistStats = [1,2,3];
plotPDR_Expt2;
title('Expt2');

pct = pct+1;
subplot(spx,spy,pct);
expIDX = '2';
plotlistCondition = [1,3,5,7];
plotlistStats = [2,5];
plotPDR_Expt2;
title('Expt2');

pct = pct+1;
subplot(spx,spy,pct);
expIDX = '2';
plotlistCondition = [5,6,7,8,9];
plotlistStats = [4,5,6];
plotPDR_Expt2;
title('Expt2');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
message_sigtime = [message_sigtime; '************* Expt3A Active ***************'];
darkgreen = [8 129 54];
lightgreen = [192 212 149];
lightblue = [147 173 215];
darkblue = [45 80 158];
lightred = [242 150 154];
darkred = [231 57 30];
colourmap = [darkgreen;lightgreen;darkblue;lightblue;darkred;lightred]/255;
colourmap2 = [darkgreen;darkblue;darkred]/255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pct = pct+1+2;
subplot(spx,spy,pct);
yrange = [-.5 1.25];
expIDX = '3A';
plotPDR_Expt1;
title('Expt3A');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pct = pct+1;
subplot(spx,spy,pct);
yrange = [-.5 1.25];
expIDX = '3A';
plotPDR_Expt1;

% rt_plot = SZ_6_Behav_Active_Expt3A;
% for k = 1:numel(rt_plot)
%     x0 = [rt_plot(k) rt_plot(k)];
%     y0 = ylim;
%     line(x0,y0,'Color',colourmap2(k,:),'LineStyle','--');
% end
title('Expt3A - from response');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
message_sigtime = [message_sigtime; '************* Expt3B Active (delay response) ***************'];
pct = pct+1+1;
subplot(spx,spy,pct);

yrange = [-.2 .8];
expIDX = '3B';
experiment = [expIDX 'hi_PDRatTrans_noBC'];
plotPDR_Expt1;
title('Expt3B');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
message_sigtime = [message_sigtime; '************* Expt4A RAND20-REG1 ***************'];

yrange = [-.2 .5];
tw_stats = [0 3];
trange = [-1 3];

darkgreen = [8 129 54];
lightgreen = [192 212 149];
lightblue = [147 173 215];
darkblue = [45 80 158];
lightred = [242 150 154];
darkred = [231 57 30];

hotpink = [255 0 255];
darkpurple = [69 39 119];
lightpurple = [147 112 219];

colourmap = [lightblue;darkblue;lightred;darkred;hotpink;lightgreen;darkgreen;]/255;
colourmap2 = [darkblue;darkred;hotpink;darkgreen]/255;
condsName = {'REG10','REG10-RAND20','RAND20','RAND20-REG10','RAND20-REG1','CONT','STEP'};

pct = pct+1+2;
subplot(spx,spy,pct);
expIDX = '4A';
ntrial = 16;
plotPDR_Expt4A;
title('Expt4A');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
message_sigtime = [message_sigtime; '************* Expt4B RAND20-REGx ***************'];

yrange = [-.2 .5];
tw_stats = [0 3];
trange = [-1 3];


condsIndex = 1:5;
condsName = {'RAND20','RAND20-REG1','RAND20-REG2','RAND20-REG5','RAND20-REG10'};
statpairs = [2,1;3,1;4,1;5,1; 2,3; 2,4; 2,5; 3,4];
cmap1 = [255 153 153; 255 0 255; 147,112,219; 75,0,130; 255 0 0]/255;

hotpink = [255 0 255];
darkpurple = [69 39 119];
lightpurple = [147 112 219];

lightred = [242 150 154];
darkred = [231 57 30];

colourmap = [lightred; hotpink; lightpurple; darkpurple; darkred]/255;
colourmap2 = [hotpink; lightpurple; darkpurple; darkred; 0 0 0; 64 64 64; 128 128 128; 224 224 224]/255;


pct = pct+1+2;
subplot(spx,spy,pct);
expIDX = '4Bg1';
plotPDR_Expt4B;
pupilmeanA = pupilmean;
sublistA = sublist;
title('Expt4B Group A');

pct = pct+1;
subplot(spx,spy,pct);
expIDX = '4Bg2';
plotPDR_Expt4B;
pupilmeanB = pupilmean;
sublistB = sublist;
title('Expt4B Group B');

pct = pct+1;
subplot(spx,spy,pct);
pupilmean = [pupilmeanA;pupilmeanB];
sublist = [sublistA;sublistB];

% Compute mean and std
cond_mean = squeeze(nanmean(pupilmean,1));
cond_std = [];
for k = 1:length(condsName)
    tmp = [];
    cond_std(k,:) = nanstd(squeeze(pupilmean(:,k,:)))/sqrt(length(sublist));
end

% find peak and time for each PDR mean, and save to diplay at the end

for k = 1:length(condsName)
    [peaky,peakx] = max(cond_mean(k,:));
    peakx = timeaxis(peakx);
    mes = ['Peak of ' condsName{k} ' : x = ' num2str(peakx) ', y = ' num2str(peaky)];
    message_sigtime = [message_sigtime; mes];
end

%% Start to plot
hold on;
for k = 1:length(condsName)
    a = cond_mean(k,:)';
    b = cond_std(k,:)';
    
    curve1 = a+b;
    curve2 = flipud(a-b);
    X = [timeaxis'; flipud(timeaxis')];
    Y = [curve1; curve2];
    figfill = fill(X,Y,colourmap(k,:),'edgecolor','none','facealpha',0.2);
    %         set(get(get(figfill,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
end

for k = 1:length(condsName)
    a = cond_mean(k,:);
    plot(timeaxis,a,'LineWidth',3,'Color',colourmap(k,:));
end
% set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
hold off;

if ~isempty(yrange); ylim(yrange); end

clusterstats = 1;

if clusterstats
    
    MM = [];
    xx = pupilmean;
    aa = ylim; Ystat = aa(1); aa = (aa(2)-aa(1))/30; Ystat=Ystat+aa;
    
    
    statstime_start = tw_stats(1); statstime_end = tw_stats(2);
    timepos = timeaxis(fFindClosestPosition(timeaxis,trange(1)):fFindClosestPosition(timeaxis,trange(2)));
    statstime=find((timepos>statstime_start) & (timepos<=statstime_end));
    
    for k = 1:size(statpairs,1)
        
        cond1 = squeeze(xx(:,statpairs(k,1),statstime))';
        cond2 = squeeze(xx(:,statpairs(k,2),statstime))';
        
        cfg = [];
        cfg.statistic        = 'ft_statfun_depsamplesT';
        cfg.numrandomization = 1000;
        cfg.correctm         = 'cluster';
        cfg.method           = 'montecarlo';
        cfg.tail             = 0;
        cfg.alpha            = 0.05;
        cfg.clusteralpha     = 0.05;
        cfg.clusterstatistic = 'maxsize';
        cfg.design           = [1:numel(sublist) 1:numel(sublist) % subject number
            ones(1,numel(sublist)) 2*ones(1,numel(sublist))];  % condition number
        cfg.uvar = 1;        % "subject" is unit of observation
        cfg.ivar = 2;        % "condition" is the independent variable
        cfg.dimord = 'time';
        %             cfg.dim=[1,numel(timeaxis)];
        cfg.dim = [1,numel(statstime)];
        cfg.connectivity =1;
        
        stat = ft_statistics_montecarlo(cfg, [cond1 cond2],cfg.design);
        
        disp([condsName{statpairs(k,1)} ' > ' condsName{statpairs(k,2)} ' : ' num2str(min(stat.prob))]);
        MM(k)=min(stat.prob);
        
        % Find indices of significant clusters
        pos=[]; neg=[];
        if isfield(stat,'posclusters')
            if ~isempty(stat.posclusters)
                pos_cluster_pvals = [stat.posclusters(:).prob];
                pos_signif_clust = find(pos_cluster_pvals < cfg.alpha);
                poss = ismember(stat.posclusterslabelmat, pos_signif_clust);
                if size(find(diff([0; poss])==-1),1) ~= size(find(diff([0; poss])==1),1)
                    bb = [find(diff([0; poss])==-1); length(poss)];
                    pos = [find(diff([0; poss])==1) bb];
                else
                    pos = [find(diff([0; poss])==1) find(diff([0; poss])==-1)];
                end
            end
        end
        if isfield(stat,'negclusters')
            if ~isempty(stat.negclusters)
                neg_cluster_pvals = [stat.negclusters(:).prob];
                neg_signif_clust = find(neg_cluster_pvals <cfg.alpha);
                negs = ismember(stat.negclusterslabelmat, neg_signif_clust);
                if size(find(diff([0; negs])==-1),1) ~= size(find(diff([0; negs])==1),1)
                    bb = [find(diff([0; negs])==-1); length(negs)];
                    neg = [find(diff([0; negs])==1) bb];
                else
                    neg = [find(diff([0; negs])==1) find(diff([0; negs])==-1)];
                end
                %                     neg = [find(diff([0; negs])==1) find(diff([0; negs])==-1)];
            end
        end
        
        hold on;
        for i = 1:size(pos,1)
            sigtime = [timeaxis(fFindClosestPosition(timeaxis,statstime_start)+pos(i,1)) timeaxis(fFindClosestPosition(timeaxis,statstime_start)+pos(i,2))];
            message_sigtime = [message_sigtime; [condsName{statpairs(k,1)} ' > ' condsName{statpairs(k,2)} ' : ' num2str(sigtime(1)) '~' num2str(sigtime(2)) ' p=' num2str(min(stat.prob))]];
            l = line(sigtime,[Ystat Ystat],'LineWidth',5,'Color',colourmap2(k,:));hold on
            %                 set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end
        for i = 1:size(neg,1)
            sigtime = [timeaxis(fFindClosestPosition(timeaxis,statstime_start)+neg(i,1)) timeaxis(fFindClosestPosition(timeaxis,statstime_start)+neg(i,2))];
            message_sigtime = [message_sigtime; [condsName{statpairs(k,1)} ' < ' condsName{statpairs(k,2)} ' : ' num2str(sigtime(1)) '~' num2str(sigtime(2))] ' p=' num2str(min(stat.prob))];
            l = line(sigtime,[Ystat Ystat], 'LineWidth',5,'Color',colourmap2(k,:)); hold on
            %                 set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end
        Ystat = Ystat+aa;
        
    end
    hold off;
    legend off;
    
    disp('------------------------------------------------------');
    disp(message_sigtime);
    
    disp('------------------------------------------------------');
    disp(MM);
end
xlim(trange);
title('Expt4B Group A+B');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 28]);
filename = ['fig_PDR_ALLexpt'];
saveas(gcf,[filename,'.png']);
