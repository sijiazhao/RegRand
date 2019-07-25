clear;close all;

expIDX = '1'; samfreq = 1000; % [Hz]
statspair = [1,2;4,3;6,5];
cl=[1,2,4,3,6,5];
cmap1 = [0 0.5 0;0.7647 0.8392 0.6078;0.5843 0.6863 0.8431; 0 0 1; 0.9020 0.7255 0.7216;1 0 0;0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0];
figname0 = 'fig3_expt1_PDevent';

% expIDX = '4'; samfreq = 1000; % [Hz]
% statspair = [2,1;3,1;4,1;5,1];
% cl = [2 3 4 5 1]; %order to plot
% figname0 = 'figs3_expt4B_PDevent';

run_mode = {'PD','PC'};
run_thrs = [75,300]; colourmap_thrs = [198 179 119; 0 0 0]/225;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
twindow = 500; % [ms] time window for counting
filename = ['PDevent_expt' expIDX];
load([filename '.mat']);
condsName = design.conds;
numcond = size(P,2);

plotoverlay = 0; plottrange = [-2 0]; % overlap STEP, REG-RAND, RAND-REG for pre-transition
ifcollapse = 0;
trange = [-2 2]; %plot time range
yrange = [];
ifttest = 1;

nfig = 0; %figure idex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PD events ~~ Figure 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:numel(run_mode)
    nfig = nfig+1;
    figure(nfig);clf;
    mode = run_mode{k};
    for j = 1:numel(run_thrs)
        thrs = run_thrs(j);
        load(['pdEvent_expt',expIDX,'_', mode, '_',num2str(thrs),'.mat'],'PM');
        
        A = PM.pdEvent;
        numcond = size(A,2);
        
        ccount = 0;
        for cond = cl % The condi tion order need to be changed for rasterplot: change, no-change, change, no-change...
            
            ccount = ccount+1;
            
            % >>> specific to expt4 <<<
            if strcmp(filename(1:2),'14')
                disp('Assign specific subplot layout for EXPT14');
                if cond == 1 %RAND20 has 2400 trials = 4 times of other 4 conditions
                    s = subplot(8,1,5:8);
                else
                    s = subplot(8,1,cond-1);
                end
            else
                s = subplot(numcond,1,ccount);
            end
            % <<< specific to expt4 >>>
            
            Y=[]; Ytoplot=[]; N=0; M=0;
            
            hold on;
            for subj = 1:numel(sublist)
                spikes = logical(A{subj,cond});
                nTrials = size(spikes,1);
                nTimeBins = size(spikes,2);
                
                [y,y2] = find(spikes);
                y = y+N;
                y2 = timeaxis(y2);
                
                s=scatter(y2,y,8,'o','Filled');
                s.MarkerFaceColor = colourmap_thrs(j,:);
                
                Y = [Y;y];
                Ytoplot = [Ytoplot;y2'];
                N = N+nTrials;
                M = nTimeBins;
            end
            hold off;
            ylim([1 N]);
            xlim(trange);
            
            title(condsName{cond});
            set(gca, 'XTick', []);
            
        end
    end
    
    ylabel('Trials');
    xlabel('Time from transition [s]');
    set(gca, 'XTick', [-2 -1.5 -1 -.5 0 .5 1 1.5 2]);
    
    figname = [figname0 '_' run_mode{k} '_raster'];
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 10]);
    print(figname,'-dpng','-r0');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PD events ~~ Figure 2: rate
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['pdRate_expt',expIDX,'.mat'],'Nonoverlapping','RunningMean','ConvolvedRate');
run_compute = {'Nonoverlapping','RunningMean','ConvolvedRate'};

for k = 1:numel(run_mode)
    nfig = nfig+1;
    figure(nfig);clf;
    
    mode = run_mode{k};
    for r = 3 %ConvolvedRate
        
        eval(['data = ',run_compute{r},';']);
        
        for j = 1:numel(run_thrs)
            thrs = run_thrs(j);
            disp(['..Plotting ',run_compute{r}, ' for threshold=' num2str(thrs)]);
            
            R = data{k,j};
            
            for iPair = 1:size(statspair,1)
                subplot(size(statspair,1),1,iPair);
                
                hold on;
                for cond = statspair(iPair,:)
                    tmp = squeeze(R(cond,:,:));
                    Ytoplot = nanmean(tmp);
                    
                    switch run_compute{r}
                        case 'Nonoverlapping'
                            err = nanstd(tmp)/sqrt(size(tmp,1));
                            
                            if numel(timeaxis)~=numel(Ytoplot)
                                twindow2 = twindow/(1000/samfreq); %twindow is set based on samfreq = 1000Hz
                                
                                timeaxis1 = downsample(timeaxis,floor((numel(timeaxis)-1)/(numel(Ytoplot)-1)));
                                timeaxis1 = timeaxis1+twindow2/samfreq/2; %shift the time axis to locate the result in the middle of time window
                            end
                            plottrange_pos = [fFindClosestPosition(timeaxis1,(trange(1)+0.25)):fFindClosestPosition(timeaxis1,trange(2)-0.25)];
                            errorbar(timeaxis1(plottrange_pos),Ytoplot(plottrange_pos),err(plottrange_pos),'color',cmap1(cond,:),'LineWidth',2);
                            
                        otherwise
                            
                            plot(timeaxis,Ytoplot,'color',cmap1(cond,:),'LineWidth',2);
                    end
                end
                clear tmp;
                
                xlim(trange);
                
                switch run_compute{r}
                    case 'Nonoverlapping'
                        yrange = [0,1];
                    case 'RunningMean'
                        yrange = [0,0.005];
                    case 'ConvolvedRate'
                        yrange = [0 0.4];
                end
                
                if ~isempty(yrange), ylim(yrange);end
                
                xl = xlim; yl = ylim;
                switch run_compute{r}
                    case 'Nonoverlapping'
                        for t = 1:numel(timeaxis1)
                            tmp1 = squeeze(R(statspair(iPair,1),:,t));
                            tmp2 = squeeze(R(statspair(iPair,2),:,t));
                            if timeaxis1(t) < xl(2)
                                try
                                    [h,p]=ttest2(tmp1,tmp2,'alpha',0.01);
                                    % add the text to the plot
                                    if h
                                        %                 if timeaxis(t)<xlim(2)
                                        x3 = timeaxis1(t);
                                        if j == 1, y3=yl(2); elseif j==2, y3=yl(1); end
                                        %                                         y3 = yl(j);%-0.1*(yl(2)-yl(1));
                                        %                 txt3 = [num2str(x3) ' = ' sprintf('%0.2f',p)];
                                        txt3 = [sprintf('%0.2f',p)];
                                        text(x3,y3,txt3,'FontSize',8,'HorizontalAlignment','center');
                                        %                 end
                                    end
                                end
                            end
                        end
                    otherwise
                        if j == 1, Ystat = yl(2); yl = (yl(2)-yl(1))/30; Ystat=Ystat-yl;
                        elseif j==2, Ystat = yl(1); yl = (yl(2)-yl(1))/30; Ystat=Ystat+yl; end
                        
                        cond1 = squeeze(R(statspair(iPair,1),:,:))';
                        cond2 = squeeze(R(statspair(iPair,2),:,:))';
                        
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
                        cfg.dim=[1,numel(timeaxis)];
                        %                     cfg.dim = [1,numel(statstime)];
                        cfg.connectivity = 1;
                        
                        stat = ft_statistics_montecarlo(cfg, [cond1 cond2],cfg.design);
                        
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
                        
                        % _____ PLOT______
                        hold on;
                        sig = [pos;neg];
                        if ~isempty(sig)
                            for i = 1:size(sig,1)
                                sigtime = [timeaxis(sig(i,1)),timeaxis(sig(i,2))];
                                l = line(sigtime,[Ystat Ystat],'LineWidth',5,'Color','k');hold on
                            end
                        end
                        Ystat = Ystat+yl;
                end
                
                legend off;
                
                disp('------------------------------------------------------');
                
                if iPair == 1
                    ylabel([run_compute{r} mode ' rate [/s]']);
                elseif iPair == size(statspair,1)
                    xlabel('Time from transition');
                end
            end
            
        end
        hold off;
        
        figname = [figname0 '_' run_mode{k} '_' run_compute{r}];
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 10]);
        print(figname,'-dpng','-r0');
    end
end
