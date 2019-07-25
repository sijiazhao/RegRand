filename = ['pupil_expt' expIDX '_PDRatTrans'];
load(filename);

%% **** Down sample for boostrap analysis ************************
ds = 1;
rate = 20; %[Hz]
if ds
    timeaxis = downsample(timeaxis,rate);
    for subj = 1:size(P,1)
        for k = 1: size(P,2)
            p = P{subj,k};
            p_new = NaN(size(p,1),numel(timeaxis));
            for trial= 1:size(p,1)
                Pn = p(trial,:);
                if any(isnan(Pn));end
                Pn = downsample(Pn,rate);
                p_new(trial,:) = Pn;
            end
            P{subj,k} = p_new;
        end
    end
end

%% Start regression
tw_bc = [-1,0];
mywindow_baseline = [find(timeaxis == tw_bc(1)):find(timeaxis == tw_bc(2))];
mywindow_post = [find(timeaxis == 0):length(timeaxis)];

condlist = condsName;

valBaseline = [];

for subj = 1:size(P,1)
    for k = 1: size(P,2)
        p = P{subj,k};
        for t = 1:size(p,1)
            pn = p(t,:);
            vb = nanmean(pn(mywindow_baseline));
            valBaseline(subj,k,t) = vb;
        end
    end
end

P0 = P;
beta = [];
for s = 1:numel(sublist)
    for k = 1:length(condsName)
        p = P{s,k};
        p0 = p; % this is for beta = 1;
        p1 = p; % beta varies along time
        for t = 1:size(p,1)
            pn = p(t,:);
            vb = nanmean(pn(mywindow_baseline));
            valBaseline(s,k,t) = vb;
        end
        
        VB = squeeze(valBaseline(s,k,:));    
        
        for z = 1:size(p,2)
            N = size(p,1); % ntrial
            VB = VB(1:N); %THIS IS HARD CODED!
            x = reshape(VB,N,1);
            y = reshape(p(:,z),N,1);
            
            I = find(isnan(y));
            y(I) = [];
            x(I) = [];
            
            I = find(isnan(x));
            x(I) = [];
            y(I) = [];
            
            r = corrcoef(x,y); % Corr coeff is the off-diagonal (1,2) element
            r = r(1,2);  % Sample regression coefficient
            
            sigx = std(x);
            sigy = std(y);
            a1 = r*sigy/sigx;   % Regression line slope
            beta(s,k,z) = a1;
            p1(:,z) = p(:,z)-a1*VB;
            p0(:,z) = p0(:,z) - VB;
        end
        
        P{s,k} = p1;
        P0{s,k} = p0;
    end
end

%% **** Compute mean for each condition each subject
pupilmean = [];
for subj = 1:size(P,1)
    for k = 1:size(P,2)
        p = P{subj,k};
        pupilmean(subj,k,:) = nanmean(p,1);
    end
end

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
        
        % _____ PLOT______
        hold on;

        for i = 1:size(pos,1)
            sigtime = [timeaxis(fFindClosestPosition(timeaxis,statstime_start)+pos(i,1)) timeaxis(fFindClosestPosition(timeaxis,statstime_start)+pos(i,2))];
%             message_sigtime = [message_sigtime; [condsName{statpairs(k,1)} ' > ' condsName{statpairs(k,2)} ' : ' num2str(sigtime(1)) '~' num2str(sigtime(2)) ' p=' num2str(min(stat.prob))]];
            l = line(sigtime,[Ystat Ystat],'LineWidth',5,'Color',colourmap2(k,:));hold on
            %                 set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end
        for i = 1:size(neg,1)
            sigtime = [timeaxis(fFindClosestPosition(timeaxis,statstime_start)+neg(i,1)) timeaxis(fFindClosestPosition(timeaxis,statstime_start)+neg(i,2))];
%             message_sigtime = [message_sigtime; [condsName{statpairs(k,1)} ' < ' condsName{statpairs(k,2)} ' : ' num2str(sigtime(1)) '~' num2str(sigtime(2))] ' p=' num2str(min(stat.prob))];
            l = line(sigtime,[Ystat Ystat], 'LineWidth',5,'Color',colourmap2(k,:)); hold on
            %                 set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end
        Ystat = Ystat+aa;
        
    end
    hold off;
    legend off;
    
    disp('------------------------------------------------------');
%     disp(message_sigtime);
    
    disp('------------------------------------------------------');
    disp(MM);
end
xlim(trange);
