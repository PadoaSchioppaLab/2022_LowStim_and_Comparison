% FigSummary_plot_generator_allcurrent.m
% This script uses complete_dataset.xlsx to generate a figure to summarize
% all the effects

% Authors: Sebastien Ballesta & Weikang Shi

% Copyright: Camillo Padoa-Schioppa lab, Washington Univerty in St. Louis

close all
clearvars

% different definition of hysteresis
def_hyst = 0; % '0': hysteresis and bias are defined as a2; '1': defined as a2/a1; '2': defined as 2rho*a2/a1

% only <15uA histogram
only15uA = 0;

% % load data and info files
% % if use complete_dataset.xlsx
[~,~,Exp1_info] = xlsread('session_info.xlsx', 1); 
[~,~,Exp1_data] = xlsread('session_detailed_data.xlsx', 1); 
info=cell2table(Exp1_info(2:end,:),'VariableNames',Exp1_info(1,:));
alldata=cell2table(Exp1_data(2:end,:),'VariableNames',Exp1_data(1,:));


% define analysis type and figure axis labels
% different current levels
limlow1  = 1;
limhigh1 = 16;
limlow2  = 25;
limhigh2 = 50;
limhigh3 = 100;

% pval_thres = 0.05;

% parsup = [0 limlow1  limlow2  limhigh2 limhigh3];  % lower boundary to choose target current level
% parinf = [1 limhigh1 limhigh2 limhigh3 Inf     ];  % upper boundary to choose target current level
% currnames = { '0 \muA'; '\leq15\muA'; '25\muA'; '50\muA'; '\geq100\muA'};
%
parsup = [limlow1  limlow2  limhigh2 limhigh3];  % lower boundary to choose target current level
parinf = [limhigh1 limhigh2 limhigh3 Inf     ];  % upper boundary to choose target current level
currnames = {'\leq15\muA'; '25\muA'; '50\muA'; '125\muA'};
%
% parsup = [limlow1];  % lower boundary to choose target current level
% parinf = [limhigh1];  % upper boundary to choose target current level
% currnames = {'\leq15\muA'};
%
ncurrs = length(currnames);
% 
% get new fitting
[info_sepeOrder] = probitLink_regression_seperateOrder(alldata, info);
[info_sameOrder] = probitLink_regression_sameOrder(alldata, info);

% for range dependent effect
AR = info_sameOrder.relative_value_ON-info_sameOrder.relative_value_OFF; % relative value difference
mrho = info_sameOrder.relative_value;                                    % mean relative value
RR = info_sameOrder.rangeA-info_sameOrder.rangeB;                        % value range difference
isplit = info_sameOrder.isplit;                                          % index for good sessions fulfilling incomplete separation
isat = info_sameOrder.isat;                                              % index for good sessions fulfilling saturation
stimcurr = info_sameOrder.StimCurrent;                                   % stimulation current level
stimlat = info_sameOrder.StimLateral;                                    % uni or bilateral micro-stimulation

% for other analysis
% Rho_OFFON       = [ info_sepeOrder.relative_value_OFF info_sepeOrder.relative_value_ON ];
% Steepness_OFFON = [ info_sepeOrder.steepness_OFF info_sepeOrder.steepness_ON ];
Rho_OFFON       = [ info.relative_value_OFF info.relative_value_ON ];
Steepness_OFFON = [ info.steepness_OFF info.steepness_ON ];
if def_hyst == 0
%     OrdBias_OFFON  = [ info_sepeOrder.order_bias_OFF info_sepeOrder.order_bias_ON ];
    OrdBias_OFFON  = [ info.order_bias_OFF info.order_bias_ON ];
elseif def_hyst == 1
    OrdBias_OFFON  = [ info_sepeOrder.order_bias_OFF info_sepeOrder.order_bias_ON ]./[ info_sepeOrder.steepness_OFF info_sepeOrder.steepness_ON ];
elseif def_hyst == 2
    OrdBias_OFFON = [ info_sepeOrder.relative_value_OFF info_sepeOrder.relative_value_ON ].*...
            [ info_sepeOrder.order_bias_OFF info_sepeOrder.order_bias_ON ]./[ info_sepeOrder.steepness_OFF info_sepeOrder.steepness_ON ];
end

anaTypes = {'range-dep bias'; 'order bias(offer1)'; 'order bias(offer2)'; 'steepness(offer1)'; 'steepness(offer2)'};
anaColrs = [ 1.0 0.5 0.0; 0.6 1.0 1.0; 0.2 0.8 0.6;  0.8 0.8 0.8; 0.6 0.0 0.6];
monkeys  = {'both'}; 
nmonkeys = length(monkeys);

for imon = 1:nmonkeys
monkeyname = monkeys{imon};

% initialize
summaryData = [];
summaryEror = [];
summaryPval = [];
summaryN    = [];
if ismember(monkeyname,'both')
    mask_mon = logical(ones(size(info_sepeOrder.monkey)));
    pval_thres = 0.05/5;
else
    mask_mon = ismember(info_sepeOrder.monkey,monkeyname);
    pval_thres = 0.05/5;
end
for icur = 1:ncurrs    
    %
    % range dependent effects
    mask = info_sameOrder.StimCurrent>=parsup(icur) & info_sameOrder.StimCurrent<parinf(icur) & mask_mon;
    limout = 3;
    K = limout;
    OUT_RR = (RR>(mean(RR)+std(RR)*K)) | (RR<(mean(RR)-std(RR)*K));  
    OUT_AR = (AR>(mean(AR)+std(AR)*K)) | (AR<(mean(AR)-std(AR)*K));  
    OUT = OUT_RR | OUT_AR;
    mask = logical(~OUT.*mask);
    [Rpe,Ppe]=corr(RR(mask),AR(mask),'Type','Pearson');
    [Rsp,Psp]=corr(RR(mask),AR(mask),'Type','Spearman');
    summaryData(icur,1) = Rpe; 
    summaryEror(icur,1) = sqrt((1-Rpe.^2)./(sum(mask)-2));
    summaryPval(icur,1) = Ppe;
    % summaryData(icur,1) = Rsp; 
    % summaryEror(icur,1) = sqrt((1-Rsp.^2)./(sum(mask)-2));
    % summaryPval(icur,1) = Psp;
    summaryN(icur,1) = sum(mask);
    
    %
    % order bias (offer1)
    if parinf(icur) == 1
        mask = info_sepeOrder.StimCurrent>=parsup(icur) & info_sepeOrder.StimCurrent<parinf(icur) & mask_mon;
    else
        mask = info_sepeOrder.StimOffer(:,1)==1 & info_sepeOrder.StimCurrent>=parsup(icur) & info_sepeOrder.StimCurrent<parinf(icur) & mask_mon;
    end
    OrdBias_OFFON_icur = OrdBias_OFFON(mask,:);
    summaryData(icur,2) = nanmean(OrdBias_OFFON_icur(:,2)-OrdBias_OFFON_icur(:,1));   
    summaryEror(icur,2) = nanstd(OrdBias_OFFON_icur(:,2)-OrdBias_OFFON_icur(:,1))./sqrt(size(OrdBias_OFFON_icur,1));
    summaryPval(icur,2) = signrank(OrdBias_OFFON_icur(:,2)-OrdBias_OFFON_icur(:,1));
    % [~,summaryPval(icur,2)] = ttest(OrdBias_OFFON_icur(:,2)-OrdBias_OFFON_icur(:,1));
    summaryN(icur,2) = size(OrdBias_OFFON_icur,1);
    
    %
    % order bias (offer2)
    if parinf(icur) == 1
        mask = info_sepeOrder.StimCurrent>=parsup(icur) & info_sepeOrder.StimCurrent<parinf(icur) & mask_mon;
    else
        mask = info_sepeOrder.StimOffer(:,1)==2 & info_sepeOrder.StimCurrent>=parsup(icur) & info_sepeOrder.StimCurrent<parinf(icur) & mask_mon;
    end
    OrdBias_OFFON_icur = OrdBias_OFFON(mask,:);
    summaryData(icur,3) = nanmean(OrdBias_OFFON_icur(:,2)-OrdBias_OFFON_icur(:,1));    
    summaryEror(icur,3) = nanstd(OrdBias_OFFON_icur(:,2)-OrdBias_OFFON_icur(:,1))./sqrt(size(OrdBias_OFFON_icur,1));
    summaryPval(icur,3) = signrank(OrdBias_OFFON_icur(:,2)-OrdBias_OFFON_icur(:,1));    
    % [~,summaryPval(icur,3)] = ttest(OrdBias_OFFON_icur(:,2)-OrdBias_OFFON_icur(:,1)); 
    summaryN(icur,3) = size(OrdBias_OFFON_icur,1);
    
    %
    % steepness (offer1)
    if parinf(icur) == 1
        mask = info_sepeOrder.StimCurrent>=parsup(icur) & info_sepeOrder.StimCurrent<parinf(icur) & mask_mon;
    else
        mask = info_sepeOrder.StimOffer(:,1)==1 & info_sepeOrder.StimCurrent>=parsup(icur) & info_sepeOrder.StimCurrent<parinf(icur) & mask_mon;
    end
    Steepness_OFFON_icur = Steepness_OFFON(mask,:);
    summaryData(icur,4) = nanmean(Steepness_OFFON_icur(:,2)-Steepness_OFFON_icur(:,1));
    summaryEror(icur,4) = nanstd(Steepness_OFFON_icur(:,2)-Steepness_OFFON_icur(:,1))./sqrt(size(Steepness_OFFON_icur,1));
    summaryPval(icur,4) = signrank(Steepness_OFFON_icur(:,2)-Steepness_OFFON_icur(:,1));
    % [~,summaryPval(icur,4)] = ttest(Steepness_OFFON_icur(:,2)-Steepness_OFFON_icur(:,1));
    summaryN(icur,4) = size(Steepness_OFFON_icur,1);
    
    %
    % steepness (offer2)
    if parinf(icur) == 1
        mask = info_sepeOrder.StimCurrent>=parsup(icur) & info_sepeOrder.StimCurrent<parinf(icur) & mask_mon;
    else
        mask = info_sepeOrder.StimOffer(:,1)==2 & info_sepeOrder.StimCurrent>=parsup(icur) & info_sepeOrder.StimCurrent<parinf(icur) & mask_mon;
    end
    Steepness_OFFON_icur = Steepness_OFFON(mask,:);
    summaryData(icur,5) = nanmean(Steepness_OFFON_icur(:,2)-Steepness_OFFON_icur(:,1));
    summaryEror(icur,5) = nanstd(Steepness_OFFON_icur(:,2)-Steepness_OFFON_icur(:,1))./sqrt(size(Steepness_OFFON_icur,1));
    summaryPval(icur,5) = signrank(Steepness_OFFON_icur(:,2)-Steepness_OFFON_icur(:,1));
    % [~,summaryPval(icur,5)] = ttest(Steepness_OFFON_icur(:,2)-Steepness_OFFON_icur(:,1));
    summaryN(icur,5) = size(Steepness_OFFON_icur,1);
end

% %
% normalize to the maximum
max_summaryData = max(abs(summaryData));
% order1 steepness and order2 steepness normalize to the same max
max_steepness = max(max_summaryData(1,4:5));
max_summaryData(1,4:5) = [max_steepness, max_steepness];
% order1 order bias and order2 order bias normalize to the same max
max_orderbias = max(max_summaryData(1,2:3));
max_summaryData(1,2:3) = [max_orderbias, max_orderbias];
% normalize
summaryData_norm = summaryData./(ones(numel(parsup),1)*max_summaryData);
summaryEror_norm = summaryEror./(ones(numel(parsup),1)*max_summaryData);
% reverse effect to make positive mean excepted effect direction
summaryData_norm_rev = summaryData_norm;
summaryData_norm_rev(:,3) = -summaryData_norm_rev(:,3);
summaryData_norm_rev(:,4) = -summaryData_norm_rev(:,4);
summaryData_norm_rev(:,5) = -summaryData_norm_rev(:,5);

%
% plot
if ~only15uA
    %
    %
    % figure 1 
    figure
    set(gcf,'position',[110 65 750 450], 'PaperPositionMode','auto');
    colormap(anaColrs);
    bplot = bar(summaryData_norm);
    % %
    % add error bar
    hold on
    % Find the number of groups and the number of bars in each group
    ngroups = size(summaryData_norm, 1);
    nbars = size(summaryData_norm, 2);
    % Calculate the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        x = x';
        errorbar(x, summaryData_norm(:,i), summaryEror_norm(:,i), 'k', 'linestyle', 'none','linewidth',1);
        % use * to represent significance
        iPval = summaryPval(:,i);
        iPval_thres = iPval<pval_thres;
        try
            plot(x(iPval_thres),ones(sum(iPval_thres),1)*1.44,'k*','MarkerSize',5);   
        end
    end
    % %
    hold off
    xticks([1:numel(parsup)]);
    xticklabels(currnames);
    axis([0,numel(parsup)+1,-1.68,1.68]);
    box off
    legend(anaTypes,'Location','northwest');
    set(gca,'FontSize',10)
    ylabel('normalized effect size');
    title({['monkey ',monkeyname],['normalized to maximal effects']});
    %
    % 
    % figure 2 
    % 
    figure
    set(gcf,'position',[110 65 750 450], 'PaperPositionMode','auto');
    colormap(anaColrs);
    bar(summaryData_norm_rev);
    % %
    % add error bar
    hold on
    % Find the number of groups and the number of bars in each group
    ngroups = size(summaryData_norm_rev, 1);
    nbars = size(summaryData_norm_rev, 2);
    % Calculate the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, summaryData_norm_rev(:,i), summaryEror_norm(:,i), 'k', 'linestyle', 'none','linewidth',1);
        % use * to represent significance
        iPval = summaryPval(:,i);
        iPval_thres = iPval<pval_thres;
        try
            plot(x(iPval_thres),ones(sum(iPval_thres),1)*1.44,'k*','MarkerSize',5);   
        end
    end
    % %
    xticks([1:numel(parsup)]);
    xticklabels(currnames);
    axis([0,numel(parsup)+1,-1.28,1.68]);
    box off
    legend(anaTypes,'Location','northwest');
    set(gca,'FontSize',10)
    ylabel('normalized effect size');
    title({['monkey ',monkeyname];['normalized to maximal effects'];['reverse effect to keep positive as expected']});


elseif only15uA
    %
    %
    % figure 1 
    figure
    set(gcf,'position',[110 65 450 450], 'PaperPositionMode','auto');
    colormap(anaColrs);
    nbars = length(summaryData_norm(1,:));
    bplot = bar(1:nbars,diag(summaryData_norm(1,:)),'stacked');
    hold on
    errorbar(1:nbars, summaryData_norm(1,:), summaryEror_norm(1,:), 'k', 'linestyle', 'none','linewidth',1);
    iPval = summaryPval(1,:);
    iPval_thres = find(iPval<pval_thres==1);
    plot(iPval_thres,ones(sum(iPval_thres),1)*1.44,'k*','MarkerSize',5); 
    hold off
    xticks([1:nbars]);
    xticklabels([]);
    xlabel('\leq15\muA');
    axis([0,nbars+1,-1.68,1.68]);
    box off
    legend(anaTypes,'Location','northwest');
    set(gca,'FontSize',10)
    ylabel('normalized effect size');
    title({['monkey ',monkeyname],['normalized to maximal effects']});
    %
    % 
    % figure 2 
    % 
    figure
    set(gcf,'position',[110 65 450 450], 'PaperPositionMode','auto');
    colormap(anaColrs);
    nbars = length(summaryData_norm_rev(1,:));
    bplot = bar(1:nbars,diag(summaryData_norm_rev(1,:)),'stacked');
    hold on
    errorbar(1:nbars, summaryData_norm_rev(1,:), summaryEror_norm(1,:), 'k', 'linestyle', 'none','linewidth',1);
    iPval = summaryPval(1,:);
    iPval_thres = find(iPval<pval_thres==1);
    plot(iPval_thres,ones(sum(iPval_thres),1)*1.44,'k*','MarkerSize',5); 
    hold off
    xticks([1:nbars]);
    xticklabels([]);
    xlabel('\leq15\muA');
    axis([0,nbars+1,-1.28,1.68]);
    box off
    legend(anaTypes,'Location','northwest');
    set(gca,'FontSize',10)
    ylabel('normalized effect size');
    title({['monkey ',monkeyname];['normalized to maximal effects'];['reverse effect to keep positive as expected']});

end % if else only15uA
end

% % %
% function 1
% % %
function [info] = probitLink_regression_seperateOrder(alldata, info)
%
% This function is to do logistic regression with "probit" link function
% in this regression:
% X = b0 + b1 * log(#B/#A) + b2 *(delta_AB - delta_BA)
% regression is done on stimOFF and stimON trials seperately
%
nsessions = size(info,1);
%
for isession = 1:nsessions
	session = info.session{isession};
	% load raw data from alldata
	ind_trials = ismember(alldata.Session, ['ST',session]);
	QtyA = alldata.QtyA(ind_trials,:);
	QtyB = alldata.QtyB(ind_trials,:);
	stim = alldata.Stim(ind_trials,:);
	chosenID = alldata.ChosenID(ind_trials,:);
	order = -alldata.Order(ind_trials,:);
	%
	% remove forced choices
	ii = ~prod([QtyA,QtyB],2)==0;
	%
	iOF = stim==-1;
	iON = stim==1;
	%
	logB_A = log(abs(QtyB)./abs(QtyA));
	%
	% logistic regression
	% for stimOFF
	jj = ii & iOF;
	% [~ , ~, stats] = glmfit([logB_A(jj,1),order(jj,1)], chosenID(jj,1), 'binomial', 'link','probit', 'constant','on');
	mdl = fitglm([logB_A(jj,1),order(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
	betaOFF = mdl.Coefficients.Estimate;
	rhoOFF = exp(-betaOFF(1)/betaOFF(2));
	steepnessOFF = betaOFF(2);
	orderbiasOFF = betaOFF(3);
	
	% for stimON
	jj = ii & iON;
	% [~ , ~, stats] = glmfit([logB_A(jj,1),order(jj,1)], chosenID(jj,1), 'binomial', 'link','probit', 'constant','on');
	mdl = fitglm([logB_A(jj,1),order(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
	betaON = mdl.Coefficients.Estimate;
	rhoON = exp(-betaON(1)/betaON(2));
	steepnessON = betaON(2);
	orderbiasON = betaON(3);
		
	% save data
	relative_value_OFF(isession,:) = rhoOFF;
	relative_value_ON(isession,:) = rhoON;
	steepness_OFF(isession,:) = steepnessOFF;
	steepness_ON(isession,:) = steepnessON;
	order_bias_OFF(isession,:) = orderbiasOFF;
	order_bias_ON(isession,:) = orderbiasON;
	
end % for isession

% save data to info
info.relative_value_OFF = relative_value_OFF;
info.relative_value_ON = relative_value_ON;
info.steepness_OFF = steepness_OFF;
info.steepness_ON = steepness_ON;
info.order_bias_OFF = order_bias_OFF;
info.order_bias_ON = order_bias_ON;
end


% % %
% function 2
% % %
function [info] = probitLink_regression_sameOrder(alldata, info)
%
% This function is to do logistic regression with "probit" link function
% in this regression:
% X = b0 + b1 * log(#B/#A)
% regression is done on stimOFF and stimON trials seperately
%
nsessions = size(info,1);
%
for isession = 1:nsessions
    session = info.session{isession};
    % load raw data from alldata
    ind_trials = ismember(alldata.Session, ['ST',session]);
    QtyA = alldata.QtyA(ind_trials,:);
    QtyB = alldata.QtyB(ind_trials,:);
    stim = alldata.Stim(ind_trials,:);
    chosenID = alldata.ChosenID(ind_trials,:);
    order = -alldata.Order(ind_trials,:);
    %
    % remove forced choices
    ii = ~prod([QtyA,QtyB],2)==0;
    %
    iOF = stim==-1;
    iON = stim==1;
    
    % logistic regression - based on single trials
    % for stimOFF
    logB_A = log(abs(QtyB)./abs(QtyA));
    jj = ii & iOF;
    mdl = fitglm([logB_A(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
    betaOFF = mdl.Coefficients.Estimate;
    rhoOFF = exp(-betaOFF(1)/betaOFF(2));
    steepnessOFF = betaOFF(2);
    % for stimON
    jj = ii & iON;
    mdl = fitglm([logB_A(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
    betaON = mdl.Coefficients.Estimate;
    rhoON = exp(-betaON(1)/betaON(2));
    steepnessON = betaON(2); 
               
    % save data
    relative_value(isession,:) = (rhoOFF + rhoON)/2;
    relative_value_OFF(isession,:) = rhoOFF;
    relative_value_ON(isession,:) = rhoON;
    steepness_OFF(isession,:) = steepnessOFF;
    steepness_ON(isession,:)  = steepnessON;
    rangeA_all(isession,:) = max(QtyA)*(rhoOFF + rhoON)/2;
    rangeB_all(isession,:) = max(QtyB);
    
    
    %%%% FILTER SESSIONS WITH NO SPLIT OR NOT SATURATED OFFERTYPES
    %%%%% CREATE QUALITY FILTER % 2 SPLITS TRIALS + GOOD SATURATION
    %%%%% (>90%)
    nlow_split = .1; nhigh_split = .9;
    nlow_satur = .1; nhigh_satur = .9;
    ordermerge = zeros(size(order)); % AB BA trials are merged
    jj = ii & iOF;
    [tableOFF, binosizeOFF] = gettable(QtyA(jj,1),QtyB(jj,1),ordermerge(jj,1),chosenID(jj,1));
    jj = ii & iON;
    [tableON,  binosizeON ] = gettable(QtyA(jj,1),QtyB(jj,1),ordermerge(jj,1),chosenID(jj,1));    
    
    if 1 % FILTER ON ALL TRIALS
        % % stim OFF
        % % based on BAratio
        % uni_logBA = unique(tableOFF.offer);
        % nuni = size(uni_logBA,1);
        % BpercOFF = [];
        % for iuni = 1:nuni
        %     ind = ismember(tableOFF.offer,uni_logBA(iuni),'rows');
        %     BpercOFF(iuni,1) = sum(tableOFF.choice(ind))./sum(binosizeOFF(ind));
        % end
        % based on offer type
        BpercOFF = tableOFF.choice./binosizeOFF;
        t1OFF = sum(BpercOFF>nlow_split & BpercOFF<nhigh_split); % SPLIT TRIALS
        t2OFF = sum(BpercOFF<=nlow_satur)>=1 | sum(BpercOFF>=nhigh_satur)>=1; % SATURATED TRIALS
        % t2OFF = sum(BpercOFF<=nlow_satur)>=1; % SATURATED TRIALS
        % % stim ON
        % % based on BAratio
        % uni_logBA = unique(tableON.offer);
        % nuni = size(uni_logBA,1);
        % BpercON = [];
        % for iuni = 1:nuni
        %     ind = ismember(tableON.offer,uni_logBA(iuni),'rows');
        %     BpercON(iuni,1) = sum(tableON.choice(ind))./sum(binosizeON(ind));
        % end
        % based on offer type
        BpercON = tableON.choice./binosizeON;
        t1ON = sum(BpercON>nlow_split & BpercON<nhigh_split); % SPLIT TRIALS
        t2ON = sum(BpercON<=nlow_satur)>=1 | sum(BpercON>=nhigh_satur)>=1; % SATURATED TRIALS
        % t2ON = sum(BpercON<=nlow_satur)>=1; % SATURATED TRIALS
        %
        isplit(isession,1) = t1OFF>=1 & t1ON>=1 ;
        isat(isession,1) = t2OFF>=1 & t2ON>=1 ;        
    end
    
end % for isession

% save data to info
info.relative_value_OFF = relative_value_OFF;
info.relative_value_ON = relative_value_ON;
info.relative_value = relative_value;
info.steepnessOFF = steepness_OFF;
info.steepnessON = steepness_ON;
info.rangeA = rangeA_all;
info.rangeB = rangeB_all;
info.isat=isat;
info.isplit=isplit;
end


% % %
% function 3
% % %
function [table02,binosize] = gettable(QtyA,QtyB,order,chosenID)

alltrials = [QtyA,QtyB,order];
[trialtypes,~,~] = unique(alltrials,'rows');
ntrialtypes = size(trialtypes,1);

table01 = [];
for itrialtype = 1:ntrialtypes
    trialtype = trialtypes(itrialtype,:);
    table01(itrialtype,1:3) = trialtype;
    %
    ind_trialtypes = ismember(alltrials,trialtype,'rows');
    ntrials = sum(ind_trialtypes);
    nchBtrials = sum(chosenID(ind_trialtypes,1));
    %
    table01(itrialtype,4) = nchBtrials;
    table01(itrialtype,5) = ntrials;
end

logB_A = log(table01(:,2)./table01(:,1));
order2 = table01(:,3);
choice2 = table01(:,4);
binosize = table01(:,5);

table02=table(logB_A,order2,choice2,'VariableNames',{'offer','order','choice'});

end

