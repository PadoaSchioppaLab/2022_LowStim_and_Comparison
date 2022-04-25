% Fig1and2_plot_generator_allcurrent_RT_errate.m
% Scatter plot to show RT and error rate effects

% Authors: Sebastien Ballesta & Weikang Shi

% Copyright: Camillo Padoa-Schioppa lab, Washington Univerty in St. Louis

close all
clearvars

% % load data and info files
% % if use complete_dataset.xlsx
[~,~,Exp1_info] = xlsread('session_info.xlsx', 1); 
[~,~,Exp1_data] = xlsread('session_detailed_data.xlsx', 1); 
info=cell2table(Exp1_info(2:end,:),'VariableNames',Exp1_info(1,:));
alldata=cell2table(Exp1_data(2:end,:),'VariableNames',Exp1_data(1,:));


% example session
plotexamplesession = 0;
examplesession = '171204d'; % '171204d', '171205c'
exampletitle = 'offer2, \leq15\muA'; 

% monkey target
monkeyname = ''; % 'G', 'J', '' (both)

% %
% % FIGURE 1
% % scatter plot: raction time; error rates
%
% define analysis type and figure axis labels
%  ylab_f3={ '0\muA offer1'; '0\muA offer2';};
ylab_f3={ '\leq15\muA offer1'; '\leq15\muA offer2';};
% ylab_f3={ '50\muA offer1'; '50\muA offer2';};
% ylab_f3={ '\geq100\muA offer1'; '\geq100\muA offer2';};
% current levels - figure 1 only analyzes <=15uA
limlow1  = 1;
limhigh1 = 16;
limlow2  = 26;
limhigh2 = 51;
limhigh3 = 99;

% %
f(1)=figure(1);
% set(gcf,'position',[110 65 1250 1250], 'PaperPositionMode','auto')
set(gcf,'position',[110 65 1450 750], 'PaperPositionMode','auto')
hold on;
% % 0uA
% paroff=[  0 0  ]; % offer1 or offer2
% parinf=[  limlow1 limlow1 ];     % use  1uA as upper boundary to choose target sessions
% parsup=[  -1  -1  ];             % use -1uA as lower boundary to choose target sessions       
% <=15uA
paroff=[  1 2  ]; % offer1 or offer2
parinf=[  limhigh1 limhigh1 ];     % use 16uA as upper boundary to choose target sessions
parsup=[  limlow1  limlow1  ];     % use 1uA  as lower boundary to choose target sessions       
% % 50uA
% paroff=[  1 2  ]; % offer1 or offer2
% parinf=[  limhigh2 limhigh2 ];     % use 51uA as upper boundary to choose target sessions
% parsup=[  limlow2  limlow2  ];     % use 26uA as lower boundary to choose target sessions   
% % >=100uA
% paroff=[  1 2  ]; % offer1 or offer2
% parinf=[  inf inf ];                 % use inf as upper boundary to choose target sessions
% parsup=[  limhigh3  limhigh3  ];     % use 100uA  as lower boundary to choose target sessions      

nplot=numel(paroff);
%
for k=1:2 % 2: two analysis target - raction time; error rate
	for n=1:numel(paroff)
		% mask contains index for target sessions
        if ~isempty(monkeyname)
            mask_monkey = ismember(info.monkey,monkeyname);
        else
            mask_monkey = logical(ones(size(info.monkey)));
        end
		mask = (info.StimOffer(:,1)==paroff(n) & info.StimCurrent>parsup(n) & info.StimCurrent<=parinf(n)) ;
        mask = mask & mask_monkey; 
        if plotexamplesession
            ind_exam = ismember(info.session(mask),examplesession);
        end
        % load data accordingly
		try
			switch k				
                case 1
					tit={'Error rates'};
					data=[ info.errRate_OFF info.errRate_ON ];
                case 2
					tit={'RTime'};
					data=[ info.RT_OFF info.RT_ON ];               
			end
		catch
			[info,~] = probitLink_regression_RT_errate(alldata, info);
			switch k
                case 1
					tit={'Error rates'};
					data=[ info.errRate_OFF info.errRate_ON ];
                case 2
					tit={'RTime'};
					data=[ info.RT_OFF info.RT_ON ];             
			end
		end
		% find the limits for two axes for plotting
		fact=1;
		for i=1:size(data,1)
			limi(i,1)=min(min(data))*fact; limi(i,2)=max(max(data))*fact;
		end
		% only looks at target sessions
		data=data(mask,:);
		% plot
		ax{n,k}=subplot(nplot,6,k+(n-1)*6);
%         ax{n,k}=subplot(nplot,2,n+(k-1)*2);
		hold on;
		% add ellipse
        try
		eli=error_ellipse(cov(data),'mu',nanmean(data),'conf',.9);
		eli.LineWidth=2; eli.Color=[0.7 0.7 0.7];
        end
        % add data points
		plot(data(:,1),data(:,2),'ko','MarkerSize',8);
        if plotexamplesession
            plot(data(ind_exam,1),data(ind_exam,2),'ko','MarkerSize',8,'MarkerFaceColor',[0.4 0.6 0.4]);
        end
        try
		plot(limi(k,:),limi(k,:),'k--','HandleVisibility','off','LineWidth',1);
		axis ([limi(k,:), limi(k,:)])
        end
        axis square
		xlabel('Stim OFF');
		title(ylab_f3{n});
		ylabel({tit{1} 'Stim ON'});
	end
end
%
% statistical analysis
for k=1:2
	for n=1:numel(paroff)
        if ~isempty(monkeyname)
            mask_monkey = ismember(info.monkey,monkeyname);
        else
            mask_monkey = logical(ones(size(info.monkey)));
        end
		mask=(info.StimOffer(:,1)==paroff(n) & info.StimCurrent>=parsup(n) & info.StimCurrent<parinf(n));
        mask = mask & mask_monkey; 
		switch k		
            case 1
                tit={'Error rates'};
                data=[ info.errRate_OFF info.errRate_ON ];
            case 2
                tit={'RTime'};
                data=[ info.RT_OFF info.RT_ON ];           
		end
		% only looks at target sessions
		data=data(mask,:);
		%
		[~, p_ABOF(k,n)]=ttest(data(:,1),data(:,2));
		[ p(k,n) ] = signrank(data(:,1),data(:,2));
		pos=ax{n,k}.Position; AX{n,k}=axes('Parent',f(1),'Position',[pos(1) pos(2)-.1 pos(3) .1 ]);
		text(0,0.8,['  n sessions= ' num2str(sum(mask))],'Color','k');
		if p(k,n)>.05
			text(0,.6,['   Wilcoxon: p = ', num2str(p(k,n),'%.4f') ])
		else
			text(0,.6,['   Wilcoxon: p = ', num2str(p(k,n),'%.4f') ], 'FontWeight', 'bold')
		end
		if p_ABOF(k,n)>.05
			text(0,.4 ,['  ttest: p = ', num2str(p_ABOF(k,n),'%.4f') ])
		else
			text(0,.4, ['   ttest: p = ', num2str(p_ABOF(k,n),'%.4f') ], 'FontWeight', 'bold')
		end
		AX{n,k}.Visible ='off';
	end
end
if ~isempty(monkeyname)
    suptitle({['Figure 1'];['monkey ',monkeyname]})
else
    suptitle({['Figure 1'];['monkey both']})
end

         
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % %
% function 1
% % %
function [info, table01] = probitLink_regression_RT_errate(alldata, info)
%
% This function is to do logistic regression with "probit" link function
% in this regression:
% X = b0 + b1 * log(#B/#A) + b2 *(A_right - A_left)
% regression is done on stimOFF and stimON trials seperately
%
nsessions = size(info,1);
%
for isession = 1:nsessions
	session = info.session{isession};
    monkey = info.monkey{isession};
    if isequal(monkey,'G')
        monkey_full = 'Gervinho';
    elseif isequal(monkey,'J')
        monkey_full = 'Juan';
    end
    
    stimlateral = info.StimLateral(isession); % 1 unilateral; 2 bilateral
    
	% load raw data from alldata
	ind_trials = ismember(alldata.Session, ['ST',session]);
	QtyA = alldata.QtyA(ind_trials,:);
	QtyB = alldata.QtyB(ind_trials,:);
	stim = alldata.Stim(ind_trials,:);
	chosenID = alldata.ChosenID(ind_trials,:);
	order = -alldata.Order(ind_trials,:);
    sides =  alldata.Side(ind_trials,:);
    [table01] = get_table01(QtyA, QtyB, chosenID);
	%
	% remove forced choices
	ii = ~prod([QtyA,QtyB],2)==0;
	%
	iOF = stim==-1;
	iON = stim== 1;
	      
    % Reaction time 
    RTs_all = alldata.Rxtime(ind_trials,:);
    % for stimOFF
	jj = ii & iOF;
    RTOFF = nanmean(RTs_all(jj,:));
    % for stimON
	jj = ii & iON;
    RTON  = nanmean(RTs_all(jj,:));
    
    
    % Error rate (percentage correct)
    disp(['error rate: session ',monkey,session]);
    clear parsession
    try
        filename = ['Z:\Backups\Analysis6_old_finalBackup_210509\Experiments\Stimulation\Sessions\ST_',monkey,session,'.m'];   
        run(filename);
        cellname = num2str([parsession.clusters(1,1)*10+parsession.clusters(1,2)]);
        filename = ['Z:\Backups\Analysis6_old_finalBackup_210509\Experiments\Stimulation\Data\',monkey_full,'\',session,'\',monkey,session,cellname,'_data.mat'];
        load(filename);
    catch
        addpath C:\Experiments\SOstim\Sessions\Juan
        addpath C:\Experiments\SOstim\Data
        filename = ['ST_',monkey,session,'.m'];   
        run(filename);
        cellname = num2str([parsession.clusters(1,1)*10+parsession.clusters(1,2)]);
        filename = ['C:\Experiments\SOstim\Data\',monkey_full,'\',session(1:end-1),'\',monkey,session,cellname,'_data.mat'];
        load(filename);
    end
    % cut off trials when the error rate is too high
    errtrialnum = [trialRecordError.trialNumber]';
    trutrialnum = [trialRecord.trialNumber]';
    alltrialnum = [[errtrialnum; trutrialnum], [ones(size(errtrialnum)); zeros(size(trutrialnum))]];
    alltrialnum = sortrows(alltrialnum,1);
    windowSize = round(size(alltrialnum,1)/20); 
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    alltrialnum(:,3) = filter(b,a,alltrialnum(:,2));
    cutoff_trialnum = alltrialnum(find(alltrialnum(:,3)>0.6,1),1);
    if isempty(cutoff_trialnum) % |1        
        ind_good_tru = logical(ones(size(trutrialnum)));
        ind_good_err = logical(ones(size(errtrialnum)));
    else
        ind_good_tru = trutrialnum < cutoff_trialnum;
        ind_good_err = errtrialnum < cutoff_trialnum;
    end
    %
    trialRecord_02 = trialRecord(ind_good_tru);
    trialRecordError_02 = trialRecordError(ind_good_err);
    ncorTrials = length(trialRecord_02);
    nerrTrials = length(trialRecordError_02);
    errtypes = [trialRecordError_02.errortype]';
    stim_cor = [];
    stim_err = [];
    ind_forced_cor = [];
    ind_forced_err = [];
    for itrial_cor = 1:ncorTrials
        stim_cor(itrial_cor,:) = abs(trialRecord_02(itrial_cor).currentOffers(1).stim); % 1: OFF 2:ON
        ind_forced_cor(itrial_cor,:) = trialRecord_02(itrial_cor).currentOffers(1).quantity == 0 | trialRecord_02(itrial_cor).currentOffers(2).quantity == 0;       
    end
    for itrial_err = 1:nerrTrials
        stim_err(itrial_err,:) = abs(trialRecordError_02(itrial_err).currentOffers(1).stim); % 1: OFF 2:ON
        ind_forced_err(itrial_err,:) = trialRecordError_02(itrial_err).currentOffers(1).quantity == 0 | trialRecordError_02(itrial_err).currentOffers(2).quantity == 0;
    end
    % remove forced choice
    stim_cor(logical(ind_forced_cor)) = [];
    stim_err(logical(ind_forced_err)) = [];
    errtypes(logical(ind_forced_err)) = [];
    %
    errRateOFF = sum(stim_err==1)./(sum(stim_err==1) + sum(stim_cor==1));
    errRateON  = sum(stim_err==2)./(sum(stim_err==2) + sum(stim_cor==2));
    %
    err1RateOFF = sum(stim_err==1 & errtypes ==1)./(sum(stim_err==1) + sum(stim_cor==1));
    err1RateON  = sum(stim_err==2 & errtypes ==1)./(sum(stim_err==2) + sum(stim_cor==2));
    %
    err2RateOFF = sum(stim_err==1 & errtypes ==2)./(sum(stim_err==1) + sum(stim_cor==1));
    err2RateON  = sum(stim_err==2 & errtypes ==2)./(sum(stim_err==2) + sum(stim_cor==2));
    %
    err3RateOFF = sum(stim_err==1 & errtypes ==3)./(sum(stim_err==1) + sum(stim_cor==1));
    err3RateON  = sum(stim_err==2 & errtypes ==3)./(sum(stim_err==2) + sum(stim_cor==2));
    %
    err4RateOFF = sum(stim_err==1 & errtypes ==4)./(sum(stim_err==1) + sum(stim_cor==1));
    err4RateON  = sum(stim_err==2 & errtypes ==4)./(sum(stim_err==2) + sum(stim_cor==2));
    %
    err5RateOFF = sum(stim_err==1 & errtypes ==5)./(sum(stim_err==1) + sum(stim_cor==1));
    err5RateON  = sum(stim_err==2 & errtypes ==5)./(sum(stim_err==2) + sum(stim_cor==2));
    %
    err6RateOFF = sum(stim_err==1 & errtypes ==6)./(sum(stim_err==1) + sum(stim_cor==1));
    err6RateON  = sum(stim_err==2 & errtypes ==6)./(sum(stim_err==2) + sum(stim_cor==2));
    
	% save data
    RT_OFF(isession,:) = RTOFF;
	RT_ON(isession,:) = RTON;
    errRate_OFF(isession,:) = errRateOFF;
    errRate_ON(isession,:) = errRateON;
    err1Rate_OFF(isession,:) = err1RateOFF;
    err1Rate_ON(isession,:) = err1RateON;
    err2Rate_OFF(isession,:) = err2RateOFF;
    err2Rate_ON(isession,:) = err2RateON;
    err3Rate_OFF(isession,:) = err3RateOFF;
    err3Rate_ON(isession,:) = err3RateON;
    err4Rate_OFF(isession,:) = err4RateOFF;
    err4Rate_ON(isession,:) = err4RateON;
    err5Rate_OFF(isession,:) = err5RateOFF;
    err5Rate_ON(isession,:) = err5RateON;
    err6Rate_OFF(isession,:) = err6RateOFF;
    err6Rate_ON(isession,:) = err6RateON;
    
end % for isession

% save data to info
info.RT_OFF = RT_OFF;
info.RT_ON = RT_ON;
info.errRate_OFF = errRate_OFF;
info.errRate_ON = errRate_ON;
info.err1Rate_OFF = err1Rate_OFF;
info.err1Rate_ON = err1Rate_ON;
info.err2Rate_OFF = err2Rate_OFF;
info.err2Rate_ON = err2Rate_ON;
info.err3Rate_OFF = err3Rate_OFF;
info.err3Rate_ON = err3Rate_ON;
info.err4Rate_OFF = err4Rate_OFF;
info.err4Rate_ON = err4Rate_ON;
info.err5Rate_OFF = err5Rate_OFF;
info.err5Rate_ON = err5Rate_ON;
info.err6Rate_OFF = err6Rate_OFF;
info.err6Rate_ON = err6Rate_ON;
end


% % %
% function 2
% % %
function [table01] = get_table01(QtyA, QtyB, ChosenID)
%
pairtrials = [QtyA, QtyB, ChosenID];
[offer, ~, groups] = unique(abs(pairtrials(:,1:2)),'rows');
noffs = size(offer,1);
choiz = pairtrials(:,3);
%
perc_B = nan(noffs,1);
Ntrials = nan(noffs,1);
for i = 1:noffs
	ind = find(groups== i);
	perc_B(i,1) = length(find(choiz(ind)== 1))/length(ind);	%perc of ch b
	Ntrials(i,1) = length(ind);
end
table01 = [offer, perc_B, Ntrials];
%
%sort choices
eps = 0.001;
aux = abs(table01) + eps;
[~, jnd] = sort(aux(:,2)./aux(:,1));
table01 = table01(jnd,:);

end
