% plot_generator_allcurrent.m
% Figure 1 - Exp1, primary analysis: rho steepness and order bias, <= 15 uA
% Figure 2 - Exp1, primary analysis, all current levels

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

%
% highlight new sessions
monkey = 'Juan'; % 'Gervinho'; 'Juan'; 'both'
eval(['sessionlist_',monkey]);
newsessions = [];
for isess = 1:length(sessions)
    newsessions{isess,1} = sessions{isess}(2:end);
end
plotnewsessions = 0;

% example session
plotexamplesession = 1;
examplesession = {'180102b', '171205c'}; % '171205c' % for offer2
exampletitle = {'offer1, \leq15\muA', 'offer2, \leq15\muA'}; 

% monkey target
monkeyname = ''; % 'G', 'J', '' (both)

doremoveoutlier = 0;
% different definition of hysteresis
def_hyst = 0; % '0': hysteresis and bias are defined as a2; '1': defined as a2/a1; '2': defined as 2rho*a2/a1

% %
% % FIGURE 1
% %

% define analysis type and figure axis labels
%  ylab_f1={ '0\muA offer1'; '0\muA offer2';};
ylab_f1={ '\leq15\muA offer1'; '\leq15\muA offer2';};
% ylab_f1={ '50\muA offer1'; '50\muA offer2';};
% ylab_f1={ '\geq100\muA offer1'; '\geq100\muA offer2';};
% current levels - figure 1 only analyzes <=15uA
limlow1  = 1;
limhigh1 = 16;
limlow2  = 26;
limhigh2 = 51;
limhigh3 = 99;

% %
f(1)=figure(1);
set(gcf,'position',[110 65 1250 1250], 'PaperPositionMode','auto')
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
%
nplot=numel(paroff);
%
for k=1:3 % 3: three analysis target - relative value, steepness and order bias
% for k=2:2 % 3: three analysis target - relative value, steepness and order bias
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
        if plotnewsessions
           ind_newsess =  ismember(info.session(mask),newsessions);
        end
        % load data accordingly
		try
			switch k
				case 1
					tit={'Rho'};
					data=[info.relative_value_OFF info.relative_value_ON];
				case 3
					tit={'Steepness'};
					data=[info.steepness_OFF info.steepness_ON ];
				case 2
					tit={'Order bias'};
                    if def_hyst == 0
                        data=[info.order_bias_OFF info.order_bias_ON ];
                    elseif def_hyst == 1
                        data=[info.order_bias_OFF info.order_bias_ON ]./[info.steepness_OFF info.steepness_ON ];
                    elseif def_hyst == 2
                        data=2.*[info.relative_value_OFF info.relative_value_ON ].*...
                                [info.order_bias_OFF info.order_bias_ON ]./[info.steepness_OFF info.steepness_ON ];
                    end
			end
		catch
			[info,table01] = probitLink_regression(alldata, info);
			switch k
				case 1
					tit={'Rho'};
					data=[ info.relative_value_OFF info.relative_value_ON ];
				case 3
					tit={'Steepness'};
					data=[info.steepness_OFF info.steepness_ON ];
				case 2
					tit={'Order bias'};
                    if def_hyst == 0
                        data=[info.order_bias_OFF info.order_bias_ON ];
                    elseif def_hyst == 1
                        data=[info.order_bias_OFF info.order_bias_ON ]./[info.steepness_OFF info.steepness_ON ];
                    elseif def_hyst == 2
                        data=2.*[info.relative_value_OFF info.relative_value_ON ].*...
                                [info.order_bias_OFF info.order_bias_ON ]./[info.steepness_OFF info.steepness_ON ];
                    end
			end
        end       
		% find the limits for two axes for plotting
		fact=1;
		for i=1:size(data,1)
			limi(i,1)=min(min(data))*fact; limi(i,2)=max(max(data))*fact;
		end
		% only looks at target sessions
		data=data(mask,:);
        %
        % remove outlier
        if doremoveoutlier
            K = 3;
            dataXY = [data(:,1);data(:,2)];
            OUT_XX = (data(:,1)>(mean(dataXY)+std(dataXY)*K)) | (data(:,1)<(mean(dataXY)-std(dataXY)*K));  
            OUT_YY = (data(:,2)>(mean(dataXY)+std(dataXY)*K)) | (data(:,2)<(mean(dataXY)-std(dataXY)*K));  
            %
            OUT = OUT_XX | OUT_YY;
            %
            data = data(~OUT,:);        
        end
        %
		% plot
		ax{n,k}=subplot(nplot,3,k+(n-1)*3);
		hold on;
		% add ellipse
		eli=error_ellipse(cov(data),'mu',nanmean(data),'conf',.9);
		eli.LineWidth=2; eli.Color=[0.7 0.7 0.7];
		% add data points
		plot(data(:,1),data(:,2),'ko','MarkerSize',8);
        if plotexamplesession
            plot(data(ind_exam,1),data(ind_exam,2),'ko','MarkerSize',8,'MarkerFaceColor',[0.4 0.6 0.4]);
        end
        if plotnewsessions
            plot(data(ind_newsess,1),data(ind_newsess,2),'ko','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.0]);
        end
        plot(limi(k,:),limi(k,:),'k--','HandleVisibility','off','LineWidth',1);
		axis ([limi(k,:), limi(k,:)])
		axis square
		xlabel('Stim OFF');
		title(ylab_f1{n});
		ylabel({tit{1} 'Stim ON'});
	end
end
%
% statistical analysis
for k=1:3
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
				tit={'Rho'};
				data=[ info.relative_value_OFF info.relative_value_ON ];
			case 3
				tit={'Steepness'};
				data=[info.steepness_OFF info.steepness_ON ];
			case 2
				tit={'Order bias'};
                if def_hyst == 0
                    data=[info.order_bias_OFF info.order_bias_ON ];
                elseif def_hyst == 1
                    data=[info.order_bias_OFF info.order_bias_ON ]./[info.steepness_OFF info.steepness_ON ];
                elseif def_hyst == 2
                    data=2.*[info.relative_value_OFF info.relative_value_ON ].*...
                            [info.order_bias_OFF info.order_bias_ON ]./[info.steepness_OFF info.steepness_ON ];
                end
		end
		% only looks at target sessions
		data=data(mask,:);
		%
        % remove outlier
        if doremoveoutlier
            K = 3;
            dataXY = [data(:,1);data(:,2)];
            OUT_XX = (data(:,1)>(mean(dataXY)+std(dataXY)*K)) | (data(:,1)<(mean(dataXY)-std(dataXY)*K));  
            OUT_YY = (data(:,2)>(mean(dataXY)+std(dataXY)*K)) | (data(:,2)<(mean(dataXY)-std(dataXY)*K));              
            %
            OUT = OUT_XX | OUT_YY;
            %
            data = data(~OUT,:);        
        end
        %
		[~, p_ABOF(k,n)]=ttest(data(:,1),data(:,2));
		[ p(k,n) ] = signrank(data(:,1),data(:,2));
		pos=ax{n,k}.Position; AX{n,k}=axes('Parent',f,'Position',[pos(1) pos(2)-.1 pos(3) .1 ]);
		text(0,0.8,['  n sessions= ' num2str(size(data,1))],'Color','k');
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


% %
% % FIGURE 2
% %

% define analysis type and figure axis labels
ylab_f2={ 'NoStim'; '\leq15muA'; '\leq15muA offer2'; '25\muA'; '25\muA Offer2'; '50\muA'; '50\muA Offer2'; '\geq100\muA'; '\geq100\muA offer2';};
% different current levels
limlow1  = 1;
limhigh1 = 16;
limlow2  = 25;
limhigh2 = 50;
limhigh3 = 100;
%
paroff=[0 1 2 1 2 1 2 1 2]; % 0: non stimulation control; 1: offer1 stimulation; 2: offer2 stimulation
parsup=[0 limlow1  limlow1  limlow2  limlow2  limhigh2 limhigh2 limhigh3 limhigh3];  % lower boundary to choose target current level
parinf=[1 limhigh1 limhigh1 limhigh2 limhigh2 limhigh3 limhigh3 Inf Inf ];         % upper boundary to choose target current level
%
for k=1:3 % 3: three analysis target - relative value, steepness and order bias
	for n=1:numel(paroff)
		% mask contains index for target sessions
        if ~isempty(monkeyname)
            mask_monkey = ismember(info.monkey,monkeyname);
        else
            mask_monkey = logical(ones(size(info.monkey)));
        end
		mask=(info.StimOffer(:,1)==paroff(n) & info.StimCurrent>=parsup(n) & info.StimCurrent<parinf(n)) ;
        mask = mask & mask_monkey; 
		try
			switch k
				case 1
					tit={'Rho'};
					data=[ info.relative_value_OFF info.relative_value_ON ];
				case 3
					tit={'Steepness'};
					data=[info.steepness_OFF info.steepness_ON ];
				case 2
					tit={'Order bias'};
                    if def_hyst == 0
                        data=[info.order_bias_OFF info.order_bias_ON ];
                    elseif def_hyst == 1
                        data=[info.order_bias_OFF info.order_bias_ON ]./[info.steepness_OFF info.steepness_ON ];
                    elseif def_hyst == 2
                        data=2.*[info.relative_value_OFF info.relative_value_ON ].*...
                                [info.order_bias_OFF info.order_bias_ON ]./[info.steepness_OFF info.steepness_ON ];
                    end
			end
		catch
			[info,table01] = probitLink_regression(alldata, info);
			switch k
				case 1
					tit={'Rho'};
					data=[ info.relative_value_OFF info.relative_value_ON ];
				case 3
					tit={'Steepness'};
					data=[info.steepness_OFF info.steepness_ON ];
				case 2
					tit={'Order bias'};
                    if def_hyst == 0
                        data=[info.order_bias_OFF info.order_bias_ON ];
                    elseif def_hyst == 1
                        data=[info.order_bias_OFF info.order_bias_ON ]./[info.steepness_OFF info.steepness_ON ];
                    elseif def_hyst == 2
                        data=2.*[info.relative_value_OFF info.relative_value_ON ].*...
                                [info.order_bias_OFF info.order_bias_ON ]./[info.steepness_OFF info.steepness_ON ];
                    end
			end
		end
		% only looks at target sessions
		data=data(mask,:);
		mdata(k,n) = mean(diff(data'));
		stdata(k,n) =std(diff(data'))./sqrt(sum(mask));
		[~,pval_ttest(k,n)] = ttest(diff(data'));
		pval_wil(k,n) = signrank(diff(data'));
		if isnan(mdata(k,n))
			keyboard
		end
		
	end
end
%
% plot
tit={'Relative value \rho' , 'Steepness \eta' , 'Order bias \epsilon'};
f(2) = figure(2);
set(gcf,'position',[110 65 1450 850], 'PaperPositionMode','auto')
for k=1:3
	mdata_off1=mdata(k,paroff==1 | paroff==0);
	mdata_off2=mdata(k,paroff==2 | paroff==0);
	stddata_off1=stdata(k,paroff==1 | paroff==0);
	stddata_off2=stdata(k,paroff==2 | paroff==0);
	pttest_off1=pval_ttest(k,paroff==1 | paroff==0);
	pttest_off2=pval_ttest(k,paroff==2 | paroff==0);
	pwil_off1=pval_wil(k,paroff==1 | paroff==0);
	pwil_off2=pval_wil(k,paroff==2 | paroff==0);
	%
	subplot(1,3,k)
	plot([0 numel(parsup)/2+1], [0 0],'k--','HandleVisibility','off'); hold on
	h=errorbar([1:numel(parsup)/2+1],mdata_off1,stddata_off1,'o','Linewidth',1,'HandleVisibility','off');
	h.MarkerEdgeColor = [0.0 0.2 0.8];
	h.MarkerFaceColor = [0.0 0.2 0.8];
	h.MarkerSize = 10;
	h = errorbar([1:numel(parsup)/2+1],mdata_off2,stddata_off2,'o','Linewidth',1,'HandleVisibility','off');
	h.MarkerEdgeColor = [1.0 0.4 0.0];
	h.MarkerFaceColor = [1.0 0.4 0.0];
	h.MarkerSize = 10;
	pl1 = plot([1:numel(parsup)/2+1],mdata_off1,'-','Linewidth',1);
	pl1.Color = [0.0 0.2 0.8];
	pl2 = plot([1:numel(parsup)/2+1],mdata_off2,'-','Linewidth',1);
	pl2.Color = [1.0 0.4 0.0];
	title(tit{k});
	legend([ pl1 pl2],{'Offer1', 'Offer2'},'Location','northwest');
	set(gca,'Xtick',[1:numel(parsup)/2+1],'XtickLabel',{'NoStim', '\leq15\muA', '25\muA', '50\muA', '\geq100\muA', ' '});
	xlabel('current level (\muA)')
    if ~isempty(monkeyname)
        pthres = 0.05;
    else
        pthres = 0.01;
    end
	switch k
		case 1
			ylim([-.5 .5])
			stars = [1:numel(parsup)/2+1];
			ind_signi1=pwil_off1<pthres;
			plot(stars(ind_signi1),ones(size(stars(ind_signi1))).*0.35,'k*');
			ind_signi2=pwil_off2<pthres;
			plot(stars(ind_signi2),ones(size(stars(ind_signi2))).*-0.35,'k*');
			%
			ylabel('\rhostimON - \rhostimOFF');
		case 3
			ylim([-1.5 1.5])
			stars = [1:numel(parsup)/2+1];
			ind_signi1=pwil_off1<pthres;
			plot(stars(ind_signi1),ones(size(stars(ind_signi1))).*1.1,'k*');
			ind_signi2=pwil_off2<pthres;
			plot(stars(ind_signi2),ones(size(stars(ind_signi2))).*-1.1,'k*');
			%
			ylabel('\etastimON - \etastimOFF');
		case 2
			ylim([-.5 .5])
			stars = [1:numel(parsup)/2+1];
			ind_signi1=pwil_off1<pthres;
			plot(stars(ind_signi1),ones(size(stars(ind_signi1))).*0.47,'k*');
			ind_signi2=pwil_off2<pthres;
			plot(stars(ind_signi2),ones(size(stars(ind_signi2))).*-0.47,'k*');
			%
			ylabel('\epsilonstimON - \epsilonstimOFF');
	end
	xlim([0.5 5.5])
	axis square
	box off
end
if ~isempty(monkeyname)
    suptitle({['Figure 2'];['monkey ',monkeyname]})
else
    suptitle({['Figure 2'];['monkey both']})
end


% %
% % FIGURE 3: example session
% %
if plotexamplesession
try
    ind_examsess = ismember(info.session, examplesession);
    betaOFF = info.beta_OFF(ind_examsess,:);
    betaON  = info.beta_ON(ind_examsess,:);
    table01s_ABOF = {table01.ABOFF{ind_examsess}};
    table01s_ABON = {table01.ABON{ind_examsess}};
    table01s_BAOF = {table01.BAOFF{ind_examsess}};
    table01s_BAON = {table01.BAON{ind_examsess}};
    rho_OFF = info.relative_value_OFF(ind_examsess);
    rho_ON  = info.relative_value_ON(ind_examsess);
    steepness_OFF = info.steepness_OFF(ind_examsess);
    steepness_ON  = info.steepness_ON(ind_examsess);
    orderbias_OFF = info.order_bias_OFF(ind_examsess);
    orderbias_ON = info.order_bias_ON(ind_examsess);
catch
    [info, table01] = probitLink_regression(alldata, info);
    betaOFF = info.beta_OFF(ind_examsess,:);
    betaON  = info.beta_ON(ind_examsess,:);
    table01s_ABOF = {table01.ABOFF{ind_examsess}};
    table01s_ABON = {table01.ABON{ind_examsess}};
    table01s_BAOF = {table01.BAOFF{ind_examsess}};
    table01s_BAON = {table01.BAON{ind_examsess}};
    rho_OFF = info.relative_value_OFF(ind_examsess);
    rho_ON  = info.relative_value_ON(ind_examsess);
    steepness_OFF = info.steepness_OFF(ind_examsess);
    steepness_ON  = info.steepness_ON(ind_examsess);
    orderbias_OFF = info.order_bias_OFF(ind_examsess);
    orderbias_ON = info.order_bias_ON(ind_examsess);
end
nexamplesession = sum(ind_examsess);
for iexamsess = 1: nexamplesession
    f(2+iexamsess) = figure(2+iexamsess);
    set(gcf,'position',[110 65 850 850], 'PaperPositionMode','auto')
    set(gca,'fontsize',9, 'fontweight','normal')
    %
    table01_ABOF = table01s_ABOF{iexamsess};
    table01_ABON = table01s_ABON{iexamsess};
    table01_BAOF = table01s_BAOF{iexamsess};
    table01_BAON = table01s_BAON{iexamsess};
    %
    % AB stimOFF 
    p_ABOF = plot_psyphy_spec(table01_ABOF, [betaOFF(iexamsess,1)+betaOFF(iexamsess,3),betaOFF(iexamsess,2)], '^', [1.0 0.2 0.2], 10, 4);
    % BA stimOFF 
    p_BAOF = plot_psyphy_spec(table01_BAOF, [betaOFF(iexamsess,1)-betaOFF(iexamsess,3),betaOFF(iexamsess,2)], '<', [0.0 1.0 1.0], 10, 4);
    % AB stimON
    p_ABON = plot_psyphy_spec(table01_ABON, [betaON(iexamsess,1)+betaON(iexamsess,3),betaON(iexamsess,2)], 'V', [0.6 0.0 0.0], 10, 4); 
    % BA stimON
    p_BAON = plot_psyphy_spec(table01_BAON, [betaON(iexamsess,1)-betaON(iexamsess,3),betaON(iexamsess,2)], '>', [0.0 0.0 0.6], 10, 4);
    %   
    lgd = legend([p_ABOF, p_BAOF, p_ABON, p_BAON],'AB stimOFF','BA stimOFF', 'AB stimON','BA stimON', 'Location','NorthWest');
    lgd.FontSize = 14;           
    xtickangle(45);
    axis square
    xlabel('log(qB/qA)');
    ylabel('B choices (%)');
    title({exampletitle{iexamsess}; examplesession{iexamsess}});
    set(gca,'FontSize',18);
    text(1.0,0.5,{['\rho_{OFF} = ',num2str(rho_OFF(iexamsess),'%.2f')];...
                  ['\rho_{ON}  = ',num2str(rho_ON(iexamsess),'%.2f')];...
                  ['\eta_{OFF} = ',num2str(steepness_OFF(iexamsess),'%.2f')];...
                  ['\eta_{ON}  = ',num2str(steepness_ON(iexamsess),'%.2f')];...
                  ['\epsilon_{OFF} = ',num2str(orderbias_OFF(iexamsess),'%.2f')];...
                  ['\epsilon_{ON}  = ',num2str(orderbias_ON(iexamsess),'%.2f')];...
                 },'fontsize', 15, 'color', 'k');  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % %
% function 1
% % %
function [info, table01] = probitLink_regression(alldata, info)
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
	
    % make and save table01
    iAB = order == 1;
    iBA = order ==-1;
    table01.ABOFF{isession,:} = get_table01(QtyA(iOF & iAB,:), QtyB(iOF & iAB,:), chosenID(iOF & iAB,:));
	table01.BAOFF{isession,:} = get_table01(QtyA(iOF & iBA,:), QtyB(iOF & iBA,:), chosenID(iOF & iBA,:));
    table01.ABON{isession,:}  = get_table01(QtyA(iON & iAB,:), QtyB(iON & iAB,:), chosenID(iON & iAB,:));
	table01.BAON{isession,:}  = get_table01(QtyA(iON & iBA,:), QtyB(iON & iBA,:), chosenID(iON & iBA,:));
    
	% save data
	relative_value_OFF(isession,:) = rhoOFF;
	relative_value_ON(isession,:) = rhoON;
	steepness_OFF(isession,:) = steepnessOFF;
	steepness_ON(isession,:) = steepnessON;
	order_bias_OFF(isession,:) = orderbiasOFF;
	order_bias_ON(isession,:) = orderbiasON;
    beta_OFF(isession,:) = betaOFF;
    beta_ON(isession,:)  = betaON;
    
end % for isession

% save data to info
info.relative_value_OFF = relative_value_OFF;
info.relative_value_ON = relative_value_ON;
info.steepness_OFF = steepness_OFF;
info.steepness_ON = steepness_ON;
info.order_bias_OFF = order_bias_OFF;
info.order_bias_ON = order_bias_ON;
info.beta_OFF = beta_OFF;
info.beta_ON  = beta_ON;
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


% % %
% function 3
% % %
function [po] = plot_psyphy_spec(table01, beta, dots, clr, dotssize, linewidth)

hold on

table01_all = table01;

table02 = table01(table01(:,1) & table01(:,2),:);
xx = log(table02(:,2)./table02(:,1));
yy = table02(:,3);

%domain for plot
x = min(xx)-log(2)/0.75 : .05 : max(xx)+log(2)/0.75;
%plot fitted sigmoid
y_choicefit = 1./(1+exp(-(beta(1)+beta(2)*x)));
plot(x,y_choicefit,'-','color',clr,'linewidth',linewidth);
%plot choice pattern
po = plot(xx, yy, dots, 'color',clr, 'markersize',dotssize, 'MarkerFaceColor',clr); %,'markerfacecolor',clr)

%add forced choices
forcedAtab = table01_all( table01_all(:,1) & ~table01_all(:,2),:);
nfA = size(forcedAtab,1);
xx_forcedA = min(xx)-log(2)*(1:nfA);	xx_forcedA = sort(xx_forcedA)';
plot(xx_forcedA, forcedAtab(:,3)', dots, 'color',clr, 'markersize',dotssize, 'MarkerFaceColor',clr)

forcedBtab = table01_all(~table01_all(:,1) &  table01_all(:,2),:);
nfB = size(forcedBtab,1);
xx_forcedB = max(xx)+log(2)*(1:nfB);	xx_forcedB = sort(xx_forcedB)';
plot(xx_forcedB, forcedBtab(:,3)', dots, 'color',clr, 'markersize',dotssize, 'MarkerFaceColor',clr)

%cosmetics
[~,ind,~] = unique(xx);	%remove doubles in xx
xxx = xx(ind);

%xlabels
xlab = cell(1,nfA);
for ifA = 1:nfA
	xlab{ifA} = [num2str(forcedAtab(ifA,2)),':',num2str(forcedAtab(ifA,1))];
end
for i = 1:size(xxx,1)
	xlab{nfA+i} = [num2str(table02(ind(i),2)),':',num2str(table02(ind(i),1))];
end
for ifB = 1:nfB
	xlab{nfA+i+ifB} = [num2str(forcedBtab(ifB,2)),':',num2str(forcedBtab(ifB,1))];
end

%add forced choices
if ~isempty(xx_forcedB) & ~isempty(xx_forcedA) 
    xxx = [xx_forcedA;xxx;xx_forcedB];    
    set(gca,'xlim',[min(xxx)-log(2)/5 max(xxx)+log(2)/5])
    set(gca,'xtick',xxx,'xticklabel',xlab)
    set(gca,'ylim',[0,1]);
    set(gca,'ytick',0:.25:1,'yticklabel',0:25:100)
end

end
