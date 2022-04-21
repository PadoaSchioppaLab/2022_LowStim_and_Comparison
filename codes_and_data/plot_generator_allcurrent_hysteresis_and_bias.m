% plot_generator_allcurrent_hysteresis_and_bias.m
% Scatter plot to show RT and error rate effects

% Authors: Sebastien Ballesta & Weikang Shi

% Copyright: Camillo Padoa-Schioppa lab, Washington Univerty in St. Louis

close all
clearvars

doremoveoutlier = 0;
% different definition of hysteresis
def_hyst = 0; % '0': hysteresis and bias are defined as a2; '1': defined as a2/a1; '2': defined as 2rho*a2/a1

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
% % scatter plot: position choice hysteresis; juice choice hysteresis
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
set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
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
for k=1:4 % 4: three analysis target - juice choice hysteresis; position choice hysteresis; order choice hysteresis
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
                    tit={'Side bias'};
					data=[ info.sidebias_OFF info.sidebias_ON ];
                case 2
					tit={'Choice hysteresis: juice type'};
					data=[ info.chhy_OFF info.chhy_ON ];
				case 3
					tit={'Choice hysteresis: side'};
					data=[ info.pohy_OFF info.pohy_ON ];
				case 4
					tit={'Choice hysteresis: order'};
					data=[ info.orhy_OFF info.orhy_ON ];
			end
		catch
			[info] = probitLink_regression_hysteresis_and_bias(alldata, info, def_hyst);
			switch k
				 case 1
                    tit={'Side bias'};
					data=[ info.sidebias_OFF info.sidebias_ON ];
                case 2
					tit={'Choice hysteresis: juice type'};
					data=[ info.chhy_OFF info.chhy_ON ];
				case 3
					tit={'Choice hysteresis: side'};
					data=[ info.pohy_OFF info.pohy_ON ];
				case 4
					tit={'Choice hysteresis: order'};
					data=[ info.orhy_OFF info.orhy_ON ];
			end
		end
		% find the limits for two axes for plotting
		fact=1;
        limi = [];
		for i=1:size(data,1)
			% limi(i,1)=min(min(data(mask,:)))*fact; limi(i,2)=max(max(data(mask,:)))*fact;
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
            OUT = OUT_XX | OUT_YY;
            data = data(~OUT,:);        
        end
		% plot
% 		ax{n,k}=subplot(nplot,4,k+(n-1)*4);
        ax{n,k}=subplot(4,nplot,n+(k-1)*nplot);
		hold on;
		% add ellipse
		eli=error_ellipse(cov(data),'mu',nanmean(data),'conf',.9);
		eli.LineWidth=2; eli.Color=[0.7 0.7 0.7];
		% add data points
		plot(data(:,1),data(:,2),'ko','MarkerSize',8);
        if plotexamplesession
            plot(data(ind_exam,1),data(ind_exam,2),'ko','MarkerSize',8,'MarkerFaceColor',[0.4 0.6 0.4]);
        end
		plot(limi(k,:),limi(k,:),'k--','HandleVisibility','off','LineWidth',1);
        plot([0 0],limi(k,:),'k:','HandleVisibility','off','LineWidth',1);
        plot(limi(k,:),[0 0],'k:','HandleVisibility','off','LineWidth',1);
		axis ([limi(k,:), limi(k,:)])
		axis square
		xlabel('Stim OFF');
		title(ylab_f3{n});
		ylabel({tit{1} 'Stim ON'});
	end
end
%
% statistical analysis
for k=1:4
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
                    tit={'Side bias'};
					data=[ info.sidebias_OFF info.sidebias_ON ];
                case 2
					tit={'Choice hysteresis: juice type'};
					data=[ info.chhy_OFF info.chhy_ON ];
				case 3
					tit={'Choice hysteresis: side'};
					data=[ info.pohy_OFF info.pohy_ON ];
				case 4
					tit={'Choice hysteresis: order'};
					data=[ info.orhy_OFF info.orhy_ON ];
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
            OUT = OUT_XX | OUT_YY;
            data = data(~OUT,:);        
        end
		%
		[~, p_ABOF(k,n)]=ttest(data(:,1),data(:,2));
		[ p(k,n) ] = signrank(data(:,1),data(:,2));
        [~, p_tt_xx(k,n)]=ttest(data(:,1),0);
		[ p_wil_xx(k,n) ] = signrank(data(:,1),0);
        [~, p_tt_yy(k,n)]=ttest(data(:,2),0);
		[ p_wil_yy(k,n) ] = signrank(data(:,2),0);
        
		pos=ax{n,k}.Position; AX{n,k}=axes('Parent',f(1),'Position',[pos(1) pos(2)-.1 pos(3) .1 ]);
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
        % text(0,.2, ['x axis differs from 0: p(Wilcoxon) = ',num2str(p_wil_xx(k,n),'%.2f'),'; p(ttest) = ',num2str(p_tt_xx(k,n),'%.2f')])
        % text(0,.0, ['y axis differs from 0: p(Wilcoxon) = ',num2str(p_wil_yy(k,n),'%.2f'),'; p(ttest) = ',num2str(p_tt_yy(k,n),'%.2f')])
        
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
function [info] = probitLink_regression_hysteresis_and_bias(alldata, info, def_hyst)
%
% This function is to do logistic regression with "probit" link function
% in this regression:
% X = b0 + b1 * log(#B/#A) + b2 *effects
% regression is done on stimOFF and stimON trials seperately
%
nsessions = size(info,1);
%
for isession = 1:nsessions
	session = info.session{isession};
    
    stimlateral = info.StimLateral(isession); % 1 unilateral; 2 bilateral
    
	% load raw data from alldata
	ind_trials = ismember(alldata.Session, ['ST',session]);
	QtyA = alldata.QtyA(ind_trials,:);
	QtyB = alldata.QtyB(ind_trials,:);
	stim = alldata.Stim(ind_trials,:);
	chosenID = alldata.ChosenID(ind_trials,:);
	order = -alldata.Order(ind_trials,:);
    sides =  alldata.Side(ind_trials,:); % -1 A left; 1 A right
    
    % % %
	% remove forced choices
	ii = ~prod([QtyA,QtyB],2)==0;
	%
	iOF = stim==-1;
	iON = stim==1;
	%
	logB_A = log(abs(QtyB)./abs(QtyA));
    % % %
    
    
    % % % 
    % % %
    % side bias
    % effects = A_right - A_left
	% logistic regression
	% for stimOFF
	jj = ii & iOF;
	% [~ , ~, stats] = glmfit([logB_A(jj,1),sides(jj,1)], chosenID(jj,1), 'binomial', 'link','probit', 'constant','on');
	mdl = fitglm([logB_A(jj,1), sides(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
	betaOFF = mdl.Coefficients.Estimate;
	rhoOFF = exp(-betaOFF(1)/betaOFF(2));
	steepnessOFF = betaOFF(2);
    if def_hyst == 0
        sidebiasOFF = betaOFF(3);
    elseif def_hyst == 1
        sidebiasOFF = betaOFF(3)./betaOFF(2);
    elseif def_hyst == 2
        sidebiasOFF = 2*rhoOFF*betaOFF(3)./betaOFF(2);
    end
	if stimlateral == 2
        sidebiasOFF = -sidebiasOFF;
    end
    %      
	% for stimON
	jj = ii & iON;
	% [~ , ~, stats] = glmfit([logB_A(jj,1),sides(jj,1)], chosenID(jj,1), 'binomial', 'link','probit', 'constant','on');
	mdl = fitglm([logB_A(jj,1),sides(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
	betaON = mdl.Coefficients.Estimate;
	rhoON = exp(-betaON(1)/betaON(2));
	steepnessON = betaON(2);
	if def_hyst == 0
        sidebiasON = betaON(3);
    elseif def_hyst == 1
        sidebiasON = betaON(3)./betaON(2);
    elseif def_hyst == 2
        sidebiasON = 2*rhoON*betaON(3)./betaON(2);
    end
	if stimlateral == 2
        sidebiasON = -sidebiasON;
    end
    

    % % % 
    % % %
    % choice hysteresis: juice type 
    % define ChHyst - choice of previous trial: -1 for A, 1 for B, 0 for else
    trialNum = alldata.trialNumber(ind_trials,:);
    indA = chosenID==0;
    trA  = trialNum(indA,1);
    indB = chosenID==1;
    trB  = trialNum(indB,1);
    ChHyst = zeros(size(trialNum)); % choice of previous trial: -1 for A, 1 for B, 0 for else
    ChHyst(ismember(trialNum,trA+1)) = -1; 
    ChHyst(ismember(trialNum,trB+1)) = 1; 
	%
	% logistic regression
	% for stimOFF
	jj = ii & iOF;
	mdl = fitglm([logB_A(jj,1), ChHyst(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
	betaOFF = mdl.Coefficients.Estimate;
	rhoOFF = exp(-betaOFF(1)/betaOFF(2));
	steepnessOFF = betaOFF(2);
	if def_hyst == 0
        ChHystOFF = betaOFF(3);
    elseif def_hyst == 1
        ChHystOFF = betaOFF(3)./betaOFF(2);
    elseif def_hyst == 2
        ChHystOFF = 2*rhoOFF*betaOFF(3)./betaOFF(2);
    end
    %      
	% for stimON
	jj = ii & iON;
	mdl = fitglm([logB_A(jj,1), ChHyst(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
	betaON = mdl.Coefficients.Estimate;
	rhoON = exp(-betaON(1)/betaON(2));
	steepnessON = betaON(2);
	if def_hyst == 0
        ChHystON = betaON(3);
    elseif def_hyst == 1
        ChHystON = betaON(3)./betaON(2);
    elseif def_hyst == 2
        ChHystON = 2*rhoON*betaON(3)./betaON(2);
    end
    

    % % % 
    % % %
    % choice hysteresis: side
    % define PoHyst - -1 for juice A in the same location as last chosen side, 1 for juice B in the same location as last chosen side, 0 for else
    trialNum = alldata.trialNumber(ind_trials,:);
    ind_chosenL = (chosenID*2-1).*sides == 1;
    tr_chosenL  = trialNum(ind_chosenL,1); % trialnum of choosing left
    ind_chosenR = (chosenID*2-1).*sides == -1; 
    tr_chosenR  = trialNum(ind_chosenR,1); % trialnum of choosing right
    trialNum_AR = trialNum(sides== 1,:);
    trialNum_AL = trialNum(sides==-1,:);
    PoHyst = zeros(size(trialNum));
    PoHyst(ismember(trialNum_AL,tr_chosenL+1)) = -1; 
    PoHyst(ismember(trialNum_AR,tr_chosenR+1)) = -1; 
    PoHyst(ismember(trialNum_AL,tr_chosenR+1)) = 1; 
    PoHyst(ismember(trialNum_AR,tr_chosenL+1)) = 1; 
	%
	% logistic regression
	% for stimOFF
	jj = ii & iOF;
	mdl = fitglm([logB_A(jj,1), PoHyst(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
	betaOFF = mdl.Coefficients.Estimate;
	rhoOFF = exp(-betaOFF(1)/betaOFF(2));
	steepnessOFF = betaOFF(2);
	if def_hyst == 0
        PoHystOFF = betaOFF(3);
    elseif def_hyst == 1
        PoHystOFF = betaOFF(3)./betaOFF(2);
    elseif def_hyst == 2
        PoHystOFF = 2*rhoOFF*betaOFF(3)./betaOFF(2);
    end
    %
	% for stimON
	jj = ii & iON;
	mdl = fitglm([logB_A(jj,1), PoHyst(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
	betaON = mdl.Coefficients.Estimate;
	rhoON = exp(-betaON(1)/betaON(2));
	steepnessON = betaON(2);
	if def_hyst == 0
        PoHystON = betaON(3);
    elseif def_hyst == 1
        PoHystON = betaON(3)./betaON(2);
    elseif def_hyst == 2
        PoHystON = 2*rhoON*betaON(3)./betaON(2);
    end
    
    
    % % % 
    % % %
    % choice hysteresis: order
    % define OrHyst - -1 for juice A in the same location as last chosen order, 1 for juice B in the same location as last chosen order, 0 for else
    trialNum = alldata.trialNumber(ind_trials,:);
    ind_chosen1 = (chosenID*2-1).*order == -1;
    tr_chosen1  = trialNum(ind_chosen1,1);
    ind_chosen2 = (chosenID*2-1).*order == 1;
    tr_chosen2  = trialNum(ind_chosen2,1);
    trialNum_A1 = trialNum(order== 1,:);
    trialNum_A2 = trialNum(order==-1,:);
    OrHyst = zeros(size(trialNum));
    OrHyst(ismember(trialNum_A1,tr_chosen1+1)) = -1; 
    OrHyst(ismember(trialNum_A2,tr_chosen2+1)) = -1; 
    OrHyst(ismember(trialNum_A1,tr_chosen2+1)) =  1; 
    OrHyst(ismember(trialNum_A2,tr_chosen1+1)) =  1; 
	%
	% logistic regression
	% for stimOFF
	jj = ii & iOF;
	mdl = fitglm([logB_A(jj,1), OrHyst(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
	betaOFF = mdl.Coefficients.Estimate;
	rhoOFF = exp(-betaOFF(1)/betaOFF(2));
	steepnessOFF = betaOFF(2);
	if def_hyst == 0
        OrHystOFF = betaOFF(3);
    elseif def_hyst == 1
        OrHystOFF = betaOFF(3)./betaOFF(2);
    elseif def_hyst == 2
        OrHystOFF = 2*rhoOFF*betaOFF(3)./betaOFF(2);
    end  
    %
	% for stimON
	jj = ii & iON;
	mdl = fitglm([logB_A(jj,1), OrHyst(jj,1)], chosenID(jj,1) , 'Distribution','binomial' , 'link' , 'probit');
	betaON = mdl.Coefficients.Estimate;
	rhoON = exp(-betaON(1)/betaON(2));
	steepnessON = betaON(2);
	if def_hyst == 0
        OrHystON = betaON(3);
    elseif def_hyst == 1
        OrHystON = betaON(3)./betaON(2);
    elseif def_hyst == 2
        OrHystON = 2*rhoON*betaON(3)./betaON(2);
    end
    
    
    % % % 
	% save data
    sidebias_OFF(isession,:) = sidebiasOFF;
	sidebias_ON(isession,:) = sidebiasON;
	chhy_OFF(isession,:) = ChHystOFF;
	chhy_ON(isession,:) = ChHystON;
    pohy_OFF(isession,:) = PoHystOFF;
	pohy_ON(isession,:) = PoHystON;
    orhy_OFF(isession,:) = OrHystOFF;
	orhy_ON(isession,:) = OrHystON;
    
end % for isession

% % % 
% save data to info
info.sidebias_OFF = sidebias_OFF;
info.sidebias_ON = sidebias_ON;
info.chhy_OFF = chhy_OFF;
info.chhy_ON = chhy_ON;
info.pohy_OFF = pohy_OFF;
info.pohy_ON = pohy_ON;
info.orhy_OFF = orhy_OFF;
info.orhy_ON = orhy_ON;

end


% % %
% function 3
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


