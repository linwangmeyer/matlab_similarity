%% To produce the plots showing in Figures 2a, 2b and 2c

%% add path
addpath fieldtrip-20150923;
addpath scripts;
addpath scripts/functions;
ft_defaults;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spatial RSA: calculate the between-pair and within-pair correlations for one subject i
%% load the paired data (trial x sensor x time) for subject i, with the corresponding trials of data_1 and data_2 predicting the same word   
i=10;
load(strcat('/PairData_lpfilter30Hz_sub', sprintf('%02d', i)));
Rsubbtw_1 = PearsonCorrSpat_shift_BtwPair(data_1, data_1);
Rsubbtw_2 = PearsonCorrSpat_shift_BtwPair(data_2, data_2);
Rsubbtw_3 = PearsonCorrSpat_shift_BtwPair(data_1, data_2);
Rsubbtw(i,:) = (Rsubbtw_1 + Rsubbtw_2 + Rsubbtw_3)./3;
Rsubwin(i,:) = PearsonCorrSpat_WinPair(data_1, data_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistical analysis to the spatial similarity R-values averaged within the defined time window
load data/spatRSA; %load within-pair and between-pair spatial similarity time series

%% to find out the time window that shows larger R values than the R=0.04 in the -1 - 0s time window
timeWind = linspace(-2.0,2.0,4800);
StartTime = find(abs(timeWind-(-1))==min(abs(timeWind-(-1))));
EndTime = find(abs(timeWind-0)==min(abs(timeWind-0)));

comb_data = (Rsubbtw+Rsubwin)./2; %combine the within-pair and between-pair time series
comb_data_avg = mean(comb_data,1); % group-averaged combined time series

time_ROI_ind = find(comb_data_avg(1,[StartTime:EndTime])>=0.04);
time_ROI_start = timeWind(StartTime+time_ROI_ind(1));
time_ROI_end = timeWind(StartTime+time_ROI_ind(end));

%% pairwise comparison between the within-pair and between-pair R values within the time window
[min_val StartSample] = min(abs(timeWind-time_ROI_start));
[max_val EndSample] = min(abs(timeWind-time_ROI_end));
within = mean(Rsubwin(:,[StartSample:EndSample]),2);
between = mean(Rsubbtw(:,[StartSample:EndSample]),2);
[h p ci stats] = ttest(within,between)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot time course of the spatial similarity time series for Figures 2a, 2b and 2c
Rwin_avg = mean(Rsubwin,1); % group-averaged within-pair time series
Rbtw_avg = mean(Rsubbtw,1); % group-averaged between-pair series

timeWind = linspace(-2.0,2.0,4800); %% define the whole time window
figure(2);
subplot(3,3,1)
plot(timeWind, comb_data_avg,'k','LineWidth',2);
axis([-2 1.0 0 0.12]);
hold on;
plot(timeWind, ones(1,length(timeWind))*0.04,'k','LineWidth',2);
title ('spatial RSA time course for combined conditions')

subplot(3,3,2)
plot(timeWind, Rwin_avg,'r','LineWidth',2);
axis([-2 1.0 0 0.12]);
hold on;
plot(timeWind, Rbtw_avg,'b','LineWidth',2);
axis([-2 1.0 0 0.12]);
hold on;
plot(timeWind, ones(1,length(timeWind))*0.04,'k','LineWidth',2);
title ('spatial RSA time course for separate conditions')

%% plot the scatter plot
subplot(3,3,3)
scatter(within,between,'filled','d')
axis([0 0.14 0 0.14]);
line(xlim,ylim,'Color','k')
title ('scatter plot of averaged R values')
