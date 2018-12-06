%% To produce the plots showing in Figure 2d, cross-temporal spatial RSA

%% add path
addpath /fieldtrip-20150923;
addpath /scripts;
addpath /scripts/functions;
ft_defaults;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spatial RSA: calculate the between-pair (randomize trials) and within-pair correlations for one subject i

%% load the paired data (trial x sensor x time) for subject i, with the corresponding trials of data_1 and data_2 predicting the same word
i=10;
load(strcat('/PairData_lpfilter30Hz_downsample_sub', sprintf('%02d', i))); %data were down sampled to 300Hz

%% define time window of interest: -1 to 0 s
timeWind = linspace(-2.0,2.0,1200); %sampling rate: 300Hz
startWind = find(timeWind-(-1)==min(abs(timeWind-(-1))));
endWind = find(timeWind-0==min(abs(timeWind-0)));

for timepoint = startWind : endWind
	Select_data_1 = data_1(:,:,timepoint);
	Select_data_2 = data_2 (:,:,[startWind:endWind]);
	Select_data_2_rand = Select_data_2(randperm(size(Select_data_2,1)),:,:); %randomize the trl info (the 1st dimension)
	
	for indtrl = 1:size(Select_data_1,1)        
		trldata_1 = squeeze(Select_data_1(indtrl,:))'; %sensor
		trldata_2 = squeeze(Select_data_2(indtrl,:,:)); %sensor x time
		trldata_2_rand = squeeze(Select_data_2_rand(indtrl,:,:));
		R(indtrl,:) = Correlation(trldata_1,trldata_2); %within-pair correlation: trial x time
		R_rand(indtrl,:) = Correlation(trldata_2,trldata_1_rand); %between-pair correlation: trial x time
	end
	win(timepoint-startWind+1,:) = mean(R,1); %average the R-values across trials, time x time matrix of R-values
	btw(timepoint-startWind+1,:) = mean(R_rand,1);	%average the R-values across trials, time x time matrix of R-values
	clear R* Select* indtrl
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% all subjects' cross-temporal similarity R-values: subject x time x time
load data/spatRSA_crosstime; %load within-pair and between-pair spatial similarity matrices

%% apply a Gaussian filter of 40 ms
data_1 = Rsubwin;
data_2 = Rsubbtw;

% Build a Gaussian filter
smooth = 40;
gaussFilter = gausswin(smooth); %40ms time window
gaussFilter = gaussFilter / sum(gaussFilter);

for j=1:size(data_1,1)
    for indRow = 1:size(data_1,2)
        data_1_extract = squeeze(data_1(j,indRow,:)); %extract single trial time series
        data_1_cal = conv(data_1_extract,gaussFilter); %convert
        conv_data_1(j,indRow,:) = data_1_cal(smooth/2:end-smooth/2); %shift

        data_2_extract = squeeze(data_2(j,indRow,:)); %extract single trial time series
        data_2_cal = conv(data_2_extract,gaussFilter); %convert
        conv_data_2(j,indRow,:) = data_2_cal(smooth/2:end-smooth/2); %shift    
    end
    clear data_1_* data_2_*
    
    for indCol = 1:size(conv_data_1,2)
        data_1_extract = squeeze(conv_data_1(j,:,indCol)); %extract single trial time series
        data_1_cal = conv(data_1_extract,gaussFilter); %convert
        conv2_data_1(j,:,indCol) = data_1_cal(smooth/2:end-smooth/2); %shift

        data_2_extract = squeeze(conv_data_2(j,:,indCol)); %extract single trial time series
        data_2_cal = conv(data_2_extract,gaussFilter); %convert
        conv2_data_2(j,:,indCol) = data_2_cal(smooth/2:end-smooth/2); %shift  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% permutation test

%% calculate the uncorrected t-values between the two conditions
data_dif = conv2_data_1 - conv2_data_2;
[h,p,ci,stats] = ttest(data_dif,zeros(size(data_dif,1),size(data_dif,2),size(data_dif,3)),'alpha',0.05,'dim',1);    
hvalues = squeeze(h);
tvalues = squeeze(stats.tstat);

%% find all the clusters and calculate the sum of the t-values
[L,n] = bwlabel(hvalues,4);
for k = 1:n
    if size(find(L==k),1) >= 1 & ~isnan(hvalues(find(L==k)))
        cluster(k) = sum(tvalues(find(L==k)));
    end
end
clear h p ci stats k

%% permutations
permN = 1000; %number of permutation
n_subs = 26; %number of subjects
alph = 0.05; %sig threshold
tail = 0; %two tails
[Perm,tmx_ptile] = twoDpermute(data_dif,permN,n_subs,alph,tail);

%% compare the observed clusters with the distribution, two tails
plotFig = L;
%pos cluster
indClusterPos = find(cluster>=tmx_ptile(2));
for m = 1:length(indClusterPos)
    highlight = find(L==indClusterPos(m));
    plotFig(highlight) = cluster(indClusterPos(m));
    pvaluePos(m) = mean(abs(Perm) > cluster(indClusterPos(m)));
end

%neg cluster
indClusterNeg = find(cluster<=tmx_ptile(1));
for m = 1:length(indClusterNeg)
    highlight = find(L==indClusterNeg(m));
    plotFig(highlight) = cluster(indClusterNeg(m));
    pvalueNeg(m) = mean(-abs(Perm) < cluster(indClusterNeg(m)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the group-averaged cross-temporal spatial similarity R-values as well as the masked difference, figure 2d
load data/stat_spatRSA_crosstime;
Rwin_avg = squeeze(mean(conv2_data_1,1)); %within-pair, smoothed R-values
Rbtw_avg = squeeze(mean(conv2_data_2,1)); %between-pair, smoothed R-values

timeWind=linspace(-1,0,300);

figure(2);
subplot(3,3,1)
clims = [-0.1 0.1];
imagesc(timeWind,fliplr(timeWind),flipud(Rwin_avg),clims)
axis xy;
colorbar
xlim(ylim)
subplot(3,3,2)
clims = [-0.1 0.1];
imagesc(timeWind,fliplr(timeWind),flipud(Rbtw_avg),clims)
axis xy;
colorbar
xlim(ylim)

%% mask the difference R-values with the statistical cluster
dif_mask = squeeze(mean(data_dif,1));
indClusterPos = find(cluster>=tmx_ptile(2));
for m = 1:length(indClusterPos)
    highlight = find(L==indClusterPos(m));
    noHighlight = setdiff([1:300*300],highlight);
    dif_mask(noHighlight) = 0;
end
subplot(3,3,3)
clims = [-0.02 0.02];
imagesc(timeWind,fliplr(timeWind),flipud(dif_mask),clims)
axis xy;
colorbar
xlim(ylim)