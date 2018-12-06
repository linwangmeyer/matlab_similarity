%% calculate the corrections with only spatial info, Pearson correlation, within pairs
function [win] = PearsonCorrSpat_WinPair(Select_data_1,Select_data_2)
NrTrl = min([size(Select_data_1,1);size(Select_data_2,1)]);%in case the number of trials is different
for indtrl = 1:NrTrl       
    trldata_1 = squeeze(Select_data_1(indtrl,:,:));
    trldata_2 = squeeze(Select_data_2(indtrl,:,:));
    R(indtrl,:) = Correlation(trldata_1,trldata_2);
end
win = mean(R,1);
