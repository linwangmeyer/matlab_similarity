%% calculate the corrections with only temporal info, Pearson correlation, within pairs
function [win] = PearsonCorrTemp_WinPair(Select_data_1,Select_data_2)
for indtrl = 1:size(Select_data_1,1)        
    trldata_1 = squeeze(Select_data_1(indtrl,:,:));
    trldata_2 = squeeze(Select_data_2(indtrl,:,:));
    R(indtrl,:) = Correlation(trldata_1',trldata_2');
end
win = mean(R,1);