%% calculate the corrections with only spatial info, Pearson correlation, use the min number of trls
function [btw_mean, btw_SE] = PearsonCorrSpat_shift_BtwPair(Select_data_1,Select_data_2)

NrTrl = min([size(Select_data_1,1);size(Select_data_2,1)]);%in case the number of trials is different
reshape_data_2 = reshape(Select_data_2,size(Select_data_2,1),size(Select_data_2,2)*size(Select_data_2,3));

for indshift = 1:NrTrl-1
    shift_data_2 = circshift(reshape_data_2,indshift);
    data_2_new = reshape(shift_data_2,size(Select_data_2,1),size(Select_data_2,2),size(Select_data_2,3));
    
    for indtrl = 1:NrTrl     
        trldata_1 = squeeze(Select_data_1(indtrl,:,:));
        trldata_2_shift = squeeze(data_2_new(indtrl,:,:));
        R_shift(indtrl,:) = Correlation(trldata_1,trldata_2_shift);
    end
    Rtrl_shift_mean(indshift,:) = mean(R_shift,1);
    Rtrl_shift_SE(indshift,:) = std(R_shift,1)/sqrt(NrTrl);
end
btw_mean = mean(Rtrl_shift_mean,1);
btw_SE = mean(Rtrl_shift_SE,1);



