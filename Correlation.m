%% calculate the correlation between multidimensional data a and b

function r = Correlation(a,b)
az = bsxfun(@minus, a, mean(a,1));%subtract the mean over the first dimension
bz = bsxfun(@minus, b, mean(b,1));
% Standard Pearson correlation coefficient formula
a2 = az .^ 2; %self variance
b2 = bz .^ 2; %self variance
ab = bsxfun(@times,az,bz); %covariance
r = sum(ab, 1) ./ sqrt(bsxfun(@times,sum(a2, 1),sum(b2, 1)));