function voraValue = voraValue(cie_cmf_pre_2, cie_cmf)
%% reference: Poorvi L. Vora and H. Joel Trussell. Measure of goodness of a set of color scanning filters. Journal of the Optical Society of America-A, vol. 10, no. 7, pp. 1499-1508, July 1993.

alpha = 3; % in this case, alpha = 3

voraValue = trace(cie_cmf_pre_2*pinv(cie_cmf_pre_2)*cie_cmf*pinv(cie_cmf)) / alpha;

end