close all;
SNR = 0:5:25;
%[Ren] = optimal_size(2,1,'Gauss',SNR);
%[Ren1] = optimal_size(2,1,'hamming',SNR);
[Ren1] = optimal_size(2,3,'Gauss',SNR);