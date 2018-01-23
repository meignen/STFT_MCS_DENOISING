nbreal = 10;
snr_demod =zeros(nbreal,2,11);
snr_demod1 =zeros(nbreal,2,11);
snr_demod2 =zeros(nbreal,2,11);
snr_direct =zeros(nbreal,2,11);
snr_direct1 =zeros(nbreal,2,11);
snr_direct2 =zeros(nbreal,2,11);

SNR_modes_0 = zeros(nbreal,2);
SNR_modes1_0 = zeros(nbreal,2);
SNR_modes2_0 = zeros(nbreal,2);
inter_moy_0  = zeros(nbreal,2);
inter_moy1_0  = zeros(nbreal,2);
inter_moy2_0  = zeros(nbreal,2);
coeff_util_0  = zeros(nbreal,2);
coeff_util1_0  = zeros(nbreal,2);
coeff_util2_0  = zeros(nbreal,2);
SNR_modes_10 = zeros(nbreal,2);
SNR_modes1_10 = zeros(nbreal,2);
SNR_modes2_10 = zeros(nbreal,2);
inter_moy_10  = zeros(nbreal,2);
inter_moy1_10  = zeros(nbreal,2);
inter_moy2_10  = zeros(nbreal,2);
coeff_util_10  = zeros(nbreal,2);
coeff_util1_10  = zeros(nbreal,2);
coeff_util2_10  = zeros(nbreal,2);
 
N =4096;
for k=1:nbreal
 k
 %mode reconstruction using the demodulation procedure
 [snr_direct(k,:,:),snr_demod(k,:,:)]   = Recons_demod(1,0,1);% 0 dB
 [snr_direct2(k,:,:),snr_demod2(k,:,:)] = Recons_demod(1,10,1);% 10 dB
   
 %mode reconstruction based downsampled STFT using optimized Hamming filter
 %with different downsampling values

 %0 dB
 [SNR_modes_0(k,:),inter_moy_0(k,:),coeff_util_0(k,:)]    = reconstruct_modes(1,'hamming',0,64);
 [SNR_modes1_0(k,:),inter_moy1_0(k,:),coeff_util1_0(k,:)] = reconstruct_modes(1,'hamming',0,32);
 [SNR_modes2_0(k,:),inter_moy2_0(k,:),coeff_util2_0(k,:)] = reconstruct_modes(1,'hamming',0,16);
     
 %10 dB
 [SNR_modes_10(k,:),inter_moy_10(k,:),coeff_util_10(k,:)]    = reconstruct_modes(1,'hamming',10,64);
 [SNR_modes1_10(k,:),inter_moy1_10(k,:),coeff_util1_10(k,:)] = reconstruct_modes(1,'hamming',10,32);
 [SNR_modes2_10(k,:),inter_moy2_10(k,:),coeff_util2_10(k,:)] = reconstruct_modes(1,'hamming',10,16);
end

d =0:10;

close all

%SNR 0 mode 1
X1   = zeros(2,11);
X1(:,:) = mean(snr_demod);
X2 = zeros(1,2);
X3 = zeros(1,2);
X4 = zeros(1,2);
X2(:) = mean(SNR_modes_0);
X3(:) = mean(SNR_modes1_0);
X4(:) = mean(SNR_modes2_0);
XX = [X2;X3;X4];

coeff0 = zeros(1,2);
coeff1 = zeros(1,2);
coeff2 = zeros(1,2);
coeff0(:) = mean(coeff_util_0)/N;
coeff1(:) = mean(coeff_util1_0)/N;
coeff2(:) = mean(coeff_util2_0)/N;

coeff = [coeff0;coeff1;coeff2];
coeff = sum(coeff,2);
figure() 
nr =2;
plot(nr*(2*d+1),X1(1,:),coeff,XX(:,1),'--');
hold on;  
%SNR 0 mode 2
plot(nr*(2*d+1),X1(2,:),'-.',coeff,XX(:,2),':');

%SNR 10 mode 1
X1   = zeros(2,11);
X1(:,:) = mean(snr_demod2);
X2 = zeros(1,2);
X3 = zeros(1,2);
X4 = zeros(1,2);
X2(:) = mean(SNR_modes_10);
X3(:) = mean(SNR_modes1_10);
X4(:) = mean(SNR_modes2_10);
XX = [X2;X3;X4];

coeff0 = zeros(1,2);
coeff1 = zeros(1,2);
coeff2 = zeros(1,2);
coeff0(:) = mean(coeff_util_10)/N;
coeff1(:) = mean(coeff_util1_10)/N;
coeff2(:) = mean(coeff_util2_10)/N;

coeff = [coeff0;coeff1;coeff2];
coeff = sum(coeff,2);

plot(nr*(2*d+1),X1(1,:),coeff,XX(:,1),'--');
%SNR 10 mode 2
plot(nr*(2*d+1),X1(2,:),'-.',coeff,XX(:,2),':');
hold off;