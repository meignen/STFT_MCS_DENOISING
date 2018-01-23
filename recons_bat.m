function [X,XX,coeffX,Y,YY,coeffY,Z,ZZ,coeffZ]=recons_bat(nr)
% Tests the computation of basin of attraction on the real bat signal
nbreal = 20;

snr_demod =zeros(nbreal,11);
snr_demod1 =zeros(nbreal,11);
snr_demod2 =zeros(nbreal,11);
snr_direct =zeros(nbreal,11);
snr_direct1 =zeros(nbreal,11);
snr_direct2 =zeros(nbreal,11);

SNR_modes_0 = zeros(nbreal,1);
SNR_modes1_0 = zeros(nbreal,1);
SNR_modes2_0 = zeros(nbreal,1);
inter_moy_0  = zeros(nbreal,nr);
inter_moy1_0  = zeros(nbreal,nr);
inter_moy2_0  = zeros(nbreal,nr);
coeff_util_0  = zeros(nbreal,nr);
coeff_util1_0  = zeros(nbreal,nr);
coeff_util2_0  = zeros(nbreal,nr);
SNR_modes_5 = zeros(nbreal,1);
SNR_modes1_5 = zeros(nbreal,1);
SNR_modes2_5 = zeros(nbreal,1);
inter_moy_5  = zeros(nbreal,nr);
inter_moy1_5  = zeros(nbreal,nr);
inter_moy2_5  = zeros(nbreal,nr);
coeff_util_5  = zeros(nbreal,nr);
coeff_util1_5  = zeros(nbreal,nr);
coeff_util2_5  = zeros(nbreal,nr);
SNR_modes_10 = zeros(nbreal,1);
SNR_modes1_10 = zeros(nbreal,1);
SNR_modes2_10 = zeros(nbreal,1);
inter_moy_10  = zeros(nbreal,nr);
inter_moy1_10  = zeros(nbreal,nr);
inter_moy2_10  = zeros(nbreal,nr);
coeff_util_10  = zeros(nbreal,nr);
coeff_util1_10  = zeros(nbreal,nr);
coeff_util2_10  = zeros(nbreal,nr);
load -ascii batsig.txt
s = batsig;
s = s(145:end)';
N = length(s); 

for k=1:nbreal
 k
 %mode reconstruction using the demodulation procedure
 
 [snr_direct(k,:),snr_demod(k,:)]   = Recons_demod_nr(3,0,1,nr);% 0 dB
 [snr_direct1(k,:),snr_demod1(k,:)] = Recons_demod_nr(3,5,1,nr);% 5 dB
 [snr_direct2(k,:),snr_demod2(k,:)] = Recons_demod_nr(3,10,1,nr);% 10 dB
 
 %mode reconstruction based downsampled STFT using optimized Hamming filter
 %with different downsampling values
 
  %0 dB
  [SNR_modes_0(k),inter_moy_0(k,:),coeff_util_0(k,:)]    = reconstruct_modes_nr(3,'hamming',0,32,nr);
  [SNR_modes1_0(k),inter_moy1_0(k,:),coeff_util1_0(k,:)] = reconstruct_modes_nr(3,'hamming',0,16,nr);
  [SNR_modes2_0(k),inter_moy2_0(k,:),coeff_util2_0(k,:)] = reconstruct_modes_nr(3,'hamming',0,8,nr);
    
  %5 dB
  [SNR_modes_5(k),inter_moy_5(k,:),coeff_util_5(k,:)]    = reconstruct_modes_nr(3,'hamming',5,32,nr);
  [SNR_modes1_5(k),inter_moy1_5(k,:),coeff_util1_5(k,:)] = reconstruct_modes_nr(3,'hamming',5,16,nr);
  [SNR_modes2_5(k),inter_moy2_5(k,:),coeff_util2_5(k,:)] = reconstruct_modes_nr(3,'hamming',5,8,nr);
 
  %10 dB
  [SNR_modes_10(k),inter_moy_10(k,:),coeff_util_10(k,:)]    = reconstruct_modes_nr(3,'hamming',10,32,nr);
  [SNR_modes1_10(k),inter_moy1_10(k,:),coeff_util1_10(k,:)] = reconstruct_modes_nr(3,'hamming',10,16,nr);
  [SNR_modes2_10(k),inter_moy2_10(k,:),coeff_util2_10(k,:)] = reconstruct_modes_nr(3,'hamming',10,8,nr);
end

%SNR 0 
X = mean(snr_demod);

X2 = mean(SNR_modes_0);
X3 = mean(SNR_modes1_0);
X4 = mean(SNR_modes2_0);
XX = [X2;X3;X4];

coeff0 = zeros(1,nr);
coeff1 = zeros(1,nr);
coeff2 = zeros(1,nr);
coeff0(:) = mean(coeff_util_0)/N;
coeff1(:) = mean(coeff_util1_0)/N;
coeff2(:) = mean(coeff_util2_0)/N;

coeff = [coeff0;coeff1;coeff2];
coeffX = sum(coeff,2);

%SNR 5
Y = mean(snr_demod1);
X2 = mean(SNR_modes_5);
X3 = mean(SNR_modes1_5);
X4 = mean(SNR_modes2_5);
YY = [X2;X3;X4];

coeff0 = mean(coeff_util_5)/N;
coeff1 = mean(coeff_util1_5)/N;
coeff2 = mean(coeff_util2_5)/N;

coeff = [coeff0;coeff1;coeff2];
coeffY = sum(coeff,2);

%SNR 10 
Z = mean(snr_demod2);
X2 = mean(SNR_modes_10);
X3 = mean(SNR_modes1_10);
X4 = mean(SNR_modes2_10);
ZZ = [X2;X3;X4];

coeff0 = mean(coeff_util_10)/N;
coeff1 = mean(coeff_util1_10)/N;
coeff2 = mean(coeff_util2_10)/N;

coeff = [coeff0;coeff1;coeff2];
coeffZ = sum(coeff,2);
end
