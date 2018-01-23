nbreal = 5;
snr_direct =zeros(nbreal,2,11);
snr_direct1 =zeros(nbreal,2,11);
snr_direct2 =zeros(nbreal,2,11);
snr_demod =zeros(nbreal,2,11);
snr_demod1 =zeros(nbreal,2,11);
snr_demod2 =zeros(nbreal,2,11);

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

for k=1:nbreal
 k
 %mode reconstruction using the demodulation procedure
 [snr_direct(k,:,:),snr_demod(k,:,:)]   = Recons_demod(1,0,1);% 0 dB
 [snr_direct1(k,:,:),snr_demod1(k,:,:)] = Recons_demod(1,5,1);% 5 dB
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
figure()
close all
%SNR 0 mode 1
 X =zeros(2,11);
 X1 =zeros(2,11);
 X2 =zeros(1,2);
 X3 =zeros(1,2);
 X4 =zeros(1,2);
 X(:,:)= mean(snr_direct);
 X1(:,:) = mean(snr_demod);
 X2(:) = mean(SNR_modes_0);
 X3(:) = mean(SNR_modes1_0);
 X4(:) = mean(SNR_modes2_0);

 plot(d,X(1,:),d,X1(1,:),'--',...
      d,X2(1)*ones(1,length(d)),'-d',...
      d,X3(1)*ones(1,length(d)),'-<',d,X4(1)*ones(1,length(d)),'->');
 
 figure()
 %SNR 0 mode 2
 plot(d,X(2,:),d,X1(2,:),'--',...
      d,X2(2)*ones(1,length(d)),'-d',...
      d,X3(2)*ones(1,length(d)),'-<',d,X4(2)*ones(1,length(d)),'->');

 figure()
 X =zeros(2,11);
 X1 =zeros(2,11);
 X2 =zeros(1,2);
 X3 =zeros(1,2);
 X4 =zeros(1,2);

 X(:,:) = mean(snr_direct2);
 X1(:,:) = mean(snr_demod2);
 X2(:) = mean(SNR_modes_10);
 X3(:) = mean(SNR_modes1_10);
 X4(:) = mean(SNR_modes2_10);

 %SNR 10 mode 1
  plot(d,X(1,:),d,X1(1,:),'--',...
      d,X2(1)*ones(1,length(d)),'-d',...
      d,X3(1)*ones(1,length(d)),'-<',d,X3(1)*ones(1,length(d)),'->');
 
 figure()
 %SNR 10 mode 2
 plot(d,X(2,:),d,X1(2,:),'--',...
      d,X2(2)*ones(1,length(d)),'-d',...
      d,X3(2)*ones(1,length(d)),'-<',d,X4(2)*ones(1,length(d)),'->');

