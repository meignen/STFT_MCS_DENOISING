% we compute the ouput SNR as function of the shift parameter
% for the three different techniques and for the two different type of modes
% the standard deviation of the noise is assumed to be known

[SNR_out1_10,tfr_out1] = test_down_three_case_noise(1,0,'hamming',32,10);
[SNR_out2_10,tfr_out2] = test_down_three_case_noise(2,0,'hamming',32,10);
[SNR_out3_10,tfr_out3] = test_down_three_case_noise(3,0,'hamming',32,10);

[SNR_out01_10,tfr_out01]  = test_down_three_case_noise(1,1,'hamming',32,10);
[SNR_out02_10,tfr_out02] = test_down_three_case_noise(2,1,'hamming',32,10);
[SNR_out03_10,tfr_out03] = test_down_three_case_noise(3,1,'hamming',32,10);
close all
figure()
B = size(tfr_out1);
imagesc(1:B(2),(0:B(1)/2-1)/B(1),abs(tfr_out1(1:B(1)/2,:)));
set(gca,'ydir','normal');
axis square

figure()
B = size(tfr_out01);
imagesc(1:B(2),(0:B(1)/2-1)/B(1),abs(tfr_out01(1:B(1)/2,:)));
set(gca,'ydir','normal');
axis square

[SNR_out1_5,tfr_out1]  = test_down_three_case_noise(1,0,'hamming',32,5);
[SNR_out2_5,tfr_out2] = test_down_three_case_noise(2,0,'hamming',32,5);
[SNR_out3_5,tfr_out3] = test_down_three_case_noise(3,0,'hamming',32,5);

[SNR_out01_5,tfr_out01]  = test_down_three_case_noise(1,1,'hamming',32,5);
[SNR_out02_5,tfr_out02] = test_down_three_case_noise(2,1,'hamming',32,5);
[SNR_out03_5,tfr_out03] = test_down_three_case_noise(3,1,'hamming',32,5);


figure()
%plot of the linear chirp results
plot(0:31,SNR_out1_5,'-',0:31,SNR_out2_5,'--',0:31,SNR_out3_5,'-.',... 
     0:31,SNR_out1_10,'-d',0:31,SNR_out2_10,'-*',0:31,SNR_out3_10,'-s','LineWidth',2,'MarkerSize',10)
 
figure()
%plot of the cosine phase modes
plot(0:31,SNR_out01_5,'-',0:31,SNR_out02_5,'--',0:31,SNR_out03_5,'-.',... 
     0:31,SNR_out01_10,'-d',0:31,SNR_out02_10,'-*',0:31,SNR_out03_10,'-s','LineWidth',2,'MarkerSize',10)
