 SNR = [0 10];
 nbreal = 5;
 %Gaussian window
 [SNR_modesG_0,downsamp]  = test_reconstruct_three_case_modes_noshift(2,1,1,'Gauss',SNR(1),5);
 [SNR_modesG1_0,downsamp] = test_reconstruct_three_case_modes_noshift(2,2,1,'Gauss',SNR(1),5);

 %Hamming Window 
 [SNR_modesH_0,downsamp]  = test_reconstruct_three_case_modes_noshift(2,1,1,'Hamming',SNR(1),5);

 %Gaussian window
 [SNR_modesG_10,downsamp]  = test_reconstruct_three_case_modes_noshift(2,1,1,'Gauss',SNR(2),5);
 [SNR_modesG1_10,downsamp] = test_reconstruct_three_case_modes_noshift(2,2,1,'Gauss',SNR(2),5);

 %Hamming Window 
 [SNR_modesH_10,downsamp]  = test_reconstruct_three_case_modes_noshift(2,1,1,'Hamming',SNR(2),5);

 figure()

 % first mode
 plot(downsamp,SNR_modesG_0(:,1),downsamp,SNR_modesG1_0(:,1),'--',...
      downsamp,SNR_modesH_0(:,1),'-.',downsamp,SNR_modesG_10(:,1),'-s',...
      downsamp,SNR_modesG1_10(:,1),'-d',downsamp,SNR_modesH_10(:,1),'-c')
 
 figure 
 
 plot(downsamp,SNR_modesG_0(:,2),downsamp,SNR_modesG1_0(:,2),'--',...
      downsamp,SNR_modesH_0(:,2),'-.',downsamp,SNR_modesG_10(:,2),'-s',...
      downsamp,SNR_modesG1_10(:,2),'-d',downsamp,SNR_modesH_10(:,2),'-c');
