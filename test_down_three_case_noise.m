function [snr_out,tfr_out]= test_down_three_case_noise(cas,sig,window,downsamp,SNR)
 
 %INPUT
 %cas      : 1, first STFT technique
 %           2, second STFT technique (filter with unit energy on its support)
 %           3, third STFT technique (filter L1 normed on its support) 
 %sig      : if sig = 1, multicomponent signal, else linear chirp.
 %window   : the type of window used hamming or Gaussian 
 %downsamp : a downsampling parameter between 1 and Lh
 %SNR      : the input signal to noise ratio
 
 %OUTPUT
 %snr_out  : the output SNR depending on the shift factor (averaged over 30
 %           realizations of the noise)
 %tfr_out  : illustration of the tfr for shift parameter equal to 0 
 
 N = 4096;
 if sig == 1
  t = (0:N-1)/N;
  a  = 2;
  s1 = a.*exp(2*pi*1i*(1000*t+60*cos(3*pi*t)));
  s2 = a.*exp(2*pi*1i*(400*t+30*cos(3*pi*t)));
  s  = s1+s2;
  s = s(:);
 else 
  s  = fmlin(N,.05,.25);
 end
 
 Nfft = 256;
 %we build the filter h

 hlength=floor(Nfft/4);
 hlength=hlength+1-rem(hlength,2);%the length of the filter has to be odd
 h = tftb_window(hlength,window); 
 
 [hrow,hcol]=size(h); 
 Lh=(hrow-1)/2;
 
 %we add some gaussian white noise to the original signal s
 nbreal = 30;
 snr_out_accu = zeros(nbreal,downsamp);
 for k = 1:nbreal, 
  n     = randn(N,1)+1i*randn(N,1);
  [s1]  = sigmerge(s,n,SNR);
  sigma = std(s1-s);
  s1    = s1(:);
 
  shift = 0;
  recons_signal = zeros(N,downsamp);  
 
  while shift < downsamp
  
   [tfr,norm2h] = tfrstft_three_case_down(s1,Nfft,cas,h,Lh,downsamp,shift);
   
   A = size(tfr);
   tfr_thresh = zeros(A);
   if shift == 0
     tfr_out = tfr;
   end    
   for icol = 1:A(2), 
    B = tfr(:,icol);
    tfr_thresh(:,icol) = B.*(abs(B)> 3*sigma*norm2h(icol));
   end
  
   recons_signal(:,shift+1) = itfrstft_three_case_down(tfr_thresh,cas,N,h,shift);
   snr_out_accu(k,shift+1)  = snr(s,s-recons_signal(:,shift+1));
   shift = shift + 1;
  end
 end
 snr_out = mean(snr_out_accu);
end


  
 