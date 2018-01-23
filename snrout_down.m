function [snr_out,nbcoeff]= snrout_down(sig,window,downsamp,SNR)
 
 %INPUT
 %sig      : if sig = 1, multicomponent signal, else linear chirp.
 %window   : the type of window used hamming or Gaussian 
 %downsamp : a downsampling parameter between 1 and Lh
 %SNR      : the input signal to noise ratio
 
 %OUTPUT
 %snr_out  : the output SNR (averaged over 20 realizations of the noise)
 %nbcoeff  : number of coefficients kept 
 
 N = 4096;
 if sig == 1
  t = (0:N-1)/N;
  a  = 2;
  s1 = a.*exp(2*pi*1i*(1000*t+60*cos(3*pi*t)));
  s2 = a.*exp(2*pi*1i*(400*t+30*cos(3*pi*t)));
  s  = s1+s2;
  s = s(:);
  Lg = 161;
  sigma_opt = 0.15;
 elseif sig == 2 
  s  = fmlin(N,.05,.25);
  Lg = 161;
  sigma_opt = 0.15;
 else
  load -ascii batsig.txt
  s = batsig;
  s = s(145:end)';
  s = hilbert(s);
  s = s(:);
  N = length(s);
  if (N == 2^(floor(log2(N))))
   Nfft = N;
  else
   Nfft = 2^(floor(log2(N))+1);
  end
  sigma_opt = 0.13;   
  Lg = 65;
 end
 
 if (sig <= 2)
  Nfft = 256;
 end
 %we build the filter h
 s    = s(:);
 
 %we build the filter h
 if strcmp(window,'hamming')
  hlength=floor(Lg);%optimal window determined by Rényi entropy
  hlength=hlength+1-rem(hlength,2);%the length of the filter has to be odd
  h = tftb_window(hlength,window);
  [hrow,hcol]=size(h); 
  Lh=(hrow-1)/2;
 else
  %the window is the Gaussian window    
  prec = 10^(-3);
  L =  sigma_opt*Nfft;
  Lh = floor(L*sqrt(-log(prec)/pi))+1;
  h = amgauss(2*Lh+1,Lh+1,L); 
 end
 
 %we add some gaussian white noise to the original signal s
 nbreal = 20;
 snr_out_accu = zeros(nbreal,1);
 nbcoeff_accu = zeros(nbreal,1);

 for k = 1:nbreal, 
  n     = randn(N,1)+1i*randn(N,1);
  [sn]  = sigmerge(s,n,SNR);
  sn    = sn(:);
 
  [tfr,norm2h] = tfrstft_three_case_down(sn,Nfft,2,h,Lh,downsamp,0);
 
  %estimation of the noise level
  Y2 = real(tfr);
  gamma_estime = median(abs(Y2(:)))/0.6745;      
  
  A = size(tfr);
  tfr_thresh = zeros(A);
      
  for icol = 1:A(2), 
   B = tfr(:,icol);
   tfr_thresh(:,icol) = B.*(abs(B)> 3*gamma_estime);
  end
  
  X = abs(tfr) > 3*gamma_estime;
  
   nbcoeff_accu(k) = sum(sum(X))/N;
   recons_signal = itfrstft_three_case_down(tfr_thresh,2,N,h,0);
   snr_out_accu(k) = snr(s,s-recons_signal);
  end
 snr_out  = mean(snr_out_accu);
 nbcoeff  = mean(nbcoeff_accu);
end