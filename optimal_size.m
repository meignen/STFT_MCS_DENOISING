function [Ren] = optimal_size(cas,sig,window,SNR)
 %INPUT
 %cas        : 1, first STFT technique
 %             2, second STFT technique (filter with unit energy on its support)
 %             3, third STFT technique (filter L1 normed on its support)  
 %sig        : type of studied signal
 %window     : choice for COLA(R) window
 %SNR        : SNR corresponding to the added noise
 
 %OUTPUT
 %Ren       : optimal window size as a function depending on the SNR 
 
 N = 4096;
 if sig == 1
  t  = (0:N-1)/N;
  a  = 2;
  s1 = a.*exp(2*pi*1i*(400*t+30*cos(3*pi*t)));
  s2 = a.*exp(2*pi*1i*(1000*t+60*cos(3*pi*t)));
  s  = s1+s2;
  s = s(:);
 elseif sig == 2 
  t  = (0:N-1)/N;   
  s  = fmlin(N,.05,0.25);
 else
  load -ascii batsig.txt
  s = batsig;
  s = s(150:end)';
  s = hilbert(s);
  s = s(:);
  N = length(s);
  t  = (0:N-1)/N;
  if (N == 2^(floor(log2(N))))
   Nfft = N;
  else
   Nfft = 2^(floor(log2(N))+1);
  end
 end
 
 if (sig <= 2)
  Nfft = 512;
 end 
 
 %we build the filter h
 
 if strcmp(window,'hamming')
  taille = 1:0.5:10;
  Ren = zeros(length(SNR),length(taille));
  for k = 1:length(SNR),
   n    = randn(N,1)+1i*randn(N,1);
   [sn] = sigmerge(s,n,SNR(k));
   sn   = sn(:);
   for p = 1:length(taille),
    
    %computation of the filter
    hlength=floor(Nfft/taille(p));
    hlength=hlength+1-rem(hlength,2);%the length of the filter has to be odd
    h = tftb_window(hlength,window);
    %h(end) = 0;
    [hrow,hcol] = size(h); 
    Lh = (hrow-1)/2;
    
    [tfr,norm2h] = tfrstft_three_case_down(sn,Nfft,cas,h,Lh,1,0);
    %computation of the Rényi entropy
    fs = 0:Nfft-1;
    Ren(k,p) = renyi_entropy(abs(tfr),t,fs',3);  
   end 
  end
  figure()
  hold on;
  hlength=floor(Nfft./taille);
  hlength=hlength+1-rem(hlength,2);
  
  for k = 1:length(SNR)
   plot(hlength,Ren(k,:))   
  end
  hold off; 
 else
   %the window is the Gaussian window    
   sigma_w = 0.1:0.5:8;
   Ren = zeros(length(SNR),length(sigma_w));
   for k = 1:length(SNR),
    n    = randn(N,1)+1i*randn(N,1);
    [sn] = sigmerge(s,n,SNR(k));
    sn   = sn(:);
    for p = 1:length(sigma_w),
     
     %computation of the filter    
     prec = 10^(-3);
     L =  sqrt(Nfft)*sigma_w(p);
     Lh = floor(sigma_w(p)*sqrt(-Nfft*log(prec)/pi))+1;
     h = amgauss(2*Lh+1,Lh+1,L);
     [tfr,norm2h] = tfrstft_three_case_down(sn,Nfft,cas,h,Lh,1,0);
     %computation of the Rényi entropy
     fs = 0:Nfft-1;
     Ren(k,p) = renyi_entropy(abs(tfr),t,fs',3);
    end 
   end
   figure()
   hold on;
   for k = 1:length(SNR)
    plot(sigma_w/sqrt(Nfft),Ren(k,:))   
    [A,I] = min(Ren(k,:));
    SNR(k)
    sigma_w(I(1))/sqrt(Nfft)
   end
   hold off; 
  end
end  