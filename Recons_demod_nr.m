function [snr_direct,snr_demod] = Recons_demod_nr(sig,SNR,fact,nr)
 
 N = 4096;
 if (sig == 1)
  t  = (0:N-1)/N;
  a  = 2;
  s1 = a.*exp(2*pi*1i*(400*t+30*cos(3*pi*t)));
  s2 = a.*exp(2*pi*1i*(1000*t+60*cos(3*pi*t)));
  %phi_theo = [ 400-90*pi*sin(3*pi*t); 1000-180*pi*sin(3*pi*t)]; 
  s  = s1+s2;
  s = s(:);
  signal = [s1; s2];
  nr = 2;
  clwin = 10;
  sigma_opt = 0.15;
 elseif (sig == 2)
  t  = (0:N-1)/N;   
  s  = fmlin(N,.05,0.25);
  signal = transpose(s);
  nr = 1;
  clwin = 10;
  sigma_opt = 0.15;
 else
  load -ascii batsig.txt
  s = batsig;
  s = s(145:end)';
  s = hilbert(s);
  s = s(:);
  N = length(s);
  t  = (0:N-1)/N;
  if (N == 2^(floor(log2(N))))
   Nfilt = N;
  else
   Nfilt = 2^(floor(log2(N))+1);
  end
  sigma_opt = 0.13;   
  Lg = 65;
  clwin = 10;
  %nr = 3;    
 end 

 if (sig <= 2)
  Nfilt  = 512;
 end
 
 %the window is the Gaussian window    
 
 prec  = 10^(-3);
 L     =  sigma_opt*Nfilt;
 Lh    = floor(L*sqrt(-log(prec)/pi))+1;
 h     = amgauss(2*Lh+1,Lh+1,L); 
 n     = randn(N,1)+1i*randn(N,1);
 [sn]  = sigmerge(s,n,SNR);
 
 %estimation of the threshold gamma
 %we use zero-padding to improve frequency resolution
 
 Nfft4 = fact*Nfilt;
 [tfr,norm2h] = tfrstft_three_case_down(sn,Nfft4,1,h,Lh,1,0); 
 Y2 = real(tfr);
 gamma = median(abs(Y2(:)))/0.6745; 
 
 %computation of the STFT, SST and VSST, and of the ridges associated with
 %VSST
 [STFT,SST,VSST] = sst2_new(sn,sigma_opt,Nfilt,Nfft4,3*gamma);
 
 [Cs]  = exridge_mult(VSST,nr,0,0,clwin*Nfft4/Nfilt);
%  imagesc(abs(VSST));
%  set(gca,'ydir','normal');
%  hold on;
%  B = size(VSST);
%  plot(1:B(2),Cs(1,:)-1,1:B(2),Cs(2,:)-1);
%  hold off; 
%  pause
 %computation of the phase of the modes and of the demodulated signals
 
 [sp1_s,integ1] = demod_multi(sn,VSST,Nfft4,t,N,nr,0,0,clwin*Nfft4/Nfilt);
 index = 1:N; 
 
 %signal reconstruction
 d = 0:1:10;
 if (sig <= 2)
  snr_direct = zeros(nr,length(d));
  snr_demod  = zeros(nr,length(d)); %using demodulated signal
 else
  snr_direct = zeros(1,length(d));
  snr_demod  = zeros(1,length(d)); %using demodulated signal
 end

 if (sig <= 2)
  sign1  = zeros(nr,length(s));
 else
  sign1  = zeros(nr,length(d),length(s));
  imf    = zeros(nr,length(d),length(s));
 end
 
 for p = 1:nr,
  [STFT1,SST1,VSST_sd_1] = sst2_new(sp1_s(p,:),sigma_opt,Nfilt,Nfilt,3*gamma);
  VSST_sd_int  = zeros(size(VSST_sd_1));
  VSST_sd_int(round(100*(Nfilt/N))-clwin:round(100*(Nfilt/N))+clwin,:) =... 
        VSST_sd_1(round(100*(Nfilt/N))-clwin:round(100*(Nfilt/N))+clwin,:);
  [C, Es] = exridge(VSST_sd_int,0,0,clwin);
  
  if (sig <= 2)
   for  k = 1:length(d),
      
    sign1(p,:) = recmodes(VSST,Cs(p,:),d(k));   
    snr_direct(p,k) =  snr(signal(p,index),sign1(p,index)-signal(p,index));
   
    imf        = recmodes(VSST_sd_1,C,d(k));
    imf        = imf.*exp(2*1i*pi*(integ1(p,:)-100.*t));
    snr_demod(p,k)  = snr(signal(p,index),imf(index)-signal(p,index));
   end
  else
   for  k = 1:length(d),
    sign1(p,k,:) = recmodes(VSST,Cs(p,:),d(k));  
    imf1   = recmodes(VSST_sd_1,C,d(k));
    imf(p,k,:)   = imf1.*exp(2*1i*pi*(integ1(p,:)-100.*t));
   end
  end
 end 
 if (sig == 3)
  sign_direct = zeros(length(d),length(s));
  sign_direct(:,:) = sum(sign1);
  sign_demod  = zeros(length(d),length(s));
  sign_demod(:,:) = sum(sign1);
  
  X = zeros(length(s),1);
  for k=1:length(d)
   X(:) = sign_direct(k,index);
   snr_direct(k) =  snr(s(index),X-s(index));
    X(:) = sign_demod(k,index);
   snr_demod(k)  = snr(s(index),X-s(index));  
  end
 end 