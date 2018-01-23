function [SNR_modes,inter_moy,coeff_util] = reconstruct_modes(sig,window,SNR,downsamp)
 %INPUT
 %sig        : type of studied signal
 %window     : window used 
 %SNR        : SNR corresponding to the added noise
 %downsamp   : downsmpling factor
 
 %OUTPUT
 %SNR_modes  : SNR associated with signal reconstruction
 %inter_moy  : average number of coefficients used at considered time instants
 %coeff_util : total number of coefficients used for mode reconstruction
 
 N = 4096;
 if sig == 1
  t  = (0:N-1)/N;
  a  = 2;
  s1 = a.*exp(2*pi*1i*(400*t+30*cos(3*pi*t)));
  s2 = a.*exp(2*pi*1i*(1000*t+60*cos(3*pi*t)));
  s  = s1+s2;
  s = s(:);
  Sacc = [s1(:) s2(:)];
  nr = 2;
  Lg = 161;
  sigma_opt = 0.15;
  clwin = 10;
 elseif sig == 2 
  s  = fmlin(N,.05,0.25);
  Sacc = s;
  nr = 1;
  Lg = 161;
  sigma_opt = 0.15;
  clwin = 10;
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
  clwin = 10;
  nr = 3;
 end
 
 if (sig <= 2)
  Nfft = 512;
 end
 
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
 
 n = randn(N,1)+1i*randn(N,1);
 [sn] = sigmerge(s,n,SNR);
 sn   = sn(:);
  
 [tfr,norm2h] = tfrstft_three_case_down(sn,Nfft,2,h,Lh,downsamp,0); 
 Abstfr = abs(tfr);
 %estimation of the noise level
 Y2 = real(tfr);
 gamma_estime = median(abs(Y2(:)))/0.6745;
     
 %ridge extraction
 [Cs] = exridge_mult(tfr,nr,0,0,clwin);
 Cs = Cs';
 modes = zeros(N,nr);
 B = size(tfr);  
 index = 1:N;

 if (sig <= 2) 
  SNR_modes = zeros(1,nr);
 end
 
 interval  = zeros(B(2),nr);
 for j=1:nr,
  %construction of the TF mask
  tfr_int = zeros(B);  
  for r = 1:B(2),  
   
   val = 3*gamma_estime; %trheshold for the transfrom depending on the noise level
  
   if (Abstfr(Cs(r,j),r) > val)
    k1 = 0;
    k2 = 0;
     
    eta1 = - 1;
    while (eta1 < 0)&&(Abstfr(Cs(r,j)-min(k1,Cs(r,j)-1),r) > val)
     if (k1 ~= Cs(r,j)-1)
      k1 = k1+1;
     else
      eta1 = k1;   
     end
    end
    if (eta1 < 0)
     eta1 = k1-1;
    end
     
    eta2 = -1;
    while (eta2 < 0) && (Abstfr(Cs(r,j)+min(k2,B(1)-Cs(r,j)),r) > val)
     if (k2 ~= B(1)-Cs(r,j))
      k2 = k2+1;   
     else
      eta2 = k2;   
     end
    end
    if (eta2 < 0)
     eta2 = k2;
    end
    interval(r,j) =eta2+eta1+1;
    tfr_int(max(1,Cs(r,j)-eta1):min(B(1),Cs(r,j)+eta2),r) = tfr(max(1,Cs(r,j)-eta1):min(B(1),Cs(r,j)+eta2),r);
   end
  end
  modes(:,j) = itfrstft_three_case_down(tfr_int,2,N,h,0);
  if (sig <= 2)
   SNR_modes(j) = snr(Sacc(index,j),modes(index,j)-Sacc(index,j));
  end
 end
 
 if (sig > 2)
  %size(modes)
  modes = transpose(modes);
  sigrec = sum(modes);
  sigrec = transpose(sigrec);
  SNR_modes = snr(s(index),sigrec(index)-s(index));
 end
 inter_moy   = mean(interval);
 coeff_util  = sum(interval);  
  
 