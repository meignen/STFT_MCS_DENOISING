function [SNR_modes,downsamp] = test_reconstruct_three_case_modes_noshift(cas,method,sig,window,SNR,nbreal)
 %INPUT
 %cas        : 1, first STFT technique
 %             2, second STFT technique (filter with unit energy on its support)
 %             3, third STFT technique (filter L1 normed on its support)  
 %sig        : type of studied signal
 %window     : Gaussian or Hamming window
 %downsamp   : downsampling factor between 1 and Lh
 %method     : if method = 1, the proposed method of the paper 
 %             if method = 2, when the window is Gaussian, we use first
 %             order approximation to recontruct
 %SNR        : SNR corresponding to the added noise
 %nbreal     : number of realizations
 
 %OUTPUT
 %SNR_modes  : SNR associated with signal reconstruction as a function of
 %             the shift parameter 
 %tfr_out    : the time frequency representation associated with
 %             downsampling factor downsamp  
 
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
  clwin = 10;
 else 
  s  = fmlin(N,.05,0.25);
  Sacc = s;
  nr = 1;
  clwin = 10;
 end

 Nfft = 512;
 s    = s(:);
 
 %we build the filter h
 if strcmp(window,'hamming')
  hlength=floor(161);%optimal window determined by Rényi entropy
  hlength=hlength+1-rem(hlength,2);%the length of the filter has to be odd
  h = tftb_window(hlength,window);
  [hrow,hcol]=size(h); 
  Lh=(hrow-1)/2;
 else
  %the window is the Gaussian window    
  prec = 10^(-3);
  L =  0.15*Nfft;
  Lh = floor(L*sqrt(-log(prec)/pi))+1;
  h = amgauss(2*Lh+1,Lh+1,L); 
 end
 
 downsamp  = [1 4 8 16 32 64];
 SNR_modes = zeros(length(downsamp),nr);
 
 for nb = 1:nbreal
  nb
  n = randn(N,1)+1i*randn(N,1);
  [sn] = sigmerge(s,n,SNR);
  sn   = sn(:);
  
  if (method == 1)
   for p = 1:length(downsamp),
  
    SNR_int = zeros(1,nr);
    [tfr,norm2h] = tfrstft_three_case_down(sn,Nfft,cas,h,Lh,downsamp(p),0); 
   
     %estimation of the noise level
     Y2 = real(tfr);
     gamma_estime = median(abs(Y2(:)))/0.6745;
     
     %ridge extraction
     [Cs] = exridge_mult(tfr,nr,0,0,clwin);
     Cs = Cs';
     Abstfr = abs(tfr);   
   
     modes = zeros(N,nr);
     B = size(tfr);
   
     for j=1:nr,
   
      %construction of the TF mask
      tfr_int = zeros(B);
    
      for r = 1:B(2),  
       val = 3*gamma_estime; %trheshold for the transfrom depending on the noise level
       k1 = 1;
       k2 = 1;
     
       eta1 = - 1;
       while (eta1 < 0)&&(Abstfr(Cs(r,j)-min(k1,Cs(r,j)-1),r) > val)
        if (k1 ~= Cs(r,j)-1)
         k1 = k1+1;
        else
         eta1 = k1;   
        end
       end
       if (eta1 < 0)
        eta1 = k1 - 1;
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
        eta2 = k2-1;
       end    
       tfr_int(max(1,Cs(r,j)-eta1):min(B(1),Cs(r,j)+eta2),r) = tfr(max(1,Cs(r,j)-eta1):min(B(1),Cs(r,j)+eta2),r);
      end
    
      modes(:,j) = itfrstft_three_case_down(tfr_int,cas,N,h,0);
      SNR_int(j) = snr(Sacc(:,j),modes(:,j)-Sacc(:,j));
     end
     
    if nb == 1
     SNR_modes(p,:) = SNR_int;
    else
     SNR_modes(p,:) = SNR_modes(p,:) + SNR_int;
    end
   end
  elseif (method == 2) 
   for p = 1:length(downsamp),   
    
    SNR_int = zeros(1,nr);
    [tfr,norm2h] = tfrstft_three_case_down(sn,Nfft,cas,h,Lh,downsamp(p),0); 
   
    %estimation of the noise level
    Y2 = real(tfr);
    gamma_estime = median(abs(Y2(:)))/0.6745;
     
    %ridge extraction
    [Cs] = exridge_mult(tfr,nr,0,0,clwin);
    Cs = Cs';
    Abstfr = abs(tfr);   
   
    modes = zeros(N,nr);
    B = size(tfr);
   
    for j=1:nr,
   
     %construction of the TF mask
     tfr_int = zeros(B);
    
     for r = 1:B(2),  
      val = 3*gamma_estime; %trheshold for the transfrom depending on the noise level
      eta = round(1/0.15*sqrt(-1/pi*log(val/(Abstfr(Cs(r,j),r)))));        
      tfr_int(max(1,Cs(r,j)-eta):min(B(1),Cs(r,j)+eta),r) = tfr(max(1,Cs(r,j)-eta):min(B(1),Cs(r,j)+eta),r);
     end
     
     modes(:,j) = itfrstft_three_case_down(tfr_int,cas,N,h,0);
     SNR_int(j) = snr(Sacc(:,j),modes(:,j)-Sacc(:,j));
    end
 
    if nb == 1
     SNR_modes(p,:) = SNR_int;
    else
     SNR_modes(p,:) = SNR_modes(p,:) + SNR_int;
    end
   end   
  end
 end
 SNR_modes = SNR_modes/nbreal;
end
