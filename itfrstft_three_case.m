function [x] = itfrstft_three_case(tfr,cas)
 %tfr  : STFT of signal x
 %N  : number of frequency bins
 %cas: if 1, no assumption on h except it does not vanish at 0 
 %     if 2, use a filter with unit energy 
 %     if 3, use a filter with unit mean
 
 %x : restored signal
 
 [N,xrow] = size(tfr);
 
 hlength=floor(N/4);
 hlength=hlength+1-rem(hlength,2);%the length of the filter has to be odd
 h=tftb_window(hlength,'HANNING');
 
 [hrow,hcol]=size(h); 
 Lh=(hrow-1)/2; 
 
 if cas == 3, 
  h = h/sum(h);
 end 
 if cas == 2,
  h=h/norm(h);
 end
 
 if (cas == 1)
  %case without periodizing
  x = zeros(xrow,1);   
  for icol=1:xrow,
   x(icol) = 1/h(Lh+1)*mean(tfr(:,icol));
  end
 end
 
 if (cas == 2)
  %case with periodization and with unit energy filters   
  tau =-Lh:Lh;
  x = zeros(xrow,1);    
  for icol = 1:xrow,    
   if (icol > Lh)&& (icol <= xrow-Lh)
    x(icol) = mean((tfr(:,icol-tau).*exp(2*1i*pi*(0:N-1)'*tau/N))*h(Lh+1+tau));
   else
    x(icol) = mean((tfr(:,1+rem((icol-1)-tau+xrow,xrow)).*exp(2*1i*pi*(0:N-1)'*tau/N))*h(Lh+1+tau));     
   end 
  end
 end
 
 if (cas == 3)
  %case with periodization and with unit mean filters   
  tau =-Lh:Lh;
  x = zeros(xrow,1);
  for icol = 1:xrow,    
   if (icol > Lh) && (icol <= xrow -Lh)
    x(icol) = mean(sum(tfr(:,icol-tau).*exp(2*1i*pi*(0:N-1)'*tau/N),2));
   else
    x(icol) = mean(sum(tfr(:,1+rem((icol-1)-tau+xrow,xrow)).*exp(2*1i*pi*(0:N-1)'*tau/N),2));     
   end
  end  
 end
end