function [tfr,norm2h] = tfrstft_three_case_down(x,N,cas,h,Lh,downsamp,shift)
 
 %x        : signal
 %N        : number of frequency bins
 %cas      : if 1, no assumption on h except it does not vanish at 0 
 %           if 2, use a filter with unit energy 
 %           if 3, use a filter with unit mean 
 %h        : the filter h used
 %Lh       : the filter is of length 2Lh+1
 %downsamp : the downsampling value between 1 and Lh.
 %shift    : the shift parameter, the value used for reconstruction are
 %           m*downsamp+shift (if downsamp = 1, shift =0)
 
 %tfr      : short time Fourier transform
 %norm2h   : the L2 norm of the filter on it support
 
 [xrow,xcol] = size(x);
 
 t = 1+shift:downsamp:xrow; %the time instant, we consider the time instant shitfed by a factor shift.
  
 tfr= zeros (N,length(t)) ;
 if (cas == 1)
  %case without periodizing   
  trans  = zeros(1,length(t));
  norm2h = zeros(1,length(t));
  
  for icol=1:length(t),
   tau = -min([Lh,t(icol)-1]):min([Lh,xrow-t(icol)]); 
   tfr(1:length(tau),icol) = x(t(icol)+tau,1).*h(Lh+1+tau);
   trans(icol)  = tau(1);
   norm2h(icol) = norm(h(Lh+1+tau)); %we compute the L2 norm of the filter on its support 
  end
  tfr=fft(tfr,N); 
  A = exp(-2/N*pi*1i*(0:N-1)'*trans);
  tfr = tfr.*A;
 end
 
 if (cas == 2)||(cas == 3)
  %case with periodization    
  tau = -Lh:Lh;
  for icol = 1:length(t), 
   if (t(icol) > Lh) && (t(icol) <= xrow-Lh)
    tfr(1:length(tau),icol) = x(t(icol)+tau,1).*h(Lh+1+tau); 
   else
    tfr(1:length(tau),icol) = x(1+rem((t(icol)-1)+tau+xrow,xrow),1).*h(Lh+1+tau);
   end
   norm2h(icol) = norm(h(Lh+1+tau));%computation of the L2 norm of h 
  end
   tfr = fft(tfr,N);
   trans = Lh*ones(1,length(t)); 
   A = exp(2/N*pi*1i*(0:N-1)'*trans);
   tfr = tfr.*A;
 end
end 
