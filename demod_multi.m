function [sp1_s,integ1] = demod_multi(s,VSST_s,Nfft1,t,N,nr,lambda,beta,jump)   
 
 %computation of the different ridges
 [Cs2] = exridge_mult(VSST_s,nr,lambda,beta,jump);
  
 integ1 = zeros(size(Cs2));
 sp1_s  =  zeros(nr,length(s));
 
 for k = 1:nr
  %numerical integration of R_vsst0 and demodulated signal
  integ1(k,:) = cumtrapz(t,N/Nfft1*(Cs2(k,:)-1));
  sp1_s(k,:)  = transpose(s).*exp(-2*1i*pi*(integ1(k,:)-100.*t));
 end 