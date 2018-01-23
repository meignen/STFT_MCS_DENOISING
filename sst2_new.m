function [STFT,SST,VSST1] = sst2_new(s,sigma,Nfilt,Nfft,gamma)

 %% sst2_new : computes the STFT of a signal and different versions of synchrosqueezing
 %
 % INPUTS:   
 %   s: real or complex signal
 %   sigma: the variance of the Gaussian window
 %   Nfft: number of frequency bins
 %   gamma: threshold on the STFT for reassignment
 % OUTPUTS:   
 %   STFT : the short-time Fourier transform
 %   SST  : standard synchrosqueezing
 %   VSST1: vertical second-order synchrosqueezing [1]
 % REFERENCES
 % [1] Behera, R., Meignen, S., & Oberlin, T. (2015). Theoretical Analysis
 % of the Second-order Synchrosqueezing Transform. To appear in ACHA
 
 s = s(:);
           
 ft   = 1:Nfft;
 bt   = 1:length(s);
 nb   = length(bt);
 neta = length(ft);
 
 prec = 10^(-3);
 L =  sigma*Nfilt;
 l = floor(L*sqrt(-log(prec)/pi))+1;
 g = amgauss(2*l+1,l+1,L);
 
 % Window definition
  
  n   = (0:2*l)'-l;
  t0  = n/Nfilt;
  t0  = t0(:);
  a   = pi/sigma^2;
  gp  = -2*a*t0.*g; 
  gpp = (-2*a+4*a^2*t0.^2).*g; % g''

 % Initialization
 STFT  = zeros(neta,nb);
 SST   = zeros(neta,nb);
 VSST1 = zeros(neta,nb);
 
 omega  = zeros(neta,nb);
 tau    = zeros(neta,nb);
 omega2 = zeros(neta,nb);
 phipp  = zeros(neta,nb);
             
 %% Computes STFT and reassignment operators

 for b=1:nb
 	% STFT, window g  
 	time_inst = -min([l,b-1]):min([l,nb-b]);
    tmp = fft(s(bt(b)+time_inst).*g(l+time_inst+1),Nfft);
 	vg  = tmp(ft);
     
 	% STFT, window xg           
 	tmp = fft(s(bt(b)+time_inst).*(time_inst)'/Nfilt.*g(l+time_inst+1),Nfft);
 	vxg = tmp(ft);
       
    % operator Lx (dtau)
	
    tau(:,b) = vxg./vg;
 	
    % STFT, window gp
 	tmp = fft(s(bt(b)+time_inst).*gp(l+time_inst+1),Nfft);
 	vgp = tmp(ft);
    
    % operator omega
    omega(:,b) = Nfilt/Nfft*(ft-1)'-real(vgp/2/1i/pi./vg);    
 	
    
    % STFT, window gpp
 	tmp  = fft(s(bt(b)+time_inst).*gpp(l+time_inst+1),Nfft);
 	vgpp = tmp(ft);
       
    %STFT, windox xgp
 	tmp  = fft(s(bt(b)+time_inst).*(time_inst)'/Nfilt.*gp(l+time_inst+1),Nfft);
 	vxgp = tmp(ft);
    
       
 	%computation of the two different omega 
        
    phipp(:,b) = 1/2/1i/pi*(vgpp.*vg-vgp.^2)./(vxg.*vgp-vxgp.*vg);
       
    %new omega2
    omega2(:,b) = omega(:,b) - real(phipp(:,b)).*real(tau(:,b))...
                              + imag(phipp(:,b)).*imag(tau(:,b)); 

	% Storing STFT       
    STFT(:,b) = vg.*exp(-2/Nfft*pi*1i*(ft-1)'*time_inst(1))/Nfft*g(l+1);%renormalized so that it fits with recmodes
 end
  
 %% reassignment step
 for b=1:nb
    for eta=1:neta
        if abs(STFT(eta,b))> gamma/Nfft
           k = 1+round(Nfft/Nfilt*omega(eta,b));
            if (k >= 1) && (k <= neta)
             % original reassignment
             SST(k,b) = SST(k,b) + STFT(eta,b);
            end
            %reassignment using new omega2
            k = 1+round(Nfft/Nfilt*omega2(eta,b));
            if k>=1 && k<=neta
                % second-order Vertical reassignment: VSST
                VSST1(k,b) = VSST1(k,b) + STFT(eta,b);
            end 
        end
    end
 end
end