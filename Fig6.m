
[X,XX,coeffX,Y,YY,coeffY,Z,ZZ,coeffZ]=recons_bat(3);
[X1,XX1,coeffX1,Y1,YY1,coeffY1,Z1,ZZ1,coeffZ1]=recons_bat(4);

d=0:10;
%SNR 0
figure()  
[snr_out,nbcoeff]= snrout_down(3,'hamming',32,0);
[snr_out1,nbcoeff1]= snrout_down(3,'hamming',16,0);
[snr_out2,nbcoeff2]= snrout_down(3,'hamming',8,0);
coeff = [nbcoeff nbcoeff1 nbcoeff2];
snrout = [snr_out snr_out1 snr_out2];
plot(3*(2*d+1),X,coeffX,XX,'--',4*(2*d+1),X1,':',coeffX1,XX1,'-.',coeff,snrout,'-*');

%SNR 10
figure()
[snr_out,nbcoeff]= snrout_down(3,'hamming',32,10);
[snr_out1,nbcoeff1]= snrout_down(3,'hamming',16,10);
[snr_out2,nbcoeff2]= snrout_down(3,'hamming',8,10);
coeff = [nbcoeff nbcoeff1 nbcoeff2];
snrout = [snr_out snr_out1 snr_out2];
plot(3*(2*d+1),Z,coeffZ,ZZ,'--',4*(2*d+1),Z1,':',coeffZ1,ZZ1,'-.',coeff,snrout,'-*');