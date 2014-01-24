clear

n = 0:314;
x = cos(pi/8*n)+cos(pi/2*n)+sin(3*pi/4*n);%+cos(2*pi*n)+cos(4*pi*n);
% Design an FIR equiripple bandpass filter to remove the lowest and highest discrete-time sinusoids.

d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1/4,3/8,5/8,6/8,60,1,60);
Hd = design(d,'equiripple');
% Apply the filter to the discrete-time signal.

y = filter(Hd,x);
freq = 0:(4*pi)/length(x):pi;
xdft = fft(x);
ydft = fft(y);
plot(abs(xdft(1:length(x)/2+1)));
hold on;
plot(abs(ydft(1:length(x)/2+1)),'r','linewidth',2);
legend('Original Signal','Bandpass Signal');
% figure
% plot(n,x);