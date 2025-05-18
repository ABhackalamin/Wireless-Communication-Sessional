clear all;
close all;
% Test with synthetically generated sinusoidal wave
f=1000;% Frequency of the audio signal
Fs =4000; % Sampling rate is 4000 samples per second.
t = [1/Fs:1/Fs:1];% total time for simulation=0.05 second.
% Number of samples=4000
Am=1.0;
signal = Am*sin(2*pi*1000*t); % Original signal
figure(1);
plot(t(1:200),signal(1:200), 'r') % red color
set(gca,'ytick',[ -1.0   0 1.0 ])
title('A segment of synthetically generated sinusiodal wavform')
grid on
xlabel( 'time(sec)');
ylabel( 'Amplitude(volt)');
maximumvalue=max(signal);
minimumvalue=min(signal);
interval=(maximumvalue-minimumvalue)/255; % interval: 
partition = [minimumvalue:interval:maximumvalue]; % -1:0.0078:1
codebook = [(minimumvalue-interval):interval:maximumvalue]; % -1.0078:0.0078:1
[index,quants,distor] = quantiz(signal,partition,codebook);
% Convertion of deci  into binary from least  to most significant
indxtrn=index';
for i=1:4000
matrix(i,1:1:8)=bitget(uint8(indxtrn(i)),1:1:8);
end,
% matrix is of  4000 rows X 8 columns
% matrixtps is a matrix of 8 rows X4000 columns
matrixtps=matrix';
% Baseband is produced, it has 32000 bits
baseband=reshape(matrixtps,4000*8,1);
Tb=1/32000;
% bit rate 32 kbps
time=[0:Tb:1];
figure(2);
stairs(time(1:500),baseband(1:500), 'g') % green color
title(' A segment of baseband signal')
xlabel('Time(sec)')
ylabel('Binary value')
set(gca,'ytick',[0  1 ])
axis([0,time(500),0,1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Serial the data for the next step.
input_to_Convolutional_encoder = baseband'; % 1 X  32000
%Now, the binary converted data is sent to th Convolutional encoder.
t=poly2trellis(7, [171 133]);
%Channel coding
code = convenc(input_to_Convolutional_encoder,t); % 1 x 64000
%Interleaving
st2 = 4831;
data_interleave = randintrlv(code,st2); % Interleave, 1 row x 64000 columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Quadrature amplitude modulation
M=4;
k=log2(M);
baseband=double(baseband);
% bit to symbol mapping
symbol=bi2de(reshape(data_interleave,k,length(data_interleave)/k).','left-msb');
Quadrature_amplitude_modulated_data = qammod(symbol,M);
% demodulation of Quadrature amplitude data
Quadrature_amplitude_demodulated_data = qamdemod(Quadrature_amplitude_modulated_data,M);
[number,ratio]= symerr(symbol,Quadrature_amplitude_demodulated_data) % symbol error
% symbol to  bit mapping
% 2-bit symbol to Binary bit mapping
Retrieved_bit = de2bi(Quadrature_amplitude_demodulated_data,'left-msb');
Retrieved_bit=Retrieved_bit';
Retrieved_bit=reshape(Retrieved_bit, 64000,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deinterleaving
errors = zeros(size(Retrieved_bit));
inter_err = bitxor(Retrieved_bit,errors); % Include burst error.
data_deinterleave=randdeintrlv(inter_err,st2); 
%Convolutional Decoding
tblen=3;
decodx= vitdec(data_deinterleave,t,tblen,'cont','hard'); %
N3=length(decodx);
NN=N3/8;
decod2(1:(N3-3))=decodx(tblen+1:end);
decod2(N3)=decodx(1);
decod2=decod2' ; % 32000 X 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseband=double(baseband);
[number,ratio]= biterr(decod2,baseband)
convert=reshape(decod2,8,4000); % First reshaping and then transposing
matrixtps=double(matrixtps);
[number,ratio]= biterr(convert,matrixtps)
convert=convert' ; % 4000 rows X 8 columns
% binary to decimally converted value
intconv=bi2de(convert); % converted into interger values(0-255) of 4000 samples
% intconv is 4000 rows X 1 column
[number,ratio]= biterr(intconv,index');
sample_value=minimumvalue +intconv.*interval;
figure(3)
subplot(2,1,1)
plot(time(1:100),signal(1:100), 'b') % blue color
set(gca,'ytick',[ -1.0   0 1.0 ])
axis([0,time(100),-1,1])
title('Graph for a segment of recoded Audio signal')
xlabel('Time(sec)')
ylabel('Amplitude')
grid on
subplot(2,1,2)
plot(time(1:100),sample_value(1:100), 'm') % magenta color
axis([0,time(100),-1,1])
set(gca,'ytick',[ -1.0   0 1.0 ])
title('Graph for a segment of retrieved Audio signal')
xlabel('Time(sec)')
ylabel('Amplitude')
grid on
