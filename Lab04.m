clear all;
close all;
xbit=[1 0 1 1 0 1 0 0 0 1 1 0];
% Initial reference bit is assumed to be 1
% Binary bit stream is in 0 and 1 : 12 bits
% NOT of Exclusive OR operation
difencod(1)=~(1-xbit(1));
for i=2:length(xbit)
    difencod(i)=~(difencod(i-1)-xbit(i));
end

% Differential Encoded binary bit stream
xbit(1)=1-~(difencod(1));
for i=2:length(xbit)
    xbit(i)=difencod(i-1)-~(difencod(i));
    if(xbit(i)==-1)
        xbit(i)=1;
    end
end

% Inphase unipolar bit stream 
% from differentially encoded baseband
for i=1:2:(length(difencod)-1)
    inp(i)=difencod(i);
    inp(i+1)=inp(i);
end

% Quadrature unipolar bit stream 
% from differentially encoded baseband
for i=2:2:(length(difencod))
    qp(i)=difencod(i);
    qp(i-1)=qp(i);
end

% Inphase bipolar NRZ bit stream 
for i=1:(length(inp))
    if(inp(i)== 1)
        it(i)=1;
    elseif(inp(i)==0)
        it(i)=-1;
    end
end

% Quadrature bipolar NRZ bit stream 
for i=1:(length(qp))
    if(qp(i)== 1)
        qt(i)=1;
    elseif(qp(i)==0)
        qt(i)=-1;
    end
end

% Raised Cosine Filter used
filtorder = 40; % Filter order
nsamp=4;
delay = filtorder/(nsamp*2);
rolloff = 0.5; % Rolloff factor of filter
rrcfilter = rcosine(1,nsamp,'fir/normal',rolloff,delay);

% Plot impulse response.
figure(1); 
impz(rrcfilter,1);
grid on
%title(' Impulse response of Raised Cosine Filter');

%% Transmitted Signal
% Upsample and apply  raised cosine filter.
itx = rcosflt(it,1,nsamp,'filter',rrcfilter);

Drate=64000; % Bit rate
T=1/Drate;
Ts=T/nsamp;
time=0:Ts:(length(itx)-1)*Ts;

figure(2); 
plot(time,itx)
%title(' Low pass filtered InPhase Component');
xlabel('Time(sec)');
ylabel('Amplitude(volt)');
grid on

tme=Ts:Ts:(length(itx)-1)*Ts+Ts;
qtx = rcosflt(qt,1,nsamp,'filter',rrcfilter);

figure(3); 
plot(tme,qtx)
title(' Low pass filtered Quadrature Component');
xlabel('Time(sec)');
ylabel('Amplitude(volt)');
grid on

fc=900*100000; % 900MHz Carrier frequency chosen
dd=2*pi*fc*time';
ddd=2*pi*fc*tme';

% One bit or 1/2 of symbol delay consideration in OQPSK
delay(1:nsamp)=0.0;
delay((nsamp+1):length(qtx))=qtx(1:(length(qtx)-nsamp));

half=filtorder/2;
mt=(cos(dd)).*itx+(sin(ddd)).*delay';

figure(4); 
plot(time,mt)
xlabel('Time(sec)');
ylabel('Amplitude(volt)');
title(' Differentially encoded OQPSK modulated signal');
grid on

snr=10;
% Signal-to-noise ratio per sample is assumed to be 10
madd=awgn(mt,snr);

figure(5); 
plot(time,madd)
grid on
xlabel('Time(sec)');
ylabel('Amplitude(volt)');
%title(' Differentially encoded OQPSK modulated signal with added white noise');

cscomp=mt.*(cos(dd));
sincomp=mt.*(sin(ddd));

figure;
plot(time,cscomp)
grid on
xlabel('Time(sec)');
ylabel('Amplitude(volt)');

lpfin = rcosflt(cscomp,1,nsamp,'filter',rrcfilter);
lpfqu = rcosflt(sincomp,1,nsamp,'filter',rrcfilter);

tmx=0:Ts:(length(lpfin)-1)*Ts;
tmy=Ts:Ts:(length(lpfqu)-1)*Ts+Ts;

figure(5); 
plot(tmx,lpfin)
grid on
xlabel('Time(sec)');
ylabel('Amplitude');

figure(6); 
plot(tmy,lpfqu)
grid on
xlabel('Time(sec)');
ylabel('Amplitude(volt)');

% Initial checking for I and Q channel bit stream
itxx=itx(half:nsamp:length(xbit)*nsamp+half-1);
for i=1:1:length(itxx)
    if(itxx(i)> 0)
        chk1(i)=1;
    elseif(itxx(i)< 0)
        chk1(i)=-1;
    end
end

ityy=qtx(half:nsamp:length(xbit)*nsamp+half-1);
for i=1:1:length(ityy)
    if(ityy(i)> 0)
        chk2(i)=1;
    elseif(ityy(i)< 0)
        chk2(i)=-1;
    end
end

disp('I channel bit stream checking')
distortion = sum((it-chk1).^2)/length(chk1); % Mean square error
disp(distortion)

disp('Q channel bit stream checking')
distortion = sum((qt-chk2).^2)/length(chk2); % Mean square error
disp(distortion)

% Differentially decoded bit stream from I and Q channels
for i=1:2:(length(xbit)-1)
    dfd(i)=chk1(i);
end
for i=2:2:(length(xbit))
    dfd(i)=chk2(i);
end
for i=1:(length(xbit))
    if(dfd(i)== 1)
        dfdecod(i)=1;
    elseif(dfd(i)==-1)
        dfdecod(i)=0;
    end
end

detected(1)=1-~(dfdecod(1));
for i=2:length(xbit)
    detected(i)=dfdecod(i-1)-(~dfdecod(i));
    if(detected(i)==-1)
        detected(i)=1;
    end
end

disp('Distortion between transmitted and received NRZ bit stream')
distortion = sum((xbit-detected).^2)/length(detected); % Mean square error
disp(distortion)

tmx=0:(1/64000):(1/64000).*(length(xbit)-1);

figure(7);
subplot(211)
stairs(tmx, xbit, 'b','LineWidth',1.5)  % blue color for transmitted
set(gca,'ytick',[0 1])
grid on
xlabel('Time(sec)');
ylabel('Binary value');
title(' Transmitted bit stream ');

subplot(212)
stairs(tmx, detected, 'r','LineWidth',1.5) % red color for received
set(gca,'ytick',[0 1])
grid on
xlabel('Time(sec)');
ylabel('Binary value');
title(' Received bit stream ');
