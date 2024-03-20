% B {Pailtjorpe, U Sydney. sixNMclusterVx.m, (24/6/21 - 24)
 % code to set up stimuli for NM models; do FFT, test S[V] forms, etc.
% Wave & Pulse density set up; FFT analysis; Autocorrln fn (24/6/21)
%  #1.3 Waves: set up: options for waveform
%  #1.4 Pulse density fn
%  #1.6 Load spike rate sequence - from sim oof 10^2-10^3 LIF neurons

%  #1.4.2 >>  & add repeat pulse[s]: periodic

%   2. sigmoid switching fn. : v -> pulse rate
%  #3.0 FFT of signals
%  #3.3 Autocorrln
%   Appx - various: tests of codes & functions, S[v(t)] etc

%% time axis & IC (0 mV)
x0=0; xf=5.0; % t range (sec)
h = 0.001;   % t step (1 ms; [& 0.5, 0.25 ms, for debug])
x = (x0:h:xf);  % t domain array (ms)
time =x/(1000);  %time =x/(1000*h); % time (sec) as row vec

%% 1.3 set up Wave train/ packet type (theta, alpha, beta, gamma)
  % cf. code: VisMarmosetBrainWaveV1.m
fprintf(' >> Add types of wave train: sin wave @ 5.3 Hz; start at 4.0 s ')
 %fprintf('    multi modes: k= 1, 2, 4 ..') 
Vyt =zeros(1,length(x)); % "x" is time vector - set in main NM code
% parameters
Vwm= 1.0  % amplitude
lambda = 41.2; % mm  %lamda = 0.03 (m) ~ brain size (marmoset ~30mm)

% iterate freq;   k /lambda
   %wfreq = 5.0/10000 % debug
wfreq = 5.3/10000 % %wfreq = 5.9/10000; %     theta [#2 or 5]: 5.3 or 5.9Hz
% wfreq = 9.9/10000  % 10.7/10000; %  alpha [for #4
% wfreq = 16.0/10000 % (beta) % was 16.1
   %wfreq = 13.0/10000 % off resonance (midway betw alpha & beta) 
   %wfreq = 7.9/10000 % off resonancemid(theta, alpha)%wfreq = 16.0/10000 % beta
%wfreq = 32.1/10000 % (gamma) wave frequency, Hz (s^-1) {was 32.1, 32.3 
  %wfreq = wfreq/2 % nb abs effectively doubles freq of envelope - for rectified
  %wfreq = wfreq*2 %  doubles freq, for half wave rectifcation
   %wfreq = wfreq*3
phi=0; % phase delay (lag, in radians) % choose phase to line up start of sin wave
         %  & antinode near RoI % for alpha: S[V~erfc], fo, fo/2, 2fo
dphi = 0  %dphi= -175 % dphi =0    % (in t-steps: 0.1ms) extra phase diff to align with delta-y osc
%pdhi=dphi*pi/180; % in rad 
Wspeed = wfreq*lambda; % (m/s or mm/ms); here 2 m/s
istart=40001  % istart=1 %-100 %  istart=140001 % turn wave on at 1,4 sec; aim near V(t)=0
% istop= 55000 +605 %  then stop: eg after 2 sec [6s-390steps] or [5s+700]; also adjust to V~0
istop = xf/h; % full t span
% row vec x is the t-line axis [0 :1ms steps: end (sec)]
 tline = x*10000.0; % convert ms scale to sec, to match pulses & DE solver
  %max(tline)
 ypt =16.0; % y[A-P] coord (mm), near antinode
 Vyt0 = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*tline -phi +dphi); % base case;  row vec
fprintf('+ wave type: 1) cont. 2). envelope.  3) modulated.  4) packet.  5) on/off: \n')
 
% > > > > > > > > > > > > > 
 WaveType = 1  % is: 1 or 2 or 3 or 4 or 5.
 
% choose which to calculate:
 switch WaveType
    case 1
 % 1.3a  continuous wave [& on/off ?]
fprintf('  + add continuous wave train at f-i, above x 1; at 4-20 s  \n')
Vyt =zeros(1,length(x));
 %Vyt = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*tline -phi +dphi); % now row vec
  %figure; plot(Vyt, 'r-');  title('stimulus  wave'); xlabel('t (msec) ')
  % V osc at fixed y, vs t (ie "x"), in 2D
 %Vyt(istart:end) = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*(tline(1:end-istart+1) -phi +dphi) ); % start at yy sec
Vyt(istart:istop) = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*(tline(1:istop-istart+1) -phi +dphi) ); % start at yy sec, then stop

    case 2
% 1.3b  enveloped simulus wave train [envelope of wave sizes, at const rate]
% envelope of the wave; at double period (ie half cycle)
fprintf('  + add enveloped wave train: gamma x 1 @ theta modulation  \n')
xb=1/5.9;  % theta % xb=1; %period of burst modulation (sec); 1/f-theta
%xb=1/16.1; %beta
 xb=2*xb; % for f/2 (rectified)
  %xb=xf; % full time period
 Wenvelope=1.0*sin(2*pi*x/(2*xb));
  figure; plot(Wenvelope); title('envelope') % debug
 Vyt=Wenvelope.*Vyt0;
  %figure; plot(Vyt, 'b-'); title('stimulus enveloped wave train'); xlabel('t (msec) ')

    case 3
% 1.3c  modulated simulus wave train [envelope of wave sizes, at const rate]
% envelope of the wave; at double period (ie half cycle)
fprintf('  + add modulated wave train: 32.8 Hz x 1 mV  \n')
xb=1/5.9;  % xb=0.1; %period of burst modulation (sec)
  %xb=xf; % full time period
 Wenvelope=1.0*sin(2*pi*x/(2*xb));
  %figure; plot(Wenvelope) % debug
 Vyt=Wenvelope.*Vyt0;
  %figure; plot(Vyt, 'b-');  title('stimulus modulated wave bursts'); xlabel('t (msec) ')
   % appears to be same as "2" ?

    case 4 
% 1.3d  "beta burst" simulus wave packet [envelope of wave sizes, at const rate]
% envelope of the wave; say 1 sec duration (ie half cycle); start at t=1 sec.
fprintf('  + add single wave packet: theta, 5.9 Hz x 1 mV; @ 4s  \n')

 Vytburst=zeros(1,10000); % set up row vec (1sec @ 0.1 ms)
 %xb=1/5.9; %theta %xb=1; %period of burst modulation (sec)
 xb=1/16.1; %beta
 xb=2*xb; % for beta/2 (rectified)
 Wenvelope=1.0*sin(2*pi*x(1:10000)/(2*xb)); % half wave in 1 sec t window
  %figure; plot(Wenvelope) % debug
  
 Vytburst=Vyt0(10001:20000).*Wenvelope; % nb. here Wenvelope has 1k/10k elements
 %Vytburst=Vyt0.*Wenvelope;
  %figure; plot(Vytburst, 'b-'); title('stimulus single wave burst'); xlabel('t (msec) ')
 Vyt(40001:50000)=Vytburst;  
   figure; plot(x,Vyt); xlim([3 6])
   title('Wave: gamma on theta, burst')
  
     case 5 
% 1.3e  on/off sequence of waves
  % sq wave envelope of the wave; a la TMS protocols
fprintf('  + add sequence of wave: 5Hz sq wv modn of pos 16.25 Hz x 10 mV; @ 4s : one burst  \n')
  Vyt =zeros(1,length(x));
  Vyt(istart:end) = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*(tline(1:end-istart+1) -phi +dphi) ); % start at yy sec
  %Vyt(istart:istop) = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*(tline(1:istop-istart+1) -phi +dphi) ); % start at yy sec, then stop
 % Monophasic TMS: pos pulse only (cf. dA/dt of coils in Laasko('13)  % induces max 1 ?
  Vyt(find(Vyt<0))=0; % figure; plot(Vyt)
  mod_freq = 5.00; % modulation frequency,  Hz (s^-1) % nb. also need to change npulse, to get n sec long stimulus
  % set up modulating [pos only] sq wave
  Ips=zeros(1,length(x));
    %tstart= 1; % tstart   set above 
   %tdelay= round(1/(h*mod_freq)) % period: needs to be integer, for indexing; allow h varies
 twidth=2000/2; % width of half pulse (dt steps; eg 0.1 ms)
 dphi= 0 % phase delay (dt steps; ms) of 1st pulse [& pulse train]
         % or increment (ramp) for ea pulse; or for 2nd pulse; etc.
  npulse=10 % pulse train % set length of pulses: eg 1 sec
  for i=0:npulse-1
   Ips(istart+ 2*i*twidth: istart +2*i*twidth +twidth) =1.0; % infill pulse
  end
    %figure; plot(x, Ips); xlim([4 6]); ylim([0 1.1]) % debug
    %fprintf('  + add sequence of wave: 5Hz sq wv modulation of const 1 mV')
   
  Vyt=Ips.*Vyt; %max(Vyt)
  figure; plot(x, Vyt, 'r-'); xlabel('t (0.1 msec) '); xlim([4 6]);
  
    otherwise
     fprintf('/n >>>>>  choose wave type! ')
 end
  clear  dphis tdelay twidth tline mod_freq npulse Ips lambda W* wfreq

% 1.4  Gather stimuli: to pass to DE solver
 % pran= [pran, 0]; % need N+1 elements (eg. 5001)
 figure;  subplot(2,1,1);  title('stimulus pulse train'); hold on
 subplot(2,1,2); plot(x, pran); hold on; plot(x, Vyt, 'r-'); xlabel('t (0.1 msec) ')
   %plot(x, Vyt, 'r-'); xlabel('t (sec) ')
  %xlim([0 20000])
 title(' waves + noise ');

 %V3 = Vyt; % accumulate multiple modes

 clear xb Vyt0 Wenvelope tline Vytburst dphi* phi W* Vwm wfreq lambda ypt
% >>>>>>>>>>>>>> 

%% 1.3a set up 2nd Wave train
fprintf('  + add 2nd wave, gamma; on at 10s \n')
%wfreq = 5.3/10000 %wfreq = 5.9/10000; %   theta [#2 or 5]: 5.3 or 5.9Hz
% wfreq = 9.9/10000  % 10.7/10000; %  alpha [for #4
% wfreq = 16.0/10000 % (beta) % was 16.1
wfreq = 32.1/10000 % (gamma) wave frequency, Hz (s^-1) {was 32.1
phi=0; % phase delay (lag, in radians) % choose phase to line up start of sin wave
         %  & antinode near RoI % for alpha: S[V~erfc], fo, fo/2, 2fo
dphi = 0 %500  %dphi= -175 % dphi =0    % (in t-steps: 0.1ms) extra phase diff to align with delta-y osc
%pdhi=dphi*pi/180; % in rad 
% parameters
Vwm= 1.0  % amplitude
lambda = 41.2; % mm  %lamda = 0.03 (m) ~ brain size (marmoset ~30mm)
Wspeed = wfreq*lambda; % (m/s or mm/ms); here 2 m/s
istart= 100001  % istart=1 %-100 %  istart=140001 % turn wave on at 1,4 sec; aim near V(t)=0
istop = xf/h; % full t span
 % row vec x is the t-line axis [0.1ms steps: end (sec)]
 tline = x*10000.0; % convert ms scale to sec, to match pulses & DE solver
 ypt =16.0; % y[A-P] coord (mm), near antinode
 %Vyt0 = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*tline -phi +dphi); % base case;  row vec
 % 1.3a  continuous wave [& on/off ?]
fprintf('  + add 2nd wave train at f-i, above x 1; at 10.0- s  \n')
Vyt2 =zeros(1,length(x));
Vyt2(istart+dphi:istop-dphi) = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*(tline(1:istop-istart -2*dphi +1) )); % start at yy sec, then stop
% Pulse1(istart+dphi: istart+dphi+2000) = Pulse(1:2000+1); % correct shift 

figure; plot(Vyt2); xlim([99000 102000]) % debug

clear Vyt0 tline  phi W* Vwm wfreq lambda ypt % dphi

%% 1.3b set up Wave train (@ f/2, f, f*2)
  % cf. code: VisMarmosetBrainWaveV1.m
fprintf(' >> Add continuous wave train: cos wave @ 32.1 Hz x 1 mV; start at 4 s ')
 %fprintf('    multi modes: k= 1, 2, 4 ..') 
Vyt =zeros(1,length(x)); % "x" is time vector - set in main NM code
% parameters
Vwm= 1.0 % mV % amplitude
lambda = 41.2; % mm  %lamda = 0.03 (m) ~ brain size (marmoset ~30mm)

% iterate freq;   k /lambda
%wfreq = 16.1/10000 % theta, alpha, beta
wfreq = 32.1/10000 % 3*beta/2; or 2*gamma
 %wfreq = 32.1/10000 % (gamma) wave frequency (0.1 ms^-1)  %freq = 32.0; %  Hz (s^-1)
wfreq = wfreq/2 % nb abs effectively doubles freq of envelope
   %wfreq = wfreq*3  % 3rd overtone
phi=0; % phase delay (lag, in radians) % choose phase to line up start of sin wave
         %  & antinode near RoI
          % for alpha: S[V~erfc], fo, fo/2, 2fo
dphi = -180  % phase diff to align with sin, using cos as base wave
pdhi=dphi*pi/180; % in rad 
Wspeed = wfreq*lambda; % (m/s or mm/ms); here 2 m/s
istart=40001  % istart=1 %-100 %  istart=140001 % turn wave on at 1,4 sec; aim near V(t)=0
 %istop= 55000 +605 %  then stop: eg after 2 sec [6s-390steps] or [5s+700]; also adjust to V~0
istop = xf/h; % full t span
 % row vec x is the t-line axis [0 :1ms steps: end (sec)]
 tline = x*10000.0; % convert ms scale to sec, to match pulses & DE solver
 ypt =16.0; % at fixed point, y[A-P] coord (mm), near antinode
 Vyt0 = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*tline -phi +dphi); % base case;  row vec
%  continuous wave - in time
Vyt =zeros(1,length(x));
 %Vyt = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*tline -phi +dphi); % now row vec
  %figure; plot(Vyt, 'r-');  title('stimulus  wave'); xlabel('t (msec) ')
  % V osc at fixed y, vs t (ie "x"), in 2D
Vyt(istart:istop) = Vwm*cos(2*pi*ypt/lambda - 2*pi*wfreq*(tline(1:istop-istart+1) -phi -dphi) ); % start at yy sec, then stop

 figure;  subplot(2,1,1);  title('stimulus wave train'); hold on
 plot(x, Vyt, 'r-'); xlim([4, 5]); grid on
 subplot(2,1,2); plot(x, pran); hold on; plot(x, Vyt, 'r-'); xlabel('t (sec) ')
   %plot(x, Vyt, 'r-'); xlabel('t (sec) ')
  %xlim([0 20000])
 title(' waves + noise ');

 clear xb Vyt0 Wenvelope tline Vytburst dphi* W* Vwm wfreq lambda ypt
%% 1.3c  Multi modes wave
 figure;  subplot(2,1,1);  title('indiv modes'); hold on; xlim([0 1])
 plot(x, V1); plot(x, V2); plot(x, V3);
 subplot(2,1,2);  hold on; xlabel('t (sec) '); xlim([0 1])
 plot(x, V1+V2+V3);
  %plot(x, V1+V2/2+V3/4);
  title('sum 3 harmonics')
  
  Vyt =V1+V2+V3; % for sim code input

%% 1.4 + Pulse density fn, as stimulus input (orig.) - First one
 % following JansenZR('93) & JR('95) / export Pulse1
   % params
   fprintf('  + add 1st [narrow/wide] pulse density stimulus \n')
nep =7; % exponent
q=0.26;  %0.26  %q=0.5; % scale
w=0.004  %w=0.005 % w=0.002; % w=0.005 % 
 %w=0.01; % rate time scale (sec; or 5 ms)
% a) start at t=0
Pulse = q*(x/w).^nep.*exp(-x/w); % nb x is vec of times in sec
 % max is 375
  %figure;  subplot(2,1,1); plot(x, Pulse); title('stimulus pulse'); hold on
  %xlabel('t (sec) ')
  %subplot(2,1,2); plot(x, Pulse); axis([0 0.2 0 400])

 % 1.4.1 single pulse @ start time: 
  % use use of "(x-tstart)" causes divg to -inf
  fprintf('  + 1st single pulse \n')
istart = 40001 % at 4.0 sec; (dt 1 0r 0.1 msec, in array) % check x(10001) or  x(40001)
%istart = 83001 % at 7sec; 8.3s, etc

dphi= 250 %0  %420 %dphi= 337    % phase delay (in 0.1ms t steps) : so aligns w lfp(t) ~ need trial& error
Pulse1 =zeros(1,length(x)); 
Pulse1(istart+dphi: istart+dphi+2000) = Pulse(1:2000+1); % load 1st pulse, wt 1

% Option: add const backgound stim; eg p-bar = 155(s^-1), ref.  Malagarriga('19)
 %Pulse1 = Pulse1+50.0;  %Pulse1 = Pulse1+155.0; 
   % decay:  goes to 10^-7 at 1k steps10^-72 at 5k steps
 figure;  subplot(2,1,1); plot(x, Pulse); title('basic stimulus pulse'); hold on
 xlabel('t (sec) ');  % xlim([0 0.5]) % check pulse width
 subplot(2,1,2); plot(x, Pulse1); title ('delayed  pulse') %axis([0.8 1.2 0 400]); 
  xlim([4.0 4.3]); grid on; xlabel('t (sec) ');
  
 % Integrated "energy" content of pulse[s]
PulseIntegrl = sum(Pulse)*h; % sum across rows % trapezoidal rule 
PulseIntegrl=PulseIntegrl'

%Pulse2 =zeros(1,length(x)); % default: no 2nd pulse
 clear   idelay PulseIntegrl % nPulse istart dphi
 
 %% nb v5d code uses 1 row for Pd
  % v6d code uses 6 rows, for 6NM
 
 %% 1.4.01 + sharper Pulse density fn, - First one
 % following JansenZR('93) & JR('95) / export Pulse1
   % params
   fprintf('  + add [narrow, sharp] pulse density stimulus \n')
nep =7 % exponent
 %q=0.26;  %0.26  %
q=0.8; % scale
w=0.001 % w=0.005; %
 %w=0.01; % rate time scale (sec; or 5 ms)
% a) start at t=0
Pulse = q*(x/w).^nep.*exp(-x/w); % nb x is vec of times in sec
 % max is 375
  %figure;  subplot(2,1,1); plot(x, Pulse); title('stimulus pulse'); hold on
  %xlabel('t (sec) ')
  %subplot(2,1,2); plot(x, Pulse); axis([0 0.2 0 400])

 % 1.4.1 single pulse @ start time: 
  % use use of "(x-tstart)" causes divg to -inf
  fprintf('  + 1st single pulse \n')
istart = 40001 % at 4.0 sec; (dt 1 0r 0.1 msec, in array) % check x(10001) or  x(40001)
%istart = 83001 % at 7sec; 8.3s, etc

dphi= 0 %dphi= 337    % phase delay (in 0.1ms t steps) : so aligns w lfp(t) ~ need trial& error
Pulse1 =zeros(1,length(x)); 
Pulse1(istart+dphi: istart+dphi+2000) = Pulse(1:2000+1); % load 1st pulse, wt 1
 figure;  subplot(2,1,1); plot(x, Pulse); title('delayed stimulus pulse'); hold on
 xlabel('t (sec) ');  % xlim([0 0.5]) % check pulse width
 subplot(2,1,2); plot(x, Pulse1); title ('delayed  pulse') %axis([0.8 1.2 0 400]); 
  xlim([4.0 4.2]); grid on
  
 % Integrated "energy" content of pulse[s]
PulseIntegrl = sum(Pulse)*h; % sum across rows % trapezoidal rule 
PulseIntegrl=PulseIntegrl'

%Pulse2 =zeros(1,length(x)); % default: no 2nd pulse
 clear   idelay PulseIntegrl % nPulse istart dphi 
%% 1.4.02 two pulses 1) @ later start time, eg after 0.2 s:
%idelay=7000 %  /for 0.7s intervals: cf data on in-links to NM#2, 6
idelay=60000 %  /for delay to 10 sec 
wt2 = 1.0;  % relative weight of 2nd pulse [cf. Marm. link wt data]
nPulse=1  % echo
% add 2nd pulse; same phase delay
fprintf('  two pulses \n')
Pulse1(istart+dphi+idelay:istart+dphi+idelay+1000) = wt2*Pulse(1:1000+1); % fill in "later"
%
figure;  subplot(2,1,1); plot(x, Pulse1); title('1st stimulus stream, 2 pulses'); hold on
 xlabel('t (sec) '); xlim([3.9 11])
 subplot(2,1,2); plot(x, Pulse1); xlim([3.9 11]) 
 title ('delayed  pulses')
 PulseIntegrl = sum(Pulse)*h; % sum across rows % trapezoidal rule 
 PulseIntegrl=PulseIntegrl'

 clear wt2 PulseIntegrl 

%% 1.4.1 set up, add 2nd Pulse density fn, at later time
fprintf('  + 2nd single pulse stream \n')
   fprintf('  use 1st pulse functional form \n')
nep =7; % exponent
q=0.26;  %0.26  %q=0.5; % scale
w=0.004  %w=0.005 % w=0.002; % w=0.005 % 
 %w=0.01; % rate time scale (sec; or 5 ms)
 dphi= 0 
% a) start at t=0
Pulse = q*(x/w).^nep.*exp(-x/w); % nb x is vec of times in sec
  %istart2 = istart+7000; % at 0.7s later (eg. for NM@2, after NM#6
%istart2 = 100001 % 2nd PUlse, 10s later; use params as above
istart2 = 100001 % at 10 s [1st sim on at 4 sec
 %dphi= 0  % same phase delay (in 0.1ms t steps) : so aligns w lfp(t) ~ need trial& error
  % use same phase delay as for 1st pulse:
wt2 = 1.0;  % relative weight of 2nd pulse / or set in sim code
Pulse2 =zeros(1,length(x)); 
Pulse2(istart2+dphi: istart2+dphi+1000) = wt2*Pulse(1:1000+1); % load 2nd pulse [params as above]
  figure; subplot(2,1,1); plot(x, Pulse); hold on; %xlim([3.5 5]) % xlim([9.5 11]) % xlim([3.5 5])
title('assemble Pulse fns: Pulse train 1 ') 
subplot(2,1,2); plot(x, Pulse2); title(' Pulse 2')
xlabel('t (sec) '); %xlim([9.5 12])
 clear istart2 wt2 nep q w wt2 dphi %Pulseistart2

 %% 1.4.2 >>  ++ add repeat 1st or 2nd pulse[s]: periodic
 % set up 1st pulse above (at #1.4)
  fprintf('  + N 2nd pulses, regular freq. theta \n')
 %idelay=1000 %  /for 5 Hz rept'n, ie @ 200 ms intervals; ie no gap
 %stim_freq=1/(idelay*h)  % h is dt (0.1 ms)
 stim_freq= 16.3  %6.4 %32.0  %5.74  % (Hz)
 idelay=int32(1/(stim_freq*h)) % need integer index (32 bit)
                                % alt:  gap between pulses; 
 nPulse= 5  % echo, total#
 for i=1:nPulse-1; 
  %Pulse1(istart+dphi+idelay*i:istart+dphi+idelay*i+1000) = Pulse(1:1000+1);
 Pulse2(istart2+dphi+idelay*i:istart2+dphi+idelay*i+1000) = Pulse(1:1000+1); % 2nd pulse
 end
% Option: add const backgound stim; eg p-bar = 155(s^-1), ref.  Malagarriga('19)
% Pulse1 = Pulse1+50.0;  %Pulse1 = Pulse1+155.0;
 
figure; plot(x, Pulse1); % xlim([3.8 4.5])
  % check against last lfp(t) calc
  %figure; plot(lfp(tstart:tstart+4000), 'b');  hold on;
  %plot(Pulse1(tstart:tstart+4000)/10, 'k') % stim % reduce scale [of pulse rate]

  % Integrated "energy" content of pulse[s]
PulseIntegrl = sum(Pulse1)*h; % sum across rows % trapezoidal rule 
PulseIntegrl=PulseIntegrl'

%Pulse2 =zeros(1,length(x)); % default: no 2nd pulse
 clear   idelay PulseIntegrl % nPulse istart dphi
 clear istart dphi idelay nPulse  PulseIntegrl % Pulse
 
%% 1.4.3 Pulse density fn, tighten up
 % following JanserZR('93) & JR('95)
% params
   fprintf('  + add tighter pulse density stimulus \n')
nep =7; % exponent
q=0.5  %0.26  %q=0.5; % scale
w=0.003; % rate time scale (sec; or 5 ms)
% a) start at t=0
Pulse = q*(x/w).^nep.*exp(-x/w); % nb x is vec of times in sec
 % max is 375
 %figure;  subplot(2,1,1); plot(x, Pulse); title('stimulus pulse'); hold on
 %xlabel('t (sec) ')
 %subplot(2,1,2); plot(x, Pulse); axis([0 0.2 0 400])

 % 1.4.1 single pulse @ start time: 
  % use use of "(x-tstart)" causes divg to -inf
  fprintf('  + single pulse \n')
istart = 1000 % (msec, in array)
dphi= -20  % phase delay
Pulse1 =zeros(1,length(x)); 
Pulse1(istart+dphi: istart+dphi+1000) = Pulse(1:1000+1);
 figure;  subplot(2,1,1); plot(x, Pulse1); title('delayed stimulus pulse'); hold on
 xlabel('t (sec) ')
 subplot(2,1,2); plot(x, Pulse1); %axis([0.8 1.2 0 400]); title ('delayed  pulse')
  
%% 1.4.x2 add N pulses 1) @ start time; 2) repeat after 0.2 s (@ 1.3s:
istart = 1000; % (msec, in array)
iwidth=100; % assume fixed
%lfp_stats_atStim= [ mean(lfp_sample), (max(lfp_sample)-min(lfp_sample))] %idelay=200; % for 3.3 Hz rept'n, ie @ 300 ms intervals
idelay=100 % for 5 Hz rept'n, ie @ 200 ms intervals
PulseN=Pulse1; % 1st pulse
nPulse=4  % to add
% add extra pulses
fprintf('  four pulses \n')
for i=1:nPulse
 PulseN(istart+dphi+(iwidth+idelay)*i :istart+dphi+(iwidth+idelay)*i +1000) = Pulse(1:1000+1); % shift along 
end
%
figure;  subplot(2,1,1); plot(x, PulseN); title('delayed stimulus, 5 pulses'); hold on
 xlabel('t (sec) ')
 subplot(2,1,2); plot(x, PulseN); axis([0.9 2.1 0 400]); title('delayed  pulses')
 
%% 1.5 N pulses, irregular sequence, via ISI timings
 fprintf('  + N pulses, irregular ISI pattern: \n')
 % set up 1st pulse above (at #1.4)
 % get t0 & ISI timings - eg. from plot of y-e(1) peaks in sim code at #3.6
 % istart = 40000+t0; (step of h ~ 0.1 ms)
  % use the observed ISI (llcations in 0.1ms steps) beat pattern:
 idelayISI= [1657 1679 1606 1602 1590 1584]; % for 6NM, alpha, S~erfc(1/v, w); p =100;
 AvStim_freq=1/(mean(idelayISI)*h)  %  (Hz)  h is dt (0.1 ms)
 %idelay=int32(tbd?) % need integer index (32 bit)
                                % alt:  gap between pulses; 
 nPulse= 6  % echo, total#
 for i=1:nPulse-1;
  idelay = idelayISI(i); 
  Pulse1(istart+dphi+idelay*i:istart+dphi+idelay*i+1000) = Pulse(1:1000+1);
 end 
 figure; subplot(2,1,1); plot(x, Pulse1); 
 subplot(2,1,2); plot(x, Pulse1); xlim([3.8 5.5]); grid on; title('Pulses at ISI pattern:')

 % clear AvStim_freq idelay nPulse

%% 1.6 Load spike rate sequence, binned - from LIF sim, x modulation
 % cf COBN (Mazzoni) c codes for LIF (e+i) neurons
 fprintf('\n load LIF Firing rate per excit neuron (in 10ms bins) \n')
 % spikes in 10 ms bins: so need to "expand" for 0.1ms resolution here
 % nb raw eFR is for COBN.eNnrn excit neurons, in 0.1ms bins: v small & noisy
   %efrStim50=csvread('LIFefrStim50.csv' ); % read saved file
  %eSR = efrStim50; % get that sim output
eSpikes = zeros(length(x),1); % Array, in 0.1ms (h) steps
for i=0:199  % (10, 20 sec @ 0.01s, 10ms steps: ie 1000 Dt or h)
    i1= i*1000+1; i2 = i*1000+1000;
     %[i1 i2]  % debug
       %eSR(i) = sum(eFR(i1:i2))/0.01; % for net_COBN.eNnrn excit neurons
    eSpikes(i1:i2) = eSR(i+1);  % per excit neuron
     % iSR(i)  % per inhib neuron
end
  %Pd = eSR; % load pulse fn
 % figure; subplot(2,1,1); plot(eSR);  % debug
 %subplot(2,1,2); plot(x, eSpikes);      % debug

  clear efrStim50 

%% 1.6.1 Load spike rate sequence, raw - from LIF sim
fprintf(' load, add in LIF spike train stim [0.1 ms bins] 20sec: \n')
 efrStim=csvread('efrStim6_300.csv' ); % read saved file, from raw LIF sim (M=300, noise 6)
  %efrStim=csvread('efrStim6_300_40s.csv' ); % 40s sim
 efrStim = [efrStim', 0]; % match length to local variables  
  %figure; plot(x, efrStim);      % debug
 
 % envelope of the wave; at double period (ie half cycle)
%fprintf('  + envelope spike train: x alpha 9.6 Hz   \n')
fprintf('  + envelope spike train: x gamma 32.0 Hz   \n')
%fprintf('  + envelope spike train: x 1 sec osc   \n')
Pulse1=zeros(1,length(x)); % set up pulse
 %xb=1/6.25;  % f-theta modulation
 %xb=1/9.6;  % f-alpha modulation
 %xb=1/14.3;  % f-beta modulation
xb=1/32.0;  % f-gamma modulation {32.8 for 4NM, 32.0 for 6NM}
%xb=5.0; %xb=1; %xb=1/6.25;  %xb=1; %period of burst modulation (sec); 
  %xb=xf; % full time period
dphi=40001+180;
 Wenvelope=1.0*sin(2*pi*x/(2*xb) -dphi);
  %figure; plot(Wenvelope); title('envelope') % debug
tmp=abs(Wenvelope).*efrStim; % pos only [rectified] - full time; or on at 4s
Pulse1(40001:end)= 1.0*tmp(40001:end); % 4 sec : full range t; x 1
 %Pulse1(40001:50000) = tmp(40001:50000); % 1sec burst [4:5] sec
%Pulse1(40001:50000) = tmp(40001:50000); % 4 sec burst [4:8] sec
figure; plot(x, Pulse1); title(' Pulse: LIF, modulated x gamma ')
  %figure; plot(Pulse1); title(' Pulse: LIF, modulated x gamma ')
  % title('stimulus enveloped wave train'); xlabel('t (msec) ')
  figure; plot(Wenvelope); xlim([3.9e4 4.1e4]); grid on % debug, check phase
%clear Wenvelope xb tmp 


 %% 1.9 synaptic response fn, for Rate -> V   operation
 % following JansenZR('93) & JR('95) % "opposite" of S[v]
   fprintf('  examine synaptic response kernel fn \n')
% time axis & params
x0=0; xf=5.0; % t range (sec)
h = 0.001;   % t step (1 ms; [&  0.1 ms]
xx = (x0:h:xf);  % t domain array (ms)
c=1.0;  %0.26  %q=0.5; % scale
w=0.010; % dendritic time scale ( 10 ms)
% a) start at t=0
SR = c*(xx/w).*exp(-xx/w); % nb x is vec of times in sec
 % max is 375
 figure;   plot(xx, SR); title('synaptic response kernel fn'); hold on
 xlabel('t (sec) '); axis([0 0.5 0 0.5])
  clear c w x0 xf h xx SR
  
 %% 2. sigmoid switching fn. : v -> pulse rate
  % C2*vm/(1+exp(r*(v0-C1*y(1,10)))) % from JR12a.m, etc
v0=6;    % (mV) midpoint (switching trheshold V) of sigmoid function, for V -> pulse conversion 
vm=5;    % ie. 2*2.5 (s^-1) max firing rate for sigmoid fn
r=0.5;   % activation rate;  Steepness of sigmoid function; default (width/2 is 3ms)
  %r=0.56   % as used by Goodfellow('11), Friston
vv= (-10:1:50); % voltage range % nb. flattens quickly
  %vv=pran; % U or G ran noise
PulseRate=zeros(1,length(vv)); % set up array
P1=zeros(1,length(vv)); % alt form
for i=1:length(vv)
    PulseRate(i)= vm/(1+exp(r*(v0-vv(i))));
    P1(i)= vm/(1+exp(r*(v0-vv(i)))) -vm/2;
end

figure; plot(vv, PulseRate); hold on
plot(vv, P1, 'k--');
plot([v0 v0], [0 vm], 'g--');  % mark switching midpoint
plot([0 15], [vm/2 vm/2], 'g--'); title(' Sigmoid: rate vs V');
xlabel('V (mV) '); ylabel('pulse rate (s^-1) ');

clear PuseRate P1

%% 2.1 S fn for y-vec(1)
 % as calc by NM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
onerow=ones(1,length(y)); % needs vector of 1's for sigmoid (over time span)
% for delta-y(1)
figure; plot(A*a*(0  +(vm./(1+exp(r*(v0*onerow -y(2,:)+y(3,:))) ))) ) % zero stim
  % in range [300 1100]  : v large! nb. A*a = 325
% components
figure; subplot(2,1,1); plot(y(2,:) ); title (' y1(t) '); xlabel('t (ms) ');
 ylabel('y (mV) '); axis([0 5000 0 max(y(2,:))+10]);
subplot(2,1,2); plot(vm./(1+exp(r*(v0*onerow -y(2,:) )) )) % 
axis([0 5000 0 vm+1]);
title (' S[y(t)] '); xlabel('t (ms) '); ylabel('pulse rate (s^-1) ');

% for y-0(1)
%figure; plot(A*a*(0  +(vm./(1+exp(r*(v0*onerow -y(1,:))) ))) ) % in 237, 244] 
%% 2.1a S fn for y-vec(2)
 % as calc byNM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
onerow=ones(1,length(y)); % needs vector of 1's for sigmoid (over time span)
% for delta-y(1)
figure; plot(A*a*(0  +(vm./(1+exp(r*(v0*onerow -y(8,:)+y(9,:))) ))) ) % zero stim
  % in range [300 1100]  : v large! nb. A*a = 325
% components
figure; subplot(2,1,1); plot(y(8,:) ); title (' y2(t) '); xlabel('t (ms) ');
 ylabel('y (mV) '); axis([0 5000 0 max(y(8,:))+10]);
subplot(2,1,2); plot(vm./(1+exp(r*(v0*onerow -y(8,:) )) )) % 
axis([0 5000 0 vm+1]);
title (' S[y(t)] '); xlabel('t (ms) '); ylabel('pulse rate (s^-1) ');

%% 2.2 S fn for ran noise
 % as calc byNM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
onerow=ones(1,length(y)); % needs vector of 1's for sigmoid (over time span)
pran=10.0*rand(1,length(y));
figure; subplot(2,1,1); plot(pran ); title (' y1(t) '); xlabel('t (ms) ');
 ylabel('y (mV) '); %axis([0 5000 0 max(y(2,:))+10]);
subplot(2,1,2); plot(vm./(1+exp(r*(v0*onerow -pran )) )) % 
 %axis([0 5000 0 vm+1]);
title (' S[y(t)] '); xlabel('t (ms) '); ylabel('pulse rate (s^-1) ');

%% 2.3   S[v(t)] for wave train or packets
 % as calc byNM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
 % Ip set up in driver code (twoN, sixCluster, etc.
onerow=ones(1,length(x)); % needs vector of 1's for sigmoid (over time span)
% get Ip and/ot Vy(t):
figure; subplot(2,1,1); plot(Vyt ); title (' Vy(t), gamma 32.1 Hz wave x 10 mV '); xlabel('t (ms) ');
 ylabel('Vy (mV) '); %axis([0 5000 0 max(y(2,:))+10]);
subplot(2,1,2); plot(vm./(1+exp(r*(v0*onerow -Vyt )) )) % 
 %axis([0 5000 0 vm+1]);
title (' S[Vy(t)] '); xlabel('t steps (0.1ms) '); ylabel('pulse rate (s^-1) ');
% 2.3a - plot vs t(ms)
figure; subplot(2,1,1); plot(x, Vyt ); title (' Vy(t), gamma 32.8 Hz wave x 5 mV '); 
xlabel('t (sec) '); grid on
 ylabel('Vy (mV) '); %axis([0 5000 0 max(y(2,:))+10]);
 PulseOut =vm./(1+exp(r*(v0*onerow -Vyt )) ); 
subplot(2,1,2); plot(x, vm./(1+exp(r*(v0*onerow -Vyt )) )) % 
 %axis([0 5000 0 vm+1]);
title (' S[Vy(t)] '); xlabel('t (sec) '); ylabel('pulse rate (s^-1) '); 
grid on; xlim([4 5]) % xlim([4 4.3]) %fast
 % compare/mix LIF output
figure; subplot(2,1,1); plot(x, efrStim); title('LIF output')
subplot(2,1,2); plot(x, (PulseOut+efrStim), 'b'); xlabel('t (sec) '); 

clear PulseOut

%% 2.3.1 S for Ip:
onerow=ones(1,length(x)); % needs vector of 1's for sigmoid (over time span)
figure; subplot(2,1,1); plot(Ip ); title (' Ip(t), 15 Hz pulse train of 10.0 mV '); xlabel('t (ms) ');
ylabel('Ip (mV) '); %axis([0 5000 0 max(y(2,:))+10]);
 %subplot(2,1,2); plot( vm./(1+exp(r*(v0*onerow -Ip )) )) % raw S[v]
subplot(2,1,2); plot( kIn2*vm./(1+exp(r*(v0*onerow -Ip )) )) % with const for #Links
 %axis([0 5000 0 vm+1]);
title (' kIn *S[Ip(t)] '); xlabel('t (ms) '); ylabel('pulse rate (s^-1) ');

%% 2.4 S for pulse spike train 
 % as calc byNM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
 % Ip set up in driver code (twoN, sixCluster, etc.
 r=0.56
onerow=ones(1,length(x)); % needs vector of 1's for sigmoid (over time span)
% get Ip
figure; subplot(2,1,1); plot(Ip ); title (' Ip(t), no stim '); xlabel('t (ms) ');
 ylabel('Vy (mV) '); %axis([0 5000 0 max(y(2,:))+10]);
subplot(2,1,2); plot(vm./(1+exp(r*(v0*onerow -Ip )) )) % 
 %axis([0 5000 0 vm+1]);
title (' S[Ip(t)] '); xlabel('t (ms) '); ylabel('pulse rate (s^-1) ');

%% 2.5 S for Pulse density fn 
 % as calc byNM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
 % Ip set up in driver code (twoN, sixCluster, etc.
onerow=ones(1,length(x)); % needs vector of 1's for sigmoid (over time span)
% get Ip
r=1.3
figure; subplot(2,1,1); plot(x, Pulse1 ); title (' Ip(t), pulse density fn '); xlabel('t (ms) ');
 ylabel('Vy (mV) '); %axis([0 5000 0 max(y(2,:))+10]);
subplot(2,1,2); plot(x, vm./(1+exp(r*(v0*onerow -Pulse1 )) )) % 
 %axis([0 5000 0 vm+1]);
title (' S[Ip(t)] '); xlabel('t (ms) '); ylabel('pulse rate (s^-1) ');

%% 2.5  plot S[v] curve
% x(time) & Pd set in main code
r0=0.5; % etc % v0 =6 (mV) % default
v=[-10:1: 50]; % voltage array
onerow=ones(1,length(v)); % needs vector of 1's for sigmoid (over time span)

figure;  title (' S[v] curve ');
plot(vm./(1+exp(r0*(v0*onerow - v )) ), 'k-'); hold on %
 %r0=0.45; plot(vm./(1+exp(r0*(v0*onerow - v )) ), 'r-') 
 axis([0 50 0 vm+1]);
title (' S[v] curve '); xlabel('V (mV) '); ylabel('pulse rate (s^-1) ');
  %legend('r: 0.3','  0.5', '  0.75', '  0.9')
  legend( 'r= 0.5')
  %text(60, 1.5, 'Qmax = 5 Hz');  text(60, 1.0, 'Vo = 6 mV');
 
%% 2.5  plot basic S[v] curve; series expansion
% x(time) & Pd set in main code
r0=0.5; % etc % v0 =6 (mV) % default
v=[-10:1: 50]; % voltage array
onerow=ones(1,length(v)); % needs vector of 1's for sigmoid (over time span)

figure;  title (' S[v] curve ');
plot(vm./(1+exp(r0*(v0*onerow - v )) ), 'k-'); hold on 
 %r0=0.45; plot(vm./(1+exp(r0*(v0*onerow - v )) ), 'r-') 
 axis([0 50 0 vm+1]);
title (' S[v] curve '); xlabel('V (mV) '); ylabel('pulse rate (s^-1) ');
  %legend('r: 0.3','  0.5', '  0.75', '  0.9')
  legend( 'r= 0.5')
  %text(60, 1.5, 'Qmax = 5 Hz');  text(60, 1.0, 'Vo = 6 mV');m*r
% Lin approx
s1 = 0.5*vm*(1.0-r0*(v0*onerow - v)/2 );
plot(s1, 'b--');
s3 = 0.5*vm*(1.0-r0*(v0*onerow - v)/2 - (r0*(v0*onerow - v)).^3/24 );
plot(s3, 'b');
clear s1 s3 onerow 
%% 2.5 simple plot S[v] curve; series expansion
% x(time) & Pd set in main code
r0=0.5; % etc % v0 =6 (mV) % default
v=[-10:1: 50]; % voltage array
onerow=ones(1,length(v)); % needs vector of 1's for sigmoid (over time span)

figure;  title (' S[v] curve ');
plot(1./(1+ exp(-r0*v ) ), 'k-'); hold on 
 %r0=0.45; plot(vm./(1+exp(r0*(v0*onerow - v )) ), 'r-') 
 axis([0 50 0 1.1]);
title ('simple S[v] curve '); xlabel('V (mV) '); ylabel('pulse rate (s^-1) ');
  %legend('r: 0.3','  0.5', '  0.75', '  0.9')
  legend( 'r= 0.5')
  %text(60, 1.5, 'Qmax = 5 Hz');  text(60, 1.0, 'Vo = 6 mV');m*r
% Lin approx
s1 = (1.0+ r0*v/2 )/2;
plot(s1, 'b--');
%
s3 = (1.0+ r0*v/2 + (r0*v).^3/24 )/2;
plot(s3, 'b');
s5 = (1.0+ r0*v/2 + (r0*v).^3/24 + (r0*v).^5/240 )/2;
plot(s5, 'r');
clear s1 s3 s5 onerow 

%% 2.5.1 simpler still plot S[v] curve; series expansion
% x(time) & Pd set in main code
v=[-10:1: 50]; % voltage array
onerow=ones(1,length(v)); % needs vector of 1's for sigmoid (over time span)

figure;  title (' S[v] curve ');
plot(1./(1+ exp(-v ) ), 'k-'); hold on 
 %r0=0.45; plot(vm./(1+exp(r0*(v0*onerow - v )) ), 'r-') 
 axis([0 50 0 1.1]);
title ('simple S[v] curve '); xlabel('V (mV) '); ylabel('pulse rate (s^-1) ');
  %legend('r: 0.3','  0.5', '  0.75', '  0.9')
  legend( 'r= 0.5')
  %text(60, 1.5, 'Qmax = 5 Hz');  text(60, 1.0, 'Vo = 6 mV');m*r
% Lin approx
s1 = (1.0+ v/2 )/2;
plot(s1, 'b--');
%
s3 = (1.0+ v/2 + (v).^3/24 )/2;
plot(s3, 'b');
s5 = (1.0+ v/2 - (v).^3/24 + (v).^5/240 )/2;
plot(s5, 'r');
clear s1 s3 s5 onerow 

%% 1/sqrt(N)
Nvec=[1 10 50 100 200];
figure; plot(Nvec, 1./sqrt(Nvec)); title(' 1/sqrt(N)')
figure; plot(Nvec, sqrt(Nvec)); title(' sqrt(N)')

%% 3.0 FFT Spectrum of signals (6/7/21)
 % cf. testRanSpectrum.m
  fprintf('  >> FFT of signal: peaks < 100 Hz ')
% close all; clear all
 %t0=0; tf=5; % t range (sec ?)  % these are set in sim codes
fs=1/h; % sampling frequency; Nyquist = 1/2h
%  pick the Signal: check, to load full LFP(t), at #3.0.1 below
% LFP:
   % signal = lfp(1:end); t=x(1:end); disp('  LFP (all t):'); % all 
% signal = lfp(40001:end); t=x(40001:end); disp('  LFP: 4-20s'); % after transients;  stim [settles down]
    % signal = LFP(40001:end); t=x(40001:end); disp('  LFP, from LIF sim: 4-20s'); % nb. def'n
   %signal = lfp(40001:92001); t=x(40001:92001); disp('  LFP(transition 4:9.2 s):'); % 4 - 9.2 s
   % signal = lfp(92002:end); t=x(92002:end); disp('  LFP (in s/s > 9.2 s):'); % ??after transition, s/s > 9.6s
   % signal = lfp(17002:end); t=x(97002:end); disp('  LFP (in s/s > 9.7s):');
   %signal = lfp0(40001:end); t=x(40001:end); disp('  LFP[y0]:'); % <y0 s> after transients
     %signal = lfp(120001:end); t=x(120001:end); disp('  LFP later:'); % after
     %signal = lfpB(20001:end); t=x(20001:end); disp('  LFP, Area B:'); % 
     %signal = lfpA(40001:end); t=x(40001:end); % for areas A, B
     %signal = lfpB(40001:end); t=x(40001:end); % for areas A, B
   %signal = Vit; t=tt; % Alt. eg cf testFftA... code, sample waveform
  %signal= Ve(40001:end);  t=x(40001:end); disp('LIF Ve:  excit potential [4:20s]:'); % e pot'l

% other sim outputs:
   % ye-1
   %signal= y(2, 40001:end);  t=x(40001:end);  disp('  y-e(1):'); % y-e(1): cf. ISI calc in sim code, at #4.1.2
    % signal= y(3, 40001:end);  t=x(40001:end); % y-i(1)
 % y0: pyrm
  % signal = y(1,40001:end); t=x(40001:end);  disp('  y0(1)-pyrm:'); %' NM-1-pyrm, output:
  % signal = y(19,40001:end); t=x(40001:end);  disp('  y0(4)-pyrm:'); %' NM-1-pyrm, output:
  % y3: i-fast Z('10)
  %signal = y(4,40001:end); t=x(40001:end);  disp('  y(3 - fast i):'); %' NM-1-fast-i, output:

% dy(1):
% signal = (y(2,40001:end)-y(3,40001:end)); t=x(40001:end);  disp('  dy(1):'); %' NM-1, output: dy: y1- y2 '); 
   %signal = (y(2,40001:end)-y(3,40001:end)-y(4,40001:end)); t=x(40001:end);  disp('  dy(1 -5popn):'); %' NM-1, output: dy: y1- y2 -y3 '); 
           % omit transients (eg t < 0.5sec); stim typically at 4s
% dy(2), dy(6), etc: for 6NM {3 popn} models
 signal =  (y(8,40001:end) -y(9,40001:end)); t=x(40001:end);  disp('  dy(2):') % NM#2:  dy: y1- y2 ');
   %signal =  (y(8,80001:end) -y(9,80001:end)); t=x(40001:end);  disp('  dy(2, t>8s):') % NM#2:  dy: y1- y2 ');
% signal =  (y(14,40001:end) -y(15,40001:end)); t=x(40001:end);  disp('  dy(3):') % NM#3: 
% signal =  (y(20,40001:end) -y(21,40001:end)); t=x(40001:end); disp('  dy(4): \n');% NM-4
% signal =  (y(26,40001:end) -y(27,40001:end)); t=x(40001:end); disp('  dy(5): \n'); % NM-5
% signal =  (y(32,40001:end) -y(33,40001:end)); t=x(40001:end); disp('  dy(6): \n'); % NM-6 % dy(6):
      %signal =  (y(32,1:end) -y(33,1:end)); t=x(1:end); disp('  dy(6): '); % NM-6 % dy(6); longer
 % y-e(2):
 %signal =  y(8,40001:end)'; t=x(40001:end); disp('  y-e(2): \n');% NM-2,  node-2: 1st input; y-e excit popn
  % y-e(6):
    % signal =  y(32,40001:end)';  %  node-6: 2nd input; y-e excit popn
    % signal =  (y(38,40001:end) -y(39,40001:end)); t=x(40001:end); % NM-7
    %  signal = ( y(44,(40001):end) - y(45,(40001):end) ); t=x(40001:end); % NM#8, A8b
% Others:
% signal =  (y(14,40001:end) -y(15,40001:end)); t=x(40001:end); disp('  dy-3 ')  % dy(3)
 % signal =  (y(20,40001:end)' -y(21,40001:end)'); t=x(40001:end); disp('  dy-4 ') % dy(4)
% signal =  (y(26,30000:end)' -y(27,30000:end)'); t=x(30000:end);  disp('  dy-5 ')% dy(5)
% Pulse [sequence]:
   %signal = Pulse; t=x;
% signal = Pulse -mean(Pulse); t=x; % if biased, enfore zero mean: educes dc(f=0) signal
% Multi mode waves
   %signal =V1+V2+V3; t=x;
 % DE driving forces - cf DebugCodes2.m   
 % signal = Drive; t= x(800:tend); disp(' DE nett driving force: ') % of 5-popn DEs

   % ran stim: cf testAMIRA.m code
    %signal = pran(10:end); t=tt(10:end); % ? elim transients
    %signal = pran; t=tt; % ran noise
     %signal = thetary; t=x(1:length(thetary)); % av rotn angle of dipoles
figure; plot(signal); title(' signal(t) ')
%
L=length(t); % length of signal
% FFT spectrum, raw
Y=fft(signal); Yr=real(Y)/L; Yi=imag(Y)/L;% was complex; has neg+pos freq.
%figure; plot(Yr); title('Re{ fft(signal) } raw')
% figure; stem(P2r); xlim([0 100]); hold on; stem(P2i, 'r') % debug
 % figure; stem(Yr); xlim([0 100]); hold on; stem(Yi, 'r') % correct
% figure; plot(angle(Y)*180/pi ); % phase angle of cmpl FFT components
%axis([0 20000 -Inf Inf])

% pos-side spectrum
P2r=Yr(1:round(L/2)); P2i=Yi(1:round(L/2)); 
% max freq is 1/2h : 500 Hz (for 1 ms time steps)
f=fs*(0:round(L/2-1) )/L; % [pos] freq steps are 1/Lh : nb. df varies with L
  % figure; plot(f,P2r);  xlabel(' f (Hz)'); axis([0 50 -5 5]) % for (1)
  %title('6NM, alpha, 0 Ii, 0.1mV noise, tune r & Gains: lfp spectrum ');
    %title('6NM, 20pc Ii, 1 pulse, noise, 8Hz wave: spectrum [Re] of lfp ');
    %title('6NM, 20pc Ii, 1 pulse, noise, 8Hz wave: spectrum [Re] of delta-y(1): A-10 ');
     % axis([0 50 -2 6]) % focus on low frequencies of  lfp
% get main spectral components - by hand / or use peaks(), below
  % mainf=find(abs(P2r) >1); % main peak locations (pos & neg)
  %peaks_main = f(mainf) % in Hz %peaks_medium = f(nextf) % in Hz
  % nextf=find( (abs(P2r) > 0.5) & (abs(P2r) < 1) );

%single sided FFT
P2=abs(Y/L); P1=P2(1:round((L/2)+1)); P1(2:end-1)=2.0*P1(2:end-1); % pos side; f(o) x1
figure; subplot(3,1,1); stem(P2r);  xlim([0 200]) % debug
title('fft: Real'); ylabel('fft ampl-Re ')
subplot(3,1,2); stem(P2i); title('fft Imag');  xlim([0 200])
subplot(3,1,3); stem(P1);  xlim([0 200]); title('fft: single sided'); 
 clear P2

% phase shifts & spectrum
ampl=sqrt(P2r.^2 + P2i.^2);
     %phase=atan(P2r./P2i)*180/pi; % in deg / wrong?
sigPhase= atan(P2i./P2r)*180/pi; % in deg [up from x axis
 %sigPhase1= angle(Y(1:round(L/2)))*180/pi; % of cmplx fft,  in deg
  % figure; plot(angle(Y)*180/pi ); xlim([0 100]) % phase angle 
  % hold on; plot(sigPhase1, 'r.', 'MarkerSize', 15); % ebiug: ok
%
figwidth = 560; figheight = 688; figposition = [500, 100, figwidth, figheight]; % tall
figure('position',figposition, 'units','pixels');  %figure; % default
subplot(3,1,1); plot(f,ampl); axis([0 50 0 0.5])
title('fft ampl: lfp spectrum');
ylabel('fft amplitude ')
subplot(3,1,2); plot(f,sigPhase); title('fft phase');  
axis([0 50 -100 100])
spectruma=ampl.^2; % check normalisation below? nb. half length
  %sum(spectruma) % is 41?
subplot(3,1,3); plot(f,spectruma); title('fft spectrum'); xlabel(' f (Hz)');
axis([0 50 0 0.5])

% PSD:  Spectrum only: (ampl^2); cf #3.01 tests, below
Y = Y(1:L/2+1); % Y is cmplx FFT
psdx = (1/(fs*L)) * abs(Y).^2;  % Re^2 + Im^2
psdx(2:end-1) = 2*psdx(2:end-1); % nb f(o) unique; x2 incl here
 % nb 10^3 mV / Volt
%
figure; plot(f,1e-6*psdx(1:round((L-1)/2+1))); 
 %figure; plot(f,1e-6*psdx(1:round((L-1)/2))); % alt ?
title('signal: fft spectrum:  PSD'); xlabel(' f (Hz)'); %axis([0 100 0 0.1])
ylabel('Power/Frequency  (V^{2}/Hz)'); xlim([0 100]);  % ylim([0 0.1e-4]) % to showup minor peaks
 % text(70, 1.6e-5, 'Adj = 0.5, 0.5');  text(70, 2.5e-5, '0.1 mV G ran noise, AR-1')
PSDscale = max(psdx) % debug
 %figure; plot(f,psdx(1:(L-1)/2+1)), xlim([0 150]) % extra
%
% >> Signal proc Tbox  <<
  % figure; hist(widths) % debug
   %[mainPeaks, locs, widths, proms]= findpeaks(ampl(1:1000)); % < 50, 100Hz
   %[mainPeaks, locs, widths, proms]= findpeaks(ampl); % all f peaks:  a)
   %[mainPeaks, locs, widths, proms]= findpeaks(psdx(1:5000),'MinPeakheight',5e-15); % for scale of psdx
%[mainPeaks, locs, widths, proms]= findpeaks(psdx(1:5000)); % b), out to 500Hz
%[mainPeaks, locs, widths, proms]= findpeaks(psdx); % default search for peaks
[mainPeaks, locs, widths, proms]=findpeaks(psdx(1:5000),'MinPeakheight',1e-2); % for dy-1-2
%[mainPeaks, locs, widths, proms]=findpeaks(psdx(1:5000),'MinPeakheight',1e-9); % for y-0

  %[mainPeaks, locs, widths, proms]= findpeaks(psdx(1:5000),'MinPeakheight',6500); % for noisy neurons Pulse fn
%[mainPeaks, locs, widths, proms]= findpeaks(psdx(1:5000),'MinPeakheight',0.3); % for wave
% [mainPeaks, locs, widths, proms]= findpeaks(psdx(1:5000),'MinPeakheight',6.0,'MinPeakwidth',1.0);
%[mainPeaks, locs, widths, proms]= findpeaks(psdx(1:5000), 'MinPeakDistance', 5); % min 5 Hz apart
 %
 %figure;  % Tests of other features
 %findpeaks(psdx(1:5000), 'MinPeakDistance', 5, 'Annotate', 'extents'); xlim(([0 100]))
 %figure; plot(f,ampl); hold on; axis([0 200 0 1]); % a)
figure; plot(f,psdx(1:(L-1)/2+1)); hold on; xlim([0 100]) %axis([0 200 0 1]); % b)
%ylim([0 psdx(2)]); % avoid dominace of 1st element 
 %ylim([0 5e4]); title('signal spectrum');
%title('stoch WC popn: fft spectrum:  PSD');
  %title('6NM alpha, 0.1mV noise, lfp spectrum');
  %title('2NM alpha, 1mV AR-1 noise, lfp spectrum');
%title('1+1 NM (alpha & beta), 1mV G noise: lfp spectrum');
title('2+1+2+1 NM (theta, alpha, beta, gamma), 1x G noise: lfp spectrum');
  %title('2NM alpha, 0.1mV G noise, AR-1(0.7): lfp spectrum');
   xlabel(' f (Hz)'); ylabel('PSD: Power/Freq (mV^{2}/Hz)'); 
%
fprintf('\n peaks in spectrum \n [indx  f(Hz), promince, width, phase, ampl(i), psd(i)]: ')
PeakLocs=[]; iP=0; % save in this vec, use counter
 for i=1:length(locs)  %100  % nb. sometimes there are 200? : check
 %for i=2:length(locs)  % omit peak at 0 Hz
     if proms(i) > 1e-4   % 0.5 %4e-3  %0.05  %0.1 % only plot peaks w signif prominence
        PeakLocs= [PeakLocs, locs(i)]; iP=iP+1; % save these, increment counter
         %[i PeakLocs(iP)] % debug
     plot(f(locs(i)), mainPeaks(i), 'ro') % plot peaks & F, at their locations
    %fprintf('\n  %4.0f   %4.2f   %4.2f   %4.2f   %4.2f',locs(i), f(locs(i)), proms(i), sigPhase(locs(i)), P2r(locs(i)))
     fprintf('\n  %4.0f   %4.2f   %5.2f   %4.2f   %4.2f   %5.3f    %6.4f',locs(i), f(locs(i)), proms(i), widths(i), sigPhase(locs(i)), ampl(locs(i)), psdx(locs(i))) 
     sigPhase=mod(sigPhase, 90); % appear to be the samw
     %fprintf('\n  %4.0f   %4.2f   %4.2f   %4.2f   %4.2f',locs(i), f(locs(i)), proms(i), sigPhase(locs(i)), ampl(locs(i))) 
   
     end
 end
 fprintf('\n ') % leave a new line%% old:  at peaks:
 hold off;
 Npeaks= length(PeakLocs)
   %mainPeaks = find( psdx>20); % need to check heights
  mainf= f(find( psdx>20));
  mainAmpl=P2r(find( psdx>20));
  disp('  Phase(deg) at main peaks: ')
  mainPhase =sigPhase(find( psdx>20))' % in deg
  at_f= f(find( psdx>20))
  
% Other plots:
figure; %semilogy(f, psdx(1:(L-1)/2+1))
plot(f, log10(psdx(1:(L-1)/2+1)) )
grid on; xlim([0 100]); %title('2NM, alpha,G(ran) AR-1: log10 PSD Using FFT'); 
title('6 NM, 1mV G noise: lfp spectrum'); 
xlabel('Frequency (Hz)'); ylabel('log10 Power/Frequency  (mV^{2}/Hz)')
%
figure; loglog(f, psdx); ylabel('ln PSD, V^{2}/Hz'); xlabel('ln f (Hz)');
% POWER
psdx(1)=0; % elim large DC value
disp('tot Power {f, t};  via fft') 
TotPower_signal=sum(signal.^2)*h/L % integral(signal^2 * dt) units: (mV^2/sec)
TotPower_fft=sum(psdx)/L  % here in mV^2/Hz
  % nb. value of integral ~ length of signal {can diverge at inf time
%
 %[sum(psdx(1:801))/(L*TotPower_fft)  sum(psdx(802:2001))/(L*TotPower_fft) ] % for h 0.1ms & T:20s
 %[sum(psdx(1:801))/(L*TotPower_fft)*h/800  sum(psdx(802:2001))/(L*TotPower_fft)*h/(1199) ] % normalised
 %frn_gam=sum(psdx(802:1601))/(L*TotPower_fft); % 40+hz for [0:20]sec
 frn_gam=sum(psdx(482:1601))/(L*TotPower_fft); % 30+hz  [0:20]sec
 %disp('Power fraction in theta band: [4:8]Hz:')
 frn_thet=sum(psdx(64:128))/(L*TotPower_fft);  % check f(i) ok for [0:20 sec sim]
 %disp('Power fraction in alpha band: [8:12]Hz:')
 frn_alph=sum(psdx(129:193))/(L*TotPower_fft); 
 %disp('Power fraction in beta band: [12:30]Hz:')
 frn_bet=sum(psdx(193:481))/(L*TotPower_fft);  % check f(i) ok
  %frn_sharpRip=sum(psdx(1602:2401))/(L*TotPower_fft) % 100:150 Hz { v small
 fprintf('\n theta   alpha   beta  gamma   fractions of power: ')
 fprintf('\n  %4.3f   %4.2f  %4.3f  %6.4f \n',frn_thet, frn_alph, frn_bet, frn_gam)
%
  %powest=1e-6*powest;
  %[powest, idx] = max(psdx); % max peak power estm
 
% PSD windowed, 
% calculate the noise windowed PSD
winlen = 2*fs;
window = hanning(winlen, 'periodic'); % DSP hann window
noverlap = winlen/2; nfft = winlen;
%
[Pxx, fx] = pwelch(signal, window, noverlap, nfft, fs, 'onesided');
figure; plot(fx, Pxx);  xlim([0 300]);  title('windowed PSD') % drops rapidly > 10 Hz
 
PxxdB = 10*log10(Pxx); % dB
figure; semilogx(fx, PxxdB, 'r', 'LineWidth', 1.5)
%f=f';  %psdx=psdx';  % col vec easier to navigate in GUI

%
clear zeroCross posSlope tstart sample PhaseEstm testM VatZero firstZero widths locs mainPeaks at_f
clear Y* peaks* mainPeaks P2r P2i P1 spectruma t signal ampl sigPhase* Peak*  frn_ %psdx zL
clear  main* widths proms iP  locs Npeaks fx Pxx* win* nover* nfft powest PSDscale frn_* % TotPower*

%% 3.01  cf. mathworks code: cf testARIMA.m 
   %signal = pran;
signal = lfp(40001:end); t=x(40001:end); % for NM, after transients;
N = length(signal);  % same as L, above
fs=1/h; % sampling frequency (10k Hz); Nyquist = 1/2h (5k Hz)
h = 0.0001;  % dt step (0.1 or 0.5ms for d's;  was 1 ms)
xdft = fft(signal); % complex
 %figure; plot(xdft); 
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;  % Re^2 + Im^2
psdx(2:end-1) = 2*psdx(2:end-1); % nb f(o) unique; x2 incl here
freq = 0:fs/N:fs/2;  % nb. same as f, above

figure; plot(freq, 1e-6*psdx(1:(N-1)/2+1))
grid on; title('G(ran) AR-1:  PSD Using FFT'); xlim([0 200])
xlabel('Frequency (Hz)'); ylabel('Power/Frequency  (V^{2}/Hz)')
[powest, idx] = max(psdx) % max power estm
freq(idx)
 %MaxPwrAtFreqhz =  [ powest  freq(idx) ]
%
figure; semilogy(freq, psdx(1:(N-1)/2+1))
grid on; title('G(ran) AR-1: log PSD Using FFT'); xlim([0 200])
xlabel('Frequency (Hz)'); ylabel('Power/Frequency  (mV^{2}/Hz)')
 sum(psdx)/N % Tot power
 TotPower_signal=sum(signal.^2)*h/N % integral(signal^2 * dt) units: (V^2/sec)

 clear xdft psdx freq MaxPwrAtFreqhz N idx powest 

%% 3.0.1 Use [full 1:T]  LFP
fprintf('\n > use LFP for full t \n')
% net lfp output: sum: detla-y(1) + detla-y(2), for the 8 NM - incl transients
 lfp= (y(2,:) - y(3,:) + y(8,:) - y(9,:) +y(14,:) - y(15,:) ...
     + y(20,:) - y(21,:)+ y(26,:) - y(27,:) +y(32,:) - y(33,:) ...
     + y(38,:) - y(39,:) + y(44,:) - y(45,:) )/8; % av of 8 nodes
  %figure; plot(lfp); title('JR model: net lfp output ')
  figure; plot(x, lfp); title('8-star cluster, JR model: net lfp(t) output '); xlabel('t (sec) ')

%% 3..0.2 Reassemble  waveform as cos(wt + phi)
 % nb. FFT coeffs: a = ampl*cos(phi), b= ampl*sin(phi)
 % x vec: time [0.1, 0.5 ms steps]  set in sixCluster code, etc
 % t is sampled part: eg skip transients
  %sigPhase= - sigPhase;  % lag ?
    %x1 =x((tstart+10000):end); % [2:end] - avoids transients
SigAprox=zeros(1,length(t))+ampl(1)/2; % set up row vec i-th cpnt of waveform; DC term
onerow=ones(1,length(t));
for i = 1:length(PeakLocs)
 iloc=PeakLocs(i); % 1st peak, from #3.0
 SigAprox = SigAprox + ampl(iloc)*cos( 2*pi*f(iloc)*t + (pi/180)*sigPhase(iloc)*onerow ); % cos(wt+phi [in rad])
  %SigAprox = SigAprox + ampl(iloc)*sin( 2*pi*f(iloc)*t + (pi/180)*sigPhase(iloc)*onerow ); % cos(wt+phi [in rad])
end

figure; 
subplot(2,1,1); plot(t,signal);  xlim([1.5 3.5]); grid on; hold on
plot([2.0 2.5], [mean(signal) mean(signal)], 'r')  % cf zero crossing, at mean
title('orignal LFP(t)'); %ylabel('fft amplitude ')
subplot(2,1,2); plot(t,SigAprox(:), 'k-');  xlim([1.5 3.5]); grid on
hold on; plot([2.0 2.5], [mean(SigAprox) mean(SigAprox)], 'r') 
%axis([ 1 2 -Inf Inf]); 
xlabel(' t (sec) '); title('LFP: fft components; 1st + harmonics')
MeanSignals =[mean(signal) mean(SigAprox)]

clear SigAprox MeanSignals onerow iloc % SigAprox
  
%% 3.1 Reassemble (progressively plot) the waveform
 % x vec: time [0.5/0.1 ms steps]  set in sixCluster code
 %phase=-phase;
Vit=zeros(1,length(x)); % set up row vec i-th cpnt of waveform
onerow=ones(1,length(x));
%
figure; subplot(2,1,1); hold on; title('LFP: progressive sum of fft cpnts'); xlabel('t (sec) ')
for i=1:Npeaks
    Vit(:) = Vit(:)+ P2r(PeakLocs(i))*cos( 2*pi*f(PeakLocs(i))*x + (pi/180)*phase(PeakLocs(i)) *onerow )'; % cost(wt+phi [in rad])
    plot(x, Vit(:) );
 % partial sum
  drawnow
  pause
end
  %axis([ 1 4 -Inf Inf]); xlabel(' t (ms) '); % expanded time scale
% & cf orig lfp(t)
subplot(2,1,2); plot(x, lfp);title('8-NM, alpha; orig LFP(t) '); xlabel('t (sec) ')
 
%% 3.2 Reassemble  real cpnts of the waveform
 % x vec: time [0.5 ms steps]  set in sixCluster code
 phase=-phase;
 %x1 =x((tstart+10000):end); % [2:end] - avoids transients
 x1=x;
Vit=zeros(6,length(x1)); % set up row vec i-th cpnt of waveform
onerow=ones(1,length(x1));
iloc=90; % 1st peak, from #3.0
Vit(1,:) = P2r(iloc)*cos( 2*pi*f(iloc)*x1 + (pi/180)*phase(iloc)*onerow ); % cost(wt+phi [in rad])
figure; plot(x1,Vit(1,:), 'k-')
axis([ 1 2 -Inf Inf]); xlabel(' t (ms) '); title('lfp: 5 fft components; 1st harm')
hold on
iloc=78; % 2nd peak, from #3.0
Vit(2,:) = P2r(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(2,:), 'b-')
iloc=71; % 3rd peak, from #3.0
Vit(3,:) = P2r(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(3,:), 'r-')
iloc=103; % 4th peak, from #3.0
Vit(4,:) = P2r(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(4,:), '--', 'Color', [0.5 0.5 0.5])
 iloc=66; % 5th peak, from #3.0 / skip others ...
Vit(5,:) = P2r(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(5,:), 'k:')
 iloc=210; % 1st harmonic
Vit(6,:) = P2r(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(6,:), 'k--')
hold off

% partial sum
figure; subplot(2,1,1); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:) );
title('lfp: sum 3 fft cpnts')
subplot(2,1,2); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:)+Vit(4,:)+Vit(5,:)+Vit(6,:) ); hold on
axis([ 2 4 -Inf Inf]); xlabel(' t (ms) '); 
% & cf orig lfp(t)
figure; subplot(2,1,1);  plot(x1, lfp);title('8-NM, alpha; orig lfp(t) '); xlabel('t (sec) ')
subplot(2,1,2); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:) );
title('lfp: sum 3 fft cpnts'); xlabel('t (sec) ')
%

% partial sum
figure; subplot(2,1,1); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:)+Vit(4,:)+Vit(5,:)+Vit(6,:) );
title('lfp: sum 6 fft cpnts ')
subplot(2,1,2); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:)+Vit(4,:)+Vit(5,:)+Vit(6,:) ); hold on
axis([ 2 4 -Inf Inf]); xlabel(' t (ms) '); % expanded time scale
% & cf orig lfp(t)
figure; subplot(2,1,1);  plot(x1, lfp);title('8-NM, alpha; orig lfp(t) '); xlabel('t (sec) ')
subplot(2,1,2); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:)+Vit(4,:)+Vit(5,:)+Vit(6,:) );
title('lfp: sum 6 fft cpnts'); xlabel('t (sec) ')

%% 3.2a Reassemble  real cpnts of the waveform as sin(wt + phi)
 % x vec: time [0.5 ms steps]  set in sixCluster code
 phase=-phase;
 %x1 =x((tstart+10000):end); % [2:end] - avoids transients
 x1=x;
Vit=zeros(6,length(x1)); % set up row vec i-th cpnt of waveform
onerow=ones(1,length(x1));
iloc=60; % 1st peak, from #3.0
Vit(1,:) = P2r(1)/2 + P2r(iloc)*sin( 2*pi*f(iloc)*x1 + (pi/180)*phase(iloc)*onerow ); % cost(wt+phi [in rad])
figure; plot(x1,Vit(1,:), 'k-')
%axis([ 1 2 -Inf Inf]); 
xlabel(' t (ms) '); title('lfp: 5 fft components; 1st harm')
hold on

%
%iloc=78; % 2nd peak, from #3.0
Vit(2,:) = P2r(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(2,:), 'b-')
iloc=71; % 3rd peak, from #3.0
Vit(3,:) = P2r(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(3,:), 'r-')
iloc=103; % 4th peak, from #3.0
Vit(4,:) = P2r(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(4,:), '--', 'Color', [0.5 0.5 0.5])
 iloc=66; % 5th peak, from #3.0 / skip others ...
Vit(5,:) = P2r(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(5,:), 'k:')
 iloc=210; % 1st harmonic
Vit(6,:) = P2r(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(6,:), 'k--')
hold off

% partial sum
figure; subplot(2,1,1); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:) );
title('lfp: sum 3 fft cpnts')
subplot(2,1,2); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:)+Vit(4,:)+Vit(5,:)+Vit(6,:) ); hold on
axis([ 2 4 -Inf Inf]); xlabel(' t (ms) '); 
% & cf orig lfp(t)
figure; subplot(2,1,1);  plot(x1, lfp);title('8-NM, alpha; orig lfp(t) '); xlabel('t (sec) ')
subplot(2,1,2); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:) );
title('lfp: sum 3 fft cpnts'); xlabel('t (sec) ')
%

% partial sum
figure; subplot(2,1,1); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:)+Vit(4,:)+Vit(5,:)+Vit(6,:) );
title('lfp: sum 6 fft cpnts ')
subplot(2,1,2); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:)+Vit(4,:)+Vit(5,:)+Vit(6,:) ); hold on
axis([ 2 4 -Inf Inf]); xlabel(' t (ms) '); % expanded time scale
% & cf orig lfp(t)
figure; subplot(2,1,1);  plot(x1, lfp);title('8-NM, alpha; orig lfp(t) '); xlabel('t (sec) ')
subplot(2,1,2); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:)+Vit(4,:)+Vit(5,:)+Vit(6,:) );
title('lfp: sum 6 fft cpnts'); xlabel('t (sec) ')

%% 3.2b Reassemble amplitudes of the waveform
 % x vec: time [0.5 ms steps]  set in sixCluster code
 %phase=-phase;
 %x1 =x((tstart+10000):end); % [2:end] - avoids transients
 x1=x;
Vit=zeros(6,length(x1)); % set up row vec i-th cpnt of waveform
onerow=ones(1,length(x1));
iloc=90; % 1st peak, from #3.0
Vit(1,:) = ampl(iloc)*cos( 2*pi*f(iloc)*x1 + (pi/180)*phase(iloc)*onerow ); % cost(wt+phi [in rad])
figure; plot(x1,Vit(1,:), 'k-')
axis([ 1 2 -Inf Inf]); xlabel(' t (ms) '); title('lfp: 5 fft components; 1st harm')
hold on
iloc=78; % 2nd peak, from #3.0
Vit(2,:) = ampl(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(2,:), 'b-')
iloc=71; % 3rd peak, from #3.0
Vit(3,:) = ampl(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(3,:), 'r-')
iloc=103; % 4th peak, from #3.0
Vit(4,:) = ampl(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(4,:), '--', 'Color', [0.5 0.5 0.5])
 iloc=66; % 5th peak, from #3.0 / skip others ...
Vit(5,:) = ampl(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(5,:), 'k:')
 iloc=210; % 1st harmonic
Vit(6,:) = ampl(iloc)*cos( 2*pi*f(iloc)*x + (pi/180)*phase(iloc)*onerow );
 plot(x1,Vit(6,:), 'k--')
hold off

% partial sum
figure; subplot(2,1,1); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:) );
title('lfp: sum 3 fft cpnts')
subplot(2,1,2); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:)+Vit(4,:)+Vit(5,:)+Vit(6,:) ); hold on
axis([ 2 4 -Inf Inf]); xlabel(' t (ms) '); 
% & cf orig lfp(t)
figure; subplot(2,1,1);  plot(x1, lfp);title('8-NM, alpha; orig lfp(t) '); xlabel('t (sec) ')
subplot(2,1,2); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:) );
title('lfp: sum 3 fft cpnts'); xlabel('t (sec) ')
%

% partial sum
figure; subplot(2,1,1); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:)+Vit(4,:)+Vit(5,:)+Vit(6,:) );
title('lfp: sum 6 fft cpnts ')
subplot(2,1,2); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:)+Vit(4,:)+Vit(5,:)+Vit(6,:) ); hold on
axis([ 2 4 -Inf Inf]); xlabel(' t (ms) '); % expanded time scale
% & cf orig lfp(t)
figure; subplot(2,1,1);  plot(x1, lfp);title('8-NM, alpha; orig lfp(t) '); xlabel('t (sec) ')
subplot(2,1,2); plot(x1, Vit(1,:)+Vit(2,:)+Vit(3,:)+Vit(4,:)+Vit(5,:)+Vit(6,:) );
title('lfp: sum 6 fft cpnts'); xlabel('t (sec) ')

%% Dev. 3.0a phase delay wrt t= 2.0s : locate zero-crossings; plot to check
tstart=20000; % t step to being sample (@ 0.1ms, 0.5 etc)
 %sample=lfp(tstart:(tstart+2000)); sample=(sample - mean(sample));
tstart=1; % for Vit
sample=Vit(tstart:(tstart+5000-1)); sample=(sample - mean(sample));
figure; plot(sample); hold on; title('signal sample'); xlabel('t (0.1ms)')
VatZero = sample(1)
firstZero=  find(abs(sample) <= 0.02, 1) % find 1st zero; as check
zeroCross=find(abs(sample) <= 0.02)
 plot(zeroCross, sample(zeroCross), 'o'); grid on
 fprintf('\n ZeroCross at: %5.4f  %5.4f %5.4f (sec)',h*zeroCross(1:3))
 fprintf('\n       Signal: %5.3f  %5.3f %5.3f \n',sample(zeroCross(1:3))) 
%posSlope = (lfp(zeroCross(2))- lfp(zeroCross(1)) >0) % check inscreasing signal
posSlope = (Vit(zeroCross(2))- Vit(zeroCross(1)) >0) % check inscreasing signal
PhaseEstm= h*zeroCross(1)*f(locs(1))*360   % ie. (delta-t/T)*360 deg

% test metric, at 1st zero crossing ~ phase lag
testM=0;
 % ampl(1)/2  % ?? needed? / use a(i)=P2r or [ampl & phase]
for i = 1:length(PeakLocs)
 iloc=PeakLocs(i); % 1st peak, from #3.0
 %testM = testM + ampl(iloc)*cos( 2*pi*f(iloc)*zeroCross(1) - (pi/180)*sigPhase(iloc)); % cos(wt+phi [in rad])
  testM = testM + P2r(iloc)*cos( 2*pi*f(iloc)*zeroCross(1) - (pi/180)*sigPhase(iloc)); % cos(wt+phi [in rad])
end
testM

%% Dev 3.0b work on the peaks of wave(t)
 % sample set in #3.0a
 % examine plot, to "guess" heights & separations
 % these are peaks of V(t)
[SamplePeaks, PeakLocs, widths, proms]= findpeaks(sample, ...
    'MinPeakDistance',500 ,'MinPeakheight',0.5); %
length(PeakLocs)
PeakTimes=PeakLocs*h  %(1:5)
 % t(locs) % 24 peaks for 5 sec(@0.1)
shortT = mean(diff(PeakLocs))*h  % short period
shortT_freq= [shortT 1/shortT]

figure; plot(sample); hold on; title('signal sample'); xlabel('t (0.1ms)')
plot(PeakLocs, sample(PeakLocs), 'o'); grid on
for i=1:length(PeakLocs)
 plot([PeakLocs(i) PeakLocs(i)] , [min(sample) max(sample)], 'r');
end
 
clear PeakTimes shortT*

%% 3.3 test waveform
 x1=x; % time domain; set in 8NM code
Vit=zeros(1,length(x)); % set up row vec i-th cpnt of waveform
onerow=ones(1,length(x));
dphi=15; % phase increment(deg)
fn= [ 7.4 8.5 9.1 14.9 16.3];   % freq (Hz)
coeff =[0.8 0.8 0.45 0.08 0.08];  % ampl for f
for i =1:length(fn)
Vit(:) = Vit(:)+ coeff(i)*cos( 2*pi*(i*fn(i))*x + (pi/180)*(i*dphi)*onerow )'; % (1/i)*cost(wt+phi [in rad])
end

figure; plot(x, Vit)
 % axis([0 2 -4 4])
 % lfp = Vit;  for fft & xcorr codes
 
%% 3.1 play the sound of NM waveform
 % assumes data in [-1 1]: use soundsc to rescale
 signal=lfp(tstart:end); %skip tranient
 signal= signal-mean(signal); % zero mean
 level=max(max(signal), abs(min(signal)) ); signal= signal/level;
 for i=1:50
    %soundsc(signal, fs)  % at sampling freq. (in Hz: here 1 kHz)
    sound(signal, fs);
 end
   % gets echo ... later ??
 %player=audioplayer(lfp, fs);% ? fails
 %play(audioplayer)
 
 clear signal level 
 
%% Appx 1: test sounds in matlab - load eg. chirp
load chirp % birds y(t): 13,129 pts; Fs 8192 (x 1.6)
signal=y; 
 sound(signal, Fs)  % play, at Fs - for 2.1 sec
 
 % player=audioplayer(signal, Fs); % ? fails
 % play(audioplayer)
% time domain
L=length(signal);
h = 1/Fs   % t step (1 ms) /sampling frequency
T=L*h
t = [0:h:T-h];  % t domain  %
figure; plot(t, signal); xlabel('t (sec)'); title('eg. bird chirp: waveform')

% FFT spectrum, raw
Y=fft(signal); Yr=real(Y)/L; Yi=imag(Y)/L;% was complex; has neg+pos freq. 
% pos-side spectrum
P2r=Yr(1:round(L/2)); % P2i=Yi(1:round(L/2)); 
% max freq is 1/2h : 500 Hz (for 1 ms time steps)
f=Fs*(0:(L/2 ))/L; % [pos] freq steps are 1/Lh : 0.2 Hz 
figure; plot(f, P2r);  xlabel(' f (Hz)');  title('eg. bird chirp: spectrum [Re] of lfp');
title('eg. bird chirp: spectrum')

clear y signal Fs f L h T t Y Yr Yi P2r
 
%% 3.3 Autocorrelation; find periods
% dt: h set in 8NM & 6NM codes: 1 or 0.1 ms
fprintf('\n AutoCorrelation fn (from 4s), Periods of:  LFP, or y-e(1), dy(1,2) \n')
fs=1/h; % sampling frequency
%t = x(1:end-1);  % t domain - default - cf x vector in NM codes

%  pick the signal ** needs zero mean
% signal = lfp(40000:end) -mean(lfp(40000:end)); % LFP:  same length as t
% signal = Vyt(40000:end); % Wave form;  alt.
    % signal = Pulse1(1:end-1); % Pulse fn;  x1, x10
% signal= y(2, 40000:end)-mean(y(2, 40000:end));  t=x(40000:end)-4.0; % y-e(1): cf. ISI calc
%signal= (y(2, 40000:end)-y(3, 40000:end)) -mean(y(2, 40000:end)-y(3,40000:end)); % dy(1)
signal= (y(8, 40000:end)-y(9, 40000:end)) -mean(y(8, 40000:end)-y(9,40000:end));  disp('\n dy-2:')  % dy(2)
t=x(40000:end)-4.0; %t=x(40000:end)-4.0; % 
L=length(t); % length of signal
fncorr=xcorr(signal);
  % nb double sided; need to flip for -t:
  tneg= -fliplr(t); dblt=[tneg(1:end), t(2:end)]; % nb count t=0 only once
figure; plot(dblt, fncorr); %title('8NM: dy(2): AutoCorrlnFn(t); Peaks ') % vs t in (sec) : need same length
title('6 NM: LFP: AutoCorrlnFn(t); Peaks ')
hold on; xlabel('Lag t (sec)'); %axis([0 10 -Inf Inf]) % pos time
% xlim([0 7]);  %xlim([-2 2])

% find peaks: read plot to get scale
middle= (length(fncorr)+1)/2;
[CorrPeaks, locs, widths, proms]= findpeaks(fncorr(middle:end)); % use [0:tmax]
% locate each type of peak
%figure;
xlim([0 Inf])
plot(t(locs), CorrPeaks,'.r', 'MarkerSize', 18); hold on % nb t shift [sample > 3 sec]
 
% dblt(locs) % 85 peaks!
short = mean(diff(locs))/fs;  % short period
shortT_freq= [short 1/short]
 % next find longer period[s]
 % LFP: 1.25e5;  Pulses 1e7;  waves 0.1
 [MedPeaks, MedLocs] =findpeaks(fncorr(middle:end), ...
    'MinPeakDistance',ceil(short)*fs,'MinPeakheight',1e6); % fast/ medium peaks 1M & 10M 
 medium = mean(diff(MedLocs))/fs; % medium period
 mediumT_freq= [medium 1/medium]
 % next find longer period[s]
 [MainPeaks, MainLocs] =findpeaks(fncorr(middle:end), ...
    'MinPeakDistance',ceil(medium)*fs,'MinPeakheight',2e5); % slow peaks % 500k & 1M
 long = mean(diff(MainLocs))/fs; % long period

Periods = t(MainLocs) % (sec)
freqs = 1./Periods % (Hz)

 % locate medium peaks
plot(t(MedLocs), MedPeaks +0.0,'.b', 'MarkerSize', 16);
plot(t(MedLocs), MedPeaks +0.0,'sb', 'MarkerSize', 14);
plot(t(MainLocs), MainPeaks +0.0,'sk', 'MarkerSize', 14);
%legend(' ', 'fast period', ' ', 'medium period', 'slow period'); 
% title('8NM: dy(1): Peaks of AutoCorrlnFn ') 
 %text(5, -4e5, 'cut at 1.25 e5') % etc

% pos side signal & auto-correlation:
 middle= (length(fncorr)+1)/2;
 figure;  subplot(2,1,1); plot(t,signal); title('signal (sample)'); hold on
  xlabel('t (sec + 4) '); xlim([0 30]) %xlim([4 5]);
 subplot(2,1,2); plot(t, fncorr(middle:end)/max(fncorr(middle:end)), 'k-'); xlabel('t (sec) ')
  grid on; xlim([0 30])  % xlim([0 1]);
  title(' Auto-correlation fn.');
 
figure; semilogy(Periods)
 clear fncorr tneg middle short* middle medium* long* MainPeaks MainLocs MedLocs Periods
 clear signal t CoorPeaks Medpeaks freqs 
 
%% Appx 2.0  test erf, erfc etc - Re: S[v]
% raw fn
%tline=[0.1:0.1:5.0];
 tline=[-2.0:0.1:10.0]; % incl neg argm
tmp1=erfc(1./tline);
figure; plot(tline, tmp1, 'k-')
xlabel('t'); ylabel l('erfc(1/t)'); grid on
  hold on; plot(tline, erf(tline), 'b--' )
legend('erfc(1/t)', 'erf(t)')
%
figure; plot(tline, erf(tline), 'k-' ); hold on
plot(tline, erf(1./tline), 'b--' );
 %legend('erf(t)', 'erf(1/t)')
%figure; 
plot(tline, 1-erf(tline), 'k--' ); %hold on
plot(tline, 1-erf(1./tline), 'r-' );
legend('erf(t)', 'erf(1/t)',  '1-erfc(t)', '1-erf(1/t)')

clear tline tmp1
 
%% 2.0.1 test erf(x) - erf(a)
% raw fn
%tline=[0.1:0.1:5.0];
 tline= [-2.0:0.1:20.0]; % incl neg argm
 onerow=ones(1,length(tline)); % needs vector of 1's for sigmoid 
 a = 2;  r0=0.3; v0=6; % for a =r*V0
figure; plot(tline, erf(tline), 'k-' ); hold on
tmp1= 1+(erf(r0*tline-a) - erf(a) )/2;
plot(tline, tmp1, 'k--' ); %legend('erf(x)', '1+[erf(a*x)-erf(a)]/2');

plot(tline, 1.0./(1+exp(r0*(v0*onerow - tline )) ), 'b--'); % S[v(t)]
legend('erf(x)', '1+[erf(a*x)-erf(a)]/2', 'S[v]');
text(15, 0.2, 'a = 2; r =0.3'); plot([-2 20], [0.5 0.5], 'b:') % 50% activation
 title('Sigmoid and erf responses'); xlabel('V input (mV) '); 
 ylabel('pulse rate output (s^-1) '); ylim([0 1])
 
clear tline onerow r0 v0 a tmp1

%%  2.1 plot S[v] sigmoid curve - vs erf form : from #2.5, above
r0=0.4; v0=6; % defaults for JR; etc
v=[-10:1: 50]; % voltage array
onerow=ones(1,length(v)); % needs vector of 1's for sigmoid (over time span)
figure;  plot(v, 1.0./(1+exp(r0*(v0*onerow - v )) ), 'k-');  hold on;  % 
 xlabel('V input (mV) '); ylabel('pulse rate output (s^-1) ');
 title(' basic S[v] curve ');
  %legend('r: 0.3') %legend('r: 0.3','  0.5', '  0.75', '  0.9')
  text(20, 0.1, 'r = 0.4; Vo = 6 mV'); xlim([-10, 40]); ylim([-0.05, 1.0])
  %legend('r= 0.3', 'r= 0.4','r= 0.5','r= 0.56')
  %text(60, 1.5, 'Qmax = 5 Hs');    
% Appx 2.1.1 vs. erf form
   % erf(r0*v0) % const offset
 %r0=0.3;
 %r1=r0/(sqrt(2));
  S1= (erf(r0*(v -v0*onerow)) +erf(r0*v0) )/2;
 plot(v, S1, 'k--' );  % xlim([-10, 40]); ylim([-0.05 1]);
  legend(' S[v]', 'erf(v)');  title('Sigmoid and erf responses')
% new r
r0=0.3;
plot(v, 1.0./(1+exp(r0*(v0*onerow - v )) ), 'r-')
S1= (erf(r0*(v -v0*onerow)) +erf(r0*v0) )/2;
plot(v, S1, 'r--' ); 
text(20, 0.05, 'r = 0.3; Vo = 6 mV', 'Color', 'r'); 

clear onerow S1 r0

%% Appx 2.1.2 Sigmoid vs. erfc(1/v) wt form
   % erf(r0*v0) % const offset
vpos=[-2:1: 30];
 %vpos=[-5:0.1: 5];
 %vpos=[-5.1:1: 50]; % avoids divergence at 0 
r0=0.5; v0= 6;  % now fixed threshold; r=1/(sqrt(2)*sigma) [of wts]
% vs. wt/ erfc form
   % erf(r0*v0) % const offset
onerowp=ones(1,length(vpos));
 % r1=1.0; % new guess?
% hi wt
 wbar= 3.23; % {max} mean wt/link, for this pop'n
S2= 1 -erfc(r0*(wbar*onerowp -v0./vpos))/2;
 maxS_hiWt = S2(end)
 figure; hold on; plot(vpos, S2, 'r--' ); %ylim([0 1]);  %xlim([-10, 40]); %ylim([0 1]);
% unit
 wbar= 1.0; % {max} mean wt/link, for this pop'n
S2= 1 -erfc(r0*(wbar*onerowp -v0./vpos))/2;
 maxS_unitWt = S2(end)
 plot(vpos, S2, 'g-' );
% lo wt 
wbar=0.48; % {min} mean wt/link
S2= 1 -erfc(r0*(wbar*onerowp -v0./vpos))/2;  %erfc(r0*(v0./vpos -wbar*onerowp))/2;
  maxS_loWt = S2(end)
  plot(vpos, S2, 'b--' );
% cf. sigmoid  
plot(vpos, 1.0./(1+exp(r0*(v0*onerowp - vpos )) ), 'k-'); % basic sigmoid S[v] form
  title('S[v] via synaptic wt distribution; erfc (1/V)')
  legend('wbar = 3.23', 'wbar = 1.0', 'wbar = 0.48', 'S[v], sigmoid') 
 text(30, 0.1, 'r = 0.5; Vo = 6 mV');  
  plot([0 50], [0.5 0.5], 'b:') % 50% activation

  clear r0 r1 onerow onerowp v vpos S1 S2 maxS*

%% Appx 2.1.2a  mean & std for. erfc(1/v) wt form
% cf. Mouse Retina data for wt-Deg - fit to LogNorm distr'n 
Wtmean = exp(6.4 + 0.8^2/2)
WtSTD = exp(0.8^2 -1)*exp(2*6.4 + 0.8^2)
  clear Wtmean WtSTD 
  
 %% Appx 2.1.2b Test alt form of S[v]~ erfc(1/v): - use 1-erf
   % erf(r0*v0) % const offset (cf. integral)
%vpos=[0.1:1: 50];
vpos=[-2:0.1: 20];
 %vpos=[-20:1: 50];
r0=0.4; v0= 6;  % now fixed threshold; r=1/(sqrt(2)*sigma) [of wts]
% vs. wt/ erfc form
   % erf(r0*v0) % const offset
onerowp=ones(1,length(vpos));
 wbar= 1; %3.23; % {max} mean wt/link, for this pop'n
 %S2= 1 -erfc(r0*(wbar*onerowp -v0./vpos))/2;
 argm=r0*(wbar*onerowp -v0./vpos);
 figure; 
 subplot(3,1,1); plot(vpos, argm, 'k' ); legend('argument')
 title('S[V - for N(wt) via erf]')
 
 subplot(3,1,2);  % plot(vpos, erf(argm), 'b' ); legend('S[V - erf form]')
S2= (1 + erf(r0*(wbar*onerowp -v0./vpos)) )/2; % alt form
 plot(vpos, S2, 'b' ); hold on
 plot(vpos, 1.0./(1+exp(r0*(v0*onerowp - vpos )) ), 'k-'); % basic sigmoid S[v] form
 legend('S[V - erf form]', 'Sigmoid'); text(6, 0.2, 'r = 0.4; wbar = 1')
  %maxS = S2(end)
  %figure; hold on; 
   
 subplot(3,1,3); 
    %plot(vpos, S2, 'r--' );  hold on; %ylim([0 1.2]); %xlim([-10, 40]);   
 % other limit of params  
wbar= 0.48; % {max} mean wt/link, for this pop'n
S2= (1 + erf(r0*(wbar*onerowp -v0./vpos)) )/2; % erf form
maxS_erf = S2(end)
 plot(vpos, S2, 'b--' ); hold on
S2= 1 -erfc(r0*(wbar*onerowp -v0./vpos))/2;  % erfc form
 maxS_erfc = S2(end)
 plot(vpos, S2, 'g--' ); 
 plot(vpos, 1.0./(1+exp(r0*(v0*onerowp - vpos )) ), 'k-'); % basic sigmoid S[v] form
 legend('S[v] - erf form', '   erfc form', 'Sigmoid')
 text(6, 0.2, 'r = 0.4; wbar = 0.48')
 clear r0 r1 onerow onerowp v vpos S1 S2 maxS*
 
  %% Appx 2.1.2c vs. erfc(1/v) wt. Test erf form for code
   % erfc(r0*v0/v(t)) % for N(wt) form
   % nb. neg argm seems ok for erfc
vpos=[0.01:0.5:50];
%vpos=[-2:1: 50]; % test neg v - nb. unusual behav'r of erfc(-)
r0=0.6; v0= 6;  % now fixed threshold; r=1/(sqrt(2)*sigma) [of wts]
wbar= 0.48; % {max} mean wt/link, for this pop'n
%wbar=0.48; % {min} mean wt/link
%wbar= 6.73; % if use > 1
wbar=1.0;
 %wbar=0.15; % min in alt form
onerowp=ones(1,length(vpos));
argm=r0*(wbar*onerowp -v0./vpos); % argument of erfc form of S[w, v]
     %stepfn= zeros(1,length(vpos)); stepfn(find(argm>0))=1.0; % build heaviside step fn
      % nb. that form is "1" for neg argm: only working at "0"
     %argm=stepfn.*argm; % apply "heaviside step fn"; exclude neg argm of erfc
     %[i dyi argm]  % debug
%  OR, just use erfc as is 
figure; subplot(3,1,1)
 plot(vpos, argm, 'k' ); legend('argument'); title('test S[w; v] erfc form - all argm')
 subplot(3,1,2);  plot(vpos, erf(argm), 'b' ); legend('erf(argm)')
 %S2= 1 -erfc(argm)/2; % w argm filtered 
  %S2= 1 -erfc(r0*(wbar*onerowp -v0./vpos))/2; % orig form
  %S2= (1 + erf(r0*(wbar*onerowp -v0./vpos)) )/2; % alt form
 S2= (1 + erf(argm) )/2; % alt w argm filtered
  %maxS = S2(end) 
 subplot(3,1,3); plot(vpos, S2, 'r--' ); %ylim([0 1]);  %xlim([-10, 40]); %ylim([0 1]);
 hold on;  xlabel('V(t)  (mV)'); text(30, 0.15, 'r = 0.6; vo = 6 mV;  wbar = 1.0 ')
 plot(vpos, 1.0./(1+exp(r0*(v0*onerowp - vpos)) ), 'k-'); % basic sigmoid S[v] form
 %plot(vpos, 1.0./(1+exp(r0*(v0*onerowp - wbar*vpos )) ), 'k-'); % sigm form,modified by wbar
 legend('S[v]', 'Sigmoid'); ylim([0 1]);
 
 % clear argm r0 r1 onerow onerowp  vpos S1 S2 maxS stepfn
 
 %% Appx 2.1.2d. test arg & components parts - for erfc(1/V)
 vpos=[0.1:1: 50]; r0=0.5; v0=6;
 figure; plot(vpos, v0./vpos, 'k--'); hold on; title(' 1/V(t) form')  
 wbar=0.48; % {min} mean wt/link
  plot(vpos, (v0./vpos -wbar), 'b-');
 wbar= 3.23; % {max} mean wt/link, for this pop'n
  plot(vpos, (v0./vpos -wbar), 'r-');
 text(40, -2, 'v0 =6 mV'); legend('wbar = 0', 'wbar = 0.48', 'wbar = 3.23')
 plot([0 50], [0 0], 'k:') % 0 axis
 plot([0 10], [1 1], 'b:') % "knee"
 % test components parts - for erfc(-1/V)
 r0=0.4; v0=6;  % r0=2; % test
 figure; plot(vpos, -v0./vpos, 'k--'); hold on; title(' [V0/V(t) -wbar] argument')  
 wbar=0.48; % {min} mean wt/link
  argm= r0*(-v0./vpos +wbar);
  plot(vpos, (-v0./vpos +wbar), 'b-');
   plot(vpos, (-v0./vpos +wbar), 'b.', 'MarkerSize', 14);
  VatZero= vpos(find((argm-0)<= 1e-2 ))
 wbar= 3.23; % {max} mean wt/link, for this pop'n
  argm= r0*(-v0./vpos +wbar);
  plot(vpos, (-v0./vpos + wbar), 'r-');
  plot(vpos, (-v0./vpos + wbar), 'r.', 'MarkerSize', 14);
  VatZero= vpos(find((argm-0)<= 1e-2 ))
 text(40, -2, 'vo =6 mV, ro =2'); legend('wbar = 0', 'wbar = 0.48', '', 'wbar = 3.23')
 plot([0 50], [0 0], 'k:') % 0 axis
 plot([0 10], [-1 -1], 'b:') % "knee"
 
 clear argm wbar VatZero 

 %% 2.1.23. test code for NM integration loop:
 fprintf('\n >> test logic code of integ over erfc ')
 %vpos=[0.1:1: 50];
 vpos=[-1: 0.5: 10];
 Scurve=zeros(1,length(vpos));
r0=0.05; v0= 6;  % now fixed threshold; r=1/(sqrt(2)*sigma) [of wts]
wbar= 3.23; % {max} mean wt/link, for this pop'n
 %wbar=0.48; % {min} mean wt/link
 %wbar=0.15; % min in alt form
 for i = 1:length(vpos)
    dyi=vpos(i);  % was delta-y of nn.
    dytest= sign(int32(dyi))   
    switch dytest
        case -1
            argm=0; In12tmp= 0; % avoid difficulties of erfc(neg)
        case 0
            argm=0; In12tmp= 0; % avoid divergence of 1/v
        case 1
           argm=r0*(wbar -v0/dyi); In12tmp= (1.0-erfc(argm)/2); % erfc form for wt; avoid v=0
    end
       %[i dyi argm]  % debug
        % In12tmp= (1.0-erfc(argm)/2);  % erfc form for wt; avoid v=0?
        %In12tmp= vm*(1.0-erfc(r0*(wbar -v0/dyi))/2);  % orign 
      %[i dyi argm In12tmp]  % debug
     Scurve(i) = In12tmp;
 end 
 
  figure; plot(vpos, Scurve);
  
  clear argm i dyi vpos vpos
  
 %% Appx 2.1.4 vs. erfc(1/v) form
   % erf(r0*v0) % const offset
vpos=[0.1:1: 50]; % assumes v > 0 !!
 %vpos=[-20:1: 50]; % debug : exclude neg argm.
r0=0.3
v0=6;  % now fixed threshold; r=1/(sqrt(2)*sigma) [of wts]
onerowp=ones(1,length(vpos));
 
wbar= 3.23; % {max} mean wt/link, for this pop'n
S2= 1.0 -erfc(r0*(wbar*onerowp -v0./vpos))/2;
 maxS = S2(end)
 figure; hold on; plot(vpos, S2, 'r--' ); %ylim([0 1]);  %xlim([-10, 40]); %ylim([0 1]);
 S3=1.0./(1+exp(r0*(v0*onerowp -wbar*vpos)) );
  plot(vpos, S3, 'k-')
  
wbar=1.0; %wbar=6.73; %wbar=0.48; % {min} mean wt/link
S2= 1.0 -erfc(r0*(wbar*onerowp -v0./vpos))/2; % raw form
  maxS = S2(end)
  plot(vpos, S2, 'b--' );
  title('S[v(t); w] via synaptic wt distribution; erfc[-1/V(t)]')
S3=1.0./(1+exp(r0*(v0*onerowp -wbar*vpos)) );
  plot(vpos, S3, 'k-')
  %legend('wbar = 3.23', 'Sigmoid', 'wbar = 0.48', 'Sigmoid'); %legend('mu = 3.23', 'mu = 0.48'); 
  legend('wbar = 6.73', 'Sigmoid', 'wbar = 1.0', 'Sigmoid');
text(40, 0.1, 'r = 0.3; Vo = 6 mV');  %  label params
  plot([0 50], [0.5 0.5], 'b:') % 50% activation
  xlabel('V input (mV) '); ylabel('pulse rate output (s^-1) ');
  
clear r0 r1 onerow onerowp v vpos S1 S2 S3 maxS
 
 %% Appx 2.1.3 Sigmoid fn form: Sigmoid + Weight (wbar) forms; plot S[v] curves
 % vthreshold form
r0=0.5; v0=6; % % "lower" limit,  etc
v=[-10:1: 50]; % voltage array
onerow=ones(1,length(v)); % needs vector of 1's for sigmoid (over time span)
figure;  title(' basic S[v] curve ');
plot(v, 1.0./(1+exp(r0*(v0*onerow - v )) ), 'k-'); hold on 
  %r0=0.4; %  ..  & "upper" limit
  %plot(v, 1.0./(1+exp(r0*(v0*onerow - v )) ), 'k--');
 xlabel('V input (mV) '); ylabel('pulse rate output (s^-1) ');
  %legend('r: 0.3') %legend('r: 0.3','  0.5', '  0.75', '  0.9')
  text(20, 0.1, 'r = 0.5; Vo = 6 mV'); xlim([-20, 60])
   title(' JR model: S[ ..(wbar*V(t) - Vo)]')
 
 % wt form
  wbar=0.48; %wbar=0.16; %  % ("lower" limit) mean wt/link of neurn pop'n
  S3=1.0./(1+exp(r0*(v0*onerow - wbar*v )) ); % is correct
  %figure
  plot(v, S3, 'b-')
  
  wbar=3.23; %wbar=1.0; %  % upper limit (mean wt/link of neurn pop'n
  S3=1.0./(1+exp(r0*(v0*onerow - wbar*v )) );
  plot(v, S3, 'r-')
  wbar=6.73; %wbar=1.0; %  % upper limit (mean wt/link of neurn pop'n
  S3=1.0./(1+exp(r0*(v0*onerow -wbar*v)) );
  plot(v, S3, 'r--')
   legend('r= 0.5  wbar=1', ' wbar= 0.48', ' wbar= 3.23', ' wbar= 6.73')
  %legend('r= 0.5, wbar=1',' wbar= 0.48', ' wbar= 2.73')
  plot([0 60], [0.5 0.5], 'b:') % 50% activation
  clear onerow wbar S3 r0
  
%% Appx 2.1.4 S S[V; wbar] for [ y-vec(1)]
 % as calc byNM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
onerow=ones(1,length(y)); % needs vector of 1's for sigmoid (over time span)
r0=0.3; wbar=0.48; %0.52; 
% components
figure; subplot(2,1,1); plot(y(2,:) ); title (' y-e: y1(t) '); xlabel('t (0.1 ms) ');
 ylabel('y (mV) '); axis([0 50000 0 max(y(2,:))+10]);
 subplot(2,1,2); plot(vm./(1+exp(r0*(v0*onerow -y(2,:) )) ), '--'); hold on % orig form
 %subplot(2,1,2); 
 plot(vm./(1+exp(r0*(v0*onerow - wbar*y(2,:) )) ) ); legend('Sigmoid', 'S[V-wt form ]')
 axis([0 50000 0 vm+1]); text(30000, 2, 'r=0.45; w = 1')
title (' S[y(t)-wt] '); xlabel('t (0.1 ms) '); ylabel('pulse rate (s^-1) ');

% Appx 2.1.4a S fn for y-vec(2) / 2nd plot
 % as calc byNM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
onerow=ones(1,length(y)); % needs vector of 1's for sigmoid (over time span)
% components
figure
subplot(2,1,1); plot(y(8,:) ); title ('y-i: y2(t) '); xlabel('t (0.1 ms) ');
 ylabel('y (mV) '); axis([0 50000 0 max(y(8,:))+10]);
 subplot(2,1,2); 
 plot(vm./(1+exp(r0*(v0*onerow -y(8,:) )) ), '--'); hold on % orig form
  %subplot(2,1,2); 
plot(vm./(1+exp(r0*(v0*onerow - wbar*y(2,:) )) )) % wt form
 legend('Sigmoid', 'S[V-wt form ]')
axis([0 50000 0 vm+1]); text(30000, 2, 'r= 0.45; w = 1')
title (' S[y(t)-wt] '); xlabel('t (0.1 ms) '); ylabel('pulse rate (s^-1) ');


%% Appx 2.1.5 S[] fn for erfc[1/dy-vec(1)] "rough" form
 % as calc by NM_JR12 codes
 fprintf('\n   plot V(t) & S[v] via erfc(1/v : N(w)) - raw form of erfc \n')
 % resultant firing rate, due to v's : ~ A*a*S(v)
  % needs heaviside step fn: " Y = zeros(size(X),'like',X);  Y(X > 0) = 1; "
r0 =0.4; v0=6.0; vm =5; % params. of S[V}, cf. sim codes
wbar= 1; %wbar= 0.52  %wbar= 3.23;  %0.48  % for node-1 in 6-star cluster
% components  % dyi=y(2,:);  % V(t)
dyi=(y(2,:)-y(3,:)); % choose  y-e or dy
onerow=ones(1,length(dyi)); % needs vector of 1's for sigmoid (over time span)
stepfn= zeros(1,length(dyi)); 

 stepfn(dyi>0)=1.0; % apply 
  %dyi=stepfn.*dyi; % XXXapply heaviside step fn to dy; exclude neg dy
figure; subplot(2,1,1); plot(dyi ); %title (' y-e-1(t), S[erfc(1/v)]');   %
title (' dy1(t), S[erfc(1/v)] : 2NM, alpha'); 
 xlabel('t (0.1 ms) '); ylabel('y (mV) '); xlim([0 10000]);
%
S2= vm*(1.0 -erfc(r0*(wbar*onerow -v0./dyi))/2);
S2=stepfn.*S2; % heaviside step fn; exclude neg dy
 subplot(2,1,2); 
 plot(vm./(1+exp(r0*(v0*onerow - dyi )) ), '--'); hold on % orig sigmoid form - y-e
 %plot(vm./(1+exp(r0*(v0*onerow -(y(2,:)-y(3,:)) )) ), '--'); hold on % orig sigmoid form - dy
 %subplot(2,1,2); 
 plot(S2); xlim([0 10000]); subplot(2,1,2); ylim([0 6]); 
text(30000, 0.5, 'wbar = 1; r = 0.4'); legend('Sigmoid')
title (' S[y(t)-wt] '); xlabel('t (0.1 ms) '); ylabel('pulse rate (s^-1) ');

%% Appx 2.1.5a S[] fn for erfc[- dy-vec(2)] : raw form
 % as calc by NM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
r0 =0.4; % params. cf. sim codes
wbar=1; %wbar=3.23;
onerow=ones(1,length(y)); % needs vector of 1's for sigmoid (over time span)
stepfn= zeros(1,length(y));
% components
dyi=(y(8,:)-y(9,:)); stepfn(dyi>0)=1.0; 
dyi=stepfn.*dyi; % apply heaviside step fn; exclude neg dy
figure; subplot(2,1,1); plot(dyi ); title (' dy2(t) '); xlabel('t (0.1 ms) ');
 ylabel('y (mV) '); xlim([0 20000]);
S2= 1.0 -erfc(r0*(wbar*onerow -v0./dyi))/2;
 % subplot(2,1,2); plot(vm./(1+exp(r*(v0*onerow -y(8,:) )) )) % orig form
subplot(2,1,2); plot(S2) % wt form
xlim([0 20000]); ylim([0 1.1])
text(3500, 0.7, 'wbar = 1; r= 0.4')
title (' S[y(t)-wt] '); xlabel('t (0.1 ms) '); ylabel('pulse rate (s^-1) ');

clear dyi S2 stepfn onerow 

%% Appx 2.1.5b S[] fn for erfc[- dy-vec(2)] : stepfn form
 % as calc by NM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
r0 =0.4; % params. cf. sim codes
wbar=1; %wbar=3.23;
onerow=ones(1,length(y)); % needs vector of 1's for sigmoid (over time span)
stepfn= zeros(1,length(y));
% components
dyi=(y(8,:)-y(9,:)); argm= r0*(wbar*onerow -v0./dyi);
stepfn(argm>0)= 1.0; % need to exclude neg argm
 %Nzeros = length(find(stepfn==0) ) % debug
stepfn(dyi<0.01) = 0.0; % and exclude y=0 divg. nb. not "precisely" zero
 %Nzeros = length(find(stepfn==0) )
 %dyi=stepfn.*dyi; % apply heaviside step fn; exclude neg dy
figure; subplot(3,1,1); plot(dyi ); title (' dy2(t) '); xlabel('t (0.1 ms) ');
 ylabel('y (mV) '); xlim([0 2e4]);

 subplot(3,1,2); plot(stepfn); title('StepFn'); xlim([0 2e4]); ylim([0 1.1]); %
   %plot(argm); title('argm'); xlim([0 2e4]); %ylim([0 1.1]); % alt.
   %S2= 1.0 -erfc(r0*(wbar*onerow -v0./dyi))/2;
S2= stepfn.*(1.0 -erfc(argm)/2 ); % modulated by step fn
S2=S2.*stepfn; % filter "problems"
subplot(3,1,3); plot(S2) % wt form
 ylim([0 1.1]);  xlim([0 2e4]); hold on;
 plot(1.0./(1+exp(r0*(v0*onerow - dyi )) ), 'k-'); % basic sigmoid S[v] form
text(3500, 0.7, 'wbar = 1; r= 0.4')
title (' S[y(t)-wt: erfc(1/V] '); xlabel('t (0.1 ms) '); ylabel('pulse rate (s^-1) ');
legend('erfc form', 'sigmoid')

clear dyi S2 stepfn onerow Nzeros argm

%% Appx 2.1.5b.1 S[] fn for erfc[- dy-vec(2)] : stepfn form : test code
 % as calc by NM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
r0 =0.4; % params. cf. sim codes
wbar=1; %wbar=3.23;
onerow=ones(1,length(y)); % needs vector of 1's for sigmoid (over time span)
stepfn= zeros(1,length(y));
% components
dyi=(y(8,:)-y(9,:)); argm= r0*(wbar*onerow -v0./dyi);
stepfn(argm>0 & dyi>0.01)= 1.0; % need to exclude neg argm & V~0 :: OK
   % and exclude y=0 divg. nb. not "precisely" zero
Nzeros = length(find(stepfn==0) ) % debug
 %Nzeros = length(find(stepfn==0) )
 %dyi=stepfn.*dyi; % apply heaviside step fn; exclude neg dy
figure; subplot(3,1,1); plot(dyi ); title (' dy2(t) '); xlabel('t (0.1 ms) ');
 ylabel('y (mV) '); xlim([0 2e4]);

 subplot(3,1,2); plot(stepfn); title('StepFn'); xlim([0 2e4]); ylim([0 1.1]); %
   %plot(argm); title('argm'); xlim([0 2e4]); %ylim([0 1.1]); % alt.
   %S2= 1.0 -erfc(r0*(wbar*onerow -v0./dyi))/2;
S2= stepfn.*(1.0 -erfc(argm)/2 ); % modulated by step fn
S2=S2.*stepfn; % filter "problems"
subplot(3,1,3); plot(S2) % wt form
 ylim([0 1.1]);  xlim([0 2e4]); hold on;
 plot(1.0./(1+exp(r0*(v0*onerow - dyi )) ), 'k-'); % basic sigmoid S[v] form
text(3500, 0.7, 'wbar = 1; r= 0.4')
title (' S[y(t)-wt: erfc(1/V] '); xlabel('t (0.1 ms) '); ylabel('pulse rate (s^-1) ');
legend('erfc form', 'sigmoid')

clear dyi S2 stepfn onerow Nzeros argm


%% Appx 2.1.5bXXX S[] fn for erfc[- dy-vec(2)] : "more careful"form {switch}
 % as calc by NM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
r0 =0.4; wbar=1;  % params. cf. sim codes
 v0=6; % % "lower" limit,  etc
v=[-10:1: 20]; % voltage array
dy1=y(2,:)-y(3,:); % delta-y 
figure; subplot(2,1,1); plot(dy1 ); title (' dy1(t) '); xlabel('t (0.1 ms) ');
 ylabel('y (mV) '); xlim([0 10000]);
%
%
Sarray=zeros(1,length(dy1));
%figure
 subplot(2,1,2); hold on  % iterate thru S[v(t-i)]
for i = 1:length(v)
     %dyi=dy1(i);  %
    dyi=y(2,i)-y(3,i); % delta-y
    argm=r0*(wbar -v0/dyi);
    dytest= sign(int32(argm));  %dytest= sign(int32(dyi));   
       switch dytest
          case -1
            argm=0; S2= 0; % avoid difficulties of erfc(neg)
          case 0
            argm=0; S2= 0; % avoid divergence of 1/v
          case 1
           argm=r0*(wbar -v0/dyi); S2= (1.0-erfc(argm)/2); % erfc form for wt; avoid v=0
       end % logic
         % S2 % debug
       Sarray(i)=S2;
    % plot(i, S2,' -.')   
end % v loop
  plot(Sarray)
xlim([0 10000]);       
text(8000, 0.8, 'wbar = 1; r= 0.4')

% clear dy1 dyi dytest argm S2 Sarray

%% Appx 2.1.5cXXX S[] fn for erfc[- dy-vec(2)] : "more careful"form {H-stepfn}
 % as calc by NM_JR12 codes
 % resultant firing rate, due to v's : ~ A*a*S(v)
r0 =0.4; wbar=1;  % params. cf. sim codes
 v0=6; % % "lower" limit,  etc
v=[-10:1: 20]; % voltage array
dy1=y(2,:)-y(3,:); % delta-y 
figure; subplot(2,1,1); plot(dy1 ); title (' dy1(t) '); xlabel('t (0.1 ms) ');
 ylabel('y (mV) '); xlim([0 10000]);
%
%
Sarray=zeros(1,length(dy1)); Sarray=zeros(1,length(dy1));
%figure
 subplot(2,1,2); hold on  % iterate thru S[v(t-i)]
for i = 1:length(v)
     %dyi=dy1(i);  %
    dyi=dy1(i);  %y(2,i)-y(3,i); % delta-y
    argm=r0*(wbar -v0/dyi);
    argtest= sign(int32(argm));  %dytest= sign(int32(dyi));   
        S2= 0; % avoid difficulties of erfc(neg) & divergence of 1/v
        if argtest>0
           %argm=r0*(wbar -v0/dyi); 
           S2= (1.0-erfc(argm)/2); % erfc form for wt; avoid v=0
        end % logic
         % S2 % debug
       Sarray(i)=S2;
    % plot(i, S2,' -.')   
end % v loop
  plot(Sarray)
xlim([0 10000]);       
text(8000, 0.8, 'wbar = 1; r= 0.4')

% clear dy1 dyi dytest argm S2 Sarray

%%  2.1.7 plot S[v] curve - vs erf[logN(Vthr)] form : from #2.5, above
r0=0.05; v0= 6; % defaults for JR; etc
v=[0.01:1: 70]; % voltage array: needs >=0
onerow=ones(1,length(v)); % needs vector of 1's for sigmoid (over time span)
figure;  plot(v, 1.0./(1+exp(0.4*(v0*onerow - v )) ), 'k-');  hold on;  % fix at r=0.4 
 xlabel('V input (mV) '); ylabel('pulse rate output (s^-1) ');
  text(20, 0.1, 'r = 0.05; Vo = 6 mV'); %xlim([-10, 40]); ylim([-0.05, 1.0])
  %text(60, 1.5, 'Qmax = 5 Hz');    

  S1= (erf(r0*(v -v0*onerow)) +erf(r0*v0) )/2; % for N(Vthr)
 plot(v, S1, 'k:' );  % xlim([-10, 40]); ylim([-0.05 1]);
  %legend(' S[v]', 'erf(v)');  title('Sigmoid and erf responses')
  %pause
  S2= 1*(onerow + erf(r0*(log(v) -v0*onerow)) )/2; % for logN(Vthr)
 plot(v, S2, 'b--' );  % xlim([-10, 40]); ylim([-0.05 1]);
  legend(' S[v; r=0.4]', 'erf(v): N', ' erf(lnV): logN');  
  title('Sigmoid and erf responses: for N(.) & logN(.)')
  
%% tests
  %figure; plot(log(v))
  hold on; plot(log(v) -v0*onerow, 'b--')
  %% new r
r0=0.3;
plot(v, 1.0./(1+exp(r0*(v0*onerow - v )) ), 'r-')
S1= (erf(r0*(v -v0*onerow)) +erf(r0*v0) )/2;
plot(v, S1, 'r--' );  
text(20, 0.05, 'r = 0.3; Vo = 6 mV', 'Color', 'r'); 
 
clear onerow S1 r0

%

%% Appx 2.2 S[v] curve, threshold form : cf wbar varies xxx wrong
  % from #2.1, above
r0=0.5; v0=6; % defaults for JR; etc
v=[-10:1: 60]; % voltage array
onerow=ones(1,length(v)); % needs vector of 1's for sigmoid (over time span)
figure;  title(' S[v] curves ');
plot(v, 1.0./(1+exp(-r0*(v - v0*onerow )) ), 'k-'); hold on
 xlabel('V input (mV) '); ylabel('pulse rate output (s^-1) ');
  %legend('r: 0.3') %legend('r: 0.3','  0.5', '  0.75', '  0.9')
  text(20, 0.1, 'r = 0.5; Vo = 6 mV'); %xlim([-10, 40])
  title(' JR model: S[ ..(V - wbar*Vo)]')
  %text(60, 1.5, 'Qmax = 5 Hs');  
  
 % wt form
  wbar=0.48; % %  % ("lower" limit) mean wt/link of neurn pop'n
  S3=1.0./(1+exp(-r0*(v  -v0*wbar*onerow)) ); % Not correct
  plot(v, S3, 'b-')
  
  wbar=2.73; %wbar=1.0; %  % upper limit (mean wt/link of neurn pop'n
  S3=1.0./(1+exp(-r0*(v -v0*wbar*onerow)) );
  plot(v, S3, 'r-')
  wbar=6.73; %wbar=1.0; %  % upper limit (mean wt/link of neurn pop'n
  S3=1.0./(1+exp(-r0*(v -v0*wbar*onerow)) );
  plot(v, S3, 'r--')
  legend('r= 0.5  wbar=1', ' wbar= 0.48', ' wbar= 2.23', ' wbar= 6.73')
  %legend('r= 0.3', 'r= 0.5 wbar=1',' wbar= 0.16', ' wbar= 1.0')
  plot([0 60], [0.5 0.5], 'b:') % 50% activation
  clear onerow wbar S3 r0
  
%% raw dy's
figure;
dyi=(y(2,:)-y(3,:));  
plot(dyi, 'k--'); hold on
dyi=(y(8,:)-y(9,:)); plot(dyi, 'b--');
xlim([0 5000])

 %% Appx 3.0 test trig ienties
 theta =[0: 0.1: 3*pi]; dphi=pi/2;
 figure; subplot(3,1,1); plot(theta, sin(theta), 'k--'); hold on; title(' sin(x, 2x)')
 plot(theta, sin(2* theta), 'b'); plot(theta, abs(sin(theta)), 'r--'); grid on
  legend('sin(x)', 'sin(2x)', 'abs sin(x)')
 subplot(3,1,2); plot(theta, 2*sin(theta).*cos(theta), 'k'); grid on
 subplot(3,1,3); plot(theta, 2*sin(theta).*cos(theta), 'k'); title(' 2sin(x)*cos(x)')
 xlabel('theta (rad)')
 
 figure; subplot(3,1,1); plot(theta, cos(theta), 'k--'); hold on; title(' cos(x, 2x)')
 plot(theta, cos(2* theta), 'b'); plot(theta, abs(cos(theta)), 'r--'); grid on
  legend('cos(x)', 'cos(2x)', 'abs cos(x)')
 subplot(3,1,2);  plot(theta, cos(2* theta -dphi), 'b-'); hold on
   plot(theta, abs(cos(2*theta -dphi)), 'r--'); grid on; legend('cos(2x -pi/2)', 'abs cos(2x -pi/2)')
 subplot(3,1,3); plot(theta, 1.0-sin(theta).^2, 'k'); title(' 1 - sin^2(x)')
 xlabel('theta (rad)')
 
 clear theta
 
 
 