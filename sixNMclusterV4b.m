% B {Pailtjorpe, U Sydney. sixNMclusterVx.m, (25/7/21 - 24)
 %Jansen-Rit model, 6 node cluster at A-10 (pFC), classic sigmoid form,
    % same f for each NM. (almost identical nodes)
  % Runga-Kutta integration use RK4 on JR, fn JR12(.., NNinputs), with stim. forcing
% Dependencies:  JR12a.m function for RK4 DE integrator; WaveTrains.m for stimuli set up 
              %  is RHS of DE: function dydx as JR12: incl Ip, p forcing, NNcoupling
  % Data file:  AdjMarmoset_Rescale2b.csv % link wts(116 x 116 array)
     % V1. dev.(23/5/21)  %  based on sixNMloopRK4_JR12V2b.m
      % Updade, V2: use sigmoid filter on coupling  (12/6/21)
         % V2b. modulate stim Ip+ran to selected nodes only (7/5/21)
      % V3 use JR12a(t, y, Ip, p, I12, a, b, A, B, C) : pass params to tune freq band;
      % V3a. pass Stim, pran via Sigmoid: pulse->V transfer (4/7/21)
      % V4. pass r to fn J12a.m; vary at ea NM
      % V4.1 vary Gain[s] ~ syn strength, per NM; pass Cfac (C1, C3) to fn J12a.m (24/7)
      % V4a. add signal delays, via d(i,j)
      % V4b. uses S[v] sigmoid form;  vlocal effects on dij ** THis version
% Cells (sections)
% 1. ADj set up
% 2. parameters for WC/JR model
% 3. set up stimuli
% 4.0  RK4 DE solver, 1st order, for n variables
% 5.0  various PLOT OUTPUT

% clear all; close all

% 1.0 Set up
addpath('/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain/Models'); % for other model codes
fprintf('\n 6 Neural Masses, JR/RK4, star cluster; wt Adj; S(v): coupling, Pulse, ran \n')
fprintf('\n  alpha band;  Sigmoid; vary r & Gain; +  stim; code: v4b \n')
 % Adj, for coupling nodes
nn=6  %nn=6 % number of neural masses (nodes)
 % examine star cluster about Out-hub A-10  (marmoset logbook, p58, 66 ("v1")
  % n-node, Marmoset star clusters: at 7ms: 6 nodes; @7.5 & 8.1 ms: 8 nodes

  %%  1.0.1 Approx  Adj2, .cf Appx. below: this is similar to Adj2, with rounding
  % fprintf('\n   use approx Adj (rounded)  \n')
  % Adj= [0.0,  9.0, 0.5, 1.0, 4.6, 8.5; ...  % #1 is out hub;
  %      0.25, 0.0, 0.5, 0.0, 0.0, 0.5; ...  % cf Workbook 3, p18a for list of Links
  %      0.5, 4.9, 0.0, 2.0, 0.25, 1.0; ...  % scale wt wrt to 1k (LNe)
  %      0.75, 0.0, 1.2, 0.0, 1.6, 0.0; ...
  %      0.70, 0.0, 0.0, 0.1, 0.0, 0.0; ...  %  
  %      0.40, 1.5, 1.0, 0.0, 0.1, 0.0];     %  
  % cf. Appx 2, below for check of Adj against orig Anew(LNe)
  % A(4,2), A(5,6) was missing; A(6,3) was too big; others are close
  
%% 1.0.2 unit Adj
fprintf('  test Adj=1:  \n')
Adj=ones(nn);
 % elim diagonals
 for i=1:nn
     Adj(i,i)=0;
 end
 Adj
 length(find(Adj(:))) % 30 links
 % 1.02 zero Adj
  %fprintf('  test Adj=1:  \n')
  %Adj=zeros(nn);

%% 1.0.3 Zero Adj: no coupling
  fprintf('  test Adj=0: no coupling  \n')
  Adj=zeros(nn);
  
%% 1.0.3  d(ij) array - cf. Appx. 3.2.1  set d in bins of 0.5 mm (approx)
 fprintf('  d(ij) in 0.5mm bins:  \n')
Ad6A = [...
     0     3.5   3.0   2.5   2.5   2.5; ...
     3.5   0     1.5   3.5   4.0   2.5; ...
     3.0   1.5   0     2.0   3.0   2.5; ...
     2.5   3.5   2.0   0     2.0   3.5; ...
     2.5   4.0   0     2.0   0     3.0; ...
     2.5   2.5   2.5   3.5   3.0   0];
%% 1.0.4  d(ij) array - cf. Appx. 4.1.2  set d in bins of 0.1 mm (approx)
fprintf('  d(ij) in 0.1mm bins:  \n')
Ad6A = [...
     0     3.4  2.9  2.4   2.4  2.5, 4.3, 5.0; ...
     3.4   0    1.7  3.4   4.1  2.7, 0.0, 4.7; ... % 0 if no link present
     2.9   1.7  0    1.9   2.8  2.6, 0.0, 3.2; ...
     2.4   3.4  1.9  0     1.8  3.4, 3.2, 2.8; ...
     2.4   4.1  0    1.8   0    2.8, 2.1, 3.2; ...
     2.5   2.7  2.6  3.4   2.8  0.0  3.6, 4.8; ...
     4.3,  4.9, 3.6, 3.2,  2.1, 3.6, 0.0, 2.4; ...
     5.0,  4.7, 3.2, 2.8,  3.2, 4.8, 2.4, 0.0];   % checks ok.
 %% 1.0.5  zero dij
 fprintf('  d(ij) =0: no signal delays  \n')
 Ad6A = zeros(nn);

%% 1.1 Adj from Marmoset (LNe) data
 % cf.Appx 2. in orig order: cf. logbook, p58
  fprintf('  use marmoset LNe Adj(2):  \n')
clusterlist=[2, 26, 25, 48, 32, 3]; % as plotted in x-y plane, viewed from Ant.
Anew= csvread('AdjMarmoset_Rescale2b.csv');  % based on LNe 
 Aclust = Anew(clusterlist, clusterlist);
Adj2=Aclust/1000  % scale used for 6 NM model
Adj=Adj2;   % for main code
nLinks= length(find(Adj2(:))) % 29 entries
clear nLinks Anew Aclust Adj2 clusterlist
 %Adj=zeros(nn); % test indep NM.
 %Adj=ones(nn); % test unit coupled NM. nb. diagonals (self links) present

% nb. set approx dij in Ad6A, above & at Appx 3.2.1

%% 2.0 Parameters: tune freq. band
 fprintf('  : cf Marmoset out hub, A10 (3-D cluster, 6); beta band; v-local \n')
 %nLinks = length(find(Adj(:))) % check # Links: was 22 [nb. some v small & omitted], now 29
 %clear nLinks
% time axis & IC (0 mV)
x0=0; xf=20.0; % t range (sec)
h = 0.0001;   % t step (1 ms; [& 0.5, 0.25 ms, for debug])
x = [x0:h:xf];  % t domain array
 %time =x/(1000);  %time =x/(1000*h); % time (sec) as row vec
tstart=10001; %tstart=2000; %1000; % for plots (@ 1, 0.5, 0.1 ms time steps)
%y = zeros(6,length(x)); % row vec for each variable
yic1=[0.0; 0.0; 0.0; 0.0; 0.0; 0.0]; % IC at t=0 for 1 node (6 DEs),  col vec, 6 rows
yic=yic1; % add additional nodes
for in=1:nn-1
 yic=[yic; yic1]; % add 6 rows per extra node; 12, 18 DEs for 2, 3 ... nodes
end
 clear yic1
%
% Model parameters (orig JR) : needs to be in the fn subroutine!
v0=6;    % (mV) midpoint (switching trheshold V) of sigmoid function, for V -> pulse conversion 
vm=5;    % (s^-1) max firing rate for sigmoid fn
 %r=0.3;   % activation rate;  Steepness of sigmoid function 
  % rvec= [0.3 0.3 0.3 0.3 0.3 0.3];  % pass to indiv nodes, via fn J12a.m
%rvec= [0.4 0.4 0.4 0.4 0.4 0.4] % test all alpha, or beta (low)
   %rvec= [0.45 0.45 0.45 0.45 0.45 0.45] % test - in beta
rvec= [0.5 0.5 0.5 0.5 0.5 0.5] % test: all beta: most NM dont osc! (at 0.56)
 %rvec= [0.56 0.40 0.49 0.46 0.40 0.49] % set r ~ sqrt(vol); cf w/b#3, p80: tune in beta: ok
   %rvec= [0.56 0.40 0.47 0.43 0.38 0.47] % revised r ~ sqrt(N); cf w/b#4, p42: tune in beta band
  %rvec= [0.42 0.30 0.37 0.34 0.30 0.37] % ditto: tune in alpha band: corrected (23/7/21)
%rvec= [0.45 0.32 0.38 0.35 0.30 0.38] % revised: tune in alpha band: corrected (18/9/21)
        % cf w/b#4, p42.
% tuned to alpha/ beta band
C=300 %C=250 %C=135   % basic (scale) connection strength, within NM % defauly J-R value         
R12=100  % basic NM - NM coupling coeff (ie. link wt) (cf. Goodfellow)
C1=C;     % Pyramidal to excitatory interneuron connection 
C2=0.8*C; % (108)   Excit interneuron to pyramidal connection [feedback]
C3=0.25*C; % (33.75)  Pyramidal to inhib interneuron connection [weaker feedback]
C4=0.25*C;

 % tune feedback gain (baed on WtIn, kIn)
CfacVec= [1.0 1.0 1.0 1.0 1.0 1.0]; % pass to fn J12a.m
 %CfacVec= [6.33; 1.0; 7.0; 4.67; 2.67; 1.33]; % follows wtIn/kIn for 6x6(wb3, p88, 90): Poor
%CfacVec= [4.6; 3.07; 1.03; 1.52; 1.0; 1.92];  % ?recalc for 8x8, wb#3, p101 (4/8/21)
  %CfacVec= [6.39; 4.26; 1.0; 2.11; 1.39; 2.66; 1.76; 2.76];  % ?recalc again (wb#3, p104 & XL)
CfacVec    % Assign in J12a:C1=C1*Cfac;  C3=C3*Cfac; 

% signal delay: set by array Ad6a, set in App3.2, based on d(i,j) % d12=1.5-4;  % eg. mm / ms
 %Ad6A=zeros(nn); % debug; no delays (base case)
A=3.25;  B=22;  % Max PSP amplitude for pyramidal (A), inhib (B) & excit (kA * A) interneurons
a=100; b=50;  % (s^-1) / default JR; set to osc at 8Hz (taus: 20, 20 mS)
%
taue= 10.0; taui= 10; % set Time Const (excit, inhib) (ms)  :  alpha, beta bands
    % move freq band; keep a*A & b*B const
a=1000/taue; b=1000/taui; % rate const (s^-1)
 %A=3.25*100/a; % rescale, to keep a*A const. dont tune B (lower C)?
 %B=22.0*50/b;  % ?dont tune B (more inhib): need lower C then.
[taue taui] % debug / use default A,B
[a b]
[A B C]
[round(a*A) round(b*B)]  %[a*A b*B]
ratioIE =  round(b*B)/round(a*A)
tstart = 40001;  % stim at 4.0 s (0.1 ms steps0

clear ratio*


% >>>>>>>>>>>>>>

%% 3.0 Zero pulse, ran or wave :: Base case
fprintf('\n >> zero inputs: 1 sd G ran noise \n')
Vyt =zeros(1,length(x));
pran=zeros(1,length(x)); 
Ip=zeros(1,length(x));  
Ip6=zeros(1,length(x)); IpOn=zeros(nn,1); IpOnB=zeros(nn,1); % need IpOn also
Pd=zeros(1,length(x)); Pd2=zeros(1,length(x)); %  pulse density fn - cf #3.6
% noise (BC)
  %pran = 0.1*rand(1,length(x)); % low noise;  Uniform ran noise:
pran = 1.0*randn(1,length(x)); % Gaussian noise, pos/neg 
figure;  plot(pran, 'Color', [0.9 0.9 0.9]); % ylim([min(pran)-1  max(pran)+1]); xlim([2e4 3e4])
 %  text(2.06e4, 0.3, 'U ran noise, \mu = 0.1 mV')
 % text(2.06e4, 0.3, 'G ran noise, \mu = 0.1 mV') 

%% 3.1 set up Stimulus pulse (spike) train
  %fprintf('  * stimulate node 2&6 (A32V, A11 ) both, equally:  \n')
 %fprintf('  * stimulate node 6 (A11 ) only:  \n')
 %fprintf('  * stimulate node 2 (A32V ) only:  \n')
%fprintf('  Pulse train stimulus - 10 mV; to NM #2 amd/or #6  \n')
stim_freq=16  % nb. also need to change npulse, to get 1 sec long stimulus
Ip=zeros(1,length(x)); % 1st input (NM#2):excit curents impulse train: input to pyramidal population
Ip6=zeros(1,length(x)); % pulse train for 2nd input (?dphi shifted), to node 6
tstart= 1000; % wait #steps past start up transients; adjust to 1 sec.
 %tdelay= round(1000/stim_freq)  % period: needs to be integer, for indexing; assumes h = 1ms
tdelay= round(1/(h*stim_freq)) % period: needs to be integer, for indexing; allow h varies
twidth=25 % width of pulse (dt steps; eg ms)
dphi= -10 % phase delay (dt steps; ms) of 1st pulse [& pulse train]
         % or increment (ramp) for ea pulse; or for 2nd pulse; etc.
npulse=5 % pulse train % set length of pulses: eg 1 sec
for i=0:npulse-1
   %Ip(tstart+i*tdelay+1)=1.0; Ip(tstart+i*tdelay+2)= -1.0; % symm. sharp dbl spike, pos/neg
   %Ip(tstart+i*tdelay+1)= -1.0; Ip(tstart+i*tdelay+2)= 1.0; % symm. sharp dbl spike, neg/pos
 Ip(tstart+ i*tdelay +dphi +1)=1.0; % simple phase shift; sharp pulse
 Ip6(tstart+ i*tdelay +dphi +1)=1.0; % asymm. sharp single unit spike, phase shifted, to Node 6
   %Ip(tstart+ i*tdelay +dphi +1: tstart+ i*tdelay +dphi +1+twidth-1)=1.0; % wider pulse

end
% other inputs; set scale
 %Ip=zeros(1,length(x)); Ip6=zeros(1,length(x)); % debug: test zero excit'n. : just follows the noise!
amplP=10.0  %mV % nb.max rate of 5Hz at ~25mV
 Ip=amplP*Ip; Ip6=amplP*Ip6; %Ii=60*r*I; % excit & inhib [set Ii in ext fn]
figure; plot(Ip); title('stimulus pulse train ')
 % p=zeros(1,length(x)); % debug; no ran noise
% set input channels, IP() at #3.7 below

 %% 3.2 set up random Noise input - filtered spectrum
fprintf('  + add rand noise (bandpass) \n')
fs=1/h; % sampling frequency
L=length(x); % length of signal
%pran = 1.0.*randn(1,L); % Gaussian (0,0.1): Excit'y input (av 0.1 mV)
pran = 0.1.*rand(1,L); % Uniform ran noise: contains all frequencies!
  %figure; plot(pran); title('test: raw ran stimulus ') 
% FFT spectrum
Y=fft(pran);
P2=abs(Y/L); P1=P2(1:floor(L/2+1)); P1(2:end-1)=2*P1(2:end-1); % 2-sided spectrum
f=fs*(0:(L/2))/L;   % ?? check L/2+1 is integer?
  %figure; plot(f,P1); title('spectrum of U ran noise'); xlabel(' f (Hz)')
   %title('spectrum of "100Hz" U ran noise')
% contains most frequencies !! 
% band pass filter
Pfilter=P1(501: 1501);  % 100 Hz - 300 Hz;  cf f vec.; 1k entries
pran100=fft(Pfilter); % nb. complex, want real part [by default]
pran100=real(pran100); pran100(1)=0; % 1st entry is spurious?  nb. row vec.
pran100=abs(pran100); % pos values
%figure; plot(pran100); title('filtered U ran noise')
% [mean(pran100) std(pran100)]
pran100=0.1*pran100/mean(pran100); % renormalise
  % pran100=(0.1/std(pran100))*pran100; % alt. renormalise
 %figure; plot(pran100); title('filtered U ran noise, renormalised ')
 %std(pran100)
tstart= 1000; tend=xf/h+1; % ms scale
% fill array to full time
pran=zeros(1,tend+1);
for i=1:xf-1
    pran(1000*i: 1000*i+1000 ) = pran100; % insert 1k+1 filtered blocks
end
pran=pran(1:end-1); % was +1 too long?
 %pran= [pran, 0]; % need N+1 elements (eg. 5001)
 % figure; plot(pran); title(' U ran noise stimulus (100:300Hz band) ')
clear P1 P2 Pfilter pran100 Y f L

% >>>>>>>>>>>>>>

%% 3.3 set up travelling Wave :: or use WaveTrain.m code now
  % cf. code: VisMarmosetBrainWaveV1.m
 fprintf('  + add wave train: 2 mV  \n')
Vwm= 10.0 % mV
istart=200001 % turn wave on at 1 sec; or add phase shift
istop= 250001 % turn wave off after 5 sec
lambda = 30.0; % mm  %lamda = 0.03 (m) ~ brain size (marmoset)
  % Set the travelling wave Frequency:
wfreq = 15.8/1000 % wave frequency (ms^-1)  %freq = 33.0; %  Hz (s^-1)
 % phiRad=PhiMS*2*pi/period
phi=2.1; % phase delay (lag, in radians) 
  % choose phase to line up start of sin wave & max near y(A10) location
  % additional phase shift to aling with dy(t) calc: cf. plots at #4.1, below
dphi= 425 % (in t-steps: 0.1ms) extra phase diff to align with delta-y osc
pdhi=dphi*pi/180; % in rad
Wspeed = wfreq*lambda % (m/s or mm/ms); here 0.5/1/2 m/s
 %t = 0; % ms
 % row vec x is the t-line
 %tline= [0: 1 : 1000]; %t axis [0 :1ms steps: 1 sec]
 tline = x*1000.0; % convert ms to sec, to match pulses & DE solver: for 1 msec steps
 ypt =16.0; % peak near this coord (A-P)
Vyt =zeros(1,length(x));
 %Vyt(istart:end) = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*(tline(1:end-istart+1) -phi +dphi) ); % start at 1 sec, continue
Vyt(istart:istop) = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*(tline(1:istop-istart+1) -phi +dphi) ); % start at 1 sec, then stop
%Vyt(istart:end) = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*(tline(1:end-istart+1) -phi+dphi)); %   & dphi in t-steps
 %Vyt = Vwm*sin(2*pi*ypt/lambda - 2*pi*wfreq*tline -phi); % all t;  now row vec

figure;  plot(Vyt); title('stimulus + waves + noise'); hold on
plot(Vyt, 'r-'); xlabel('t (msec) ')
 %plot(tline,Vyt, 'r-','LineWidth',1)  % V osc at fixed y vs t )ie "x"), in 2D
 

 %% 3.5 Constant stim - use Pd; not filtered by S[v]
% explore stable pts and limit cycles.
fprintf('\n ++const stimulus: at 0.1 sec ')
tstart=40001; % for plots
StimLevel= 155
Pd=zeros(1,length(x));  tstart0= 1001; %tstart0= 20001; % of this stim; cf separate for Pd
Pd(tstart0:end)=StimLevel; % && turn on the stim!!

%Ip6=zeros(1,length(x)); % pulse train/fn for 2nd input (?dphi shifted) - node 6

clear StimLevel PulseHeight tstart0

%% 3.6 Load 1st Pulse density function[s]
 % a la JR('93 & '95); set up in code WaveTrain.m @ #1.4
 fprintf(' Pulse shape stim +155(@0.5s); & 0.1mV ran noise: \n')
 Pd= 1.0*Pulse1; % 1 or 2 pulses;  vary amplitude (ie. pulse density
 tstart=40001;  % cft wave... code, for pulse setup
  % Pd= 1.0*PulseN; % N pulses, in train
 %Pd = Pd+100.0;  % add const stim (all t)
 Pd((1001):end) = Pd((1001):end)+155.0;  % add const stim , with pulse (not before)
  PulseHeight= max(Pd)
  %IpOn=zeros(nn,1); % switch Pulse stim on/off to ea node
  %IpOn(2)=0; % 1st input: stimulate Node 2 only (A32V: strongest input links)
  %IpOn(6)=1; %0.33; % 2nd input, to node 6 (A11), lower wts, 
    Ip6=zeros(1,length(x)); % pulse train/fn for 2nd input (?dphi shifted), to node 6
 Pd2=zeros(1, length(x)); % unless set & used, below
  % pran = 0.1*rand(1,length(x)); % low noise;  Uniform ran noise:
  clear PulseHeight 

 %% 3.6.1 Load 2nd Pulse density function[s]
 fprintf(' 2nd Pulse shape stim: \n')
 Pd2= 1.0*Pulse2; % 1 or 2 pulses;  vary amplitude (ie. pulse density
 Pd2 = Pd2+155.0;   % add const stim
 
%% 3.7 Input channels
 % selected nodes to get the stimulus inputs: Ip, Pd 
 % IpOn=ones(nn,1); % 1st stimulus On for all nodes
  IpOn=zeros(nn,1); % switch 1st stim on/off to ea node
 % stimulate Node 2 (A-32V): strongest, more wt, fewer links)
 % IpOn(2)=1; % 1st input: more wt; fewer #Links*0.61; % 1st or 2nd input?
  IpOn(6)=1;  %0.64; % 2nd input?, wt*0.64, to node 6 (A-11), lower wts, more Links  
 
 % stimulate Out Hub [1]A-10
 % IpOn(1)=1; % 0.4;  % A-10 <- TPO, TE3,2,11, Aud (weak), Opt, TEO, V4

 % stimulate other nodes:
 %IpOn(5)=1;  %0.8 % [5] A46D <- MO-pre
 %IpOn(4)=1;  %0.8 % [4] A9 <- Cing/RSP  
 IpOn' % check
  % direct 2nd pulse, if present?
 IpOnB=zeros(nn,1); % switch 1st stim on/off to ea node
 %IpOnB(2)=1; % 2nd Pulse to 1st input:
  %IpOnB(6)=1;  % 2nd Pulse to 2nd input (A11)
 IpOnB' % check
 
 %% 3.8 Gather stimuli & Plot, to check stim
 %pran= [pran100, 0] % need N=1 elements (5001)
 %Vyt=Vytburst;  % to pass to DE solver
 
 figure;  plot(pran, 'Color', [0.9 0.9 0.9]);  hold on
 plot(Pd+Pd2 + Ip, 'b'); title('pulse + Vstimulus + waves + noise');
  plot(Vyt, 'r-'); xlabel('t steps (0.1 msec) ')
  ylim([0 max(Pd)+20])
    if max(Pd+Pd2) <= 3
     ylim([-1 3])  %ylim([0 1])  
   else
     ylim([-1 max(Pd +Pd2)+20])
    end
   
% >>>>>>>>>>>>>> 
%% 4.0  RK4 DE solver, 1st order, for n variables
 % nb. matlab has default double precision: needed here
   % check timulus
   %figure;  plot(Ip+pran); title('stimulus + waves + noise'); hold on; plot(Vyt, 'r-'); xlabel('t (msec) ')
% reset arrays
y=zeros(nn*6,length(x)); % set up array for solutions; 2,3, etc nodes (6 rows each)
 %y=double(y); % enforce double precision (64 bit) : should be default
deltay=zeros(nn,length(x)); % set up array for y3-y2; 6 nodes (1 rows ea)
 %deltay=double(deltay);
y(:,1)= yic ; % y's are row vec, one for ea DE; IC is a col vec, for i=1 or t=0.   
In12=0;
  
% RK4 DE solution step:  iterate forward in time (here x); scan nodes(in) & linked NN (j)
  % fn JR1s incl forcing fn. as arguments, at this t-step; 
  % inputs: pran + p(i,j)+ [NM-NM couplings] {= sum(y1-y2) ; ext stim: Ip: to node 1 only
  %NOT NOW  onerow=ones(1,length(x)); % needs vector of 1's for sigmoid (over time span)
 %fprintf('\n   i  j wt(ij) dr(ij) index r Cfac ') % debug- header
tic
        % start "late" to allow for delay back ref, to L [allow 5ms w 0.5 steps]
for i = 50:(length(x)-1) % scan time steps & start at 1, 10, 100 [for delay DE]
    for in =1:nn % scan nodes in cluster (NM)
         j11=(in-1)*6+1; j16=(in-1)*6+6;  % working on this node #in
         r=rvec(in); Cfac=CfacVec(in); % for indiv NM : slope of S[v]; Gain C1,3
        %fprintf(' ... working on node  %4.0f \n', in)
    % Node #in - gets the ext inputs (Ip, Ii, p and In12 from linked nodes (i <- j)
     % and Vwave adds to In12 (a la lfp, at soma)
       jlist=find(Adj(:,in));  % linked nodes in (col) to node #in
       jnn=length(jlist); % #NN of 1
       for j=1:jnn % scan NN list of node #in
           if isempty(jlist)    
               continue % no NN, In12 stays at 0
           else
           jj=jlist(j); j11=(jj-1)*6+1; j16=(jj-1)*6+6; % for nodes 2, 3 ..
             %[in jj j11+1 j11+2]  % debug: check correct indices for delta-y inputs
             % will use: ((j-1)*6) +2 & ((j-1)*6 +3)
            % signal delay assumes v = 1 ms/ (mm/ms)
           vlocal=1.0; % unmyeln axon vel ~ signal delay ~ 0.1-1.0 m/s 
           %deltai=int16(Ad6A(jj,in)/(1000*h*vlocal)); % adjust index, for signal delay time from nn, jj->in / as integer
          deltai=int32(Ad6A(jj,in)/(1000*h*vlocal)); % check discontinuities
           %deltai=Ad6A(jj,in)/(1000*h); % Old: adjust index, for signal delay time from nn, jj->in
           In12= R12*Adj(jj,in)*vm./(1 +exp(r*(v0-y(j11+1,i-deltai)+y(j11+2,i-deltai)))); % sigmoid filter coupling, vector
           In12=In12+Vyt(i); % add travelling/standing wave here, now
    %fprintf('\n  %2.0f %2.0f %5.3f %4.1f %4.0f %4.2f %4.2f ', in, jj, Adj(jj,in), Ad6A(jj,in), (i-deltai), r, Cfac )
           end % check for NN  
       end % scan linked NNs; gather inputs (in <- jj) from linked NNs
           %fprintf('\n  ') % for debug
           %In21=zeros(12,length(x)); % debug, no input from other node[s]
           %[jj j11 j16]  % debug:
        j11=(in-1)*6+1; j16=(in-1)*6+6;  % update node #in
          % nb. y7:12 in Node-2; y13:18 is Node-3, etc.
        % Stimulus pulses to selected node (eg 2, 6), via IpOn switch
        % add travelling wave Voltage, at y=16mm (Ant), at this time(i), to _all_ NM
        %Istim=Vyt(i)+IpOn(in)*Ip(i); %pstim=IpOn(in)*pran(i); % turn stim on only for selected nodes - set above, at IC.
         % gather ran & wave train together: to all NM
        Istim=IpOn(in)*Ip(i); pstim=pran(i) +Vyt(i); % Ip to selected node; Vyt & ran noise to ALL nodes
         % apply S(v) to pran, Vyt & Ip in fn J12a; 
         pstim=vm/(1+exp(r*(v0-pstim))); %pstim=vm/(1+exp(r*(v0-pran(i)))); 
         Pdstim=IpOn(in)*Pd(i); Istim=vm/(1+exp(r*(v0-Istim)));  % switch 1st stim pulse on; filter Istim
         Pdstim=Pdstim+IpOnB(in)*Pd2(i);% switch 2nd stim pulse on, if present?;         
        if in == 6
             Istim=Vyt(i) +IpOn(in)*Ip6(i); pstim=pran(i);  % stim stream Ip to node #6; 
             Istim=vm/(1+exp(r*(v0-Istim))); %pstim=vm/(1+exp(r*(v0-pran(i))));  % V-> rate / Pd already set
        end
        k_1 = JR12a( x(i), y(j11:j16,i), Pdstim, Istim, pstim, In12, a, b, A, B, C, Cfac, r);  % these k's should be col vec, one for ea DE
        k_2 = JR12a( x(i)+0.5*h, y(j11:j16,i)+0.5*h*k_1, Pdstim, Istim, pstim, In12, a, b, A, B, C, Cfac, r); % seem to need to force col vec ?
        k_3 = JR12a( (x(i) +0.5*h), (y(j11:j16, i) +0.5*h*k_2), Pdstim, Istim, pstim, In12, a, b, A, B, C, Cfac, r);
        k_4 = JR12a( (x(i)+h), (y(j11:j16,i) +k_3*h), Pdstim, Istim, pstim, In12, a, b, A, B, C, Cfac, r);
        y(j11:j16,i+1) = y(j11:j16,i) + (k_1 +2*k_2 +2*k_3 +k_4)*(h/6); % load y1:y6 at this t-step      
    end %  scan nodes
end  % scan t-steps
 %fprintf('\n  ') % for debug
toc %timing
 % clear i in deltai In12 j*  Istim pstim
 
%% 5.0 PLOT OUTPUT 
tstart=4001; % gets lost sometimes ??
% Node 1 output [this is centre of star cluster, and Out-Hub (A-10)
  %figure; plot(y(1:3,:)'); title('JR/ node-1: y0, y1, y2')  % nb. need cols
  %title('A-10:  JR/ node-1, outputs:'); 
  %legend('y0', 'y1', 'y2'); % orig notation [not array index]
 figure; plot(y(2,:)'-y(3,:)'); 
 title('A-10:  JR/ node-1, output: detla-y: y1(1)- y2(1) '); 
    % skip transients
    %figure; plot(y(2,(tstart-100):end) - y(3,(tstart-100):end));  
     %title('A-10:  JR node 1: output, delta-y: y1(1) - y1(1)') 
     %text(3000, -1.5, 'stimulate #2, 6 '); text(3000, -2, ' 10Hz stimulus, x10 pulses ');
 figure; subplot(2,1,1); plot(x, y(2,:)'); hold on; plot(x, y(3,:)'); 
 title('JR/ node-1: output: y1, y2 '); legend('y1(1)', 'y2(1)'); hold off
 subplot(2,1,2); plot(x, y(2,:)'); hold on; plot(x, y(3,:)'); 
 axis([1.5 2.5 10 70]);  grid on; xlabel('t (s) ')
 figure; plot(x, (y(2,:)-y(3,:)) ); title('6-NM, node-1, delta-y '); xlabel('t (sec) ')
  
% Node 2 outputs [this has the primary inputs]:
 %figure; plot(y(8,(tstart-100):end) - y(9,(tstart-100):end));  title('JR node-2: output, y2(2) - y1(2)')
  %figure; plot(y(7:12,:)'); title('JR/ node-2: output, y1, y2, .. y6')  % nb. need cols
  %legend('y0', 'y1', 'y2', 'y3', 'y4', 'y5'); % 
 figure; subplot(2,1,1); plot(x,y(8,:)'); hold on; plot(x,y(9,:)'); 
 title('JR/ node-2: output: y1, y2 '); legend('y1(2)', 'y2(2)'); hold off
 subplot(2,1,2); plot(x,y(8,:)'); hold on; plot(x,y(9,:)'); 
 axis([1.5 2.5 20 100]);  grid on; xlabel('t (s) ')
     % other outputs
     %figure; plot(y(6,(tstart-100):end));  title('JR model, Node-1 output, y5  ')
     %figure; plot(y(12,(tstart-100):end));  title('JR model, Node-2 output, y5  ')
     % axis([0 5 -500 0]);  axis([0 20  -10000 10000]);
     figure;  plot(y(8,(tstart-800):end) - y(9,(tstart-800):end)); title('Node 2, delta-y ')
     xlabel('t (x0.1 ms) ')
 
 % Node 3 outputs:
 %figure; plot(y(14,(tstart-800):end) - y(15,(tstart-800):end));
 %title('JR model: node-3 output, y2(3) - y1(3)')   
% Node 4 outputs:
 %figure; plot(y(20,(tstart-800):end) - y(21,(tstart-800):end));
  %title('JR model: node-4 output, y2(4) - y1(4)') 
% Node 5 outputs:
 %figure; plot(y(26,(tstart-800):end) - y(27,(tstart-800):end));
 %title('JR model: node-5 output, y2(5) - y1(5)')
 %figure; subplot(2,1,1); plot(x,y(26,:)'); hold on; plot(x,y(27,:)'); 
 %title('JR/ node-5: output: y1, y2 '); legend('y1(5)', 'y2(5)'); hold off
 %subplot(2,1,2); plot(x,y(26,:)'); hold on; plot(x,y(27,:)'); 
 %axis([1.5 2.5 0 50]);  grid on; xlabel('t (s) ') 

% Node 6 outputs:
 tstart= 40001; % debug
 figure; plot(y(32,(tstart-100):end) - y(33,(tstart-100):end));
  title('JR model: node-6 output, y1(6) - y2(6)')
 figure; plot(y(32,:)'); hold on; plot(y(33,:)'); 
  title('JR/ node-6: output: y1, y2 '); legend('y1(6)', 'y2(6)');
 figure; subplot(2,1,1); plot(x,y(32,:)'); hold on; plot(x,y(33,:)'); 
  title('JR/ node-6: output: y1, y2 '); legend('y1(6)', 'y2(6)'); hold off
  subplot(2,1,2); plot(x,y(32,:)'); hold on; plot(x,y(33,:)');
  plot(x, Pd/10) % stim
  xlim([3.9 4.4]); grid on; xlabel('t (s) '); legend('ye(6)', 'yi(6)', 'Stim') 
 
% all 6 nodes together:
tstart= 40001; % debug
 figure; 
 plot(x((tstart-800):end), y(2,(tstart-800):end) - y(3,(tstart-800):end), 'k', ...
     'LineWidth', 0.75);  hold on % node #1, A-10 
 xlim([4 4.2]); xlabel('t (s)'); xlabel('t (ms)');  %pause
 plot(x((tstart-800):end) ,(y(8,(tstart-800):end) - y(9,(tstart-800):end) ), ...
     'LineWidth', 0.75);   % #2, A-32V
 %pause  % - to closely examine waveforms
 plot(x((tstart-800):end) ,y(14,(tstart-800):end) - y(15,(tstart-800):end), ...
     'LineWidth', 0.75); % #3 A-32
 %pause
 plot(x((tstart-800):end) ,y(20,(tstart-800):end) - y(21,(tstart-800):end), ...
     'LineWidth', 0.75); % #4, A-9
 %pause
 plot(x((tstart-800):end) ,y(26,(tstart-800):end) -   y(27,(tstart-800):end), ...
     'LineWidth', 0.75); % #5, A-46D
 %pause  
 plot(x((tstart-800):end) ,y(32,(tstart-800):end) - y(33,(tstart-800):end), ...
     'LineWidth', 0.75); % #6, A-11
 xlabel('t (s)')
 plot(x((tstart-800):end) ,Pd((tstart-800):end)+pran((tstart-800):end), 'Color', [0.9 0.9 0.9]);  
 xlabel('t (s)'); ylabel('dy(t) (mV)')
 %plot(x((tstart-800):end) ,Ip6((tstart-800):end), 'g--') % stimulii, on same scale
 %plot(x((tstart-800):end) ,Vyt((tstart-800):end), 'r--') % stimulii: wave
  legend('dy(1)', 'dy(2)', 'dy(3)', 'dy(4)', ...
      'dy(5)','dy(6)', 'stimulus' ); %text(3000, 0.3, 'wt(1-2) = 100')
  % axis([0 5000 -20 20]);  %axis([1 2 -25 25]); 
%
% mean and range of detla-y, overall: 6 nodes together:  cf. #4.3 below, for calc at 17:19sec 
 d_yMean=zeros(nn,1); d_yAmpl=zeros(nn,1); 
 d_yMean(1)=mean(y(2,(tstart+2000):end) - y(3,(tstart+2000):end));   % node #1, A-10 
 d_yAmpl(1)=max(y(2,(tstart+2000):end) - y(3,(tstart+2000):end)) ... % sample after stim settles
     - min(y(2,(tstart+2000):end) - y(3,(tstart+2000):end));
 d_yMean(2)=mean(y(8,(tstart+2000):end) - y(9,(tstart+2000):end));   % #2, A-32V
 d_yAmpl(2)=max(y(8,(tstart+2000):end) - y(9,(tstart+2000):end)) ...
     - min (y(8,(tstart+2000):end) - y(9,(tstart+2000):end));
 d_yMean(3)=mean(y(14,(tstart+2000):end) - y(15,(tstart+2000):end)); % #3 A-32
 d_yAmpl(3)=max(y(14,(tstart+2000):end) - y(15,(tstart+2000):end)) ...
      -min(y(14,(tstart+2000):end) - y(15,(tstart+2000):end));
 d_yMean(4)=mean(y(20,(tstart+2000):end) - y(21,(tstart+2000):end)); % #4, A-9
 d_yAmpl(4)=max(y(20,(tstart+2000):end) - y(21,(tstart+2000):end)) ...
      - min(y(20,(tstart+2000):end) - y(21,(tstart+2000):end));
 d_yMean(5)=mean(y(26,(tstart+2000):end) - y(27,(tstart+2000):end)); % #5, A-46D
 d_yAmpl(5)=max(y(26,(tstart+2000):end) - y(27,(tstart+2000):end)) ...
      - min(y(26,(tstart+2000):end) - y(27,(tstart+2000):end));
 d_yMean(6)=mean(y(32,(tstart+2000):end) - y(33,(tstart+2000):end)); % #6, A-11
 d_yAmpl(6)=max(y(32,(tstart+2000):end) - y(33,(tstart+2000):end)) ...
      -min(y(32,(tstart+2000):end) - y(33,(tstart+2000):end));
 figure; stem(d_yMean, '+');  hold on; stem(d_yAmpl, ':diamond');
  % hs(1).Marker='+'; hs(2).Marker='diamond'; 
 axis([0 7 -15 60 ]); xlabel('NM #'); ylabel('delta-y12(i)  (mV)'); title('dy(i) overall')
  legend('mean', 'p-p ampl.', 'Location', 'northeast') % place out of the way
   % text(5.5, 23, 'Adj = 1');  % text(5.5, 10, 'alpha; C x vector')
 
   [d_yMean d_yAmpl] % test/ debug
%   %figure; plot(rvec, d_yAmpl, '.', 'MarkerSize', 20); hold on
   % plot(rvec, d_yMean, '+', 'MarkerSize', 5);
 
% 3 key nodes together:
 figure;   hold on 
 plot(y(8,(tstart-800):end) - y(9,(tstart-800):end), 'r') % ', 'Color', [0 0.3 0]);  % #2, A-32V : Input 1st
 plot(y(32,(tstart-800):end) - y(33,(tstart-800):end),'g'); % #6, A-11 : Input 2nd
 plot(Pd((tstart-800):end), 'k');  %plot(Ip6((tstart-800):end), 'k--') % stimulii, on same scale
 plot(y(2,(tstart-800):end) - y(3,(tstart-800):end),'b'); % #1, A-10 : Output
  title('6-star cluster, JR model: In/Out nodes delta-y ')
  legend('delta-y(1): A-10', 'delta-y(2)', ...
         'delta-y(6)', 'stimulus' );
%   
% LFP LFP:  net lfp output: sum: detla-y(1) + detla-y(2), for the 6 NM
 lfp= (y(2,:) - y(3,:) + y(8,:) - y(9,:) +y(14,:) - y(15,:) ...
     + y(20,:) - y(21,:)+ y(26,:) - y(27,:) +y(32,:) - y(33,:) )/6; % av of 6 nodes
% figure; plot(lfp); title('JR model: net lfp output ')
  figure; plot(x, lfp); title('6-star cluster, JR model: net lfp output '); xlabel('t (sec) ')
  ylabel('LFT (t)  (mV)'); % hold on; plot(20,-1, 'b^'); plot(20,-1.1, 'b^'); % wave starts
   % text(15, -3.5, 'wave ampl = 10 mV') % plot(4,-4, 'k^'); % Pulse starts
    % text(2.0, -120.0, '6-cluster: ran noise input ') 
    % text(4000, -2.0, 'stimulate #2, 6 ')
    % axis([500 5000 -0.1 0.1])  % axis([500 5000 -2.0 0])
   % text(2.5, -4, '16 Hz wave + 100 Hz pulses to A32V')
%   
% Group of 4 close together: net lfp output: sum: detla-y(2:5) 
 lfp25= ( y(8,:) - y(9,:) +y(14,:) - y(15,:) ...
     + y(20,:) - y(21,:)+ y(26,:) - y(27,:) )/4; % av of 4 nodes
  %figure; plot(lfp); title('JR model: net lfp output ')
%figure; plot(x, lfp25); title('4NM of 6, local lfp output '); xlabel('t (sec) ')
  % text(2.0, -120.0, '4NM of 6: ran noise input ') 
  % text(4000, -2.0, 'stimulate #2, 6 ')
         % axis([500 5000 -0.1 0.1])  % axis([500 5000 -2.0 0])
   % text(2.5, -4, '16 Hz wave + 100 Hz pulses to A32V')
  lfp_stats_2to5= [mean(lfp25(tstart:end))  (max(lfp25(tstart:end)) - min(lfp25(tstart:end)) )]
  clear lfp25  lfp_stats_2to5
  
% Grouped summary: resize figure
  %figwidth = 1024; figheight = 896; figposition = [100, 100, figwidth, figheight]; % large
  figwidth = 560; figheight = 704; figposition = [500, 100, figwidth, figheight]; % tall
 figure('position',figposition, 'units','pixels');  %figure; % default
 subplot(5,1,1);
 plot(y(2,:) - y(3,:),'b'); % #1, A-10 : Output
 tlim=length(y(2,:) - y(3,:))-1; axis([0,tlim, -Inf, Inf]) % keep vetical as is
  title('6-star cluster, NM#1, (JR model: A-10: Out node delta-y ')
 subplot(5,1,2);
 plot(y(8,:) - y(9,:), 'r') % ', 'Color', [0 0.3 0]);  % #2, A-32V : Input 1st
  title('     NM#2, A-32V: 1st In node delta-y '); axis([0,tlim, -Inf, Inf])
 subplot(5,1,3);
 plot(y(32,:) - y(33,:),'g'); % #6, A-11 : Input 2nd
  title('     NM#6, A-11: 2nd In node delta-y '); axis([0,tlim, -Inf, Inf])
  subplot(5,1,4);  plot(lfp,'k'); title('     lfp output ');  %  net output
  xlim([0,tlim]); % xlim([0 3.0e5]); %axis([0,tlim, -Inf, Inf])
 subplot(5,1,5); 
 plot((Pd)*IpOn(2) +pran+Vyt, 'g--'); hold on; 
  %plot((Pd2)*IpOnB(2), 'g--'); % incl 2nd pulse fn
  plot((Pd)*IpOn(6) +pran+Vyt, 'r--') % stimulii, on same scale
  %plot((Pd2)*IpOnB(6), 'r--'); xlim([0,tlim]); %xlim([0 3.0e5])
  % plot((Ip+Pd)*IpOn(1) +pran+Vyt, 'b--') % if needed?
  legend( 'stimulus-2' , 'stimulus-6' ); title('     Inputs: stimulus->2, ->6 ');  
   xlabel('t steps (0.1 msec) '); %axis([0,tlim, -Inf, Inf])
    % text(5500, 50, 'Stim 300 to both #2, #6')
%
 lfp_sample=lfp(40001:50000); % during stimulus [0.5s sample]
   %lfp_sample=lfp(1100:2900); % debug, nb. wider sample for half t-step
 lfp_stats_atStim= [ mean(lfp_sample), (max(lfp_sample)-min(lfp_sample))]
 % figure; plot(lfp_sample - mean(lfp_sample)); title('JR, 6xNM: sum (lfp -mean), during stim. ')

 lfp_sample=lfp(60001:70000); % after stim, in s/s
  %lfp_sample=lfp(5000:6000); % debug, nb. wider sample for half t-step
 lfp_stats_atSS= [ mean(lfp_sample), (max(lfp_sample)-min(lfp_sample))]
 
 lfp_sample=lfp(70001:80000);   % much later
 lfp_stats_atLaterSS= [ mean(lfp_sample), (max(lfp_sample)-min(lfp_sample))]
 
 lfp_sample=lfp(80001:90000);   % much later
 lfp_stats_at8to9secSS= [ mean(lfp_sample), (max(lfp_sample)-min(lfp_sample))]
% 
  if xf >15 % for longer sim (eg 20 sec)
   lfp_sample=lfp(180000:190000);   % V. much later
   lfp_stats_at18secLaterStill= [ mean(lfp_sample), (max(lfp_sample)-min(lfp_sample))]
  end
  if xf >25 % for v longer sim (eg 30 sec)
   lfp_sample=lfp(280000:290000);   % V. much later
   lfp_stats_at18secLaterStill= [ mean(lfp_sample), (max(lfp_sample)-min(lfp_sample))]
  end
  
 clear jlist i j in jj j11 j16 k*  In12 deltay lfp_* d_* dy6* % Ip pran yic

% > > > > >

%% joint plot of dy-i & lfp
figure; subplot(2,1,1); plot(x, lfp); title('6-NM cluster, JR model: net lfp output ');
xlabel('t (sec) ');  ylabel('LFT (t)  (mV)');
% all 6 nodes together:
tstart= 40001; % debug
 subplot(2,1,2); plot(x((tstart-800):end), y(2,(tstart-800):end) - y(3,(tstart-800):end), 'k', ...
     'LineWidth', 0.75);  hold on % node #1, A-10 
 xlim([4.0 5.0]); xlabel('t (s)'); xlabel('t (ms)'); 
 plot(x((tstart-800):end) ,(y(8,(tstart-800):end) - y(9,(tstart-800):end) ), ...
     'LineWidth', 0.75);   % #2,
 plot(x((tstart-800):end) ,y(14,(tstart-800):end) - y(15,(tstart-800):end), ...
     'LineWidth', 0.75); % #3
 plot(x((tstart-800):end) ,y(20,(tstart-800):end) - y(21,(tstart-800):end), 'b', ...
     'LineWidth', 0.75); % #4, A9
 plot(x((tstart-800):end) ,y(26,(tstart-800):end) -   y(27,(tstart-800):end), ...
     'LineWidth', 0.75); % #5, A-46D  
 plot(x((tstart-800):end) ,y(32,(tstart-800):end) - y(33,(tstart-800):end), ...
     'LineWidth', 0.75); % #6, A-11 
 plot(x((tstart-800):end) ,Pd((tstart-800):end)+pran((tstart-800):end), 'Color', [0.8 0.8 0.8]);  
  %plot(x((tstart-800):end) ,Vyt((tstart-800):end), 'r--') % stimulii: wave
  legend('dy(1)', 'dy(2)', 'dy(3)', 'dy(4)', 'dy(5)','dy(6)', 'stimulus' ); 

  
%% 4.1 Align waveforms, by phase: 1st Pulse
% a) at NM#2 [the input]
dy2_sample=y(9, 40000:60000) - y(8, 40000:60000); % NM#2, delta-y 
figure; plot(dy2_sample); hold on; xlim([0 2000]) %axis([0.5 2.5 -50 50])
dym=mean(dy2_sample); plot([0.5 2.5], [dym dym], 'g--'); xlabel('t (ms)')
%
title('JR node #2, 1t pulse: output delta-y: y1(2) - y2(2)') % e - i 
plot(Pd(40000:60000)/10, 'k') % stim % reduce scale [of pulse rate]
plot(Ip(40000:60000), 'k') ; plot(Vyt(40000:60000)+dym, 'Color', [1 0 1]) % stim V pulses, wave
legend('delta-y(2)', 'Pulse', 'wave'); 

% a.1)
dy6_sample=y(32, 40000:60000) - y(33, 40000:60000); % NM#6, delta-y 
figure; plot(dy6_sample); hold on; xlim([0 2000]) %axis([0.5 2.5 -50 50])
dym=mean(dy6_sample); plot([0.5 2.5], [dym dym], 'g--'); xlabel('t (ms)')
title('JR node #6: output delta-y: y1(6) - y2(6)') % e - i 
plot(Pd(40000:60000)/10, 'k') % stim % reduce scale [of pulse rate]
plot(Ip(40000:60000), 'k') ; plot(Vyt(40000:60000)+dym, 'Color', [1 0 1]) % stim V pulses, wave
legend('delta-y(6)', 'Pulse', 'wave')

% b) lfp [the output]
%lfp_sample=lfp(1:1200); % all NM#2, nett lfp
lfp_sample=lfp(40000:60000);
figure; plot( lfp_sample); hold on; xlim([0 2000]) %axis([600 1200 -5 45])
dym=mean(lfp_sample); plot([3.500 6.0], [dym dym], 'g--')
title('JR 6 NM: align lfp w stim') 
plot(Pd(40000:60000)/10, 'k') % stim % reduce scale
plot(Ip(40000:60000), 'k');  %plot(Ip(1:1200), 'k'); % 
plot(Vyt(40000:60000)+dym, 'Color', [1 0 1]) % wave; align axes
legend('lfp', 'Pulse', 'wave'); %axis([0.5 2.5 -10 30]); 
xlabel('t steps (0.1ms)')

% c) align Pulse with y-e(2, 6)
dym=min(y(8, 40000:60000));
figure; subplot(2,1,1); plot(y(8, 40000:60000)); hold on
plot(dym+Pd(40000:60000)/10, 'k'); legend('y-e(2)', 'Stim');
xlim([0 4000]); grid on; title('Align Pulse (@ 4s) with y-e and y-i, NM #2, #6 ')
dym=min(y(32, 40000:60000));
subplot(2,1,2); plot(y(32, 40000:60000)); hold on; grid on
plot(dym+Pd(40000:60000)/10, 'k'); legend('y-e(6)', 'Stim')
xlim([0 4000]); xlabel('t steps (0.1ms)')

 clear dy2_* dy6_ dym 
%% 4.1.1 Align waveforms, and 2nd Pulse by phase - at later times (20s)
% a) at NM#2 [the input]
dy2_sample=y(9, 190000:210000) - y(8, 190000:210000); % NM#2, delta-y 
figure; plot(dy2_sample); hold on; xlim([9000 12000]) %axis([0.5 2.5 -50 50])
dym=mean(dy2_sample); plot([19.0 21.0], [dym dym], 'g--'); xlabel('t (ms)')
%
title('JR node #2, 2nd pulse: output delta-y: y1(2) - y2(2)') % e - i 
plot(Pd2(190000:210000)/10, 'k') % stim % reduce scale [of pulse rate]
plot(Ip(190000:210000), 'k') ; plot(Vyt(190000:210000)+dym, 'Color', [1 0 1]) % stim V pulses, wave
legend('delta-y(2)', 'Pulse', 'wave'); 

% a.1)
dy6_sample=y(32, 190000:210000) - y(33, 190000:210000); % NM#6, delta-y 
figure; plot(dy6_sample); hold on; xlim([9000 12000]) %axis([0.5 2.5 -50 50])
dym=mean(dy6_sample); plot([19.0 21.0], [dym dym], 'g--'); xlabel('t (0.1ms)')
title('JR node #6, 2nd pulse: output delta-y: y1(6) - y2(6)') % e - i 
plot(Pd2(190000:210000)/10, 'k') % stim % reduce scale [of pulse rate]
plot(Ip(190000:210000), 'k') ; plot(Vyt(190000:210000)+dym, 'Color', [1 0 1]) % stim V pulses, wave
legend('delta-y(6)', 'Pulse', 'wave')

% b) lfp [the output]
%lfp_sample=lfp(1:1200); % all NM#2, nett lfp
lfp_sample=lfp(190000:210000);
figure; plot( lfp_sample); hold on; xlim([9000 12000]) %axis([600 1200 -5 45])
dym=mean(lfp_sample); plot([3.500 6.0], [dym dym], 'g--')
title('JR 6 NM LFP, 2nd pulse: align w stim') 
plot(Pd2(190000:210000)/10, 'k') % stim % reduce scale
plot(Ip(190000:210000), 'k');  %plot(Ip(1:1200), 'k'); % 
plot(Vyt(190000:210000)+dym, 'Color', [1 0 1]) % wave; align axes
legend('lfp', 'Pulse', 'wave'); %axis([0.5 2.5 -10 30]); 
xlabel('t steps (0.1ms)')

 clear dy2_* dy6_ dym 
 
 %% 4.2 Phase plot (3 y's for  NM #1)

figure;  subplot(3,1,1); plot(y(1,:), y(4,:)); 
 %title('dy0/dt vs y0(t);   zero input, C= 135, tau(15, 20), 1e-12 ran noise'); hold on
title('dy0/dt vs y0(t); C= 135, tau(15, 20),  0.1mV noise, stim = 300 '); hold on

ylabel('dy0/dt' ); xlabel('y0(t)')
 subplot(3,1,2); plot(y(2,:), y(5,:)); title ('dy1/dt vs y1(t)')
   ylabel('dy1/dt' ); xlabel('y1(t)')
 subplot(3,1,3); plot(y(3,:), y(6,:)); title ('dy2/dt vs y2(t)')
   ylabel('dy2/dt' ); xlabel('y2(t)'); hold off

%% 4.2.1 state space plot (3 y's for 1 NM)

figure; plot(y(2,:), y(3,:)); 
title('y1(t) vs y2(t); step stim, C= 250, tau(12, 16 ms), 10 mV G ran noise'); hold on
xlabel('y1(t) '); ylabel('y2(t) ');
figure; plot3(y(2,:), y(3,:), y(1,:)); 
title('y0(t) vs y1(t) vs y2(t); no stim, C= 300, tau(10, 10)  ');
xlabel('y1(t) '); ylabel('y2(t) '); zlabel('y0(t) '); grid on
 %text(40, 10, 0.02, 'single pulse, peak 195')
 %text(40, 10, 0.01, ' 0 inhib, Ii ')

  %% 4.2.2 state space plot (3 y's for 1 NM): at stim & at s/s
 % during stim; 100ms after transients
figwidth = 560; figheight = 688; figposition = [500, 100, figwidth, figheight]; % tall
figure('position',figposition, 'units','pixels');  %figure; % default
 subplot(2,1,1); plot3(y(2,100:1000), y(3,100:1000), y(1,100:1000)); hold on;
 title('6NM, before stim; y0-y1-y2')
fp1=mean(y(2,100:1000)); fp2=mean(y(3,100:1000)); fp3=mean(y(1,100:1000)); 
plot3(fp1, fp2, fp3, 'k.', 'MarkerSize', 20)
FixedPoint_stim=[fp1 fp2 fp3]
title('y0(t) vs y1(t) vs y2(t) at stim;  C= 135, tau(15, 20), 8Hz wave  ');
xlabel('y1(t) '); ylabel('y2(t) '); zlabel('y0(t) '); grid on
 %thisaxis=axis; %axis manual % freeze
% figure; 
subplot(2,1,2); plot3(y(2,2201:end), y(3,2201:end), y(1,2201:end)); hold on
plot3(fp1, fp2, fp3, 'k.', 'MarkerSize', 20) % replot marker
title('y0(t) vs y1(t) vs y2(t);  in steady state  '); 
xlabel('y1(t) '); ylabel('y2(t) '); zlabel('y0(t) '); grid on
fp1=mean(y(2,2201:end)); fp2=mean(y(3,2201:end)); fp3=mean(y(1,2201:end)); 
plot3(fp1, fp2, fp3, 'b.', 'MarkerSize', 25);
 %axis(thisaxis); % same scales
FixedPoint_ss=[fp1 fp2 fp3]
  clear fp* FixedPoint* thisaxis

%% 4.3   dy(i) at later times: 17-19s, just before wave arrives
% mean and range of detla-y, overall: 6 nodes together:
 d_yMean=zeros(nn,1); d_yAmpl=zeros(nn,1); 
 d_yMean(1)=mean(y(2,(170001):190000) - y(3,(170001):190000));   % node #1, A-10 
 d_yAmpl(1)=max(y(2,(170001):190000) - y(3,(170001):190000)) ... % sample after stim settles
     - min(y(2,(170001):190000) - y(3,(170001):190000));
 d_yMean(2)=mean(y(8,(170001):190000) - y(9,(170001):190000) );   % #2, A-32V
 d_yAmpl(2)=max(y(8,(170001):190000) - y(9,(170001):190000)) ...
     - min (y(8,(170001):190000) - y(9,(170001):190000));
 d_yMean(3)=mean(y(14,(170001):190000) - y(15,(170001):190000)); % #3 A-32
 d_yAmpl(3)=max(y(14,(170001):190000) - y(15,(170001):190000)) ...
      -min(y(14,(170001):190000) - y(15,(170001):190000));
 d_yMean(4)=mean(y(20,(170001):190000) - y(21,(170001):190000)); % #4, A-9
 d_yAmpl(4)=max(y(20,(170001):190000) - y(21,(170001):190000)) ...
      - min(y(20,(170001):190000) - y(21,(170001):190000));
 d_yMean(5)=mean(y(26,(170001):190000) - y(27,(170001):190000)); % #5, A-46D
 d_yAmpl(5)=max(y(26,(170001):190000) - y(27,(170001):190000)) ...
      - min(y(26,(170001):190000) - y(27,(170001):190000));
 d_yMean(6)=mean(y(32,(170001):190000) - y(33,(190000):190000)); % #6, A-11
 d_yAmpl(6)=max(y(32,(170001):190000) - y(33,(170001):190000)) ...
      -min(y(32,(170001):190000) - y(33,(190000):190000));
 figure; stem(d_yMean, '+');  hold on; stem(d_yAmpl, ':diamond');
  % hs(1).Marker='+'; hs(2).Marker='diamond'; 
 axis([0 7 -15 60 ]); xlabel('NM #'); ylabel('delta-y12(i)  (mV)'); title('dy(i) at 17-19s')
  legend('mean', 'p-p ampl.', 'Location', 'northeast') % place out of the way
   % text(5.5, 23, 'Adj = 1');  % text(5.5, 10, 'alpha; C x vector')
 
   [d_yMean d_yAmpl] % test/ debug
   clear d_*
%

 %% Appx B. Period; 
  % Lorenz map (later?)
% plot z(n+1) vs z(n)  % ie. find local max
% 1st find period - ie. zero crossings, to within dt
figure; plot(y(1,1000:1200)); hold on
title('y(1,t) zero crossings ')
ym=mean(y(1,1000:1200)); % "zero" line
yzeros=find( abs(y(1,1000:1200) - ym) < 1.5e-3) % need to guess tol ?
T= (yzeros(2)-yzeros(1))*2 % ie. 2 half waves
f = 1000/T % Hz
%
plot([0 250], [ym ym], 'g-'); 
for j=1:length(yzeros)
plot(yzeros(j), y(1,1000+yzeros(j)), 'sr');  
end

%clear ym yzeros T f

%% Appx B.1 Period of y1 ie y(0);  at s/s
% plot z(n+1) vs z(n)  % ie. find local max
% 1st find period - ie. zero crossings, to within dt
Inode=5
dysample= (y((Inode-1)*6+2,3000:3200) -y((Inode-1)*6+3,3000:3200) ); % delta-y: y2- y3 for this node; 201 pts
 % y(8,(tstart-800):end) - y(9,(tstart-800):end) % Node#2
figure; plot(dysample); hold on  % nb. (Inode-1)*6+1 is y0(Inode)
title('NM#6, delta-y(t) zero crossings ')
ym=mean(dysample) % "zero" line
plot([0 200], [ym ym], 'g-');  % mark mean ("zero") line
yzeros=find( abs(dysample - ym) < 2.0e-1) % need to guess tol ?
T= (yzeros(2)-yzeros(1))*2; % ie. 2 half waves
f = 1000/T; % Hz
PeriodAndFreq = [T f]
%
plot([0 250], [ym ym], 'g-'); 
for j=1:length(yzeros)
plot(yzeros(j), dysample(yzeros(j)), 'sr');
end
PeakToPeak = max(dysample) - min(dysample);
OscAt = [ym PeakToPeak]
% clear ym yzeros T f Inode dysample PeriodAndFreq PeakToPeak


 %% Appx 1. check Adj setup - use marmoset data directly
 % 6 nodes are:  2, A10; 3, A11; 26, A32V; 25, A32; 48, A9; A32, A46D
 clusterlist=[2, 3, 26, 25, 48, 32]; % in y-coord order
   % nb. log p66: uses x-y plot clockwise order: 2, 26, 25, 32, 3, 48. 
 
 Anew= csvread('AdjMarmoset_Rescale2b.csv');  % based on LNe
 
 Aclust = Anew(clusterlist, clusterlist);
 Aclust/1000  % as used for 6 NM model
 % enter by hand, rounded
 Adj1 = [0.0,    8.54,    9.1,    0.39,   0.84   4.56; ...
         0.42,   0.0,     1.46,   0.03,   0.44,  0.063; ...
         0.23,   0.596,   0.0,    0.325,  0.09,  0.003; ...
         0.53,   1.14,    4.88,   0.0,    2.2,   0.23; ....
         0.76,   1.28,    0.68,   1.19,   0.0,   1.64; ....                  
         0.68,   0.17,    0.005,  0.0,    0.114, 0.0]
     
length(find(Adj1(:))) % 29 entries
 %Aclust = Anew(clusterlist, clusterlist);
 
 %% Appx 1.1 [total] Deg
clusterlist=[2, 26, 25, 48, 32, 3]; % 6*
DegInTot = sum((Anew(:,clusterlist)))'

kOutTot=length(find(Anew(clusterlist,:)))
DegWtOut = sum((Anew(clusterlist,:)))

%% Appx 2. in orig order: cf. logbook, p58
  fprintf('  use marmoset LNe Adj(2):  \n')
clusterlist=[2, 26, 25, 48, 32, 3]; % as plotted in x-y plane, viewed from Ant.
Anew= csvread('AdjMarmoset_Rescale2b.csv');  % based on LNe 
 Aclust = Anew(clusterlist, clusterlist);
Adj2=Aclust/1000  % as used for 6 NM model
%Adj=Adj2;   % for main code
length(find(Adj2(:))) % 29 entries

 %clear Anew Aclust
 
 %% Appx 3. Links & wts In/Out nodes
 %  Adj from Marmoset (LNe) data
 % cf.Appx 2. in orig order: cf. logbook, p58
 fprintf('  use marmoset LNe Adj(2): net In links & wt \n')
clusterlist=[2, 26, 25, 48, 32, 3]; % 6* as plotted in x-y plane, viewed from Ant.
Anew= csvread('AdjMarmoset_Rescale2b.csv');  % based on LNe 
% #Links: un-wt-Deg, less int links ~ net ext links
 kIn_6=length(find(Anew(:,3))) - length(find(Anew(clusterlist,3)))  % within cluster
 kOut_2=length(find(Anew(26,:))) -length(find(Anew(26, clusterlist)))
 kOut_6=length(find(Anew(3,:))) - length(find(Anew(3, clusterlist)))
 
 kIn_1 =length(find(Anew(:,2))) -length(find(Anew(clusterlist,2))) 
 kOut_1=length(find(Anew(2,:))) - length(find(Anew(2, clusterlist)))
 
% wt-Deg
 DegIn_2=sum((Anew(:,26))) -sum((Anew(clusterlist,26))) % less int links
 DegIn_6=sum((Anew(:,3))) - sum((Anew(clusterlist,3)))  % within cluster
 DegOut_2=sum((Anew(26,:))) -sum((Anew(26, clusterlist)))
 DegOut_6=sum((Anew(3,:))) - sum((Anew(3, clusterlist)))
 % Out Hub
 DegIn_1 =sum((Anew(:,2))) -sum((Anew(clusterlist,2))) 
 DegOut_1=sum((Anew(2,:))) - sum((Anew(2, clusterlist)))
 %Other key Inputs
 kIn_5=length(find(Anew(:,32))) - length(find(Anew(clusterlist,32)))  % within cluster
 DegIn_5 =sum((Anew(:,32))) -sum((Anew(clusterlist,32)))
 
 kIn_4=length(find(Anew(:,48))) - length(find(Anew(clusterlist,48)))  % within cluster
 DegIn_4 =sum((Anew(:,48))) -sum((Anew(clusterlist,48)))
 
 kIn_3=length(find(Anew(:,25))) - length(find(Anew(clusterlist,25)))  % within cluster
 DegIn_3 =sum((Anew(:,25))) -sum((Anew(clusterlist,25)))
 
  clear k* Deg* Anew clusterlist
  
  %% Appx 3.1 Int Links only: kIn, wtDegIn
  % use Adj6x6 set in #1.1 above
  fprintf('  Adj internal to cluster: In / Out links & wt \n')
  kIn=zeros(nn,1); DegIn=zeros(nn,1); kOut=zeros(nn,1); DegOut=zeros(nn,1);
  fprintf('\n In: ')
  ratios=[];
  for i= 1:nn
     kIn(i) =length(find(Adj(:,i)));  DegIn(i) =sum((Adj(:,i)));
      %[kIn(i)  DegIn(i)] % debug
     fprintf('\n kIn, DegIn, wt/k: %3.0f  %6.3f  %6.3f', kIn(i),  DegIn(i), DegIn(i)/kIn(i) )
     ratios =[ratios, DegIn(i)/kIn(i)];
  end
  fprintf('\n')
  [mean(ratios) std(ratios)] % mean & std dev
  fprintf('\n Out: ')
  for i= 1:nn
      kOut(i) =length(find(Adj(i,:)));  DegOut(i) =sum((Adj(i,:)));
      %[kOut(i)  DegOut(i)]
     fprintf('\n kOut, DegOut, wt/k: %3.0f  %6.3f  %6.3f', kOut(i),  DegOut(i),  DegOut(i)/kOut(i))
  end
  fprintf('\n')
 % 
  fprintf('\n InWt/InLinks: ')
  for i= 1:nn
     %[(DegIn(i)/kIn(i)) + DegOut(i)/kOut(i) ]
      fprintf('\n kOut, DegOut, wt/k: %3.0f  %6.3f  %6.3f', kIn(i),  DegIn(i),  DegIn(i)/kIn(i))
  end
  fprintf('\n Total v & links, In+Out: ')
  for i= 1:nn
     [(DegOut(i)+DegIn(i))/(kIn(i)+ kOut(i))]
  end
  
  clear kIn DegIn kOut DegOut
   
%% Appx 3.2 Link distances of 6 nodes
 %  d(i,j)
 addpath(genpath('/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain')); % Data files;  include sub-directories 
fprintf('  use marmoset d(i,j) from Atlas \n')
clusterlist=[2, 26, 25, 48, 32, 3]; % as plotted in x-y plane, viewed from Ant.
Adist= csvread('MarmosetDistPairs.csv');  % based on Atlas; saved file 
% csvwrite( 'MarmosetDistPairs.csv',Adist); % 
%
figure; hist(Adist(:), 100); hold on  % ~ skewedGaussian, peaks at 7, 11, 13 mm
Adist6 =Adist(clusterlist,clusterlist);
figure; hist(Adist6(:), 50); title('Marmoset: 6 NM cluster d(i,j)');  xlabel('d(i,j)  (mm) ')
clear Adist 

%% Appx 3.2.1  set d in bin of 0.5 mm (approx)
Ad6A = [...
     0     3.5   3.0   2.5   2.5   2.5; ...
     3.5   0     1.5   3.5   4.0   2.5; ...
     3.0   1.5   0     2.0   3.0   2.5; ...
     2.5   3.5   2.0   0     2.0   3.5; ...
     2.5   4.0   0     2.0   0     3.0; ...
     2.5   2.5   2.5   3.5   3.0   0];
 
%% Appx 3.2.2  set d in bin of 0.1 mm (approx)
Ad6A = [...
     0     3.4   2.9   2.4   2.4   2.5; ...
     3.4   0     1.7   3.4   4.1   2.7; ...
     2.9   1.7   0     1.9   2.8   2.6; ...
     2.4   3.4   1.9   0     1.8   3.5; ...
     2.4   4.1   0     2.0   0     3.4; ...
     2.5   2.7   2.6   3.4   2.8   0];
 
%% Appx 4. 8-9-10 NM Adj setup - use marmoset data directly
 % 2+6 nodes are:  2, A10; 3, A11; 26, A32V; 25, A32; 48, A9; A32, A46D,
 % + A8aD, A8b; then maybe later: A47L, A6DR
     % cf Logbook, p66 ; workbook#3, p17
 clusterlist=[2, 3, 26, 25, 48, 32, 45, 47]; % in y-coord order (A-P)
   % nb. log p66: uses x-y plot clockwise order: 2, 26, 25, 32, 3, 48. 
 Anew= csvread('AdjMarmoset_Rescale2b.csv');  % based on LNe
 
 Aclust = Anew(clusterlist, clusterlist);
 Aclust/1000  % as used for 6/8 NM model
 length(find(Aclust(:))) % 53 entries
 Aclust
 
 %% A4.01 Adj (8x8) estm. enter by hand, rounded
 Adj1 = [0.0,   8.54,  9.1,    0.39,  0.84, 4.56, 0.045, 0.23; ...
         0.42,  0.0,   1.46,   0.03,  0.44, 0.063, 0.033 0.13, ; ...
         0.23,  0.596,  0.0,   0.325, 0.09,  0.003, 0.0, 0.004; ...
         0.53,  1.14,  4.88,   0.0,   2.2,   0.23,  0.0, 2.01; ....
         0.76,  1.28,  0.68,   1.19,  0.0,  1.64, 0.029, 1.81; ....                  
         0.68,  0.17,  0.005,  0.0,   0.114, 0.0, 0.13, 0.145; ...
         0.114, 0.51,  0.001,  0.006, 1.04,  2.73, 0.0  1.426; ...
         0.05,  0.54,  0.01,   1.12,  3.38,  0.74, 0.10, 0.0]; % checks ok
     
length(find(Adj1(:))) % 29 entries
 Aclust = Anew(clusterlist, clusterlist)
 
clear Anew Aclust

%% Appx 4.1 Link distances of 8 nodes
 %  d(i,j)
 addpath(genpath('/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain')); % Data files;  include sub-directories 
fprintf('  use marmoset d(i,j) from Atlas \n')
clusterlist=[2, 26, 25, 48, 32, 3, 45, 47]; % as plotted in x-y plane, viewed from Ant.
Adist= csvread('MarmosetDistPairs.csv');  % based on Atlas; saved file 
% csvwrite( 'MarmosetDistPairs.csv',Adist); % 
%%
%figure; hist(Adist(:), 100); hold on  % ~ skewedGaussian, peaks at 7, 11, 13 mm
Adist8 =Adist(clusterlist,clusterlist);
figure; hist(Adist8(:), 50); title('Marmoset: 8 NM cluster d(i,j)');  xlabel('d(i,j)  (mm) ')
clear Adist 
%% Appx 4.1.2  set d in bin of 0.1 mm (approx)
Ad6A = [...
     0     3.4  2.9  2.4   2.4  2.5, 4.3, 5.0; ...
     3.4   0    1.7  3.4   4.1  2.7, 0.0, 4.7; ... % 0 if no link present
     2.9   1.7  0    1.9   2.8  2.6, 0.0, 3.2; ...
     2.4   3.4  1.9  0     1.8  3.4, 3.2, 2.8; ...
     2.4   4.1  0    1.8   0    2.8, 2.1, 3.2; ...
     2.5   2.7  2.6  3.4   2.8  0.0  3.6, 4.8; ...
     4.3,  4.9, 3.6, 3.2,  2.1, 3.6, 0.0, 2.4; ...
     5.0,  4.7, 3.2, 2.8,  3.2, 4.8, 2.4, 0.0];   % checks ok.
 
%%  Appx 5.  cluster CM & boundary; 3D plot 
% read coords in MarmosetReadData.m code :  NodeCoord(1:3, :) : n x 3 array
NodeCoord=csvread('CoordsMarmoset.csv'); % nb. some absent: NaN 
clusterlist=[2, 26, 25, 48, 32, 3]; % as plotted in x-y plane, viewed from Ant.
XYZclust = NodeCoord(clusterlist, :) 
CMcluster  = mean(XYZclust(:, 1:3))

figure; hold on
set(gcf,'Renderer','OpenGL');
set(gca, 'Color', [0.7 0.7 0.7]) % background gray
 %axis vis3d  % freeze aspect ratio (during rot)
axis equal  % get aspect ratio correct
 %set(gca, 'Projection', 'perspective');  % orthographic vs perspective view 
theta=230; grid on;  % Start view, for p-FC (RH, front  : choose these as needed
view(theta, 20); % to begin rotation scan
%axis([-11 11 -10 20 0 20]) % full
axis([-11 11 0 20 0 20]) % front half only
camlight(-120, 30); camlight(105,50); camlight(-40, 20);  % light from above & side
camlight(45, 60) % & above, rear 
floorgridOnly % % add the grids: floor; draw 5 mm grid & marker on floor
transp=0.5;

for i=1:6
 xi=XYZclust(i,1); yi=XYZclust(i,2); zi=XYZclust(i,3); % col-1 is ID# 
 ii=clusterlist(i); % get correct ID#
   plot3(xi, yi, zi, 's', 'MarkerFaceColor', nodecolors(i,:),'MarkerEdgeColor', nodecolors(ii,:),'MarkerSize', 10) % square, Schematic coord
   text( (xi-0.1), (yi+0.15), (zi+0.1), Acrn{ii},'FontWeight','Bold','FontSize',11); % Label the pt.
      % ,'Fontname','Times New Roman'
end

plot3(CMcluster(1), CMcluster(2), CMcluster(3), 'h', 'MarkerFaceColor','g','MarkerEdgeColor', 'g','MarkerSize', 10) % square, Schematic coord
text( (CMcluster(1)-0.1), (CMcluster(2)+0.15), (CMcluster(3)+0.1), 'CM','FontWeight','Bold','FontSize',11); % Label the pt.


%% Appx 6. Cort'l Col'm dimensions - from sixClusterWavesV4b.m 
 % from  code sixClusterWavesV4b.m; cf w/b 6, p93 & Marm voxel data, bv vol.
fprintf('\n Geom.  A-8b')
 %vol = 23.92; % A-10 (mm^3)  %vol = 12.0; % A-11 (mm^3)  %vol = 11.4; % A-32
 %vol = 2.68; % A-32V (mm^3)  %vol = 3.29; % A-46D  %vol = 7.74; % A-9  
 %vol = 13.0; % A-8aD [NM #7
vol = 17.83; % A-8b   [NM #8
rsphr = (3*vol/(4*pi))^0.333  % (mm)

 %thick= 2.0; % A-10. cort thickness (mm) [ref Atapor Cereb Ctx'`19, Fig 10]
 %thick= 1.4; % A-11  %thick= 1.6; % A-32  %thick= 1.0; % A-32V  %thick= 1.5; % A-46D
 thick= 1.8; % A-9, A8b;  %thick= 1.5; % A-8aD
rcolm = sqrt(vol/(pi*thick))
Acol = pi*rcolm^2
PerimCol=2*pi*rcolm

clear vol rsphr thick rcolm Acol PerimCol

%% & for other Source Hubs
 fprintf('\n Geom.  TPO') % #103 TPO
  vol = 16.71; % TPO
rsphr = (3*vol/(4*pi))^0.333  % (mm)
 thick= 1.8; % est from Attapour, Fig 10
rcolm = sqrt(vol/(pi*thick))
Acol = pi*rcolm^2
PerimCol=2*pi*rcolm 


 