% B {Pailtjorpe, U Sydney. sixNMclusterVx.m, (28/10/21 - 24) 
% Jansen-Rit model, 6 node cluster at A-10, Runga-Kutta integration
  % 6 node cluster at A10 (pFC), JR NM model;
   %  use RK4 integrator fn JR12(.., NNinputs), with stim. forcing
% Dependencies:  JRv5d.m for DE integrator; WaveTrains.m for stimuli set up 
 % Data file: AdjMarmoset_Rescale2b.csv  % 116x116 link wts (based on LNe) 
    % V1. dev.(23/5/21)  %  based on sixNMloopRK4_JR12V2b.m
    % Updade, V2: use sigmoid filter on coupling  (12/6/21)
    % V2b. modulate stim Ip+ran to selected nodes only (7/5/21)
    % V3 use JR12a(t, y, Ip, p, I12, a, b, A, B, C) : pass params to tune freq band;
    % V3a. pass Stim, pran via Sigmoid: pulse->V transfer (4/7/21)
    % V4. pass r to fn J12a.m; vary at ea NM
    % V4.1 vary Gain[s] ~ syn strength, per NM; pass Cfac (C1, C3) to fn J12a.m (24/7)
    % V4a. add signal delays, via d(i,j)
    % V5. Sigmoid S[v], erfc based on Distribution(wt/link) {in RK4 subr (29/9/21) % cf. workbook#4, p59
    % V5d [] - erfc((1/V)) form for N(wt) distrn; uses "if", JRv5d.m 
       % also used for perturbed NM params & link weights (#1.0.2) (./9/22)
       % nb. signal delay assumes v = 1 ms/ (mm/ms) @ #4.0
       % add restart from last sim: write/ read y(50 steps)
       
% clear all; close all

% 1.0 Set up
addpath('/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain/Models'); % for other model codes
fprintf('\n 6 Neural Masses, JR/RK4, star cluster; v5d; wt Adj; S(v * wt) form \n')
fprintf('\n  alpha band; stim: 0 or 100pc e or i; S[erfc], r & wbar +  stim \n')
 % Adj, for coupling nodes
nn=6  %nn=6 % number of neural masses (nodes)
 % examine star cluster about Out-hub A-10  (marmoset logbook, p58, 66 ("v1")
  % n-node, Marmoset star clusters: at 7ms: 6 nodes; @7.5 & 8.1 ms: 8 nodes

  %%  1.0.1 Approx  Adj2, .cf Appx. below: this is similar to Adj2, with rounding
   fprintf('\n   use approx Adj (rounded)  \n')
   Adj= [0.0,  9.1, 0.4, 0.8, 4.6, 8.5; ...  % #1 is out hub;
        0.25, 0.0, 0.4, 0.1, 0.003, 0.6; ...  % cf Workbook 3, p18a for list of Links
        0.5, 4.9, 0.0, 2.2, 0.23, 1.1; ...  % scale wt wrt to 1k (LNe)
        0.75, 0.7, 1.2, 0.0, 1.6, 1.3; ...
        0.70, 0.005, 0.0, 0.1, 0.0, 0.2; ...  %  
        0.4, 1.5, 0.03, 0.44, 0.1, 0.0];     %  
  % cf. Appx 2, below for check of Adj against orig Anew(LNe)
  % A(4,2), A(5,6) was missing; A(6,3) was too big; others are close
 
   %%  1.0.2 Approx  Adj2, & adjustments
   fprintf('\n   use perturbed approx Adj (rounded)  \n')
   Adj=[0.0,  9.1, 0.0, 0.8, 4.6, 8.5; ...  % #1 is out hub;
        0.0, 0.0, 0.0, 0.0, 0.00, 0.6; ...  % cf Workbook 3, p18a for list of Links
        0.5, 4.9, 0.0, 2.2, 0.23, 1.1; ...  % scale wt wrt to 1k (LNe)
        0.75, 0.7, 1.2, 0.0, 1.6, 1.3; ...
        0.70, 0.00, 0.0, 0.0, 0.0, 0.0; ...  %  
        0.0, 1.5, 0.0, 0.0, 0.0, 0.0];     %  2-5 & 6-3 to 0; then 2-4, 5-2, 6-5 to 0
   Adj                                    % then 1-3, 2-3, 5-4, 6-1, 6-4 to 0
   length(find(Adj(:))) % 27 entries
%% 1.0.2 unit Adj & adjustments
fprintf('  test Adj=1:  \n')
Adj=ones(nn);
 % elim diagonals
 for i=1:nn
     Adj(i,i)=0;
 end
  % Adj
 %length(find(Adj(:))) % 30 links
 % 1.02 zero Adj
  %fprintf('  test Adj=1:  \n')
  %Adj=zeros(nn);

  fprintf(' + strengthen some links; elim. & weaken others:  \n')
  Adj(1,2)= 9; Adj(1,6)= 8;  % selectively strenthen into 2,6
  Adj(1,5)= 4; Adj(3,2)= 4; Adj(3,4)= 2; % & into 5, 
  Adj(6,3)=0; Adj(5,3)=0;  
  Adj(2,4)=0.1; Adj(2,5)=0.1; Adj(3,5)=0.1; Adj(5,2)=0.1; Adj(5,4)=0.1; Adj(5,6)=0.1; % NM#5 has weak in- & out-links 
  Adj(6,5)=0.1;  % for closer matching to Marm LNe2 data
  Adj(1,3)= 0.5; Adj(2,1)= 0.5; Adj(2,3)= 0.5; Adj(2,6)= 0.5; Adj(3,1)= 0.5; Adj(5,1)= 0.5; 
  Adj(6,1)= 0.5; Adj(6,4)= 0.5; Adj(4,2)= 0.5;
  Adj(4,5)= 1.5; Adj(6,2)= 1.5; % selectively strenthen into 2,5,6 
  Adj
  nLinks =length(find(Adj(:)))
  TotLinkWts =sum(Adj(:))
  clear nLinks TotLinkWts 
  
%% 1.0.3 Zero Adj: no coupling
  fprintf('  test Adj=0: no coupling  \n')
  Adj=zeros(nn);
  
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
Ad6A= Ad6A(1:6, 1:6)

 %% 1.0.5  zero dij
 fprintf('  d(ij) =0: no signal delays  \n')
 Ad6A = zeros(nn);

%% 1.0.9  random Adj (U ran < 1)
fprintf('  U. random Adj; wt <1  \n')
rng('shuffle') % seed rand, to vary order...
 Adj=rand(nn);
 % elim diagonals
 for i=1:nn
     Adj(i,i)=0;
 end
 Adj
 mean(Adj(:))
 AdjRan=Adj; % save, for below
 
%% 1.1 Adj from Marmoset (LNe) data
 % cf.Appx 2. in orig order: cf. logbook, p58
  fprintf('  use marmoset LNe Adj(2):  \n')
clusterlist=[2, 26, 25, 48, 32, 3]; % as plotted in x-y plane, viewed from Ant.
Anew= csvread('AdjMarmoset_Rescale2b.csv');  % based on LNe 
 Aclust = Anew(clusterlist, clusterlist);
Adj2=Aclust/1000;  % scale used for 6 NM model
Adj =Adj2;   % for main code
 %fprintf(' + selectively tweak some links:  \n')
  %Adj(6,3) =0; Adj(5,2) = 0.0;% were v weak links originally
     % Adj(5,1) = 0.5; % ??
     %Adj(1,2) = 7; %Adj(1,6) = 6; % weaken strongest links %Adj(1,5) = 3;
     %Adj(5,1) = 0.5; % tweak "minor" links %Adj(2,6) = 0.4; %Adj(4,5) = 0.5;
       %Adj(2,1) = 0.5; Adj(3,2) = 4; Adj(4,5) = 1; % 1)tweak group of "minor, peripheral" links
         % Adj(5,6) = 0.13; % 2) alt
       %Adj(5,6) = 0.26; Adj(3,4) = 1.65; Adj(2,6) = 0.4; % 2a)tweak group of links
       %Adj(2,6) = 0.6; Adj(5,1) = 0.5; Adj(4,5) = 1; %Adj(4,3) = 0.1; % 3, 3a, 3b)tweak group of links
Adj
nLinks= length(find(Adj(:))) % 29 entries [originally]
%TotLinkWts =sum(Adj(:))

clear nLinks Anew Aclust  clusterlist Adj2 Tot*

% nb. set approx dij in Ad6A, above & at Appx 3.2.1

%% 2.0 Parameters: tune freq. band
 fprintf('  : cf Marmoset out hub, A10 (3-D cluster, 6); beta band; wt x \n')
 %nLinks = length(find(Adj(:))) % check # Links: was 22 [nb. some v small & omitted], now 29
 %clear nLinks
% time axis & IC (0 mV)
x0=0; xf= 20.0; % t range (sec)
h = 0.0001;   % t step (1 ms; [& 0.5, 0.25 ms, for debug])
x = [x0:h:xf];  % t domain array
 %time =x/(1000);  %time =x/(1000*h); % time (sec) as row vec
tstart=4001; %tstart=2000; %1000; % for plots (@ 1, 0.5, 0.1 ms time steps)
%y = zeros(6,length(x)); % row vec for each variable
%yic1=[0.0; 0.0; 0.0; 0.0; 0.0; 0.0]; % IC at t=0 for 1 node (6 DEs),  col vec, 6 rows
%yic1=[1.0; 2.0; 1.0; 0; 0; 0]; % need non-zero IC
yic1=[0.5; 1.0; 0.5; 0; 0; 0]; % smaller, non-zero IC
yic=yic1; % add additional nodes
for in=1:nn-1
 yic=[yic; yic1]; % add 6 rows per extra node; 12, 18 DEs for 2, 3 ... nodes
end
 clear yic1
%
% Model parameters (orig JR) : needs to be in the fn subroutine!
v0=6;    % (mV) midpoint (switching trheshold V) of sigmoid function, for V -> pulse conversion 
vm=5;    % (s^-1) max firing rate for sigmoid fn

 % activation rate;  Steepness of sigmoid function 
 %rvec= [0.3 0.3 0.3 0.3 0.3 0.3]  % pass to indiv nodes, via fn J12a.m
%rvec= [0.4 0.4 0.4 0.4 0.4 0.4] % test all alpha, or beta (low)
 %rvec= [0.45 0.45 0.45 0.45 0.45 0.45] % test - in beta
rvec= [0.5 0.5 0.5 0.5 0.5 0.5] % test: all beta: most NM dont osc!
 %rvec= [0.55 0.55 0.55 0.55 0.55 0.55] % test:
  %rvec= [0.6 0.6 0.6 0.6 0.6 0.6] % test - for wt form
  %rvec= [0.7 0.7 0.7 0.7 0.7 0.7] % test - for wt form
  %rvec= [0.56 0.40 0.49 0.46 0.40 0.49] % set r ~ sqrt(vol); cf w/b#3, p80: tune in beta: ok
 %rvec= [0.56 0.40 0.47 0.43 0.38 0.47] % revised r ~ sqrt(N); cf w/b#4,p42 & 5,p72: tune in BETA band
   %rvec= [0.42 0.30 0.37 0.34 0.30 0.37] % ditto: tune in alpha band: corrected (23/7/21)
%rvec= [0.45 0.32 0.38 0.35 0.30 0.38] % tune in ALPHA (or theta) band: wb 5, p70 (18/9/21)
    %rvec= ones(1, nn)*mean(rvec) %-0.05  % test: all = mean [& +/- sigma]

% weight/synp link, for Sigmoid response fn:
  %wvec=[0.52 3.23 0.48 0.74 1.30 2.35] % raw wt-In/#link-In: av in each NM: wb4,p59a
  %wvec=[1.08 6.73 1.00 1.54 2.71 4.90] % ditto: rescale so >=1
  %wvec=[2.6 1.74 1.14 0.92 3.58 1.41] % based on Av[wt-In/#link-In & wt-Out/#link-Out] [wb5],p102a
        % nb ^ w(5) in error, should be 0.77  ??
%wvec=[3.4 2.3 1.5 1.2 1.0 1.8] % rescaled Av[w-in/k & w-out/k] wb5,p102a; [8], p64
%wvec=[4.7 2.3 1.8 1.1 1.3 2.4] % trial 8.3) use max{in, out, av} -wb[8], p65
wvec=[2.6 3.2 2 1.5 2 3] % trial 8.3) mixed -wb[8], p66 & adjust

%wvec= ones(1, nn)*mean(wvec) -2.35  % test: all = mean [& +/- sigma]
  %wvec(2)=5.0; wvec(6)=4.0; wvec(5)= 1.0; wvec  % perturb largest one-by-one
    %wvec=[0.40 2.50 0.37 0.57 1.01 1.82] % 1st revise, to reset #2
    %wvec=[0.32 2.00 0.30 0.46 0.81 1.46] % 2nd revise
    %wvec=[0.16 1.00 0.15 0.23 0.4 0.73] % alt form ?

% tuned to alpha/ beta band
%C=350  % gamma;  C= 135  % 
 C= 250   % basic (scale) connection strength, within NM % defauly J-R value         
%C= 300 % beta band [5], p72
 %C= 135 % alpha band -  cf tunimg 3NM,[wb 7], p15 / 6NM p57
R12=200;   % Adj/5;  basic NM - NM coupling coeff (ie. link wt) (cf. Goodfellow)
C1=C;     % Pyramidal to excitatory interneuron connection 
C2=0.8*C; % (108)   Excit interneuron to pyramidal connection [feedback]
C3=0.25*C; % (33.75)  Pyramidal to inhib interneuron connection [weaker feedback]
C4=0.25*C;
% tune feedback gain (baed on WtIn, kIn): same as wvec above!
CfacVec= [1.0; 1.0; 1.0; 1.0; 1.0; 1.0]; % pass to fn J12a.m
 %CfacVec= [6.33; 1.0; 7.0; 4.67; 2.67; 1.33]; % follows wtIn/kIn for 6x6(wb3, p88, 90): Poor
 %CfacVec= [4.6; 3.07; 1.03; 1.52; 1.0; 1.92];  % ?recalc for 8x8, wb#3, p101 (4/8/21)
  %CfacVec= [6.39; 4.26; 1.0; 2.11; 1.39; 2.66; 1.76; 2.76];  % ?recalc again (wb#3, p104 & XL)
  %CfacVec= [4.67; 3.1; 2.6; 2.19; 1.81; 2.44; 1.0; 1.65];  % ?recalc for 8x8, wb#3, p104 (6/8/21)
%CfacVec=wvec;
%CfacVec' % Assign in J12x.m:C1=C1*Cfac;  C3=C3*Cfac; 
 % signal delay: set by array Ad6a, set above, based on d(i,j) % d12~1.5-4;  % eg. mm / ms
% 
A=3.25;  B=22;  % Max PSP amplitude for pyramidal (A), inhib (B) & excit (kA * A) interneurons
  %a=100; b=50;  % (s^-1) / default JR; set to osc at 8Hz (tau's: 20, 20 ms)
% Time Const. - equal for each node
% taue= 10; taui= 10; % set time const (excit, inhib) (ms) :  BETA band, wb5, p72
%taue= 25; taui= 25; %  (ms)  :  ALPHA band : results in theta!
%taue= 16; taui= 17; % also beta, ok.
   %taue= 20; taui= 20; %  (ms)  :  ALPHA - another test
taue= 10; taui= 10; %  (ms)  :  BETA band 
% taue= 5; taui= 5; % gamma band
    % move freq band; keep a*A & b*B const
a=1000/taue; b=1000/taui; % rate const (s^-1)
 %A=3.25*100/a; % rescale, to keep a*A const. dont tune B (lower C)?
 %B=22.0*50/b;  % ?dont tune B (more inhib): need lower C then.
[taue taui] % debug / use default A,B
[a b]
[A B C]
[round(a*A) round(b*B)]  %[a*A b*B] % "damping" coeff
 round(b*B)/round(a*A)
 
tstart = 40001;  % stim at 4.0 s (0.1 ms steps0

% >>>>>>>>>>>>>>
%% 3.0 Zero stim + ran noise :: Base case
 fprintf('\n >> Zero inputs:  1 x G ran (pos, neg) noise \n')
% fprintf('\n >> Zero inputs: 0.1mV Gaussian ran noise \n')
Vyt =zeros(1,length(x));
pran=zeros(1,length(x)); 
Ip=zeros(1,length(x));  
IpOn=zeros(nn,1); IpOnB=zeros(nn,1); % need IpOn also
Pd=zeros(1,length(x)); Pd2=zeros(1,length(x)); %  pulse density fn - cf #3.6
 % Pd=zeros(6,length(x)); % need 6 x stim later
 %pran = 0.1*rand(1,length(x)); % low noise;  Uniform, pos ran noise:
 pran = 1.0*randn(1,length(x)); % Gaussian noise, pos/neg

%% 3.1 set up Stimulus pulse (spike) train
  %fprintf('  * stimulate node 2&6 (A32V, A11 ) both, equally:  \n')
 %fprintf('  * stimulate node 6 (A11 ) only:  \n')
 %fprintf('  * stimulate node 2 (A32V ) only:  \n')
%fprintf('  Pulse train stimulus - 10 mV; to NM #2 amd/or #6  \n')
stim_freq=16  % nb. also need to change npulse, to get 1 sec long stimulus
Ip=zeros(1,length(x)); % 1st input (NM#2):excit curents impulse train: input to pyramidal population
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
   %Ip(tstart+ i*tdelay +dphi +1: tstart+ i*tdelay +dphi +1+twidth-1)=1.0; % wider pulse

end
% other inputs; set scale
 %Ip=zeros(1,length(x)); % debug: test zero excit'n. : just follows the noise!
amplP=10.0  %mV % nb.max rate of 5Hz at ~25mV
 Ip=amplP*Ip;  %Ii=60*r*I; % excit & inhib [set Ii in ext fn]
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
clear P1 P2 Pfilter pran100 Y f

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
%tstart=40001; % for plots
StimLevel= 100
Pd=zeros(1,length(x));  tstart0= 2001;  % of this stim; avoid transients
Pd(tstart0:end)=StimLevel; % && turn on the stim!!

% AND - set Input Channels: at #3.7:

clear StimLevel PulseHeight tstart0

 %% 3.5.1 Modulated stim rate - use Pd x |Vyt|, full rectification; not filtered by S[v]
% modulate by Vyt fn - set up in WaveTrains.m
%fprintf('\n ++ modulated stimulus - via  Vyt(f/2): at 4.0 sec ')
fprintf('\n ++ modulated stimulus - const +  Vyt(f): at 4.0 sec ')
tstart=40001; % start stim & for plots
StimLevel= 300.0
Pd=zeros(1,length(x)); %tstart0= tstart; %tstart0= 1001; %tstart0= 20001; %  separate for Pd
 %tend=50001; % for short stim [4:5 sec] etc
 %Pd(tstart:tend)=StimLevel*abs(Vyt(tstart:tend)); % short, const stim
Pd(tstart:end)=StimLevel*abs(Vyt(tstart:end)); % modulated stim! pos. via abs(), so 2*f
%Pd(tstart:end)= StimLevel + 0.1*StimLevel*Vyt(tstart:end); % modulated stim! pos. via abs(), so 2*f

StimStats= [mean(Pd) std(Pd) max(Pd)]
PulseIntegrl = sum(Pd)*h % sum across rows % trapezoidal rule
 test=Pd(tstart:end); %test=Pd(tstart:tend);
 StimStats1sec= [mean(test) std(test) max(test)]

 Vyt =zeros(1,length(x)); % reset Vyt: no V wave
 
 % PdSave=Pd; % for reuse
 % Pd=zeros(1,length(x));  tstart=40001; tend=47001; % reset
 %  Pd(tstart:tend)=PdSave(tstart:tend); % short, modulated stim
 
 InI =0  % controls stim to Excit/Inhib popn:  apply stim to Excit popn (default)
 % InI =1 % apply stim to Inhib 
 
clear StimLevel PulseHeight tstart0 tend PulseIntegrl StimStats test


%% 3.6 Load 1st Pulse density function[s]
 % a la JR('93 & '95); set up in code WaveTrain.m @ #1.4
 fprintf(' Pulse shape stim +0 (@4.0 s); & 1 mV G ran noise: \n')
 Pd= 2.0*Pulse1; % 1 or 2 pulses;  vary amplitude (ie. pulse density fn.
 tstart=40001;  % cft wave... code, for pulse setup
  % Pd= 1.0*PulseN; % N pulses, in train
 %Pd = Pd+100.0;  % add const stim (all t)
 Pd((1001):end) = Pd((1001):end)+ 0.0;  % add const stim , with pulse (not at t=0)
  PulseHeight= max(Pd)
  PulseIntegrl = sum(Pulse)*h; % sum across rows % trapezoidal rule 
  PulseIntegrl=PulseIntegrl'
  %IpOn=zeros(nn,1); % switch Pulse stim on/off to ea node
  %IpOn(2)=0; % 1st input: stimulate Node 2 only (A32V: strongest input links)
  %IpOn(6)=1; %0.33; % 2nd input, to node 6 (A11), lower wts, 
 Pd2=zeros(1, length(x)); % unless set & used, below
  % pran = 0.1*rand(1,length(x)); % low noise;  Uniform ran noise:
 InI =0  % controls stim to Excit/Inhib popn:  apply stim to Excit popn
  %InI =1 % apply stim to Inhib 
  clear PulseHeight PulseIntegrl 

 %% 3.6.1 Load 2nd Pulse density function[s]
 fprintf(' 2nd Pulse shape stim: \n')
 Pd2= 1.0*Pulse2; % 1 or 2 pulses;  vary amplitude, width (ie. pulse density fn)
 Pd2 = Pd2+ 0.0;   % add extra const stim
 
%% 3.7 Input channels
 % selected nodes to get the stimulus inputs: Ip, Pd 
  IpOn=zeros(nn,1); % switch 1st stim on/off to ea node
 % IpOn=ones(nn,1); % 1st stimulus On for all nodes
 % stimulate Node 2 (A-32V): strongest, more wt, fewer links)
  IpOn(2)=1; % 1st input: more wt; fewer #Links*0.61; % 1st or 2nd input?
  IpOn(6)=1;  %0.64; % 2nd input?, wt*0.64, to node 6 (A-11), lower wts, more Links  
 
 % stimulate Out Hub [1]A-10
 % IpOn(1)=1; % 0.4;  % A-10 <- TPO, TE3,2,11, Aud (weak), Opt, TEO, V4

 % stimulate other nodes:
 % IpOn(5)=1;  %0.8 % [5] A46D <- MO-pre
 %IpOn(4)=1;  %0.8 % [4] A9 <- Cing/RSP  
 IpOn' % check
 
 % direct 2nd pulse, if present?
 IpOnB=zeros(nn,1); % switch 1st stim on/off to ea node
 %IpOnB(2)=1; % 2nd Pulse to 1st input:
  %IpOnB(6)=1;  % 2nd Pulse to 2nd input (A11)
 % IpOnB' % check
 
 %% 3.8 Gather stimuli & PLOT, to check stim inputs
 %pran= [pran100, 0] % need N=1 elements (5001)
 %Vyt=Vytburst;  % to pass to DE solver
 
 figure;  plot(pran, 'Color', [0.9 0.5 0.5]);  hold on
 plot(Pd+Pd2 , 'b'); title('pulse + Vstimulus + waves + noise');
  plot(Vyt, 'g--'); xlabel('t steps (0.1 msec) ')
  ylim([-10 max(Pd)+20])
text(1.5e5, (max(Pd)+7), strcat('Stim To  :  ', num2str(IpOn')) )
  
% >>>>>>>>>>>>>> 

%% 3.9 Optional:: read y(last 50 steps) for restart
 % cf. #5.0.3 for save file
  fprintf(' \n > read saved 50 t-steps for sim restart \n')
y=zeros(nn*6,length(x)); % reset array
y_end= csvread('y_6NMerfc_a.csv'); % was saved for restart {case of NM#3 does osc}
y(1:36,1:50) = y_end;
 size(y_end) % debug
 % figure; plot(y(2, 1:50)) % debug
  % y_6NMerfc_a = y_end; % Here NM#3 does osc. : full y is too big
clear y_end

%% 4.0  RK4 DE solver, 1st order, for n variables
 % nb. matlab has default double precision: needed here
   % check timulus
   %figure;  plot(Ip+pran); title('stimulus + waves + noise'); hold on; plot(Vyt, 'r-'); xlabel('t (msec) ')
% reset arrays
fprintf(' \n > start sim from y-IC ...  ')
y=zeros(nn*6,length(x)); % set up array for solutions; 2,3, etc nodes (6 rows each)
 %y=double(y); % enforce double precision (64 bit) : should be default
deltay=zeros(nn,length(x)); % set up array for y3-y2; 6 nodes (1 rows ea)
y(:,1)= yic ; % y's are row vec, one for ea DE; IC is a col vec, for i=1 or t=0. 
  for i=1:50  % need (IC > 0) to escape "trap" at zero
   y(:,i)= yic ; % extend IC - for delay DE "late start" - eg to 5 ms (0.005 s)
  end
  
  % 4.0.1 optional restart :: edit out
  %if restart from last sim, enter here: with y(50 steps) read in at #3.9
   % fprintf(' \n > restart last sim \n')
   % deltay=zeros(nn,length(x)); % set up array for y3-y2; 6 nodes (1 rows ea)
    % use initial y(i, t) from # 3.9
    
In12=0;  % i-j interaction; 
Istim =0; % default: cf InI - stim to e only; vary below, as specified
   % signal delay assumes v = 1 m/s (mm/ms)
   vlocal=1.0; % unmyeln axon vel ~ signal delay ~ 0.1-1.0 m/s 
  
% RK4 DE solution step:  iterate forward in time (here x); scan nodes(in) & linked NN (j)
  % fn JR1s incl forcing fn. as arguments, at this t-step; 
  % inputs: pran + p(i,j)+ [NM-NM couplings] {= sum(y1-y2) ; ext stim Pd etc
 %fprintf('\n   i  j wt(ij) dr(ij) index r Cfac ') % debug- header
tic
        % start "late" to allow for delay back ref, to L [allow 0.5/5 ms w 0.1 ms steps]
for i = 50:(length(x)-1) % scan time steps & start at 1, 10, 100 [for delay DE]
    for in =1:nn % scan nodes in cluster (NM)
         j11=(in-1)*6+1; j16=(in-1)*6+6;  % working on this node #in
         r=rvec(in); Cfac=CfacVec(in); wbar=wvec(in); % for indiv NM: slope of S[v]; Gain C1,3; wt/link
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
           %deltai=int16(Ad6A(jj,in)/(1000*h*vlocal)); % adjust index, for signal delay time from nn, jj->in / as integer
          deltai=int32(Ad6A(jj,in)/(1000*h*vlocal)); % check discontinuities
          dyi= y(j11+1,i-deltai)-y(j11+2,i-deltai); % delta-y of nn.
          argm=r*(wbar -v0/dyi);  In12=0; % need to exclude neg argm & V~0 :: OK
          if argm >0 & dyi>0.01 
           In12= R12*Adj(jj,in)*vm*(1.0-erfc(argm)/2); % erfc form for wt; avoid v=0 : modulated by step fn
          end
          In12=In12+Vyt(i); % add travelling/standing wave here, apply direct to soma
    %fprintf('\n  %2.0f %2.0f %5.3f %4.1f %4.0f %4.2f %4.2f ', in, jj, Adj(jj,in), Ad6A(jj,in), (i-deltai), r, Cfac )
           end % check for NN  
       end % scan linked NNs; gather inputs (in <- jj) from linked NNs
           %fprintf('\n  ') % for debug
           %[jj j11 j16]  % debug:
        j11=(in-1)*6+1; j16=(in-1)*6+6;  % update node #in
          % nb. y7:12 in Node-2; y13:18 is Node-3, etc.
        % Stimulus pulses to selected node (eg 2, 6), via IpOn switch
        % add travelling wave Voltage, at y=16mm (Ant), at this time(i), to _all_ NM
        %Istim=Vyt(i)+IpOn(in)*Ip(i); %pstim=IpOn(in)*pran(i); % turn stim on only for selected nodes - set above, at IC.
           % ?? tbd:  Istim=vm/(1+exp(r*(v0-wbar*Istim)));  % filter Istim pulses
        % gather ran, pulse[s] & wave train together: to all NM
        pstim=pran(i) +Vyt(i); %  Vyt & ran noise to ALL nodes
         Pdstim=IpOn(in)*Pd(i); % switch 1st stim pulse on, to selected NM;
         Pdstim=Pdstim +IpOnB(in)*Pd2(i); % also switch 2nd stim pulse on, if present?;         
        k_1 = JR12v5d( x(i), y(j11:j16,i), Pdstim, Istim, pstim, In12, a, b, A, B, C, Cfac, r, wbar);  % these k's should be col vec, one for ea DE
        k_2 = JR12v5d( x(i)+0.5*h, y(j11:j16,i)+0.5*h*k_1, Pdstim, Istim, pstim, In12, a, b, A, B, C, Cfac, r, wbar); % seem to need to force col vec ?
        k_3 = JR12v5d( (x(i) +0.5*h), (y(j11:j16, i) +0.5*h*k_2), Pdstim, Istim, pstim, In12, a, b, A, B, C, Cfac, r, wbar);
        k_4 = JR12v5d( (x(i)+h), (y(j11:j16,i) +k_3*h), Pdstim, Istim, pstim, In12, a, b, A, B, C, Cfac, r, wbar);
        y(j11:j16,i+1) = y(j11:j16,i) + (k_1 +2*k_2 +2*k_3 +k_4)*(h/6); % load y1:y6 at this t-step      
    end %  scan nodes
end  % scan t-steps
fprintf('\n  done! ')
toc %timing
 clear i in deltai In12 j*  Istim pstim Pdstim stepfn dyi 
 
%% 5.0 PLOT OUTPUT 
tstart=40001; % if needed? : start of stim.
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
 axis([1.5 2.5 10 100]);  grid on; xlabel('t (s) ')
  % figure; plot(x, (y(2,:)-y(3,:)) ); title('6-NM, node-1, delta-y '); xlabel('t (sec) ')
  
% Node 2 outputs [this has the primary inputs]:
 %figure; plot(y(8,(tstart-100):end) - y(9,(tstart-100):end));  title('JR node-2: output, y2(2) - y1(2)')
  %figure; plot(y(7:12,:)'); title('JR/ node-2: output, y1, y2, .. y6')  % nb. need cols
  %legend('y0', 'y1', 'y2', 'y3', 'y4', 'y5'); % 
 figure; subplot(2,1,1); plot(x,y(8,:)'); hold on; plot(x,y(9,:)'); 
 title('JR/ node-2: output: y1, y2 '); legend('y1(2)', 'y2(2)'); hold off
 subplot(2,1,2); plot(x,y(8,:)'); hold on; plot(x,y(9,:)'); 
 axis([1.5 2.5 0 100]);  grid on; xlabel('t (s) ')
     % other outputs
     %figure; plot(y(6,(tstart-100):end));  title('JR model, Node-1 output, y5  ')
     %figure; plot(y(12,(tstart-100):end));  title('JR model, Node-2 output, y5  ')
     % axis([0 5 -500 0]);  axis([0 20  -10000 10000]);
     figure;  plot(y(8,(tstart-800):end) - y(9,(tstart-800):end)); title('Node 2, delta-y ')
     xlabel('t (x0.1 ms) ')
 
 % Node 3 outputs:
 figure; plot( y(14,1:end) - y(15,1:end) ) ;
 title('JR model: node-3 output, y2(3) - y1(3)')   
% Node 4 outputs:
 %figure; plot(y(20,(tstart-800):end) - y(21,(tstart-800):end));
  %title('JR model: node-4 output, y2(4) - y1(4)') 
% Node 5 outputs:
 figure; plot(y(26,(tstart-800):end) - y(27,(tstart-800):end));
 title('JR model: node-5 output, y2(5) - y1(5)')
 figure; subplot(2,1,1); plot(x,y(26,:)'); hold on; plot(x,y(27,:)'); 
 title('JR/ node-5: output: y1, y2 '); legend('y1(5)', 'y2(5)'); hold off
 subplot(2,1,2); plot(x,y(26,:)'); hold on; plot(x,y(27,:)'); 
  grid on; xlabel('t (s) '); %axis([1.5 2.5 0 50]); 

% Node 6 outputs:
 tstart= 40001; % debug
 figure; plot(y(32,(tstart-100):end) - y(33,(tstart-100):end));
 title('JR model: node-6 output, y1(6) - y2(6)')
 figure; plot(y(32,:)'); hold on; plot(y(33,:)'); 
 title('JR/ node-6: output: y1, y2 '); legend('y1(6)', 'y2(6)');

% all 6 nodes together:
tstart= 40001; % debug
 figure; plot(x((tstart-800):end), y(2,(tstart-800):end) - y(3,(tstart-800):end), 'k' );  hold on % node #1, A-10 
 xlim([4 4.2]); xlabel('t (s)'); xlabel('t (ms)');  %pause
 plot(x((tstart-800):end) ,(y(8,(tstart-800):end) - y(9,(tstart-800):end) ));   % #2, A-32V
 %pause  % - to closely examine waveforms
 plot(x((tstart-800):end) ,y(14,(tstart-800):end) - y(15,(tstart-800):end)); % #3 A32
 %pause
 plot(x((tstart-800):end) ,y(20,(tstart-800):end) - y(21,(tstart-800):end)); % #4, A9
 %pause
 plot(x((tstart-800):end) ,y(26,(tstart-800):end) -   y(27,(tstart-800):end)); % #5, A-xx
 %pause  
 plot(x((tstart-800):end) ,y(32,(tstart-800):end) - y(33,(tstart-800):end)); % #6, A11
 xlabel('t (s)')
 %plot(x((tstart-800):end) ,Ip((tstart-800):end), 'g');  
 xlabel('t (s)')
 plot(x((tstart-800):end) ,Vyt((tstart-800):end), 'r--') % stimulii: wave
  legend('delta-y(1)', 'delta-y(2)', 'delta-y(3)', 'delta-y(4)', ...
      'delta-y(5)','delta-y(6)', 'stimulus-2' , 'wave' ); %text(3000, 0.3, 'wt(1-2) = 100')
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
 %plot(Ip((tstart-800):end), 'k');  % stimulii, on same scale
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
%   
% Group of 4 close together: net lfp output: sum: detla-y(2:5) 
  %lfp25= ( y(8,:) - y(9,:) +y(14,:) - y(15,:) ...
  %  + y(20,:) - y(21,:)+ y(26,:) - y(27,:) )/4; % av of 4 nodes
  %lfp_stats_2to5= [mean(lfp25(tstart:end))  (max(lfp25(tstart:end)) - min(lfp25(tstart:end)) )]
  %clear lfp25  lfp_stats_2to5
  
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
  xlim([0 2.0e5]); %axis([0,tlim, -Inf, Inf])
 subplot(5,1,5); 
 plot((Pd)*IpOn(2) +pran+Vyt, 'g--'); hold on; 
  plot((Pd2)*IpOnB(2), 'g--'); % incl 2nd pulse fn
  plot((Pd)*IpOn(6) +pran+Vyt, 'r--') % stimulii, on same scale
  plot((Pd2)*IpOnB(6), 'r--'); 
  xlim([0 2.0e5])
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

%% 5.0.1  mean and range of detla-y, After stim: 6 nodes  10:20sec 
% check dy3 before stim
 %dy3LT4 = max( y(14,1:35000) - y(15,1:35000) ) - min( y(14,1:35000) - y(15,1:35000) ) 
 %figure; plot ( y(14,1:45000) - y(15,1:45000) ); title(' NM # 3 - check state')
 % well after stim:
tstart = 100000;
fprintf('\n > dy(i) stats > 10 sec \n')
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
 d_yAmpl(5)= max(y(26,(tstart+2000):end) - y(27,(tstart+2000):end)) ...
      - min(y(26,(tstart+2000):end) - y(27,(tstart+2000):end) );
 d_yMean(6)=mean(y(32,(tstart+2000):end) - y(33,(tstart+2000):end)); % #6, A-11
 d_yAmpl(6)=max(y(32,(tstart+2000):end) - y(33,(tstart+2000):end)) ...
      -min(y(32,(tstart+2000):end) - y(33,(tstart+2000):end));
%
  figure; stem(d_yMean, '+');  hold on; stem(d_yAmpl, ':diamond');
  % hs(1).Marker='+'; hs(2).Marker='diamond'; 
 axis([0 7 -15 60 ]); xlabel('NM #'); ylabel('delta-y12(i)  (mV)'); title('dy(i) 10 : 20 sec')
  legend('mean', 'p-p ampl.', 'Location', 'northeast') % place out of the way
   % text(5.5, 23, 'Adj = 1');  % text(5.5, 10, 'alpha; C x vector')
 
 dy_mean_ampl = [d_yMean d_yAmpl] % test/ debug
 clear dy3LT4 dy_mean_ampl d_y* 
 
 %% 5.0.2 LFP LFP:  net lfp output: sum: detla-y(1) + detla-y(2), for the 6 NM
 fprintf('\ LFP, after stim: \n')
 lfp_after= (y(2,60001:end) - y(3,60001:end) + y(8,60001:end) - y(9,60001:end) +y(14,60001:end) - y(15,60001:end) ...
     + y(20,60001:end) - y(21,60001:end)+ y(26,60001:end) - y(27,60001:end) +y(32,60001:end) - y(33,60001:end) )/6; % av of 6 nodes
% figure; plot(lfp); title('JR model: net lfp output ')
  figure; plot(x(60001:end), lfp_after); title('6-star cluster, JR model: net lfp output, after stim (4 s) '); xlabel('t (sec) ')
  ylabel('LFT (t)  (mV)'); % hold on; plot(20,-1, 'b^'); plot(20,-1.1, 'b^'); % wave starts
lfp_mean=mean(lfp_after); lfp_ampl = max(lfp_after) - min(lfp_after);
lfp_mean_ampl = [lfp_mean lfp_ampl] % test/ debug  
  clear dy_mean* lfp_*

 %% 5.0.3 Save y vec for restarts
  % eg  for "states" a, b for ran resamples:  wb7, p191
  fprintf(' \n > save last 50 t-steps for sim restart \n')
clear y_end % to be sure
iend= length(lfp)-50+1; % want 50 elements
 y_end= y(1:36,iend:end);
 size(y_end) % debug    
  % y_6NMerfc_a = y_end; % Here NM#3 does osc. : full y is too big
csvwrite('y_6NMerfc_a.csv', y_end); % save for restart
clear y_end iend %y_6NMerfc_a

%% 5.0.3a LFP at end, for the 6 NM
 fprintf('\ LFP, last 50 steps: \n')
 istart= length(lfp)-50;
 lfp_end= (y(2,istart:end) - y(3,istart:end) + y(8,istart:end) - y(9,istart:end) +y(14,istart:end) - y(15,istart:end) ...
     + y(20,istart:end) - y(21,istart:end)+ y(26,istart:end) - y(27,istart:end) +y(32,istart:end) - y(33,istart:end) )/6; % av of 6 nodes
% figure; plot(lfp); title('JR model: net lfp output ')
  figure; plot(x(istart:end), lfp_end); title('6-star cluster, JR model: net lfp output, at end (20 s) '); xlabel('t (sec) ')
  ylabel('LFT (t)  (mV)'); % hold on; plot(20,-1, 'b^'); plot(20,-1.1, 'b^'); % wave starts
lfp_mean=mean(lfp_end); lfp_ampl = max(lfp_end) - min(lfp_end);
lfp_mean_ampl = [lfp_mean lfp_ampl] % test/ debug  
  clear dy_mean* lfp_* istart 

%% 5.2 Transitions, Pulse, Wave arrives:  PLOT all 6 nodes together:
tstart= 30001; % debug
 figure; plot(x((tstart-1000):end), y(2,(tstart-1000):end) - y(3,(tstart-1000):end), 'k', ...
     'LineWidth', 1.0);  hold on; xlim([3.5 5]); % node #1, A-10 
  %xlabel('t (sec)'); xlabel('t (ms)'); 
  ylabel('dy(t)  (mV)');  
 pause % - to closely examine indiv. waveforms
 plot(x((tstart-1000):end) ,(y(8,(tstart-1000):end) - y(9,(tstart-1000):end) ), 'b', ...
     'LineWidth', 1.0);   % #2, A-32V
 pause  % - to closely examine waveforms
 plot(x((tstart-1000):end) ,y(14,(tstart-1000):end) - y(15,(tstart-1000):end), ...
     'LineWidth', 0.75, 'Color', [0, 0.7, 0.7]); % #3 A-xx
 pause
 plot(x((tstart-1000):end) ,y(20,(tstart-1000):end) - y(21,(tstart-1000):end), ...
     'LineWidth', 0.75,'Color', [0.7, 0.7, 0]); % #4, A9
 pause
 plot(x((tstart-1000):end) ,y(26,(tstart-1000):end) -   y(27,(tstart-1000):end), 'm',...
     'LineWidth', 1.0); % #5, A46D
 pause
 plot(x((tstart-1000):end) ,y(32,(tstart-1000):end) - y(33,(tstart-1000):end), 'g', ...
     'LineWidth', 1.0); % #6, A11
 pause
  xlabel('t (sec)')
 hold off; hold on   % scan the [switched] 6 pulse streams stimulii
 plot(x((tstart-1000):end) , IpOn(1)*Pd(1, (tstart-1000):end)-20, 'k-', 'LineWidth', 0.65) % force color sequence
 %plot(x((tstart-1000):end) , IpOn(2)*Pd(2, (tstart-1000):end)-20, 'b-', 'LineWidth', 0.65) % highlight hubs
 %plot(x((tstart-1000):end) , IpOn(3)*Pd(3, (tstart-1000):end)-20, '-', 'LineWidth', 0.65, ...
 %    'Color', [0, 0.7, 0.7])
 %plot(x((tstart-1000):end) , IpOn(4)*Pd(4, (tstart-1000):end)-20, '-', 'LineWidth', 0.65, ...
 %    'Color', [0.7, 0.7, 0]) %
 %plot(x((tstart-1000):end), IpOn(5)*Pd(5, (tstart-1000):end)-20, 'm-', 'LineWidth', 0.65)
 %plot(x((tstart-1000):end), IpOn(6)*Pd(6, (tstart-1000):end)-20, 'g-', 'LineWidth', 0.65) 
 %plot(x((tstart-1000):end) ,Vyt((tstart-1000):end), 'r--', 'LineWidth', 1.0) % stimulii: wave 
 %pause 
 legend('delta-y(1)', 'delta-y(2)', 'delta-y(3)', 'delta-y(4)', ...
      'delta-y(5)','delta-y(6)', 'Pulses-1', ' wave'); % ...
        %'Pulses-2','Pulses-3', 'Pulses-4', 'Pulses-5', 'Pulses-6', ' wave'); %text(3000, 0.3, 'wt(1-2) = 100')
  % axis([0 5000 -20 20]);  %axis([1 2 -25 25]); xlim([4 4.2]);
  xlim([3.5 5]); ylim([-75 100]); % xlim([3.9 5]); % short wave %xlim([4.5 6]); % early
  % xlim([8.2 9.5]);  % later
  % xlim([7.8 9.5]);  % covers 2nd NN volley
 title('6 NM delta-y(i; t) - random Adj, noise')
  %title('6 NM delta-y(i; t) - unit Adj, noise')
  % title('8 NM delta-y(i; t) - V2 pulses arrive')  % pulse arrive 
  % xlim([13 16]);  
  %title('6NM delta-y(j; t) - wave arrives') % wave arrives
   % text(4.1,35, '10mV, 2.8Hz wave x 1.5s') % label, info
 
%% 4.1 Align waveforms, by phase: 1st Pulse
% a) at NM#2 [the input]
dy2_sample=y(9, 40000:60000) - y(8, 40000:60000); % NM#2, delta-y 
figure; plot(dy2_sample); hold on; xlim([0 2000]) %axis([0.5 2.5 -50 50])
dym=mean(dy2_sample); plot([0.5 2.5], [dym dym], 'g--'); xlabel('t (ms)')
%
title('JR node #2, 1t pulse: output delta-y: y1(2) - y2(2)') % e - i 
plot(Pd(40000:60000)/10, 'k') % stim % reduce scale [of pulse rate]
plot(Ip(40000:60000), 'k') ; plot(5*Vyt(40000:60000)+dym, 'Color', [1 0 1]) % stim V pulses, wave
legend('delta-y(2)', 'mean', 'Pulse', 'Spikes', 'wave'); 

% a.1)
dy6_sample=y(32, 40000:60000) - y(33, 40000:60000); % NM#6, delta-y 
figure; plot(dy6_sample); hold on; xlim([0 2000]) %axis([0.5 2.5 -50 50])
dym=mean(dy6_sample); plot([0.5 2.5], [dym dym], 'g--'); xlabel('t (ms)')
title('JR node #6: output delta-y: y1(6) - y2(6)') % e - i 
plot(Pd(40000:60000)/10, 'k') % stim % reduce scale [of pulse rate]
plot(Ip(40000:60000), 'k') ; plot(5*Vyt(40000:60000)+dym, 'Color', [1 0 1]) % stim V pulses, wave
legend('delta-y(6)', 'mean', 'Pulse', 'Spikes', 'wave')

% b) lfp [the output]
%lfp_sample=lfp(1:1200); % all NM#2, nett lfp
lfp_sample=lfp(40000:60000);
figure; plot( lfp_sample); hold on; xlim([0 2000]) %axis([600 1200 -5 45])
dym=mean(lfp_sample); plot([3.500 6.0], [dym dym], 'g--')
title('JR 6 NM: align lfp w stim') 
plot(Pd(40000:60000)/10, 'k') % stim % reduce scale
plot(Ip(40000:60000), 'k');  %plot(Ip(1:1200), 'k'); % 
plot(5*Vyt(40000:60000)+dym, 'Color', [1 0 1]) % wave; align axes
legend('lfp', 'Pulse', 'wave'); %axis([0.5 2.5 -10 30]); 
xlabel('t steps (0.1ms)')

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

  %% 4.1.2 ISI.  {was 3.6}  ISI for y-e(2);  dy(2) / y-e(1);  dy(1)
   % fprintf('\n > ISI for y-i(1):  '); 
   %fprintf('\n > ISI for y-e(1): from 3.0 sec '); 
 %fprintf('\n > ISI for dy(1):  from 3.0 sec ');
 fprintf('\n > ISI for LFP (all NM, t): from 3.0 sec ');
 %signal= y(2, 40000:end);  % y-e(1)  <:: edit, to pick signal
% signal= y(2, :);  % y-e(1) all
   % signal= y(2, 30000:100000); 
% signal= y(8, 30000:end);     % y-e(2) <:: long time
  % signal= y(8, 35000:100000);  % y-e(2)
   % signal= y(8, :); % check all t
   %signal= y(32, 30000:100000);  % y-e(6)
  % signal= y(32, 30000:end);      % <<::
% signal= y(3, 30000:100000);  % y-i(1)
 % signal= y(3, 35000:end); 
   % signal= y(3, :); % check all t
% signal= y(9, 35000:100000);  % y-i(2)
% signal= y(33, 35000:100000);  % y-i(6)

% signal= y(2, 40000:end)-y(3, 40000:end);  % <:: dy(1), for NM#1 
  % signal= y(2,:)-y(3, :);  % <:: all dy(1), for NM#1
    %signal= y(2, 1000:end)-y(3, 1000:end);  % dy(1), for 1/2NM
% signal= y(8,30000:end)-y(9,30000:end);  % dy(2), for NM#2
  % signal= y(32,35000:100000)-y(33,35000:100000);  % dy(6), for NM#6

 signal= lfp(30000:end);  % <:: lfp(t), av over all NM

%  % these are peaks of V(t) - max
[SamplePeaks, PeakLocs]= findpeaks(signal, ...
    'MinPeakDistance',500 ,'MinPeakheight',-10.0); % for alpha band: 1000; beta 500 / 8.0
% length(PeakLocs)
PeakTimes=PeakLocs*h;  %(1:5)
ISI = diff(PeakLocs)*h;
 ISI  %ISI(1:20)
 Range = [min(diff(PeakLocs)*h) max(diff(PeakLocs)*h)]
 figure; hist(ISI); title('6NM, ISI istribution'); xlabel(' ISI (sec) ')
 figure; hist(1./ISI); title('6NM, equiv freq. distribution'); xlabel(' f (Hz) ')
 
% 3 plots together:
 figwidth = 560; figheight = 704; figposition = [500, 100, figwidth, figheight]; % tall
figure('position',figposition, 'units','pixels');  % larger figure; 
subplot(3,1,1); stem(ISI); % xlim([0 50]); %ylim([0.1 0.2]);  % for 6NM, alpha
 ylabel('ISI (sec)'); title('6MN cluster; node #1: ISI times (sec)')
subplot(3,1,2); stem(ISI);  xlabel('peak # (from 3.0 sec) '); % ylim([0.11 0.22]); % alpha 
ylim([0.07 0.26]); % beta
  grid on; ylabel('ISI (sec)'); %xlim([40 length(ISI)]); %xlim([80 160]); %
  % text(7, 0.205, 'v', 'Color', 'r')
 % subplot(3,1,1); text(20, 0.177, 'NM #1, y-e(1);  10 pulses @ 6.35Hz')    % '5 pulses @ 3Hz')
subplot(3,1,3); hold on; title('pulse train'); % text(2.5, 0.25, '6 NM, 0.1mV U ran noise')
 t0 = PeakLocs(1)*h; tt=t0; % (sec)
for i =1:length(PeakTimes)-1
   line([tt tt], [0.0 0.2 ], 'LineWidth', 1.0 );
   tt=tt+ISI(i); % accumulate ISI's to get this time
end
grid on; xlabel(' t (sec +3.0)'); ylim([0 0.3]); xlim([0 4]);  
%subplot(3,1,3); text(1.0, 0.21, 'v', 'Color', 'r')  % mark stimulus time
%                                                   [@3.5 + 0.5 s;  or  @1.0+3.0 s]
%text(3.0, 0.21, 'v', 'Color', 'r')

% show ISI(t) & Pulse  Train together
figwidth = 1300; figheight = 350; figposition = [100, 380, figwidth, figheight]; % wide
figure('position',figposition, 'units','pixels');  % larger figure; 
subplot(2,1,1); hold on; %stem(ISI);
ylabel('ISI (sec)'); title('6MN cluster; node #1:  ISI times (sec)')
 t0 = PeakLocs(1)*h; tt=t0; % (sec)
for i =1:length(PeakTimes)-1
   line([tt tt], [0.0 ISI(i) ], 'LineWidth', 1.0 );
   plot(tt, ISI(i), 'ob')
   tt=tt+ISI(i); % accumulate ISI's to get this time
end
text(1.0, 0.5, 'v', 'Color', 'r')  % mark stimulus time, ie 4s (nb. t -3s)
 %text(6.0, 0.5, 'v', 'Color', 'r') 
 subplot(2,1,2); hold on; title('pulse train'); % text(2.5, 0.25, '6 NM, 0.1mV U ran noise')
 t0 = PeakLocs(1)*h; tt=t0; % (sec)
for i =1:length(PeakTimes)-1
   line([tt tt], [0.0 0.2 ], 'LineWidth', 1.0, 'Color','k' );
   tt=tt+ISI(i); % accumulate ISI's to get this time
end
grid on; xlabel(' t (sec +3.0)', 'FontWeight', 'bold'); ylim([0 0.3]); %xlim([0 4]);  
text(1.0, 0.21, 'v', 'Color', 'r')  % mark stimulus time
% text(6.0, 0.22, 'v', 'Color', 'r') % marker: wave starts
 % text(6.0, 0.235, 'wave arrives', 'Color', 'r') 
 % subplot(2,1,2); text(1.0, 0.24, '+1 Pulse -> #2]', 'Color', 'r')  % describe stimulu2
%
 figure; plot(signal); hold on; %title('6NM: y-e(1), dy(1);  ISI')
 title('signal V(t)'); xlabel('t (3s + 0.1ms)'); ylabel('V(t)  (mV) ');
plot(PeakLocs, signal(PeakLocs), 'ro'); grid on; % mark peaks - to be sure all id'd
%
Npeaks=length(PeakLocs)
figure; hold on
for i=1:Npeaks-1
plot(ISI(i), signal(PeakLocs(i)), 'ko');
 %drawnow
 %pause(0.007) % animate: jumps L/R, U/D
end
grid on; title('6MN; node #1: y-e(1) vs ISI'); xlabel(' ISI (sec)'); ylabel('y-e(1)(t) (mV) ');
 
 shortT = mean(diff(PeakLocs))*h  % short period
 StdDev = std(diff(PeakLocs))*h
 AvshortT_freq= [shortT 1/shortT]

% show ISI(t) & V(t) Signal together
figwidth = 1300; figheight = 350; figposition = [100, 380, figwidth, figheight]; % wide
figure('position',figposition, 'units','pixels');  % larger figure; 
subplot(2,1,1); hold on; xlim([0 17]); %stem(ISI);
ylabel('ISI (sec)'); title('6MN cluster; node #1:  ISI times (sec)')
text(1.0, 0.25, 'v', 'Color', 'r')  % mark stimulus time
 t0 = PeakLocs(1)*h; tt=t0; % (sec)
for i =1:length(PeakTimes)-1
   line([tt tt], [0.0 ISI(i) ], 'LineWidth', 1.0 );
   plot(tt, ISI(i), 'ob')
   tt=tt+ISI(i); % accumulate ISI's to get this time
end
text(1.0, 0.5, 'v', 'Color', 'r'); grid on  % mark stimulus time, ie 4s (nb. t -3s)
 %text(6.0, 0.5, 'v', 'Color', 'r')
 
 subplot(2,1,2); hold on; title('V waveform'); % text(2.5, 0.25, '6 NM, 0.1mV U ran noise')
 plot(x(30000:end)-4.0 ,signal); hold on; xlim([0 17]); xlabel(' t (sec +3.0)'); % nb. adjust t axis
 grid on
for i =2:length(PeakLocs)
   plot(x(PeakLocs(i)), signal(PeakLocs(i)), 'or') % highlight peaks
end
 
% save "code" in 10 ms bins, as "binary" array - for later classification
%Code=zeros(1, int32(length(x)/100)); % now 3k long, for 30 sec [not 300k]
%tmp=int32(PeakLocs/100); % place in bins of 100*(0.1 ms) long
%Code(tmp)=1;
 % Code_3Pulses=Code; csvwrite('Code_3Pulses.csv', Code_3Pulses); % save
 % ISI_3Pulses=ISI; csvwrite('ISI_3Pulses.csv', ISI_3Pulses); % for analysis

 clear tmp Code to Peak* Range Npeaks SamplePeaks shortT StdDev AvshortT_freq tt signal

 %% show "f"(t) & Pulse  Train together
figwidth = 1300; figheight = 350; figposition = [100, 380, figwidth, figheight]; % wide
figure('position',figposition, 'units','pixels');  % larger figure; 
subplot(2,1,1); hold on; %stem(ISI);
ylabel('f (Hz)'); title('6MN cluster; node #1:  f (Hz) - local to each spike')
 t0 = PeakLocs(1)*h; tt=t0; % (sec)
for i =1:length(PeakTimes)-1
   line([tt tt], [0.0 1/ISI(i) ], 'LineWidth', 1.0 );
   plot(tt, 1/ISI(i), 'ob')
   tt=tt+ISI(i); % accumulate ISI's to get this time
end
%text(1.0, 0.5, 'v', 'Color', 'r')  % mark stimulus time (nb. t -3s)
 text(4.05, 11, 'v', 'Color', 'r') 
 subplot(2,1,2); hold on; title('pulse train'); % text(2.5, 0.25, '6 NM, 0.1mV U ran noise')
 t0 = PeakLocs(1)*h; tt=t0; % (sec)
for i =1:length(PeakTimes)-1
   line([tt tt], [0.0 0.2 ], 'LineWidth', 1.0, 'Color','k' );
   tt=tt+ISI(i); % accumulate ISI's to get this time
end
grid on; xlabel(' t (sec +3.0)', 'FontWeight', 'bold'); ylim([0 0.3]); %xlim([0 4]);  
  %text(1.0, 0.21, 'v', 'Color', 'r')  % mark stimulus time
% text(6.0, 0.22, 'v', 'Color', 'r') % marker: wave starts
 % text(6.0, 0.235, 'wave arrives', 'Color', 'r') 
 %
 
 %% 4.2 Phase plot (3 y's for  NM #1)

figure;  subplot(3,1,1); plot(y(1,:), y(4,:)); 
 %title('dy0/dt vs y0(t);   zero input, C= 135, tau(15, 20), 1e-12 ran noise'); hold on
title('dy0/dt vs y0(t); C= 135, tau(15, 20),  0.1mV noise, stim = 300 '); hold on

ylabel('dy0/dt' ); xlabel('y0(t)')
 subplot(3,1,2); plot(y(2,:), y(5,:)); title ('dy1/dt vs y1(t)')
   ylabel('dy1/dt' ); xlabel('y1(t)')
 subplot(3,1,3); plot(y(3,:), y(6,:)); title ('dy2/dt vs y2(t)')
   ylabel('dy2/dt' ); xlabel('y2(t)'); hold off

%% 4.2a  Phase plot (3 y's for  NM #2 )

figure;  subplot(3,1,1); plot(y(7,:), y(10,:)); 
 %title('dy0/dt vs y0(t);   zero input, C= 135, tau(15, 20), 1e-12 ran noise'); hold on
title('NM#2:  dy0/dt vs y0(t); C= 300, tau(10, 20),  0.1mV G noise, S[v-wt] '); hold on

ylabel('dy0/dt' ); xlabel('y0(t)')
 subplot(3,1,2); plot(y(8,:), y(11,:)); title ('dy1/dt vs y1(t)')
   ylabel('dy1/dt' ); xlabel('y1(t)'); 
 subplot(3,1,3); plot(y(9,:), y(12,:)); title ('dy2/dt vs y2(t)')
   ylabel('dy2/dt' ); xlabel('y2(t)'); hold off
   
%% 4.2.1 state space plot (3 y's for 1 NM)

figure; plot(y(2,:), y(3,:)); 
title('y1(t) vs y2(t); step stim, C= 300, tau(10, 10), 0.1mV noise  '); hold on
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

  %% 4.2d Phase plot y1'(t) vs y1(t) {ie y-e}  for NM#1 : animated
figure; hold on; title ('y-e (1):  dy1/dt(t) vs y1(t)')
%title ('dy1/dt vs y1(t)') %axis([10 60 10 60]); 
ylabel('dy1/dt(t)' ); xlabel('y1(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ -0 40 -600 1000]); % initiation
%for i = 1:5:5000  % nb transients < 2s
 for i = 5000:5:18000  % at s/s
   %for i = 10000:5:15000  % etc
  %for i = 100000:5:120000  % after 10 s [transition]
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
 plot(y(2,i), y(5,i), 'k.', 'MarkerSize', 10);
 pause(0.0005); drawnow
end
fprintf('\n  plot done!  ') 

 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
%ylabel('dy1/dt' ); xlabel('y1(t)');
hold off

%% 4.2e Phase plot y2'(t) vs y2(t) {ie y-i}  for NM#1 : animated
figure; hold on; title ('y-i (1):  dy2/dt(t) vs y2(t)')
ylabel('dy2/dt(t)' ); xlabel('y21(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ 0 40 -500 1000]); % initiation
%for i = 1:5:5000 % startup transients
  for i = 5000:5:10000  % at s/s 
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(y(3,i), y(6,i), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fprintf('\n  done!  ') 

 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
hold off

%% 4.2f Phase plot y1'(t) vs y1(t) {ie y-e}  for NM#2 : animated
figure; hold on; title ('y-e (2):  dy1/dt(t) vs y1(t)')
%title ('dy1/dt vs y1(t)') %axis([10 60 10 60]); 
ylabel('dy1/dt(t)' ); xlabel('y1(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ -0 40 -600 1000]); % initiation
%for i = 1:5:5000  % nb transients < 2s
 for i = 5000:5:10000  % at s/s 
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(y(8,i), y(11,i), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fprintf('\n  plot done!  ') 

 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
%ylabel('dy1/dt' ); xlabel('y1(t)');
hold off

%% 4.2g Phase plot y2'(t) vs y2(t) {ie y-i}  for NM#2 : animated
figure; hold on; title ('y-i (2):  dy2/dt(t) vs y2(t)')
ylabel('dy2/dt(t)' ); xlabel('y21(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ 0 40 -500 1000]); % initiation
for i = 1:5:5000
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(y(9,i), y(12,i), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fprintf('\n  done!  ') 

 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
hold off

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

%% 4.4 Phase plot LFP(all): av[y1-y2)'(t)] vs av[y1(t)-y2(t)] {ie y-e - y-i} for NM#1-6 : animated
fprintf('\n  > LFP[e-i] av. for all NM  \n')
 % LFP LFP:  net lfp output: sum: detla-y(1) + detla-y(2), etc
 lfp= (y(2,:) - y(3,:) + y(8,:) - y(9,:) +y(14,:) - y(15,:) ... % calcabove
     + y(20,:) - y(21,:) + y(26,:) - y(27,:) +y(32,:) - y(33,:) )/6; % av of 6 nodes
 lfpdot= (y(5,:) - y(6,:) + y(11,:) - y(12,:) +y(17,:) - y(18,:) ...
     + y(23,:) - y(24,:) + y(29,:) - y(30,:) +y(35,:) - y(36,:) )/6;
figure; hold on; title ('LFP (all):  dy1-2/dt(t) vs y1-2(t)')
ylabel('dy1-2/dt(t)' ); xlabel('y1-2(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ 0 40 -500 1000]); % initiation
for i = 35000:5:40000 % time windows : before stim
plot(lfp(i), lfpdot(i), 'g.', 'MarkerSize', 10); %  nb order: plot (x,y
pause(0.0005)
drawnow
end
fp1=mean(lfpdot(35000:40000)); fp2=mean(lfp(35000:40000)); % before stim
plot(fp2, fp1, 'g.', 'MarkerSize', 25);
FixedPoint_before=[fp1 fp2]
% stim arrives at 4 s
for i = 40000:5:45000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(lfp(i), lfpdot(i), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
% another color, for later
for i = 45000:5:50000 % time windows
 %for i = 60000:5:65000 % later
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(  lfp(i), lfpdot(i), 'b.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
%fp1=mean(lfpdot(40000:end)); fp2=mean(lfp(40000:end)); % overall, after stim
fp1=mean(lfpdot(40000:50000)); fp2=mean(lfp(40000:50000)); % just for plotted orbit
plot(fp2, fp1, 'b.', 'MarkerSize', 25);
  %text(57, -2100, '14 to 16  sec'); %text(27, 12, '0.1 sec')
 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
 % text(5, 1200, '4:5 sec')
hold off
FixedPoint_after=[fp1 fp2]
fprintf('\n  av:LFP [e-i] done!  \n')   
 clear fp* FixedPoint*  %lfpdot 
 
 %% 4.4.1 add later orbit
 hold on
 % another color, for later
for i = 110000:5:120000 % time windows
 %for i = 60000:5:65000 % later
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot( lfp(i), lfpdot(i), 'r.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fp1=mean(lfpdot(110000:120000)); %fp2=mean(lfp(40000:end)); % overall, after stim
fp2=mean(lfp(110000:120000)); % just for plotted orbit
plot(fp2, fp1, 'r.', 'MarkerSize', 25);
FixedPoint_ss=[fp1 fp2] % later, in beating phase
fprintf('\n  extra done!  \n')
hold off
clear fp* FixedPoint*  lfpdot

%% Appx. 1 coords of NM (2D xz)
 addpath(genpath('/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain')); % Data files;  include sub-directories 
fprintf('  use marmoset coords from Atlas; cluster in x-z plane \n')
clusterlist=[2, 26, 25, 48, 32, 3, 45, 47]; % as plotted in x-y plane, viewed from Ant.
 % from Marm..ReadDataV1.m #1.1b  coords - prev calc from Atlas volumes; cf Appx. A.1 below
  % from 3D vol: 'atlas_segmentation.nii'
NodeCoord=csvread('CoordsMarmoset.csv'); % 116x3 [x y z]; nb. some absent: NaN 
ClusterCoord= NodeCoord(clusterlist,:);
 % diff in xz plane
dXZ=zeros(5,4);  % [x z r theta] for ea point [coord of NM]; r wrt NM(1): A10
dXZ(:,1)= ClusterCoord(2:6,1) - ClusterCoord(1,1); % dx(1-j)
dXZ(:,2)= ClusterCoord(2:6,3) - ClusterCoord(1,3);  % dz(1-j)
dXZ(:,3)= sqrt(dXZ(:,1).^2 + dXZ(:,2).^2 ); % r coord
 %dXZ(:,4)= atan(dXZ(:,2)./dXZ(:,1)); % rad, simple form
dXZ(:,4)= atan2(dXZ(:,2), dXZ(:,1)); % rad, 4 quadrant form; theta in xz plane
 %dXZ(:,1)= -dXZ(:,1); % flip x, for view from Antr (front)
   % dXZ(:,4)*(180/pi) % debug
        % -46.4075, 7.9610, 60.0953, 115.7795,  -138.1053
%
figure; plot(dXZ(:,1), dXZ(:,2), '.', 'MarkerSize', 25); hold on
plot(0, 0,'bo', 'MarkerSize', 12); % A10, central
xlabel('dx (mm)  M-L'); ylabel('dz (mm)  V-D');
line([ 0 1.5], [0 0]);  line([0 0], [0 2]); % x z  axis
%
text( (0+0.05), (0+0.1), num2str(1),'FontSize',9); % central pt, A10
for i=1:5
 text( (dXZ(i,1)+0.05), (dXZ(i,2)+0.1), num2str(i+1),'FontSize',9); % number pts, shift slightly
end

%% Appx 1.1 plot dy in radial vectors (2D xz)
 % @ ea time step
 % dV calc at #5.1 above
  %[min(dV(1,:))  max(dV(1,:))] % scale
  lcolr = char(5); lcolr(1)='b'; lcolr(2)='c'; lcolr(3)='g'; lcolr(4)='m'; lcolr(5)='r';
  thetary = [];
 figure; hold on; %
 axis([-100 100 -100 100]); xlabel('dV * dx (mV mm)  M-L'); ylabel('dV * dz (mV mm)  V-D');
 title('6 NM dipoles; vector sum vs t')
  %axis([-50 50 -50 75]);  %axis([-10 10 -10 10]); 
    %plot(dV(1,:),'b');  ylabel(' dV(1-2)  (mV)')  % 1-2, A10 - A32v: 
 % polar coords
 %for it=1:100:20000 %length(x) % tsteps (0.1 ms)
 for it=35000:100:50000 % before/after stim @ 4s
  xAv=0; yAv=0; % form net vector sum of ea dV; restet at ea t step
  for i = 1:5 % which dV
   rp =dXZ(i,3); thetap = dXZ(i,4); %rp =dXZ(i,6); % [r, theta] coords
   if dV(1,it) < 0
      thetap =thetap + pi; % "neg" direction for dipole
   end % if
   thisx=rp*cos(thetap)*dV(i,it); thisy=rp*sin(thetap)*dV(i,it);
   % [rp*cos(thetap) rp*sin(thetap)]; [rp*cos(thetap)*dV(1,i) rp*sin(thetap)*dV(1,i)] % debug
   hline(i)=line([0 thisx], [0 thisy], 'color', lcolr(i),'LineWidth',1);
   %line([0 rp*cos(thetap)*dV(1,it)], [0 rp*sin(thetap)*dV(1,it)], 'LineWidth',1) % ??cycle thru default colors
   drawnow
   xAv=xAv+thisx; yAv=yAv+thisy; % accumulate vectors
  end % scan nodes (i)
   thetaAv=atan2(yAv, xAv); % rad, 4 quadrant form
   thetary = [thetary, thetaAv]; % acccumulate net circ angle(t)
  hline0=line([0 xAv], [0 yAv], 'color', 'k', 'LineWidth',1.5);
  timestamp=num2str(it*h, 4); timestamp=horzcat(timestamp, ' sec'); % dt is h
  htext =text(60, 80, timestamp);
  pause(0.1)
  delete(hline0);  delete(hline); delete(htext); % so can draw a fresh line & txt
  patch([0; xAv; 0], [0; yAv; 0], 'k', 'EdgeColor', ...
               [0.1 0.1 0.1], 'LineSmoothing','on', 'EdgeAlpha', 0.2, 'FaceColor', 'none'); % accum gray trace
 end % scan t steps (it)
 line([0 xAv], [0 yAv], 'color', 'k', 'LineWidth',1.5); % leave last one
 htext =text(60, 80, timestamp); % final time
 FinishAt =timestamp
 
 % gather circular vector [angle] path {no graphics}
 clear thetary; thetary = [];
 for it=1:length(dV) % use _all_ tsteps (0.1 ms)
  xAv=0; yAv=0; % form net vector sum of ea dV; restet at ea t step
  for i = 1:5 % which dV
   rp =dXZ(i,3); thetap = dXZ(i,4); % [r, theta] coords
   if dV(1,it) < 0
      thetap =thetap + pi; % "neg" direction for dipole
   end % if
   thisx=rp*cos(thetap)*dV(i,it); thisy=rp*sin(thetap)*dV(i,it);
   xAv=xAv+thisx; yAv=yAv+thisy; % accumulate vectors
  end % scan nodes (i)
   thetaAv=atan2(yAv, xAv); % rad, 4 quadrant form
   thetary = [thetary, thetaAv]; % acccumulate net circ angle(t)
 end % scan t steps (it)
 
 figure; plot(thetary); xlabel('t step (0.1 ms)'); ylabel('dipole angle(t) (rad)')
 title('6 NM, net dipole angle') % xlim([0 2.0e4])
 
  clear lcolr xAv yAv rp thetap thisx thisy hline* htext timestamp thetaAv FinishAt
 
%% Appx. 1.2 coords of NM (3D xyz)
 addpath(genpath('/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain')); % Data files;  include sub-directories 
fprintf('  use marmoset coords from Atlas; cluster in x-z plane \n')
clusterlist=[2, 26, 25, 48, 32, 3]; %, 45, 47]; % 6, 8; as plotted in x-y plane, viewed from Ant.
 % from Marm..ReadDataV1.m #1.1b  coords - prev calc from Atlas volumes; cf Appx. A.1 below
  % from 3D vol: 'atlas_segmentation.nii'
NodeCoord=csvread('CoordsMarmoset.csv'); % 116x3 [x y z]; nb. some absent: NaN 
ClusterCoord= NodeCoord(clusterlist,:);
 % diff in xz plane
dXYZ=zeros(5,6);  % [x y z r theta] for ea point [coord of NM]; r wrt NM(1): A10
dXYZ(:,1)= ClusterCoord(2:6,1) - ClusterCoord(1,1); % dx(1-j)
dXYZ(:,2)= ClusterCoord(2:6,2) - ClusterCoord(1,2); % dy(1-j)
dXYZ(:,3)= ClusterCoord(2:6,3) - ClusterCoord(1,3); % dz(1-j)
[dXYZ(:,4), dXZ(:,5), dXZ(:,6)] = cart2sph(dXYZ(:,1), dXYZ(:,2), dXYZ(:,3)); % order: [az, elev, r]
 % dXZ(:,2)= -dXZ(:,2); % flip y; point to Antr
  % az was theta, above - in xz plane
 %dXZ(:,1)= -dXZ(:,1); % flip x, for view from Antr (front)
   % dXZ(:,4)*(180/pi) % debug
        % -46.4075, 7.9610, 60.0953, 115.7795,  -138.1053
%
figure; plot3(dXYZ(:,1), dXYZ(:,2), dXYZ(:,3),'.', 'MarkerSize', 25); hold on
plot(0, 0,'bo', 'MarkerSize', 12); hold on % A10, central
plot(0, 0,'bo', 'MarkerSize', 6);
xlabel('dx (mm)  M-L'); ylabel('dy (mm)  A-P'); zlabel('dz (mm)  V-D');
line([ 0 1.5], [0 0], [0 0]); line([0 0], [0 0], [0 2]); % x z  axis
grid on; view(245, 30);
%
text( (0+0.05), (0+0.05), (0+0.1), num2str(1),'FontSize',9); % central pt, A10
for i=1:5
 text( (dXYZ(i,1)+0.05), (dXYZ(i,2)+0.05), (dXYZ(i,3)+0.1), num2str(i+1),'FontSize',9); % number pts, shift slightly
end 

%% Appx 1.3 plot dV(1,j) in radial vectors (3D xyz), 2D cut plane
 % @ ea time step
 % dV calc at #5.1 above
  %[min(dV(1,:))  max(dV(1,:))] % scale
  lcolr = char(5); lcolr(1)='b'; lcolr(2)='c'; lcolr(3)='g'; lcolr(4)='m'; lcolr(5)='r';
  azary = []; elevary = [];
 figure; hold on; axis([-150 200 10 100]); grid on; % 2D
 %axis([-100 100 -100 100 0 2e3]); grid on; view(245, 30); % 3D
  xlabel('dV * dx (mV mm)  M-L'); ylabel('dV * dy (mV mm)  A-P');
   %zlabel('dV * dz (mV mm)  V-D'); 
  title('6 NM dipoles in horiz. x-y plane; vector sum vs t')
  %axis([-50 50 -50 75]);  %axis([-10 10 -10 10]); 
    %plot(dV(1,:),'b');  ylabel(' dV(1-2)  (mV)')  % 1-2, A10 - A32v: 
 % polar coords
 %for it=1:100:20000 %length(x) % tsteps (0.1 ms)
 for it=1:100:20000 % before/after stim @ 4s
  xAv=0; yAv=0; zAv=0;% form net vector sum of ea dV; reset at ea t step
  for i = 1:5 % which dV
   azp =dXZ(i,4); elevp =dXZ(i,4); rp=dXZ(i,6); % [az, elev(wasTheta), r] coords
   if dV(1,it) < 0
      azp =azp + pi; elevp =elevp + pi; % "neg" direction for dipole
   end % if
   [thisx, thisy, thisz] = sph2cart(azp, elevp, rp*dV(i,it)); % rescale radius by dV
    % debug
   hline(i)=line([0 thisx], [0 thisy], 'color', lcolr(i),'LineWidth',1); % in 2D, at this t
   %line([0 rp*cos(thetap)*dV(1,it)], [0 rp*sin(thetap)*dV(1,it)], 'LineWidth',1) % ??cycle thru default colors
   drawnow
   xAv=xAv+thisx; yAv=yAv+thisy; zAv=zAv+thisz; % accumulate vectors
  end % scan nodes (i)
   [azAv elevAv, rAv] =cart2sph(xAv, yAv, zAv);
   azary = [azary, azAv]; elevary = [elevary, elevAv]; % acccumulate net circ angle(t)
  hline0=line([0 xAv], [0 yAv], [0 zAv], 'color', 'k', 'LineWidth',1.5);
  timestamp=num2str(it*h, 4); timestamp=horzcat(timestamp, ' sec'); % dt is h
  htext =text(60, 150, timestamp);
  pause(0.1)
  delete(hline0);  delete(hline); delete(htext); % so can draw a fresh line & txt
  patch([0; xAv; 0], [0; yAv; 0], 'k', 'EdgeColor', ...
               [0.1 0.1 0.1], 'LineSmoothing','on', 'EdgeAlpha', 0.2, 'FaceColor', 'none'); % accum gray trace
 end % scan t steps (it)
 line([0 xAv], [0 yAv],  'color', 'k', 'LineWidth',1.5); % leave last one
 htext =text(60, 150, timestamp); % final time
 FinishAt =timestamp
 
 %% gather circular vector [angle] path {no graphics}
 clear thetary; thetary = [];
 for it=1:length(dV) % use _all_ tsteps (0.1 ms)
  xAv=0; yAv=0; % form net vector sum of ea dV; restet at ea t step
  for i = 1:5 % which dV
   rp =dXZ(i,3); thetap = dXZ(i,4); % [r, theta] coords
   if dV(1,it) < 0
      thetap =thetap + pi; % "neg" direction for dipole
   end % if
   thisx=rp*cos(thetap)*dV(i,it); thisy=rp*sin(thetap)*dV(i,it);
   xAv=xAv+thisx; yAv=yAv+thisy; % accumulate vectors
  end % scan nodes (i)
     thetaAv=atan2(yAv, xAv); % rad, 4 quadrant form
   thetary = [thetary, thetaAv]; % acccumulate net circ angle(t)
 end % scan t steps (it)
 
 figure; plot(thetary); xlabel('t step (0.1 ms)'); ylabel('dipole angle(t) (rad)')
 title('6 NM, net dipole angle') % xlim([0 2.0e4])
 
  clear lcolr xAv yAv rp thetap thisx thisy hline* htext timestamp thetaAv FinishAt
  
 %% Appx 1.3a plot dV(1,j) in radial vectors (3D xyz)
 % @ ea time step
 % dV calc at #5.1 above
  %[min(dV(1,:))  max(dV(1,:))] % scale
  lcolr = char(5); lcolr(1)='b'; lcolr(2)='c'; lcolr(3)='g'; lcolr(4)='m'; lcolr(5)='r';
  azary = []; elevary = [];
 figure; hold on; %axis([-150 200 10 100]); grid on; % 2D
 axis([-150 200 -10 100 -100 100]); 
 grid on; view(20, 30); % 3D
  xlabel('dV * dx (mV mm)  M-L'); ylabel('dV * dy (mV mm)  A-P');
  zlabel('dV * dz (mV mm)  V-D'); 
  title('6 NM dipoles in horiz. x-y plane; vector sum vs t')
  %axis([-50 50 -50 75]);  %axis([-10 10 -10 10]); 
    %plot(dV(1,:),'b');  ylabel(' dV(1-2)  (mV)')  % 1-2, A10 - A32v: 
 % polar coords
 %for it=1:100:20000 %length(x) % tsteps (0.1 ms)
 for it=35000:100:50000 % before/after stim @ 4s
  xAv=0; yAv=0; zAv=0;% form net vector sum of ea dV; reset at ea t step
  for i = 1:5 % which dV
   azp =dXZ(i,4); elevp =dXZ(i,4); rp=dXZ(i,6); % [az, elev(wasTheta), r] coords
   if dV(1,it) < 0
      azp =azp + pi; elevp =elevp + pi; % "neg" direction for dipole
   end % if
   [thisx, thisy, thisz] = sph2cart(azp, elevp, rp*dV(i,it)); % rescale radius by dV
    % debug
   hline(i)=line([0 thisx], [0 thisy], 'color', lcolr(i),'LineWidth',1); % in 2D, at this t
   %line([0 rp*cos(thetap)*dV(1,it)], [0 rp*sin(thetap)*dV(1,it)], 'LineWidth',1) % ??cycle thru default colors
   drawnow
   xAv=xAv+thisx; yAv=yAv+thisy; zAv=zAv+thisz; % accumulate vectors
  end % scan nodes (i)
   [azAv elevAv, rAv] =cart2sph(xAv, yAv, zAv);
   azary = [azary, azAv]; elevary = [elevary, elevAv]; % acccumulate net circ angle(t)
  hline0=line([0 xAv], [0 yAv], [0 zAv], 'color', 'k', 'LineWidth',1.5);
  timestamp=num2str(it*h, 4); timestamp=horzcat(timestamp, ' sec'); % dt is h
  htext =text(60, 70, 70, timestamp);
  pause(0.1)
  delete(hline0);  delete(hline); delete(htext); % so can draw a fresh line & txt
  patch([0; xAv; 0], [0; yAv; 0], [0; zAv; 0], 'k', 'EdgeColor', ...
               [0.1 0.1 0.1], 'LineSmoothing','on', 'EdgeAlpha', 0.2, 'FaceColor', 'none'); % accum gray trace
 end % scan t steps (it)
 line([0 xAv], [0 yAv], [0 zAv], 'color', 'k', 'LineWidth',1.5); % av dipole direction
 text(xAv/2+5, yAv/2+2, zAv/2+2, 'Av dipole');
 htext =text(60, 70, 70, timestamp); % final time
 FinishAt =timestamp
  % 
 line([ 0 -100], [0 0], [0 0], 'color', 'b'); line([0 0], [0 50], [0 0], 'color', 'b');
 line([0 0], [0 0], [0 100], 'color', 'b'); % x y z  axis
 plot3(0, 0, 0,'bo', 'MarkerSize', 12); % mark A10, central ("apex") node
 plot3(0, 0, 0,'bo', 'MarkerSize', 6);
 
 
 %% Appx 2.0 Dispersion reln for 1 NM, excitn
  % excit popn (tau-e ~ 10ms etc) - cf wb6, p 125a etc
  
  C =100
   %omegae =sqrt(100 - 600*C +400*C^2); % rad/s  ;  ao ~ 1 x ae: ampl of osc
   %omegae =sqrt(100 - 6*C +4*C^2) % rad/s  ;  ao ~ 10^-2 x ae
  %omegae =sqrt(100 - 60*C +40*C^2) % rad/s  ;  ao ~ 10^-1 x ae : about right?
  omegae =sqrt(4e4 - 12*C +16*C^2) % rad/s  ;  ao/ae ~ 2/50 (mV) : from sims
   %argm= - 600*C +400*C^2   %argm= - 6*C +4*C^2 % argm= - 60*C +40*C^2
  argm= - 12*C +16*C^2;
  fe = omegae/(2*pi) % s^-1, Hz
  
  clear argm omegae fe 

 % Appx 2.1 Dispersion reln for 1 NM, inhib
  % inhib popn (tau-i ~ 20ms etc) - cf wb6, p 125x etc
  C =10
  
   omegae =sqrt(2.5e3 - 12*C +16*C^2) % rad/s  ;  ao/ai ~ 2/50 (mV) : from sims