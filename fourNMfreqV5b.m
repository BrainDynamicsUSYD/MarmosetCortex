% B {Pailtjorpe, U Sydney. fourNMfreqVx.m, (14/4/23)  
 % Jansen-Rit model tests, 4 node cluster: 3 alpha + 1 beta
 %  here, 1: A11, 2: A46D, 3: A9, 4: A10. dbl triangle layout 
    % based on twoNMfreqV5.m; wb6, p166, 177 (Feb '23)
%Dependencies: JR12v5.m - for RK4 DE solver, incl. pass r & wbar to fn J12; vary at ea NM
    % nb. V4.1 vary Gain[s] ~ syn strength, per NM; pass Cfac; add signal delays, via d(i,j)
         % V4c. use S[ wbar*V(t) -Vo)] -  uses JR12v5.m * : equivalent to V5 codes
         % V5. assign freq bands to marm anat areas; update (5/4/23), wb7, p57    
         % V5b. use S[~1/V, erfc (N(w)]; matches twoNM..v5b & sixCluster..v5b (15/4/23)
                % & JR12v5d.m : uses tests for erfc; has 0.2*excit stim to inhib popn
% Cells (sections)
% 1. ADj set up
% 2. parameters for WC/JR model
% 3. set up stimuli
% 4.0  RK4 DE solver, 1st order, for n variables
% 5.0  various PLOT OUTPUT
         
% 1.0 Set up
addpath('/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain/Models'); % for other model codes
fprintf('\n 4NM: 1+1+1+1; tune theta, alpha, beta, gamma bands; + delays; S[1/V - erfc], V5b code \n')
fprintf('\n      0 pc stim -> i popn \n')
 % Adj, for coupling nodes
nn=4; % number of neural masses (nodes)

%% 1.01 Adj: zeros (all)
 fprintf('  Adj(ij) = 0: no coupling  \n')
Adj= zeros(nn) % no coupling: indep NM 

%% 1.02 Adj: ones (all)
 fprintf('  Adj(ij) = 1: unit coupling, all-all  \n')
Adj= [0.0, 1.0, 1.0, 1.0; ...
       1.0, 0.0, 1.0, 1.0; ...
       1.0, 1.0, 0.0, 1.0; ...
       1.0, 1.0, 1.0, 0.0] % unit coupling, symm

%% 1.03 Adj: weak (all)
 fprintf('  Adj(ij) = 0.1: weak coupling, all-all  \n')
Adj=  [0.0, 0.1, 0.1, 0.1; ...
       0.1, 0.0, 0.1, 0.1; ...
       1.0, 1.0, 0.0, 0.1; ...
       0.1, 0.1, 0.1, 0.0] % unit coupling, symm
det(Adj)
%% 1.03.1  Adj: medium (all)
 fprintf('  Adj(ij) = 0.1: weak coupling, all-all  \n')
Adj= [0.0, 0.5, 0.5, 0.5; ...
       0.5, 0.0, 0.5, 0.5; ...
       0.5, 0.5, 0.0, 0.5; ...
       0.5, 0.5, 0.5, 0.0] % unit coupling, symm
det(Adj)
%% 1.03.2  Adj: weak (all)
 fprintf('  Adj(ij) = 0.1: weak coupling, #1 isolated  \n')
Adj= [0.0, 0.0, 0.0, 0.0; ...
      0.0, 0.0, 0.1, 0.0; ...
      0.0, 0.1, 0.0, 0.1; ...
      0.0, 0.0, 0.1, 0.0] % unit coupling, symm
det(Adj)

%% 1.03.3  Adj: weak (dbl triangle)
 fprintf('  Adj(ij) = 0.1: weak coupling, dbl triang;e  \n')
Adj= [0.0, 0.1, 0.0, 0.1; ...
      0.1, 0.0, 0.1, 0.1; ...
      0.0, 0.1, 0.0, 0.1; ...
      0.1, 0.1, 0.1, 0.0] % unit coupling, symm
det(Adj)
%% 1.03.4  Adj: ramp up wts (dbl triangle)
 fprintf('  Adj(ij) = up: weak coupling, dbl triang;e  \n')
Adj= [0.0, 0.1, 0.0, 0.1; ...
      0.1, 0.0, 0.1, 0.1; ...
      0.0, 1.0, 0.0, 0.1; ...
      3.0, 2.0, 2.0, 0.0] % unit coupling, symm; need (3-4) to 1 
det(Adj)
%% 1.03.4a  Adj: ramp up wts & tune feedbk (dbl triangle)
 fprintf('  Adj(ij) = varies: medium coupling, weak feedback, dbl triangle  \n')
Adj= [0.0, 0.1, 0.0, 0.4; ...
      0.2, 0.0, 0.1, 0.1; ...
      0.0, 1.0, 0.0, 0.5; ...
      3.0, 2.0, 2.0, 0.0] % medium coupling, weaker feedbk
det(Adj)

%% 1.03.4b  Adj: ramp up wts & tune feedbk (dbl triangle)
 fprintf('  Adj(ij) = varies: medium+ coupling, follow Marm, dbl triangle  \n')
Adj= [0.0, 0.1, 0.0, 0.4; ...
      0.2, 0.0, 0.1, 0.5; ...
      1.0, 1.0, 0.0, 0.5; ...  % nb. A(3,1) is new link added here
      5.0, 4.0, 1.0, 0.0] % medium coupling, weaker feedbk
det(Adj)

%% 1.04 Adj: medium (beta star)
 fprintf('  Adj(ij) = 0.5: beta star coupling  \n')
Adj= [0.0, 0.0, 0.0, 0.5; ...
      0.0, 0.0, 0.0, 0.5; ...
      0.0, 0.0, 0.0, 0.5; ...
      0.5, 0.5, 0.5, 0.0] % medium coupling, symm 
det(Adj)
%% 1.04 Adj: stronger (beta star)
 fprintf('  Adj(ij) = 1.0: beta star coupling  \n')
Adj= [0.0, 0.0, 0.0, 1.0; ...
      0.0, 0.0, 0.0, 1.0; ...
      0.0, 0.0, 0.0, 1.0; ...
      1.0, 1.0, 1.0, 0.0] % medium coupling, symm 
det(Adj)
%% 1.04.1 Adj: weak ( "square" closed loop)
 fprintf('  Adj(ij) = 0.1/0.1: "square" coupling  \n')
Adj= [0.0, 0.1, 0.0, 0.1; ...
      0.1, 0.0, 0.1, 0.0; ...
      0.0, 0.1, 0.0, 0.1; ...
      0.1, 0.0, 0.1, 0.0] % weak coupling, symm 
det(Adj)
%% 1.04.1a Adj: medium (Linear chain coupling; not closed loop)
 fprintf('  Adj(ij) = 0.1: linear coupling  \n')
Adj= [0.0, 0.1, 0.0, 0.0; ...
      0.1, 0.0, 0.1, 0.0; ...
      0.0, 0.1, 0.0, 0.1; ...
      0.0, 0.0, 0.1, 0.0] % weak coupling, symm 
det(Adj)

%% 1.04.2 Adj: medium ( chain: "square" closed loop)
 fprintf('  Adj(ij) = 0.5: closed loop, linear coupling  \n')
Adj= [0.0, 0.5, 0.0, 0.5; ...
      0.5, 0.0, 0.5, 0.0; ...
      0.0, 0.5, 0.0, 0.5; ...
      0.5, 0.0, 0.5, 0.0] % medium coupling, symm 
det(Adj)
%% 1.04.2b Adj: medium ( chain: "square" closed loop)
 fprintf('  Adj(ij) = 0.1/0.5/1.0: closed loop, asym coupling  \n')
Adj= [0.0, 0.1, 0.0, 0.5; ...
      0.1, 0.0, 0.5, 0.0; ...
      0.0, 1.0, 0.0, 0.5; ...
      1.0, 0.0, 0.5, 0.0] % medium coupling, symm 
det(Adj)
%% 1.04.3 Adj: medium (Linear / dbl alpha star)
 fprintf('  Adj(ij) = 0.1: linear coupling  \n')
Adj= [0.0, 0.1, 0.1, 0.0; ...
      0.1, 0.0, 0.0, 0.1; ...
      0.1, 0.0, 0.0, 0.1; ...
      0.0, 0.1, 0.1, 0.0] % medium coupling, symm 
det(Adj)
%% 1.04.3a Adj: medium (Linear ring / dbl alpha )
 fprintf('  Adj(ij) = 0.1: linear coupling  \n')
Adj= [0.0, 0.1, 0.1, 0.0; ...
      0.1, 0.0, 0.0, 0.1; ...
      0.1, 0.0, 0.0, 0.0; ...
      0.0, 0.1, 0.0, 0.0] % medium coupling, symm 
det(Adj)
%% 1.04.4 Adj: medium, asymm (beta star)
 fprintf('  Adj(ij) = 0.5/0.1: star coupling  \n')
Adj= [0.0, 0.5, 0.1, 0.5; ...
      0.1, 0.0, 0.5, 0.5; ...
      0.5, 0.1, 0.0, 0.5; ...
      0.1, 0.1, 0.1, 0.0] % medium coupling, asymm + edges
 det(Adj)
%% 1.04.4a  Adj: weak, symm (beta star)
 fprintf('  Adj(ij) = 0.1: star coupling to beta \n')
Adj= [0.0, 0.0, 0.1, 0.0; ...
      0.0, 0.0, 0.1, 0.0; ...
      1.0, 1.0, 0.0, 0.1; ... % stronger feedbk to alpha
      0.0, 0.0, 0.1, 0.0] % weak coupling, symm
 det(Adj) 
%% 1.04.4b  Adj: weak, symm (beta dbl star)
 fprintf('  Adj(ij) = 0.1: star coupling to dbl beta \n')
Adj= [0.0, 0.1, 0.1, 0.0; ...
      1.0, 0.0, 0.5, 0.1; ...
      0.5, 0.5, 0.0, 0.1; ... % x//stronger feedbk to alpha
      0.0, 0.1, 0.1, 0.0] % weak coupling, symm
 det(Adj) 
 %% 1.04.4c  Adj: weak, symm (alpha dbl star)
 fprintf('  Adj(ij) = 0.1: star coupling to dbl alpha \n')
Adj= [0.0, 0.1, 0.5, 0.1; ...
      0.1, 0.0, 0.1, 0.1; ...
      0.5, 0.5, 0.0, 0.0; ... % x//stronger feedbk to alpha
      0.1, 0.1, 0.0, 0.0] % weak coupling, symm
 det(Adj) 
%% 1.05 Adj: asymm links 
 fprintf('  Adj(ij): asymm coupling  \n')
Adj= [0.0, 1.0, 1.0, 1.0; ...
       1.0, 0.0, 1.0, 1.0; ...
       1.0, 1.0, 0.0, 1.0; ...
       1.0, 1.0, 1.0, 0.0] % unit coupling, symm
%% 1.06 Adj: ones (beta star)
 fprintf('  Adj(ij) = 0.5: ring coupling  \n')
Adj= [0.0, 0.5, 0.0, 0.5; ...
      0.5, 0.0, 0.5, 0.0; ...
      0.0, 0.5, 0.0, 0.5; ...
      0.5, 0.0, 0.0, 0.0] % medium coupling, symm 
det(Adj) 
%% 1.07 Adj:  (beta star cluster; asymm)
 % cf. Marm pFC(6)
 fprintf('  Adj(ij) = 1.0/0.1: beta star coupling; asymm  \n')
Adj= [0.0, 1.0, 1.0, 0.1; ... % weak feedback; strong feedfwd
      0.5, 0.0, 0.1, 0.1; ... % weak around edges
      0.1, 0.1, 0.0, 0.1; ... % all weak
      1.0, 1.0, 5.0, 0.0] % medium coupling, strong feedfwd from beta 
det(Adj)

%% 1.0.9  zero delays:  dij
 fprintf('  d(ij) =0: no signal delays  \n')
 Ad6A = zeros(nn);
%% 1.0.10  delays:  dij (mm / ms) from marm data 
 % cf # 1.0.4, of sixCluster codes:  d(ij) array -  set d in bins of 0.1 mm (approx)
fprintf('  d(ij) in 0.1mm bins: from Marm. \n')
 %Ad6A=Adist4; % cf appx 1.x below
Ad6A = [...
     0     2.8000  3.4000  2.5000; ...
    2.8000 0       1.8000  2.4000; ...
    3.4000 1.8000   0      2.4000; ...
    2.5000 2.4000   2.4000 0] 

%% 1.1 Adj from Marmoset (LNe 6x6) data
 % cf.Appx 2. in orig order: cf. logbook, p58
fprintf('  use marmoset LNe Adj(2): ordered here \n')
%clusterlist=[2, 26, 25, 48, 32, 3]; % as plotted in x-y plane, viewed from Ant.
clusterlist=[3, 32, 48, 2]   % 4 NM dbl triangle, as explored here; cf. logbook, p67
Anew= csvread('AdjMarmoset_Rescale2b.csv');  % based on LNe 
 Aclust = Anew(clusterlist, clusterlist);
Adj2=Aclust/1000  % scale used for 6 NM model
Adj=Adj2;   % for main code
nLinks= length(find(Adj2(:))) % 29 entries
%det(Adj)
clear nLinks Anew Aclust clusterlist Adj2 
 %Adj=zeros(nn); % test indep NM.
 %Adj=ones(nn); % test unit coupled NM. nb. diagonals (self links) present

% nb. set approx dij in Ad6A, above & at Appx 3.2.1

%% 2.0 Parameters: tune freq. band
fprintf('\n 1+1+1+1 NM;  dbl triangle \n')
% time axis & IC (0 mV)
x0=0; xf= 20.0; % t range (sec)
h = 0.0001;   % t step (0.1 ms; [& 0.5, 0.25 ms, for debug])
x = [x0:h:xf];  % t domain array
 %time =x/(1000);  %time =x/(1000*h); % time (sec) as row vec
tstart=10001; %tstart=2000; %1000; % for plots (@ 1, 0.5, 0.1 ms time steps)
%y = zeros(6,length(x)); % row vec for each variable
yic1=[1.0; 2.0; 1.0; 0.0; 1.0; 0.0]; % IC at t=0 for 1 node (6 DEs),  col vec, 6 rows
yic=5*yic1; % add additional nodes; and stronger?
for in=1:nn-1
 yic=[yic; yic1]; % add 6 rows per extra node; 12, 18 DEs for 2, 3 ... nodes
end
 clear yic1
%
% Model parameters (orig JR) : pass to fn subroutine! cf wb7, p66
v0=6;    % (mV) midpoint (switching trheshold V) of sigmoid function, for V -> pulse conversion 
vm=5;    % (s^-1) max firing rate for sigmoid fn
  %rvec= [0.5 0.43 0.5 0.56] % v1. theta, alpha, beta, gamma; mixed
rvec= [0.5 0.5 0.5 0.5] % v2. beta, theta, alpha, gamma; revised (wb7, p57)
   %rvec= [0.5 0.5 0.5 0.5] % default
   %rvec= [0.45 0.32 0.38 0.35 0.30 0.38] % for pFC; revised: alpha band: wb#4, p42 (18/9/21)
% feedback: tuned to bands
   %Cvec = [135 135 300 500]  % orig for alpha / beta /  gamma
  % Cvec = [180 200 250 350]  % test orig
  % Cvec = [180 200 250 350]  % orig  beta/ theta /alpha / gamma [more marginal osc]
Cvec = [250 200 220 220]  % retuned  beta/ theta /alpha / gamma 
R12=100;   % NM - NM coupling coeff (ie. link wt) (cf. Goodfellow('11))
           % base coupling: 100 
 % tune feedback gain (based on WtIn, kIn)
CfacVec= [1.0 1.0 1.0 1.0]; % pass to fn J12a.m integrator
    % nb. Assign in J12: C1=C1*Cfac;  C3=C3*Cfac; 
% weight/synp link, for Sigmoid response fn:
  %wvec= [1.0 1.0 1.0 1.0] % default case: basic S[V] sigmoid
  %wvec= [1.5 1.5 3.0 5.0] % larger w for  larger r
wvec= [3.0 1.5 1.5 3.0] % revised: for S[v ~w] & S[1/V,erfc]
   %wvec=[0.52 3.23 0.48 0.74 1.30 2.35] % for pFC; raw wt/link: av in each NM: wb4,p59a
    %wvec=[1.08 6.73 1.00 1.54 2.71 4.90] % ditto: rescale>=1
% signal delay: set by array Ad6a, above; cf. App3.2, based on d(i,j) % d12=1.5-4;  % eg. mm / ms
 %Ad6A=zeros(nn); % debug; no delays (base case)
% Rate params - pass to JR12a fn.
 A= 3.25;  B =22;  % Max PSP amplitude for pyramidal (A), inhib (B) & excit (kA * A) interneurons
   %a=100; b=50;  % (s^-1) / default JR; set to osc at 8Hz (taus: 20, 20 mS)
    %taue= 15; taui= 20; % set Time Const (excit, inhib) (ms)  :  alpha, beta bands
    % retune in DE loop:  keep a*A & ? b*B const
       %a=1000/taue; b=1000/taui; % rate const (s^-1)
     %tauevec=[10.0 20.0 15.0 5.0] % revised (early/4/23) for S[V ~ w]
     %tauivec=[10.0 22.0 16.0 5.0] % for NM 1-4
tauevec=[10.0 25.0 16.0 5.0] % revised again (15/4/23): wb7, p68
tauivec=[10.0 25.0 17.0 5.0] % for NM 1-4
 %[1000*A/tauevec(1)  1000*B/tauivec(1)] % check "driving forces"

% Coupling const for S[v]: cf. Steigler(2011)
 kIn2=1; kIn6=1; 
 
% >>>>>>>>>>>>>>
%% 3.0 Zero stim -pulse, ran or wave :: Base case
fprintf('\n >> zero inputs: 1 mV G ran noise \n')
Vyt =zeros(1,length(x));
pran=zeros(1,length(x)); 
 %Ip=zeros(1,length(x));  %Ip6=zeros(1,length(x)); % not used here
IpOn=zeros(nn,1);  % no stim inputs; need IpOn also
Pd=zeros(1,length(x)); % Pd2=zeros(1,length(x)); %  pulse density fn - cf #3.6
 %pran = 0.1*rand(1,length(x)); % low noise;  Uniform ran noise:
pran = 1*randn(1,length(x)); % G ran noise (+-);
  %  tba ??
InI =0  % controls stim to Excit/Inhib popn:  default - apply stim to Excit popn
%InI =1 % apply stim to Inhib
% 
 %% 3.5 Constant stim rate - use Pd; not filtered by S[v]
% explore stable pts and limit cycles.
fprintf('\n ++const stimulus: at 4.0 sec ')
tstart=40001; % for plots
StimLevel= 150
Pd=zeros(1,length(x)); tstart0= tstart; %tstart0= 1001; %tstart0= 20001; %  separate for Pd
  %tend=100001; % for short stim [4:5sec]
  % Pd(tstart0:tend)=StimLevel; % short, const stim
Pd(tstart0:end)=StimLevel; % && turn on the const stim!!


%Ip6=zeros(1,length(x)); % pulse train/fn for 2nd input (?dphi shifted) - node 6

clear StimLevel PulseHeight tstart0

%% 3.6 Load 1st Pulse density function[s]
 % a la JR('93 & '95); set up in code WaveTrain.m @ #1.4
 fprintf(' Pulse shape stim(@ 4s) + 0 (@0.5s); & 1 mV ran noise: \n')
 Pd= 100.0*Pulse1; % 1 or 2 pulses;  vary amplitude (ie. pulse density
 tstart=40001;  % cf. wave... code, for pulse setup
  % Pd= 1.0*PulseN; % N pulses, in train
 %Pd = Pd+100.0;  % add const stim (all t)
 Pd((1001):end) = Pd((1001):end)+ 0.0;  % add const stim , with pulse (not before)
  PulseStats= [mean(Pd) std(Pd) max(Pd)]
  PulseStatsAtStim = [mean(Pd(40001:50000)) std(Pd(40001:50000)) ]
  PulseIntegrl = sum(Pd)*h; % sum across rows % trapezoidal rule 
  PulseIntegrl=PulseIntegrl'

 %Pulse2 =zeros(1,length(x)); % default: no 2nd pulse
 InI =0  % controls stim to Excit/Inhib popn:  apply stim to Excit popn (default)
  %InI =1 % apply stim to Inhib 
  
  clear PulseHeight PulseIntegrl PulseStats* 

 %% 3.6.1 Load 2nd Pulse density function[s]
 fprintf(' add 2nd Pulse shape stim: \n')
 %Pd2= 1.0*Pulse2; % 1 or 2 pulses;  vary amplitude (ie. pulse density
  %Pd2 = Pd2+155.0;   % add const stim
 Pd = Pd + 1.0*Pulse2; % add 2nd pulse
 IpOnB =zeros(nn,1); % switch 2nd stim on/off to ea node - tbd
 
 %% 3.6.2 Load Spike train fn - eg. from from LIF sim
  % Load eSR (10ms bins) into eSpikes array in WaveTrains #1.6; from COBN.m .c
 fprintf(' add in LIF spike train stim x frn: \n')
 efrStim=csvread('efrStim6_300.csv' ); % read saved file, from LIF sim (M=300, noise 6)
  %efrStim=csvread('efrStim6_300_40s.csv' ); % 40s sim
 efrStim = [efrStim', 0]; % match length to local variables
   %efrStim = [efrStim', efrStim', 0]; % for 40s
 
 Pd= Pd + 10.0*efrStim; % shape matches lfp, x, y etc  [col vec]
   %Pd =[Pd, 0]; % ? adjust length to match (t =0 "problem")
  % Pd= 1.0*eSpikes'; % alt form
  % Pd= 0.042*eSpikes'; % also: frn stim ~ wt frn
   %Pd = Pd+155.0;   % add const bias stim - optional
 
  % another Alt. use use COBN output eFR (in 0.1ms bins)
     %Pd=zeros(1,length(x)); % correct length
     %Pd(2:end) = eFR; % matches lengths
 PulseHeight= max(Pd)
 clear PulseHeight
 
 %% 3.6.2a Load current LIF output
  eFR =[double(eFR); 0.0]; % match length
  Pd= Pd + 10.0*eFR'; % resampled LIF output
  
%% 3.7 Input channels
 % selected nodes to get the stimulus inputs: Ip, Pd 
 % IpOn=ones(nn,1); % 1st stimulus On for ALL nodes
  IpOn=zeros(nn,1); % switch 1st stim on/off to ea node
 % stimulate Node 2 (A-32V): strongest, more wt, fewer links)
 IpOn(1)=1; % 1st beta; A11.  cf wb7, p68
% IpOn(2)=1;  % [2] theta {A46D
% IpOn(3)=1;  %  [3] alpha { A9 
% IpOn(4)=1;  % gamma  {A10 assigned }
  StimTo =IpOn' % check
  % direct 2nd pulse, if present?
    %IpOnB=zeros(nn,1); % switch 1st stim on/off to ea node
 clear StimTo
 
 %% 3.8 Gather stimuli & Plot, to check stim
 figure;  plot(pran, 'Color', [0.9 0.9 0.9]);  hold on
 plot(Pd, 'b'); title('pulse + Vstimulus + waves + noise');
  plot(Vyt, 'g--'); xlabel('t steps (0.1 msec) ')
   if max(Pd+pran+ Vyt) < 5
     ylim([-2 2.2])  %ylim([0 1])  
   else
     ylim([-2 (max(max(Pd +pran)) +20)])
   end
  % xlim([3.9e4 5e4]) % focus on time of stim
  confirm_Excit0Inhib1 = InI
   clear confirm_*
   
% >>>>>>>>>>>>>> 
%% 4.0  RK4 DE solver, 1st order, for n variables
% nb. matlab has default double precision: needed here
tic
% check stimulus
   %figure;  plot(Ip+pran); title('stimulus + waves + noise'); hold on; plot(Vyt, 'r-'); xlabel('t (msec) ')
% reset arrays
y=zeros(nn*6,length(x)); % set up array for solutions; 2,3, etc nodes (6 rows each)
 %y=double(y); % enforce double precision (64 bit) : should be default
deltay=zeros(nn,length(x)); % set up array for y3-y2; 6 DEs (1 row ea)
 %deltay=double(deltay);
y(:,1)= yic ; % y's are row vec, one for ea DE; IC is a col vec, for i=1 or t=0. 
  for i=1:50  % may need to escape "trap" at zero
   y(:,i)= yic ; % extend IC - for delay DE "late start" - eg to ~1-5 ms
  end
In12=0;  % reset here; and for ea  NM
% RK4 DE solution step:  iterate forward in time (here x); scan nodes(in) & linked NN (j)
  % fn JR12 incl forcing fn. as arguments, at this t-step; 
  % inputs: pran + p(i,j)+ [NM-NM couplings] {= sum(y1-y2) ; 
       %fprintf('\n   i  j wt(ij) dr(ij) index r Cfac ') % debug- header
        % start "late" to allow for delay back ref, to L [allow 5ms w 0.1 steps]
for i = 50:(length(x)-1) % scan time steps & start at 1, 10, 100 [for delay DE]
    for in =1:nn % scan nodes in cluster (NM)
         j11=(in-1)*6+1; j16=(in-1)*6+6;  % working on this node #in
         r=rvec(in); C=Cvec(in); Cfac=CfacVec(in); wbar=wvec(in); % for indiv NM : slope of S[v]; Gain C1,3
         a=1000.0/tauevec(in); b= 1000.0/tauivec(in); % time const for ea. NM (nb. ms)
          %A=3.25; A=3.25*100/a; % ? rescale, to keep a*A const. Optional?
          %[(a*A) (b*B)]  % debug
          %fprintf(' ... working on node  %4.0f \n', in)
    % Node #in - gets the ext inputs (Ip, Ii, p and In12 from linked nodes (i <- j)
     % and Vwave adds to In12 (a la lfp, at soma)
       jlist=find(Adj(:,in));  % linked nodes in (col) to node #in
       jnn=length(jlist); % #NN of in
       for j=1:jnn % scan NN list of node #in
           In12=0; %  initialise for ea NM
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
           dyi= y(j11+1,i-deltai)-y(j11+2,i-deltai); % delta-y of nn. JR case (input to pyrm)
            %dyi= y(j11,i-deltai); %  ?test: use y0 of nn. pyrm for i-j interactions
           argm=r*(wbar -v0/dyi); In12=0; 
          if argm >0 & dyi>0.01  % avoid difficulties of erfc(neg, 0)
            In12= R12*Adj(jj,in)*vm*(1.0-erfc(argm)/2); % erfc form for wt; avoid v=0
          end 
            %[i dyi argm In12]  % debug  
            % In12= R12*Adj(jj,in)*vm*(1.0-erfc(r*(wbar -v0/dyi))/2);  % erfc form for wt; avoid v=0?
          % workin on node in:  weight by In link, from j to i: jj -> in 
    %fprintf('\n  %2.0f %2.0f %5.3f %4.1f %4.0f %4.2f %4.2f ', in, jj, Adj(jj,in), Ad6A(jj,in), (i-deltai), r, Cfac )
           end % check for NN  
       end % scan linked NNs; gather inputs (in <- jj) from linked NNs
         %In12=In12+Vyt(i); % add travelling/standing wave direct, here, now
           %fprintf('\n  ') % for debug
           %In21=zeros(12,length(x)); % debug, no input from other node[s]
           %[jj j11 j16]  % debug:
        j11=(in-1)*6+1; j16=(in-1)*6+6;  % update node #in
          % nb. y7:12 in Node-2; y13:18 is Node-3, etc.
        % Stimulus pulses to selected node (eg 2, 6), via IpOn switch
        % add travelling wave Voltage, at y=16mm (Ant), at this time(i), to _all_ NM
          %Istim=Vyt(i)+IpOn(in)*Ip(i); %pstim=IpOn(in)*pran(i); % turn stim on only for selected nodes - set above, at IC.
         % gather ran & wave train together: to all NM
         %Istim=IpOn(in)*Ip(i); % use if spike train Ip present
         pstim=pran(i); %  & ran noise to ALL nodes
         Pdstim=IpOn(in)*Pd(i);   % switch 1st stim pulse on; 
            % nb/ here InI replaces Ip, and is 0; feed-fwd,bk calc in JR fn.
        k_1 = JR12v5d( x(i), y(j11:j16,i), Pdstim, InI, pstim, In12, a, b, A, B, C, Cfac, r, wbar);  % these k's should be col vec, one for ea DE
        k_2 = JR12v5d( x(i)+0.5*h, y(j11:j16,i)+0.5*h*k_1, Pdstim, InI, pstim, In12, a, b, A, B, C, Cfac, r, wbar); % seem to need to force col vec ?
        k_3 = JR12v5d( (x(i) +0.5*h), (y(j11:j16, i) +0.5*h*k_2), Pdstim, InI, pstim, In12, a, b, A, B, C, Cfac, r, wbar);
        k_4 = JR12v5d( (x(i)+h), (y(j11:j16,i) +k_3*h), Pdstim, InI, pstim, In12, a, b, A, B, C, Cfac, r, wbar);
        y(j11:j16,i+1) = y(j11:j16,i) + (k_1 +2*k_2 +2*k_3 +k_4)*(h/6); % load y1:y6 at this t-step      
    end %  scan nodes
end  % scan t-steps
fprintf('\n  calc done/ ') % for debug
toc %timing
 clear argm i in deltai In12 j* k_ a b Istim pstim pdstim vlocal 
 
%% 5.0 PLOT OUTPUT 
% Node 1 output [this is centre of star cluster, and Out-Hub (A-10)
  tstart=40001; % 4 sec
  %figure; plot(y(1:3,:)'); title('JR/ node-1: y0, y1, y2')  % nb. need cols
  %title('A-10:  JR/ node-1, outputs:'); 
  %legend('y0', 'y1', 'y2'); % orig notation [not array index]
figure; plot(x, y(1,:));  title('4 NM, node-1, pyrm: output: y-e-pyrm '); xlabel('t (s) ')
%figure; plot(y(2,:)'-y(3,:)'); 
% title('A-10:  JR/ node-1, output: detla-y: y1(1)- y2(1) '); 
 figure; subplot(3,1,1); plot(x, y(2,:)'); hold on; plot(x, y(3,:)'); 
 title('4NM, node-1: output: y-e, y-i '); 
 legend('y-e(1)', 'y-i(1)', 'Location', 'SouthEast'); hold off
 subplot(3,1,2); plot(x, y(2,:)'); hold on; plot(x, y(3,:)');   
 xlim([3.9 5.4]); 
 subplot(3,1,3); plot(x, y(2,:)'); hold on; plot(x, y(3,:)'); 
 grid on; xlabel('t (s) '); xlim([13.6 15.6]);
 figure; plot(x, (y(2,:)-y(3,:)) ); title('4-NM, node-1, delta-y '); 
 hold on; xlabel('t (sec) ')  % checktransition in dy1
  dy1Sample=y(2,45001:65000)-y(3,45001:65000); % 2 sec sample, 0.5s after stim
  dy1Up = mean(dy1Sample) + std(dy1Sample); 
  dy1Down = mean(dy1Sample) - std(dy1Sample);
  plot([4.5, 7], [dy1Up, dy1Up], '--r', 'LineWidth', 1.5)
  plot([4.5, 7], [dy1Down, dy1Down], '--r', 'LineWidth', 1.5); hold off
  dy1RangeSD = [ mean(dy1Sample) (dy1Up - dy1Down) ]
  dy1RangeAbs = max(dy1Sample) - min(dy1Sample)
   clear dy1Sample dy1Up dy1Range*
   
% Node 2 outputs [this has the primary inputs]:
 %figure; plot(y(8,(tstart-100):end) - y(9,(tstart-100):end));  title('JR node-2: output, y2(2) - y1(2)')
  %figure; plot(y(7:12,:)'); title('JR/ node-2: output, y1, y2, .. y6')  % nb. need cols
  %legend('y0', 'y1', 'y2', 'y3', 'y4', 'y5'); % 
 figure; subplot(2,1,1); plot(x,y(8,:)'); hold on; plot(x,y(9,:)'); 
 title('JR/ node-2: output: y1, y2 '); legend('y1(2)', 'y2(2)'); hold off
 subplot(2,1,2); plot(x,y(8,:)'); hold on; plot(x,y(9,:)'); 
  %axis([1.5 2.5 20 100]);
  xlim([3.9 5.4]); grid on; xlabel('t (s) ')
     % other outputs
     %figure; plot(y(6,(tstart-100):end));  title('JR model, Node-1 output, y5  ')
     %figure; plot(y(12,(tstart-100):end));  title('JR model, Node-2 output, y5  ')
     % axis([0 5 -500 0]);  axis([0 20  -10000 10000]);
     figure;  plot(y(8,(tstart-800):end) - y(9,(tstart-800):end)); title('Node 2, delta-y ')
     xlabel('t (x0.1 ms) ')
% Node 3 outputs:  
 figure; subplot(2,1,1); plot(x, y(14,:)'); hold on; plot(x, y(15,:)'); 
 title('JR/ node-3 output: y1, y2 '); legend('y1(3)', 'y2(3)'); hold off
 subplot(2,1,2); plot(x, y(14,:)'); hold on; plot(x, y(15,:)');  
 xlim([3.9 5.4]); grid on; xlabel('t (s) ')
 figure; plot(y(14,(tstart-800):end) - y(15,(tstart-800):end));
 title('JR model: node-3 output, dy: y2(3) - y1(3)') 
% Node 4 outputs:
 figure; plot(y(20,(tstart-800):end) - y(21,(tstart-800):end));
  title('JR model: node-4 outputs, dy: y1(4) - y2(4)')
 figure; subplot(2,1,1); plot(x, y(20,:)'); hold on; plot(x, y(21,:)'); 
 title('JR/ node-4: output: y1, y2 '); legend('y1(4)', 'y2(4)'); hold off
 subplot(2,1,2); plot(x, y(20,:)'); hold on; plot(x, y(21,:)');  
 xlim([3.9 5.4]); grid on; xlabel('t (sec) ')
 
 figure; plot(x, y(19,:));  title('4 NM, node-4, pyrm: output: y-e-pyrm '); xlabel('t (s) ')
  ylim([0 0.2]);  xlim([4 4.2]) 
% all 4 nodes together:
tstart= 40001; % debug
 figure; plot(x((tstart-800):end), y(2,(tstart-800):end) - y(3,(tstart-800):end), 'k' );  hold on % node #1, A-10 
 xlim([4.0 5.0]); xlabel('t (s)'); xlabel('t (ms)');  %pause
 plot(x((tstart-800):end) ,(y(8,(tstart-800):end) - y(9,(tstart-800):end) ));   % #2,
 %pause  % - to closely examine waveforms
 plot(x((tstart-800):end) ,y(14,(tstart-800):end) - y(15,(tstart-800):end)); % #3
 %pause
 plot(x((tstart-800):end) ,y(20,(tstart-800):end) - y(21,(tstart-800):end), 'b'); % #4,
  %plot(x((tstart-800):end) ,Ip((tstart-800):end), 'g');  
  %plot(x((tstart-800):end) ,Ip6((tstart-800):end), 'g--') % stimulii, on same scale
 plot(x((tstart-800):end) ,Vyt((tstart-800):end), 'r--') % stimulii: wave
  legend('delta-y(1)', 'delta-y(2)', 'delta-y(3)', 'delta-y(4)', ...
       'stimulus-2' , 'wave' ); %text(3000, 0.3, 'wt(1-2) = 100')
  % axis([0 5000 -20 20]);  %axis([1 2 -25 25]); 
%
% mean and range of detla-y, overall: 4 nodes together:  cf. #4.3 below, for calc at 17:19sec 
 d_yMean=zeros(nn,1); d_yAmpl=zeros(nn,1); 
 d_yMean(1)=mean(y(2,(tstart+2000):end) - y(3,(tstart+2000):end));   % node #1,  
 d_yAmpl(1)=max(y(2,(tstart+2000):end) - y(3,(tstart+2000):end)) ... % sample after stim settles
     - min(y(2,(tstart+2000):end) - y(3,(tstart+2000):end));
 d_yMean(2)=mean(y(8,(tstart+2000):end) - y(9,(tstart+2000):end));   % #2,
 d_yAmpl(2)=max(y(8,(tstart+2000):end) - y(9,(tstart+2000):end)) ...
     - min (y(8,(tstart+2000):end) - y(9,(tstart+2000):end));
 d_yMean(3)=mean(y(14,(tstart+2000):end) - y(15,(tstart+2000):end)); % #3 
 d_yAmpl(3)=max(y(14,(tstart+2000):end) - y(15,(tstart+2000):end)) ...
      -min(y(14,(tstart+2000):end) - y(15,(tstart+2000):end));
 d_yMean(4)=mean(y(20,(tstart+2000):end) - y(21,(tstart+2000):end)); % #4, 
 d_yAmpl(4)=max(y(20,(tstart+2000):end) - y(21,(tstart+2000):end)) ...
      - min(y(20,(tstart+2000):end) - y(21,(tstart+2000):end));
 disp(' dy(i): mean, amplitude (global) \n'); [d_yMean d_yAmpl] % test/ debug
% mean and range of detla-y, later, after stim
tstart=90001; % [9: end] sec
 d_yMean=zeros(nn,1); d_yAmpl=zeros(nn,1); 
 d_yMean(1)=mean(y(2,(tstart+2000):end) - y(3,(tstart+2000):end));   % node #1,  
 d_yAmpl(1)=max(y(2,(tstart+2000):end) - y(3,(tstart+2000):end)) ... % sample after stim settles
     - min(y(2,(tstart+2000):end) - y(3,(tstart+2000):end));
 d_yMean(2)=mean(y(8,(tstart+2000):end) - y(9,(tstart+2000):end));   % #2,
 d_yAmpl(2)=max(y(8,(tstart+2000):end) - y(9,(tstart+2000):end)) ...
     - min (y(8,(tstart+2000):end) - y(9,(tstart+2000):end));
 d_yMean(3)=mean(y(14,(tstart+2000):end) - y(15,(tstart+2000):end)); % #3 
 d_yAmpl(3)=max(y(14,(tstart+2000):end) - y(15,(tstart+2000):end)) ...
      -min(y(14,(tstart+2000):end) - y(15,(tstart+2000):end));
 d_yMean(4)=mean(y(20,(tstart+2000):end) - y(21,(tstart+2000):end)); % #4, 
 d_yAmpl(4)=max(y(20,(tstart+2000):end) - y(21,(tstart+2000):end)) ...
      - min(y(20,(tstart+2000):end) - y(21,(tstart+2000):end)); 
 disp(' dy(i): mean, amplitude (> 9sec) \n'); [d_yMean d_yAmpl] % test/ debug
  
  figure; stem(d_yMean, '+');  hold on; stem(d_yAmpl, ':diamond');
  % hs(1).Marker='+'; hs(2).Marker='diamond'; 
 xlim([0 7]);  %axis([0 7 -15 60 ]); 
 xlabel('NM #'); ylabel('delta-y12(i)  (mV)'); title('dy(i) overall')
  legend('mean', 'p-p ampl.', 'Location', 'northeast') % place out of the way
   % text(5.5, 23, 'Adj = 1');  % text(5.5, 10, 'alpha; C x vector')  
  %   %figure; plot(rvec, d_yAmpl, '.', 'MarkerSize', 20); hold on
   % plot(rvec, d_yMean, '+', 'MarkerSize', 5);
 
% 3 key nodes together:
 figure;   plot(y(2, (tstart-800):end)-y(3,(tstart-800):end) ); hold on; % #1 
 plot(y(8,(tstart-800):end) - y(9,(tstart-800):end), 'r') % ', 'Color', [0 0.3 0]);  % #2,
 %plot(y(32,(tstart-800):end) - y(33,(tstart-800):end),'g'); % #6, A-11 : Input 2nd
  plot(y(20,(tstart-800):end) - y(21,(tstart-800):end),'b'); % #4,
  plot(Pd((tstart-800):end), 'k');  %plot(Ip((tstart-800):end), 'k--') % stimulii, on same scale
  title('4 NM cluster, JR model:  delta-y ')
  legend('delta-y(1)', 'delta-y(2)', 'delta-y(4)', 'stimulus'  );  

%  
% LFP  a) net lfp output: sum: detla-y(1) + detla-y(2), for the 4 NM
 lfp= (y(2,:) - y(3,:) + y(8,:) - y(9,:) +y(14,:) - y(15,:) ...
     + y(20,:) - y(21,:) )/4; % av of 4 nodes
 %    b) alt. av of y0's [pyrm]
  lfp0= (y(1,:) + y(7,:) +y(13,:)  + y(19,:) )/4; % av of 4 nodes
  
% figure; plot(lfp); title('JR model: net lfp output ')
 figure; plot(x, lfp); title('4-NM cluster, JR model: net lfp output '); xlabel('t (sec) ')
  ylabel('LFT (t)  (mV)'); % hold on; plot(20,-1, 'b^'); plot(20,-1.1, 'b^'); % wave starts
 figure; plot(x, lfp0); title('4-NM cluster, JR model: net lfp[y0] output '); xlabel('t (sec) ')
  ylabel('LFT (t)  (mV)'); ylim([0 0.2]); 
%   
% a) Group of 4 close together: net lfp output: sum: delta-y(2:5) 
% lfp25= (y(2,:) - y(3,:) + y(8,:) - y(9,:) +y(14,:) - y(15,:) ...
 %    + y(20,:) - y(21,:) )/4; % av of 4 nodes
 % lfp_stats_2to5= [mean(lfp25(tstart:end))  (max(lfp25(tstart:end)) - min(lfp25(tstart:end)) )]
 % clear lfp25  lfp_stats_2to5
  
% Grouped summary: resize figure [5 by 1 plots]
  %figwidth = 1024; figheight = 896; figposition = [100, 100, figwidth, figheight]; % large
  figwidth = 560; figheight = 704; figposition = [500, 100, figwidth, figheight]; % tall
 figure('position',figposition, 'units','pixels');  %figure; % default
 subplot(5,1,1);
 plot(y(2,:) - y(3,:),'b'); % #1, A-10 : Output
 tlim=length(y(2,:) - y(3,:))-1; axis([0,tlim, -Inf, Inf]) % keep vetical as is
  title('4-NM cluster, NM#1, delta-y ')
 subplot(5,1,2);
 plot(y(8,:) - y(9,:), 'r') % ', 'Color', [0 0.3 0]);  % #2,
  title('     NM#2, delta-y '); axis([0,tlim, -Inf, Inf])
 subplot(5,1,3);
 plot(y(20,:) - y(21,:) ); %  #4,
  title('     NM#4, delta-y '); axis([0,tlim, -Inf, Inf])
  subplot(5,1,4);  plot(lfp,'k'); title('     lfp output ');  %  net output
  xlim([0,tlim]); % xlim([0 3.0e5]); %axis([0,tlim, -Inf, Inf])
 subplot(5,1,5); 
  plot((Pd)*IpOn(1) +pran+Vyt, 'g--'); hold on; % theta
  plot((Pd)*IpOn(4) +pran+Vyt, 'r--') % stimulii to gamma, on same scale
  plot(pran, 'Color', [0.9 0.9 0.9]); 
  xlim([0,tlim]); ylim([-10 (max(Pd)+10)])
  % plot((Ip+Pd)*IpOn(1) +pran+Vyt, 'b--') % if needed?
  legend( 'stimulus-1', 'stimulus-4', 'ran'  ); title('     Inputs: stimulus-> NM ');  
   xlabel('t steps (0.1 msec) '); %axis([0,tlim, -Inf, Inf])
    % text(5500, 50, 'Stim 300 to both #2, #6')
    
%   & LFP(dy12)
 lfp_sample=lfp(40001:end); % during stimulus [0.5s sample]
   %lfp_sample=lfp(1100:2900); % debug, nb. wider sample for half t-step
 lfp_stats= [ mean(lfp_sample), (max(lfp_sample)-min(lfp_sample))]
 % figure; plot(lfp_sample - mean(lfp_sample)); title('JR, 6xNM: sum (lfp -mean), during stim. ')
 
% pyrm: y0 outpute for ea NM
 figure; plot(x, y(1,:)); hold on;
 plot(x, y(7,:));   plot(x, y(13,:));   plot(x, y(19,:));
 title('4 NM cluster, JR model:  y0(i) ')
  legend('y0(1)', 'y0(2)', 'y0(3)', 'y0(4)'  ); xlabel('t (sec) ')
  ylim([-0.1 0.5]); xlim([3 5]); % xlim([0 0.5]) % nb transients at start up
%  
 lfp0_sample=lfp0(2000:40000); % during stimulus [0.5s sample]; avoid transient
 lfp0_stats_Before= [ mean(lfp0_sample), (max(lfp0_sample)-min(lfp0_sample))]
 lfp0_sample=lfp0(40001:end); % during stimulus [0.5s sample]
 lfp0_stats_After= [ mean(lfp0_sample), (max(lfp0_sample)-min(lfp0_sample))]
 
% mean and range of y0, overall: 4 nodes together, after stim (4s) 
 d_y0Mean=zeros(nn,1); d_y0Ampl=zeros(nn,1); 
 d_y0Mean(1)=mean(y(1,(tstart+2000):end) );   % node #1,  
 d_y0Ampl(1)=max(y(1,(tstart+2000):end) ) - min(y(1,(tstart+2000):end) ); % sample after stim settles
 d_y0Mean(2)=mean(y(7,(tstart+2000):end) );   % #2,
 d_y0Ampl(2)=max(y(7,(tstart+2000):end) ) - min (y(7,(tstart+2000):end) );
 d_y0Mean(3)=mean(y(13,(tstart+2000):end) ); % #3 
 d_y0Ampl(3)=max( y(13,(tstart+2000):end) ) -min(y(13,(tstart+2000):end) );
 d_y0Mean(4)=mean(y(19,(tstart+2000):end) ); % #4, 
 d_y0Ampl(4)=max(y(19,(tstart+2000):end) ) - min(y(19,(tstart+2000):end) );
 disp(' y0(i): mean, amplitude (after stim) \n'); [d_y0Mean d_y0Ampl] % 
 
% 
 % if xf >15 % for longer sim (eg 20 sec)
  % lfp_sample=lfp(180000:190000);   % V. much later
  % lfp_stats_at18secLaterStill= [ mean(lfp_sample), (max(lfp_sample)-min(lfp_sample))]
  % end
    
 clear jlist i j in jj j11 j16 k*  In12 deltay lfp_* lfp0_* d_* dy6* % tstart  Ip pran yic
% > > > > >

%% 5.0.1 compare dy12 vs y0
fprintf('\n > plot yo vs y-e, y-i: \n');
figure; plot(x, (y(1,:)*50.0 +40), 'b' ); % its small - check offsets
hold on; plot(x, y(2,:), 'g', 'LineWidth', 1.2);  
plot(x, y(3,:), 'r', 'LineWidth', 1.0 ); % for y-e & y-i
 % hold on; plot(x, (y(2,:) -y(3,:)) ) % for lfp
legend('yo *50', 'y1', 'y2')
 % legend('yo *50', 'y1 -y2')
xlim([0.2 0.4]); grid on
 % ylim([-25.0 30])
 title('4NM, dbl tri; #1, beta; yo lag varies')
%
figure; plot(x, (y(7,:)*50.0 +40), 'b' ); % its small - check offsets
hold on; plot(x, y(8,:), 'g', 'LineWidth', 1.2);  
plot(x, y(9,:), 'r', 'LineWidth', 1.0 ); % for y-e & y-i
 % hold on; plot(x, (y(2,:) -y(3,:)) ) % for lfp
legend('yo *50', 'y1', 'y2')
 % legend('yo *50', 'y1 -y2')
 %xlim([0.2 0.4]); 
xlim([0.2 1.2]); grid on; 
 % ylim([-25.0 30])
 title('4NM, dbl tri; #2, theta; yo & ye, yi')  
 %
figure; plot(x, (y(13,:)*50.0 +40), 'b' ); % its small - check offsets
hold on; plot(x, y(14,:), 'g', 'LineWidth', 1.2);  
plot(x, y(15,:), 'r', 'LineWidth', 1.0 ); % for y-e & y-i
 % hold on; plot(x, (y(2,:) -y(3,:)) ) % for lfp
legend('yo *50', 'y1', 'y2')
 % legend('yo *50', 'y1 -y2')
 %xlim([0.2 0.4]); 
xlim([0.2 1.2]); grid on; 
 % ylim([-25.0 30])
 title('4NM, dbl tri; #3, alpha; yo & ye, yi') 
%
figure; plot(x, (y(19,:)*50.0 +20), 'b' ); % its small - check offsets
hold on; plot(x, y(20,:), 'g', 'LineWidth', 1.3);  
plot(x, y(21,:), 'r', 'LineWidth', 0.9 ); % for y-e & y-i
 % hold on; plot(x, (y(2,:) -y(3,:)) ) % for lfp
legend('yo *50', 'y1', 'y2')
 % legend('yo *50', 'y1 -y2')
 %xlim([0.2 0.4]); 
xlim([0.2 1.2]); grid on; 
 % ylim([-25.0 30])
 title('4NM, dbl tri; #4, gamma; yo & ye, yi')
 
%% 5.0.2 compare lfp(dy12) vs y0
figure; plot(x, (y(1,:)*50.0 +40), 'b' ); % its small - check offsets
 %hold on; plot(x, y(2,:), 'g', 'LineWidth', 1.2);  plot(x, y(3,:), 'r' ); % for y-e & y-i
 hold on; plot(x, (y(2,:) -y(3,:)), 'k' ) % for lfp
 %legend('yo *50', 'y1', 'y2')
 legend('yo *50', 'y1 -y2')
xlim([0.2 0.4]); grid on
 % xlim([0.2 1.2]);
 % ylim([-25.0 30])
 title('4NM, dbl tri; #1, beta; yo lag varies')
 
%% 5.1 revised LFP:  net lfp output: sum: detla-y(1) + detla-y(2), for only 3 NM
  % special case: use 4NM code for only 3NM(#2-3-4)
 lfp= (y(8,:) - y(9,:) +y(14,:) - y(15,:) ... % omit NM#1, to get 3
     + y(20,:) - y(21,:) )/3; % av of 3 nodes
 lfp3_sample=lfp(40001:end); % during stimulus [0.5s sample]
   %lfp_sample=lfp(1100:2900); % debug, nb. wider sample for half t-step
 lfp3_stats= [ mean(lfp3_sample), (max(lfp3_sample)-min(lfp3_sample))]
 clear lfp3_*
 
%% 4.1 Align waveforms, by phase: 1st Pulse
% a) at NM#2 [the input]
dy2_sample=y(8, 40000:60000) - y(9, 40000:60000); % NM#2, delta-y 
figure; plot(dy2_sample); hold on; xlim([0 2000]) %axis([0.5 2.5 -50 50])
dym=mean(dy2_sample); plot([0.5 2.5], [dym dym], 'g--'); xlabel('t (0.1ms)')
%
title('JR node #2, 1t pulse: output delta-y: y1(2) - y2(2)') % e - i 
plot(Pd(40000:60000)/10, 'k') % stim % reduce scale [of pulse rate]
plot(Vyt(40000:60000)+dym, 'Color', [1 0 1]) % stim V pulses, wave
legend('delta-y(2)', 'Pulse', 'wave'); 

% a.1)
dy6_sample=y(32, 40000:60000) - y(33, 40000:60000); % NM#6, delta-y 
figure; plot(dy6_sample); hold on; xlim([0 2000]) %axis([0.5 2.5 -50 50])
dym=mean(dy6_sample); plot([0.5 2.5], [dym dym], 'g--'); xlabel('t (0.1ms)')
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

%% 4.1.2 ISI.  {was 3.6}  ISI for y-e(2);  dy(2) / y-e(1);  dy(1)
% fprintf('\n > ISI for y-i(1):  '); 
 %fprintf('\n > ISI for y-e(1): from 3.0 sec '); 
fprintf('\n > ISI for dy(1):  from 3.0 sec ');
%fprintf('\n > ISI for lfp(all NM, t): from 3.0 sec ');
 % signal= y(2, 30000:end);  % y-e(1)  <:: edit, to pick signal
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

       %signal= y(2, 35000:100000)-y(3, 35000:100000);  % short
%signal= y(2, 30000:end)-y(3, 30000:end); disp('  dy(1), beta: \n'); % for NM#1
signal= y(3, 30000:end); disp('  y-i(1), beta: \n'); % for NM#1, inhib popn
 % signal= y(2,:)-y(3, :);  % <:: all dy(1), for NM#1
    %signal= y(2, 1000:end)-y(3, 1000:end);  % dy(1), for 1/2NM
    %signal= y(8,30000:100000)-y(9,30000:100000);  % short
% signal= y(8,30000:end)-y(9,30000:end);  disp('  dy(4): \n'); %for NM#2
%signal =  (y(20,30000:end) -y(21,30000:end)); t=x(40001:end); disp('  dy(4, gamma): \n');% NM-4
% signal= y(32,35000:100000)-y(33,35000:100000);  % dy(6), for NM#6
% LFP LFP:  net lfp output: sum: detla-y(1) + detla-y(2), for the 6 NM
% lfp= (y(2,:) - y(3,:) + y(8,:) - y(9,:) +y(14,:) - y(15,:) ...
%     + y(20,:) - y(21,:)+ y(26,:) - y(27,:) +y(32,:) - y(33,:) )/6; % av of 6 nodes
 %signal= lfp(30000:end)- mean(lfp(30000:end));  % <:: lfp(t), av over all NM
   % signal= lfp(10000:end)-mean(lfp(10000:end)); % short t
 
 
%  % these are peaks of V(t) - max
[SamplePeaks, PeakLocs]= findpeaks(signal); %,'MinPeakheight',10.0); % default 
%[SamplePeaks, PeakLocs]= findpeaks(signal, ...
 %   'MinPeakDistance',100 ,'MinPeakheight',8.0); % for alpha band: 1000; beta 500 / 8.0
% length(PeakLocs)
PeakTimes=PeakLocs*h;  %(1:5)
ISI = diff(PeakLocs)*h;
 ISI  %ISI(1:20)
 Range = [min(diff(PeakLocs)*h) max(diff(PeakLocs)*h)]
 figure; hist(ISI); title('4 NM, ISI istribution'); xlabel(' ISI (sec) ')
 figure; hist(1./ISI); title('4 NM, equiv freq. distribution'); xlabel(' f (Hz) ')
 
% 3 plots together:
 figwidth = 560; figheight = 704; figposition = [500, 100, figwidth, figheight]; % tall
figure('position',figposition, 'units','pixels');  % larger figure; 
subplot(3,1,1); stem(ISI); % xlim([0 50]); %ylim([0.1 0.2]);  % for 6NM, alpha
 ylabel('ISI (sec)'); title('4 MN cluster; node #1: ISI times (sec)')
subplot(3,1,2); stem(ISI);  xlabel('peak # (from 1.0 sec) '); % ylim([0.11 0.22]); % alpha 
ylim([0.07 0.26]); % beta
  grid on; ylabel('ISI (sec)'); xlim([30 length(ISI)]); %xlim([80 160]); % <<<< ??
  % text(7, 0.205, 'v', 'Color', 'r')
 % subplot(3,1,1); text(20, 0.177, 'NM #1, y-e(1);  10 pulses @ 6.35Hz')    % '5 pulses @ 3Hz')
subplot(3,1,3); hold on; title('pulse train'); % text(2.5, 0.25, '6 NM, 0.1mV U ran noise')
 t0 = PeakLocs(1)*h; tt=t0; % (sec)
for i =1:length(PeakTimes)-1
   line([tt tt], [0.0 0.2 ], 'LineWidth', 1.0 );
   tt=tt+ISI(i); % accumulate ISI's to get this time
end
grid on; xlabel(' t (sec +1.0)'); %ylim([0 0.3]); xlim([0 4]);  
%subplot(3,1,3); text(1.0, 0.21, 'v', 'Color', 'r')  % mark stimulus time
%                                                   [@3.5 + 0.5 s;  or  @1.0+3.0 s]
%text(3.0, 0.21, 'v', 'Color', 'r')

% show ISI(t) & Pulse  Train together
figwidth = 1300; figheight = 350; figposition = [100, 380, figwidth, figheight]; % wide
figure('position',figposition, 'units','pixels');  % larger figure; 
subplot(2,1,1); hold on; %stem(ISI);
ylabel('ISI (sec)'); title('4 MN cluster; node #1:  ISI times (sec)')
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
  %text(1.0, 0.21, 'v', 'Color', 'r')  % mark stimulus time
% text(6.0, 0.22, 'v', 'Color', 'r') % marker: wave starts
 % text(6.0, 0.235, 'wave arrives', 'Color', 'r') 
 % subplot(2,1,2); text(1.0, 0.24, '+1 Pulse -> #2]', 'Color', 'r')  % describe stimulu2

%
 t0 = PeakLocs(1) % (0.1 msec)
 Npeaks=length(PeakLocs)
 shortT = mean(diff(PeakLocs))*h  % short period
 StdDev = std(diff(PeakLocs))*h
 AvshortT_freq= [shortT 1/shortT]
%
 figure; plot(signal); hold on; %title('6NM: y-e(1), dy(1);  ISI')
 title('signal V(t)'); xlabel('t (1s + 0.1ms)'); ylabel('V(t)  (mV) ');
plot(PeakLocs, signal(PeakLocs), 'ro'); grid on; % mark peaks - to be sure all id'd
%
figure; hold on
for i=1:Npeaks-1
plot(ISI(i), signal(PeakLocs(i)), 'ko');
 %drawnow
 %pause(0.007) % animate: jumps L/R, U/D
end
grid on; title('4 MN; node #1: y-e(1) vs ISI'); xlabel(' ISI (sec)'); ylabel('y-e(1)(t) (mV) ');

%ISIe = ISI; Spikese=PeakTimes;  %save for later analysis
%ISIi = ISI;  % Pick one
%ISIdy1 = ISI; SpikesDy1=PeakTimes;
%ISIdy2 = ISI; SpikesDy2=PeakTimes;

% save "code" in 10 ms bins, as "binary" array - for later classification
%Code=zeros(1, int32(length(x)/100)); % now 3k long, for 30 sec [not 300k]
%tmp=int32(PeakLocs/100); % place in bins of 100*(0.1 ms) long
%Code(tmp)=1;
 % Code_3Pulses=Code; csvwrite('Code_3Pulses.csv', Code_3Pulses); % save
 % ISI_3Pulses=ISI; csvwrite('ISI_3Pulses.csv', ISI_3Pulses); % for analysis
 clear tmp PeakLocs PeakTimes Npeaks ISI shortT signal StdDev AvshortT_freq Sample* %Code
 
 
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
figure; hold on; title ('phase plot,  y-e (1):  dy1/dt(t) vs y1(t)')
%title ('dy1/dt vs y1(t)') %axis([10 60 10 60]); 
ylabel('dy1/dt(t)' ); xlabel('y1(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ -0 40 -600 1000]); % initiation
%for i = 1:5:5000  % nb transients < 2s
 for i = 15000:5:20000  % at s/s 
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(y(2,i), y(5,i), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end 
 text(95, -600, '10 to 15  sec'); %text(27, 12, '0.1 sec')
%ylabel('dy1/dt' ); xlabel('y1(t)');
hold off
fprintf('\n  plot done!  ')

  %% 4.2e Phase plot y2'(t) vs y2(t) {ie y-i}  for NM#1 : animated
figure; hold on; title ('y-i (1):  dy2/dt(t) vs y2(t)')
ylabel('dy2/dt(t)' ); xlabel('y21(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ 0 40 -500 1000]); % initiation
for i = 140000:5:160000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(y(3,i), y(6,i), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
  text(57, -2100, '14 to 16  sec'); %text(27, 12, '0.1 sec')
fprintf('\n  done!  ') 

 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
hold off

%% Appx. a
% for 8x8 erfc(1/V) form, based on Marm data
rvec= [0.45 0.32 0.38 0.35 0.30 0.38]; % 0.50 0.56];
wvec=[1.08 6.73 1.00 1.54 2.71 4.90]; % 2.05 2.05];
figure; plot(rvec, wvec, '.', 'MarkerSize', 22)
title('Marmoset pFC(6) data')
ylim([0 7]); xlim([0.29 0.64])
xlabel('r'); ylabel('w-bar')
% alt 
figure; plot(wvec, rvec, '.', 'MarkerSize', 22)
title('Marmoset pFC(6) data')
xlim([0 7]); ylim([0.25 0.65])
ylabel('r'); xlabel('w-bar')

clear wvec rvec

%% Appx. 1.x
% get d(i j) dist from marm pFC-6 data
% cf # 1.0.4, of sixCluster codes:  d(ij) array -  set d in bins of 0.1 mm (approx)
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
 list4 = [6, 5, 4, 1]; % corrsponds to A11, A46D, A9, A10 in dbl triangle model
 Adist4 = Ad6A(list4, list4)
 clear list4
 
 
%% Appx cpnts of DE's : fast inhib popn, 
tend=10000; g =1000/taugvec(1);  %g== 167.0  %1/60Hz % prev; G is shared var
 %C=135; C5= 1.0; C6= 1.0; %0.2*C;  %C6 =0; % feedforward: i-slow -> i-fast
 %C7= 0.005*C;  C8= 1.0;  % from JR .. subroutine
C5= 50.0; C6= 0.2; C7= 20.0;  C8= 0.4; % tuning of extra C's
wbar=1.0; v0 =6.0; vf0 =v0; % defaults   %v0= -25; % test?
%
Feedback03= G*g*C7*( vm./(1+exp(r*(v0 -wbar*C5*y(1,800:tend) +wbar*C6*y(3,800:tend)  ...
             -wbar*C8*y(4,800:tend)-wbar*y(5,800:tend) ))) ); % around 0.1 sec
           % in JR code:  -wbar*C5*y(1) +wbar*C6*y(3)  +wbar*C8*y(4) -wbar*y(5)
Drive = Feedback03 -2*g*y(9, 800:tend) -g^2*y(4, 800:tend); % tot for DE
figure; 
 plot(x(800:tend) , -2*g*y(9, 800:tend), 'r--'  ); hold on % fast, damping term, y3' :: V small !
 pause 
 plot(x(800:tend) , -g^2*y(4, 800:tend), 'g-' ); pause;  % osc driving term; y3
 plot(x(800:tend), Feedback03, 'b--');  title('DE driving forces: fast inhib popn'); pause
 plot(x(800:tend), Drive, 'k'); legend( 'damping', 'driving','feedback', 'Total');
  %ylim([-1.5e4 3e4])
 
 figure; plot(x(800:tend), Feedback03, 'b');  title('fast inhib popn :feedback only')
 figure;  plot(x(800:tend), Drive, 'k'); grid on; xlabel('t (sec)'); title('total Drive')
 
 figure; plot(  (v0 -wbar*C5*y(1,800:tend) +wbar*C6*y(3,800:tend)  ...
             -wbar*C8*y(4,800:tend)-wbar*y(5,800:tend)) ); title('nett Cj *dy-s')
         
  figure; plot(  (v0 -wbar*y(1,800:tend) +wbar*y(3,800:tend)  ...
             -wbar*y(4,800:tend)-wbar*y(5,800:tend)) ); title('nett dy-s'); grid on
  %figure; plot(  (y(1,800:tend) +y(3,800:tend)  ...
   %          -y(4,800:tend)-y(5,800:tend)) ); title('basic nett dy-s')
  % 
  %vf = -mean( -wbar*y(1,800:tend) +wbar*y(3,800:tend)  ...
   %          -wbar*y(4,800:tend)-wbar*y(5,800:tend) )
  vf0 = -mean( -wbar*C5*y(1,5000:tend) +wbar*C6*y(3,5000:tend)  ...
             -wbar*C8*y(4,5000:tend)-wbar*y(5,5000:tend) ) % nb. this is s/s, _after_ the sim
         
   %sigmoid= vm./(1+exp(r*(vf0 -wbar*y(1,800:tend) +wbar*y(3,800:tend)  ...
    %         -wbar*y(4,800:tend)-wbar*y(5,800:tend) )))  ; 
   sigmoid= vm./(1+exp(r*(v0 -wbar*C5*y(1,800:tend) +wbar*C6*y(3,800:tend)  ... % weigjted
             -wbar*C8*y(4,800:tend)-wbar*y(5,800:tend) )))  ; 
         %meanDv = -mean( -wbar*y(1,800:tend) +wbar*y(3,800:tend)  ...
          %   -wbar*y(4,800:tend)-wbar*y(5,800:tend) )
         
  figure; plot(sigmoid); title('sigmoid to i-fast; at v0')
  
  clear Feedback03 meanDv sigmoid vf0 % Drive
  
  %% pyrm popn
  % int feedback
r=rvec(in); C=Cvec(in); wbar=wvec(in); % for indiv NM : slope of S[v]; Gain C1,3
 a=1000.0/tauevec(in); b= 1000.0/tauivec(in); % time const for ea. NM (nb. ms)
 A= 3.25;  B =22; 
dyi=zeros(nn,length(x)); % set up array for y3-y2;
onevec=ones(nn,length(x));
  dyi=y(2,:)-y(3,:); % delta-y : feedforward    
  argm=r*(wbar -v0/dyi); Sfn=0;
    if argm >0 & dyi>0.01
       Sfn= (1.0-erfc(argm)/2); % erfc form for wt; avoid v=0
    end          
 Feedback01=  A*a*vm*Sfn  % no I12 to y0
    -(2*a*y(4)) % damping
     -(a^2*y(1))  % ie y0''
   Drive = Feedback01 -2*g*y(9, 800:tend) -g^2*y(4, 800:tend); % tot for DE

%% excit inter n.
C1=C; C2=0.8*C; A= 3.25; 
  yi=y(1); % feedback to both e & i from y0.    
  argm5=r*(wbar -v0/(C1*yi)); Sfn5=0; % nb mean weight
    if argm5 >0 & yi>0.01
       Sfn5= (1.0-erfc(argm5)/2);  % erfc form for N(wt); avoided v=0
    end
    A*a*(Estim + p)% y0, Pulse, ran & I12 to y1 (e)
     A*a*I12 % nn interaction
    Feedback02A*a*C2*vm*Sfn5 %
    -2*a*y(5) 
    -a^2*y(2)
    
    
%% cf  %% feed bk to excit popn
in=4; j11=(in-1)*6+1;   % working on this node #in (4 -> 19]
 % working on : 6 y's @ t-i: sub vec y(j11:j16,i)
 tend =5000;  a=1000.0/tauevec(in);  C=Cvec(in); C1=C; C2=0.8*C; A= 3.25; 
 Feedback01=zeros(1,length(x)); Drive=zeros(1,length(x)); % set up array
 for i=1:tend % scan t [some] steps
    yi= y(j11,i);  %  use y0 from pyrm
    argm=r*(wbar -v0/(C1*yi)); Feedback01(i)=0; % nb mean weight
    if argm >0 & dyi>0.01  % avoid difficulties of erfc(neg, 0)
       Feedback01(i)= A*a*C2*vm*(1.0-erfc(argm)/2); % erfc form for wt; avoid v=0
    end 
 end
 figure; plot(x(800:tend), y(j11,800:tend) ); title('y0( 4)')
 figure; plot(x(800:tend), Feedback01(800:tend), 'b', 'LineWidth', 1.5);  title('y0*( 4)') % S[V] 
 pause
 %
 figure;  hold on; title('4NM, y-e DE driving forces: pyrm [excit] popn');  
 plot(x(800:tend), -2*a*y(j11+4, 800:tend), 'r--', 'LineWidth', 1.5); hold on  % ** probs!!!  slow, damping term y1'  % 
 pause 
 plot(x(800:tend) , -a^2*y(j11+1, 800:tend), 'g--', 'LineWidth', 1.5); % osc driving term; y1
 pause
 plot(x(800:tend), Feedback01(800:tend), 'b', 'LineWidth', 1.0); % of order 
 pause 
 Drive = Feedback01 -2*a*y(j11+4, :) -a^2*y(j11+1, :); % tot for DE
 plot(x(800:tend), Drive(800:tend), 'k', 'LineWidth', 1.0); grid on
 legend('Damping', 'Osc restore', 'Feedback', 'Tot Drive')
 clear A a argm C C1 C2 Feedback01 yi Drive
 
%% ext interactions, to A10 #4; eg from #2 only
  tend =5000; in=4; jj =3
  %j11=(in-1)*6+1; % working on this node #in
 j11=(jj-1)*6+1; % this is the NN [linked]
 r=rvec(in); C=Cvec(in); wbar=wvec(in); % for indiv NM : slope of S[v]; Gain C1,3
 A= 3.25; a=1000.0/tauevec(in); % time const for ea. NM (nb. ms)
    %deltai=int32(Ad6A(jj,in)/(1000*h*vlocal)); % check discontinuities           
    % dyi= y(j11+1,i-deltai)-y(j11+2,i-deltai); % default % delta-y of nn. JR case (input to pyrm)
for i=1:tend % scan t [some] steps
  dyi= y(j11,i); %dyi= y(j11,i-deltai); %  use y0 of nn. pyrm
   %dyi= y(j11+1,i)-y(j11+2,i); % default % delta-y of nn. % dij=0
  argm=r*(wbar -v0/dyi); In12(i)=0; 
    if argm >0 & dyi>0.01  % avoid difficulties of erfc(neg, 0)
       In12(i)= R12*Adj(jj,in)*vm*(1.0-erfc(argm)/2); % erfc form for wt; avoid v=0
    end
end % scan of t steps

   figure; plot(x(800:tend), A*a*In12(800:tend), 'b', 'LineWidth', 1.5);   title('Interaction (4-3) [y0]') 
    % title('Interaction (4-2) [dy 1-2]') % 
   grid on
    
  % clear a A dyi argm Sfn Sfn5 r C wbar a b In12 
   
%% feed fwd to pyrm popn
 tend = 50000; in=4; jj =3 % A =Avec(1); 
 A= 3.25; a=1000/tauevec(1); % default
 in=4; wbar=wvec(in); r=0.5; 
 j11=(jj-1)*6+1; % this is the NN [linked]
  %tmp = C2*y(2,800:tend)- C4*y(3,800:tend) -C7*y(4,800:tend); % in Z('06) *10^2 bigger
 % tmp = y(2,800:tend)- y(3,800:tend);  % orig form dy, of jj
 tmp= y(j11+1,:)-y(j11+2,:); % default % delta-y of NN % dij=0
  %tmp= y(j11,:); % y0  choice
   figure; plot(x(800:tend), tmp(800:tend)); title('y0 '); % title('  nett input, dy to pyrm')
   
  tmp1 = A*a*vm./(1+exp(r*(v0- wbar*(tmp(800:tend) )))); 
   figure; plot(x(800:tend), tmp1); title(' S[dy, nett input]  to pyrm, wbar form')
  
  clear A a tmp tmp1
 
