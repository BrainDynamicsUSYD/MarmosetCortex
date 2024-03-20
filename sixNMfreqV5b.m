% % B {Pailtjorpe, U Sydney. sixNMfreqV5b.m (27/4/23 - 24)
 % Jansen-Rit 6 node cluster at A10 (pFC): 2 theta, 1 alpha, 2 beta, 1 gamma bands
 % uses erfc(1/V) form of V-rate transformation 
% Dependencies: JR12v5d.m in RK4 DE solver loop, incl. pass r & wbar for ea NM
                % WaveTrains.m for stimuli set up 
   % Data files:  AdjMarmoset_Rescale2b.csv, 116x166 link wts, based on LNe 
   
    % V5b. use S[~1/V, erfc (N(w)]; matches twoNM..v5b & sixCluster..v5b (15/4/23)
                % & JR12v5d.m : uses tests for erfc; incl param trials 2, 6 & 'the 10'
% Cells (sections):
% 1. ADj set up
% 2. parameters for WC/JR model
% 3. set up stimuli
% 4.0  RK4 DE solver, 1st order, for n variables
% 5.0  various PLOT OUTPUT
   % Appx 3.2a,b: calc w-in,out / k-in,out for local & ext links
   % Appx 3.3: plot w-bar calc
    
% 1.0 Set up
addpath('/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain/Models'); % for other model codes
fprintf('\n 6NM: 2+1+2+1; theta, alpha, beta, gamma bands; + delays; S[1/V - erfc], V5b code \n')
fprintf('\n      100 pc stim -> e popn \n')
 % Adj, for coupling nodes
nn=6; % number of neural masses (nodes)

% Also: Appx 1, 2 - calc DE driving forces, S[V(t)]

%% 1.01 Adj: zeros (all)
 fprintf('  Adj(ij) = 0: no coupling  \n')
Adj= zeros(nn) % no coupling: indep NM 

%% 1.02 Adj: ones (all)
 fprintf('  Adj(ij) = 1: unit coupling, all-all  \n')
Adj=ones(nn); % unit coupling, symm
 % elim diagonals
 for i=1:nn
     Adj(i,i)=0;
 end
 Adj
 length(find(Adj(:))) % 30 links
   
%% 1.1 Adj from Marmoset (LNe 6x6) data
 % cf.Appx 2. in orig order: cf. logbook, p58
  fprintf('  use marmoset LNe Adj(2): & scale R12=200 (raw Adj/5) \n')
clusterlist=[2, 26, 25, 48, 32, 3]; % as plotted in x-y plane, viewed from Ant.
Anew= csvread('AdjMarmoset_Rescale2b.csv');  % based on LNe 
 Aclust = Anew(clusterlist, clusterlist);
Adj2=Aclust/1000;  % scale used for 6 NM model
Adj=Adj2  % nb use of R12*Adj {R12=100 
R12=200  %R12=100   % NM - NM coupling coeff (ie. link wt) (cf. Goodfellow('11))
          % base coupling: 100 : cf.scaling of Adj wts; applied in DE solver
          % ie. this is Adj/5 : 20% [of expt] link wts
nLinks= length(find(Adj2(:))) % 29 entries
% det(Adj)
clear nLinks Anew Aclust Adj2 clusterlist
 %Adj=zeros(nn); % test indep NM.
 %Adj=ones(nn); % test unit coupled NM. nb. diagonals (self links) present
   % isolate #2,3 [6->4NM]
   % Adj(2,:)=0; Adj(:,2)=0 Adj(3,:)=0; Adj(:,3)=0
    % nb. set approx dij in Ad6A, below & at Appx 3.2.1

%% 1.0.9  zero delays:  dij
 fprintf('  d(ij) =0: no signal delays  \n')
 Ad6A = zeros(nn);
%% 1.0.10  delays:  dij (mm / ms) from marm data 
 % cf # 1.0.4, of sixCluster codes:  d(ij) array -  set d in bins of 0.1 mm (approx)
fprintf('  d(ij) in 0.1mm bins: from Marm. \n')
 %Ad6A=Adist4; % cf appx 1.x below; data for 8x8, use 6x6
Ad6A = [...
     0     3.4  2.9  2.4   2.4  2.5, 4.3, 5.0; ...
     3.4   0    1.7  3.4   4.1  2.7, 0.0, 4.7; ... % 0 if no link present
     2.9   1.7  0    1.9   2.8  2.6, 0.0, 3.2; ...
     2.4   3.4  1.9  0     1.8  3.4, 3.2, 2.8; ...
     2.4   4.1  0    1.8   0    2.8, 2.1, 3.2; ...
     2.5   2.7  2.6  3.4   2.8  0.0  3.6, 4.8; ...
     4.3,  4.9, 3.6, 3.2,  2.1, 3.6, 0.0, 2.4; ...
     5.0,  4.7, 3.2, 2.8,  3.2, 4.8, 2.4, 0.0];   % checks ok.
 tmp =Ad6A(1:6,1:6); max(tmp(:))
 Ad6A = tmp;
 % figure; hist(tmp(:)); title('pFC(6* cluster) d(ij) distribution'); xlabel('d(i j)  (mm  / ms)')
 %figure; hist(Ad6A(:)); xlim([1 5]); % isolate "0"
  % isolate #2,3:  Ad6A(2,:)=0; Ad6A(:,2)=0 Ad6A(3,:)=0; Ad6A(:,3)=0
clear tmp

%% 2.0 Parameters: tune freq. band
fprintf('\n 2+1+2+1 NM; 6 star cluster, set parameters \n')
fprintf('\n   Tuning #6, for w-av   \n')
% time axis & IC (0 mV)
x0=0; xf= 30.0; % t range (sec)
h = 0.0001;   % t step (1 ms; [& 0.5, 0.25 ms, for debug])
x = [x0:h:xf];  % t domain array
 %time =x/(1000);  %time =x/(1000*h); % time (sec) as row vec
tstart=4001; %tstart=2000; %1000; % for plots (@ 1, 0.5, 0.1 ms time steps)
%y = zeros(6,length(x)); % row vec for each variable
yic1=[1.0; 2.0; 1.0; 1.0; 1.0; 0.0]; % IC at t=0 for 1 node (6 DEs),  col vec, 6 rows
%yic=5*yic1; % add additional nodes; and stronger?
yic= yic1;
for in=1:nn-1
 yic=[yic; yic1]; % add 6 rows per extra node; 12, 18 DEs for 2, 3 ... nodes
end
 clear yic1
%
% Model parameters (orig JR) : pass to fn subroutine! cf wb7, p66
v0=6;    % (mV) midpoint (switching threshold V) of sigmoid function, for V -> pulse conversion 
vm=5;    % ie. 2*2.5 (s^-1) max firing rate for sigmoid fn
rvec= [0.5 0.5 0.5 0.5 0.5 0.5] % v2+ revised params (wb7, p57 & 84) % default
 % & alt: from prev 6NM [single tau] models
%rvec= [0.56 0.40 0.47 0.43 0.38 0.47] % revised r ~ sqrt(N); cf w/b#4,p42 & 5,p72: tune in BETA band
   %rvec= [0.42 0.30 0.37 0.34 0.30 0.37] % ditto: tune in alpha band: corrected (23/7/21)
   %rvec= [0.45 0.32 0.38 0.35 0.30 0.38] % tune in ALPHA (or theta) band: wb 5, p70 (18/9/21)
% feedback: tuned to bands
Cvec = [220 200 200 220 200 250] % retuned  gamma/ theta/ beta/ alpha/ theta/ beta/; wb7, p86 
 %Cvec = [220 250 300 220 200 250] % retuned #2,3 to osc ok
% tune feedback gain (based on WtIn, kIn)
CfacVec= [1.0 1.0 1.0 1.0 1.0 1.0]; % pass to fn J12a.m integrator
    % nb. Assign in J12: C1=C1*Cfac;  C3=C3*Cfac; 
% weight/synp link, for Sigmoid response fn:
  %wvec= [1.0 1.0 1.0 1.0] % default case: basic S[V] sigmoid 
  %wvec= [3.0 1.5 3.0 1.5 1.5 3.0] % vary f-i: tuning #1, wb[7],p85
  %wvec= [3.0 3.0 3.0 1.5 3.0 3.0] % vary f-i: tuning #1a, wb[8],p40
  %wvec= [3.0 5.0 3.0 1.5 4.0 3.0] % vary f-i: tuning: guesses #1a, b, c, d wb[8],p40
 %wvec= [1.0 5.0 1.0 1.5 2.7 4.8] % check wt-in/k-in : "Trial 7" (14/12/23)

wvec= [3.0 5.0 3.0 1.5 3.5 3.0] % nb dominant link #16 wb[8],p42 : ** "Trial 6" **
  %wvec= [3.0 5.0 3.0 1.5 3.0 3.0] % nb change #5 wb[8],p69 : "Trial 9"
 %wvec= [3.0 5.0 3.0 1.5 5.0 3.0] % reassign for theta wb[7], p102, orig f-i tuning: "Trial 2"
 %wvec= [3.0 2.5 3.0 1.5 3.5 3.0] % variations  ~ theta [8], p37
   %wvec=[0.52 3.23 0.48 0.74 1.30 2.35] % for pFC; raw wt/link: av in each NM: wb4,p59a
    %wvec=[2.6 1.74 1.14 0.92 0.72 1.41] % based on Av[w-in/k & w-out/k] for 6x6: wb5,p102a
  %wvec=[4.24 2.84 1.86 1.5 5.84 2.30] % rescale so #4 matches: [7], p20a
  %wvec=[3.4 2.3 1.5 1.2 1.0 1.8] % [orig] rescaled Av[w-in/k & w-out/k] wb5,p102a; [8], p64
%wvec=[3.38 2.35 1.48 1.2 1.0 1.84] % 3a. expt Av[w-in/k & w-out/k] wb[8], p77; mss Table S2
    %wvec=[2.6 3.2 1.8 1.5 2.0 2.4] % mix of max or av, Trial 6.1 [8],p74 
   %wvec=[1.9 1.26 1.17 1.11 2.69 1.1] % and fm 8x8 calc
% signal delay: set by array Ad6a, above; cf. App3.2, based on d(i,j) % d12=1.5-4;  % eg. mm / ms
% Rate params - pass to JR12a fn.
 A= 3.25;  B =22;  % Max PSP amplitude for pyramidal (A), inhib (B) & excit (kA * A) interneurons
   %a=100; b=50;  % (s^-1) / default JR; set to osc at 8Hz (taus: 20, 20 mS)
    %taue= 15; taui= 20; % default Time Const (excit, inhib) (ms)  :  alpha, beta bands
    % retune in DE loop:  keep a*A & ? b*B const  %a=1000/taue; b=1000/taui; % rate const (s^-1)
tauevec=[5.0 25.0 10.0 16.0 25.0 10.0] % revised again (15/4/23): wb7, p68 & 8;  % #2 was 25
tauivec=[5.0 25.0 10.0 17.0 25.0 10.0] % for NM 1-6
Drive_1= [1000*A/tauevec(1)  1000*B/tauivec(1)] % check "driving forces"
Drive_6=  [1000*A/tauevec(6)  1000*B/tauivec(6)]
% Coupling const for S[v]: cf. Steigler(2011)
 kIn2=1; kIn6=1; 
 clear Drive_*
 
% >>>>>>>>>>>>>>
%% 3.0 Zero stim - pulse, ran or wave :: Base case
fprintf('\n >> zero inputs: 1 x G ran noise \n') % nb scale units are Hz, not "mV"
Vyt =zeros(1,length(x));
pran=zeros(1,length(x)); 
IpOn=zeros(nn,1);  % switch on/off to ea NM for stim inputs; need IpOn also
IpOnB=zeros(nn,1); % 2nd pulses, default
Pd=zeros(1,length(x));  Pd2=zeros(1,length(x)); %  pulse density fn - cf #3.6
 %pran = 0.1*rand(1,length(x)); % low noise;  Uniform ran noise:
pran = 1.0*randn(1,length(x)); % G ran noise (+-); [Hz]
  %pran =7.75*pran; % test:  sqrt(60) std, a la W & Z (4/5 popn) models
    % 
InI =0  % controls stim to Excit/Inhib popn:  default - apply stim to Excit popn
%InI =1 % apply stim to Inhib


%% 3.5 Constant stim rate - use Pd; 
 % nb.not filtered by S[v]
 % explore stable pts and limit cycles.
fprintf('\n ++const stimulus: on at 4.0 sec ')
tstart=40001; % also used for plots
StimLevel= 100
Pd=zeros(1,length(x)); tstart0= tstart; %tstart0= 1001;  %separate for Pd
 %tend=60001; % for short stim [4:5sec] etc
 %Pd(tstart0:tend)=StimLevel; % short, const stim
Pd(tstart0:end)=StimLevel; % && turn on the const stim!!
PulseIntegrl = sum(Pd)*h % sum across rows % trapezoidal rule
  %Ip6=zeros(1,length(x)); % pulse train/fn for 2nd input (?dphi shifted) - was to node 6
 InI =0  % controls stim to Excit/Inhib popn:  apply stim to Excit popn (default)
 %InI =1 % apply stim to Inhib 

clear StimLevel PulseHeight tstart0 PulseIntegrl tstart0

 %% 3.5.1 Modulated stim rate - use Pd x |Vyt|, full wave rectification; 
% modulate by Vyt fn - set up in WaveTrains.m;  not filtered by S[v]
%fprintf('\n ++ modulated stimulus - via  Vyt(f/2): at 4.0 sec ')
fprintf('\n ++ modulated stimulus - const * |Vyt(f)|: at 4.0  sec ')
 fprintf('\n    full wave rectification ( |abs| ) \n')
tstart=40001; % start stim & for plots
StimLevel= 100.0
Pd=zeros(1,length(x)); %tstart0= tstart; %tstart0= 1001; %tstart0= 20001; %  separate for Pd
 %tend=140001; % for short stim [4:5 sec] etc
 %Pd(tstart:tend)=StimLevel*abs(Vyt(tstart:tend)); % short, const*modulated stim
Pd(tstart:end)=StimLevel*abs(Vyt(tstart:end)); % modulated stim! pos. via abs(), so 2*f
  %Pd(tstart:end)= StimLevel + 0.1*StimLevel*Vyt(tstart:end); % modulated stim! pos. via abs(), so 2*f

StimStats= [mean(Pd) std(Pd) max(Pd)]
PulseIntegrl = sum(Pd)*h % sum across rows % trapezoidal rule
 %test=Pd(tstart:end); %test=Pd(tstart:tend);
 %StimStats1sec= [mean(test) std(test) max(test)]

%Vyt =zeros(1,length(x)); % reset Vyt: no V wave
 
 % PdSave=Pd; % for reuse
 % Pd=zeros(1,length(x));  tstart=40001; tend=47001; % reset
 %  Pd(tstart:tend)=PdSave(tstart:tend); % short, modulated stim
 
 InI =0  % controls stim to Excit/Inhib popn:  apply stim to Excit popn (default)
 % InI =1 % apply stim to Inhib 
 
clear StimLevel PulseHeight tstart0 tend PulseIntegrl StimStats test


%% 3.5.1a Vyt pos only, Half Wave rectification, Modulated stim rate - use Pd x Vyt
% modulate by Vyt fn - set up in WaveTrains.m
fprintf('\n ++ modulated stimulus, pos only (f): at 4.0 sec ')
 fprintf('\n    Half wave rectification (V>0 only) *2 {const. E \n')
tstart=40001; % start stim & for plots
StimLevel= 200
Pd=zeros(1,length(x)); 

Vpos=zeros(1,length(x));
tmp=find(Vyt>=0)';
Vtmp=Vyt(tmp)'; % pos only
  % figure; plot(tmp,Vtmp); xlim([3.9e4 4.4e4]) % debug
Vpos=zeros(1,length(x));
Vpos(tmp)=Vtmp(1:end)';
 % length(xpos) is  80027 for [4.0:20] sec
 % figure; plot(x, Vpos); xlim([3.9 4.4]) % check:  ok
Pd(tstart:end)=StimLevel*Vpos(tstart:end); % modulated stim! strictly pos.
 % tend=100001; % for short stim [4:14 sec] etc
 % Pd(tstart:tend)=StimLevel*Vpos(tstart:tend); % alt. short, const*modulated stim

StimStats= [mean(Pd) std(Pd) max(Pd)]
PulseIntegrl = sum(Pd)*h % sum across rows % trapezoidal rule
 % Vyt =zeros(1,length(x)); % reset Vyt: no wave
 
clear StimLevel PulseHeight tstart0 PulseIntegrl StimStats tmp* Vtmp Vpos

 %% 3.5.2 Modulated stim rate : const + Pd x Vyt, raw wave;
% modulate by Vyt fn - set up in WaveTrains.m
%fprintf('\n ++ modulated stimulus - via  Vyt(f/2): at 4.0 sec ')
fprintf('\n ++ modulated stimulus = const1 + const2 * Vyt(f): at 4.0 sec ')
 fprintf('\n    raw wave, const offset ) \n')
tstart=40001; % start stim & for plots
StimLevel= 100.0
Pd=zeros(1,length(x)); %tstart0= tstart; %tstart0= 1001; %tstart0= 20001; %  separate for Pd
 %tend=50001; % for short stim [4:5 sec] etc
 %Pd(tstart:tend)=StimLevel*abs(Vyt(tstart:tend)); % short, const stim
Pd(tstart:end)= (StimLevel+0.1)/2 +(StimLevel/2)*(Vyt(tstart:end)); % apply offset, so Pd >0 strictly
 %Pd= (StimLevel+0.1)/2 +Pd; % modulated stim! pos. via abs(), so 2*f

StimStats= [mean(Pd) std(Pd) max(Pd)]
PulseIntegrl = sum(Pd)*h % sum across rows % trapezoidal rule
 %test=Pd(tstart:end); %test=Pd(tstart:tend);
 %StimStats1sec= [mean(test) std(test) max(test)]

 %Vyt =zeros(1,length(x)); % reset Vyt: no V wave
 
 figure; plot(Pd); hold on;  xlim([3.9e4 4.4e4]); % check form
 avline = max(Pd)/2;
 line([4.0e4 4.4e4], [avline avline]); %, '--', 'Color', 'g'); % ?? 
 title('Stimulus, modulated + bias');
 
 % PdSave=Pd; % for reuse
 % Pd=zeros(1,length(x));  tstart=40001; tend=47001; % reset
 %  Pd(tstart:tend)=PdSave(tstart:tend); % short, modulated stim
 
 InI =0  % controls stim to Excit/Inhib popn:  apply stim to Excit popn (default)
 % InI =1 % apply stim to Inhib 
 
clear StimLevel PulseHeight tstart0 tend PulseIntegrl StimStats test avline

%% 3.6.1 Load 1st Pulse density function[s]
 % a la JR('93 & '95); set up in code WaveTrain.m @ #1.4
 fprintf('+ Load Pulse density stim(@ 4s) + 0 (@0.5s); & 1 mV G ran noise: \n')
 Pd= 1.0*Pulse1; % 1 or 2 pulses;  vary amplitude (ie. pulse density
 tstart=40001;  % cf. wave... code, for pulse setup
  % Pd= 1.0*PulseN; % N pulses, in train
 %Pd = Pd+100.0;  % add const stim (all t)
 Pd((1001):end) = Pd((1001):end)+ 0.0;  % add const stim , with pulse (not before)
  PulseStats= [mean(Pd) std(Pd) max(Pd)]
  %PulseStatsAtStim = [mean(Pd(40001:50000)) std(Pd(40001:50000)) ]
  PulseIntegrl = sum(Pd)*h; % sum across rows % trapezoidal rule 
  PulseIntegrl=PulseIntegrl'

 %Pulse2 =zeros(1,length(x)); % default: no 2nd pulse
 InI =0  % controls stim to Excit/Inhib popn:  apply stim to Excit popn (default)
  %InI =1 % apply stim to Inhib 
  
  clear PulseHeight PulseIntegrl PulseStats* 

 %% 3.6.1 Load 2nd Pulse density function[s]
 fprintf(' add in 2nd Pulse fn stim: \n')
 StimLevel= 1.0 %100.0
  %Pd2= 1.0*Pulse2; % 1 or 2 pulses;  vary amplitude (ie. pulse density
  %Pd2 = Pd2+155.0;   % add const stim
  Pd2 =  StimLevel*Pulse2; % add 2nd pulse
   %Pd = Pd + StimLevel*Pulse2; % alt:  add 2nd pulse to orig pulse  stream
 PulseIntegrl = sum(Pd2)*h; % sum across rows % trapezoidal rule 
  PulseIntegrl=PulseIntegrl'
 IpOnB =zeros(nn,1); % switch 2nd stim on/off to ea node - tbd
 % IpOnB=ones(nn,1); % 2nd stimulus On for ALL nodes
  IpOnB(1)=1; % A-10 gamma; [Out hub
 % IpOnB(2)=1;  % [node 2] A-32V  theta  [In hub
 % IpOnB(3)=1;  % beta [3] A32 
 %  IpOnB(4)=1;  % alpha {A9}
 % IpOnB(5)=1;  % theta {A46-D
%  IpOnB(6)=1;  % beta  {A11 orig } [another In hub
  IpOnB' % read out 
  clear PulseIntegrl
  
%% 3.6.2 Load 2nd wave train; Full wave retification
fprintf(' add in 2nd modulated wave stim at: 5.0:5.1 s \n')
Pd2=zeros(1,length(x)); % new array for 2nd stim
 StimLevel= 100.0 %100.0
 tstart=100001; % start stim & for plots
 tend=  101001; % for short stim % 100 ms, et
 Pd2(tstart:tend)=StimLevel*abs(Vyt2(tstart:tend)); % short, modulated stim
 %Pd(tstart:end)=StimLevel*abs(Vyt(tstart:end)); % modulated stim! pos. via abs(), so 2*f
%Pd(tstart:end)= StimLevel + 0.1*StimLevel*Vyt(tstart:end); % modulated stim! pos. via abs(), so 2*f

StimStats= [mean(Pd2) std(Pd2) max(Pd2)]
PulseIntegrl = sum(Pd2)*h % sum across rows % trapezoidal rule
 test=Pd2(tstart:end); %test=Pd(tstart:tend);
 StimStats1sec= [mean(test) std(test) max(test)]

%Vyt2 =zeros(1,length(x)); % reset Vyt2: no V 2nd wave left
% set targets 
 IpOnB =zeros(nn,1); % switch 2nd stim on/off to ea node - tbd
 % IpOnB=ones(nn,1); % 2nd stimulus On for ALL nodes
 % IpOnB(1)=1; % A-10 gamma; [Out hub
  IpOnB(2)=1;  % [node 2] A-32V  theta  [In hub
 % IpOnB(3)=1;  % beta [3] A32 
  % IpOnB(4)=1;  % alpha {A9}
 % IpOnB(5)=1;  % theta {A46-D
  IpOnB(6)=1;  % beta  {A11 orig } [another In hub
  IpOnB' % read out 
 
clear StimLevel tstart0 tend PulseIntegrl StimStats   

%% 3.6.3 Load 2nd wave train; Half wave retification
fprintf(' add in 2nd modulated wave, half wave, stim at: 10.0: 10.1 s \n')
Pd2=zeros(1,length(x)); % new array for 2nd stim
 StimLevel= 200.0 %100.0
 dphi = 200 % set here
 tstart=100001+dphi; % start stim [t steps of 0.1 ms] 
  tend= 101001+dphi; % for short stim % 100 ms, etc. Nb. need to pass dphi
  tstart2= 102001+dphi; tend2= 103001+dphi; % 2nd of 100ms pulse
  tstart3= 104001+dphi; tend3= 105001+dphi; % 3rd of 100ms pulse
 Vpos=zeros(1,length(x));
tmp=find(Vyt2>=0)'; % use 2nd wave loaded
Vtmp=Vyt2(tmp)'; % pos only
  % figure; plot(tmp,Vtmp); xlim([3.9e4 4.4e4]) % debug
Vpos(tmp)=Vtmp(1:end)';
 % length(xpos) is  80027 for [4.0:20] sec
 % figure; plot(x, Vpos); xlim([3.9 4.4]) % check:  ok
 %Pd2(tstart:end)=StimLevel*Vpos(tstart:end); % Basic option; continuous from tstart
 Pd2(tstart:tend)=StimLevel*Vpos(tstart:tend); % short, modulated stim strictly pos.
  Pd2(tstart2:tend2)=StimLevel*Vpos(tstart2:tend2); % 2nd burst
  Pd2(tstart3:tend3)=StimLevel*Vpos(tstart3:tend3); % 3rd burst
   % %Pd2(tstart:end)=StimLevel*abs(Vyt(tstart:end)); % modulated stim! pos. via abs(), so 2*f
   %Pd2(tstart:end)= StimLevel + 0.1*StimLevel*Vyt(tstart:end); % biased, modulated stim! pos. via abs(), so 2*f

StimStats= [mean(Pd2) std(Pd2) max(Pd2)]
PulseIntegrl = sum(Pd2)*h % sum across rows % trapezoidal rule
 %test=Pd2(tstart:end); %test=Pd(tstart:tend);
 %StimStats1sec= [mean(test) std(test) max(test)]

%Vyt2 =zeros(1,length(x)); % reset Vyt2: no V 2nd wave left
% set targets 
 IpOnB =zeros(nn,1); % switch 2nd stim on/off to ea node - tbd
  IpOnB=ones(nn,1); % 2nd stimulus On for ALL nodes
  IpOnB(1)=1; % A-10 gamma; [Out hub
 % IpOnB(2)=1;  % [node 2] A-32V  theta  [In hub
 % IpOnB(3)=1;  % beta [3] A32 
  % IpOnB(4)=1;  % alpha {A9}
 % IpOnB(5)=1;  % theta {A46-D
 % IpOnB(6)=1;  % beta  {A11 orig } [another In hub
  TargetNodes=IpOnB' % read out 
 
clear  tstart0 tend PulseIntegrl StimStats tmp  Vtmp TargetNodes % StimLevel Vpos 
  
 %% 3.7 Load Spike train fn - eg. from from LIF sim
  % Load eSR (10ms bins) into eSpikes array in WaveTrains #1.6; from COBN.m .c
 fprintf(' add in LIF spike train stim x frn: \n')
 Pd=zeros(1,length(x)); 
 %efrStim=csvread('efrStim6_300.csv' ); % read saved file, from LIF sim (M=300, noise 6)
 efrStim=csvread('efrStim2_100.csv'); % N=100, frn 0,2, noise 2 ~ that[&beta] dominated
  %efrStim=csvread('efrStim6_300_40s.csv' ); % 40s sim
 efrStim = [efrStim', 0]; % match length to local variables
   %efrStim = [efrStim', efrStim', 0]; % for 40s
 Pd= Pd + 400.0*efrStim; % x scale; shape matches arrays lfp, x, y etc  [col vec]
 
  % another Alt. use use COBN output eFR (in 0.1ms bins)
     %Pd=zeros(1,length(x)); % correct length
     %Pd(2:end) = eFR; % matches lengths
 PulseHeight= max(Pd)
 PulseIntegrl = sum(Pd)*h % sum across rows % trapezoidal rule
[mean(Pd) std(Pd) max(Pd)]

 clear PulseHeight PulseIntegrl
 
 %% % 3.7.1 Load Spike train fn from from LIF sim x modulate f-band
  % envelope of the wave; at double period (ie half cycle)
fprintf('  + envelope spike train  x100: x gamma 32.1 Hz, at 4.0 sec   \n')
%fprintf('  + envelope spike train: x alpha 9.6 Hz   \n')
%fprintf(' LIF spike train x envelope @ gamma 32.0 Hz   \n')
%fprintf('  + envelope spike train: x 1 sec osc   \n')
Pd=zeros(1,length(x)); 
efrStim=csvread('efrStim6_300.csv' ); % read saved file, from LIF sim (N=300, noise 6) ~gamma
%efrStim=csvread('efrStim2_100.csv'); % N=100, frn 0,2, noise 2 ~ that[&beta] dominated
  %efrStim=csvread('efrStim6_300_40s.csv' ); % 40s sim
 efrStim = [efrStim', 0]; % match length to local variables
   %efrStim = [efrStim', efrStim', 0]; % for 40s 
Pd=zeros(1,length(x)); 
  % Pd= 100.0*efrStim; % x scale;  shape matches lfp, x, y etc  [col vec]
  % Pd(40001:end) = 250.0*efrStim(40001:end);  % 4 sec : or full range t;  
 %xb=1/5.9;  % f-theta modulation: full rectification
 % xb=2.0*xb; % full rectification 
 %xb=1/9.6;  % f-alpha modulation
 %xb=1/14.3;  % f-beta modulation
xb=1/32.1;  % f-gamma modulation {32.8 for 4NM, 32.1 for 6NM}
 %xb=2.0*xb; % full rectification - ie. appliy f/2
   %xb=5.0; %xb=1; %xb=1/6.25;  %xb=1; %period of burst modulation (sec); 
dphi=0;  %40001+ 0;
 Wenvelope=1.0*sin(2*pi*x/(2*xb) -dphi);
  figure; plot(abs(Wenvelope)); title('envelope'); xlim([5e4 5.6e4]) % debug
       % Note - to set scale here:
tmp= abs(Wenvelope).*efrStim*400.0; % pos only [rectified] - full time; or on at 4s
   %figure; plot(tmp); %xlim([5e4 5.6e4]) % debug
%Pd=tmp; % full time span
Pd(40001:end)=tmp(40001:end); % start at 4.0 sec
  %  Pulse1(40001:50000) = tmp(40001:50000); % 1sec burst [4:5] sec
   %Pd((1001):end) = Pd((1001):end)+ 0.0;  % add const stim , with pulse (not before)
    figure; plot(Pd); xlim([5e4 5.6e4]) % debug
      % figure; plot(x,Pd); % debug
PulseStats= [mean(Pd) std(Pd) max(Pd)]
  %PulseStatsAtStim = [mean(Pd(40001:50000)) std(Pd(40001:50000)) ]
PulseIntegrl = sum(Pd)*h; % sum across rows % trapezoidal rule 
PulseIntegrl=PulseIntegrl'
 InI =0  % controls stim to Excit/Inhib popn:  apply stim to Excit popn (default)
  %InI =1 % apply stim to Inhib 

clear tmp Pulse1 PulseHeight PulseIntegrl PulseStats* Wenvelope 
 
 %% 3.7.1a Load current sim LIF output
  eFR =[double(eFR); 0.0]; % match length
  Pd= Pd + 10.0*eFR'; % resampled LIF output
  
%% 3.7 Input channels [of 1st pulse, wave]
 % selected nodes to get the stimulus inputs: Pd, etc < 1st stim 
 %IpOn=ones(nn,1); % 1st stimulus On for ALL nodes
 IpOn=zeros(nn,1); % switch 1st stim on/off to ea node
 % stimulate Node 2 (A-32V): strongest, more wt, fewer links)
 % IpOn(1)=1; % A-10 gamma;  cf wb7, p85
  IpOn(2)=1;  % [2] A-32V  theta 
 % IpOn(3)=1;  % beta [3] A32 
 %  IpOn(4)=1;  % alpha {A9}
 % IpOn(5)=1;  % theta {A46-D
 IpOn(6)=1;  % beta  {A11 orig }
  StimTo_1 =IpOn' % check
 % StimTo_2 =IpOnB' 
    % direct 2nd pulse, if present? : set at #3.6.1 above
   % IpOnB=zeros(nn,1); % switch 2nd stim on/off to ea node; if needed
      %IpOnB(2)=1; % 2nd Pulse to NM#2 (1st in Hub A32-V):
      %IpOnB(6)=1;  % 2nd Pulse to 2nd input (A11)
      %IpOnB' % check  
 clear StimTo*
 
 %% 3.8 Gather stimuli & Plot, to check stim
 figure;  plot(pran, 'Color', [0.9 0.9 0.9]);  hold on
 plot(Pd, 'b'); title('stimulus pulse fn + wave + noise');
 plot(Pd2, 'k');
  plot(Vyt, 'g--'); 
 xlabel('t steps (0.1 msec) ')
   if max(Pd+pran+ Vyt) < 5
     ylim([-2 2.2])  %ylim([0 1])  
   else
     ylim([-2 (max(max(Pd +pran)) +20)])
   end
     % xlim([3.9e4 5e4]) % focus on time of stim
   confirmTo_Excit0Inhib1 = InI
   text(1.05e5, (max(Pd)+9), strcat('Stim To  :   ', num2str(IpOn')) )
   text(1.05e5, (max(Pd)+5), strcat('      and  :   ', num2str(IpOnB')) )
    % xlim([0.5e5 0.6e5]) % to examine details
    % xlim([0 2e5]) %  stop at 20s 
   clear confirmTo_*
% >>>>>>>>>>>>>> 

%% 4.0  RK4 DE solver, 1st order, for 6n variables
% nb. matlab has default double precision: needed here
fprintf('\n > start sim, DE solver:  ')
 % signal delay, usually assumes v = 1 m/s; (mm/ms)
 vlocal= 1.0;  % unmyeln axon vel (~ signal delay) ~ 0.1-1.0 m/s 
tic
% reset arrays
y=zeros(nn*6,length(x)); % set up array for solutions; nn nodes (x 6 rows each)
 %y=double(y); % enforce double precision (64 bit) : should be default
deltay=zeros(nn,length(x)); % set up array for ye-yi; 6 DEs (1 row ea)
 %deltay=double(deltay);
y(:,1)= yic ; % y's are row vec, one for ea DE; IC is a col vec, for i=1 or t=0. 
  for i=1:50  % may need to escape "trap" at zero
   y(:,i)= yic; % extend IC - for delay DE "late start" - eg to ~1-5 ms
  end
In12=0;  % reset here; and for ea  NM
% RK4 DE solution step:  iterate forward in time (here x); scan nodes(in) & linked NN (j)
  % fn JR12 incl forcing fn. as arguments, at this t-step; 
  % inputs: pran + p(i,j)+ [NM-NM couplings] {= sum(y1-y2) ; 
       %fprintf('\n   i  j wt(ij) dr(ij) index r Cfac ') % debug- header
        % start "late" to allow for delay back ref, to L [allow 5ms w 0.1ms steps]
for i = 50:(length(x)-1) % scan time steps & start at 1, 10, 100 [for delay DE]
    for in =1:nn % scan nodes in cluster (NM)
         j11=(in-1)*6+1; j16=(in-1)*6+6;  % working on this node #in
         r=rvec(in); C=Cvec(in); Cfac=CfacVec(in); wbar=wvec(in); % for indiv NM : slope of S[v]; Gain C1,3
         a=1000.0/tauevec(in); b= 1000.0/tauivec(in); % time const for ea. NM (nb. ms units)
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
            % signal delay assumes v = 1 ms/ (mm/ms): set above
            % vlocal=1.0; % unmyeln axon vel ~ signal delay ~ 0.1-1.0 m/s 
           %deltai=int16(Ad6A(jj,in)/(1000*h*vlocal)); % adjust index, for signal delay time from nn, jj->in / as integer
           deltai=int32(Ad6A(jj,in)/(1000*h*vlocal)); % check discontinuities
           %deltai=Ad6A(jj,in)/(1000*h); % Old: adjust index, for signal delay time from nn, jj->in
           dyi= y(j11+1,i-deltai)-y(j11+2,i-deltai); % delta-y of nn.
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
         Pdstim=Pdstim+IpOnB(in)*Pd2(i);% switch 2nd stim pulse on, if present?; 
          %pstim=pstim+  vm./(1 +exp(r*(v0- wbar*Vyt(i))) ); % wave stim, to both e & i; via p-ran
            % nb/ here InI replaces Ip, and is 0; feed-fwd,bk calc in JR fn.
        k_1 = JR12v5d( x(i), y(j11:j16,i), Pdstim, InI, pstim, In12, a, b, A, B, C, Cfac, r, wbar);  % these k's should be col vec, one for ea DE
        k_2 = JR12v5d( x(i)+0.5*h, y(j11:j16,i)+0.5*h*k_1, Pdstim, InI, pstim, In12, a, b, A, B, C, Cfac, r, wbar); % seem to need to force col vec ?
        k_3 = JR12v5d( (x(i) +0.5*h), (y(j11:j16, i) +0.5*h*k_2), Pdstim, InI, pstim, In12, a, b, A, B, C, Cfac, r, wbar);
        k_4 = JR12v5d( (x(i)+h), (y(j11:j16,i) +k_3*h), Pdstim, InI, pstim, In12, a, b, A, B, C, Cfac, r, wbar);
        y(j11:j16,i+1) = y(j11:j16,i) + (k_1 +2*k_2 +2*k_3 +k_4)*(h/6); % load y1:y6 at this t-step      
    end %  scan nodes
end  % scan t-steps
fprintf('\n  calc done/  ') % for debug
toc %timing
 clear argm i in deltai In12 j* k_ a b Istim pstim pdstim vlocal 
 
%% 5.0 PLOT OUTPUT 
% Node 1 output [this is centre of star cluster, and Out-Hub (A-10)
  tstart=40001; % 4 sec
  %figure; plot(y(1:3,:)'); title('JR/ node-1: y0, y1, y2')  % nb. need cols
  %title('A-10:  JR/ node-1, outputs:'); 
  %legend('y0', 'y1', 'y2'); % orig notation [not array index]
    %figure; plot(x, y(1,:));  title('6NM, node-1, pyrm: output: y-e-pyrm '); xlabel('t (s) ')
% title('A-10:  JR/ node-1, output: detla-y: y1(1)- y2(1) '); 
 figure; subplot(3,1,1); plot(x, y(2,:)'); hold on; plot(x, y(3,:)'); 
 title('6NM, node-1: output: y-e, y-i '); 
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
   clear dy1Sample dy1Up dy1Down dy1Range*
   
% Node 2 outputs [this has the primary inputs]:
  figure;  plot(x, (y(8,(1):end) - y(9,(1):end)) ); 
  title('6NM Node 2, A32V;  delta-y ')
  xlabel('t (sec) ');  ylabel('dy2(t) = y-e - y-i  (mV) ')
 % text(3.5, 50, 'On'); text(4, 48, '|'); text(4, 46, '|')
 %  text(11, 40, 'stim 100Hz * gamma, full wave')
 
      %text(4.0, 50, 'v', 'FontSize', 12); text(4.05, 52, '|', 'FontSize', 12);
     %text(4.4, 50, 'Gamma modulated stimulus on')
 %
 % text(15e4, 45, '6NM, S[erfc], 100Hz x (theta/2)'); text(4e4, -37,'\^'); text(4e3, -37, '|')
 figure; subplot(2,1,1); plot(x,y(8,:)'); hold on; plot(x,y(9,:)'); 
 title(' node-2, A32V: output: y-e, y-i '); 
 legend('y-e (2)', 'y-i (2)', 'Location', 'SouthEast'); 
 ylabel('y-e, y-i  (mV)' ); hold off
 subplot(2,1,2); plot(x,y(8,:)'); hold on; plot(x,y(9,:)'); 
  plot(x, Pd/2); plot(x, Pd2/2, 'k'); % add stim
  legend('y-e (2)', 'y-i (2)', 'Stim/2') 
  %axis([1.5 2.5 20 100]);
  xlim([3.9 5.4]); grid on; xlabel('t (s) ')
     %
     % other outputs
     %figure; plot(y(6,(tstart-100):end));  title('JR model, Node-1 output, y5  ')
     %figure; plot(y(12,(tstart-100):end));  title('JR model, Node-2 output, y5  ')
     % axis([0 5 -500 0]);  axis([0 20  -10000 10000]);     

 % Node 3 outputs:  
 figure; plot(y(14,1:end) - y(15,1:end));
 title(' node-3 output, dy: y1(3) - y2(3)')
 figure; subplot(2,1,1); plot(x, y(14,:)'); hold on; plot(x, y(15,:)'); 
 title(' node-3 output: y1, y2 '); legend('y1(3)', 'y2(3)'); hold off
 subplot(2,1,2); plot(x, y(14,:)'); hold on; plot(x, y(15,:)');  
 xlim([3.9 5.4]); grid on; xlabel('t (s) ')  
% Node 4 outputs:
 figure; plot(y(20,1:end) - y(21,1:end));
  title('JR model: node-4 outputs, dy: y1(4) - y2(4)')
 figure; subplot(2,1,1); plot(x, y(20,:)'); hold on; plot(x, y(21,:)'); 
 title(' node-4: output: y1, y2 '); legend('y1(4)', 'y2(4)'); hold off
 subplot(2,1,2); plot(x, y(20,:)'); hold on; plot(x, y(21,:)');  
 xlim([3.9 5.4]); grid on; xlabel('t (s) ')
 % Node 5 outputs:
 figure; plot(y(26,1:end) - y(27,1:end));
 title('node-5 output, y1(5) - y2(5)')
 figure; subplot(2,1,1); plot(x,y(26,:)'); hold on; plot(x,y(27,:)'); 
 title('node-5: output: y1, y2 '); legend('y1(5)', 'y2(5)'); hold off
 subplot(2,1,2); plot(x,y(26,:)'); hold on; plot(x,y(27,:)'); 
 axis([3.5 5 20 130]);  grid on; xlabel('t (s) ') 
 % Node 6 outputs:
 tstart= 40001; % debug
 figure; plot(y(32,1:end) - y(33,1:end));
  title(' node-6 output, y1(6) - y2(6)')
 figure; plot(y(32,:)'); hold on; plot(y(33,:)'); 
  title(' node-6: output: y1, y2 '); legend('y1(6)', 'y2(6)');
 figure; subplot(2,1,1); plot(x,y(32,:)'); hold on; plot(x,y(33,:)'); 
  title(' node-6: output: y1, y2 '); legend('y1(6)', 'y2(6)'); hold off
  subplot(2,1,2); plot(x,y(32,:)'); hold on; plot(x,y(33,:)');
  plot(x, Pd/10) % stim
  xlim([3.9 4.4]); grid on; xlabel('t (s) '); legend('y-e (6)', 'y-i (6)', 'Stim') 
% 
% all 6 nodes together:
tstart= 30001; % debug
 figure; plot(x((tstart-800):end), y(2,(tstart-800):end) - y(3,(tstart-800):end), 'k' );  hold on % node #1, A-10 
 xlim([4.0 5.0]); xlabel('t (s)'); xlabel('t (ms)');  %pause
 plot(x((tstart-800):end) ,(y(8,(tstart-800):end) - y(9,(tstart-800):end) ));   % #2, A-32V
 %pause  % - to closely examine waveforms
 plot(x((tstart-800):end) ,y(14,(tstart-800):end) - y(15,(tstart-800):end)); % #3, A-32
 %pause
 plot(x((tstart-800):end) ,y(20,(tstart-800):end) - y(21,(tstart-800):end), 'b'); % #4, A9
 %pause
 plot(x((tstart-800):end) ,y(26,(tstart-800):end) -   y(27,(tstart-800):end)); % #5, A-46D
 %pause  
 plot(x((tstart-800):end) ,y(32,(tstart-800):end) - y(33,(tstart-800):end)); % #6, A-11
 plot(x((tstart-800):end) ,Vyt((tstart-800):end), 'r--') % stimulii: wave
  legend('delta-y(1)', 'delta-y(2)', 'delta-y(3)', 'delta-y(4)', ...
       'delta-y(5)','delta-y(6)', 'stimulus-2' , 'wave' ); %text(3000, 0.3, 'wt(1-2) = 100')
  % axis([0 5000 -20 20]);  %axis([1 2 -25 25]); 
%
% mean and range of detla-y, overall: 6 nodes together:  
tstart= 10001; tend = 40000-100;% after transients & before stim
 d_yMean=zeros(nn,1); d_yAmpl=zeros(nn,1); 
 d_yMean(1)=mean(y(2,(tstart+2000):tend) - y(3,(tstart+2000):tend));   % node #1,  
 d_yAmpl(1)=max(y(2,(tstart+2000):tend) - y(3,(tstart+2000):tend)) ... % sample before/after stim settles
     - min(y(2,(tstart+2000):tend) - y(3,(tstart+2000):tend));
 d_yMean(2)=mean(y(8,(tstart+2000):tend) - y(9,(tstart+2000):tend));   % #2,
 d_yAmpl(2)=max(y(8,(tstart+2000):tend) - y(9,(tstart+2000):tend)) ...
     - min (y(8,(tstart+2000):tend) - y(9,(tstart+2000):tend));
 d_yMean(3)=mean(y(14,(tstart+2000):tend) - y(15,(tstart+2000):tend)); % #3 
 d_yAmpl(3)=max(y(14,(tstart+2000):tend) - y(15,(tstart+2000):tend)) ...
      -min(y(14,(tstart+2000):tend) - y(15,(tstart+2000):tend));
 d_yMean(4)=mean(y(20,(tstart+2000):tend) - y(21,(tstart+2000):tend)); % #4, 
 d_yAmpl(4)=max(y(20,(tstart+2000):tend) - y(21,(tstart+2000):tend)) ...
      - min(y(20,(tstart+2000):tend) - y(21,(tstart+2000):tend));
 d_yMean(5)=mean(y(26,(tstart+2000):tend) - y(27,(tstart+2000):tend)); % #5, A-46D
 d_yAmpl(5)=max(y(26,(tstart+2000):tend) - y(27,(tstart+2000):tend)) ...
      - min(y(26,(tstart+2000):tend) - y(27,(tstart+2000):tend));
 d_yMean(6)=mean(y(32,(tstart+2000):tend) - y(33,(tstart+2000):tend)); % #6, A-11
 d_yAmpl(6)=max(y(32,(tstart+2000):tend) - y(33,(tstart+2000):tend)) ...
      -min(y(32,(tstart+2000):tend) - y(33,(tstart+2000):tend));   
 disp(' dy(i): mean, amplitude (before) \n'); [d_yMean d_yAmpl] % test/ debug

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
   d_yMean(5)=mean(y(26,(tstart+2000):end) - y(27,(tstart+2000):end)); % #5, A-46D
 d_yAmpl(5)=max(y(26,(tstart+2000):end) - y(27,(tstart+2000):end)) ...
      - min(y(26,(tstart+2000):end) - y(27,(tstart+2000):end));
 d_yMean(6)=mean(y(32,(tstart+2000):end) - y(33,(tstart+2000):end)); % #6, A-11
 d_yAmpl(6)=max(y(32,(tstart+2000):end) - y(33,(tstart+2000):end)) ...
      -min(y(32,(tstart+2000):end) - y(33,(tstart+2000):end));  
 disp(' dy(i): mean, amplitude (> 9sec) \n'); [d_yMean d_yAmpl] % test/ debug
  figure; stem(d_yMean, '+');  hold on; stem(d_yAmpl, ':diamond');
   % hs(1).Marker='+'; hs(2).Marker='diamond'; 
 xlim([0 7]);  %axis([0 7 -15 60 ]); 
 xlabel('NM #'); ylabel('delta-y12(i)  (mV)'); title('dy(i) after stimulus')
  legend('mean', 'p-p ampl.', 'Location', 'northeast') % place out of the way
   % text(5.5, 23, 'Adj = 1');  % text(5.5, 10, 'alpha; C x vector')  
  %   %figure; plot(rvec, d_yAmpl, '.', 'MarkerSize', 20); hold on
   % plot(rvec, d_yMean, '+', 'MarkerSize', 5);
%   
 disp(' >> dy(2) after 16.5 sec:  ')
 tstart=165001; % [9: end] sec 
 d_yMean(2)=mean(y(8,(tstart+2000):end) - y(9,(tstart+2000):end));   % #2,
 d_yAmpl(2)=max(y(8,(tstart+2000):end) - y(9,(tstart+2000):end)) ...
     - min (y(8,(tstart+2000):end) - y(9,(tstart+2000):end));
 [d_yMean(2) d_yAmpl(2)] 
% 
% 3 key nodes together:
 figure;   plot(y(2, (tstart-800):end)-y(3,(tstart-800):end) ); hold on; % #1 
 plot(y(8,(tstart-800):end) - y(9,(tstart-800):end), 'r') % ', 'Color', [0 0.3 0]);  % #2,
 %plot(y(32,(tstart-800):end) - y(33,(tstart-800):end),'g'); % #6, A-11 : Input 2nd
  plot(y(20,(tstart-800):end) - y(21,(tstart-800):end),'b'); % #4,
  plot(Pd((tstart-800):end), 'k');  % stimulii, on same scale
  title('6 NM cluster, JR model:  delta-y ')
  legend('delta-y(1)', 'delta-y(2)', 'delta-y(4)', 'stimulus'  );
%   
% LFP LFP:  net lfp output: sum: detla-y(1) + detla-y(2), for the 4 NM
 lfp= (y(2,:) - y(3,:) + y(8,:) - y(9,:) +y(14,:) - y(15,:) ...
     + y(20,:) - y(21,:) + y(26,:) - y(27,:) +y(32,:) - y(33,:) )/6; % av of 6 nodes
% figure; plot(lfp); title('JR model: net lfp output ')
  figure; plot(x, lfp); title('6-NM cluster, JR model: net lfp output '); xlabel('t (sec) ')
  ylabel('LFT (t)  (mV)'); % hold on; plot(20,-1, 'b^'); plot(20,-1.1, 'b^'); % wave starts
   % text(15, -3.5, 'wave ampl = 10 mV') % plot(4,-4, 'k^'); % Pulse starts
    % text(2.0, -120.0, '6-cluster: ran noise input ') 
    % text(4000, -2.0, 'stimulate #2, 6 ')
    % axis([500 5000 -0.1 0.1])  % axis([500 5000 -2.0 0])
   % text(2.5, -4, '16 Hz wave + 100 Hz pulses to A32V')
  %lfp_stats_2to5= [mean(lfp25(tstart:end))  (max(lfp25(tstart:end)) - min(lfp25(tstart:end)) )]
  clear lfp25  lfp_stats_2to5
  
% Grouped summary: resize figure [5 by 1 plots]
  %figwidth = 1024; figheight = 896; figposition = [100, 100, figwidth, figheight]; % large
  figwidth = 560; figheight = 704; figposition = [500, 100, figwidth, figheight]; % tall
 figure('position',figposition, 'units','pixels');  %figure; % default
 subplot(5,1,1);
 plot(y(2,:) - y(3,:),'b'); % #1, A-10 : Output
 tlim=length(y(2,:) - y(3,:))-1; axis([0,tlim, -Inf, Inf]) % keep vetical as is
  title('6-NM cluster, NM#1, delta-y ')
 subplot(5,1,2);
 plot(y(8,:) - y(9,:), 'r') % ', 'Color', [0 0.3 0]);  % #2,
  title('     NM#2, delta-y '); axis([0,tlim, -Inf, Inf])
 subplot(5,1,3);
      %plot(y(20,:) - y(21,:) ); %  #4,
  plot(y(32,1:end) - y(33,1:end));  % #6
  title('     NM#6, delta-y '); axis([0,tlim, -Inf, Inf])
  subplot(5,1,4);  plot(lfp,'k'); title('     lfp output ');  %  net output
  xlim([0,tlim]); % xlim([0 3.0e5]); %axis([0,tlim, -Inf, Inf])
 subplot(5,1,5); 
  plot((Pd)*IpOn(1) +pran+Vyt, 'g--'); hold on; % theta
  plot((Pd)*IpOn(6) +pran+Vyt, 'r--') % stimulii to gamma, on same scale
  plot(pran, 'Color', [0.9 0.9 0.9]); 
  xlim([0,tlim]); ylim([-10 (max(Pd)+10)])
  legend( 'stimulus-1', 'stimulus-6', 'ran'  ); title('     Inputs: stimulus-> NM ');  
   xlabel('t steps (0.1 msec) '); %axis([0,tlim, -Inf, Inf])
    % text(5500, 50, 'Stim 300 to both #2, #6')
    
%   & LFP:
 lfp_sample=lfp(40001:end); % during stimulus [0.5s sample]
   %lfp_sample=lfp(1100:2900); % debug, nb. wider sample for half t-step
 lfp_stats= [ mean(lfp_sample), (max(lfp_sample)-min(lfp_sample))]
 % figure; plot(lfp_sample - mean(lfp_sample)); title('JR, 6xNM: sum (lfp -mean), during stim. ')
% 
  %if xf >15 % for longer sim (eg 20 sec)
  % lfp_sample=lfp(180000:190000);   % V. much later
  % lfp_stats_at18secLaterStill= [ mean(lfp_sample), (max(lfp_sample)-min(lfp_sample))]
  %end
%  
lfp0= (y(1,:) + y(7,:) +y(13,:)  + y(19,:) + y(25,:)  +y(31,:) )/6; % av of 6 nodes
 figure; plot(x, lfp0); title('6-NM cluster, JR model: net lfp[y0] output '); xlabel('t (sec) ')
  ylabel('LFT (t)  (mV)'); ylim([0 0.2]);
lfp0_sample=lfp0(40001:end); % during stimulus [0.5s sample]
% lfp0_stats= [ mean(lfp0_sample), (max(lfp0_sample)-min(lfp0_sample))]
 
 % Node 2 outputs [this has the primary inputs]:
 % figure;  plot(x, y(7,(1):end) ); 
 % title('6NM Node 2, A32V;  y-0 ')
 % xlabel('t (sec) ');  ylabel('y-0[2](t) (mV) ')

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
 d_y0Mean(5)=mean(y(25,(tstart+2000):end) ); % #4, 
 d_y0Ampl(5)=max(y(25,(tstart+2000):end) ) - min(y(25,(tstart+2000):end) );
 d_y0Mean(6)=mean(y(31,(tstart+2000):end) ); % #4, 
 d_y0Ampl(6)=max(y(31,(tstart+2000):end) ) - min(y(31,(tstart+2000):end) );
 %disp(' y0(i): mean, amplitude (after stim) \n'); [d_y0Mean d_y0Ampl] % 
 
 clear jlist i j in jj j11 j16 k*  In12 deltay lfp_*  d_yM* dy6* lfp0* % tstart pran yic 
% > > > > >

%% 5.0a   1 x replot dy2
figure; title('6NM Node #2, A32V;  delta-y '); hold on
plot(x, (y(8,(1):end) - y(9,(1):end)) ); 
xlabel('t (sec) ');  ylabel('dy2(t) = y-e - y-i  (mV) ')
 %text(3.5, 50, 'On'); text(4, 48, '|'); text(4, 46, '|')
 % text(15, 50, 'gamma stim > #2 & 6'); text(15.5, 47, 'full wave')
% text(14,38,'|'); text(14,40,'|'); text(13.8, 42, 'Off') % label off
%% extra overlay
  figure; title('6NM Node #2, A32V;  delta-y '); hold on
  plot(x(150001:end), (y(8,150001:end) - y(9,150001:end) ), 'Color', [0.3 0.3 0.3] ); % from 15s
  xlabel('t (sec) ');  ylabel('dy2(t) = y-e - y-i  (mV) ')
%% 5.0b Labels on, off
text(4,47,'|'); text(4,48,'|'); text(3.8, 50, 'On') % label on
  %text(10,47,'|'); text(10,48,'|'); text(9.8, 50, 'Off') % label off
  text(10,44,'|'); text(10,45,'|'); text(9.8, 48, 'Reset pulse') % 
% text(14, 47, '6NM + beta stim -> #2,6'); % explain stim
% text(15, 44, 'Half wave'); %text(15, 37, 'full wave, biased');
 
 %% 5.0b : 2 multi plots of dy-2
 % Node 2 outputs [this has the primary inputs]:
  figure;   subplot(2,1,1); %title('6NM Node #2, A32V;  delta-y ')
  ylabel('dy2(t) = y-e - y-i  (mV) ')
  subplot(2,1,2); plot(x, (y(8,(1):end) - y(9,(1):end)) ); 
  xlabel('t (sec) ');  ylabel('dy2(t) = y-e - y-i  (mV) ')
   % text(3.5, 50, 'On'); text(4, 48, '|'); text(4, 46, '|')
   %  text(11, 40, 'stim 100Hz * gamma, full wave')
    %text(19, 57, 'd'); text(13, 50, 'gamma')
 
%% 5.0c : 4 multi plots of dy-2
% Node 2 outputs [this has the primary inputs]:
  figure;   subplot(2,2,1); title('6NM Node #2, A32V;  delta-y ')
  ylabel('dy2(t) = y-e - y-i  (mV) ')
  subplot(2,2,4); %  gamma
  %subplot(2,2,2); %  beta
  plot(x, (y(8,(1):end) - y(9,(1):end)) ); 
  xlabel('t (sec)');  ylabel('dy2(t) = y-e - y-i  (mV) ')
  text(3.5, 50, 'On'); text(4, 48, '|'); text(4, 46, '|')
 %  text(11, 40, 'stim 100Hz * gamma, full wave')
    text(13, 50, 'gamma');  % text(19, 57, 'd');
 
 subplot(2,2,1); text(13, 50, 'theta'); hold on
  plot(x, (y(8,(1):end) - y(9,(1):end)) ); ylabel('dy2(t) = y-e - y-i  (mV) ')
 
  subplot(2,2,2); text(13, 45, 'alpha'); hold on
  plot(x, (y(8,(1):end) - y(9,(1):end)) ); ylabel('dy2(t) = y-e - y-i  (mV) ');
  
 subplot(2,2,3); text(13, 53, 'beta'); hold on
  plot(x, (y(8,(1):end) - y(9,(1):end)) ); ylabel('dy2(t) = y-e - y-i  (mV) '); xlabel('t (sec)');

% publication quality print
% hh =gcf
% print(hh, 'Fig3.tiff', '-dtiff' , '-r300')

%% 5.0.1  mean and range of detla-y, After stim: 6 NM  10:20sec 
% check dy3 before stim
 %dy3LT4 = max( y(14,1:35000) - y(15,1:35000) ) - min( y(14,1:35000) - y(15,1:35000) ) 
 %figure; plot ( y(14,1:45000) - y(15,1:45000) ); title(' NM # 3 - check state')
 % well after stim:
tstart = 110000;
fprintf('\n > dy(i) stats > 11 sec \n')
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
 axis([0 7 -15 80 ]); xlabel('NM #'); ylabel('delta-y12(i)  (mV)'); title('dy(i) 10 : 20 sec')
  legend('mean', 'p-p ampl.', 'Location', 'northeast') % place out of the way
   % text(5.5, 23, 'Adj = 1');  % text(5.5, 10, 'alpha; C x vector')
 
 dy_mean_ampl = [d_yMean d_yAmpl] % test/ debug
 lfp_sample=lfp(tstart:end); % during stimulus [0.5s sample]
   %lfp_sample=lfp(1100:2900); % debug, nb. wider sample for half t-step
 lfp_stats_10s= [ mean(lfp_sample), (max(lfp_sample)-min(lfp_sample))]
 clear dy3LT4 dy_mean_ampl d_y* lfp_st*
 
%% 5.0.2 compare dy12 vs y0
fprintf('\n > plot yo vs y-e, y-i: \n');
figure; plot(x, (y(1,:)*50.0 +20), 'b' ); % its small - check offsets
hold on; plot(x, y(2,:), 'g', 'LineWidth', 1.2);  
plot(x, y(3,:), 'r', 'LineWidth', 1.0 ); % for y-e & y-i
 % hold on; plot(x, (y(2,:) -y(3,:)) ) % for lfp
legend('yo *50', 'y1', 'y2'); grid on
 % legend('yo *50', 'y1 -y2')
xlim([3.5 5]);  % xlim([0.2 0.4]); grid on
 % ylim([-25.0 30])
 title('4NM, dbl tri; #1, beta; yo lag varies')
%
figure; plot(x, (y(7,:)*50.0 +40), 'b' ); % y0, its small - check offsets
hold on; plot(x, y(8,:), 'g', 'LineWidth', 1.2);  
plot(x, y(9,:), 'r', 'LineWidth', 1.0 ); % for y-e & y-i
 % hold on; plot(x, (y(2,:) -y(3,:)) ) % for lfp
legend('yo *50', 'y1', 'y2'); grid on;
 % legend('yo *50', 'y1 -y2')
xlim([3.5 5]);  % xlim([0.2 1.2]); grid on; 
 % ylim([-25.0 30])
 title('4NM, dbl tri; #2, theta; yo & ye, yi')  
 %
figure; plot(x, (y(13,:)*50.0 +40), 'b' ); % its small - check offsets
hold on; plot(x, y(14,:), 'g', 'LineWidth', 1.2);  
plot(x, y(15,:), 'r', 'LineWidth', 1.0 ); % for y-e & y-i
 % hold on; plot(x, (y(2,:) -y(3,:)) ) % for lfp
legend('yo *50', 'y1', 'y2'); grid on;
 % legend('yo *50', 'y1 -y2') 
xlim([3.5 5]);  % xlim([0.2 1.2]); grid on; 
 % ylim([-25.0 30])
 title('4NM, dbl tri; #3, alpha; yo & ye, yi') 
%
figure; plot(x, (y(19,:)*50.0 +20), 'b' ); % its small - check offsets
hold on; plot(x, y(20,:), 'g', 'LineWidth', 1.3);  
plot(x, y(21,:), 'r', 'LineWidth', 0.9 ); % for y-e & y-i
 % hold on; plot(x, (y(2,:) -y(3,:)) ) % for lfp
legend('yo *50', 'y1', 'y2'); grid on;
 % legend('yo *50', 'y1 -y2')
 %xlim([0.2 0.4]); 
xlim([3.5 5]);  %xlim([0.2 1.2]); % span transition at 4s
 % ylim([-25.0 30])
 title('.. NM, #4, xx; yo & ye, yi')
%
figure; plot(x, (y(25,:)*50.0 +20), 'b' ); % its small - check offsets
hold on; plot(x, y(26,:), 'g', 'LineWidth', 1.3);  
plot(x, y(27,:), 'r', 'LineWidth', 0.9 ); % for y-e & y-i
 % hold on; plot(x, (y(2,:) -y(3,:)) ) % for lfp
legend('yo *50', 'y1', 'y2'); grid on;
 % legend('yo *50', 'y1 -y2')
 %xlim([0.2 0.4]); 
xlim([3.5 5]);  %xlim([0.2 1.2]); % span transition at 4s
 % ylim([-25.0 30])
 title('.. NM, #5, xx; yo & ye, yi')

 %% 5.0.3 Transitions: dy-2;  ln(dy) plot
% Node 2 outputs [this os the primary input Hub]:
  figure;  plot(x, (y(8,(1):end) - y(9,(1):end)) ); title('6NM Node-2, delta-y p-p amplitude')
  ylabel('dy-2 = y-e(2) - y-i(2) (mV) '); xlabel('t (sec) ');
  figure;  semilogy(x, (y(8,1:end) - y(9,1:end)) );
  ylim([1 100]); ylabel('ln[ dy-2 = y-e(2) - y-i(2)] (mV) '); xlabel('t (sec) ');
  title('6NM Node-2, log delta-y p-p amplitude')
   % text(12,50, '6NM + p=100 * 2alpha -> NM #2, 6')
   % overall 
  d_yMean(2)=mean(y(8,(tstart+2000):end) - y(9,(tstart+2000):end));   % #2,
  d_yAmpl(2)=max(y(8,(tstart+2000):end) - y(9,(tstart+2000):end)) ...
     - min (y(8,(tstart+2000):end) - y(9,(tstart+2000):end));
 
%% 5.1 revised LFP:  compare dy's, cf 4 NM
  % special case: use 4NM code for only 
  fprintf('\n  lfp : sample 4 nodes only  ')
 lfp= (y(2,:) - y(3,:)  ... % omit NM#2 & 3
     + y(20,:) - y(21,:) + y(26,:) - y(27,:) +y(32,:) - y(33,:) )/4; % av of 4 nodes: 1,4,5,6
 lfp4_sample=lfp(40001:end); % during stimulus [0.5s sample]
   %lfp_sample=lfp(1100:2900); % debug, nb. wider sample for half t-step
 lfp4_stats= [ mean(lfp4_sample), (max(lfp4_sample)-min(lfp4_sample))]
 % check dy6: cf dy1
 figure; subplot(3,1,1); plot(x, y(2,:) - y(3,:)); title('6NM: dy-1');  xlim([5 6]);
  %text(5.7, -8,'v-local = 10 m/s')
 subplot(3,1,2); plot(x, y(8,:) - y(9,:)); title('6NM: dy-2; [-- dy-1]'); 
  xlim([5 6]); hold on; grid on
  plot(x, y(2,:) - y(3,:), 'k--');
 subplot(3,1,3); plot(x, y(14,:) - y(15,:)); title('6NM: dy-3; [-- dy-1]'); hold on;
  xlabel('t (sec)');  xlim([5 6]);  grid on
  plot(x, y(2,:) - y(3,:), 'k--');
 clear lfp4_*

%% 4.1 Align waveforms, by phase: 1st Pulse
% a.0) at NM#1 
dy1_sample=y(2, 40000:120000) - y(3, 40000:120000); % NM#1, delta-y 
figure; plot(dy1_sample); hold on; title('dy(1)'); xlim([10000 20000]) %axis([0.5 2.5 -50 50])
  plot(Vyt(40000:120000), 'g'); % same t window
dym=mean(dy1_sample); plot([0.5 2.5], [dym dym], 'b--'); xlabel('t (0.1ms)')
%
title('JR node #1, 1t pulse: output delta-y: y1(2) - y2(2)') % e - i 
plot(Pd(40000:120000)/10, 'g') % 1st stim [eg beta] % reduce scale [of pulse rate]
plot(Pd2(40000:120000)/10, 'k') % 2nd stim [eg gamma]
plot(Vyt(40000:120000)+dym, 'Color', [1 0 1]) % stim V pulses, wave
legend('delta-y(1)', 'Pulse1', 'Pulse2','wave'); 

%% 4.1 a) Align waveforms,at NM#2 [the input] ~ 4:8 sec
dy2_sample=y(8, 40000:120000) - y(9, 40000:120000); % NM#2, delta-y : ~4:8s
%dy2_sample=y(8, 80000:120000) - y(9, 80000:120000); %
figure; plot(dy2_sample,'LineWidth', 1.0); 
  hold on; title('dy(2)'); %xlim([10000 20000])  %xlim([0 2000]) %axis([0.5 2.5 -50 50])
dym=mean(dy2_sample); plot([0.5 2.5], [dym dym], 'g--'); xlabel('t (0.1ms)')
title('JR node #2, 1t pulse: output delta-y: y1(2) - y2(2)') % e - i 
plot(Pd(40000:120000)/10, 'b') % stim % reduce scale [of pulse rate]
plot(Pd2(40000:120000)/10, 'r') % stim % reduce scale [of pulse rate]
plot(Vyt(40000:120000)+dym, 'Color', 'b') % stim V pulses, wave
plot(Vyt2(40000:120000)+dym, 'Color', 'r') % s
legend('delta-y(2)', 'Pulse1', 'Pulse2', 'wave1', 'wave2'); 
xlim([50000 70000]) % to capture gamma burst at 10 = 4+6 s

%% 4.1 b) Align waveforms,at NM#2 [the input] ~ around 10 sec
%dy2_sample=y(8, 40000:80000) - y(9, 40000:80000); % NM#2, delta-y : ~4:8s
istart= 70000;
dy2_sample=y(8, istart:120000) - y(9, istart:120000); %
figure; plot(dy2_sample,'LineWidth', 1.0); 
  hold on; title('dy(2)'); %xlim([10000 20000])  %xlim([0 2000]) %axis([0.5 2.5 -50 50])
dym=mean(dy2_sample); plot([0.5 2.5], [dym dym], 'g--'); xlabel('t (0.1ms)')
title('JR node #2, 1t pulse: output delta-y: y1(2) - y2(2)') % e - i 
plot(Pd(istart:120000)/10, 'b') % stim % reduce scale [of pulse rate]
plot(Pd2(istart:120000)/10, 'r') % stim % reduce scale [of pulse rate]
plot(Vyt(istart:120000)+dym, 'Color', 'b') % stim V pulses, wave
plot(Vyt2(istart:120000)+dym, 'Color', 'r') % s
legend('delta-y(2)', 'Pulse1', 'Pulse2', 'wave1', 'wave2'); 
% xlim([25000 30000]) % to capture gamma burst
 clear istart 

%% a.1)
dy6_sample=y(32, 40000:60000) - y(33, 40000:60000); % NM#6, delta-y 
figure; plot(dy6_sample); hold on; xlim([0 2000]) %axis([0.5 2.5 -50 50])
dym=mean(dy6_sample); plot([0.5 2.5], [dym dym], 'g--'); xlabel('t (0.1ms)')
title('JR node #6: output delta-y: y1(6) - y2(6)') % e - i 
plot(Pd(40000:60000)/10, 'k') % stim % reduce scale [of pulse rate]
 plot(Vyt(40000:60000)+dym, 'Color', [1 0 1]) % stim V pulses, wave
legend('delta-y(6)', 'Pulse', 'wave')

% b) lfp [the output]
%lfp_sample=lfp(1:1200); % all NM#2, nett lfp
lfp_sample=lfp(40000:60000);
figure; plot( lfp_sample); hold on; xlim([0 2000]) %axis([600 1200 -5 45])
dym=mean(lfp_sample); plot([3.500 6.0], [dym dym], 'g--')
title('JR 6 NM: align lfp w stim') 
plot(Pd(40000:60000)/10, 'k') % stim % reduce scale
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

% d) align Pulse with y-e,i(5)
dym=min(y(26, 40000:60000));
figure; subplot(2,1,1); plot(y(26, 40000:60000)); hold on
plot(dym+Pd(40000:60000)/10, 'k'); legend('y-e(5)', 'Stim');
xlim([0 4000]); grid on; title('Align Pulse (@ 4s) with y-e and y-i, NM # 5 ')
dym=min(y(27, 40000:60000));
subplot(2,1,2); plot(y(27, 40000:60000)); hold on; grid on
plot(dym+Pd(40000:60000)/10, 'k'); legend('y-i(5)', 'Stim')
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
 plot(Vyt(190000:210000)+dym, 'Color', [1 0 1]) % stim V pulses, wave
legend('delta-y(2)', 'Pulse', 'wave'); 

% a.1)
dy6_sample=y(32, 190000:210000) - y(33, 190000:210000); % NM#6, delta-y 
figure; plot(dy6_sample); hold on; xlim([9000 12000]) %axis([0.5 2.5 -50 50])
dym=mean(dy6_sample); plot([19.0 21.0], [dym dym], 'g--'); xlabel('t (0.1ms)')
title('JR node #6, 2nd pulse: output delta-y: y1(6) - y2(6)') % e - i 
plot(Pd2(190000:210000)/10, 'k') % stim % reduce scale [of pulse rate]
 plot(Vyt(190000:210000)+dym, 'Color', [1 0 1]) % stim V pulses, wave
legend('delta-y(6)', 'Pulse', 'wave')

% b) lfp [the output]
%lfp_sample=lfp(1:1200); % all NM#2, nett lfp
lfp_sample=lfp(190000:210000);
figure; plot( lfp_sample); hold on; xlim([9000 12000]) %axis([600 1200 -5 45])
dym=mean(lfp_sample); plot([3.500 6.0], [dym dym], 'g--')
title('JR 6 NM LFP, 2nd pulse: align w stim') 
plot(Pd2(190000:210000)/10, 'k') % stim % reduce scale 
plot(Vyt(190000:210000)+dym, 'Color', [1 0 1]) % wave; align axes
legend('lfp', 'Pulse', 'wave'); %axis([0.5 2.5 -10 30]); 
xlabel('t steps (0.1ms)')

 clear dy2_* dy6_ dym 

%% 4.1.2 ISI.  {was 3.6}  ISI for y-e(2);  dy(2) / y-e(1);  dy(1)
% fprintf('\n > ISI for y-i(1):  '); 
 %fprintf('\n > ISI for y-e(1): from 3.0 sec '); 
fprintf('\n > ISI for dy(5):  from 3.0 sec ');
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
%signal= y(2, 30000:end)-y(3, 30000:end); disp('  dy(1), gamma: \n'); % for NM#1
%signal= y(3, 30000:end); disp('  y-i(1), beta: \n'); % for NM#1, inhib popn
 % signal= y(2,:)-y(3, :);  % <:: all dy(1), for NM#1
    %signal= y(2, 1000:end)-y(3, 1000:end);  % dy(1), for 1/2NM
    %signal= y(8,30000:100000)-y(9,30000:100000);  % short
% signal= y(8,30000:end)-y(9,30000:end);  disp('  dy(2) theta: \n'); %for NM#2
%signal =  (y(20,30000:end) -y(21,30000:end)); t=x(40001:end); disp('  dy(4, gamma): \n');% NM-4
signal =   y(26,30000:end) - y(27,30000:end);  t=x(40001:end); disp('  dy(5, theta): \n');% NM-5
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
 figure; hist(ISI); title('6 NM, ISI distribution'); xlabel(' ISI (sec) ')
 figure; hist(1./ISI); title('6 NM, equiv freq. distribution'); xlabel(' f (Hz) ')
 
% 3 plots together:
 figwidth = 560; figheight = 704; figposition = [500, 100, figwidth, figheight]; % tall
figure('position',figposition, 'units','pixels');  % larger figure; 
subplot(3,1,1); stem(ISI); % xlim([0 50]); %ylim([0.1 0.2]);  % for 6NM, alpha
 ylabel('ISI (sec)'); title('6 MN cluster; node #1: ISI times (sec)')
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
figwidth = 1300; figheight = 450; figposition = [100, 380, figwidth, figheight]; % wide
figure('position',figposition, 'units','pixels');  % larger figure; 
subplot(3,1,1); hold on; %stem(ISI);
 ylabel('ISI (sec)'); title('6 MN cluster; node #5:  ISI times (sec)')
 t0 = PeakLocs(1)*h; tt=t0; % (sec)
for i =1:length(PeakTimes)-1
   line([tt tt], [0.0 ISI(i) ], 'LineWidth', 1.0 );
   plot(tt, ISI(i), 'ob')
   tt=tt+ISI(i); % accumulate ISI's to get this time
end
 text(1.0, 0.18, 'v', 'Color', 'r')  % mark stimulus time, ie 4s (nb. t -3s)
 %text(6.0, 0.5, 'v', 'Color', 'r') 
 subplot(3,1,2); hold on; title('pulse train'); % text(2.5, 0.25, '6 NM, 0.1mV U ran noise')
 t0 = PeakLocs(1)*h; tt=t0; % (sec)
for i =1:length(PeakTimes)-1
   line([tt tt], [0.0 0.2 ], 'LineWidth', 1.0, 'Color','k' );
   tt=tt+ISI(i); % accumulate ISI's to get this time
end
grid on; xlabel(' t (sec +3.0)', 'FontWeight', 'bold'); ylim([0 0.3]); %xlim([0 4]);
  text(1.0, 0.25, 'v', 'Color', 'r') 
  %text(1.0, 0.21, 'v', 'Color', 'r')  % mark stimulus time
    % text(6.0, 0.22, 'v', 'Color', 'r') % marker: wave starts
    % text(6.0, 0.235, 'wave arrives', 'Color', 'r') 
    % subplot(2,1,2); text(1.0, 0.24, '+1 Pulse -> #2]', 'Color', 'r')  % describe stimulu2
subplot(3,1,3); hold on; 
 ylabel('ISI (sec)'); title('detail of ISI times (sec)')
 t0 = PeakLocs(1)*h; tt=t0; % (sec)
for i =1:length(PeakTimes)-1
   line([tt tt], [0.0 ISI(i) ], 'LineWidth', 1.0 );
   plot(tt, ISI(i), 'ob')
   tt=tt+ISI(i); % accumulate ISI's to get this time
end
 xlabel(' t (sec +3.0)', 'FontWeight', 'bold'); 
 xlim([2.5 5.5]); % zoom in - after stim {here at [3]+1  sec}
 text(1.0, 0.28, 'v', 'Color', 'r')  % mark stimulus time, ie 4s (nb. t -3s)
 %text(6.0, 0.5, 'v', 'Color', 'r') 
 
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
grid on; title('6 MN; node #1: y-e(1) vs ISI'); xlabel(' ISI (sec)'); ylabel('y-e(1)(t) (mV) ');

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
 clear tmp PeakLocs PeakTimes Npeaks ISI shortT signal StdDev AvshortT_freq Sample* Range %Code
 
 
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
fprintf('\n  NM#1-e plot done!  ')

  %% 4.2e Phase plot y2'(t) vs y2(t) {ie y-i}  for NM#1 : animated
figure; hold on; title ('y-i (1):  dy2/dt(t) vs y2(t)')
ylabel('dy2/dt(t)' ); xlabel('y21(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ 0 40 -500 1000]); % initiation
for i = 40000:5:60000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(y(3,i), y(6,i), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
  %text(57, -2100, '14 to 16  sec'); %text(27, 12, '0.1 sec')
fprintf('\n  #1-[e-i] done!  ') 

 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
hold off


 %% 4.2e Phase plot y2'(t) vs y2(t) {ie y-i}  for NM # 2 : animated
figure; hold on; title ('y-i (2):  dy2/dt(t) vs y2(t)')
ylabel('dy2/dt(t)' ); xlabel('y2(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ 0 40 -500 1000]); % initiation
for i = 35000:10:55000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(y(9,i), y(12,i), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
 % text(115, 700, '3.5 : 5.5  sec'); %text(27, 12, '0.1 sec')
fprintf('\n  dy-2 done!  ') 

 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
hold off

%% 4.2f Phase plot y1'(t) vs y1(t) {ie y-e}  for NM # 5 : animated
disp('  y-i(5), theta: ');% NM-5
figure; hold on; title ('y-e (5):  dy1 /dt(t) vs y1 (t)')
ylabel('dy-e/dt(t)' ); xlabel('y-e(t)'); 
  %axis([ 30 60 -1000 2000]); % alpha, for phase plot
  %axis([ 0 40 -500 1000]); % initiation
istop=55000;
for i = 35000:10:istop % time window
     % for i = 5000:5:9000 % s/s
     %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); % NM #1
     %plot(y(9,i), y(12,i), 'k.', 'MarkerSize', 10); % NM #2
plot(y(26,i), y(29,i), 'k.', 'MarkerSize', 10); % NM #5 "y2" - ie y-i
timestamp=num2str(i*h, 5); timestamp=horzcat( timestamp, ' sec');
 htext =text(85, -800, timestamp);
 pause(0.001); drawnow
 delete(htext);
pause(0.0005)
drawnow
end
  text(y(27,istop)+1, y(30,istop)+20, timestamp); % at marker
   %text(57, -2100, '14 to 16  sec'); %text(27, 12, '0.1 sec')
   % axis([5 10 -40 40]) % zoom in to tight orbit
fprintf('\n  done!  ') 

 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
hold off
 clear istop 
 
%% 4.2g Phase plot y2'(t) vs y2(t) {ie y-i}  for NM # 5 : animated
disp('  y-i(5), theta: ');% NM-5
figure; hold on; title ('y-i (5):  dy2 /dt(t) vs y2 (t)')
ylabel('dy-i/dt(t)' ); xlabel('y-i(t)'); 
  %axis([ 30 60 -1000 2000]); % alpha, for phase plot
  %axis([ 0 40 -500 1000]); % initiation
istop=55000;
for i = 35000:10:istop % time window
     % for i = 5000:5:9000 % s/s
     %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); % NM #1
     %plot(y(9,i), y(12,i), 'k.', 'MarkerSize', 10); % NM #2
plot(y(27,i), y(30,i), 'k.', 'MarkerSize', 10); % NM #5 "y2" - ie y-i
timestamp=num2str(i*h, 5); timestamp=horzcat( timestamp, ' sec');
 htext =text(85, -800, timestamp);
 pause(0.001); drawnow
 delete(htext);
pause(0.0005)
drawnow
end
  text(y(27,istop)+1, y(30,istop)+20, timestamp); % at marker
   %text(57, -2100, '14 to 16  sec'); %text(27, 12, '0.1 sec')
   % axis([5 10 -40 40]) % zoom in to tight orbit
fprintf('\n  done!  ') 

 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
hold off
 clear istop 

%% 4.3a Phase plot LFP(1): (y1-y2)'(t) vs [y1(t)-y2(t)] {ie y-i}  for NM#1 : animated
figure; hold on; title ('LFP (1):  dy1-2/dt(t) vs y1-2(t)')
ylabel('dy1-2/dt(t)' ); xlabel('y1-2(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ 0 40 -500 1000]); % initiation
for i = 40000:5:45000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot((y(2,i)-y(3,i)), (y(5,i)-y(6,i)), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
% another color, for later
for i = 45000:5:50000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot((y(2,i)-y(3,i)), (y(5,i)-y(6,i)), 'b.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
  %text(57, -2100, '14 to 16  sec'); %text(27, 12, '0.1 sec')
fprintf('\n  #1-LFP[e-i] done!  ') 
 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
hold off

%% 4.3b Phase plot LFP(2): (y1-y2)'(t) vs [y1(t)-y2(t)] {ie y-e - y-i}  for NM#2 : animated
figure; hold on; title ('LFP (2):  dy1-2/dt(t) vs y1-2(t)')
ylabel('dy1-2/dt(t)' ); xlabel('y1-2(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ 0 40 -500 1000]); % initiation
for i = 40000:10:45000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot((y(8,i)-y(9,i)), (y(11,i)-y(12,i)), 'k.', 'MarkerSize', 10);
 timestamp=num2str(i*h, 5); timestamp=horzcat( timestamp, ' sec');
 htext =text(-30, -1500, timestamp);
 pause(0.001); drawnow
 delete(htext);
pause(0.0005)
drawnow
end
fprintf('\n  4:4.5s done! ')
fp1=mean( y(11,40000:45000)-y(12,40000:45000) ); fp2=mean(y(8,40000:45000)-y(9,40000:45000)); %
fp_before = fp2
plot(fp2, fp1, 'k.', 'MarkerSize', 25); % plot y(x coord) & ydot(y coord)
% another color, for later, during stim
for i = 66000:10:75000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot((y(8,i)-y(9,i)), (y(11,i)-y(12,i)), 'b.', 'MarkerSize', 10);
 timestamp=num2str(i*h, 5); timestamp=horzcat( timestamp, ' sec');
 htext =text(-30, -1500, timestamp);
 pause(0.001); drawnow
 delete(htext);
pause(0.0005)
drawnow
end
htext =text(-30, -1500, timestamp); % mark final time
fprintf('\n  6.6:7s done! ')
  text(25, 2700, '4.0:4.5 sec')
  text(30, 2500, '6.6:7.5 sec', 'color', 'b');
fprintf('\n  #2-LFP[e-i] done!  ') 
fp1=mean( y(11,66000:75000)-y(12,66000:75000) ); fp2=mean(y(8,66000:75000)-y(9,66000:75000)); %
fp_atStim = fp2
plot(fp2, fp1, 'b.', 'MarkerSize', 25);% plot y(x coord) & ydot(y coord)
 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
hold off
clear fp* htext timestamp 

 %% 4.3b.1 add later orbit, s/s [ongoing stim
 hold on
 % another color, for later
for i = 110000:10:120000 % time windows
 %for i = 60000:5:65000 % later
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot( (y(8,i)-y(9,i)), (y(11,i)-y(12,i)), 'r.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fp1=mean( y(11,110000:120000)-y(12,110000:120000) ); fp2=mean(y(8,110000:120000)-y(9,110000:120000)); %
fp_after = fp2
plot(fp2, fp1, 'r.', 'MarkerSize', 25); % plot y(x coord) & ydot(y coord)
 text(35, -1800, '11:12 sec', 'color', 'r');
fprintf('\n  extra done [to 12s]!  \n')
hold off
clear fp* 

%% 4.3b.2 add still later orbit, recovery [after Off]
 hold on
 % another color, for later
for i = 140000:10:150000 % time windows
 %for i = 60000:5:65000 % later
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot( (y(8,i)-y(9,i)), (y(11,i)-y(12,i)), 'r.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fp1=mean( y(11,140000:150000)-y(12,140000:150000) ); fp2=mean(y(8,140000:150000)-y(9,140000:150000)); %
fp_after = fp2
plot(fp2, fp1, 'r.', 'MarkerSize', 25); % plot y(x coord) & ydot(y coord)
 text(35, -1800, '14:15 sec', 'color', 'r');
fprintf('\n  extra+ done [to 15s]!  \n')
hold off
clear fp* 

%% 4.3b.3 add still later orbit, s/s  [after recovery phase
 hold on
 % another color, for later
for i = 170000:10:180000 % time windows
 %for i = 60000:5:65000 % later
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot( (y(8,i)-y(9,i)), (y(11,i)-y(12,i)), 'r.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fp1=mean( y(11,170000:180000)-y(12,170000:180000) ); fp2=mean(y(8,170000:180000)-y(9,170000:180000)); %
fp_after = fp2
plot(fp2, fp1, 'r.', 'MarkerSize', 25); % plot y(x coord) & ydot(y coord)
 text(35, -1800, '17:18 sec', 'color', 'r');
fprintf('\n  extra+ done [to 15s]!  \n')
hold off
clear fp* 

%% 4.3c Phase plot LFP(5): (y1-y2)'(t) vs [y1(t)-y2(t)] {ie y-e - y-i}  for NM#2 : animated
figure; hold on; title ('LFP (5):  dy1-2/dt(t) vs y1-2(t)')
ylabel('dy1-2/dt(t)' ); xlabel('y1-2(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ 0 40 -500 1000]); % initiation
for i = 40000:5:45000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot((y(26,i)-y(27,i)), (y(29,i)-y(30,i)), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
% another color, for later
for i = 45000:5:50000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot((y(26,i)-y(27,i)), (y(29,i)-y(30,i)), 'b.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
  %text(57, -2100, '14 to 16  sec'); %text(27, 12, '0.1 sec')
fprintf('\n  #5-LFP[e-i] done!  ') 
 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
hold off

%% 4.3c Phase plot LFP(5): (y1-y2)'(t) vs [y1(t)-y2(t)] {ie y-e - y-i}  for NM#5 : animated
figure; hold on; title ('LFP (#5):  dy1-2/dt(t) vs y1-2(t)')
ylabel('dy1-2/dt(t)' ); xlabel('y1-2(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ 0 40 -500 1000]); % initiation
for i = 40000:5:45000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot((y(26,i)-y(27,i)), (y(29,i)-y(30,i)), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
% another color, for later
for i = 45000:5:50000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot((y(26,i)-y(27,i)), (y(29,i)-y(30,i)), 'b.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
  %text(57, -2100, '14 to 16  sec'); %text(27, 12, '0.1 sec')
  text(-5, 1500, '4.0: 4.5 sec'); text(5, 1300, '4.5: 5.0 sec', 'color', 'b')
fprintf('\n  #5-LFP[e-i] done!  ') 
 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
hold off

%% 4.3c.1 add later orbit for dy(5)
 hold on
 % another color, for later
for i = 110000:5:120000 % time windows
 %for i = 60000:5:65000 % later
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot( (y(26,i)-y(27,i)), (y(29,i)-y(30,i)), 'r.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fp1=mean(y(26, 110000:120000)-y(27, 110000:120000)); % later, after stim
fp2=mean(y(29, 110000:120000)-y(30, 110000:120000)); % just for plotted orbit
plot(fp2, fp1, 'r.', 'MarkerSize', 25);
FixedPoint_ss=[fp1 fp2] % later, in beating phase
 text(0, -1200, '11:12 sec', 'color', 'r');
fprintf('\n  extra done!  \n')
hold off

clear fp* FixedPoint* 

%% 4.3d Phase plot LFP(6): (y1-y2)'(t) vs [y1(t)-y2(t)] {ie y-e - y-i}  for NM#6 : animated
figure; hold on; title ('LFP (6):  dy1-2/dt(t) vs y1-2(t)')
ylabel('dy1-2/dt(t)' ); xlabel('y1-2(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ 0 40 -500 1000]); % initiation
for i = 40000:5:45000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot((y(32,i)-y(33,i)), (y(35,i)-y(36,i)), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
% another color, for later
for i = 45000:5:50000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot((y(32,i)-y(33,i)), (y(35,i)-y(36,i)), 'b.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
  %text(57, -2100, '14 to 16  sec'); %text(27, 12, '0.1 sec')
fprintf('\n  #6-LFP[e-i] done!  \n') 
 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
hold off

%% tmp

fp1=mean(y(2,2201:end)); fp2=mean(y(3,2201:end)); %
plot(fp1, fp2, 'b.', 'MarkerSize', 25);
FixedPoint_ss=[fp1 fp2]
 
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
for i = 35000:10:40000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(lfp(i), lfpdot(i), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fp1=mean(lfpdot(35000:40000)); fp2=mean(lfp(35000:40000)); %
fp_beforeStim = fp2
plot(fp2, fp1, 'k.', 'MarkerSize', 25);
% another color, for later
for i = 40000:10:50000 % time windows
 %for i = 60000:5:65000 % later
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(  lfp(i), lfpdot(i), 'b.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fp1=mean(lfpdot(40000:45000)); fp2=mean(lfp(40000:45000)); %
fp_atStim = fp2
plot(fp2, fp1, 'b.', 'MarkerSize', 25);
  %text(57, -2100, '14 to 16  sec'); %text(27, 12, '0.1 sec')
 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
 text(5, 1200, '3.5: 4.0 sec (before)'); text(5, 1100, '4.0: 4.5 sec (stim)', 'color', 'b')
hold off
FixedPoint_ss=[fp1 fp2]
fprintf('\n  av:LFP [e-i] done!  \n')   
 clear fp* FixedPoint*  %lfpdot 

 %% 4.4.1 add later orbit
 hold on
 % another color, for later
for i = 110000:10:120000 % time windows
 %for i = 60000:5:65000 % later
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot( lfp(i), lfpdot(i), 'r.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fp1=mean(lfpdot(110000:120000)); fp2=mean(lfp(110000:120000)); % overall, after stim
fp_afterStim = fp2
plot(fp2, fp1, 'r.', 'MarkerSize', 25); % plot (x,y) ie (lfp, lfp-dot)
FixedPoint_ss=[fp1 fp2] % later, in beating phase
 text(5, 1000, '11:12 sec (s/s)', 'color', 'r');
fprintf('\n  extra done!  \n')
hold off

clear fp* FixedPoint*  lfpdot

%% 4.5 Phase plot LFP0(all): av[y0)'(t)] vs av[y0(t)] {ie interneurons} for NM#1-6 : animated
fprintf('\n  > LFP0 av. for all NM  \n')
 % LFP0 LFP0:  net lfp0 output:
 
 lfp0= (y(1,:) + y(7,:) +y(13,:)  + y(19,:) + y(25,:)  +y(31,:) )/6; % av of 6 nodes
 lfpdot= (y(4,:) + y(10,:) +y(16,:) + y(22,:) + y(28,:) +y(34,:) )/6;
figure; hold on; title ('LFP0 (all):  dy0/dt(t) vs y0(t)')
ylabel('dy0/dt(t)' ); xlabel('y0(t)'); 
%axis([ 30 60 -1000 2000]); % alpha, for phase plot
%axis([ 0 40 -500 1000]); % initiation
for i = 35000:10:40000 % time windows
 % for i = 5000:5:9000 % s/s
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(lfp0(i), lfpdot(i), 'k.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fp1=mean(lfpdot(35000:40000)); fp2=mean(lfp0(35000:40000)); %
fp_beforeStim = fp2
plot(fp2, fp1, 'k.', 'MarkerSize', 25);
% another color, for later
for i = 40000:10:50000 % time windows
 %for i = 60000:5:65000 % later
  %plot(y(2,i), y(3,i), 'k.', 'MarkerSize', 8); 
plot(  lfp0(i), lfpdot(i), 'b.', 'MarkerSize', 10);
pause(0.0005)
drawnow
end
fp1=mean(lfpdot(40000:45000)); fp2=mean(lfp0(40000:45000)); %
fp_atStim = fp2
plot(fp2, fp1, 'b.', 'MarkerSize', 25);
  %text(57, -2100, '14 to 16  sec'); %text(27, 12, '0.1 sec')
 % text(37, 20, '0.6 sec'); text(27, 12, '0.1 sec')
 text(5, 1200, '3.5: 4.0 sec (before)'); text(5, 1100, '4.0: 4.5 sec (stim)', 'color', 'b')
hold off
FixedPoint_ss=[fp1 fp2]
fprintf('\n  av:LFP0 done!  \n')   
 clear fp* FixedPoint*  %lfpdot 


%% Appx. 0.  w - r 
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

%% Appx 1.0 driving forces

%   >   >   >   >
%% Appx 2 calc, plot S[V(t)] forms, driving forces
 fprintf('\n >>> NM #2, S[V] fn *scaling, 4:5 sec: \n')
in=5  %in=5;  % working on NM #2, 5 : target
%in=6
r=rvec(in); wbar=wvec(in); % for indiv NM 
v0=6;    % (mV) midpoint (switching threshold V) of sigmoid fn
tstart =40001; tend=50000; % tend=tstart+5;  % test %tend=50001; % examine 4:5 sec
onerow =ones(1,(tend -tstart+1)); % need vector form here
StimTot =zeros(1,(tend -tstart+1));  % StimTot =zeros(1,(tend -tstart+1));
 % accumulator, in time vector, for ea NN of #in
StimTotAll = StimTot; % to accumulate tot of all NM
jlist=find(Adj(:,in));  % linked nodes in (col) to node #in
nabors= jlist' % debug
jnn=length(jlist); % #NN of in 
for j=1:jnn % scan NN list of node #in
            %In12=0.0*onerow; %  initialise for ea NM
           StimTot =zeros(1,(tend -tstart+1));  % reset for ea. input NN 
           if isempty(jlist)    
               continue % no NN, In12 stays at 0
           else
           jj=jlist(j); j11=(jj-1)*6+1; j16=(jj-1)*6+6; % for nodes 1, 2, 3 ..
           end % check for NN
           
           for i = tstart:tstart+9000 %tend-1 % scan time steps {eg. start 4 : 5 sec
             In12=0.0; %  initialise for ea NM at each time step
             dyi= y(j11+1,i)-y(j11+2,i); % delta-y of nn. at this time step (i)
            %[jj dyi] % debug
           %argm=r*(wbar*onerow -v0./dyi); In12=0.0*onerow; % erfc form for wt; avoid v=0?
           argm=r*(wbar -v0/dyi); In12=0.0; % nb. scalar at ea  t step
             % nb need vector form here
          if argm >0 & dyi>0.01  % avoid difficulties of erfc(neg, 0)
            %In12= R12*Adj(jj,in)*vm*(1.0*onerow-erfc(argm)/2); % erfc form for wt; avoid v=0
            In12= R12*Adj(jj,in)*vm*(1.0-erfc(argm)/2);
            % [jj In12] % debug
            % sum(Adj(jlist,in))
          end
          %[jj In12] % debug
          %[i dyi argm In12]  % debug      
          % working on node in:  weight by In link, from j to i: jj -> in 
    %fprintf('\n  %2.0f %2.0f %5.3f %4.1f %4.0f %4.2f %4.2f ', in, jj, Adj(jj,in), Ad6A(jj,in), (i-deltai), r, Cfac )        
         StimTot(i-tstart+1) = StimTot(i-tstart+1) + In12; %nb counter starts at tstart (~ 40k)
         %StimTot(i)
           end % loop over time steps(i)
         %nabor = jj % debug
         %StimTot = R12*(A/tauevec(in))*StimTot; % wt sum for driving force in DE * scaling
         StimTot = (A/tauevec(in))*StimTot; % scaling; In12 already incl R12, above
         figure; plot(StimTot); % nb. t course for each NN
           title(strcat('S[V(t)] from NM  # ', num2str(jj), '  to  # ', num2str(in)) )
         StimTotAll = StimTotAll + StimTot;
       end % scan linked NNs(j); gather inputs (in <- jj) from linked NNs
 
       % figure; plot(dyi) % debug
       figure; plot(StimTotAll);  title('S[V(t)] from All NM ')    
       % figure; plot(x(tstart:tend),StimTot); 

  clear argm i j* wbar In12 nabor* onerow tstart tend % StimTot* in

 %% Appx 2.01 other driving forces in DE
 in=2  %in=5;  % working on NM #2, 5 : target
 j11=(in-1)*6+1; j16=(in-1)*6+6; % for nodes 1, 2, : find the 6 DE  variables3
 ye=y(j11+1,:); yedot =y(j11+4,:);
 damping = (2/tauevec(in))*yedot;
 figure; plot(damping);  title('DE damping term - NM #2')
  xlim([5e4 5.1e4])
 restoring = (1/tauevec(in)^2)*ye;
 figure; plot(restoring);  title('DE restoring term - NM #2')
  xlim([5e4 5.1e4]); ylabel(' [mV / s]'); xlabel('t steps (0.1 ms) - 5 s ');
 pulse = (A/tauevec(in))*Pd; % ext, modulated stim
 figure; plot(pulse);  title('ext stim pulse - NM #2')
  xlim([5e4 5.1e4]); ylabel(' [mV / s]'); xlabel('t steps (0.1 ms) - 5 s ');
  clear ye yedot j* pulse
  
%% Appx 2.02 plot S[V(t)]
   % cf WaveTrains.m  # Appx 2.1.2
tstart =40001; tend=42500;  %50001; % examine 4:5 sec
onerow =ones(1,(tend -tstart+1));
jj=1 % NM #
j11=(jj-1)*6+1; j16=(jj-1)*6+6; % for nodes 1, 2, 3,4,6 
dyi= y(j11+1,tstart:tend)-y(j11+2,tstart:tend); % delta-y of nn.
figure; subplot(2,1,1); plot(x(tstart:tend), dyi); grid on
  ylabel('dy(1) (mV)'); title('V(t)')
argm=r*(wbar*onerow -v0./dyi);
 S2=vm*(1.0*onerow-erfc(argm)/2);
  %S2= 1 -erfc(r*(wbar*onerow -v0./vpos))/2;
  %maxS_Wt5 = S2(end)
 subplot(2,1,2); plot(x(tstart:tend), S2, 'r-' ); grid on
 xlabel('t (sec)'); ylabel('stim S (Hz)'); title('S[ V(t) ]') 
 
 clear argm S2 maxS* 

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
 wbar= 5.0; % {max} mean wt/link, for this pop'n
S2= 1 -erfc(r0*(wbar*onerowp -v0./vpos))/2;
 maxS_Wt5 = S2(end)
 figure; hold on; plot(vpos, S2, 'r--' ); %ylim([0 1]);  %xlim([-10, 40]); %ylim([0 1]);
 xlim([0 50])
  wbar= 3.0; % {max} mean wt/link, for this pop'n
S2= 1 -erfc(r0*(wbar*onerowp -v0./vpos))/2;
 maxS_Wt3 = S2(end)
 plot(vpos, S2, 'b--' ); %ylim([0 1]);  %xlim([-10, 40]); %ylim([0 1]);

% unit av wt
 wbar= 1.0; % {max} mean wt/link, for this pop'n
S2= 1 -erfc(r0*(wbar*onerowp -v0./vpos))/2;
 maxS_Wt1 = S2(end)
 plot(vpos, S2, 'g --' );
% lo wt 
wbar=1.5; % {min} mean wt/link
S2= 1 -erfc(r0*(wbar*onerowp -v0./vpos))/2;  %erfc(r0*(v0./vpos -wbar*onerowp))/2;
  maxS_t1p5 = S2(end)
  plot(vpos, S2, 'b-.' );
% cf. sigmoid  
plot(vpos, 1.0./(1+exp(r0*(v0*onerowp - vpos )) ), 'k-'); % basic sigmoid S[v] form
  title('S[v] via synaptic wt distribution; erfc (1/V)')
  legend('wbar = 5.0', 'wbar = 3.0', 'wbar = 1.0', 'wbar = 1.5', 'S[v], sigmoid') 
 text(30, 0.1, 'r = 0.5; Vo = 6 mV');  
  plot([0 50], [0.5 0.5], 'b:') % 50% activation

  clear r0 r1 onerow onerowp v vpos S1 S2 maxS*

%% Appx 3.0 wt Degree (Strength) in 6* cluster
 fprintf('\n >>> Cluster Adj: Wt Degree: \n')
DegIn=sum(Adj,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
DegOut=sum(Adj,2); % i->j; row sum 
DegIn'
DegOut'

%% Appx 3.2a Local, internal Links & wts In/Out nodes, 6x6 cluster
 %  Adj from Marmoset (LNe) data / code from sixClusterWavesV4.m 
 % cf.Appx 2. in orig order: cf. logbook, p58
  fprintf('  use marmoset LNe Adj(2): net In links & wt, 8x8 \n')
clusterlist=[2, 26, 25, 48, 32, 3]; % for 6x6; as plotted in x-y plane, viewed from Ant.
%clusterlist=[2, 26, 25, 48, 32, 3, 45, 47]; % for 8x8
Anew= csvread('AdjMarmoset_Rescale2b.csv');  % based on LNe
Anew = Anew/1.0e3; % scaled as for NM calc
% #Links: un-wt-Deg, less int links ~ net ext links
fprintf('\n :: Local # links: \n')
 kIn_1 =length(find(Anew(clusterlist,2))) 
 kOut_1=length(find(Anew(2, clusterlist)))
 kIn_2=length(find(Anew(clusterlist,26)))  % # - within cluster
 kOut_2=length(find(Anew(26, clusterlist)))
 kIn_3 =length(find(Anew(clusterlist,25)))
 kOut_3=length(find(Anew(25, clusterlist)))
 kIn_4=length(find(Anew(clusterlist,48)))  % within cluster
 kOut_4=length(find(Anew(48, clusterlist)))
 kIn_5=length(find(Anew(clusterlist,32)))  % within cluster
 kOut_5=length(find(Anew(32, clusterlist)))
 kIn_6=length(find(Anew(clusterlist,3)))  %  6x6
 kOut_6=length(find(Anew(3, clusterlist)))
 kIn_7=length(find(Anew(clusterlist,45)))  %  
 kOut_7=length(find(Anew(45, clusterlist)))
 kIn_8=length(find(Anew(clusterlist,47)))  %  8x8 now
 kOut_8=length(find(Anew(47, clusterlist)))
  
% wt-Deg, Local links
% Out Hub
fprintf('\n :: Local{internal} links: wt \n')
 DegIn_1 =sum((Anew(clusterlist,2))) 
 DegOut_1=sum((Anew(2, clusterlist)))
 DegIn_2=sum((Anew(clusterlist,26))) % # wt int links
 DegOut_2=sum((Anew(26, clusterlist)))  % within cluster only
 DegIn_3 =sum((Anew(clusterlist,25)))
 DegOut_3=sum((Anew(25, clusterlist)))
 DegIn_4 =sum((Anew(clusterlist,48)))
 DegOut_4=sum((Anew(48, clusterlist)))
 DegIn_5 =sum((Anew(clusterlist,32)))
 DegOut_5=sum((Anew(32, clusterlist)))
 DegIn_6=sum((Anew(clusterlist,3)))
 DegOut_6=sum((Anew(3, clusterlist)))
 DegIn_7=sum((Anew(clusterlist,45)))
 DegOut_7=sum((Anew(45, clusterlist)))
 DegIn_8=sum((Anew(clusterlist,47)))
 DegOut_8=sum((Anew(47, clusterlist)))
 
 fprintf('\n :: ratios wtIn/kIn{local} \n')
 wtIn_k_1=DegIn_1/kIn_1
 wtOut_k_1=DegOut_1/kOut_1
 RatioAv_1 = (wtIn_k_1 +wtOut_k_1)/2
 
 wtIn_k_2=DegIn_2/kIn_2
 wtOut_k_2=DegOut_2/kOut_2
 RatioAv_2 = (wtIn_k_2 +wtOut_k_2)/2
 
 wtIn_k_3=DegIn_3/kIn_3
 wtOut_k_3=DegOut_3/kOut_3
 RatioAv_3 = (wtIn_k_3 +wtOut_k_3)/2
 
 wtIn_k_4=DegIn_4/kIn_4
 wtOut_k_4=DegOut_4/kOut_1
 RatioAv_4 = (wtIn_k_4 +wtOut_k_4)/2
 
 wtIn_k_5=DegIn_5/kIn_5
 wtOut_k_5=DegOut_1/kOut_5
 RatioAv_5 = (wtIn_k_5 +wtOut_k_5)/2
  
 wtIn_k_8=DegIn_6/kIn_6
 wtOut_k_6=DegOut_6/kOut_6
 RatioAv_6 = (wtIn_k_8 +wtOut_k_6)/2
 
 wtIn_k_7=DegIn_7/kIn_7
 wtOut_k_7=DegOut_7/kOut_7
 RatioAv_7 = (wtIn_k_7 +wtOut_k_7)/2
 
 wtIn_k_8=DegIn_8/kIn_8
 wtOut_k_8=DegOut_6/kOut_8
 RatioAv_8 = (wtIn_k_8 +wtOut_k_8)/2
 
 clear k* Deg* clusterlist wtI* wtO* Ratio* Anew  

%% Appx 3.2b Ext Links & wts In/Out nodes, 6x6 cluster
 %  Adj from Marmoset (LNe) data / code from sixClusterWavesV4.m 
 % cf.Appx 2. in orig order: cf. logbook, p58
  fprintf('  use marmoset LNe Adj(2): net In links & wt, \n')
%clusterlist=[2, 26, 25, 48, 32, 3]; % as plotted in x-y plane, viewed from Ant.

Anew= csvread('AdjMarmoset_Rescale2b.csv');  % based on LNe
nLinks = length(find(Anew(:))) % tot # links: 3474, 26% of tot possible
% #Links: un-wt-Deg, less int links ~ net ext links
fprintf('\n :: ext # links: \n')
 kIn_1 =length(find(Anew(:,2))) -length(find(Anew(clusterlist,2))) 
 kOut_1=length(find(Anew(2,:))) - length(find(Anew(2, clusterlist)))
 kIn_2=length(find(Anew(:,26))) - length(find(Anew(clusterlist,26)))  % tot - within cluster
 kOut_2=length(find(Anew(26,:))) -length(find(Anew(26, clusterlist)))
 kIn_6=length(find(Anew(:,3))) - length(find(Anew(clusterlist,3)))  % tot - within cluster
 kOut_6=length(find(Anew(3,:))) - length(find(Anew(3, clusterlist)))
  
% wt-Deg, ext links
% Out Hub
 DegIn_1 =sum((Anew(:,2))) -sum((Anew(clusterlist,2))) 
 DegOut_1=sum((Anew(2,:))) - sum((Anew(2, clusterlist)))
fprintf('\n :: ext links: wt \n')
 DegIn_2=sum((Anew(:,26))) -sum((Anew(clusterlist,26))) % less int links
 DegIn_6=sum((Anew(:,3))) - sum((Anew(clusterlist,3)))  % within cluster
 DegOut_2=sum((Anew(26,:))) -sum((Anew(26, clusterlist)))
 DegOut_6=sum((Anew(3,:))) - sum((Anew(3, clusterlist)))
 
 fprintf('\n :: ratios wtIn/kIn{Ext} \n') 
 %Other key Inputs
 kIn_5=length(find(Anew(:,32))) - length(find(Anew(clusterlist,32)))  % within cluster
 DegIn_5 =sum((Anew(:,32))) -sum((Anew(clusterlist,32)))
 
 kIn_4=length(find(Anew(:,48))) - length(find(Anew(clusterlist,48)))  % within cluster
 DegIn_4 =sum((Anew(:,48))) -sum((Anew(clusterlist,48)))
 
 kIn_3=length(find(Anew(:,25))) - length(find(Anew(clusterlist,25)))  % within cluster
 DegIn_3 =sum((Anew(:,25))) -sum((Anew(clusterlist,25)))
 
  %clear k* Deg* Anew % clusterlist
  
%% Appx 3.3 w-bar calc
% from Marm link data (LNe v2b) in 6*pFC cluster: cf. Marmoset_wbarCalc.xls
 wink = [0.522, 3.22, 0.485, 0.736, 1.3,  2.346]; % local w-in/k-in
 woutk= [4.684, 0.25, 1.80,  1.11, 0.243, 0.484]; % local w-out/k-out
 wavk = [2.60,  1.74, 1.14,  0.92, 0.77,  1.42]; % an(in , out)

figure; plot(wink,  '-^', 'MarkerSize', 12); hold on
plot(woutk,  '-o', 'MarkerSize', 12);
plot(wavk,  '-ks', 'MarkerSize', 12);
title('Marmoset 6* pFC cluster, wt/#Links'); xlim([0.5 6.6]);
ylabel('wt / #Links'); xlabel('Neural Mass #')
legend('In', 'Out', 'Av')
text(1.15, 4.68,'A-10'); text(2.15, 3.25, 'A32-V'); text(5.8, 2.6,'A-11');
hold off 

%% Appx 3.3.1 w-bar vs Vol
% from Marm link data (LNe v2b) in 6*pFC cluster: cf. Marmoset_wbarCalc.xls
 NodeVol = [23.9, 2.7, 11.4, 7.7, 1.2, 12]; % (mm^3)
 Nn = [19.1, 2.4, 7.6, 5.6, 2.6, 8.6 ]; % #neurons in Anat Area (x1e5)
 wink = [0.522, 3.22, 0.485, 0.736, 1.3,  2.346]; % local w-in/k-in
 woutk= [4.684, 0.25, 1.80,  1.11, 0.243, 0.484]; % local w-out/k-out
 wavk = [2.60,  1.74, 1.14,  0.92, 0.77,  1.42]; % an(in , out)

figure; plot(NodeVol, wink,  '-^', 'MarkerSize', 12); hold on; pause
plot(NodeVol, woutk,  '-o', 'MarkerSize', 12); pause
plot(NodeVol, wavk,  '-ks', 'MarkerSize', 12);
title('Marmoset 6* pFC cluster, wt/#Links'); %xlim([0.5 6.6]);
ylabel('wt / #Links'); xlabel('Node Vol (mm3) #')
legend('In', 'Out', 'Av')
 %text(1.15, 4.68,'A-10'); text(2.15, 3.25, 'A32-V'); text(5.8, 2.6,'A-11');
hold off 


%% Appx 6.  w's vs N-rank
 % calc in Marmoset116NodeInfoV2b.xls, tab 2 (1/24)
 % cf. NM_LIF_calc.m
     %  ID#   N-rank   w-i/k-in  w-ou/k-out Av(wi,o) Max;   nb. (*10^3)
 tmp = [ 1    8.7400    0.5500    4.6840   2.6000  4.6840; ...
         2    1.0000    3.2200    0.2500   1.7400  3.2200; ...
         3    3.7900    0.4850    1.8000   1.1400  1.8000; ...
         4    2.4500    0.7360    1.1100   0.9200  1.1100; ...
         5    1.0500    1.3000    0.2400   0.7700  1.3000; ...
         6    3.9900    2.3460    0.4840   1.4200  2.346]; ...
         
figure; plot(tmp(:, 2), tmp(:, 3), 'ob') % w-in vs N-rank 
hold on;  plot(tmp(:, 2), tmp(:, 4), '^r', 'MarkerSize', 10)
plot(tmp(:, 2), tmp(:, 5), 'sk', 'MarkerSize', 10)
plot(tmp(:, 2), tmp(:, 6), 'vg', 'MarkerSize', 10)
xlim([0.5 9.5]); xlabel('N-rank')
%
figure; subplot(2,2,1); plot(tmp(:, 2), tmp(:, 3), 'ob')
 xlim([0.5 9.5]); xlabel('N-rank'); ylabel('w-in / k-in')
 title('6 NM,  synaptic weight vs node size, N')
subplot(2,2,2); plot(tmp(:, 2), tmp(:, 4), '^r', 'MarkerSize', 10)
 xlim([0.5 9.5]); xlabel('N-rank'); ylabel('w-out / k-out')
subplot(2,2,3); plot(tmp(:, 2), tmp(:, 5), 'sk', 'MarkerSize', 10)
 xlim([0.5 9.5]); xlabel('N-rank'); ylabel('Av(w/k - in, out)')
subplot(2,2,4); plot(tmp(:, 2), tmp(:, 6), 'vg', 'MarkerSize', 10)
 xlim([0.5 9.5]); xlabel('N-rank'); ylabel('Max(w/k - in, out)')

%% Appx 3.3.2 f-i vs Vol
% from Marm link data (LNe v2b) in 6*pFC cluster: cf. Marmoset_wbarCalc.xls
 NodeVol = [23.9, 2.7, 11.4, 7.7, 1.2, 12]; % (mm^3)
 freqi = [32.1 5.4 16.1 9.9 5.9 16.1]; % bands, as assigned in Trail 2: ~ N
 figure; figure; plot(NodeVol, freqi,  '-o', 'MarkerSize', 12); hold on;
 ylabel('freq (Hz)'); xlabel('Node Vol (mm3) ')
 
 %% Appx 3.3.2 f-i vs wt/k
 figure; plot(wink, freqi,  '-^', 'MarkerSize', 12); hold on; pause
 plot(woutk, freqi, '-o', 'MarkerSize', 12); pause
 plot(wavk, freqi, '-ks', 'MarkerSize', 12);
 ylabel('freq (Hz)');
 
wmaxk= [4.7, 3.2, 1.8, 1.1, 1.3, 2.3] % max of i, o, or av
figure; plot(wmaxk, freqi,  '-^', 'MarkerSize', 12); hold on;
ylabel('freq (Hz)'); xlabel('wt-max / #Links');

figure; plot(wink, freqi,  '-^', 'MarkerSize', 12); hold on;
ylabel('freq (Hz)'); xlabel('wt-in / #Links');

 figure; plot(woutk, freqi,  '-o', 'MarkerSize', 12); hold on;
ylabel('freq (Hz)'); xlabel('wt-out / #Links');

figure; plot(wavk, freqi,  'ks', 'MarkerSize', 12); hold on;
ylabel('freq (Hz)'); xlabel('av(wt / #Links)');

wmixk= [2.6, 3.2, 2, 1.5, 2, 3] % mixture of av & max of i,o [8],p66
figure; plot(wmixk, freqi,  '-^', 'MarkerSize', 12); hold on;
ylabel('freq (Hz)'); xlabel('wt-mixture / #Links');


%% Appx 3.3.3 f-i vs Number of neurons
figure; plot(NodeVol, Nn,  '-o', 'MarkerSize', 12); hold on;
xlabel('Node Vol (mm3)'); ylabel('Number of neurons (x1e5)');
 % linear
figure; plot(Nn, freqi,  '-ks', 'MarkerSize', 12); hold on;
ylabel('freq (Hz)'); xlabel('Number of neurons (x1e4)');
 % 1/sqrt(N)
 N2= 1.0./sqrt(Nn);
figure; plot(N2, freqi,  '-ks', 'MarkerSize', 12); hold on;
ylabel('freq (Hz)'); xlabel('1 / sqrt(N)'); xlabel('Node #')

figure; plot(N2); xlim([0.5 6.6]); grid on
ylabel('1 / sqrt(N)'); xlabel('Node #')

freqr = [5.4, 32.1, 9.9, 16.1, 32.1, 9.9]; % reverse, revised order: ~ 1/sqrt(N)
figure; plot(N2, freqr,  '-ks', 'MarkerSize', 12); hold on;
ylabel('freq (Hz)'); xlabel('1 / sqrt(N)'); 

%% Appx 4.0  w-av in sim vs w-av from w(i,o)/k data
wav_calc = [2.6000 1.7400 1.1400 0.9200 0.7700 1.4200]; % calc in 116Node..xl, Table S2.
wav_rescale = wav_calc/0.77; % scal from min value
wav_max = [4.6840 3.2200 1.8000 1.1100 1.3000 2.3460]; % max of i, o: cf xl.
w_trial6 =[3, 5, 3, 1.5, 3.5, 3];  % wb[8], p43
fi = [32.1 5.3 16.0 9.9 5.9 16.1]; % from 6NM sim, Table 1, mss. (1/24)
Acrn6={'A10', 'A32V', 'A32','A9', 'A46D','A11'};

figure; hold on;
 title('6 NM, synaptic weight: sim vs expt')
 ylabel('w av(i, o), expt. rescaled'); xlabel('wav sim, Trial 6') % now reversed
 ylim([0.5 3.5]); xlim([1 5.5]);
for i = 1:6 
 plot(w_trial6(i),wav_rescale(i),  'ob', 'MarkerSize', 10)
 text( w_trial6(i)+0.05,wav_rescale(i)+0.05, Acrn6{i})
    %text(wav_rescale(4)+0.05, wav_sim(4)+0.05, 'A9') % mark outliers
    %text(wav_rescale(2)+0.05, wav_sim(2)+0.05, 'A32V')
end
%%  vs. raw w_av
figure; plot(wav_calc, wav_sim, 'ob')
 title('6 NM, synaptic weight: sim vs expt')
 xlabel('w_av(i, o), expt.'); ylabel('w_av sim, Trial 6')
 xlim([0.5 3]); ylim([1 5.5]);
 % vs. w_max(i,o)
 figure; plot(wav_max, wav_sim, 'sk', 'MarkerSize', 10)
 title('6 NM, synaptic weight: sim vs expt - max(i, o)')
 xlabel('w_max(i, o), expt.'); ylabel('w_av sim, Trial 6')
 xlim([1 5]); ylim([1 5.5]);
 
 %% f-i vs. w_each
 figure; plot(wav_sim, fi, '-ob'); hold on;
 title('6 NM, synaptic weight: f-i vs. sim & expt')
 xlabel('w_a_v(i, o), expt.'); ylabel('f (Hz)')
 xlim([0.5 5.5]); ylim([4.5 35]);
 pause
 plot(wav_rescale, fi, '-sk', 'MarkerSize', 10)
 % confused plot??
 
 %%  w_each vs. f-i
 figure; plot(fi, wav_sim, '-ob'); hold on;
 title('6 NM, synaptic weight: sim & expt vs. f-i ')
 ylabel('w-av sim, expt.(i, o)'); xlabel('f (Hz)')
 ylim([0.5 5.5]); xlim([4.5 35]);
 pause
 plot(fi, wav_rescale, '-sk', 'MarkerSize', 10)
 
 %% reversed
 figure; plot(wav_sim, wav_rescale, 'ob')
 title('6 NM, synaptic weight: expt vs sim')
 ylabel('w av(i, o), expt. rescaled'); xlabel('wav sim, Trial 6')
 ylim([0.5 3.5]); xlim([1 5.5]);
  %text(wav_sim(4)+0.05,wav_rescale(4)+0.05,  'A9')
  %text(wav_sim(2)+0.05, wav_rescale(2)+0.05,  'A32V')
 text(wav_sim(5)+0.05, wav_rescale(5)+0.05, 'A46D')
 
 
 %% 4.1  w-i,o / neuron
 % local w-i,o / # neurons in anat area
 wav_expt = [2.6000 1.7400 1.1400 0.9200 0.7700 1.4200]; % av(w/k, i, o)
   % calc in 116Node..xl, Table S2.
 win = [0.12, 6.59, 0.21, 0.62, 2.54, 1.21];  % (*10^-2)  116Node..xl, tab2
 won = [1.10, 0.51, 0.97, 0.93, 0.38, 0.25];
  %wavn=(win + won)/2;
 w_trial6 =[3, 5, 3, 1.5, 3.5, 3];  % wb[8], p43
 Acrn6={'A10', 'A32V', 'A32','A9', 'A46D','A11'};
 
 figure; plot(w_trial6, win, 'ob', 'MarkerSize', 12); hold on
  plot(w_trial6, won, 'sk', 'MarkerSize', 12)
   %plot(w_trial6, won, '^r', 'MarkerSize', 12)  % av of #2 is way too low!
 title('6 NM, synaptic weight: w/#n vs w-trial')
 xlabel('w_a_v sim, Trial 6'); ylabel('w-i,o/ #n')
 xlim([0.5 5.5]);  ylim([-0.2 7]);
  %text(w_trial6(4)+0.05, win(4)+0.05,  'A9')
  %text(w_trial6(2)+0.05, win(2)+0.05,  'A32V')
  %text(w_trial6(5)+0.05,win(5)+0.05, 'A46D')
 % labels
 for i = 1:6 
  %plot(w_trial6(i), win(i), 'ob')
  text(w_trial6(i)+0.1, win(i)+0.05, Acrn6{i})
 end
 
  % legend off
  
 %% double plot
 figure; subplot(2,1,1); plot(w_trial6, win, 'ob', 'MarkerSize', 14); %hold on
 ylabel('w-i/ #n'); title('6 NM, synaptic weight: w/#n vs w-trial'); 
 xlim([1 5.5]); ylim([-1 8]);
   % labels
   for i = 1:6 
    text(w_trial6(i)+0.1, win(i)+0.05, Acrn6{i})
   end
 subplot(2,1,2); plot(w_trial6, won, 'sk', 'MarkerSize', 10)
 xlabel('w_a_v sim, Trial 6'); ylabel('w-o/ #n')
 xlim([1 5.5]);  ylim([-0.2 3]);
   for i = 1:6 
    text(w_trial6(i)+0.1, won(i)+0.05, Acrn6{i})
   end
   
   % clear wav*
 
%% 4.1.1  w-i,o / neuron // reversed plots
 % local w-i,o / # neurons in anat area
 wav_expt = [2.6000 1.7400 1.1400 0.9200 0.7700 1.4200]; % av(w/k, i, o)
   % calc in 116Node..xl, Table S2.
 win = [0.12, 6.59, 0.21, 0.62, 2.54, 1.21];  % (*10^-2)  116Node..xl, tab2
 won = [1.10, 0.51, 0.97, 0.93, 0.38, 0.25];
 w_trial6 =[3, 5, 3, 1.5, 3.5, 3];  % wb[8], p43
 Acrn6={'A10', 'A32V', 'A32','A9', 'A46D','A11'};
 
 figure; plot(w_trial6, win, 'ob', 'MarkerSize', 12); hold on
 % plot(w_trial6, won, 'sk', 'MarkerSize', 12)
   %plot(w_trial6, won, '^r', 'MarkerSize', 12)  % av of #2 is way too low!
 title('6 NM, synaptic weight: w/#n vs w-trial')
 xlabel('w_a_v sim, Trial 6'); ylabel('w-i,o/ #n')
 xlim([0.5 5.5]);  ylim([-0.2 7]);
  %text(w_trial6(4)+0.05, win(4)+0.05,  'A9')
  %text(w_trial6(2)+0.05, win(2)+0.05,  'A32V')
  %text(w_trial6(5)+0.05,win(5)+0.05, 'A46D')
 % labels
 for i = 1:6 
  %plot(w_trial6(i), win(i), 'ob')
  text(w_trial6(i)+0.1, win(i)+0.05, Acrn6{i})
 end
 
  %% legend off
  
 %% double plot
 figure; subplot(2,1,1); plot(win, w_trial6, 'ob', 'MarkerSize', 14); %hold on
 xlabel('w-i/ #neuron'); title('6 NM, synaptic weight:  w-trial vs w/#n'); 
 ylim([1 5.5]); xlim([-0.2 8]); ylabel('w_a_v sim, Trial 6');
   % labels
   for i = 1:6 
    text(win(i)+0.05, w_trial6(i)+0.1, Acrn6{i})
   end
 subplot(2,1,2); plot( won, w_trial6,'sk', 'MarkerSize', 10)
 ylabel('w_a_v sim, Trial 6'); xlabel('w-o/ #neuron')
 ylim([1 5.5]);  xlim([-0.2 3]);
   for i = 1:6 
    text( won(i)+0.05, w_trial6(i)+0.1,Acrn6{i})
   end
   % legend off
   % clear wav*
   
 %% 4.1.2  w-av in sim vs w-av from w(i,o)/k data
wik = [0.522, 3.22, 0.485, 0.736, 1.3, 2.346];  % w-in/k-in (*10^3)
wav_calc = [2.6000 1.7400 1.1400 0.9200 0.7700 1.4200]; % av(i,o/k) calc in 116Node..xl, Table S2.
wav_rescale = wav_calc/0.77; % scal from min value
w_trial6 =[3, 5, 3, 1.5, 3.5, 3];  % wb[8], p43
fi = [32.1 5.3 16.0 9.9 5.9 16.1]; % from 6NM sim, Table 1, mss. (1/24)
Acrn6={'A10', 'A32V', 'A32','A9', 'A46D','A11'};

% double plot
 figure; subplot(2,1,1); plot(wik, w_trial6, 'ob', 'MarkerSize', 14); %hold on
 xlabel('w-in/ #links in'); title('6 NM, synaptic weight: w-trial vs w/k '); 
 ylim([1 5.5]); xlim([0 4]); ylabel('w_a_v sim, Trial 6');
   % labels
   for i = 1:6 
    text(wik(i)+0.05, w_trial6(i)+0.1, Acrn6{i})
   end
 subplot(2,1,2); plot( wav_calc, w_trial6,'sk', 'MarkerSize', 10)
 ylabel('w_a_v sim, Trial 6'); xlabel('Av(w-i,o/ #links)')
 ylim([1 5.5]);  xlim([0 4]);
   for i = 1:6 
    text( wav_calc(i)+0.05, w_trial6(i)+0.1,Acrn6{i})
   end
   % legend off
   % clear wav*
