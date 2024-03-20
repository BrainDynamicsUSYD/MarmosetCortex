function dy = JR12v5b(t, y, Pulse, Ip, p, I12, a, b, A, B, C, Cfac, r, wbar)
   % JR DE's RHS: goes w twoNM...V5a.m : for erfc(1/V) form of S[v] ~ N(wt, wbar) 
   % Version JR12a carries params as arg: a, b, A, B - to adj freq band.
   % Inputs: Pulse fn, Ip (spikes), p (ran noise), I12 (i-j coupling)
   % S[V(t)]: int feedback: S[V] ~ erfc(1/V) form 
    % nb NM-NM coupling, I12, only goes into y4 (ie dy1/dy)
  %y' % debug
 % params needed in ODE's rhs function % Model parameters (orig JR)
v0=6;  % (mV) midpoint (switching threshold V) of sigmoid fn, for V -> pulse conversion 
vm=5;    % (s^-1) max firing rate for sigmoid fn; also called 2e-0
  %r=0.3;   % activation rate;  Steepness of sigmoid function- pass as argument now
  %C=135;   % basic (scale) feedback connection strength; - pass as argument now
C1=C;     % Pyramidal to excitatory interneuron connection [strong feedback]
C2=0.8*C; % (108)   Excit interneuron to [fwd] pyramidal connection
C3=0.25*C; % (33.75)  Pyramidal to inhib interneuron connection [weaker feedback]
C4=0.25*C;
C1=C1*Cfac; % vary feedback gain (x Voltage): incrs excitn
C3=C3*Cfac;
  %C2=C2*Cfac; % vary gain (x rates): inhib
  %C4=C4*Cfac; % excit
  % defaults:
   %A=3.25; B=22; % Max PSP amplitude(mV) for pyramidal (A), inhib (B) & excit (kA * A) interneurons
   %a=50; %b=100;   % std values (s^-1)
% apply sigmoid filter [V->rate]: to Vspikes, pulse & coupling, at this t(i);
  % but not to pulse fn - which is already a firing rate
Ii=0.2*(Ip +Pulse); % ext. excit stim of inhib population, via S(v); or zero (orig JR)
%Ii=0.0; % default case

dy=zeros(6,1); % set up col vec
double(dy);
  % encodes Jansen-Rit ODEs for neural mass model;  
        dy(1) = y(4); % ie y0'
        dy(2) = y(5);    % y1'
        dy(3) = y(6);    % y2'
          dyi=y(2)-y(3); % delta-y : feedforward    
           % prev. incl "switch" - cf OLD form
          argm=r*(wbar -v0/dyi); Sfn= (1.0-erfc(argm)/2); % erfc form for wt; avoid v=0          
        dy(4) = (A*a*(0  + vm*Sfn) ...  % no I12 to y0
                -(2*a*y(4))-(a^2*y(1)) );  % ie y0''
        yi=y(1); % feedback to both e & i from y0.    
        argm=r*(wbar -v0/(C1*yi)); Sfn5= (1.0-erfc(argm)/2); % erfc form for N(wt); avoided v=0
        argm=r*(wbar -v0/(C3*yi)); Sfn6= (1.0-erfc(argm)/2);     
        dy(5) = A*a*(Pulse + p +Ip + I12 + C2*vm*Sfn5 )... % y0, Ip, ran & I12 to y1 (e)
                -2*a*y(5) -a^2*y(2);            
        dy(6) = B*b*(Ii +C4*vm*Sfn6)  ...  % y0, no Ip or p-ran to y2 (i)
                -2*b*y(6) -b^2*y(3);
        
end