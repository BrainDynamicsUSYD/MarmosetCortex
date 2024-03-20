function dy = JR12a(t, y, Pulse, Ip, p, I12, a, b, A, B, C, Cfac, r)
   % JR DE's for 2 nodes: "Ip" becomes delta-y input; keep noise for now
   % Version JR12a carries params as arg: a, b, A, B - to adj freq band.
   % Inputs: Pulse fn, Ip (spikes), p (ran noise), I12 (i-j coupling)
   % S[V(t)] for internal feedbacks: classic sigmoid form
    % nb NM-NM coupling, I12, only goes into y4 (ie dy1/dy)
   % params needed in ODE's rhs function % Model parameters (orig JR)
v0=6;  % (mV) midpoint (switching threshold V) of sigmoid fn, for V -> pulse conversion 
vm=5;    % (s^-1) max firing rate for sigmoid fn; also called 2e-0
  %r=0.3;   % activation rate;  Steepness of sigmoid function- argument now
  %r=0.56;   % as used by Goodfellow('11), Friston
  %C=135;   % basic (scale) feedback connection strength; was Default JR value 
  %C=C*Cfac; % vary All  gain (x V & rates): e & i balanced
C1=C;     % Pyramidal to excitatory interneuron connection [strong feedback]
C2=0.8*C; % (108)   Excit interneuron to [fwd] pyramidal connection
C3=0.25*C; %5.0; % (33.75)  Pyramidal to inhib interneuron connection [weaker feedback]
C4=0.25*C;
C1=C1*Cfac; % vary feedback gain (x V): inhib
C3=C3*Cfac;
  %C2=C2*Cfac; % vary gain (x rates): excit
  %C4=C4*Cfac; % excit
  % defaults:
   %A=3.25; B=22; % Max PSP amplitude(mV) for pyramidal (A), inhib (B) & excit (kA * A) interneurons
   %a=50; %b=100;   % std values (s^-1)
% apply sigmoid filter [V->rate]: to Vspikes, pulse & coupling, at this t(i);
  % but not to pulse fn - which is already a firing rate
%Ii=0.2*(Ip +Pulse); % ext. excit stim of inhib population, via S(v); or zero (orig JR)
  % ?? apply to interactions?, I12 also: kills osc!
Ii=0.0; % default case, no stim to inhib popn

dy=zeros(6,1); % set up col vec
double(dy);
  % encodes std. Jansen-Rit ODEs for neural mass model; not global params;  
     % w inputs to y0, ie y(1)(excit)
        dy(1) = y(4); % ie y0'
        dy(2) = y(5);    % y1'
        dy(3) = y(6);    % y2'
        %[a b A B] % debug
        dy(4) = (A*a*(0  +(vm/(1+exp(r*(v0-y(2)+y(3)))))) ...  % no I12 to y0
                -(2*a*y(4))-(a^2*y(1)) );  % ie y0''
        dy(5) = (A*a*(Pulse + p +Ip + I12 +(C2*vm/(1+exp(r*(v0-C1*y(1))))) ) ... % Ip, ran & I12 to y1 (e)
                -(2*a*y(5)) -(a^2*y(2)) );
        dy(6) = (B*b*(Ii +(C4*vm/(1+exp(r*(v0-C3*y(1)))))) ...  % no Ip to y2 (i): test add p
                -(2*b*y(6) )-(b^2*y(3)) );
        
end