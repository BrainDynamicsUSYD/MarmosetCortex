% bap: Marmoset brain: test calc on data  V3 (20/8/20 - 2/24)
% based on: MarmosetReadDataV1.m (27/9/19) & V2 ('21)
%    use updated data: Majka & Rosa, etc, Nature Comm (2020);
%     & calc in XL: Marmoset_raw1_2_sort.xlsx 
   % V2. read prev  calc, cf V1, saved as Anew.mat
   % V3. clean up, with comments & list data files (2/24)

% Code cells (sections)
% 0.0, 0.01, 0.02 read data, info 
% #0.1   read rawLN2.csv - data from XL workbook " ":  [ Node-Id, TotLN, LNe, LNi, Av-LNtot, Av-LNi] 
%           Calc #LN(i-j) from raw data   & calc LNav ; that is "V2b" processing - cf Logbook

% #0.3 dependencies - LN,i (in col 4) on InjVol (repeats), etc
%  0.34 LNe,i vs Target Vol ~ Fig S6  [line 250]
%  0.35 LNe vs Injn Vol ~ Fig S5      [line 329]
% #0.4 linear model fits : all repeats [line 443]

% #1.0,  calc new Adj, using LNe per target: Anew =rescale Adj(LNe)  [line 500]
%         ie.  Calc #LN(i-j) from raw LN data & FLNe; save/ read in as Anew;    
%  1.0.1 read prev. calc of LN per target & Anew
%  1.3 wt Degree (Strength)
%  1.4 un-wt Deg (k) & 1.5 check dependencies (k)
  % 1.8.3 Calc Dist(i-j); & wt vs Dist
 % 1.8.4a Plot LinkWts vs Dist for rescaled data: Fig S7 [out of order]

% 2.0  Examine Hubs: links dist (rescaled)
%  Appx. X.0, line 10097  Clustering coeff.
%        X.1 Rich club calc.

%   Data files: AdjFullMarmoset.csv, VolsMarmoset.csv, CoordsMarmoset.csv, InjnRepeats.csv, 
%   MarmosetLabels139.csv, rawLN2.csv, rawLN2.csv  

%% 0.0 >>>>   Set Up, read some data 
 %clear all; close all
 %addpath(genpath('/bap_working/MatLabfiles/MarmosetBrain')); % Data files;  include sub-directories 
  %addpath('/bap_working/MatLabfiles/MarmosetBrain/CodesOther'); % for other codes & dependencies
addpath(genpath('/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain')); % Data files;  include sub-directories 
%
fprintf('\n >> Vis Marmoset Brain: : test calcs on loaded data \n')
% nb. Afull: calc prev. (116 x 116)  % read raw A (all S, some T); no self-loops
             % [n, ~] = size(Afull)  % number of Sources, 
             % its original FLNe (fractions)
Afull= csvread('AdjFullMarmoset.csv'); % calc prev. (116 x 116)  % read raw A (all S, some T); no self-loops
[n, ~] = size(Afull) 
% and Actx - square A(55 x 55); selection of raw A (both S&T); no self-loops

% set up larger figs for later
figwidth = 1024; figheight = 896; figposition = [100, 100, figwidth, figheight];

% read other info
NodeVols=csvread('VolsMarmoset.csv'); % (mm^3 - calc from bv.img[.nii]))
 % NodeVoxels=csvread('VoxelsMarmoset.csv'); % in voxel space, if needed

 %% printing
 % save figure to file  for mss.
  hh =gcf; print(hh, 'FigS4.tif', '-dtiff' , '-r300'); % for mss.
 % print(hh, 'FigXX.tif', '-dtiff' , '-r300','-opengl'  ) % 300dpi, OpenGL renderer, .tiff,
%% 0.01  Read, InjVol, examine spread of Injection Volumes, repeat expts:
% data in Marmoset_raw_1_2_sort.xlxs
%InjVols=csvread('InjnVols.csv'); % 143 cases: many repeats at 55 target sites
InjVols=csvread('InjnRepeats.csv');
figure; plot(InjVols(:,1), InjVols(:,2), '.', 'MarkerSize', 16);
title('Marmoset Brain:  Injections volumes (repeats) at 55 areas  ')
ylabel('Injn Vol (mm ^3 ) '); xlabel('Node ID # ');

figure; plot(InjVols(:,1), InjVols(:,3), '.', 'MarkerSize', 16);
title('Marmoset Brain:  Injections volumes (average) at 55 areas  ')
ylabel('Av Injn Vol (mm ^3 ) '); xlabel('Node ID # ');

InjVolsAv=[InjVols(:,1), InjVols(:,3)];
the55=find(InjVolsAv(:,2)); % eliminate repeats
InjVolsAv55=InjVolsAv(the55,:); % & pick out the averages
 %clear InjVolsAv the55 
 % nb. no InjVol entry for the non targets [amongst the 116 nodes]
 
%% 0.02 read Node Acrn, Labels 
% read the 116 Acrn for the Adj (Inj Source)
 % reader skips blanks, so load 'na' to preserve numbering
filename = '/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain/MarmosetAcrn116.csv';
delimiter = ',';  fileID = fopen(filename,'r');  formatSpec = '%s%[^\n\r]';
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false); fclose(fileID); % these are Cells, containing strings
Acrn = dataArray{1, 1}(1:end); % these are cells; 
clearvars filename delimiter formatSpec fileID dataArray ans
Acrn{107} % check this is 'V1': ok
% 1.1a read the 139 Atlas labels
 % reader skips blanks, so load 'na' to preserve numbering
filename = '/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain/MarmosetLabels139.csv';
delimiter = ',';  fileID = fopen(filename,'r');  formatSpec = '%s%[^\n\r]';
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false); fclose(fileID); % these are Cells, containing strings
Labels = dataArray{1, 1}(1:end); % these are cells; 
clearvars filename delimiter formatSpec fileID dataArray ans
 Labels{126} % check this is 'V1': ok
 
%% 0.03 node coords - prev calc from Atlas volumes;
 % cf Appx. A.1 of MarmosetReadDataV1.m  % from 3D vol: 'atlas_segmentation.nii'
NodeCoord=csvread('CoordsMarmoset.csv'); % nb. some absent: NaN 
   % adjust z-origin, to match Paxinos atlas??
    % NodeCoord(:,3)=NodeCoord(:,3)-20; % origin at V (top) & points down
max(NodeCoord(:, 3)) % z

%% 0.1 read details of LN, (V2b)  w repeat injn
fprintf('\n Calc #LN(i-j) from raw data \n')
 % data is [ Node-Id, TotLN, LNe, LNi, Av-LNtot, Av-LNi], with ~2-9 repeats; append av
  % from tab "LN_calc" of Marmoset_raw1_2_sort.xlxs
 % calc using the average LN & LNi repeats (x2-9) for the 55 targets - in cols [5, 6]
   % as calc in XL.
rawLN=csvread('rawLN2.csv'); % 6 cols copied from XL, cf above
LNtot=sum(rawLN(:,5)); % tot of all mean(LN-repeats) over all Targets; now 905k
% nb. Col 3 not reliable: may be 1st entry of LNi, not the av.
LNav = zeros(116,3); % to store av values [FLNtot, FLNe FLNi] for ea target,
  % with 0's for non targets
count=1; % nb. need to skip over repeats - just get the av.
for i =1:143
   if rawLN(i,5) > 0 % then av for these Targets recorded
       ii=rawLN(i,1); % index recorded in col 1 (w repeats)
        %[ii rawLN1(i,6)] % debug
        % nb. LNe = LNtot - LNi (col5 - col6)
       LNav(ii,:) = [rawLN(i,5), (rawLN(i,5)-rawLN(i,6)), rawLN(i,6)]; % calc FLNe also
       count=count+1;
   end
end
count % debug - check
clear count
nodeID=[1:n]'; % for 116 nodes, col vec
LNavID = [nodeID, LNav]; % preppend ID for later use nb. 0 for non-55Targets

[sum(rawLN(:,5)) sum(rawLN(:,5)-rawLN(:,6)) sum(rawLN(:,6)) ]
 % 905k , 688k, 217k : check grand totals of LN : ok
 
%% 0.1a plots, to check
figure; plot(rawLN(:,1), rawLN(:,3), '.', 'MarkerSize', 16); 
 title('Marmoset LN,e  vs Node ID (repeats) ')  
 xlabel('Node ID# '); ylabel(' LN,e (extrinsic) ');

 %% 0.3 dependencies - LN,i (in col 4) on InjVol (repeats), etc
figure; hist(InjVols(:,2),50); title('Marmoset InjVol (repeats) distribution  ') 

figure; plot(InjVols(:,2), '.', 'MarkerSize', 16);
 title('Marmoset InjVol (repeats) vs NodeID (all) ') 
 ylabel('Injn Vol (mm ^3) '); xlabel('Node ID # ');
  % LNi vs InjVol (repeated expts)
figure; plot(InjVols(:,2), rawLN(:,4), '.', 'MarkerSize', 16); % read in at #0.1
 title('Marmoset LN,i vs InjVol (repeats, all) ')  
 xlabel('Injn Vol (mm ^3 ) '); ylabel(' LN,i (intrinsic) ');
% 
figure; plot(InjVols(:,2), (rawLN(:,4)./NodeVols(rawLN(:,1)) ), '.', 'MarkerSize', 16);
 title('Marmoset: LN,i/ NodeVol vs InjVol (repeats, all) ')  
 xlabel('Injn Vol (mm ^3 ) '); ylabel(' LN,i/ NodeVol ');
   % fairly const: ~ av 285 +/- 395 ;  fit ~ 260 + 0.3*260 * vol
   [ mean(rawLN(:,4)./NodeVols(InjVols(:,1)) )  std(rawLN(:,4)./NodeVols(InjVols(:,1)) ) ]
%
 figure; plot(InjVols(:,2), (rawLN(:,4)./InjVols(:,2) ), '.', 'MarkerSize', 18);
 title('Marmoset: LN,i/ InjVol vs InjVol (repeats, all) ')  
 xlabel('Injn Vol (mm ^3 ) '); ylabel(' LN,i/ InjVol ');
% 
 % figure; plot( NodeVols(InjVols(:,1)), (rawLN(:,4)./InjVols(:,2) ), '.', 'MarkerSize', 16);
 figure; plot( NodeVols(rawLN(:,1)), (rawLN(:,4)./InjVols(:,2) ), '.', 'MarkerSize', 16);
 title('Marmoset: LN,i/ InjVol vs NodeVol (repeats, all) ')  
 xlabel('Node Vol (mm ^3 ) '); ylabel(' LN,i/ InjVol ');
  % label 3 outliers : Node 37, row 39; Node 32, row 30  Node 108, row 129
  text( (NodeVols(rawLN(39,1))+3), (rawLN(39,4)./InjVols(39,2) ), Acrn{37},'FontSize',10); % label
  text( (NodeVols(rawLN(30,1))+3), (rawLN(30,4)./InjVols(30,2) ), Acrn{32},'FontSize',10);
  text( (NodeVols(rawLN(129,1))+3), ((rawLN(129,4)./InjVols(129,2)-0.3e5) ), Acrn{108},'FontSize',10);
   % get  corrds right?
  [ mean(rawLN(:,4)./InjVols(:,2)), std(rawLN(:,4)./InjVols(:,2)) ]
  
%% 0.31  LN,e  repeats, dependencies
 % plot all measurements (repeats) (LNe is col 3)
figure; plot(InjVols(:,2), rawLN(:,3), '.', 'MarkerSize', 14);
 title('Marmoset LN,e vs InjVol (repeats, all) ')  
 xlabel('Injn Vol (mm ^3 ) '); ylabel(' LN,e (extrinsic) ');
% 
figure; plot(InjVols(:,2), (rawLN(:,3)./NodeVols(rawLN(:,1)) ), '.', 'MarkerSize', 14);
 title('Marmoset density: LN,e/ NodeVol vs InjVol (repeats, all) ')  
 xlabel('Injn Vol (mm ^3 ) '); ylabel(' LN,i/ NodeVol ');
   % fairly const: ~ av 285 +/- 395 ;  fit ~ 260 + 0.3*260 * vol
   [ mean(rawLN(:,3)./NodeVols(rawLN(:,1)) )  std(rawLN(:,3)./NodeVols(rawLN(:,1)) ) ]
%  
figure; plot(NodeVols(rawLN(:,1)), rawLN(:,3), '.', 'MarkerSize', 14);
 title('Marmoset LN,e vs Node Vol (repeats, all) ')  
 xlabel('Node Vol (mm ^3 ) '); ylabel(' LN,e (extrinsic) ');

 %% xx 0.31.1  3-way comparisons: LN,int
figure; plot3(InjVols(:,2), NodeVols(InjVols(:,1)), rawLN1(:,3), '.', 'MarkerSize', 11); % nb Node ID in 1st col
 title('Marmoset LN,e vs InjVol & NodeVol (repeats, all) ')  
 xlabel('Injn Vol (mm3) '); ylabel('Node Vol (mm3) '); zlabel(' LN,i (extrinsic) '); grid on
% nb k1 : Links In,Out calc at #1.2 below
figure; plot3(InjVols(:,2), kIn1(InjVols(:,1)), rawLN1(:,3), '.', 'MarkerSize', 11); % nb Node ID in 1st col
 title('Marmoset LN,e vs InjVol & k-In (repeats, all) ')  
 xlabel('Injn Vol (mm3) '); ylabel('k-In '); zlabel(' LN,i (extrinsic) '); grid on
%
 figure; plot3(kIn1(InjVols(:,1)), NodeVols(InjVols(:,1)), rawLN1(:,3), '.', 'MarkerSize', 11); % nb Node ID in 1st col
 title('Marmoset LN,e vs InjVol & NodeVol (repeats, all) ')  
 xlabel('k-In '); ylabel('Node Vol (mm3) '); zlabel(' LN,i (extrinsic) '); grid on

 %% 0.32 dependencies for the av LNi, LNe
 % find the unique 55 targets (amongst repeat expts)
 % nb LNeAv here is same as LNav above
the55=find(rawLN(:,5)); % eliminates repeats
 % col 6 is av-LNi
LNeAv=rawLN(the55,5)-rawLN(the55,6); % calc LN,e =(LN,tot - LN,i); NodeId listed in col 1
LNiAv=rawLN(the55,6);
theNodeVols=NodeVols(rawLN(the55,1)); % & pick out the av vol amongst repeats, for ea ID
theInjVols=InjVolsAv(the55,2);  % and av InjVol (amongst repeats for ea node)

 figure; plot(theNodeVols, LNiAv, '.', 'MarkerSize', 11);
 title('Marmoset av LN,i vs Node Vol (55 tragets) ')  
 xlabel('Node Vol (mm3) '); ylabel(' Average LN,i (intrinsic) ');
  legend('av LN,i', 'linear fit')
  
 figure; plot(theInjVols, LNiAv, '.', 'MarkerSize', 11);
 title('Marmoset av LN,i vs Inj Vol (55 tragets) ')  
 xlabel('Inj Vol (mm3) '); ylabel(' Average LN,i (intrinsic) ');
  legend('av LN,i', 'linear fit')
  
  figure; plot(theNodeVols, LNeAv, '.', 'MarkerSize', 11);
 title('Marmoset av LN,e vs Node Vol (55 tragets) ')  
 xlabel('Node Vol (mm3) '); ylabel(' Average LN,e (extrinsic) ');
  legend('av LN,e', 'linear fit')
  
 figure; plot(theInjVols, LNeAv, '.', 'MarkerSize', 11);
 title('Marmoset av LN,e vs Node Vol (55 tragets) ')  
 xlabel('Inj Vol (mm3) '); ylabel(' Average LN,e (extrinsic) ');
  legend('av LN,e', 'linear fit')

  [ mean(LNiAv), std(LNiAv) ]
  [ mean(LNeAv), std(LNeAv) ]
  
  
  %InjAvRows=find(InjVolsAv(:,2)); % ? finds 56, not 55?
  %tmp=InjVolsAv(InjAvRows,:);
 % LN,i 
 figure; plot(theNodeVols, (LNiAv./theInjVols), '.', 'MarkerSize', 11);
 title('Marmoset av (LNi/InjVol) vs Node Vol (55 tragets) ')  
 xlabel('Node Vol (mm3) '); ylabel(' Average LNi/InjVol ');
 
 figure; plot(theInjVols, (LNiAv./theInjVols), '.', 'MarkerSize', 11);
 title('Marmoset av (LNi/InjVol) vs Inj Vol (55 tragets) ')  
 xlabel('Inj Vol (mm3) '); ylabel(' Average LNi/InjVol ');
 
 [ mean(LNiAv./theInjVols), std(LNiAv./theInjVols) ]
 % LN,e:
 figure; plot(theNodeVols, (LNeAv./theInjVols), '.', 'MarkerSize', 11);
 title('Marmoset av (LNe/InjVol) vs Inj Vol (55 tragets) ')  
 xlabel('Node Vol (mm3) '); ylabel(' Average LNe/InjVol ');
 
 figure; plot(theInjVols, (LNeAv./theInjVols), '.', 'MarkerSize', 11);
 title('Marmoset av (LNe/InjVol) vs Inj Vol (55 tragets) ')  
 xlabel('Inj Vol (mm3) '); ylabel(' Average LNe/InjVol ');
 
 [ mean(LNeAv./theInjVols), std(LNeAv./theInjVols) ]

figure; loglog(theInjVols, (LNeAv./theInjVols), '.', 'MarkerSize', 11);
 title('Marmoset av (LNe/InjVol) vs Inj Vol (log log plot) ')  
 xlabel('log Inj Vol (mm3) '); ylabel('log Av LNe/InjVol '); 
 
 %% 0.34 LNe,i vs Target Vol ~ Fig S6 
 %  LNe vs Target vol
figure; hold on
for i=1:length(rawLN)
    plot(NodeVols(rawLN(i,1)), rawLN(i,3), '.', 'MarkerSize', 12); hold on % NodeID in col 1
end
title('Marmoset cortex:  LNe vs. TargetVol  ');
xlabel('Target Vol (mm ^3 ) '); hold off

%  LNe/InjVol vs Target vol
 LNeInj=zeros(143,1); % accum
figure; hold on; xlim([0 300]); % ylim([0 12.2e5]); % fix scale
title('Marmoset cortex:  LNe/ InjVol vs. TargetVol  ');
xlabel('Target Vol (mm ^3 ) '); ylabel('Counts/(mm ^3 )');
for i=1:length(rawLN)
    plot(NodeVols(rawLN(i,1)), (rawLN(i,3)/InjVols(i)), '.', 'MarkerSize', 24); hold on % NodeID in col 1
    LNeInj(i) =(rawLN(i,3)/InjVols(i));
     ii= rawLN(i); % ?? (i,1)
     [i ii] % debug
     text( (NodeVols(rawLN(i,1))+3), LNeInj(i), Acrn{ii},'FontSize',10); % label label
      pause % debug: check single points
end  
mean(LNeInj)
outList = find(LNeInj >= 1e6)' % outliers
rawLN(outList,1)' % as Node ID's
% % highlight outliers
for i = 1:length(outList)
    ii=rawLN(outList(i),1);
 text( (NodeVols(rawLN(ii,1))+3), LNeInj(ii), Acrn{ii},'FontSize',10); % label
end
 line([0, 300], [mean(LNeInj), mean(LNeInj)], 'Color', 'k' ) % mean
 line([0, 300], [(mean(LNeInj)+3*std(LNeInj)), (mean(LNeInj)+3*std(LNeInj))], ...
     'LineStyle', '--', 'Color', 'r' )
 hold off

%  LNi vs Target vol
figure; hold on
for i=1:length(rawLN)
    plot(NodeVols(rawLN(i,1)), rawLN(i,4), '.', 'MarkerSize', 12); hold on % NodeID in col 1
end
title('Marmoset cortex:  LNi vs. TargetVol  ');
xlabel('Target Vol (mm ^3 ) '); ylabel('Counts'); 

hold off

 %  LNi/InjVol vs Target vol
 LNiInj=zeros(143,1); % accum
figure; hold on; xlim([0 300]); ylim([0 12.2e5]); % fix scale
for i=1:length(rawLN) % scan rows of raw data 
    plot(NodeVols(rawLN(i,1)), (rawLN(i,4)/InjVols(i)), '.', 'MarkerSize', 24); hold on % NodeID in col 1
    LNiInj(i) =(rawLN(i,4)/InjVols(i)); % the ratio, to use below
    %ii=rawLN(i,1);
    % [i ii] % debug
    % text( (NodeVols(rawLN(i,1))+3), LNiInj(i), Acrn{ii},'FontSize',10); % label label
    %  pause % debug: check single points
end
title('Marmoset cortex:  LNi/ InjVol vs. TargetVol  ');
xlabel('Target Vol (mm ^3 ) '); ylabel('LNi/ InjVol  (Counts/(mm ^3 ))');
  
mean(LNiInj)
outList = find(LNiInj >= 5e5)' % outliers, as position in list
nodeIDs= rawLN(outList,1)' % as Node ID's
% % highlight & label outliers
for i = 1:length(outList)
    ii= rawLN(outList(i),1); % outList(i);
 %text( (NodeVols(rawLN(ii,1))+3), LNiInj(ii), Acrn{ii},'FontSize',10); % label
 text( (NodeVols(rawLN(outList(i),1))+3), LNiInj(outList(i)), Acrn{ii},'FontSize',10); % label
end
 line([0, 300], [mean(LNiInj), mean(LNiInj)], 'Color', 'k' ) % mean
 line([0, 300], [(mean(LNiInj)+3*std(LNiInj)), (mean(LNiInj)+3*std(LNiInj))], ...
     'LineStyle', '--', 'Color', 'r' )
hold off

% hh =gcf; print(hh, 'FigS6.tif', '-dtiff' , '-r300'); % for mss. Fig S6 
 
  %plot(NodeVols(rawLN(32,1)), (rawLN(32,4)/InjVols(rawLN(32,1))), 'ro', 'MarkerSize', 12); % A46D
  %(rawLN(:,4)./NodeVols(rawLN(:,1))

   clear outlist nodeIDs

%% 0.35 LNe vs Inj Vol ~ Fig S5 
%  LNe/ vs InjVol vol
figure; hold on;  xlim([0 2]);  ylim([0 5e4]); % fix scale
title('Marmoset cortex:  LNe vs. Injection Vol  ');
xlabel('Injection Vol (mm ^3 ) '); ylabel('Counts/(mm ^3 )');
for i=1:length(rawLN)
    plot(InjVols(i), rawLN(i,3), 'k.', 'MarkerSize', 20); hold on % NodeID in col 1
     %ii= rawLN(i,1);
     %[i ii] % debug
     %text( (InjVols(i)+0.03), rawLN(i,3), Acrn{ii},'FontSize',10); % label
     % pause % debug: check single points
end  
mean(rawLN(:,3))
outList = find(rawLN(:,3) >= 4e4)' % outliers
rawLN(outList,1)' % as Node ID's
% highlight & label outliers
for i = 1:length(outList)
    ii= rawLN(outList(i),1); % outList(i);
 %text( (NodeVols(rawLN(ii,1))+3), LNiInj(ii), Acrn{ii},'FontSize',10); % label
 %text( (InjVols(i)+0.05), rawLN(i,3),  Acrn{ii},'FontSize',10); % label
end

 %line([0, 2], [mean(rawLN(:,3)), mean(rawLN(:,3))], 'Color', 'k' ) % mean
  %line([0, 2], [(mean(rawLN(:,3))+3*std(rawLN(:,3))), (mean(rawLN(:,3))+3*std(rawLN(:,3)))], ...
   %  'LineStyle', '--', 'Color', 'r' )
% hh =gcf; print(hh, 'FigS5.tif', '-dtiff' , '-r300'); % for mss. Fig S5 


 %% Test plot for A10: linear trends?  2:7  % for A-10 : poor -too scattered??
%figure; hold on;  xlim([0 1.4]);  ylim([0 2e4]); % fix scale
%title('Marmoset A10:  LNe/ InjVol vs. TargetVol  ');
%xlabel('Injection Vol (mm ^3 ) '); ylabel('Counts/(mm ^3 )');
A10list =[3 2 4 5 6 7]; % needs to be ordered in incrs InjVol - for lib fit
for ii=1:length(A10list) % for A-10 
    i = A10list(ii);
    plot(InjVols(i), rawLN(i,3), 'r.', 'MarkerSize', 24); hold on % NodeID in col 1
     ii= rawLN(i,1);
     %[i ii] % debug
     text( (InjVols(i)+0.03), rawLN(i,3), Acrn{ii},'FontSize',10); % label
      %pause % debug: check single points
end 
 line([0, 2], [7599, (7599-2*2158)]) % lin fir for A10

  %  2:7  % for A-10
% A10 lin fit : rawLN(A10list,3) to InjVols(A10list)' : both as col vec
  %ffit = fit(InjVols(A10list)', rawLN(A10list,3), 'poly1')
    % >  p1 =  -2158  (-1.646e+04, 1.215e+04), ffit(x) = p1*x + p2 
     %   p2 =   7599  (-118.9, 1.532e+04) & 95% CI 
%clear A10list
%  116:121 % for V1 
%figure; hold on;  xlim([0 1.4]);  ylim([0 2e4]); % fix scale
%title('Marmoset V1:  LNe/ InjVol vs. TargetVol  ');
%xlabel('Injection Vol (mm ^3 ) '); ylabel('Counts/(mm ^3 )');
V1list =[116 117 118 119 120 121]; % needs to be ordered in incrs InjVol - for lib fit
for ii=1:length(A10list) % for A-10 
    i = V1list(ii);
    plot(InjVols(i), rawLN(i,3), 'b.', 'MarkerSize', 24); hold on % NodeID in col 1
     ii= rawLN(i,1);
     %[i ii] % debug
     text( (InjVols(i)+0.03), rawLN(i,3)+0.2, Acrn{ii},'FontSize',10); % label
      %pause % debug: check single points
end 
  line([0, 1.4], [-167.6, (-167.6 +1.4 *1.225e+04)], 'Color', 'b') % lin fit for V1
    % ffit = fit(InjVols(V1list)', rawLN(V1list,3), 'poly1') % lin fit for V1
    %  p1 =   1.225e+04  (7385, 1.712e+04)  % ffit(x) = p1*x + p2 
    %   p2 =      -167.6  (-2757, 2422) % & 95% CI
 %hold off
 
%  57:61 % for A6Va  : too cluttered near origin, small spread
%figure; hold on;  xlim([0 1.4]);  ylim([0 2e4]); % fix scale
%title('Marmoset V1:  LNe/ InjVol vs. TargetVol  ');
%xlabel('Injection Vol (mm ^3 ) '); ylabel('Counts/(mm ^3 )');
A6Vlist =[57:61]; % needs to be ordered in incrs InjVol - for lib fit
for ii=1:length(A46list) % for A-10 
    i = A6Vlist(ii);
    plot(InjVols(i), rawLN(i,3), 'r.', 'MarkerSize', 24); hold on % NodeID in col 1
     ii= rawLN(i,1);
     %[i ii] % debug
     text( (InjVols(i)+0.03), rawLN(i,3), Acrn{ii},'FontSize',10); % label
      %pause % debug: check single points
end 
   line([0, 0.3], [581.6, (581.6 + 0.3 *5.044e+04)], 'Color', 'r') % lin fit for V1
    % ffit = fit(InjVols(A6Vlist)', rawLN(A6Vlist,3), 'poly1') % lin fit for V1
    %   p1 =   1.021e+04  (-3.003e+04, 5.044e+04) % ffit(x) = p1*x + p2
    %   p2 =       581.6  (-7116, 8279) % & 95% CI
 %hold off  

% 122:128 & 130 % for V2, 
 % omit worst oulier 
%figure; hold on;  xlim([0 1.4]);  ylim([0 2e4]); % fix scale
%title('Marmoset V1:  LNe/ InjVol vs. TargetVol  ');
%xlabel('Injection Vol (mm ^3 ) '); ylabel('Counts/(mm ^3 )');
V2list =[123;127;122;125;126;124]; % ordered by Inj Vol
 %V2list =[122 : 128]; % needs to be ordered in incrs InjVol - for lib fit
for ii=1:length(V2list) % for V2 
    i = V2list(ii);
    plot(InjVols(i), rawLN(i,3), 'm.', 'MarkerSize', 24); hold on % NodeID in col 1
     ii= rawLN(i,1);
     %[i ii] % debug
     text( (InjVols(i)+0.03), rawLN(i,3)+0.2, Acrn{ii},'FontSize',10); % label
      %pause % debug: check single points
end 
   line([0, 0.7], [ -584.3, ( -584.3 +0.7 *2.818e+04)], 'Color', 'm') % lin fit for V2
     %ffit = fit(InjVols(V2list), rawLN(V2list,3), 'poly1') % lin fit for V1
     % raw fit %  p1 =   2.818e+04  (1.338e+04, 4.299e+04) % ffit(x) = p1*x + p2
               %  p2 =       993.6  (-3473, 5461) % & 95% CI
    % Ordered fit
       %p1 =   3.136e+04  (2.046e+04, 4.225e+04)
       %p2 =      -584.3  (-4100, 2932)
 %hold off
  
 
%% 0.4 linear model fits : all repeats
LNeAll=rawLN(:,2); % all repeats of LNe; vs InjVol
[ mean(LNeAll), std(LNeAll) ]
InjVolAll=InjVols(:,2);
 % cf LNeAv % av over repeats of 55 targets
mdl=fitlm(InjVolAll, LNeAll, 'linear'); % "x, y" / stats ToolBox
mdla=fitlm(InjVolAll, LNeAll, 'constant'); % no lin term
% Alt.
 %mdlb=fitlm(InjVolAll, LNeAll, 'linear', 'Exclude', [1, 2]); % exclude obvious outliers
mdl.Coefficients
mdl.LogLikelihood
mdl.Rsquared

mdla.Coefficients
mdla.Rsquared

% LNe1
LNe1All=LNeAll./InjVolAll;
[ mean(LNe1All), std(LNe1All) ]
mdl1=fitlm(InjVolAll, LNe1All, 'linear');
mdl1.Coefficients
mdl1.LogLikelihood
mdl1.Rsquared

mdl1a=fitlm(InjVolAll, LNe1All, 'constant');
mdl1a.Coefficients
mdl1a.Rsquared

%% 0.41 linear model fits : av for 55 targets
% av calc at #0.32 above
[ mean(LNeAv), std(LNeAv) ]
 % cf LNeAv % av over repeats of 55 targets
mdl=fitlm(theInjVols, LNeAv, 'linear'); % "x, y" / stats ToolBox
mdla=fitlm(theInjVols, LNeAv, 'constant'); % no lin term
% Alt.
 %mdlb=fitlm(InjVolAll, LNeAll, 'linear', 'Exclude', [1, 2]); % exclude obvious outliers
mdl.Coefficients
mdl.LogLikelihood
mdl.Rsquared

mdla.Coefficients
mdla.Rsquared

%% LNe
LNe1All=LNeAll./InjVolAll;
[ mean(LNe1All), std(LNe1All) ]
mdl1=fitlm(InjVolAll, LNeAll, 'linear');
mdl1.Coefficients
mdl1.LogLikelihood
mdl1.Rsquared

mdl1a=fitlm(InjVolAll, LNeAll, 'constant');
mdl1a.Coefficients
mdl1a.Rsquared

% >>.    >   >   >   >

%% 1.0 trial calc new Adj, using LNe per target
fprintf('\n Calc #LN(i-j) from raw LN data & FLNe \n')
 % calc LNav at #0.1, above   checked (19/2/24)
  % LNav is LNi,e,tot : 3 cols x 116 rows (nodes)
  % and LNavID has col1 as node ID [redundant]
 %targetList=find(A(:,1)); % tests % eg j=1
 %targetFLNe=A(targetList,1);
 %Tot1=sum(targetFLNe)  % debug:  is 1.0
 % OLd: % Anew1=LNav(1,2)*targetFLNe  % ie. tot # ext LN * Frn for this link
 
Anew=zeros(n);
% loop over cols (Targets), undo the fractions : get raw #LN ext labelled
for j=1:n
   Anew(:,j)= Afull(:,j)*LNav(j,2);  % mult raw, av LNe by orig FLNe for this source-target pair
     % target is col-j 
end
sum(sum(Anew))  % is 6.8822e+05  ok
 % alt defn, with max = 1.
 % Anew1=Anew*(1/sum(sum(Anew))); % rescale to frn of grand total counts (LNe)
%
posWts=Anew(find(Anew(:)>0)); % find the non-zero entries % 3474 of
 % sum(posWts) % check: ok
figure; hist(posWts, 100); title('Marmoset rescaled LNe distribution  ')

% nb. v steep exp: need log plot
wtlog=log10(posWts);
figure; hist(wtlog,50); title('Marmoset rescaled log-10 (LNe) distribution  ')
xlabel('log10(weight, LNe) '); ylabel('counts')
  % approx gaussian
% save figure to file  for mss.
hh =gcf; print(hh, 'FigS4.tif', '-dtiff' , '-r300'); % for mss.
 % print(hh, 'FigXX.tif', '-dtiff' , '-r300','-opengl'  ) % 300dpi, OpenGL renderer, .tiff,
%
posWts1=Anew(find(Anew(:)>0)); % overall normalised
posWts=Afull(find(Afull(:)>0)); % compare to original FLNe: normalised at ea. target
figure; plot(posWts, posWts1, '.', 'MarkerSize', 14);
 title('Marmoset rescaled LNe vs original FLNe (all) ')
 %axis([0 0.6 0 0.016])
   % nb different slopes
  % altv. yields same plot
 figure; %plot(Afull(:,1), Anew(:,1), '.', 'MarkerSize', 11); % Target #1
 plot(Afull(:), Anew(:), '.', 'MarkerSize', 14);  % all targets
 title('Marmoset rescaled LNe vs original FLNe (all Nodes) ')
 xlabel('FLNe, fraction'); ylabel('LNe, counts')
 %axis([0 0.6 0 0.016])

 %% 1.0.1 read in Anew V2b LNe form)
 % save as file: Anew= csvread('AdjMarmoset_Rescale2b.csv');
 
 %% 1.0.2 LinkList:  calc from Adj: best way to get log-10(wt) distribution
LinkList=adj2edgeL(Anew); % finds 3474 links
 max(LinkList(:,3)) % max wt ! 10k
 min(LinkList(:,3)) % non-zero:  0.0340 
 % v % check
  % eliminate wt:0 
  % posWts=find(LinkList(:,3)>0); %index,  find the non-zero entries % 3474 of
  % LinkList=LinkList(posWts,:);
  %clear posWts
[nl ~] =size(LinkList) % now 3474 non-0 wts
% save for InfoMap calc: .txt file, ' ' delim
% save lists for InfoMap (C++) calc
% dlmwrite( 'MarmosetPairList.txt', LinkList,'Delimiter',' ');  % InfoMap code wants .txt , space delim

figure; hist(LinkList(:,3), 50) % appears exp
logWts = log(LinkList(:,3)); % nb. all pos.

figure; hist(logWts, 50)

%% 1.3 wt Degree (Strength)
 % calc Anew at #1.0, above
DegIn=sum(Anew,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
DegOut=sum(Anew,2); % i->j; row sum 
 %figure; hist(DegIn); title('Marmoset rescaled wt-DegIn ')
 % figure; hist(DegOut); hold on; title('Marmoset rescaled wt-DegOut '); hold off
 % figure; plot(DegOut, DegIn, '.', 'MarkerSize', 9 ); title('marmoset: rescaled DegOut vs DegIn  ')
figure; hist([DegIn, DegOut]); hold on; title('Marmoset rescaled Strength: wt-DegIn, Out '); 
legend('wt-DegIn', 'wt-DegOut');  hold off
% Fig 1 for mss.
 %figure; stem (DegIn, '.'); hold on  %  ~ "wt-DegIn"
 %stem (DegOut, 'r.');  % wt-DegOut (LNe/mm^3 measure)
 %legend('wt-DegIn', 'wt-DegOut'); title('Marmoset rescaled (LNe) wt-Deg-In,Out '); hold off 
 
% exp fit?
logDegIn=log(DegIn); logDegOut=log(DegOut);
figure; hist(logDegIn, logDegOut, 50)
pd = fitdist(DegIn,'exponential')
scale = 68/DegInFit(1)
 % [min(DegIn) max(DegIn) ] % get scale 
DegInAxis=linspace(0,  max(DegIn), 50);
DegInFit = pdf(pd, DegInAxis );
figure; plot(DegInAxis, DegInFit)
plot(DegInAxis, scale*DegInFit, 'k-')
title('Marmoset rescaled Strength: wt-DegIn '); 
legend('wt-DegIn', 'exp fit');  
%
pd = fitdist(DegOut,'exponential')
DegOutAxis=linspace(0,  max(DegOut), 50);
DegOutFit = pdf(pd, DegOutAxis );
scale = 12/DegOutFit(1)  % get max from fig
figure; hist(DegOut, 50); hold on
plot(DegOutAxis, scale*DegOutFit, 'k-')
title('Marmoset rescaled Strength: wt-DegOut '); 
legend('wt-DegOut', 'exp fit'); 
%% 1.4 un-wt Deg 
fprintf('\n Calc un-wt Degree - from Adj \n')

[sii sjj v]=find(Anew);  
Alogical=full(sparse(sii, sjj, 1));    % "mask" of selected entries, to re-form sq array of 1's.
 % here thats not square?
tmp=zeros(116,1); % need dummy entries for cols 115, 116 - which have no inputs, or were not sampled (cf. A)
Alogical=[Alogical, tmp, tmp]; % now its 116 x 116
 % nb. need to use sparse to form the array
A1 = Alogical.*1;  % pick out the selected entries; yields a "full" array
figure; spy(A1)
clear Alogical sii sjj v tmp
length(find(A1(:))) % check 3474 links: ok
 %LinkList1=adj2edgeL(A1);
 %dlmwrite( 'MarmosetPairListR1.txt',LinkList,'Delimiter',' '); % InfoMap code wants .txt , space delim

% count links (<~ 50-150)
kIn1=sum(A1,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
kOut1=sum(A1,2); % i->j; row sum 
 clear A1
figure; hist(kIn1, 50)
title('Marmoset k: DegIn of binarised links')
axis([0 110 0 10]) % nb. big # at 0: not sampled

figure; hist(kOut1, 50); hold on
title('Marmoset k: DegOut of binarised links')
 %axis([0 80 0 130]); 
figure; stem (kIn1, '.'); hold on

stem (kOut1, 'r.');
legend('k-In1', 'k-Out1')
title('Marmoset Deg-In,Out of binarised links')

figure; plot(kIn1, kOut1, '.')
title('Marmoset Deg-Out vs. In of binarised links')

%% alt calc of LinkList
LinkList=adj2edgeL(A);  %LinkList=adj2edgeL(Afull);
  % use "filtered" LinklList (of 3462 links ~0)
  % LinkList=csvread('MarmosetLinkList.csv'); % i-j-wt
[nl ~]=size(LinkList)
wt1=ones(nl);
LinkList1=[LinkList(:,1:2) wt1]; % append wt=1
A1=edgeL2adj(LinkList1);
 % nb. need to use sparse to form the array
figure; spy(A1)
%
clear wt1 LinkList LinkList1
 % figure; spy(A1) % 1948 entries: banded & v sparse
 length(find(A1(:)) )

 %% 1.5 check dependencies K (# links)
figure; plot(NodeVols, kIn1, '.', 'MarkerSize', 11);
 title('Marmoset degree (k-In) vs Node vol ') 
  % linear
  
figure; plot(NodeVols, kOut1, '.', 'MarkerSize', 11);
 title('Marmoset degree (k-Out) vs Node vol ') 
   % fairly const; 2 outliers:  V1 & V2;  also: V3 (66mm^3: not injn site)

figure; plot(LNav(:,2), kIn1, '.', 'MarkerSize', 11);
 title('Marmoset degree (wtDeg-In) vs avLN-ext ') 
  %  linear
figure; plot(LNav(:,2), kOut1, '.', 'MarkerSize', 11);
 title('Marmoset degree (k-Out) vs avLN-ext ')
  % linear; flat line, weak slope 

  % nb InjVol at 55 targets only
 the55= InjVolsAv55(:, 1); % col 1
 figure; plot(InjVolsAv55(:, 2), kIn1(InjVolsAv55(:, 1)), '.', 'MarkerSize', 11);
 title('Marmoset degree (k-In) vs InjVol ') 
 figure; plot(InjVolsAv55(:, 2), kOut1(InjVolsAv55(:, 1)), '.', 'MarkerSize', 11);
 title('Marmoset degree (k-Out) vs InjVol ') 
  % nb. Vol fairly const across k-In or -Out
    % one outlier: #48 (A9) with 1.9 mm^3 (v large Injn ?)
    
% Sampling metric: num'r & denom:
   %SM= NodeVols(InjVolsAv55(:, 1))./kIn1(InjVolsAv55(:, 1)); % tests
   %figure; plot(SM, '.', 'MarkerSize', 11); hold on
   %plot(NodeVols(InjVolsAv55(:, 2)), 'g.', 'MarkerSize', 12);
  % calc the ration
  SM= (NodeVols(InjVolsAv55(:, 1))./kIn1(InjVolsAv55(:, 1)) )./InjVolsAv55(:, 2) ;
  figure; plot(SM, '.', 'MarkerSize', 11); hold on
  title('Marmoset Sampling Metric at Inj Vol sites ') 
 % [NodeVols(InjVolsAv55(:, 1))./kIn1(InjVolsAv55(:, 1))  InjVolsAv55(:, 2) ] % debug

%% 2.0 Volume & Surface area effects, densities:
NodeR=zeros(n,1); NodeAreas=zeros(n,1); % unit mm, mm^2 
for i=1:n
   NodeR(i)=(3*NodeVols(i)/(4*pi))^0.333;
   NodeAreas(i)=4*pi*NodeR(i)^2;
end
% Link density at source & target surface areas
 % nb. #73 is v sml vol, so an outlier
 % only 55 "In" measurements
kInDen=kIn1./NodeAreas; kOutDen=kOut1./NodeAreas; % Link "fluxes"
kInDenPos=kInDen(find(kInDen));
[ mean(kInDenPos)  mean(kOutDen)] % ~ 3, 5 /mm^2
%
figure; plot(kInDen, '.', 'MarkerSize', 11); hold on
axis([0 120 0 10])
plot(kOutDen, 'g.', 'MarkerSize', 11);
legend('Links-In/Area', 'Links-Out/Area');
figure; hist([kInDen, kOutDen]);  legend('Links-In/Area', 'Links-Out/Area');
figure; hist([kInDen, kOutDen], 200);  legend('Links-In/Area', 'Links-Out/Area');
axis([0 10 0 100]) % shows more detail

 % LabelledNeurons / Vol:
 
 
 
 
 
% OLD: from V1 code:

%% 3.0 LinkList:  calc from Adj: best way to get log-10(wt) distribution
LinkList=adj2edgeL(Anew); % finds 3474 links  % V2b now
  % eliminate wt:0 
  %posWts=find(LinkList(:,3)>0); % find the non-zero entries % 3474 of
  % LinkList=LinkList(posWts,:);
  %clear posWts
[nl ~] =size(LinkList) % now 3474 non-0 wts

% nb. Acrn-i, Acrn-j, wt] from .txt file & in LinkListACrn.csv
wt=LinkList(:,3);
wtlog=log10(wt);
figure; hist(wtlog,50); title('Marmoset  log10 (LNe1(i,j)) distribution  ')

% save for InfoMap calc
dlmwrite( 'MarmosetPairListRescale.txt',LinkList,'Delimiter',' '); % InfoMap code wants .txt , space delim
  
% sort list by wt (weak : strong)
wt=LinkList(:,3); 
[wtsort, wtindex] = sort(wt);
LinkListSortWt=LinkList(wtindex,:); % now sorted by raw wt (incrs);
  % nb 2851 wt entries are 0 ?? ; min ~ 4.59e-6
clear wtsort wtindex wt


%% 1.03 LinkList:  calc from Adj: best way to get log-10(wt) distribution
LinkList=adj2edgeL(Anew); % finds 3474 links 
  % eliminate wt:0 
  %posWts=find(LinkList(:,3)>0); % find the non-zero entries % 3474 of
  % LinkList=LinkList(posWts,:);
  %clear posWts
[nl ~] =size(LinkList) % now 3474 non-0 wts

% nb. Acrn-i, Acrn-j, wt] from .txt file & in LinkListACrn.csv
wt=LinkList(:,3);
wtlog=log10(wt);
figure; hist(wtlog,50); title('Marmoset  log10 (LNe1(i,j)) distribution  ')

% save for InfoMap calc
dlmwrite( 'MarmosetPairListRescale.txt',LinkList,'Delimiter',' '); % InfoMap code wants .txt , space delim
  
% sort list by wt (weak : strong)
wt=LinkList(:,3); 
[wtsort, wtindex] = sort(wt);
LinkListSortWt=LinkList(wtindex,:); % now sorted by raw wt (incrs);
  % nb 2851 wt entries are 0 ?? ; min ~ 4.59e-6
clear wtsort wtindex wt

% match Acrn
ListAcrn=cell(3474,2);
for i=1:3474
    ListAcrn{i,1}=Acrn{LinkListSortWt(i,1)}; ListAcrn{i,2}=Acrn{LinkListSortWt(i,2)};
end


%%  1.8.3 Calc Dist; & wt vs Dist
 % NodeCoord read in at #0,03 above
 NodeCoord=csvread('CoordsMarmoset.csv'); % nb. some absent: NaN 
fprintf('\n Calc Dist(i-j) from node coords \n')
DistCol=zeros(length(LinkList),1);
for i=1:nl
    dist=( (NodeCoord(LinkList(i,1),1)- NodeCoord(LinkList(i,2),1) )^2 ...
        + (NodeCoord(LinkList(i,1),2)- NodeCoord(LinkList(i,2),2) )^2 ...
        + (NodeCoord(LinkList(i,1),3)- NodeCoord(LinkList(i,2), 3) )^2 );
    DistCol(i)=sqrt(dist);
end
clear i dist

LinkList=adj2edgeL(Anew); % finds 3474 links 
[nl ~] =size(LinkList) % now 3474 non-0 wts

LinkListDist=[LinkList DistCol];
clear DistCol LinkList

% sort list by wt
wt=LinkListDist(:,3); 
[wtsort, wtindex] = sort(wt);
LinkListSortWt=LinkListDist(wtindex,:); % now sorted by raw wt (incrs);
  % nb 2851 wt entries are 0 ?? ; min ~ 4.59e-6
clear wtsort wtindex wt

 % LinkList92=LinkListSortWt(284:end,:); %  top 92% : wt> 5e-5  
% other 
DistCol=LinkListSortWt(:,4);
wts=LinkListSortWt(:,3);
wtslog=log10(wts);
mean(wtslog) % now -2.80
 
% Next: sort list by dist [for curve fits]
DistCol=LinkListDist(:,4);
[wtsort, wtindex] = sort(DistCol);
LinkListSortDist=LinkListDist(wtindex,:); % now sorted by d(i-j) (incrs);
 % max Dist is 22.65mm; nb #5941:6325 are NaN :no Injn ??
clear wtsort wtindex DistCol
% now in dist-sorted order
DistCol=LinkListSortDist(:,4);
wts=LinkListSortDist(:,3);
wtslog=log10(wts);
 mean(wtslog) % now -7.590
  % length(find(wtslog < -100) ) % man wt: 0 > -Inf for log
  % log10(mean(wts))  % -2.79 : too many 0!!
  % wtsnz=wts(wts>0); log10(mean(wtsnz)) % 3474 non-zero elements; mean -1.80 

% eliminate wt:0 % Dist: nan
posWts=find(LinkListDist(:,3)>0); % find the non-zero entries % 3474 of
LinkListDist=LinkListDist(posWts,:);
clear posWts
goodDist=find(LinkListDist(:,4) < 30); % 3223 of 
  % max(LinkListDist(:,4)) % 21.04 mm
LinkListDist=LinkListDist(goodDist,:);
clear goodDist

figure; plot(DistCol, wtslog, '.', 'MarkerSize', 14); 
  % lin fit = 2.1 -0.095 * dist
  % save for InfoMap calc
  %dlmwrite( 'MarmosetPairList.txt',LinkList,'Delimiter',' '); % InfoMap code wants .txt , space delim
  
   %LinkListSmall=LinkListDist(:,1:3);
   %dlmwrite( 'MarmosetSmallPairList.txt',LinkListSmall,'Delimiter',' ');
%
% exp fit
 wtsln=log(wts); Distln = log(DistCol);
figure; plot(DistCol, wtsln, '.', 'MarkerSize', 14); 
 % lin fit = 4.8 - 0.22 x 
 % fitLogWt = fit(DistCol, wtsln, 'poly1')
 % cftool(DistCol, wtsln) % app interface  : R^2 = 0.16, big scatter
 %% log - log fit
 cftool(Distln, wtsln)

%% 1.8.4a Plot LinkWts vs Dist / sort dist : for V2b data, rescaled : Fig S7
 % code from MBplotWtsDist.m  % MarmosetReadDataV1.m
 
figure; % figure('position',figposition, 'units','pixels'); hold on; % big
plot(DistCol, wtslog, '.', 'MarkerSize', 18); % 16 for big plot, 
 title('Marmoset cortex, log10 (LNe) vs Link Dist (mm) ','Fontname','Times New Roman','FontSize',12,'FontWeight','Bold') % V2b
%title('Marmoset Brain, link weight vs dist ','Fontname','Times New Roman','FontSize',14,'FontWeight','Bold')
xlabel('dist (mm)','Fontname','Times','FontSize',12,'FontWeight','Bold') 
ylabel('log10 (LNe1 wt) ','Fontname','Times New Roman','FontSize',12,'FontWeight','Bold');
 
% highlight weakest wt (< 7.9e-11) - for SortByWt list !!
   % ie. log10(9.18e-6): -10.0
   % find(wtslog < -10.0)' : 34 of: lines 944, 969... 
   % find(wtslog < -10.1)' : 13 of: lines 944, 969... 
lowList=find(wtslog < -10.1);
for i=1:length(lowList)  % nb the 17 pts in wt-ordered list; need to find in dist-sorted list
       %plot(DistCol(i), wtslog(i), '.', 'Color', 'red','MarkerSize', 11)
       plot(DistCol(lowList(i)), wtslog(lowList(i)), 'o', 'Color', 'red','MarkerSize', 13)
end
clear lowList
 % xxxx linear fit, from menu: logWt = -1.924 - 0.0977*dist (ResNorm = 55.7) V2a
  % logwt:  2.072 - 0.0952 * dist(mm)  (ResNorm = 57.013) V2b : correct  
 
 hold on; legend off
% trend lines
plot([0 25], [2.072 -0.3080], 'k-', 'LineWidth', 1.5) % along the mean of the data
     % lin fit is logwt:  2.072 - 0.0952*dist
 plot([0 25], [(2.072 -2.9022) (-0.3080 -2.9072)], '--r') % -3*sigma line  % 3*sig = -2.9022
 plot([0 25], [(2.072 +2.9022) (-0.3080 +2.9072)], '--r') % +3*sigma line  % 3*sig = 2.9022
plot([0 25], [(2.072  -0.9674) (-0.3080 -0.9674)], '--r', 'LineWidth', 1.2) % -sigma line  % 1*sig = -0.9674
plot([0 25], [(2.072  +0.9674) (-0.3080 +0.9674)], '--r', 'LineWidth', 1.2) % +sigma line  % 1*sig = 0.9674
ylim([-2.5 5])
hold off
% save figure to file  for mss. % fig S7 
hh =gcf; print(hh, 'FigS7.tif', '-dtiff' , '-r300');   
 % nb. again,  all data within 3*sigma boundaries; but some closer now !
 
%% 1.9a form residuals, about log-lin trend: V2b
wtslog0=wtslog-(2.072 - 0.0952*DistCol); % % logwt:  2.072 - 0.0952 * dist(mm)  (ResNorm = 57.013) V2b  
figure('position',figposition, 'units','pixels'); 
plot(DistCol, wtslog0, '.k', 'MarkerSize', 12); hold on
title('Marmoset Brain, log10(wt-LinTrend) vs dist ')

% distribution of residuals
figure; hist(wtslog0, 50)  % ~ log ; symm
title('Marmoset Brain, distribution of  log10(wt-LNe1) residuals ')

mean(wtslog0) % 2.4e-4  % ?? 4.1442
std(wtslog0) % 0.9674

%% 1.8.4b Plot LinkWts vs Dist / sort dist : for V2 data, rescaled
 % code from MBplotWtsDist.m % MarmosetReadDataV1.m
 % read DistCol at #1.8.3 -via Node Coords; & wtslog
 % figure; 
figure('position',figposition, 'units','pixels'); hold on; % big
plot(DistCol, wtslog, '.', 'MarkerSize', 11); 
 title('Marmoset Brain, V2  log10(LNe1) vs Link Dist(mm) ','Fontname','Times New Roman','FontSize',14,'FontWeight','Bold') % V1
%title('Marmoset Brain, link weight vs dist ','Fontname','Times New Roman','FontSize',14,'FontWeight','Bold')
xlabel('dist(mm)','Fontname','Times','FontSize',14,'FontWeight','Bold') 
ylabel('log10(LNe1 wt) ','Fontname','Times New Roman','FontSize',14,'FontWeight','Bold');
 
% highlight weakest wt (< 7.9e-11) - for SortByWt list !!
   % ie. log10(9.18e-6): -10.0
   % find(wtslog < -10.0)' : 34 of: lines 944, 969... 
   % find(wtslog < -10.1)' : 13 of: lines 944, 969... 
lowList=find(wtslog < -10.1);
for i=1:length(lowList)  % nb the 17 pts in wt-ordered list; need to find in dist-sorted list
       %plot(DistCol(i), wtslog(i), '.', 'Color', 'red','MarkerSize', 11)
       plot(DistCol(lowList(i)), wtslog(lowList(i)), 'o', 'Color', 'red','MarkerSize', 8)
end
clear lowList
 % linear fit, from menu: logWt = -1.924 - 0.0977*dist (ResNorm = 55.7)
 legend off
 
% linear fit, to logwt:  -6.721 - 0.0973*dist ; mean 0 about that line; std : xxx1.312
%plot([0 25], [-6.721 -9.154], 'k-', 'LineWidth', 2, 'LineSmoothing', 'on') % along the mean of the data

%plot([0 25], [(-6.721 -3.227) (-9.154 -3.227)], '--r') % -3*sigma line  % 3*sig = 3.227
%plot([0 25], [(-6.721 +3.227) (-9.154 +3.227)], '--r') % +3*sigma line  % 3*sig = 3.227

hold off
 % nb. again,  all data within 3*sigma boundaries; but some closer now !
 
%% 1.9b form residuals, about log-lin trend: V2a
wtslog0=wtslog-(-1.924 - 0.0977*DistCol);
figure('position',figposition, 'units','pixels'); 
plot(DistCol, wtslog0, '.k', 'MarkerSize', 5); hold on
title('Marmoset Brain, log10(wt-LinTrend) vs dist ')

% distribution of residuals
figure; hist(wtslog0, 50)  % ~ log ; symm
title('Marmoset Brain, distribution of  log10(wt-LNe1) residuals ')

mean(wtslog0) % -4.7929
std(wtslog0) % 1.0756

%% 2.0 Examine Hubs: links dist (for Version2, rescaled)
fprintf('\n examine Hubs: V2  \n')
% Hubs-Out, R6: TE3 (96);          R5: A8b(47), AuML (57), V2(107).
 % Hubs-In: R6: A9 (#48), TEO(97);  R5: AuA1(52), A6DC(39).
thisNode=107;
find(Anew(thisNode,:))

listOut=find(LinkListSortDist(:,1) == thisNode); % find rows in LinkLiost
listIn=find(LinkListSortDist(:,2) == thisNode);

%figure; hist(LinkListSortDist(listIn,4), 50); title('Marmoset out-Hub AuML: InLink dist distribution  ')
%figure; hist(LinkListSortDist(listOut,4), 50); title('Marmoset out-Hub AuML: OutLink dist distribution  ')

figure; hist(LinkListSortDist(listIn,4), 50); title('Marmoset in-Hub V2: InLink dist distribution  ')
figure; hist(LinkListSortDist(listOut,4), 50); title('Marmoset in-Hub V2: OutLink dist distribution  ')
figure; hist(LinkListSortDist(listIn,3), 50); title('Marmoset in-Hub V2: InLink log-wt distribution  ')
figure; hist(LinkListSortDist(listOut,3), 50); title('Marmoset in-Hub V2: OutLink log-wt distribution  ')

%% wt-dist : for V2
figure; plot(DistCol(listIn), wtslog(listIn), '.', 'MarkerSize', 14); 
  title('Marmoset in-Hub V2: InLink log-wt vs dist  '); axis([2 25 -10 -4]) % need common
    % linear fit  logwt = -5.30 - 0.140*dist (IN)
    
  figure; plot(DistCol(listOut), wtslog(listOut), '.', 'MarkerSize', 14); 
  title('Marmoset in-Hub V2: OutLink log-wt vs dist  '); axis([2 25 -10 -4]) % nb. stronger link
   % linear fit  logwt = -5.714 - 0.143*dist  { OUT  similar slope

% common plot + fit for V2 
  figure; plot(DistCol(listIn), wtslog(listIn), '.', 'MarkerSize', 14); hold on 
  title('Marmoset in-Hub V2: InLink log-wt vs dist  '); axis([2 25 -10 -4])
  plot(DistCol(listOut), wtslog(listOut), 'or', 'MarkerSize', 7); 
  legend('In links', 'Out links')
   plot([2 25], [-5.58 -8.80], 'b-', 'LineWidth', 1, 'LineSmoothing', 'on') % lin fit, In 
   plot([2 25], [-6.0 -9.289], '--r') % lin fit, Out
 
%% wt-dist : for AuA1
figure; plot(DistCol(listIn), wtslog(listIn), '.', 'MarkerSize', 11); 
  title('Marmoset in-Hub AuA1: InLink log-wt vs dist  '); axis([2 20 -10 -5]) % need common
    % linear fit  logwt = -4.94 - 0.323*dist
  figure; plot(DistCol(listOut), wtslog(listOut), '.', 'MarkerSize', 11); 
  title('Marmoset in-Hub Auå1: OutLink log-wt vs dist  '); axis([2 20 -10 -5])
   % linear fit  logwt = -6.575 - 0.141*dist  { 50% steeper!

% common plot + fit for A9 
  figure; plot(DistCol(listIn), wtslog(listIn), '.', 'MarkerSize', 11); hold on 
  title('Marmoset in-Hub AuA1: InLink log-wt vs dist  '); axis([2 20 -10 -5])
  plot(DistCol(listOut), wtslog(listOut), 'or', 'MarkerSize', 7); 
  legend('In links', 'Out links')
  plot([2 6], [-5.586 -6.878], 'b-', 'LineWidth', 1, 'LineSmoothing', 'on') % lin fit, In 
  plot([2 16], [-6.8570 -8.831], '--r') % lin fit, Out
  
%% wt-dist : for A9
figure; plot(DistCol(listIn), wtslog(listIn), '.', 'MarkerSize', 11); 
  title('Marmoset in-Hub A9: InLink log-wt vs dist  '); axis([2 20 -10 -5]) % need common
    % linear fit  logwt = -5.825 - 0.068*dist
  figure; plot(DistCol(listOut), wtslog(listOut), '.', 'MarkerSize', 11); 
  title('Marmoset in-Hub A9: OutLink log-wt vs dist  '); axis([2 20 -10 -5])
   % linear fit  logwt = -6.373 - 0.149*dist  { 50% steeper!

% common plot + fit for A9 
  figure; plot(DistCol(listIn), wtslog(listIn), '.', 'MarkerSize', 11); hold on 
  title('Marmoset in-Hub A9: InLink log-wt vs dist  '); axis([2 20 -10 -5])
  plot(DistCol(listOut), wtslog(listOut), 'or', 'MarkerSize', 7); 
  legend('In links', 'Out links')
  plot([2 20], [-5.961 -7.185], 'b-', 'LineWidth', 1, 'LineSmoothing', 'on') % lin fit, In 
  plot([2 20], [-6.671 -9.353], '--r') % lin fit, Out
  
%% wt-dist : for TEO
figure; plot(DistCol(listIn), wtslog(listIn), '.', 'MarkerSize', 11); 
  title('Marmoset in-Hub TEO: InLink log-wt vs dist  '); axis([2 20 -10 -5]) % need common
    % linear fit  logwt = -6.2 - 0.15*dist
  figure; plot(DistCol(listOut), wtslog(listOut), '.', 'MarkerSize', 11); 
  title('Marmoset in-Hub TEO: OutLink log-wt vs dist  '); axis([2 20 -10 -5])
   % linear fit  logwt = -6.3 - 0.10*dist  { 50% steeper!

   % common plot + fit for TEO 
  figure; plot(DistCol(listIn), wtslog(listIn), '.', 'MarkerSize', 11); hold on 
  title('Marmoset in-Hub TEO: InLink log-wt vs dist  '); axis([2 20 -10 -5])
  plot(DistCol(listOut), wtslog(listOut), 'or', 'MarkerSize', 7); 
  legend('In links', 'Out links')
  plot([2 20], [-6.5 -9.2], 'b-', 'LineWidth', 1, 'LineSmoothing', 'on') % lin fit, In 
  plot([2 20], [-6.5 -8.3], '--r') % lin fit, Out
  
  
%% 2.1 Examine non hub Connectors: links dist (for Version2, rescaled)
% non hub Connectors-Out, R3: A8aV(46), A23b(14), & A47L also
 % non hub Connectors:    R3: A47L(34), A8C(44)
fprintf('\n examine no-hub Connectors (R3): A47L  \n')
thisNode=34;
find(Anew(thisNode,:))

listOut=find(LinkListSortDist(:,1) == thisNode); % find rows in LinkLiost
listIn=find(LinkListSortDist(:,2) == thisNode);

figure; hist(LinkListSortDist(listIn,4), 50); title('Marmoset in-Hub A47L: InLink dist distribution  ')
figure; hist(LinkListSortDist(listOut,4), 50); title('Marmoset in-Hub A47L: OutLink dist distribution  ')

figure; hist(LinkListSortDist(listIn,3), 50); title('Marmoset in-Hub A47L: InLink log-wt distribution  ')
figure; hist(LinkListSortDist(listOut,3), 50); title('Marmoset in-Hub A47L: OutLink log-wt distribution  ')

%% 2.1a wt-dist : for A47L
figure; plot(DistCol(listIn), wtslog(listIn), '.', 'MarkerSize', 11); 
  title('Marmoset in-Hub A47L: InLink log-wt vs dist  '); axis([2 20 -10 -5]) % need common
    % linear fit  logwt = -7.372 - 0.0220*dist || need to do this for each node!
  figure; plot(DistCol(listOut), wtslog(listOut), '.', 'MarkerSize', 11); 
  title('Marmoset in-Hub A47L: OutLink log-wt vs dist  '); axis([2 20 -10 -5])
   % linear fit  logwt = -6.856 - 0.0280*dist  { 50% steeper!
   
   
%% Appx. X.0 Clustering coeff
 % fraction of transitive triples
  % BCT (Rubinov) code needs A in [0,1]
   %maxWt=max(A(:))
   %A1=A./maxWt;
   %max(A1(:)) % check
   %C = clustering_coef_wd(A1); % clustering coefficient vector
   %figure; plot(C, '.')
   %Call=mean(C)  % av over all nodes

% Clustering coeff, Altv.  Buonova (Newman / MIT) code
[C1,C2,C] = clust_coeff(Anew); % C of In-, out-, Av over all Links
figure; plot(C, '.')
xlabel('node ID '); ylabel('Clustering Coeff ');
title('Marmoset Brain (rescaled wt), Clustering Coefficient ')
  figure; hist(C,50)
  title('Marmoset Brain, Clustering Coefficient distribution ')
   %axis([0.1 3.5 0 10])  % nb. many at "0"
   
%% X.1 Rich club calc. for binary, dir Adj
% needs A1, un-wt
% DegIn1=[DegIn1; 0; 0];; % missing 2 x 0 at end
figure; hist(DegIn1) % check hi-k nodes
find (DegIn1>80)' % 16 nodes; 9 are > 90
figure; hist(DegOut1)  % nb. about half of k-In
find (DegOut1>40)' % 16 nodes;
% uses:  degree=DegIn1+DegOut1';
%        find (degree>120)  % 15 of 
% un-wt
%[R,Nk,Ek] = rich_club_bd(A1); % scans DegIn+DegOut by default
% [R,Nk,Ek] = rich_club_bd(A1, 80); % set max k-level
Rw=rich_club_wd(Anew); % wt-dir network
figure; stem(Nk)
figure; stem(Ek)
figure; stem(R)

figure; semilogy(Nk, R, '.')

Nkud=flipud(Nk); % order by incrs tot deg
figure; loglog(Nkud, R, '.')
axis([1 1000 0.01 1])
title('Marmoset: rich club spectrum  ')


  