% B Pailthorpe, U. Sydney.  V1 (27/9/19-24)
 % Marmoset brain: Load data: connectivity; Coords;
% Dependencies:  data files: AdjMarmoset.csv, AdjFullMarmoset.csv, AdjMarmoset_Rescale2b.csv, 
      % MarmosetAcrn116.csv, MarmosetAcrn139to116.csv, CoordsMarmoset.csv, MarmosetLobes116.csv,
      % VolsMarmoset.csv, MarmosetPairListRescale2b.clu % InfoMap (c code) module ID list (.clu file).
      % Marmoset_NodeFlow, Marmoset_FlowWt % Prob'y flow over Nodes & Links (.csv):
        % for original or extra details & 3D vis.:
        % Maybe: MarmosetGrayImageSlice.csv % saved image volume data, midslice,
               % from Atlas volume (atlas.nii); nb. large files (86MB) - cf Appx.
               % atlas_segmentation.nii (26MB);  %  needs NIFTI tools (Matlab), load_nii, etc. 
      
% Do these things:
% #1,  import link data & other info: Acrn;  ??Coords (Atlas,  
%         & midslice of vol image
%       Flows, nodes & links
%       Infomap module ID, in array idx
%  1.8  d(i,j) dist matrix, from coords
%       Link List ; top 40% (line 322)
%  1.9 Investigate node vol effects on FLNe wts  (line 406)
%  & 1.3a #Links: un-wt Deg (line 525 - out of order)

%   1.8.2 Calc LinkList of pos Wts & plot hist(log10)

%  2.0  LNe,1 (normalised to 1 mm^3 target vol) data (line 617 )
%  2.1  calc LN { * denom : tot#LN - LN,i } for ea target
%  2.2  calc LNe : new adj & wts   (Line (778

%  3.0 "Symmetrise" the 116 x 55 Adj (line 692 +)

% 4.0 examine spread of Injection Volumes [line 741]
% 4.1 re-examine LN data (3 exp measurements): Denom ~ 0.2% variation:  OK

% Apendices:  % A.0 get RGB colours for areas: from Atlas files (line 958)
               %  & test & devmt codes
%   A.1  get Atlas files  { line 983
%         calc  Node volumes, CM, from bvol 3D data  
%   A.1.1 3D vol of areas ID# (in atlas_labels.txt)   {line
%   A.1.2a:  mid slice grayscale (mm scale) from image data (line ; And 1193 [clean code])
 
% >>>>
%% 1.0 Read data files & setup
  % may use matlab toolboxes for NIFTI files etc - in sub-paths
  % nb. simpler path on this comp - vs laptop
  % nb. read from InfoMap flow calc in code: PlotInfoMapV2.m 
  %     & calc flows over i-j links in MB_flows.m (at #2.)
clear all; %close all
addpath(genpath('/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain')); % Data files;  include sub-directories 
%addpath('/bap_working/MatLabfiles/MarmosetBrain/CodesOther'); % for other codes & dependencies
%
fprintf('\n >> Vis Marmoset Brain: : Load (raw) linklist data \n')
Actx=csvread('AdjMarmoset.csv');  % read square A(55 x 55); selection of raw A (both S&T); no self-loops
nsquare = size(Actx,1)  % number of nodes: Sources AND Targets

Afull= csvread('AdjFullMarmoset.csv'); % calc prev. (116 x 116)  % read raw A (all S, some T); no self-loops
[n, ~] = size(Afull)  % number of Sources, 

% set one of these Adj to A, for calc below:
figure; spyc(Actx) % 1854 links:
title('Marmoset Adj: Ctx-Ctx: 55 x 55 ')
figure; spyc(Afull) % 3474 links
title('Marmoset Adj: Ctx-Ctx: 116 x 55 ')

%plotArray(A)    % coarser grid, larger squares, but only 1379 squares?

% set up larger figs for later
figwidth = 1024; figheight = 896; figposition = [100, 100, figwidth, figheight];
 % A = Afull; % easier for other codes 
 % clear Afull
 
%% 1.01 fully connected 55x55 subnetwork (FLNe)
% count links (<~ 50-150)
n55Links = length(find(Actx(:)))
frn55 = n55Links/(55*54)
[sii sjj v]=find(Actx);  
Alogical=full(sparse(sii, sjj, 1));    % "mask" of selected entries, to re-form sq array of 1's.
 % nb. need to use sparse to form the array
A1 = Alogical.*1;  % pick out the selected entries; yields a "full" array
figure; spy(A1) % debug
clear Alogical sii sjj v
length(find(A1(:))) % check 3474 links: 
  clear sii sjj v n55Links frn55
% length(find(A1(:)) )
% count links (<~ 50-150)
DegIn1=sum(A1,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
DegOut1=sum(A1,2); % i->j; row sum 
figure; hist([DegIn1,DegOut1], 50); legend('k-In', 'k-Out')
xlabel('Degree (link count)'); ylabel('counts')
title('Marmoset, Degree of 55x55 binarised links')
% wt. degree (node strength)
DegIn55=sum(Actx,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
DegOut55=sum(Actx,2); % i->j; row sum 
figure; hist([DegIn55,DegOut55], 50); legend('wtDeg-In', 'wtDeg-Out')
xlabel('weight (FLNe)'); ylabel('counts')
title('Marmoset Deg of 55x55 frn. wt. links')
[max(DegIn55) max(DegOut55) ]
 % nb. scale doesnt make sense, since FLN is fraction amongst 116 sources (not 55).

 %% 1.1d wt Degree (FLNe)
 % unwt deg at 1.3a - far below
DegIn=sum(Afull,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
DegOut=sum(Afull,2); % i->j; row sum 
 % figure; hist(DegIn); title('Marmoset wt DegIn ')
 %figure; hist(DegOut); hold on; title('Marmoset wt DegOut '); hold off

hFig=figure; hist([DegIn, DegOut]);
legend('DegIn', 'DegOut'); title(' Marmoset brain, FLN,e: Deg-In,-Out ')

 % cf un-wt Deg (#Links) at  #1.3a, line 480
 kTot = length(find(Afull(:)))
 
 clear n55Links frn55
 
 %% 1.1a read V2b ADj, based on LNe
Anew= csvread('AdjMarmoset_Rescale2b.csv'); % read re-scaled Adj (116 x 116)
[n, ~] = size(Anew)  % number of Sources
 
%% 1.1 read Node Acrn, Labels & Coord (CM)
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

NodeID=csvread('MarmosetAcrn139to116.csv'); % for the 116 Nodes (cf Links) vs 136 Labels
 % translates the full 139 List (Atlas) to the 116 nodes (Source Injn)

 % 1.1b  coords - prev calc from Atlas volumes; cf Appx. A.1 below
        % from 3D vol: 'atlas_segmentation.nii'
NodeCoord=csvread('CoordsMarmoset.csv'); % nb. some absent: NaN 
% adjust z-origin, to match Paxinos atlas??
 % NodeCoord(:,3)=NodeCoord(:,3)-20; % origin at V (top) & points down
max(NodeCoord(:, 3)) % z

NodeVols=csvread('VolsMarmoset.csv'); % (mm^3 - calc from bv.img[.nii]))
figure; hist(NodeVols,50); title('Marmoset: Node Vol Distrn (mm^3) ')

% NodeVoxels=csvread('VoxelsMarmoset.csv'); % get file to laptop

%% 1.1c read the 116 Lobes assigned (in Paxinos Atlas)
 % reader skips blanks, so load 'na' to preserve numbering
filename = '/bap_working/MatLabfiles/MatlabFiles/MarmosetBrain/MarmosetLobes116.csv';
delimiter = ',';  fileID = fopen(filename,'r');  formatSpec = '%s%[^\n\r]';
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false); fclose(fileID); % these are Cells, containing strings
Lobes = dataArray{1, 1}(1:end); % these are cells; 
clearvars filename delimiter formatSpec fileID dataArray ans
Lobes{108} % check this is 'Vis': ok

%LobesTxt=char(Lobes); % use char array to search the 116 Lobe names - faster

%% 1.1d wt Degree (FLNe)
 % unwt deg at 1.3a - far below
DegIn=sum(Afull,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
DegOut=sum(Afull,2); % i->j; row sum 
 % figure; hist(DegIn); title('Marmoset wt DegIn ')
 %figure; hist(DegOut); hold on; title('Marmoset wt DegOut '); hold off

hFig=figure; hist([DegIn, DegOut]);
legend('DegIn', 'DegOut'); title(' Marmoset brain, FLN,e: Deg-In,-Out ')

 % cf un-wt Deg (#Links) at  #1.3a, line 480
 kIn = length(find(Afull(:)))
 

 %% 1.1d.1 wt Degree (LNe)
DegIn=sum(Anew,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
DegOut=sum(Anew,2); % i->j; row sum 
 % figure; hist(DegIn); title('Marmoset wt DegIn ')
 %figure; hist(DegOut); hold on; title('Marmoset wt DegOut '); hold off

hFig=figure; hist([DegIn, DegOut]);
legend('DegIn', 'DegOut'); title(' Marmoset brain, (LN,e): Deg-In,-Out ')

% #Links: un-wt Deg -  binarise 116x116 A:  UN-wt links
[sii sjj v]=find(Anew);  
Alogical=full(sparse(sii, sjj, 1));    % "mask" of selected entries, to re-form sq array of 1's.
 % nb. need to use sparse to form the array
A1 = Alogical.*1;  % pick out the selected entries; yields a "full" array
figure; spy(A1)
clear Alogical sii sjj v
length(find(A1(:))) % check 3474 links:
% count links (<~ 50-150)
DegIn1=sum(A1,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
DegIn1=[DegIn1; 0; 0]; % append 2 missing 0's
DegOut1=sum(A1,2); % i->j; row sum 
figure; hist(DegIn1, 50)
title('Marmoset DegIn of binarised links')
axis([0 110 0 10]) % nb. big # at 0: not sampled
%
figure; hist(DegOut1, 50); hold on
title('Marmoset DegOut of binarised links')
 %axis([0 80 0 130]); 
figure; stem (DegIn1, '.'); hold on
stem (DegOut1, 'r.');
legend('DegIn1', 'DegOut1')
title('Marmoset Deg of binarised links')
% joint plot
figure; hist([DegIn1, DegOut1], 50); xlim([2 Inf]); hold on
title('Marmoset DegIn, Out of binarised links')
legend('DegIn1', 'DegOut1')
%% node vol dependence of Deg, k
figure; subplot(2,2,1); plot(NodeVols, DegIn, '.', 'MarkerSize', 12); title('Marmoset: wt-DegIn')
subplot(2,2,2); plot(NodeVols, DegIn1, '.', 'MarkerSize', 12); title(' #Links, k-In')
subplot(2,2,3); plot(NodeVols, DegOut, '.', 'MarkerSize', 12); title(' wt-DegOut')
subplot(2,2,4); plot(NodeVols, DegOut1, '.', 'MarkerSize', 12); title(' #Links, k-Out')
% & expand scale:
figure; subplot(2,2,1); plot(NodeVols, DegIn, '.', 'MarkerSize', 12); title('Marmoset: wt-DegIn')
 axis([0 40 -Inf Inf])
subplot(2,2,2); plot(NodeVols, DegIn1, '.', 'MarkerSize', 12); title(' #Links, k-In')
axis([0 40 -Inf Inf]); 
subplot(2,2,3); plot(NodeVols, DegOut, '.', 'MarkerSize', 12); title(' wt-DegOut')
axis([0 40 -Inf Inf]); xlabel('Node Vol (mm^3) ')
subplot(2,2,4); plot(NodeVols, DegOut1, '.', 'MarkerSize', 12); title(' #Links, k-Out')
axis([0 40 -Inf Inf]); xlabel('Node Vol (mm^3) ')

% fits
figure; plot(NodeVols, DegIn1, '.', 'MarkerSize', 12); title(' #Links, k-In') 
axis([0 40 -Inf Inf]); % use GUI for linear fit; avoid bias by many 0's
nodes1=find(DegIn1 >0);
DegIn11=DegIn1(nodes1);   % need non-zero,
vols1=NodeVols(nodes1);
figure; plot(vols1, DegIn11, '.', 'MarkerSize', 12); title(' #Links, k-In') 
 axis([0 40 -Inf Inf]); %
 [mean(DegIn11) std(DegIn11)] % stsats of non-zero k
     
clear nodes1 DegIn11 vols1

% Out links vs 2D perimeter
NodeRad=(3*NodeVols/(4*pi)); NodePerim =2*pi*NodeRad.^2; %  assume cyl regions
figure; subplot(1,2,1); plot(NodePerim, DegIn1, '.', 'MarkerSize', 14); title(' #Marmoset: Links, k-In')
xlabel('Node perimeter (mm) '); axis([0 70 -Inf Inf]) % omit 2 big outliers
subplot(1,2,2); plot(NodePerim, DegOut1, '.', 'MarkerSize', 14); title(' #Links, k-Out')
xlabel('Node perimeter (mm) '); axis([0 70 -Inf Inf])

clear NodeRad NodePerim

 %% 1.1e Strongest links
 
 length(find(Afull(:)>= 0.1)) %  strongest Out links
 length(find(Afull(:)>= 0.4)) %  V strongest Out links 

 %% 1.1g  Deg, k - vol sampling  (LNe)

 figure; subplot(2,2,1); plot(NodeVols, 
 
 %% ++++>>> 1.4.1 NEEDED:  set mid-transverse(x,y)-plane slice : as reference plane in 3D plots (mm space)
% needs the Atlas volume (.nii) - cf Appx  , below
% read saved image volume data, midslice
gray_imageflip=csvread('MarmosetGrayImageSlice.csv'); % prev calc RGB; x-y flipped
gray_imageflip=uint8(gray_imageflip); % Need integer; nb image was x-y flipped
gray_imageflipR=fliplr(gray_imageflip); % reflect LHS to get RHS image
% set x-y box & desired z position of the image plane.
zpix=350;  pixsize=[0.0400 0.5000 0.0400]; % set manually, avoid reloading large image
min_x = 1; max_x = 825; min_y = 1; max_y = 63; 
min_z = 1; max_z = 550; imgzposition = zpix; % "half" height, pixel space
%  pix-mm conversion, % apply vox -> mm tranform 
origin = [411,24, 57]; % from bv.hdr.hist.originator 
pixsize=[0.0400 0.5000 0.0400]; % set manually, avoid reloading large image
min_xmm=(min_x-origin(1))*pixsize(1); max_xmm=(max_x-origin(1))*pixsize(1); 
min_ymm=(min_y-origin(2))*pixsize(2); max_ymm=(max_y-origin(2))*pixsize(2);
min_zmm=(min_z-origin(3))*pixsize(3); max_zmm=(max_z-origin(3))*pixsize(3);
imgzpositionmm=(imgzposition-origin(3))*pixsize(3);  %  centered, in mm-coords
clear origin pixsize min_x max_x min_y max_y min_z max_z imgzposition zpix
% test plot in mm space the image plane using surf.
figure; hold on; grid on  % nb. need to eliminate "box"
surf([min_xmm max_xmm],[min_ymm max_ymm],repmat(imgzpositionmm, [2 2]),...
    gray_imageflip,'facecolor','texture', 'FaceAlpha', 0.2,'edgecolor', 'none')  % 
colormap(gray); view(-30, 30)
plot3(0, 0, 20, '+k','MarkerSize', 30);
title('Marmoset Brain:  Mid Slice z= -5.72 mm ')
ylabel('y (mm)'); xlabel('x (mm)'); % nb. swap
axis equal; axis([-10 10 -10 20 0 20])
% add RHS - symmetric about mid-line
surf([min_xmm max_xmm],[min_ymm max_ymm],repmat(imgzpositionmm, [2 2]),...
    gray_imageflipR,'facecolor','texture', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.0) 
hold off


%% 1.7  ++NEEDED:  Read Module IDs - InfoMap calc; Assign colours (ipsi-lateral modules)
% use A; etc, needs LinkList - calc from A116x52, A52x52, Aunwt52x52, etc
% needs idx: InfoMap module ID list = output of InfoMap calc, in .clu file 
  %idx=dlmread('MarmosetPairListRescale.clu'); % V2 of data: rescale by LNe/mm^3 
%idx=dlmread('MarmosetPairList.clu');  % wt >=1; V1, based on FLNe
idx=dlmread('MarmosetPairListRescale2b.clu');
  % idx=dlmread('MarmosetPairList55.clu');  %  55 x 55 Ctx only
  % aggregate few orphans, based on linkages, tbd; 7 modules #8-14 orphans
% idx(find(idx>=7))=7; % orphans  
idx(find(idx>=8))=8; % orphans: "last" module[s] {might be to #1 - majority of in/out links
  %idx(51)=8; % reassign Apir as orphan [has no detected links in or out]
% & set Node, Module colors
 %nm=max(idx) % # modules to use (from InfoMap)
fprintf(' >>>>  use idx:  Infomap of links,  Assign Colors    \n')
%  %fprintf(' >>>>  use idx:  Infomap of links, bottom 10pc,  Assign Colors    \n')
nodecolors=zeros(n,3); modulecolors=zeros(8,3);
% ** update these (9/1/18) - to match modules members to Allen Atlas colours
 % cf.   http://atlas.brain-map.org/ for assigned colors
% >>> need to revise colours to match marmoset atlas
for i = 1:n
        if idx(i) == 1 
            %nodecolors(i,:)=[1, 102/255, 216/255]; % for mouse: bright pink:  MidBrain;
           nodecolors(i,:)=[0.9412	1.0000	0.2314]; %use average of Atlas Node Colors
            %modulecolors(1,:)=[1, 102/255, 216/255];% for mouse
           modulecolors(1,:)= [0.9412	1.0000	0.2314]; % adjust for more contrast
                  %[0.9294	0.818	0.2839]; %use average of Atlas Node Colors
        elseif idx(i) == 2 
             %nodecolors(i,:)  =[2/255, 162/255, 88/255];   % mid-green: Vis
            nodecolors(i,:)  =[0.7037, 0.9495, 0.2913];  %use average of Atlas Node Colors
             %modulecolors(2,:)=[2/255, 162/255, 88/255];  % for mouse
            modulecolors(2,:)=[0.7037, 0.9495, 0.2913];  %use average of Atlas Node Colors
        elseif idx(i) == 3
             %nodecolors(i,:)  = [250/255, 20/255, 20/255];  % red-ish/pink  % DlpFC
            nodecolors(i,:)  = [1.0000	0.3255	0.0000]; % adjust for more contrast
                %[0.9001	0.505	0.2539]; %use av
             %modulecolors(3,:)= [250/255, 20/255, 20/255];
            modulecolors(3,:)= [1.0000	0.3255	0.0000]; % use average of Atlas Node Colors
        elseif idx(i) == 4
              %nodecolors(i,:)  =[121/255, 213/255, 190/255];   % :  Tan 
             nodecolors(i,:)  =[0.9592	0.4400	0.2694]; %use av
              %modulecolors(4,:)=[0, 213/255, 190/255];
             modulecolors(4,:)=[0.9592	0.4400	0.2694]; %use av
        elseif idx(i) == 5 
             %nodecolors(i,:) = [0 100/255 137/255];   % : brownish  : VL & DL pFC
            nodecolors(i,:) = [0.9131	0.7196	0.4680]; % use av
             %modulecolors(5,:)=[0 100/255 137/255];
            modulecolors(5,:)=[0.9131	0.7196	0.4680]; % use av
        elseif idx(i) == 6   
             %nodecolors(i,:)  =[30/255, 240/255, 60/255];   % Yellow: Cing
            nodecolors(i,:)  =[0.9000	0.8961	0.2319]; % use av
             % modulecolors(6,:)=[30/255, 240/255, 60/255];
            modulecolors(6,:)=[0.9000	0.8961	0.2319]; % use av
        elseif idx(i) == 7
             %nodecolors(i,:)  =[0.9,112/255, 160/255];  % yellow/orange: LTemp, Insul
            nodecolors(i,:) =[0.9508	0.7109	0.2050]; % use av
             %modulecolors(7,:)=[0.9,112/255, 160/255];
            modulecolors(7,:)=[0.9508	0.7109	0.2050]; % use av
        else   % Gray - (grp 8+ are orpans) 
            nodecolors(i,:)      =[0.5,0.5,0.5];  % others, minor modules :  assign colours later?
            modulecolors(8,:)=[0.5,0.5,0.5];    
        end
end 

%% 1.7.1 Flow:  Prob'y FLOWS over links - InfoMap
% use code PlotInfoMapV2.m to parse the InfoMap .map output file 
% read flow over Nodes & Links:
NodeFlow= csvread('Marmoset_NodeFlow'); 
FlowWt=csvread('Marmoset_FlowWt');

%% 1.7.2 Flow: sort lists 
list=FlowWt(find(FlowWt>0) ); % 1854 links
figure; hist(list, 1090)
listb=find(FlowWt(:)>1e-3); % 261 of
% make & sort list
[ii jj]=ind2sub(n,listb);
flowList=[ii jj zeros(length(ii),1)];
%
for i=1:length(ii)
    flowList(i,3)  =FlowWt(ii(i), jj(i) );
end
%
flowList3=flowList(:,3);
[degsort, degindex] = sort(flowList3); % sorted list & get new index
flowListSort=flowList(degindex,:);
flowListSort=flipud(flowListSort);  % list are largest... smaller
sum(flowListSort(:,3)) % is 0.7669, ie. 77% of all flows
    
 clear degsort degindex flowList3 
%% 1.7.3 Flow:  list top 10, 20, 30 of Deg, Flow
 %NodeFlow= csvread('MBnodeflows.csv');
[degsort, degindex] = sort(NodeFlow); % sorted list & get new index
hifi = degindex(n-19:n);      % list top 20, 30
hifi=flipud(hifi);
clear degsort degindex
 fprintf('\n >> Marmoset Brain: hi (In-) InfoFlow nodes \n')
 
 for ii=1:20 %length(hifi)
     i=hifi(ii);
     fprintf(' %4.0f %s %s %7.4e \n', i, Acrn{i}, Lobes{i}, NodeFlow(i) )
 end
 clear degsort degindex  % use hifi

%% 1.7.3 Calc flows over wt links [or read prev calc - at 1.6.1 above]
fprintf('\n >  Marmoset Brain, flows In: thru Nodes & over wt-links (V2b)  \n')
FlowWt=zeros(n,n); % allocate array
A = Anew; % v2b wts
 % length(find(A(:) ) ) % #links =  3474 ok
for ii=1:n
 %OUT: list of 1st NN
% jj = find(A(ii,:)); %  find only strong Out links; > this cut off (>0
 %jj = find(A(ii,:)>=1); %  find only  Out links
 jj = find(A(:, ii)>=1); %  find only  In links
 wts = A(ii,jj);  % check the i-{nn} Out weights
 TotOutWt=sum(wts); % sum =1 ok
 %wts = round(A(ii,jj))  % check the Out weights
 for j=1:length(jj) % loop over NN of ii
    jjj= jj(j);
    FlowWt(ii,jjj)=  NodeFlow(ii)*A(ii,jjj)/TotOutWt; % partiton out flow to this link
 end  % loop over linked nodes 
end % loop over all nodes
  sum(FlowWt(:) ) % = 1.0 ok
 % csvwrite('Marmoset_FlowWt', FlowWt);
clear jj wts TotOutWt
% checks
[sum(FlowWt(:, 3)) sum(FlowWt(3, :))] % : 0.0172    0.0174 % sum in & out
NodeFlow(3)  % 0.0174  ok 
 find( FlowWt == max(FlowWt(:)) )
 [ii jj] = ind2sub(116, 4794) % is (38, 42)

% sort flows over links:
 flowLinkList=adj2edgeL(FlowWt); % 16864 links (no self loops)
 % csvwrite('MBrawLinkList.csv', rawLinkList); % save file
% sort by wts: weak to strong
wts=flowLinkList(:,3);
[wtsort, wtindex] = sort(wts);
flowLinkListSort=flowLinkList(wtindex,:); % now sorted by raw wt (incrs)
flowListSort=flipud(flowLinkListSort); % Decrs. order

% clear wts wtsort wtindex flowLinkList flowLinkListSort 
%
fprintf('top 20 flows, NodeFlow(i:out)')
for ii=1:20 
     i=flowListSort(ii,1) ; j=flowListSort(ii,2) ;
     fprintf(' %4.0f %s %s %7.4e %4.0f %s %s \n', i, Acrn{i}, Lobes{i}, NodeFlow(i), j, Acrn{j}, Lobes{j} )
end

%
fprintf('\n top 20 flows, i -> j: OUT links \n')
for ii=1:20 
     i=flowListSort(ii,1) ; j=flowListSort(ii,2) ;
     fprintf(' %4.0f %s %s %7.6f %4.0f %s %s \n', i, Acrn{i}, Lobes{i}, flowListSort(ii,3), j, Acrn{j}, Lobes{j} )
end


% form dist matrix

%% 1.8  Calc i-j Distance Matrix - for links that exist
% needs Coords array - read up at #1.1, or here:
NodeCoord=csvread('CoordsMarmoset.csv'); % nb. some absent: NaN 
Adist=csvread( 'MarmosetDistPairs.csv'); % read prev calc 
 % prev calc from bv.img [.nii] below at A.1.2 [line ~400]
fprintf('\n Calc dist array, from calc Atlas Coords (for links >=0) \n')
Adist=zeros(n); % initial;ise array:  d(i,j) in mm - from Atlas area coords
for jjj=1:n
    for kk=1:n
        if Afull(jjj,kk)>0  % for non-zero wts only 
  xj=NodeCoord(jjj,1); yj=NodeCoord(jjj,2); zj=NodeCoord(jjj,3);
  dxjk=abs(NodeCoord(kk,1)-xj); %col vec % nb. need to use actual x-coord, not "idx"
  dyjk=abs(NodeCoord(kk,2)-yj); % 
  dzjk=abs(NodeCoord(kk,3)-zj);
  drjk=sqrt(dxjk*dxjk+ dyjk*dyjk+ dzjk*dzjk); %clear dxs dys dzs
  Adist(jjj,kk)=drjk;
        end
    end
end
clear xj yj zj dxjk dyjk dzjk drjk
length(find(Adist(:)>0)) % check 3223 

%
% csvwrite( 'MarmosetDistPairs.csv',Adist); % save file
%
figure; hist(Adist(:), 100); hold on  % ~ skewedGaussian, peaks at 7, 11, 13 mm 
title('MarmosetBrain: Link distances','Fontname','Times New Roman','FontSize',12,'FontWeight','Bold'); 
axis([0.2 22.5 0 100])
xlabel('link length (mm)','Fontname','Times','FontSize',12,'FontWeight','Bold') 
ylabel('counts','Fontname','Times New Roman','FontSize',12,'FontWeight','Bold'); hold off
  % looks log Normal 
  
%% 1.8.1 Calc LinkList of Wts (FLNe) & plot hist(log10)
fprintf('\n Calc LinkList from rescaled Adj-new (LNe) \n')
 % LinkList=csvread('MarmosetLinkList.csv'); % i-j-wt; 6325 entries, incl many 0's
 % Or, calc from Adj:
 LinkList=adj2edgeL(Anew); % Rescaled; finds 3474 links for 116 x 166 Adj
% eliminate wt:0 
posWts=find(LinkList(:,3)>0); % find the non-zero entries % 3474 of
LinkList=LinkList(posWts,:);
clear posWts
[nl ~] =size(LinkList) % now 3474 non-0 wts

% nb. Acrn-i, Acrn-j, wt] from .txt file & in LinkListACrn.csv
wt=LinkList(:,3);
figure; hist(wt,50); title('Marmoset rescaled LNe distribution  ')
wtlog=log10(wt);
figure; hist(wtlog,50)
  %axis([-5.3 0 0 150]) % omit large peak at "0")
title('Marmoset original log-10 wt (LNe) distribution  ')
xlabel('log10 (weight, LNe)'); ylabel('Counts')
 %  hh =gcf; print(hh, 'FigS3.tif', '-dtiff' , '-r300');
 

% sort list of wt
wts=LinkList(:,3);
[wtsort, wtindex] = sort(wts);
LinkListSort=LinkList(wtindex,:); % now sorted by raw wt (incrs)
clear wts wtsort wtindex

% Top 40% wt>= 0.0033 ~10^-2.48 (#2093-3474: 1382 Links)
%LinkListT40=LinkList(2093:end, :);
 %A40=edgeL2adj(LinkListT40); % no diag, Ok / missing 3 nodes?
 %S=LinkList(:, 1);  S=unique(S); % 116 ok
 %T=LinkList(:, 2);  T=unique(T); % 55, as before - 
 
 % clear LinkList 
%% 1.8.2 Alt.  Calc LinkList of Wts (FLNe) & plot hist(log10)
fprintf('\n Calc LinkList from orig Adj-full \n')
 % LinkList=csvread('MarmosetLinkList.csv'); % i-j-wt; 6325 entries, incl many 0's
 % Or, calc from Adj:
 LinkList=adj2edgeL(Afull); % finds 3474 links for 116 x 166 Adj
% eliminate wt:0 
posWts=find(LinkList(:,3)>0); % find the non-zero entries % 3474 of
LinkList=LinkList(posWts,:);
clear posWts
[nl ~] =size(LinkList) % now 3474 non-0 wts

% nb. Acrn-i, Acrn-j, wt] from .txt file & in LinkListACrn.csv
wt=LinkList(:,3);
figure; hist(wt,50); title('Marmoset original FLNe distribution  ')
wtlog=log10(wt);
figure; hist(wtlog,50)
axis([-5.3 0 0 150]) % omit large peak at "0")
title('Marmoset original log-10 (FLNe) distribution  ')
xlabel('log10 (weight, FLNe)'); ylabel('Counts')
 %  hh =gcf; print(hh, 'FigS3.tif', '-dtiff' , '-r300');
 clear LinkList

% sort list of wt
wts=LinkList(:,3);
[wtsort, wtindex] = sort(wts);
LinkListSort=LinkList(wtindex,:); % now sorted by raw wt (incrs)
clear wts wtsort wtindex

% Top 40% wt>= 0.0033 ~10^-2.48 (#2093-3474: 1382 Links)
%LinkListT40=LinkList(2093:end, :);
 %A40=edgeL2adj(LinkListT40); % no diag, Ok / missing 3 nodes?
 %S=LinkList(:, 1);  S=unique(S); % 116 ok
 %T=LinkList(:, 2);  T=unique(T); % 55, as before - 
 
%% need to do long-hand
A40=zeros(116); % set up array
for i=1:length(LinkListT40)
    A40(LinkListT40(i, 1), LinkListT40(i, 2)) = LinkListT40(i, 3);
end

%%  1.8.3 Calc Dist; & wt vs Dist
% LinkList calc at #
DistCol=zeros(length(LinkList),1);
for i=1:nl
    dist=( (NodeCoord(LinkList(i,1),1)- NodeCoord(LinkList(i,2),1) )^2 ...
        + (NodeCoord(LinkList(i,1),2)- NodeCoord(LinkList(i,2),2) )^2 ...
        + (NodeCoord(LinkList(i,1),3)- NodeCoord(LinkList(i,2), 3) )^2 );
    DistCol(i)=sqrt(dist);
end
clear i dist

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
 
% Next: sort list by dist
DistCol=LinkListDist(:,4);
[wtsort, wtindex] = sort(DistCol);
LinkListSortDist=LinkListDist(wtindex,:); % now sorted by d(i-j) (incrs);
 % max Dist is 22.65mm; nb #5941:6325 are NaN :no Injn ??
clear wtsort wtindex % DistCol {is needed for plot}
% now in dist-sorted order
DistCol=LinkListSortDist(:,4);
wts=LinkListSortDist(:,3);
wtslog=log10(wts);
 mean(wtslog) % now -2.80
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

  % save for InfoMap calc
  %dlmwrite( 'MarmosetPairList.txt',LinkList,'Delimiter',' '); % InfoMap code wants .txt , space delim
  
   %LinkListSmall=LinkListDist(:,1:3);
   %dlmwrite( 'MarmosetSmallPairList.txt',LinkListSmall,'Delimiter',' ');
%
%% 1.9 Investigate node vol effects on FLNe
LinkListVol=zeros(length(LinkList),5);
LinkListVol(:,1) = LinkList(:,1); LinkListVol(:,3) = LinkList(:,2); % Load Source, Target IDs
LinkListVol(:,5) = LinkList(:,3); % and wts
for i=1:length(LinkList)
    LinkListVol(:,2) = NodeVols(LinkList(:,1)); % Source vol
    LinkListVol(:,4) = NodeVols(LinkList(:,2)); % Target vol 
end
%
wtslog=log10(LinkListVol(:,5));
%volslog=log10(LinkListVol(:,2)) %+log10(LinkListVol(:,4)); % log(Vs*Vt)
volslog=log10(LinkListVol(:,4));
figure; hold on; %plot(LinkListVol(:,2), LinkListVol(:,5), '.')
 %plot(volslog, wtslog, '.') % log-lin plot
 %plot(LinkListVol(:,4), wtslog, '.') % log-lin plot
plot(volslog, wtslog, '.') % log-log plot
 % title('Marmoset, compare FLNe vs Vol-Source')
 %title('Marmoset, compare FLNe vs Vol-Target')
title('Marmoset, compare FLNe vs Vols-Source*Target')
xlabel('log(vol-T)   (mm^3)')  %xlabel('vol-Source (mm^3)')
ylabel('log(FLNe)')

%% 1.2.1 Plot LinkWts, LNe or FLNe vs Dist {for exp
 % calc DistCol at # 1.8.3, from LinkList [calc at #  1.8.1 (LNe) or 1.8.2 (FLNe)
 fprintf('\n Calc LinkList from rescaled Adj-new (LNe) \n')
figure; %figure('position',figposition, 'units','pixels'); hold on; % big
plot(DistCol, wtslog, '.', 'MarkerSize', 16); hold on
 title('Marmoset Brain log10(LNe) vs ListDist(mm)  ')
%title('Marmoset Brain, link weight vs dist ','Fontname','Times New Roman','FontSize',14,'FontWeight','Bold')
xlabel('dist(mm)','Fontname','Times','FontSize',14,'FontWeight','Bold') 
ylabel('log10(LNe wt) ','Fontname','Times New Roman','FontSize',14,'FontWeight','Bold');
 
% highlight weakest wt (< 9e-6)
for i=1:17  % fitst 17 pts in wt-ordered list
plot(DistCol(i), wtslog(i), '.', 'Color', 'red','MarkerSize', 16)
end
tmp=LinkListSortWt(1:17,:);
% save S-T
tmpS=cell(17,1);
for i=1:17
tmpS=Acrn{tmp(i,1)};
tmpT=Acrn{tmp(i,2)};
[tmpS tmpT]
end
clear tmp

% highlight next group of weaker, wt (< 5e-5)
for i=18:283  % fitst 17 pts in wt-ordered list
plot(DistCol(i), wtslog(i), '.', 'Color', [1, 0.65, 0],'MarkerSize', 16) % orange
end
% top 40%, by wt is at cutoff 10^0.0033 : 1.0076
%  ie. rows 2085:3474 : ie. 1390 links

% drop 10% - ie. top 90% is rows 348:3474 

%%       & Curve Fitting tool
% use Matlab Curve Fitting Toolbox
cftool % invokes GUI, import data
% ** refit - so edit below **
% linear fit, to logwt:  -2.019 - 0.0863*dist ; mean 0 about that line; std : xx
 plot([0 22], [-2.019 -3.9176], 'k-', 'LineWidth', 2, 'LineSmoothing', 'on') % along the mean of the data
% & 3*sigma CI: nb. log scale

plot([0 22], [-2.019+3*0.967 -3.9176+3*0.967], 'r--', 'LineWidth', 1.5, 'LineSmoothing', 'on') % 3*sigma
plot([0 22], [-2.019-3*0.967 -3.9176-3*0.967], 'r--', 'LineWidth', 1.5, 'LineSmoothing', 'on')
 
% get residuals - in log-space
[nll ~]=size(LinkListSortDist);
resLogWt=zeros(nll,1);
resLogWt(:)=wtslog(:) - (-2.019 - 0.0863*DistCol(:) );
 % [mean(resLogWt) std(resLogWt)] % -4.6944e-04, 0.9665
 % short < 10mm
 %  [mean(resLogWt(1:2116)) std(resLogWt(1:2116)) ]
 % short < 10mm
  % [mean(resLogWt(2117:end)) std(resLogWt(2117:end)) ]
% explore wt cut-off
% log-wt < 1e-5
  %mean(LinkListSortWt(1:17, 4)) % link dist: 13.5
% log-wt < 2e-5
  %mean(LinkListSortWt(1:65, 4)) % link dist: 11.23
% log-wt < 3e-5
  % mean(LinkListSortWt(1:129, 4)) % link dist: 10.7
% residuals - in original-space (FLNe)
  % reswt=10.^resLogWt;

  %% 1.2.1 Plot LinkWts, LNe or FLNe vs Dist {for exp fit
 % calc DistCol at # 1.8.3, from LinkList [calc at #  1.8.1 (LNe) or 1.8.2 (FLNe)
 fprintf('\n Plot Wt-Dist from rescaled Adj-new (LNe) \n')
figure; %figure('position',figposition, 'units','pixels'); hold on; % big
plot(DistCol, wtslog, '.', 'MarkerSize', 16); hold on
 title('Marmoset Brain log10(LNe) vs ListDist(mm)  ')
%title('Marmoset Brain, link weight vs dist ','Fontname','Times New Roman','FontSize',14,'FontWeight','Bold')
xlabel('dist(mm)','Fontname','Times','FontSize',14,'FontWeight','Bold') 
ylabel('log10(LNe wt) ','Fontname','Times New Roman','FontSize',14,'FontWeight','Bold');
 
%% highlight weakest wt (< 9e-6) : for FLNe
for i=1:17  % fitst 17 pts in wt-ordered list
plot(DistCol(i), wtslog(i), '.', 'Color', 'red','MarkerSize', 16)
end
tmp=LinkListSortWt(1:17,:);
% save S-T
tmpS=cell(17,1);
for i=1:17
tmpS=Acrn{tmp(i,1)};
tmpT=Acrn{tmp(i,2)};
[tmpS tmpT]
end
clear tmp

% highlight next group of weaker, wt (< 5e-5)
for i=18:283  % fitst 17 pts in wt-ordered list
plot(DistCol(i), wtslog(i), '.', 'Color', [1, 0.65, 0],'MarkerSize', 16) % orange
end
% top 40%, by wt is at cutoff 10^0.0033 : 1.0076
%  ie. rows 2085:3474 : ie. 1390 links

% drop 10% - ie. top 90% is rows 348:3474 

 %% 1.2.1a Log-Log Plot LinkWts, LNe or FLNe vs Dist {for power law fit
% nb. LinkListDist(:,3) is wt; LinkListDist(:,4) is dist
 % calc DistCol at # 1.8.2or3, from LinkList [calc at #  1.8.1 (LNe) or 1.8.2 (FLNe)
    % wtslog=log10(wts);  mean(wtslog) % now -2.80
 fprintf('\n Plot Wt-Dist from rescaled Adj-new (LNe) \n')
   % figure; plot(DistCol, wts) % debug
DistLn=log(DistCol); WtsLn =log(wts);
figure; %figure('position',figposition, 'units','pixels'); hold on; % big
  % plot(DistCol, wtslog, '.', 'MarkerSize', 16); hold on % was log 10
  %loglog(DistCol, wts); hold on
  plot(DistLn, WtsLn, '.', 'MarkerSize', 16); hold on % now ln = log_e
  %xlim([0.5 25])
 title('Marmoset Cortex log_e(LNe) vs log_e(ListDist (mm))  ')
%title('Marmoset Cortex, link weight vs dist ','Fontname','Times New Roman','FontSize',14,'FontWeight','Bold')
xlabel('log_e(dist)  (mm)','Fontname','Times','FontSize',14,'FontWeight','Bold') 
ylabel('log_e(LNe wt) ','Fontname','Times New Roman','FontSize',14,'FontWeight','Bold');
xlim([-0.2 3.25])
 % cftool: GUI for curve fitting
  line ([0 3], [6.9   (6.9 - 2.0*3) ], 'Color', 'k', ...
     'LineWidth', 1.5); % linear fit 6.9 -2.0*d(mm)
 
 %% 1.2.1b Plot LinkWts, LNe vs Dist {+ exp fit)
 % calc DistCol at # 1.8.3, from LinkList [calc at #  1.8.1 (LNe) or 1.8.2 (FLNe)
 fprintf('\n Plot Wt-Dist from rescaled Adj-new (LNe) \n')
figure; %figure('position',figposition, 'units','pixels'); hold on; % big
plot(DistCol, WtsLn, '.', 'MarkerSize', 16); hold on
 title('Marmoset Cortex: log_e(LNe) vs ListDist(mm)  ')
%title('Marmoset Cortex, link weight vs dist ','Fontname','Times New Roman','FontSize',14,'FontWeight','Bold')
xlabel('dist(mm)','Fontname','Times','FontSize',14,'FontWeight','Bold') 
ylabel('log_e(LNe, wt) ','Fontname','Times New Roman','FontSize',14,'FontWeight','Bold');
 line ([1 24], [4.77- 0.219  5.77 - 0.219*24], 'Color', 'k', ...
     'LineWidth', 1.5); % linear fit 4.77 - 0.219d
 
%% 1.3 node wt-Degree / Strength  :: nb. numbering out of order
DegIn=sum(Afull,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
DegOut=sum(Afull,2); % i->j; row sum

 %sum(Afull(:,1)) % check : these are "1" ??
figure; hist(DegIn,50)
title('Marmoset Brain, wt Degree-In, histogram','FontSize',14)

% get stats on Deg.
binrange=linspace(0,3.4E4,100); % get range from hist plot already
bincounts=histc(DegIn,binrange);
figure; bar(binrange, bincounts,'histc')
% appears to be exp fall off up to ~ 4xE7; only greater; & are 1 off 
figure; semilogy(binrange, bincounts, '.','MarkerSize',14)
hold on; title('Marmoset Brain, wt Degree-In,  log histogram','FontSize',14)
xlabel('Degree (bins)','FontSize',14) 
ylabel('log counts','FontSize',14); hold off
% clear binrange bincounts
%
[degsort, degindex] = sort(DegIn); % sorted list & get new index
hdin = degindex(end-9:end);      % list top 10, Hi Deg nodes (Deg = 113 .. 16  )
hdin=flipud(hdin);  % list are largest... smaller
 %mdi = degindex(n-43:n-20)';   % next 24, Medium deg   (Deg = 15 .. 10 )
% list top 10, 20
md=mean(DegIn)
sd=std(DegIn)
 for ii=1:length(hdin)
     i=hdin(ii);
     fprintf(' %4.0f %s %7.3f %5.2f \n', i, Acrn{i}, DegIn(i), (DegIn(i)-md)/sd )
 end
clear md sd 

%% 1.3a #Links: un-wt Deg -  binarise full 116x116 A:  UN-wt links
[sii sjj v]=find(Afull);  
Alogical=full(sparse(sii, sjj, 1));    % "mask" of selected entries, to re-form sq array of 1's.
 % nb. need to use sparse to form the array
A1 = Alogical.*1;  % pick out the selected entries; yields a "full" array
 %figure; spy(A1)
clear Alogical sii sjj v
length(find(A1(:))) % check 3474 links: 

% alt calc 
 %Afull= csvread('AdjFullMarmoset.csv');  % was calc from expt data
LinkList=adj2edgeL(Afull);  %LinkList=adj2edgeL(Afull);
  % use "filtered" LinklList (of 3462 links ~0)
  % LinkList=csvread('MarmosetLinkList.csv'); % i-j-wt
[nl ~]=size(LinkList)
wt1=ones(nl);
LinkList1=[LinkList(:,1:2) wt1]; % append wt=1
A1=edgeL2adj(LinkList1);
 % nb. need to use sparse to form the array
 %figure; spy(A1)
%
clear wt1 LinkList LinkList1
 % figure; spy(A1) % 1948 entries: banded & v sparse
 length(find(A1(:)) )
% count links (<~ 50-150)
DegIn1=sum(A1,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
 DegIn1=[DegIn1; 0; 0]; % append 2 missing 0's
DegOut1=sum(A1,2); % i->j; row sum 
figure; hist(DegIn1, 50)
title('Marmoset DegIn of binarised links')
axis([0 110 0 10]) % nb. big # at 0: not sampled
%
figure; hist(DegOut1, 50); hold on
title('Marmoset DegOut of binarised links')
 %axis([0 80 0 130]); 
figure; stem (DegIn1, '.'); hold on
stem (DegOut1, 'r.');
legend('DegIn1', 'DegOut1')
title('Marmoset Deg of binarised links')
% all
figure; hist(DegIn1, 50); xlim([2 Inf]); hold on; 
hist(DegOut1, 50, 'Color', 'r') % magenta
ylim([0 20]);

% all - joint p[lot
figure; hist([DegIn1, DegOut1], 50); hold on % Fails, diff # in & out 
title('Marmoset Distribution of in- and out-links'); legend('k-In', 'k-Out')
%axis([1 110 0  20]) % omit large peaks at 0 

%% plot & fit distrn
DegIn11=DegIn1(find(DegIn1 >0));   % need non-zero, eg for logNorm
DegOut11=DegOut1(find(DegOut1 >0)); 

figure; hist([DegIn11, DegOut11], 50); hold on % Fails, diff # in & out 
title('Marmoset Distribution of in- and out-links'); legend('k-In', 'k-Out')
axis([1 110 0  20]) % omit large peaks at 0 
 clear DegIn11 DegOut11
 
% csvwrite('Marmoset_Adj1unwt.csv', A1);  % save  binarises Adj (UN-wt Afull) 

%% 1.3a.1 probe missing In-link info
linkRatio=DegIn1./DegOut1;
 % figure; stem(linkRatio, '.'); xlabel('Node # (all nodes')
title('Marmoser (LNe) k-in/ k-out');
figure; plot(linkRatio, '.','MarkerSize',18) % lin fit dominated by 0's
xlabel('Node # (all nodes'); title('Marmoser (LNe) k-in/ k-out');
linkRatNon=linkRatio(find(linkRatio)); % get the 55 non-zero ratios
figure; plot(linkRatNon, '.','MarkerSize',18); xlabel('Node # (injection sites only'); hold on  
title('Marmoset: #linksIn/#Out (non-0)  ');
 text(11+1, 4.25, Acrn{11}); text(32+1, 4.68, Acrn{32}); text(33+1, 6.33, Acrn{33}); % label outliers
figure; hist(linkRatNon, 50);  % skewed normal

 clear linkRatio linkRatNon

%% 1.3.2  wt/link distrn
% wt/Link-In, per anatomical area (node)
WtLinkIn=DegIn./DegIn1;
 length(find(isnan(WtLinkIn)))  % check on 0 #Links:
 WtLinkIn(find(isnan(WtLinkIn)))=0; % omit the 0/0: NaN
figure; hist(WtLinkIn, 50); title('Marmoset, Distribution of Wt-Deg-In/ k-In ')
figure; %stem(WtLinkIn)
plot(WtLinkIn , '.','MarkerSize',12); title('Marmoset, Wt-Deg-In/ k-In vs  Node-ID '); hold on 
  % this plot shows av ~ xx wt/Link; 
[mean(WtLinkIn) mode(WtLinkIn) std(WtLinkIn)] % finds mode=low, since exp distrn
length(find(WtLinkIn==0)) %  61 of; only 55 targets measured
WtLinkIn1=WtLinkIn(find(WtLinkIn >0));   % need non-zero, eg for logNorm
 % typical range is xx wt/Link

 % wt/Link-Out, per area (node)
WtLinkOut=DegOut./DegOut1;
 length(find(isnan(WtLinkOut)))  % check on 0 #Links: 
 WtLinkIn(find(isnan(WtLinkOut)))=0; % omit the 0/0: NaN
figure; hist(WtLinkOut, 50); title('Marmoset, Distribution of Wt-Deg-Out/ k-Out ')
% use dfittool % title('Marmoset, Distribution of Weight-In/link-In ')

figure; %stem(WtLinkOut)
plot(WtLinkOut , '.','MarkerSize',12); title('Marmoset, Wt-Deg-Out/ k-Out vs  Node-ID '); hold on 
  % this plot shows av ~ xx wt/Link; 
[mean(WtLinkOut) mode(WtLinkOut) std(WtLinkOut)] % finds mode=xxx, since exp distrn
 % typical range is xx wt/Link


%% 1.4 form other LinkLists, with CutOffs, from Adj - for InfoMap calc
% cf noTarget & noTarget lists at A.2.1, line 785.
% for Ctx: 55 x 55 Source AND Target sites
 % cf. A.2.1, below, to find many missing Targets (116-52 of)
 %  ie. omit Cols 4-9, etc
  % nb. this discards a few Sources also: ie. Rows 4-9, ...
LinkList55=adj2edgeL(A);
length(find(A(:))) % 1854 links

dlmwrite( 'MarmosetPairList52.txt',LinkList55,'Delimiter',' '); % needs space deliminter 

% un-wt
tmp=ones(length(LinkList55),1);
LinkList55(:,3)=tmp;
clear tmp
dlmwrite( 'MarmosetPairList55u.txt',LinkList55,'Delimiter',' '); % needs space deliminter 

% check A-unwt
Au=adjLu2adj(LinkList52); % un-wt
figure; spy(Au)
title('Marmoset: Adj 55x55 un-wt')
length(find(Au(:)) ) % 1854, ok

% count links (<~ 50)
DegIn1=sum(Au,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
DegOut1=sum(Au,2); % i->j; row sum 
figure; hist(DegIn1)
title('Marmoset 55x55 DegIn of binarised links')
figure; hist(DegOut1); hold on
title('Marmoset 55x55 DegOut of binarised links')
 %axis([0 80 0 130]); 
figure; stem (DegIn1, '.'); hold on
stem (DegOut1, 'r.');
legend('DegIn1', 'DegOut1')
title('Marmoset 55x55 Deg of binarised links')


%% 2.0 LNe,1:  raw data, Test:  normalised to 1 mm^3 target vol
 % 55 injn sites (targets)
 % cf XL workings %& calc:  Marmoset_raw1_2_sort.xlx
LNe1=csvread('LNE1.csv');
wt=LNe1(:,3);
figure; hist(wt,50); title('Marmoset original LNe1 distribution  ')
wtlog=log10(wt);
figure; hist(wtlog,50)
 %axis([-5.3 0 0 150]) % omit large peak at "0")
title('Marmoset original log-10 (LNe1) distribution  ')

%% 2.1 calc LN from LNe link list & LNt,e,i data;  cf #4.1 below
fprintf('\n Calc #LN(i-j) from raw data \n')
 % data is [ Node Id, TotLN, LNe, LMi], with ~2-9 repeats
 % calc using the average LN & LNi for the 55 targets - in cols [5, 6]
  % as calc by hand in XL. see Marmoset_raw1_2_sort.xlsx : combines 2 data sets
rawLN1=csvread('rawLN1.csv');
LNtot1=sum(rawLN1(:,5)); % tot of all mean(LN-repeats); now 905k
FLNdenom = zeros(116,2);
FLNav = zeros(116,3); % to store av values [FLNtot, FLNe FLNi] for ea target
count=1; % nb. need to skip over repeats - just get the av.
for i =1:143
   if rawLN1(i,5) > 0 % then av for these Targets recorded
       ii=rawLN1(i,1); % index recorded on col 1 (w repeats)
       [ii rawLN1(i,6)] % debug
       FLNdenom(ii,1)=ii; % target index [only 55 of the 116]
       FLNdenom(ii,2)=(LNtot1-rawLN1(i,6)); % av LN, after repeats
       FLNav = [rawLN1(i,5), rawLN1(i,6), (rawLN1(i,5)-rawLN1(i,6))]; % calc FLNe also
       count=count+1;
   end
end
count % debug - check
clear count
sum(FLNdenom(:,2))/55 % to get mean
[min(FLNdenom(:,2)) sum(FLNdenom(:,2))/55 max(FLNdenom(:,2))]
% variation about mean:  (9.05-9.01)/9.01 : 0.0044;  ie. 0.4% : Trivial !
%
figure; hold on; plot(FLNdenom(:,2), '.', 'MarkerSize', 11);
 % axis([0 150 2.13e6 2.18e6])
title('Marmoset Brain:  Tot#LN - LNi (average at 55 targets)  ')
ylabel('All(LN) - LNi '); xlabel('Node ID # ');

%% 2.2 scan LinkList & undo normalisation by denom: get LN, Anew 
 % LinkList of [i j Wt] calc at #1.8.2 [line 322], above, from raw wts 
[nl ~] = size(LinkList)
extraCol=zeros(nl,1); % append to the 3-col LinkList
%LinkList= [LinkList, extraCol]; % 4th col wil be LNe(i-j); 5th LNe/Tot#LN
for i=1:nl
    LinkList(i,4)=LinkList(i,3)*FLNdenom(LinkList(i,2),2); % * denom(target)
    LinkList(i,5)=LinkList(i,3)*FLNdenom(LinkList(i,2),2)/LNtot1; % and as frn of tot#LN
end

clear extraCol 
figure; plot(LinkList(:,3), LinkList(:,5), '.', 'MarkerSize', 11); % good st line!!
% calc new Adj - need 3 cols in link list now:
LinkListNew=[LinkList(:,1:2) LinkList(:,4)];  % is [i j LN(i-j)/TotLN ]
Anew=zeros(n);
for i=1:nl
    Anew(LinkListNew(i,1), LinkListNew(i,2))=LinkListNew(i,3);
end

%% 2.3 new (LNe) wt Deg-In,Out (Strength)
DegIn=sum(Anew,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
DegOut=sum(Anew,2); % i->j; row sum 
figure; hist(DegIn); title('Marmoset rescaled wt-DegIn ')
figure; hist(DegOut); hold on; title('Marmoset rescaled wt-DegOut '); hold off
figure; plot(DegOut, DegIn, '.', 'MarkerSize', 9 ); title('marmoset: rescaled DegOut vs DegIn  ')

figure; stem (DegIn, '.'); hold on  %  ~ "wt-DegIn"
stem (DegOut, 'r.');  % wt-DegOut (LNe/mm^3 measure)
legend('wt-DegIn', 'wt-DegOut'); title('Marmoset rescaled (LNe) wt-Deg-In,Out '); hold off 

 

%% 2.4.1 check extra 61 rows (FLNe):
missingCols=[4:9, 12, 16:18, 20:24, 27:28, 33, 35:36, 49, 51:54, 58:60, 62:69 ...
   73:74, 83:88, 90:95, 98:102, 104:106,109, 115,116]; % check 61 cols : 116-55 ok
% for rows of non-targets
fprintf('\n check non Targets for links, using FLNe')
Anew61 = Afull(missingCols, :);
 nLinks61= length(find(Anew61)) % 1611 links 
 figure; spy(Anew61) % debug : 55 rows x 116 cols
 DegIn61=sum(Anew61,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
  nIn61 = length(find(DegIn61)) % is 55 : # of Injn sites, ok 
 DegOut61 = sum(Anew61,2); % i->j; row sum 
 figure; stem(DegOut61); title('Marmoset (FLNe), missing links, wt-Deg-Out')
   text(59+1, 1.33, Acrn{missingCols(59)}); text(33+1, 0.78, Acrn{missingCols(33)})
   text(14+1, 0.62, Acrn{missingCols(14)}); text(7+1, 0.66, Acrn{missingCols(7)}) % label outliers
 figure; stem(Anew61(59, :)); title('Marmoset (FLNe), missing links from V3, wt-Deg-Out')
 figure; plot(log10(Anew61(59, :)))
%%
figure; stem(DegIn61); title('Marmoset (FLNe), missing links, wt-Deg-In')
   text(61+1.5, 0.73, Acrn{61}); % is Injn
   text(55-1, 0.59+0.02, Acrn{55}) % is Injn
    text(11+1.5, 0.58, Acrn{11}); % is Injn
    text(103+1.5, 0.53, Acrn{103}) % in Injn 
%   
 figure; subplot(2,1,1); hist(DegIn61); legend('wtDegIn'); ylabel('Counts')
  subplot(2,1,2); hist(DegOut61); legend('wtDegOut'); xlabel('Link weight (FLNe)')
  title('Marmoset (FLNe),wt Deg-In, Out ~ 61 non Targets')
  %clear missingCols Anew61 nIn61 

%% 2.4.1a #Links: un-wt Deg -  for 61 non-Targets
fprintf('countf missing links ~ 61 non Targets')
[sii sjj v]=find(A61);  
Alogical=full(sparse(sii, sjj, 1));    % "mask" of selected entries, to re-form sq array of 1's.
 % nb. need to use sparse to form the array
A611 = Alogical.*1;  % pick out the selected entries; yields a "full" array
 %figure; spy(A1)
clear Alogical sii sjj v
length(find(A611(:))) % check 3474 links: 
 kIn61=sum(A611,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
 figure; stem(kIn61); %legend('k-In'); % 1611 : ok
 title('Marmoset (FLNe), missing links, k-In')
  nIn61 = sum(kIn61) % is 1611
  nSources61 = length(find(kIn61)) % is 55: = # Injn sites
 kOut61 = sum(A611,2); % i->j; row sum 
      % sum(A611(1, :)) % debug  % =  38, correct, row sum ok
 figure; stem(kOut61) ; %legend('k-Out');
  title('Marmoset (FLNe), missing links, k-Out')
 nOut61 = sum(kOut61) % 1611 : ok
 
 clear nIn* nOut* nSou*
 
%% 2.4.2 check extra 61 rows (LNe):
missingCols=[4:9, 12, 16:18, 20:24, 27:28, 33, 35:36, 49, 51:54, 58:60, 62:69 ...
   73:74, 83:88, 90:95, 98:102, 104:106,109, 115,116]; % check 61 cols : 116-55 ok
fprintf('\n check non Targets for links, using LNe')
% for rows of non-targets
Anew61 = Anew(missingCols, :);
 figure; spy(Anew61) % debug : 55 rows x 116 cols
 DegOut61 = sum(Anew61,2); % i->j; row sum 
 figure; stem(DegOut61)
 
  clear missingCols Anew61
  
%% 3.0 symmetrise the 116 x 55 T <- S Adj
% simple-minded symm. wt Adj
 %As= (tril(Afull) + triu(Afull') ); % symmetric, retain wt. / but looses too much data??
 %length(find(As)) % again 3964 links : 2* 1847 un-dir links [vs orig 3474 links]
 %figure; spy(As)  % looks wierd! data lost
% focus only on the missing columns of data
At=Afull';
missingCols=[4:9, 12, 16:18, 20:24, 27:28, 33, 35:36, 49, 51:54, 58:60, 62:69 ...
   73:74, 83:88, 90:95, 98:102, 104:106,109, 115,116]; % check 61 cols : 116-55 ok
As=Afull;
As(:, missingCols) =At(:, missingCols);
length(find(As)) % now 5065 links : nb. 1854 un-dir links in Act(55x55)  [vs orig 3474 links]
figure; spy(As)
title('Marmoset 116x116 Adj, symmetrised missing links')
clear At

% Link list: get i - j - wt list
LinkListS=adj2edgeL(As); % 

 % use this (A1-symm in InfoMap calc (C++) of modules (14/1/19)
 % dlmwrite( 'MarmosetPairListS.txt',LinkListS,'Delimiter',' '); % InfoMap code wants .txt , space delim   

%% 3.1  un-wt  binarised links 
[nl ~]=size(LinkListS)
wt1=ones(nl);
LinkList1=[LinkListS(:,1:2) wt1]; % append wt=1
A1=edgeL2adj(LinkList1);
 % nb. need to use sparse to form the array
figure; spy(A1)
%
clear wt1  LinkList1 %LinkListS
 % figure; spy(A1) % 1948 entries: banded & v sparse
 length(find(A1(:)) )
% count links (<~ 50-150)
DegIn1=sum(A1,1)'; % i<-j;  col sum (1st index): produces a row % show as Col
DegOut1=sum(A1,2); % i->j; row sum 
figure; hist(DegIn1, 50)
title('Marmoset DegIn of sym, binarised links')
axis([0 110 0 10]) % nb. big # at 0: not sampled

figure; hist(DegOut1, 50); hold on
title('Marmoset DegOut of sym, binarised links')
 %axis([0 80 0 130]); 
figure; stem (DegIn1, '.'); hold on
stem (DegOut1, 'r.');
legend('DegIn1', 'DegOut1')
title('Marmoset Deg of symmetrised, binarised links')


%% 4.0 examine spread of Injection Volumes:
% data in Marmoset_raw_1_2_sort.xlxs

InjVols=csvread('InjnVols.csv'); % 143 cases: many repeats at 55 target sites
figure;
plot(InjVols(:,1), InjVols(:,2), '.', 'MarkerSize', 11);
title('Marmoset Brain:  Injections volumes (repeats) at 55 areas  ')
ylabel('Injn Vol (mm3) '); xlabel('Node ID # ');

%% 4.1 re-examine all LN data (3 exp measurements x 143 repeats)
 % has 143 expts, on 55 target sites - with ~ 2-9 repeats
 % this csv has 2 lines of 0's at end? use onl 143 entries
rawLN=csvread('rawLN.csv'); % data is [ Node Id, TotLN, LNe, LMi], with ~2-9 repeats
 %rawLN=dlmread('rawLN.csv'); % rawLN=dlmread('rawLN.csv', ','); % trials
% calc denom for FLNe calc
FLNdenom = zeros(143,1);
LNtot=sum(rawLN(:,2)); % tot of all repeats: 2.16m
for i =1:143
   FLNdenom(i)= LNtot-rawLN(i,4); % sum all except this target LNi
end

(max(FLNdenom)-min(FLNdenom))/mean(FLNdenom) % range ~ 2% relative to mean
[min(FLNdenom) mean(FLNdenom) max(FLNdenom)]
figure; hold on; plot(FLNdenom, '.', 'MarkerSize', 11);
 % axis([0 150 2.13e6 2.18e6])
title('Marmoset Brain:  #LN - LNi in 143 expts (repeats at 55 targets)  ')
ylabel('All(LN) - LNi '); xlabel('Node ID # ');


%% test calc using the average LN & LNi for the 55 targets - in cols [5, 6]
 % calc by hand in XL.
rawLN1=csvread('rawLN1.csv');
LNtot1=sum(rawLN1(:,5)); % tot of all mean(LN-repeats); now 905k
figure; hold on
title('Marmoset Brain:  LN (repeats); Av(#LN) & Av(LNi) at 55 targets  ')
ylabel('Denom: 1/[Tot(LN) - LNi] '); xlabel('Node ID # ');
% axis([0 55 0 1.0e6])
[min(rawLN1(:,5)) max(rawLN1(:,5)) ]
for i =1:143
    plot(i, rawLN1(i,2), '.', 'MarkerSize', 11); % plot all LN repeats
   if rawLN1(i,5) > 0
       plot(i, rawLN1(i,5), 'ks', 'MarkerSize', 6); % av LN, after repeats
       plot(i, rawLN1(i,6), 'go', 'MarkerSize', 6); % av LN,i after repeats
   end
end

LNtot1=sum(rawLN1(:,5)); % tot LN of 55 target averages: 905k
%
figure; hold on
title('Marmoset Brain, Denom: TotLN-Av(#LN,) at 55 targets  ')
axis([0 55 1.0e-6 1.2e-6])
ylabel('Denom: 1/[Tot(LN) - LNi] '); xlabel('Node ID # ');
count=1; % Target node counter
for i =1:143
   if rawLN1(i,5) > 0 % then its a Target
       plot(count, (1/(LNtot1-rawLN1(i,6))), 'ks', 'MarkerSize', 6); % av LN, after repeats
       count=count+1;
   end
end

[1/(min(LNtot1-rawLN1(:,6))) 1/(mean(LNtot1-rawLN1(:,6))) 1/(max(LNtot1-rawLN1(:,6))) ]  % range
% variation about mean:  (9.05-9.035)/9.035 : 0.0017;  ie. 0.2% Trivial !


%  APPENDIX         APPENDIX         APPENDIX


%% Appx.0 get RGB colours for areas: from Atlas files
 %NodeID=csvread('MarmosetAcrn139to116.csv'); % for the 116 Nodes (cf Links) vs 136 Labels
 % translates the full 139 List (Atlas) to the 116 nodes (Source Injn)
tmp=csvread('Marmoset135LabelsRGB.csv'); % get listed ID-RGB for the 135/136 Labels
% scan all labels:
NodeColors=zeros(n,3); % array for the 116 R-G-B colours
for i=1:length(NodeID)   
    if (NodeID(i) ~= 0 )% skip over these
        %fprintf('> gets here \n')
        %[i NodeID(i)] % debug
       thisArea=NodeID(i);
        % Labels{thisArea}; % debug
       ThisColor=tmp(i,2:4);
       NodeColors(thisArea, :)= tmp(i,2:4); % 1st col is ID
       fprintf('%u %u %s %u %u %u \n', i, thisArea, Acrn{thisArea}, ThisColor )   
    end % if      
end % loop
clear tmp thisArea ThisColor 
NodeColors=NodeColors./255; % need 0-1 range
csvwrite('Marmoset116NodesRGB.csv',NodeColors); % save file

NodeColors=csvread('Marmoset116NodesRGB.csv'); % read file, prev calc


%  APPENDIX         APPENDIX         APPENDIX
%% Appx.1 get Atlas files
% need nifti Tools; code from MouseBrainV1-3.m
% Atlas in sub-dir /marmoset_brain_template
% get brain volume; nb. file needs to be in this dir
% Lavels for 139 sites: 116 nodes + general laebels in atlas_labels.txt & MarmosetLabels139.csv 

fprintf('\n >> Marmoset Brain: Atlas \n')
bv = load_nii('atlas.nii'); % this is 86MB file : loads into a struct, multiple details:
bv.hdr.hk % info
bv.hdr.dime 
   % 21 fields: bv.hdr.dime.dim [5, 825, 63, 550, 1, 3, 1, 1]
   %   & bv.hdr.dime.vox_offset [ 352 ] 
bv.hdr.hist
 % 18 sub fields:  bv.hdr.hist.qoffset_x [16.56], _y [19.50], _z [19.72] % "mm"
 % bv.hdr.hist.srow_x [-0.04, 0, 0, 16.56], _y , _z 
 % bv.hdr.hist.originator [411, 24, 57, 0, 0]
% view_nii(bv); % shows 3 x views of slices of 3D image [value: 0-255]
bvol=double(bv.img); % nb floating pt, for interpolation
  % size(bvol) % get rid of dim "1" from the 5D file
bvol=squeeze(bvol); size(bvol) %  825x63x550x3 : RGB 3D image, a 4-D array

 % bvmask=load_nii('atlas_mask.nii'); % this is 29MB file  
 % view_nii(bvmask); % shows "0" or "1" oursite, inside brain tissue

%% Appx.1.1 3D vol of areas ID# - cf list in atlas_labels.txt: load into
 % nb. 135 useful labels: find the 116 nodes amongst these
% find CM of each area
fprintf('\n >> Marmoset Brain: Atlas segmentation file \n')
bvsegm=load_nii('atlas_segmentation.nii'); % this is 29MB file also  
% view_nii(bvsegm); % bvsegm.img is 825 x 63 x 550 uint8 voxels % opens
% interactive figure
                  % origin at (411,24, 57) 
                  % bvsegm.hdr.dime.pixdim(2:4) are pixel dimensions in mm
origin = [411,24, 57]
pixsize= bvsegm.hdr.dime.pixdim(2:4) % ???

% scan all:  or explore single area, below
  % count voxels
NodeCoord=zeros(n,3); % array for x-y-z coords
NodeVol=zeros(n,1); voxVol=zeros(n,1);
for i=1:length(NodeID)   
    if (NodeID(i) ~= 0 )% skip over these
        %fprintf('> gets here \n')
        %[i NodeID(i)] % debug
       thisArea=NodeID(i);
       Labels{thisArea}; 
       ii= find(bvsegm.img(:, :, :) == i); % need to use 'Label' ID# here;  finds 340,519 indices
       % length(ii)
       [ix jy kz]= ind2sub(size(bvsegm.img), ii);
       thisCM = [mean(ix) mean(jy) mean(kz)]; % get CM in vox coords
       CMmm= (thisCM-origin).*pixsize; % CMM, centered, in mm-coords
       NodeCoord( NodeID(i), :)= CMmm;  % nb. only 116 of these
       voxVol( NodeID(i))= length(ix);
       vol = length(ix)*pixsize(1)*pixsize(2)*pixsize(3);
       NodeVol( NodeID(i))= length(ix)*pixsize(1)*pixsize(2)*pixsize(3);
       fprintf('%u %u %s %f %f %f %f \n', i, NodeID(i), Acrn{NodeID(i)}, CMmm, vol )   
    end % if
      
end % loop
clear thisArea ii ix jy kz thisCM CMmm vol
% save coords to file % nb. 5 missing in Vol & no Injn sites : yields NaN
%csvwrite('CoordsMarmoset.csv', NodeCoord);
%csvwrite('VolsMarmoset.csv', NodeVol); % (mm^3) 

% now calc node surface areas
 %NodeVol=csvread('VolsMarmoset.csv'); % read saved file
tmp=(3*NodeVol/(4*pi) ).^0.333;  % radius
NodeArea=4*pi*(tmp).^2;

%% Appx.1.1  Tests: find coords of a single area: eg V1, V2: in pixel space
 % check against online atlas browsers
 % Label ID's in MarmosetLinksNodes.xls [nb. skip of non-targets]
 
 fprintf('\n >> Marmoset Brain Atlas:  explore single area \n')
  % find CM of single area
  % bvsegm=load_nii('atlas_segmentation.nii'); % this is 29MB file also
  % load NodeID at #1.1a, above
  
 %thisArea=126 % V1 % check coords by hand against on-line InjnSite coords
 %thisArea=2  % A-10
 %thisArea=25  % A-32
 %thisArea=132  % V-5
% i = 127 % OPt(65)
fprintf('\n >> Marmoset Brain Atlas:  explore single area \n')
 %i = 124  % check problem with OPt, PIR 
 i = 31  % check problem with A45, A47L  (D & V?)
 i = 127  %  : #108, V2
thisArea=NodeID(i) % index in file "atlas_labels.txt"
%thisArea=85 % force bt hand, if problem?
if thisArea ~= 0 % check valid index?
 %Labels{i} 
        % ii= find(bvsegm.img(:, :, :) == thisArea); %finds 340,519 indices
       ii= find(bvsegm.img(:, :, :) == i); % need to use the 'Label' ID#
       [ix jy kz]= ind2sub(size(bvsegm.img), ii);
       thisCM = [mean(ix) mean(jy) mean(kz)] % get CM in vox coords
       CMmm= (thisCM-origin).*pixsize % CMM, centered, in mm-coords
       %NodeCoord( NodeID(i), :)= CMmm;  % nb. only 116 of these
       fprintf('%u %u %s %f %f %f  \n', i, NodeID(i), Acrn{NodeID(i)}, CMmm) % nb. need to lookup ID#
 % examine slices
 % figure; hist(ix, 50)
 % figure; hist(jy, 50)
 % figure; hist(kz, 50)

% plot point cloud in voxel coords
figure; hold on  %figure('position',figposition, 'units','pixels'); hold on; 
% view(-85, 40)
hold on
 %plot3(ix, jy, kz, 'b.', 'MarkerSize', 4); % dense colours
plot3(ix, jy, kz, '.', 'Color', [0 0 0.9], 'MarkerSize', 4);
plot3(thisCM(1), thisCM(2), thisCM(3), 'd', 'Color', [0 0 0.5], 'MarkerFaceColor', [0 0 0.5], 'MarkerSize', 8) % mark CM
text((thisCM(1)-150), (thisCM(2) - 5), (thisCM(3)+75), Acrn(NodeID(i)), 'FontSize', 12) % label the area (offset)
axis([0 825 0 63 0 550]) % "cube" of voxels
grid on
xlabel('x (vox)  Lat - Med ','FontSize',12)
ylabel('y (vox)  Antr - Post ','FontSize',12); 

title('Marmoset Brain:  areas  A45 & A47L, in voxel space  ');
clear ii ix jy kz
else % check was not ok:
       fprintf('\n >> problem with chosen area \n') % mistake?
end % while "ok"

%% A.1.1a  now in mm space: find coords of a single area: eg V1, V2:
 % check against online atlas browsers
 % Label ID's in MarmosetLinksNodes.xls [nb. skip of non-targets]
 
 fprintf('\n >> Marmoset Brain Atlas:  explore single area \n')
  % find CM of single area
  % bvsegm=load_nii('atlas_segmentation.nii'); % this is 29MB file also
  % load NodeID at #1.1a, above
  
 %thisArea=126 % V1 % check coords by hand against on-line InjnSite coords
 %thisArea=2  % A-10
 %thisArea=25  % A-32
 %thisArea=132  % V-5
  % i = 127 % OPt(65)
  %i = 124  % check problem with OPt, PIR 
  %i = 31  % check problem with A45, A47L  (D & V?)
 i = 127  %  : #108, V2
thisArea=NodeID(i) % index in file "atlas_labels.txt"
%thisArea=85 % force bt hand, if problem?
if thisArea ~= 0 % check valid index?
 %Labels{i} 
        % ii= find(bvsegm.img(:, :, :) == thisArea); %finds 340,519 indices
       ii= find(bvsegm.img(:, :, :) == i); % need to use the 'Label' ID#
       [ix jy kz]= ind2sub(size(bvsegm.img), ii);
       thisCM = [mean(ix) mean(jy) mean(kz)] % get CM in vox coords
       CMmm= (thisCM-origin).*pixsize % CMM, centered, in mm-coords
       %NodeCoord( NodeID(i), :)= CMmm;  % nb. only 116 of these
       fprintf('%u %u %s %f %f %f  \n', i, NodeID(i), Acrn{NodeID(i)}, CMmm) % nb. need to lookup ID#
 % examine slices
 % figure; hist(ix, 50)
 % figure; hist(jy, 50)
 % figure; hist(kz, 50)

% plot point cloud in voxel coords
figure; hold on  %figure('position',figposition, 'units','pixels'); hold on; 
% view(-85, 40)
hold on
 %plot3(ix, jy, kz, 'b.', 'MarkerSize', 4); % dense colours
plot3(ix, jy, kz, '.', 'Color', [0 0 0.9], 'MarkerSize', 4);
plot3(thisCM(1), thisCM(2), thisCM(3), 'd', 'Color', [0 0 0.5], 'MarkerFaceColor', [0 0 0.5], 'MarkerSize', 8) % mark CM
text((thisCM(1)-150), (thisCM(2) - 5), (thisCM(3)+75), Acrn(NodeID(i)), 'FontSize', 12) % label the area (offset)
axis([0 825 0 63 0 550]) % "cube" of voxels
grid on
xlabel('x (vox)  Lat - Med ','FontSize',12)
ylabel('y (vox)  Antr - Post ','FontSize',12); 

title('Marmoset Brain:  areas  A45 & A47L, in voxel space  ');
clear ii ix jy kz
else % check was not ok:
       fprintf('\n >> problem with chosen area \n') % mistake?
end % while "ok"


  %   App & Test codes
  
  
%% A.1.2 3D volume RGB plots (pixel space)
% load bv at A.1 line 255 above, needs the Atlas volume (.nii) : Marmoset brain:  [825, 63, 550]
  %bv = load_nii('atlas.nii'); % at A.1, above; this is 86MB  % a NIFTI tool
  %bvol=double(bv.img); % nb floating pt, for interpolation
bvol=bv.img; % keep uint8 form
size(bvol)
% get rid of dim "1" from the 5D file
bvol=squeeze(bvol); size(bvol) %  825x63x550x3 : RGB 3D image
% Tests:
cvol=bvol(:,:,:); % 3D array only
cvol(find(cvol(:) ==255))=0; % mask out the "white" surround
figure; hist(cvol(:),50); axis([5 255 0 7e5])
[x, y, z] =meshgrid(1:1:825, 1:1:63, 1:1:550);
  %figure('position',figposition, 'units','pixels'); hold on;
  %slice(x, y, z, bvol, 410 ,40, 200) % slice in mid-planes
  %  figure; slice(x, y, z, cvol, 50 ,Inf, Inf) % slice in x-y mid-plane only
  % v slow
% Tests: get the rgb channel
 % rgbbvol=permute(bvol, [4,1,2,3]); % get last page
 % rgbbvol=rgbbvol(1,:,:,:);

% get mid slice at z(pixels)=nnn
 % bvol=squeeze(bvol); size(bvol)  %  825x63x550x3 : RGB 3D image
 midslice=bvol(:,:,350, :); size(midslice) % keep RGB channels
 midslice=squeeze(midslice); size(midslice)  % again! need to eleim "1" dim 
  % TESTS: figure; h =imshow(midslice); % fails
  % midslice=squeeze(midslice);
   % midslice(find(midslice(:) ==255))=0; % mask out the "white" surround
  %figure; hist(midslice(:),50); title('Marmoset Mid Slice')
  % axis([5 255 0 500])
  % figure; image(midslice);  title('Marmoset Mid Slice')

%% A.1.2a  AGAIN:  mid slice from image data
 zpix=350
 midslice=bvol(:,:,zpix, :); size(midslice) % horoz slice at height z = nnn pix
 % midslice=bvol(:,:,350); size(midslice) % this one is 2D ! ok
 midslice=squeeze(midslice); size(midslice)  % 3D= 2DX[RGB], again need to elim "1" dim 
 midslice=uint8(midslice); % Need integer
  %figure; imshow(midslice) % ok
%
h=imshow(midslice);
[imdata flag] =getimage(h);  % gets R-G-B at each pixel in 2D array
  % imdata=uint8(imdata); % Need integer
  % figure; hist(imdata(:),50)
gray_image=rgb2gray(imdata); % now 2D array 
gray_imageflip=gray_image'; % need to flip x-y axes
  %figure; hist(gray_image(:),50)% now in [0 1] scale 
 %figure; hg=imshow(gray_image); title('Marmoset Mid Slice z=350 pix ') 
 % now shows gray-scale 2D 825 x 63 pic x-y swapped
% plot mid-plane transverse (x,y) slice, in pixel space  / a la MouseBr
% desired z position of the image plane.
imgzposition = zpix; % half height,  pixel space
%  corners of the x-y domain for data occurs.
min_x = 1; max_x = 825; min_y = 1; max_y = 63; 
 %min_x = 1; max_x = 63; min_y = 1; max_y = 825;  % x-y swapped?
% save image as file
csvwrite('MarmosetGrayImageSlice.csv', gray_imageflip);

% plot the image plane using surf. (pixel space)
figure; hold on; grid off %  still draws a "box" ?? / nb. x-y swapped already
surf([min_x max_x],[min_y max_y],repmat(imgzposition, [2 2]),...
    gray_imageflip,'facecolor','texture', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.0) % nb. no Edge (box)
 %colormap(gray) % default
title('Marmoset Brain:  Mid Slice z=350 pix ')
ylabel('y (pixels)'); xlabel('x (pixels)'); % nb. swap
colormap(gray(256)); grid on;  axis([0 825 0 63  0 550]) % nb confusion w x-y sawp??

% plot3(0, 0, 0, '+k','MarkerSize', 45); % Mark origin (mm)
plot3(411, 24, (550-57), '+k','MarkerSize', 30); % Mark origin (pix) / z down
 %plot3(24, 411, 57, '+k','MarkerSize', 45); % Mark origin (pix) / flip
% add RHS - symmetric
gray_imageR=flipud(gray_image);
gray_imageflipR=gray_imageR'; 
surf([min_x max_x],[min_y max_y],repmat(imgzposition, [2 2]),...
    gray_imageflipR,'facecolor','texture', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.0) 
hold off

% now convert to phys space % apply vox -> mm tranform 
% nb x-y  appear to be swapped?
min_x= -bv.hdr.hist.srow_x(4) % is (825-411)*0.04 ok
min_y= -bv.hdr.hist.srow_y(4) % nb. these are in mm!!
min_z= -bv.hdr.hist.srow_z(4) % offset / check 550-57)*0.04: 19.7200 
max_x= min_x + bv.hdr.hist.srow_x(1) % convert "width" of x-voxels to x-mm
max_y= min_y + bv.hdr.hist.srow_y(2) % & for y
max_z= min_z+bv.hdr.hist.srow_z(4)

%  pix-mm conversion, by hand??
origin = [411,24, 57] % from bv.hdr.hist.originator 
 % pixsize= bvsegm.hdr.dime.pixdim(2:4) 
pixsize=[0.0400 0.5000 0.0400] % set manually, avoide reloading large image
min_xmm=(min_x-origin(1))*pixsize(1); max_xmm=(max_x-origin(1))*pixsize(1); 
min_ymm=(min_y-origin(2))*pixsize(2); max_ymm=(max_y-origin(2))*pixsize(2);
min_zmm=(min_z-origin(3))*pixsize(3); max_zmm=(max_z-origin(3))*pixsize(3);
imgzpositionmm=(imgzposition-origin(3))*pixsize(3);  %  centered, in mm-coords
imgzpositionmm=min_z+ bv.hdr.hist.srow_z(3)*(-imgzposition)

% plot the image plane using surf.
figure; hold on; grid on
 %figure( 'units','centimeters'); hold on; grid on; %grid off %  still draws a "box" ??
%surf([min_xmm max_xmm],[min_ymm max_ymm],repmat(imgzpositionmm, [2 2]),...
%    gray_imageflip,'facecolor','texture', 'FaceAlpha', 0.2)  % has box
% eliminate box
surf([min_xmm max_xmm],[min_ymm max_ymm],repmat(imgzpositionmm, [2 2]),...
    gray_imageflip,'facecolor','texture', 'FaceAlpha', 0.2,'edgecolor', 'none')  % 
colormap(gray)
plot3(0, 0, 0, '+k','MarkerSize', 30);
title('Marmoset Brain:  Mid Slice z= -5.72 mm ')
ylabel('y (mm)'); xlabel('x (mm)'); % nb. swap
axis([-10 10 -10 20 -19 0])
% add RHS - symmetric
gray_imageR=flipud(gray_image);
gray_imageflipR=gray_imageR'; 
surf([min_xmm max_xmm],[min_ymm max_ymm],repmat(imgzpositionmm, [2 2]),...
    gray_imageflipR,'facecolor','texture', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.0) 
hold off
 %clear bvol

%  pix-mm conversion, by hand??
origin = [411,24, 57] % from bv.hdr.hist.originator 
 %origin = [411,0, 57] % reset to mid-line 
 % pixsize= bvsegm.hdr.dime.pixdim(2:4) 
pixsize=[0.0400 0.5000 0.0400] % set manually, avoid reloading large image
min_x = 1; max_x = 825; min_y = 1; max_y = 63; 
min_z = 1; max_z = 550; imgzposition = zpix; % "half" height, pixel space

min_xmm=(min_x-origin(1))*pixsize(1); max_xmm=(max_x-origin(1))*pixsize(1); 
min_ymm=(min_y-origin(2))*pixsize(2); max_ymm=(max_y-origin(2))*pixsize(2);
min_zmm=(min_z-origin(3))*pixsize(3); max_zmm=(max_z-origin(3))*pixsize(3);
imgzpositionmm=(imgzposition-origin(3))*pixsize(3);  %  centered, in mm-coords

%% Other tests - saggital (y,z) slice
sagslice=bvol(200,:,:,:); % size(sagslice) % sagital slice
sagslice=squeeze(sagslice); size(sagslice) % 63 x 550 x 3[RGB]
hs =imshow(sagslice);
[imdatas flag] =getimage(hs);  % gets R-G-B at each pixel in 2D array
  % figure; hist(imdata(:),50)
gray_images=rgb2gray(imdatas); % now 2D array 
figure; hg=imshow(gray_images); title('Marmoset Mid Slice x=200 pix ') 

 % plot mid-line vertical / saggital (y,z) slice  / a la MouseBr
% desired z position of the image plane.
figure;
imgzpositionx = 200; % half height,  pixel space
%  corners of the x-y domain for data occurs.
min_x = 1; max_x = 825; min_y = 1; max_y = 63; min_z=1; max_z=550;
 %min_x = 1; max_x = 63; min_y = 1; max_y = 825;  % x-y swapped?
% plot the image plane using surf.
figure; hold on; grid off %  still draws a "box" ??
% slice at x=imgzpositionx , in y-z plane
%surf([min_y max_y],[min_z max_z],repmat(imgzpositionx, [2 2]),...
%    gray_image,'facecolor','texture', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.0) % ng. no Edge (box)
surf(repmat(imgzpositionx, [2 2]),[min_y max_y],[min_z max_z],...
    gray_images,'facecolor','texture', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.0) % ng. no Edge (box)

 %colormap(gray) % default
title('Marmoset Brain:  Sag Slice, x=200 pix ') 
colormap(gray(256)); grid on;  axis([0 825 0 63  0 550]) % nb confusion w x-y sawp??

% plot3(0, 0, 0, '+k','MarkerSize', 45); % Mark origin (mm)
plot3(411, 24, (550-57), '+k','MarkerSize', 45); % Mark origin (pix) / z down
 %plot3(24, 411, 57, '+k','MarkerSize', 45); % Mark origin (pix) / flip
 
 
%% tests  %midslice=bvol(:,:,348);  % is 825 x 63 pixels; & 228 rows x 264 cols  (~ 0.4MB)
  %midslice(find(midslice(:) ==0))=256; % mask out the "black" surround : to white
% figure; image(midslice)  % check layout:  nb: shows x-y swap 
% colormap(gray) % default [64?]
midslice(find(midslice(:) ==255))=0; % mask out the "white" surround
figure; hist(midslice(:),50); title('Marmoset Mid Slice')
figure; image(midslice); title('Marmoset Mid Slice')

 % ??????
corslice=tmp(:,30,:); % coronal
 %squeeze(corslice); % get rid of dim "1" % fails!
corslice=reshape(corslice, 825, 1650);
corslice(find(corslice(:) ==255))=0; % mask out the "white" surround
figure; hist(corslice(:),50); title('Marmoset Coronal Slice')
figure; image(corslice); title('Marmoset Coronal Slice')
colormap(gray)  % shows 3 x 550 wide slices together!
%
sagslice=tmp(200,:,:); % sagital
sagslice=reshape(sagslice, 63, 1650); % get rid of dim "1"
sagslice(find(sagslice(:) ==255))=0; % mask out the "white" surround
figure; hist(sagslice(:),50); title('Marmoset Sagl Slice')
figure; image(sagslice); title('Marmoset Sagl Slice')   % again shows 3 x 450 wide slices together!
% colormap(gray) 


 
 %% plot mid-transverse(x,y)-plane slice / from M..B..V2.m [228 x 264 image]
% desired z position of the image plane.

imgzposition = 348; % @ ~half height,  pixel space
%  corners of the x-y domain for data occurs & convert to phys space (mm)
min_x=bv.hdr.hist.srow_x(4); min_y=bv.hdr.hist.srow_y(4); min_z=bv.hdr.hist.srow_z(4); % offset
% apply vox -> mm tranform (*2 for 1/2 grid & *4 for 1/4 grid sub-sample)
max_x= min_x + bv.hdr.hist.srow_x(1)*264; % convert "width" of x-voxels to x-mm
% y & z appear to be swapped?
max_y= min_y + bv.hdr.hist.srow_y(2)*264; % & y
imgzpositionmm=min_z+ bv.hdr.hist.srow_z(3)*imgzposition; % & z

% convert to phys space (mm - cf. 2.3) % apply vox -> mm tranform 
% nb x-y  appear to be swapped?
min_x=bv.hdr.hist.srow_x(4); min_y=bv.hdr.hist.srow_y(4); min_z=bv.hdr.hist.srow_z(4); % offset
max_x= min_x + bv.hdr.hist.srow_x(1)*228; % convert "width" of x-voxels to x-mm
max_y= min_y + bv.hdr.hist.srow_y(2)*264; % & for y
imgzpositionmm=min_z+ bv.hdr.hist.srow_z(3)*imgzposition;

% plot the image plane using surf.
figure; hold on; grid off %  still draws a "box" ??
surf([min_x max_x],[min_y max_y],repmat(imgzpositionmm, [2 2]),...
    midslice','facecolor','texture', 'FaceAlpha', 0.05)  % check plot 
colormap(gray)
clear bvol


%% A.2.0 test code to READ info from online data
   % LinkList=csvread('LinkList.csv'); % fails??
 
% **  first time: need to use GUI:  uiimport & nb. specify text
      % appears to be in [Target Source wt] order
 % assemble as [Source Target wt]  LinkList 
figure; hist(A(:), 50) % exp fall off; v long tail 
wtlog=log10(wt);
figure; hist(wtlog,50)
axis([-5.3 0 0 150]) % omit large peak at "0")
title('Marmoset log10(FLNe) distribution  ')
 
 Acrn=unique(Source); % 117 of [incl junk row]
 % delete empty row? / get 116 areas 
 tmp=cell(116,1);
 for i=1:116
     tmp{i}=Acrn{i+1};
 end
 % check alt. list
 
 AcrnT=unique(Target); % 55 of [+ junk row]
 % save this list
 dlmwrite('MarmosetAcrn.csv', Acrn);
 % find numerical index to Arcn used
 LinkList=zeros(6326,3);
% scan list & identify ACrn: need numerical index 
 for il=1:length(LinkList)
   for ia=1:116
    if strncmp(char(Source{il}) , char(Acrn{ia}), 7) == 1 % nb. up to 7 char long
       LinkList(il,1)=ia;
    end
   end
   for ia=1:116
    if strncmp(char(Target{il}) , char(Acrn{ia}), 7) == 1
       LinkList(il,2)=ia;
    end
   end
   LinkList(il,3)=wt(il); % wt: FLNe
 end
 % check last line: had NaN
 % save the list [S, T, wt]
 csvwrite('MarmosetLinkList.csv', LinkList);
 % form Adj
 Afull=edgeL2adj(LinkList);
 figure; spyc(Afull) % 3474 links : not square!
 csvwrite('AdjFullMarmoset.csv', Afull);
 
% add node Acrn as well
 LinkListAcrn=cell(6325,5);
 for i =1:6325
  LinkListAcrn{i,1}=num2cell(LinkList(i,1)); % Source
  LinkListAcrn{i,3}=num2cell(LinkList(i,2)); % Target
  LinkListAcrn{i,5}=num2cell(LinkList(i,3)); % wt
  LinkListAcrn{i,2}={Acrn{LinkList(i,1)}};
  LinkListAcrn{i,4}={Acrn{LinkList(i,2)}};
 end
% csvwrite('LinkListACrnMarmoset.csv', LinkListAcrn); % fails

%% A.2.1 nb many missing Targets - so need to get sub-set
  %  to form square Adj matrix
  % scan cols, find all "0"
noTarget=[]; yesTarget=[]; 
none=1; % counter
for i=1:n
    if Afull(:,i) ==0
        noTarget(none)=i;
    else
        yesTarget(none)=i;
        none=none+1;
    end    
end
noTarget=noTarget'; % get Col vec
yesTarget=yesTarget';
%
clear none
% save these Acrn: need to do individually?
  Acrn55=cell(length(yesTarget),1);
  out=1; % counter
  for i=1:length(yesTarget)
  Acrn55{i}=Acrn{yesTarget(out)};
  out=out+1;
 end
 clear out

 %% form square Adj of S-T list
 %Afull= csvread('AdjFullMarmoset.csv'); % calc prev. 
A = Afull(yesTarget, yesTarget); % pick the subset, from Target list: calc just below:
[n55 ~]=size(A) % finds 61 ??
csvwrite('AdjMarmoset.csv', A); % save this calc. 

%
% noTarget=find( (Afull(1,:) ==0) & (Afull(40,:) ==0) & (Afull(78,:) ==0)  & (Afull(89,:) ==0))'; % by inspection
   % appears incomplete? (11/10/19)
% hand checked: 64 missing: so 52 x 52 Adj
 AcrnNoTarget=cell(length(noTarget),1);
 in=1; % counter
 for i=1:length(noTarget)
  AcrnNoTarget{i}=Acrn{noTarget(in)};
  in=in+1;
 end
  clear in
 %  get Target list:
 ListT=zeros(116,2);
 for i=1:116
  ListT(i,1)=i;
  ListT(i,2)=1; % assume target
 end
  ListT(noTarget,2)=0; % idenify & set  No Targets
  yesTarget=find(ListT(:,2)==1); % get list of Target sites (ie. the others)

%% this Adj is actually 116 rows (S) x 55 cols (T)
% check symmetry
Aup=triu(Afull);
Alow=tril(Afull);
Aup=triu(Afull); Alow=tril(Afull);
Aout=Alow(:); Ain=Aup(:);
Aup=Aup'; % match structure
Ain=Aup(:);
% plot Aij vs Aji : to check reciprocal links?
figure; plot(Ain, Aout, '.', 'MarkerSize', 12) % very scattered!!

%
 AcrnTarget=cell(length(yesTarget),1);
 in=1; % counter
 for i=1:length(yesTarget)
  AcrnTarget{in}=Acrn{yesTarget(i)};
  in=in+1;
 end
  clear in
% 
% form square Adj of S-T list
Afull= csvread('AdjFullMarmoset.csv'); % calc prev. 
A = Afull(yesTarget, yesTarget); % pick the subset, from Target list
[n ~]=size(Afull)

figure; spyc(A) % 1758 links; does not appear symm.
csvwrite('AdjMarmoset.csv', A);

figure('position',figposition, 'units','pixels'); hold on; 
plotArray(A) % coarser grid % 1239 links: Looses many??

% check stats of Adj
max(A(:)) % FLNe 0.6055
rank(A)  % 52; so well condt.  ?? calls svd 
det(A)  % "0" ~ 1e-43
trace(A) % 0: no diag.

%% clean code to read saved image volume data, midslice: use this in Vis codes
gray_imageflip=csvread('MarmosetGrayImageSlice.csv'); % prev calc RGB; x-y flipped
gray_imageflip=uint8(gray_imageflip); % Need integer; nb image was x-y flipped
gray_imageflipR=fliplr(gray_imageflip); % reflect LHS to get RHS image
% set x-y box & desired z position of the image plane.
zpix=350;  pixsize=[0.0400 0.5000 0.0400]; % set manually, avoid reloading large image
min_x = 1; max_x = 825; min_y = 1; max_y = 63; 
min_z = 1; max_z = 550; imgzposition = zpix; % "half" height, pixel space
%  pix-mm conversion, % apply vox -> mm tranform 
origin = [411,24, 57] % from bv.hdr.hist.originator 
pixsize=[0.0400 0.5000 0.0400]; % set manually, avoid reloading large image
min_xmm=(min_x-origin(1))*pixsize(1); max_xmm=(max_x-origin(1))*pixsize(1); 
min_ymm=(min_y-origin(2))*pixsize(2); max_ymm=(max_y-origin(2))*pixsize(2);
min_zmm=(min_z-origin(3))*pixsize(3); max_zmm=(max_z-origin(3))*pixsize(3);
imgzpositionmm=(imgzposition-origin(3))*pixsize(3);  %  centered, in mm-coords
clear origin pixsize min_x max_x min_y max_y min_z max_z imgzposition zpix
% test plot in mm space the image plane using surf.
figure; hold on; grid on  % nb. need to eliminate "box"
surf([min_xmm max_xmm],[min_ymm max_ymm],repmat(imgzpositionmm, [2 2]),...
    gray_imageflip,'facecolor','texture', 'FaceAlpha', 0.2,'edgecolor', 'none')  % 
colormap(gray)
plot3(0, 0, 20, '+k','MarkerSize', 30);
title('Marmoset Brain:  Mid Slice z= -5.72 mm ')
ylabel('y (mm)'); xlabel('x (mm)'); % nb. swap
axis equal; axis([-10 10 -10 20 0 20])
% add RHS - symmetric about mid-line
surf([min_xmm max_xmm],[min_ymm max_ymm],repmat(imgzpositionmm, [2 2]),...
    gray_imageflipR,'facecolor','texture', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.0) 
hold off

%% examine spread of Injection Volumes:
% data in Marmoset_raw_1_2_sort.xlxs

InjVols=csvread('InjnVols.csv'); % 143 cases: many repeats at 55 target sites
figure;
plot(InjVols(:,1), InjVols(:,2), '.', 'MarkerSize', 11);
title('Marmoset Brain:  Injections volumes (repeats) at 55 areas  ')
ylabel('Injn Vol (mm3) '); xlabel('Node ID # ');




