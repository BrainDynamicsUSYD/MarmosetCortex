% B Pailthorpe, U Sydney:  (2019-23)
 % Apply Temporal Network Theory to Marmoset Brain connectivity (21/10/19)
   % to trace eiolving links in 3D
    % based on Mouse:  MBtemporalNetV2.m (2018)
    % V2  (23/10/19): multi-sense paths;  paths around Hubs, Connectors ('23)
% Dependencies: use MarmosetReadDataV1.m to input data files:
  % NodeCoord, Acrn, Lobes, idx (module IDs), A (Adj, link wts; also Anew)
  
% assumes v ~ 1 m/s = 1000 mm/s = 1 mu-m/ms; + tauS, synaptic delay ~ 2 ms.
   % Ipsi Lateral version ((1/2/18) draw on RHS (injection side)
     %  build i- j- k Out links search; order by dist or 'time'; draw-as-you-go
     %  accumulate array of i-j-k links; sort by t(ij), then
%                   t(ijk), to draw "shells" of neighbors
%  derived from Graph Vis of Mouse brain connectivity VisMouseBrainV7c,etc (Mar-Aug'17)
   
   % data files in  .../MatlabFiles/MarmosetBrain
   % codes in          .../MarmosetBrain/CodesOther/: std network codes; simpler path on different computers
% uses Cell mode ( Sections):
%   Sections (cells): work down these as needed
%      #1 - read data; setup; basics:  Now done in separate code: MBreadData.m
%       1.6 - set colours, by Module s ID: Use Paxinos Atlas colours 
% Here-in:
%   11. time-slice "layers" of Adj(i-j-k) Out links
%   11b.1 & Draw ...
%   12. time-slice "layers" of Adj(i-j-k) In links
%   12.1 & Draw ..

% 1. Read data files & setup (from Oh, Nature 2014): use code MBreadData.m 
  % may use matlab toolboxes for NIFTI files etc - in sub-paths

% Other Dependencies:
% ** needs array FlowWt: calc in MB_flows.m : use NodeFlow from Rosval's InfoMap C++ calc
FlowWt=zeros(n); % temporary fix {used in print outs 

% get coords, if needed (mm space -from Allen on-line viewer etc): Update missing coords
 
% > >   >   >   >   >   >   >   >
% > >   >   >   >   >   >   >   >  OUT      OUT     OUT  time ordered i-j links  +k    +k 2nd NN
%% 11b. first build  i-j & i-j-k OUT lists;  draw later
  %  use code: 10a PLOT All nodes:  explore SINGLE Region /  Module; OUT Links  : base plot 
  % use dr as proxy for dt (assume, say, v-signal = 1 m/sec = 1 mm/mSec)
  % nb. array FlowWt calc in code:  MB_flows.m (wt av per link)
  % Pick source node(ii); Cut-Off weight for links; nb. find/if >=  or <=, for Top or Botm links
tau=0.5;  % time slice size
contactTimes = [0.5:tau:25.0]; % ascending vector of all possible contact times
 % nb Marmoset has links up to 24.5 mm long
tauS=2.0;  % est of Synaptic delay (mS, or mm for v=1m/s)
distijk= zeros(1,5); % first row: set up array [i, j, k, tij, tijk]
CutWt=5.0e-5; %CutWt=78.0  % CutWt=10  % CutWt=2.5085   % cut off wt xx for link wts: top 40% / or, and wt<=xxx for botm 10% 
fprintf(' >>>>  Marmoset, DlpFC, build list of OUT Links, wt>= %6.2f [at #11b] \n', CutWt)


% picksources = these Nodes!  %  group of sensory outputs
%ThisList=[107, 52];  % Vis-1, Aud-1 
%ThisList=[29, 1, 30];  % SS / SSp
%ThisList=[13, 14];  % A2a,b (Cing)
ThisList=[2, 32, 48, 33, 45, 46, 47, 31];  % DlpFC group, in Ant-Post order

totLinks=0; irow=1;  %  accumulate count
%OUT
jj=[]; alljj=[]; dxs=[]; dys=[]; dzs=[]; drs=[]; % arrays to accumulate results
%
for i=1:length(ThisList)
    ii=ThisList(i);
fprintf(' %6.0f %s %s ; Module: %4.0f \n', ii, Acrn{ii}, Lobes{ii}, idx(ii) )
Xc=NodeCoord(ii,1); Yc=NodeCoord(ii,2);
Zc=NodeCoord(ii,3);   % coords of sphr centre (mm);  L side only so far
%OUT
   jj = find(A(ii,:)); % all links
      %jj = find(A(ii,:)>=CutWt); %  find [strong] Out links (1st NN From i- across i-th row): checked against DegOut
      %jj = find(A(ii,:)<=CutWt & A(ii,:)>1); %  find only weak Out links 
   totLinks=totLinks+1; % just count this i-j link here
    %wts = round(A(ii,jj));  % check the weights
   dxs=abs(NodeCoord(jj,1)-Xc); % L-M dist % nb. need to use actual x-coord
   dys=abs(NodeCoord(jj,2)-Yc); % A-P dist;  col vec: i-jList
   dzs=abs(NodeCoord(jj,3)-Zc); % D-V
    %drs=sqrt(dxs.^2+dys.^2+dzs.^2);
   thedr= sqrt(dxs.^2+dys.^2+dzs.^2); % ColVec
   thisdr=[ii*ones(length(jj), 1), jj', thedr]; % assemble mx3 ColVec of [i j dr]
   drs= vertcat(drs, thisdr);  % nb vert cat. of mx3 arrays
   alljj = [alljj jj]; % add new 1-st NN
end 
    jj=alljj;
   [dsort, dindex] = sort(drs(:,3)); % sorted (by dr, 3rd col) list & get new index
   dijrsort=drs(dindex,:); % reassemble the i-j-dr list
   clear dxs dys dzs thedr alljj 
    %figure; hist(drs,50) % debug
    %dindex(1:10)'
    %dsort(1:10)'
%   
prevT=0;  % for ( , ] timeSlice
 %fprintf('\n  Node#  Module  wt(i-j) dr(mm) Acrn Region  \n') % header
 %fprintf('\n  Node#  Module  wt(i-j) Flow(i-j)  dr(mm) Acrn Lobe  \n') % another header
fprintf('\n  Node-i  Node-j Module wt(i-j) dr(mm) Acrn Lobe  \n') % another header

for l=1:length(contactTimes)  % step through time slices; 
    %fprintf(' %10.0f %8.2f \n', l, contactTimes(l)) % debug
    for j=1:length(jj)
        ii=dijrsort(j,1);  %ii=dijrsort(dindex(j),1);   % need to find the i also
        jjj =dijrsort(j,2); 
        %[ii jjj] % debug
          %jjj= jj(dindex(j));  % start with closest NN;  nb dindex maps sorted drs to orig find list
       if ( dsort(j) <= contactTimes(l) & dsort(j) > prevT ) % limit to this ( , ] timeSlice only
       fprintf(' %6.0f %6.0f %7.0f %8.5f %5.1f   %s %s \n', ii, jjj, idx(jjj), A(ii,jjj), dsort(j), Acrn{jjj}, Lobes{jjj} )                
       xj=NodeCoord(jjj,1); yj=NodeCoord(jjj,2); zj=NodeCoord(jjj,3);
         %[xj, yj, zj]  % debug
         
      % % % %  % % % %
      % follow 2nd NN - OUT links
        kk = find(A(jjj,:)); % all links
         %kk = find(A(jjj,:)>=CutWt); %  find only strong Out links from the j-node:
        %kk = find(A(jjj,:)<=CutWt & A(ii,:)>1); %  find only weak Out links 
          %out=length(kk) %debug
        totLinks=totLinks+length(kk);
          %wts = round(A(jjj,kk));  % check the weightsfor j=1:2  %length(jj) % 1st few only, for now
        dxjk=abs(NodeCoord(kk,1)-NodeCoord(jjj,1)); %col vec % nb. need to use actual x-coord, not "idx"
        dyjk=abs(NodeCoord(kk,2)-NodeCoord(jjj,2)); % 
        dzjk=abs(NodeCoord(kk,3)-NodeCoord(jjj,3)); drjk=sqrt(dxjk.^2+dyjk.^2+dzjk.^2); %clear dxs dys dzs 
          %drs'
        [drjksort, drjkindex] = sort(drjk); % sorted list; closest/smalest first & get new index
        kList=kk(drjkindex);
        for k=1:length(kList) % closest first
          kkk= kList(k);   %kkk= kk(kList(k)) 
          tik=dsort(j)+tauS+drjksort(k); % tij + t-synapse + tjk : 2 steps
          distijk(irow,:)= [ii jjj kkk dsort(j) tik]; % accumulate info
          irow=irow+1;  % need a row for every k
        end  % of k-loop
      
      end  % if - time slice
      
    end % of j-loop
 
 prevT=contactTimes(l); % update startT for time slice

end % of timeSlice-loop
tikCol=distijk(:,5);
[tsort, tindex] = sort(tikCol); % sorted list & get new index
distijkSort=distijk(tindex,:); % now sort the whole array, by incrs tijk  
 clear tsort tindex tikCol
 clear i j jjj jlist k kkk  Xc Yc Zc  % prevT  % ii jj kk  
totLinks

%% 11.b.1 Draw:  Now, loop over i-j OUT links; then scan j-k list in progressive t-slices
%
fprintf(' >>>>  Marmoset,  trace OUT Links, time ordered; in MarmTemporalNetV1.m (at #11b.)  \n')
figure('position',figposition, 'units','pixels'); hold on; 
set(gcf,'Renderer','OpenGL');
set(gca, 'Color', [0.7 0.7 0.7]) % background gray
 %axis vis3d  % freeze aspect ratio (during rot)
axis equal  % get aspect ratio correct
 %set(gca, 'Projection', 'perspective');  % orthographic vs perspective view 
view(-64, 32);
  %  plot saved surface: approx to plial: at Isosuf=rface of brain-vol: 50, in mm-space
  % patch(FVmm, 'FaceColor', [0.7 0.7 0.7],'FaceAlpha', 0.03,'EdgeColor','none'); % for 1/2 res
% set mid-transverse(x,y)-plane slice : as reference in 3D plots; min_x etc set above, #1.4.1  
surf([min_xmm max_xmm],[min_ymm max_ymm],repmat(imgzpositionmm, [2 2]),...
    gray_imageflip,'facecolor','texture', 'FaceAlpha', 0.2,'edgecolor', 'none')  % mm coords
% add RHS - symmetric
surf([min_xmm max_xmm],[min_ymm max_ymm],repmat(imgzpositionmm, [2 2]),...
    gray_imageflipR,'facecolor','texture', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.0)
colormap(gray) % default [64?]
floorgridOnly % draw 5 mm grid & marker on floor

% Labels:
  %title('Marmoset Brain (mm): All areas; Out Links  ','FontSize',16,'FontWeight','bold')
title('Marmoset Brain: DlpFC: Out Links  ','FontSize',16,'FontWeight','bold')
text( 6.0, -1, 15.0, 'DlpFC: Out Links  ','FontSize',16,'FontWeight','bold'); % for zoom-in label
text( 6.0, -1, 13.0, 'module colours  ','FontSize',12,'FontWeight','bold'); %
  % text( 5.0, -1, 11.0, 'close links (d<7mm)
    % ','FontSize',12,'FontWeight','bold'); % optional
axis([-11 11 -10 20 0 20])
%plot3(0,0, 20, '+k','MarkerSize', 35); % Mark origin
line([0 0], [-1.3 1.2], [20 20], 'Color', 'k') % finish the 3D cross
line([0 0], [0 0], [18.7 21.3], 'Color', 'k')
line([-1.3; 1.3],[0; 0], [20; 20], 'Color', 'k','LineWidth',1)
text( 0, -1.6, 20.1, 'P','FontSize',16); text( 0, 1.9, 20, 'A','FontSize',16); % Label directions
text( 0.2, -0.3, 21.2, 'D','FontSize',16); text( 0, -0.3, 18.8, 'V','FontSize',16);
xlabel('x (mm)  Lat - Med ','FontSize',14)
ylabel('y (mm)  Antr - Post  ','FontSize',14)  % nb. x-y swapped above - in cols of FV.vertices
zlabel('z (mm)  Ventr - Dors ','FontSize',14)

  % axis off % remove white background - helps see colours

% Plot all nodes: coords of sphr centre (mm);  L side only so far 
for i=1:n
     plot3(NodeCoord(i,1), NodeCoord(i,2), NodeCoord(i,3), '.', 'Color', nodecolors(i,:),'MarkerSize', 20) % orig coord
        %text( (0-0.1), (NodeCoord(ii,2)-0.1), (NodeCoord(ii,3)+0.1), num2str(ii),'FontSize',8); % number the pt.
  %text( (NodeCoord(ii,1)-0.25), (NodeCoord(ii,2)-0.25), (NodeCoord(ii,3)+0.1), Acrn{ii},'FontSize',10); % L side % label node  
     %pause % debug
end
%
% + + + + + +
% overlay 2nd source (i) from here
 %  jj = find(A(ii,:)); % all links % found above
% do i-j outer Loop; draw i-j Link once only; then scan j-k List
% find min-T for first 2ndNN j-k Link:
minTk=distijkSort(1,5);  % this list ordered by t-ijk: this is the first t-ijk to look for
tau=0.5;  % time slice size % reset, to be sure
% when t-ij >= minTk, then search for j-k links:

% Plot selected Source Node[s] & their OUT Links (wt > 0, etc)
% plot a  sphere at ea. node/area  % colour by module ID
transp=0.6;  % for Alpha (line drawing): hi wt links heavier
%
for i=1:length(ThisList)
    ii=ThisList(i);
    fprintf(' %4.0f %s %s links to: \n', ii, Acrn{ii}, Lobes{ii} )
   xi=NodeCoord(ii,1); yi=NodeCoord(ii,2); zi=NodeCoord(ii,3); % Mark Sources;  ATlas coords
   plot3(xi, yi, zi, '.', 'Color', nodecolors(ii,:),'MarkerSize', 30) % schematic coord
  text( (xi-0.25), (yi-0.25), (zi+0.1), Acrn{ii},'FontSize',10);  %,'Fontname','Times New Roman','FontWeight','Bold','FontSize',10); % R side % label node   
end

% OUT :: scan the list found above: ie  kept jj-list
%jj(dindex) %debug % sorted list of j-1stNN
%tmpj=jj(dindex'); %debug

[nLinks ~] = size(distijk);
prevT=0;  % for ( , ] timeSlice
kListRow=1; % for scanning j-k (t-ijk) links
fprintf(' %4.0f %s %s \n', ii, Acrn{ii}, Lobes{ii} ) % header

[nLinks ~] = size(distijk);
prevT=0;  % for ( , ] timeSlice
kListRow=1; % for scanning j-k (t-ijk) links
fprintf('\n  Node#-i  Node#-j  Module  wt(i-j) dr(mm) Acrn Lobe  \n') % header
  % fprintf('\n        step#   t-slice  \n') % header
%prevT=contactTimes(11-1)  % Debug
for l=1:length(contactTimes)  % step through successive time slices; 
    %fprintf('\n  step, TimeWindow: %6.0f %8.2f %8.2f ', l, prevT, contactTimes(l)) % debug
    %fprintf('\n   Node#  Module  wt(i-j) dr(mm) Acrn Region  \n') % header
    for j=1:length(jj)
        ii =drs(dindex(j),1);   % need to find the i also, to draw i-j link
        jjj=drs(dindex(j),2);  % start with closest NN;  nb dindex maps sorted drs to orig find list
         %jjj= jj(dindex(j));  % start with closest NN;  nb dindex maps sorted drs to orig find list
        %fprintf(' step %2.0f work on j: %3.0f d-ij %7.1f \n', j, jjj, dsort(j) ); % debug
      if ( dsort(j) <= contactTimes(l) & dsort(j) > prevT ) % limit to this (t-l-1, t-l] timeSlice only
      fprintf(' %6.0f %6.0f %7.0f %8.5f %7.2f   %s %s \n', ii, jjj, idx(jjj), A(ii,jjj), dsort(j), Acrn{jjj}, Lobes{jjj} )             
      xi=NodeCoord(ii,1); yi=NodeCoord(ii,2); zi=NodeCoord(ii,3); % need Source coord
      xj=NodeCoord(jjj,1); yj=NodeCoord(jjj,2); zj=NodeCoord(jjj,3); % & Target
          %[xj, yj, zj] % debug 
      plot3(xj, yj, zj, '.','Color', nodecolors(jjj,:), 'MarkerSize', 25);  % plot & label this j-Node point
      text( (xj-0.25), (yj-0.25), (zj+0.1), Acrn{jjj},'FontSize',10); % Label the pt.
      text( (xj-0.25), (yj-0.25), (zj-0.15), Lobes{jjj},'FontSize',9, 'Color', [0.3 0.3 0.3]); % Label the target Region (faint)
      % colour Out i ->j link by the Target(jjj) node
      patch([xi; xj; xi], [yi; yj; yi], [zi; zj; zi], 'k', 'EdgeColor', ...
               nodecolors(jjj,:), 'LineSmoothing','on', 'EdgeAlpha', transp, 'FaceColor', 'none'); % stronger colour / alpha 0.05-0.2
      arrowhead([xi, yi, zi],[xj, yj,zj], nodecolors(jjj,:), transp);  % add an OUT (i -> j) arrowhead 85% along line
      % Strong links:
        if A(ii,jjj)>=0.1; %  highlight Stronger Out links; only label these points
          line([xi; xj], [yi; yj], [zi; zj],'Color', nodecolors(jjj,:),'LineWidth',2)
          %text( (xj-0.1), (yj-0.1), (zj), Acrn{jjj},'FontSize',10); % Label the pt.
        elseif A(ii,jjj)>=0.4; %  highlight Strongest Out links 
            %w5k=jjj
          line([xi; xj], [yi; yj], [zi; zj],'Color', nodecolors(jjj,:),'LineWidth',3)
          %text( (xj-0.1), (yj-0.1), (zj), Acrn{jjj},'FontSize',10); % Label the pt.
        end
        arrowhead([xi, yi, zi],[xj, yj,zj], nodecolors(jjj,:), transp+0.2);  % fix the arrowhead
      % % % %  % % % %
      % follow 2nd NN - OUT links: found above (so dont search again)
        % ?? if contactTimes(l) >= minTk  % any k-links to show yet?
      if dsort(j) >= minTk  % any k-links at/near this i-j time ?
           %fprintf(' +k links exist at this timeSlice \n') % debug
          for idummy=1:50 % "dummy" loop to scan j-k list
           %fprintf('      > link between:  %6.0f %5.0f t-ijk %4.2f \n', distijkSort(kListRow, 2), distijkSort(kListRow, 3), distijkSort(kListRow, 5) ); % debug
          thisJ=distijkSort(kListRow, 2); % check we have correct row?
          kkk=distijkSort(kListRow, 3); % select this k
              % ?? if (distijkSort(kListRow, 5) > dsort(j) & dsort(j) <= dsort(j)+ tau)  % check t-ijk in this slice
           if distijkSort(kListRow, 5) <= dsort(j)     % debug 
            fprintf('        & draw link:  %4.0f %4.0f %s %s \n', thisJ, kkk, Acrn{kkk}, Lobes{kkk} ); % debug
            xk=NodeCoord(kkk,1); yk=NodeCoord(kkk,2); zk=NodeCoord(kkk,3); % & draw this link:
            xj=NodeCoord(thisJ,1); yj=NodeCoord(thisJ,2); zj=NodeCoord(thisJ,3); % nb. not same "j" as above
            %[xk, yk, zk]  % debug / label the k-node
            text( (NodeCoord(kkk,1)-0.25), (NodeCoord(kkk,2)-0.25), (NodeCoord(kkk,3)+0.1), Acrn{kkk},'FontSize',10); % L side % label node  %OUT :: scan the list found above
              %text( (xk+0.1), (yk-0.15), (zk+0.1), Region{kkk}(1:4),'FontSize',10, 'Color', [0.5 0.5 0.5]); 
               %text( (CoordSize(kkk,1)-0.1), (CoordSize(kkk,2)-0.15), (CoordSize(kkk,3)-0.15), num2str(kkk),'FontSize',9); % shift slightly
            % draw j->k link;  fainter
            patch([xj; xk; xj], [yj; yk; yj], [zj; zk; zj], 'k', 'EdgeColor', ...
               nodecolors(kkk,:), 'LineSmoothing','on', 'EdgeAlpha', transp-0.2, 'FaceColor', 'none'); 
            arrowhead([xj, yj,zj],[xk, yk, zk], nodecolors(kkk,:), transp-0.2);  % add an OUT (j -> k) arrowhead
            kListRow=kListRow+1;  % increment row: go to next row in j-k list
            %pause
           else
                 % ?? kListRow=kListRow+1;  % also need to increment row here
               break % skip out & go on to next j, in outer loop
           end % if
          end % for loop over list 
      end % if - search for k at this time?
      
     end  % if - in this time slice
    
    end % j-loop
     %missing % to look up Atlas coords
     %wts = round(A(ii,missing))  % check the weights of missing Out links
     %[ii length(missing)]
     
 prevT=contactTimes(l); % update startT for time slice
 pause
end % l-timeSlice-loop

%clear i j jjj jlist k kk kkk wts  Xc Yc Zc ThisList prevT SecondNN % ii jj kk 

hold off  
  
%           Direct          Direct          Direct Links
%% 11.b.2 Draw:  Now, loop over direct i-j OUT links (only)
%
fprintf(' >>>>  Marmoset,  trace OUT Links, time ordered; in MarmTemporalNetV1.m (at #11b.)  \n')
figure('position',figposition, 'units','pixels'); hold on; 
set(gcf,'Renderer','OpenGL');
set(gca, 'Color', [0.7 0.7 0.7]) % background gray
 %axis vis3d  % freeze aspect ratio (during rot)
axis equal  % get aspect ratio correct
 %set(gca, 'Projection', 'perspective');  % orthographic vs perspective view 
view(-64, 32);
  %  plot saved surface: approx to plial: at Isosuf=rface of brain-vol: 50, in mm-space
  % patch(FVmm, 'FaceColor', [0.7 0.7 0.7],'FaceAlpha', 0.03,'EdgeColor','none'); % for 1/2 res
% set mid-transverse(x,y)-plane slice : as reference in 3D plots; min_x etc set above, #1.4.1  
surf([min_xmm max_xmm],[min_ymm max_ymm],repmat(imgzpositionmm, [2 2]),...
    gray_imageflip,'facecolor','texture', 'FaceAlpha', 0.2,'edgecolor', 'none')  % mm coords
% add RHS - symmetric
surf([min_xmm max_xmm],[min_ymm max_ymm],repmat(imgzpositionmm, [2 2]),...
    gray_imageflipR,'facecolor','texture', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.0)
colormap(gray) % default [64?]
floorgridOnly % draw 5 mm grid & marker on floor

% Labels:
  %title('Marmoset Brain (mm): All areas; Out Links  ','FontSize',16,'FontWeight','bold')
title('Marmoset Brain: DlpFC: Direct Out Links  ','FontSize',16,'FontWeight','bold')
text( 6.0, -1, 15.0, 'DlpFC: Direct Out Links  ','FontSize',16,'FontWeight','bold'); % for zoom-in label
text( 6.0, -1, 13.0, 'module colours  ','FontSize',12,'FontWeight','bold'); %
  % text( 5.0, -1, 11.0, 'close links (d<7mm)
    % ','FontSize',12,'FontWeight','bold'); % optional
axis([-11 11 -10 20 0 20])
%plot3(0,0, 20, '+k','MarkerSize', 35); % Mark origin
line([0 0], [-1.3 1.2], [20 20], 'Color', 'k') % finish the 3D cross
line([0 0], [0 0], [18.7 21.3], 'Color', 'k')
line([-1.3; 1.3],[0; 0], [20; 20], 'Color', 'k','LineWidth',1)
text( 0, -1.6, 20.1, 'P','FontSize',16); text( 0, 1.9, 20, 'A','FontSize',16); % Label directions
text( 0.2, -0.3, 21.2, 'D','FontSize',16); text( 0, -0.3, 18.8, 'V','FontSize',16);
xlabel('x (mm)  Lat - Med ','FontSize',14)
ylabel('y (mm)  Antr - Post  ','FontSize',14)  % nb. x-y swapped above - in cols of FV.vertices
zlabel('z (mm)  Ventr - Dors ','FontSize',14)

  % axis off % remove white background - helps see colours

% Plot all nodes: coords of sphr centre (mm);  L side only so far 
for i=1:n
     plot3(NodeCoord(i,1), NodeCoord(i,2), NodeCoord(i,3), '.', 'Color', nodecolors(i,:),'MarkerSize', 20) % orig coord
        %text( (0-0.1), (NodeCoord(ii,2)-0.1), (NodeCoord(ii,3)+0.1), num2str(ii),'FontSize',8); % number the pt.
  %text( (NodeCoord(ii,1)-0.25), (NodeCoord(ii,2)-0.25), (NodeCoord(ii,3)+0.1), Acrn{ii},'FontSize',10); % L side % label node  
     %pause % debug
end
%
% + + + + + +
% overlay 2nd source (i) from here
 %  jj = find(A(ii,:)); % all links % found above
% do i-j outer Loop; draw i-j Link once only; then scan j-k List
% find min-T for first 2ndNN j-k Link:
minTk=distijkSort(1,5);  % this list ordered by t-ijk: this is the first t-ijk to look for
tau=0.5;  % time slice size % reset, to be sure
% when t-ij >= minTk, then search for j-k links:

% Plot selected Source Node[s] & their OUT Links (wt > 0, etc)
% plot a  sphere at ea. node/area  % colour by module ID
transp=0.4;  % for Alpha (line drawing): hi wt links heavier
%
for i=1:length(ThisList)
    ii=ThisList(i);
    fprintf(' %4.0f %s %s links to: \n', ii, Acrn{ii}, Lobes{ii} )
   xi=NodeCoord(ii,1); yi=NodeCoord(ii,2); zi=NodeCoord(ii,3); % Mark Sources;  ATlas coords
   plot3(xi, yi, zi, '.', 'Color', nodecolors(ii,:),'MarkerSize', 30) % schematic coord
  text( (xi-0.25), (yi-0.25), (zi+0.1), Acrn{ii},'FontSize',10);  %,'Fontname','Times New Roman','FontWeight','Bold','FontSize',10); % R side % label node   
end

% OUT :: scan the list found above: ie  kept jj-list
%jj(dindex) %debug % sorted list of j-1stNN
%tmpj=jj(dindex'); %debug

[nLinks ~] = size(distijk);
prevT=0;  % for ( , ] timeSlice
kListRow=1; % for scanning j-k (t-ijk) links
fprintf(' %4.0f %s %s \n', ii, Acrn{ii}, Lobes{ii} ) % header

[nLinks ~] = size(distijk);
prevT=0;  % for ( , ] timeSlice
kListRow=1; % for scanning j-k (t-ijk) links
fprintf('\n  Node#-i  Node#-j Module wt(i-j) Flow(i-j) dr(mm) Acrn Lobe  \n') % header
  % fprintf('\n        step#   t-slice  \n') % header
%prevT=contactTimes(11-1)  % Debug
for l=1:length(contactTimes)  % step through successive time slices; 
    %fprintf('\n  step, TimeWindow: %6.0f %8.2f %8.2f ', l, prevT, contactTimes(l)) % debug
    %fprintf('\n   Node#  Module  wt(i-j) dr(mm) Acrn Region  \n') % header
    for j=1:length(jj)
        ii =drs(dindex(j),1);   % need to find the i also, to draw i-j link
        jjj=drs(dindex(j),2);  % start with closest NN;  nb dindex maps sorted drs to orig find list
         %jjj= jj(dindex(j));  % start with closest NN;  nb dindex maps sorted drs to orig find list
        %fprintf(' step %2.0f work on j: %3.0f d-ij %7.1f \n', j, jjj, dsort(j) ); % debug
      if ( dsort(j) <= contactTimes(l) & dsort(j) > prevT ) % limit to this (t-l-1, t-l] timeSlice only
    fprintf(' %6.0f %6.0f %7.0f %8.5f %8.4e %7.2f   %s %s \n', ii, jjj, idx(jjj), A(ii,jjj), FlowWt(ii,jjj), dsort(j), Acrn{jjj}, Lobes{jjj} )             
      xi=NodeCoord(ii,1); yi=NodeCoord(ii,2); zi=NodeCoord(ii,3); % need Source coord
      xj=NodeCoord(jjj,1); yj=NodeCoord(jjj,2); zj=NodeCoord(jjj,3); % & Target
          %[xj, yj, zj] % debug 
      plot3(xj, yj, zj, '.','Color', nodecolors(jjj,:), 'MarkerSize', 25);  % plot & label this j-Node point
      text( (xj-0.25), (yj-0.25), (zj+0.1), Acrn{jjj},'FontSize',10); % Label the pt.
      text( (xj-0.25), (yj-0.25), (zj-0.15), Lobes{jjj},'FontSize',9, 'Color', [0.3 0.3 0.3]); % Label the target Region (faint)
      % colour Out i ->j link by the Target(jjj) node
      patch([xi; xj; xi], [yi; yj; yi], [zi; zj; zi], 'k', 'EdgeColor', ...
               nodecolors(jjj,:), 'LineSmoothing','on', 'EdgeAlpha', transp, 'FaceColor', 'none'); % stronger colour / alpha 0.05-0.2
      arrowhead([xi, yi, zi],[xj, yj,zj], nodecolors(jjj,:), transp);  % add an OUT (i -> j) arrowhead 85% along line
      % Strong links:
        if A(ii,jjj)>=0.1; %  highlight Stronger Out links; only label these points
          line([xi; xj], [yi; yj], [zi; zj],'Color', nodecolors(jjj,:),'LineWidth',2)
          %text( (xj-0.1), (yj-0.1), (zj), Acrn{jjj},'FontSize',10); % Label the pt.
        elseif A(ii,jjj)>=0.4; %  highlight Strongest Out links 
            %w5k=jjj
          line([xi; xj], [yi; yj], [zi; zj],'Color', nodecolors(jjj,:),'LineWidth',3)
          %text( (xj-0.1), (yj-0.1), (zj), Acrn{jjj},'FontSize',10); % Label the pt.
        end
        arrowhead([xi, yi, zi],[xj, yj,zj], nodecolors(jjj,:), transp+0.2);  % fix the arrowhead
      % % % %  % % % %
      % omit 2nd NN - OUT links
      
     end  % if - in this time slice
    
    end % j-loop
     %missing % to look up Atlas coords
     %wts = round(A(ii,missing))  % check the weights of missing Out links
     %[ii length(missing)]
     
 prevT=contactTimes(l); % update startT for time slice
 drawnow
 %pause
end % l-timeSlice-loop

%clear i j jjj jlist k kk kkk wts  Xc Yc Zc ThisList prevT SecondNN % ii jj kk 
hold off  
  


% > >   >   >   >   >   >   >   >  IN      IN     IN
% > >   >   >   >   >   >   >   >  IN      IN     IN  time ordered i-j links  +k    +k 2nd NN
%% 12. first build  i-j & i-j-k IN lists;  draw later
  %  use code: 10a PLOT All nodes:  explore SINGLE Region /  Module; OUT Links  : base plot 
  % use dr as proxy for dt (assume, say, v-signal = 1 m/sec = 1 mm/mSec)
tau=0.5;  % time slice size
contactTimes = [0.5:tau:25.0];   % ascending vector of all possible contact times
tauS=2.0;  % est of Synaptic delay (mS, or mm for v=1m/s)
distijk= zeros(1,5); % first row: set up array [i, j, k, tij, tijk]
fprintf(' >>>>  MB, build list of IN Links [at #12.1] \n')

% picksources = these Nodes!  %  group of sensory outputs
%ThisList=[107, 52];  % Vis-1, Aud-1 
%ThisList=[29, 1, 30];  % SS / SSp
%ThisList=[13, 14];  % A2a,b (Cing)
ThisList=[2, 32, 48, 33, 45, 46, 47, 31];  % DlpFC group, in Ant-Post order

CutWt=5.0e-5; %CutWt=78.0 % CutWt=2.5085 %  % cut off wt for link wts: 
fprintf(' >>>>  Marmoset, DlpFC, build list of IN Links, wt>= %4.1f [at #11b] \n', CutWt)

% pick this Node!
totLinks=0; irow=1;  %  accumulate count
%OUT
jj=[]; alljj=[]; dxs=[]; dys=[]; dzs=[]; drs=[]; % arrays to accumulate results
%
for i=1:length(ThisList)
    ii=ThisList(i);
fprintf(' %4.0f %s %s %4.0f \n', ii, Acrn{ii}, Lobes{ii}, idx(ii) )
Xc=NodeCoord(ii,1); Yc=NodeCoord(ii,2);
Zc=NodeCoord(ii,3);   % coords of sphr centre (mm);  L side only so far
%IN
   jj = find(A(:,ii)>=CutWt)'; %  find [strong] In links (To "j" ~ here its ii)
   %jj = find( (A(:,ii)<=CutWt) & A(:,ii)>1 )'; %  find only weak In links
   totLinks=totLinks+1; % just count this i-j link here
   %wts = round(A(jj,ii))'  % check the IN wts
   dxs=abs(NodeCoord(jj,1)-Xc); % L-M dist % nb. need to use actual x-coord
   dys=abs(NodeCoord(jj,2)-Yc); % A-P dist;  col vec: i-jList
   dzs=abs(NodeCoord(jj,3)-Zc); % D-V
    %drs=sqrt(dxs.^2+dys.^2+dzs.^2);
   thedr= sqrt(dxs.^2+dys.^2+dzs.^2); % ColVec
   thisdr=[ii*ones(length(jj), 1), jj', thedr]; % assemble mx3 ColVec of [i j dr]
   drs= vertcat(drs, thisdr);  % nb vert cat. of mx3 arrays
   alljj = [alljj jj]; % add new 1-st NN
end   
   jj=alljj;
   [dsort, dindex] = sort(drs); % sorted (by dr, 3rd col) list & get new index
   dijrsort=drs(dindex,:); % reassemble the i-j-dr list
    %figure; hist(drs,50) % debug
    %dindex(1:10)'
    %dsort(1:10)'
   
prevT=0;  % for ( , ] timeSlice
 %fprintf('\n  Node#  Module  wt(j-i) dr(mm) Acrn Region  \n') % header
fprintf('\n  Node#  Module  wt(j-i) Flow(j-i)  dr(mm) Acrn Region  \n') % another header
for l=1:length(contactTimes)  % step through time slices; 
    %fprintf(' %10.0f %8.2f \n', l, contactTimes(l)) % debug
    for j=1:length(jj)
        jjj= jj(dindex(j));  % start with closest NN;  nb dindex maps sorted drs to orig find list
       if ( dsort(j) <= contactTimes(l) & dsort(j) > prevT ) % limit to this ( , ] timeSlice only
       %fprintf(' %6.0f %5.0f %8.0f %7.1f   %s %s \n', jjj, idx(jjj), A(jjj,ii), dsort(j), Acrn{jjj}, Region{jjj}(1:4) )
       fprintf(' %6.0f %5.0f %8.1f %9.4e %5.1f   %s %s \n', jjj, idx(jjj), A(jjj,ii), FlowWt(jjj,ii), dsort(j), Acrn{jjj}, Lobes{jjj} )
      
         xj=NodeCoord(jjj,1); yj=NodeCoord(jjj,2); zj=NodeCoord(jjj,3);
         %[xj, yj, zj]  % debug
      % % % %  % % % %
      % follow 2nd NN - IN links
        kk = find(A(:,jjj)>=CutWt)'; %  find only strong Out links from the j-node:
        %kk = find( (A(:,jjj)<=CutWt) & A(:,jjj)>1 )'; %  find only weak In links
        %out=length(kk) %debug
        totLinks=totLinks+length(kk);
        %wts = round(A(kk,jjj))'  % check the weights (in)
        dxjk=abs(NodeCoord(kk,1)-NodeCoord(jjj,1)); %col vec % nb. need to use actual x-coord, not "idx"
        dyjk=abs(NodeCoord(kk,2)-NodeCoord(jjj,2)); % 
        dzjk=abs(NodeCoord(kk,3)-NodeCoord(jjj,3)); drjk=sqrt(dxjk.^2+dyjk.^2+dzjk.^2); %clear dxs dys dzs 
        %drs'
        [drjksort, drjkindex] = sort(drjk); % sorted list; closest/smalest first & get new index
        kList=kk(drjkindex);
        for k=1:length(kList) % closest first
          kkk= kList(k);   %kkk= kk(kList(k)) 
          tik=dsort(j)+tauS+drjksort(k); % tij + t-synapse + tjk
          distijk(irow,:)= [ii jjj kkk dsort(j) tik]; % accumulate info
          irow=irow+1;  % need a row for every k
        end  % of k-loop
      
      end  % if - time slice
      
    end % of j-loop
 
 prevT=contactTimes(l); % update startT for time slice

end % of timeSlice-loop

tikCol=distijk(:,5);
[tsort, tindex] = sort(tikCol); % sorted list & get new index
distijkSort=distijk(tindex,:); % now sort the whole array, by incrs tijk  
clear tsort tindex tikCol
clear i j jjj jlist k kkk  Xc Yc Zc ThisList prevT % ii jj kk missing 
totLinks

%% 12.1 Draw:  Now, loop over i-j IN links; then scan j-k list in progressive t-slices
tikCol=distijk(:,5); % nb 2mS synaptic delay
[tsort, tindex] = sort(tikCol); % sorted list & get new index
distijkSort=distijk(tindex,:); % now sort the whole array, by incrs tijk  
clear tsort tindex tikCol
%
fprintf(' >>>>  MB, hiDegOUT: IN Links, time ordered; in MBtemporalNetV2.m (at #12.1)  \n')
figure('position',figposition, 'units','pixels'); hold on; 
set(gcf,'Renderer','OpenGL');
 %axis vis3d  % freeze aspect ratio (during rot)
axis equal  % get aspect ratio correct
 %set(gca, 'Projection', 'perspective');  % orthographic vs perspective view 
%view(-70, 40);  % 
view(110, 20);
  % use 4.1  plot saved surface: approx to plial: at Isosuf=rface of brain-vol: 50, in mm-space
  %patch(FVmm, 'FaceColor', [0.7 0.7 0.7],'FaceAlpha', 0.01,'EdgeColor','none'); % for 1/2
% set mid-transverse(x,y)-plane slice : as reference in 3D plots; min_x etc set above, #1.4.1  
surf([min_x max_x],[min_y max_y],repmat(imgzpositionmm, [2 2]),...
    midslice','facecolor','texture', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.0) % ng. no Edge (box)) 
colormap(gray) % default [64?]

title('Mouse Brain (mm): Area   VTA: In Links:top 40pc  ','FontSize',16,'FontWeight','bold')
text( 2.0, -2.5, 2, ' VTA: In Links: top 40pc  ','FontSize',16,'FontWeight','bold'); % for zoom-in label
 %text( 2.0, -2.5, 0, 'MPO: In Links: top 40pc','FontSize',16,'FontWeight','bold'); reposition, for zoom-n
%text( 1.0, -2.5, 0.25, 't < 4 mS  ','FontSize',16)  % Or, all t 
plot3(0,0,0, '+k','MarkerSize', 42); % Mark origin
line([0 0], [-0.3 0.3], [0 0], 'Color', 'k') % finish the 3D cross
line([0 0], [0 0], [-0.3 0.3], 'Color', 'k')
text( 0, -0.5, 0.1, 'P','FontSize',16); 
text( 0, 0.5, 0, 'A','FontSize',16); % Label directions
text( 0, -0.1, 0.4, 'D','FontSize',16);
text( 0, 0, -0.4, 'V','FontSize',16);
xlabel('x (mm)  Lat - Med ','FontSize',14)
ylabel('y (mm)  Antr - Post  ','FontSize',14)  % nb. x-y swapped above - in cols of FV.vertices
zlabel('z (mm)  Ventr - Dors ','FontSize',14)
axis([-5 5 -8 4  -6.7 1]); % rely on min/max val of y-mm & extend y to show features

% overlay Region Coords on plot of Surface(mm)
% Mark origin
plot3(0,0,0, '+k','MarkerSize', 45);  % Bregma
% Plot all nodes: coords of sphr centre (mm);  L side only so far 
  % place "other" points at their y-Region; and on mid-LHS {find their coords later
nocoord=0; % track missing coords :: 98 missed
for i=1:n
    if(abs(CoordSize(i,1)) <= 0.1 &  abs(CoordSize(i,3)) <= 0.1) % "others have 0 entry for x & z
        % [i CoordSize(i,2)] % check ok
     nocoord=nocoord+1;
     plot3(0, CoordSize(i,2), -3.5, '.', 'Color',nodecolors(i,:),'MarkerSize', 15)
     %text( (0-0.1), (CoordSize(i,2)-0.1), (-3.5), num2str(i),'FontSize',8); % number the pt.
     %pause % debug
    else
     plot3(CoordSize(i,1), CoordSize(i,2), CoordSize(i,3), '.', 'Color', nodecolors(i,:),'MarkerSize', 15) % orig coord
     %text( (CoordSize(i,1)-0.1), (CoordSize(i,2)-0.1), (CoordSize(i,3)-0), Acrn{i},'FontSize',10); % Label (Acronym) the pt.
    end
end

% overlay 2nd source (i) from here
% do i-j outer Loop; draw i-j Link once only; then scan j-k List
% find min-T for first 2ndNN j-k Link:
minTk=distijkSort(1,5);  % this list ordered by t-ijk: this is the first t-ijk to look for
tau=0.5;  % time slice size % reset, to be sure
% when t-ij >= minTk, then search for j-k links:

% Plot selected Node[s] & their OUT Links (>10, 100)
% plot a  sphere at CM of ea. region, in top-10,20 list  % colour by module ID
transp=0.4;  % for Alpha (line drawing): hi flow nodes (%1-10) heavier

fprintf(' %4.0f %s %s \n', ii, Acrn{ii}, Region{ii} )
Xc=CoordSize(ii,1); Yc=CoordSize(ii,2);
Zc=CoordSize(ii,3);   % coords of sphr centre (mm);  L side only so far
  %[Xc, Yc, Zc]  % debug
  if ( (abs(Xc)<=0.02) & (abs(Zc) <= 0.11) )  % check if coord is missing?
     Xc =0; Zc=-3.5;  % default (ie to "others"): place them at base (for now)
     Yc=CoordSize(ii,2);  
     missing =[missing ii]; % note: still to look up coords for these
  else
     Xc=CoordSize(ii,1); Yc=CoordSize(ii,2); Zc=CoordSize(ii,3) ; % use known coords
  end
  plot3(Xc, Yc, Zc, '.', 'Color', nodecolors(ii,:),'MarkerSize', 25) % orig coord, i-Node
     %text( (CoordSize(ii,1)+0.1), (CoordSize(ii,2)-0.2), (CoordSize(ii,3)+0.1), Region{ii}(1:4),'FontSize',10, 'Color', [0.5 0.5 0.5]); % L side % label original region
  text( (Xc+0.1), (Yc+0.1), (Zc-0), Acrn{ii},'FontSize',10); % R side % label node  
   %plot3(CoordSize(ii,1), CoordSize(ii,2), CoordSize(ii,3), '.', 'Color', nodecolors(ii,:),'MarkerSize', 20) % orig coord
   %text( (CoordSize(ii,1)+0.1), (CoordSize(ii,2)-0.2), (CoordSize(ii,3)+0.1), Region{ii}(1:4),'FontSize',10, 'Color', [0.5 0.5 0.5]); % L side, label original region
   %text( (CoordSize(ii,1)-0.1), (CoordSize(ii,2)-0.1), (CoordSize(ii,3)-0), Acrn{ii},'FontSize',10); % L side % label node  
%In :: scan the list found above
  %jj = find(A(:,ii)>=CutWt)'; %  find [strong] In links (1st NN To i): checked against DegOut
   %wts = round(A(jj,ii))';  % check the weights
   dxs=abs(CoordSize(jj,1)-Xc); % L-M dist % nb. need to use actual x-coord
   dys=abs(CoordSize(jj,2)-Yc); % A-P dist;  col vec: i-jList
   dzs=abs(CoordSize(jj,3)-Zc); % D-V
   drs=sqrt(dxs.^2+dys.^2+dzs.^2);
   [dsort, dindex] = sort(drs); % sorted list & get new index
    %figure; hist(drs,50) % debug
    %dindex(1:10)'
    %dsort(1:10)'
%jj(dindex) %debug % sorted list of j-1stNN
tmpj=jj(dindex'); %debug

[nLinks ~] = size(distijk);
prevT=0;  % for ( , ] timeSlice
kListRow=1; % for scanning j-k (t-ijk) links
fprintf('\n   Node#  Module  wt(j-i) dr(mm) Acrn Region  \n') % header
  % fprintf('\n        step#   t-slice  \n') % header
%prevT=contactTimes(11-1)  % Debug
for l=1:length(contactTimes)  % step through successive time slices; 
    %fprintf('\n  step, TimeWindow: %6.0f %8.2f %8.2f ', l, prevT, contactTimes(l)) % debug
    %fprintf('\n   Node#  Module  wt(j-i) dr(mm) Acrn Region  \n') % header
    for j=1:length(jj)
        jjj= jj(dindex(j));  % start with closest NN;  nb dindex maps sorted drs to orig find list
        %fprintf(' step %2.0f work on j: %3.0f d-ij %7.1f \n', j, jjj, dsort(j) ); % debug
      if ( dsort(j) <= contactTimes(l) & dsort(j) > prevT ) % limit to this ( , ] timeSlice only
      fprintf(' %6.0f %5.0f %8.0f %7.2f   %s %s \n', jjj, idx(jjj), A(jjj,ii), dsort(j), Acrn{jjj}, Region{jjj}(1:4) )             
         xj=CoordSize(jjj,1); yj=CoordSize(jjj,2); zj=CoordSize(jjj,3);
         %[xj, yj, zj]  % debug
          if ( (abs(xj)<=0.02) & (abs(zj) <= 0.11) )  % check if coord is missing?
            xj =0; zj=-3.5;  % default (ie to "others"): place them at base (for now)
            yj=CoordSize(jjj,2);  
            missing =[missing jjj]; % note: still to look up coords for these
          else
            xj=CoordSize(jjj,1); yj=CoordSize(jjj,2); zj=CoordSize(jjj,3); % use known coords
          end
          %[xj, yj, zj] % debug 
      plot3(xj, yj, zj, '.','Color', nodecolors(jjj,:), 'MarkerSize', 15);  % plot & label this j-Node point
      text( (xj+0.1), (yj+0.1), (zj-0), Acrn{jjj},'FontSize',10); % Label the pt.
      text( (xj+0.1), (yj-0.15), (zj+0.15), Region{jjj}(1:4),'FontSize',10, 'Color', [0.5 0.5 0.5]); % Label the target Region (faint)
      % colour Out i ->j link by the Target(jjj) node
      patch([Xc; xj; Xc], [Yc; yj; Yc], [Zc; zj; Zc], 'k', 'EdgeColor', ...
               nodecolors(jjj,:), 'LineSmoothing','on', 'EdgeAlpha', transp, 'FaceColor', 'none'); % stronger colour / alpha 0.05-0.2
      arrowhead([xj, yj,zj],[Xc, Yc, Zc], nodecolors(jjj,:), transp);  % add an In (i <- j) arrowhead 85% along line
      % Strong links:
        if A(jjj,ii)>=1000; %  highlight Stronger Out links; only label these points
          line([Xc; xj], [Yc; yj], [Zc; zj],'Color', nodecolors(ii,:),'LineWidth',2)
          %text( (xj-0.1), (yj-0.1), (zj), Acrn{jjj},'FontSize',10); % Label the pt.
             %wt1k=jjj
        elseif A(jjj,ii)>=5000; %  highlight Strongest Out links 
            %wt5k=jjj
          line([Xc; xj], [Yc; yj], [Zc; zj],'Color', nodecolors(ii,:),'LineWidth',3)
          %text( (xj-0.1), (yj-0.1), (zj), Acrn{jjj},'FontSize',10); % Label the pt.
        end
        arrowhead([xj, yj,zj], [Xc, Yc, Zc], nodecolors(jjj,:), transp);  % fix the arrowhead
      % % % %  % % % %
      % follow 2nd NN - OUT links
        % ?? if contactTimes(l) >= minTk  % any k-links to show yet?
      if dsort(j) >= minTk  % any k-links at/near this i-j time ?
           %fprintf(' +k links exist at this timeSlice \n') % debug
          for idummy=1:50 % "dummy" loop to scan j-k list
           %fprintf('      > link between:  %6.0f %5.0f t-ijk %4.2f \n', distijkSort(kListRow, 2), distijkSort(kListRow, 3), distijkSort(kListRow, 5) ); % debug
          thisJ=distijkSort(kListRow, 2); % check we have correct row?
          kkk=distijkSort(kListRow, 3); % select this k
              % ?? if (distijkSort(kListRow, 5) > dsort(j) & dsort(j) <= dsort(j)+ tau)  % check t-ijk in this slice
           if distijkSort(kListRow, 5) <= dsort(j)     % debug 
            fprintf('        & draw link:  %4.0f %4.0f %s %s \n', thisJ, kkk, Acrn{kkk}, Region{kkk}(1:4)); % debug
            xk=CoordSize(kkk,1); yk=CoordSize(kkk,2); zk=CoordSize(kkk,3); % & draw this link:
            xj=CoordSize(thisJ,1); yj=CoordSize(thisJ,2); zj=CoordSize(thisJ,3); % nb. not same "j" as above
            %[xk, yk, zk]  % debug
            text( (xk+0.1), (yk+0.1), (zk-0), Acrn{kkk},'FontSize',10);
            %text( (xk+0.1), (yk-0.15), (zk+0.1), Region{kkk}(1:4),'FontSize',10, 'Color', [0.5 0.5 0.5]); 
             %text( (CoordSize(kkk,1)-0.1), (CoordSize(kkk,2)-0.15), (CoordSize(kkk,3)-0.15), num2str(kkk),'FontSize',9); % shift slightly
            % draw j->k link
            patch([xj; xk; xj], [yj; yk; yj], [zj; zk; zj], 'k', 'EdgeColor', ...
               nodecolors(kkk,:), 'LineSmoothing','on', 'EdgeAlpha', transp, 'FaceColor', 'none'); 
            arrowhead([xk, yk, zk],[xj, yj,zj], nodecolors(kkk,:), transp);  % add an In (j <- k) arrowhead
            kListRow=kListRow+1;  % increment row: go to next row in j-k list
            %pause
           else
                 % ?? kListRow=kListRow+1;  % also need to increment row here
               break % skip out & go on to next j, in outer loop
           end % if
          end % for loop over list 
      end % if - search for k at this time?
      
     end  % if - in this time slice
    
    end % j-loop
     %missing % to look up Atlas coords
     %wts = round(A(ii,missing))  % check the weights of missing Out links
     %[ii length(missing)]
     
 prevT=contactTimes(l); % update startT for time slice
 pause
end % l-timeSlice-loop

clear i j  jjj jlist k kk kkk wts wts_missing  Xc Yc Zc ThisList prevT SecondNN % ii jj kk missing 
hold off  

  