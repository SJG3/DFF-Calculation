%% Functional network construction and cell ranking for in-vitro neuronal populations
% Zac Bowen 2021 modified SJG 2021-2023
close all force
clear 
clc
%% Manual input parameters -- worked well with the example data
expectedNeuronRadiusPix = 5; % Radius of circular ROI on each neuron
dffBaselineSize = 300; % window size IN FRAMES for baseline fluo estimation
dffBaselinePercent = 75; % bottom percentile of fluo to use for baseline estimate
minpkheight = 0.2; %3 seems to be a good value
% fps = ;
% total_time = ;
%% Load image stack -- change this code to load IMG matrix as specified below
%added by SJG
    [fileX,pathX] = uigetfile('*.tif');
    disp("Opening File:" + fileX);
    im_path = fullfile(pathX,fileX);
    info = imfinfo(im_path);
    FileTif= im_path;

    %save the Tif into an array
    tic;
    TifLink = Tiff(FileTif, 'r');
    for frame = 1:numel(info)
        TifLink.setDirectory(frame);
        IMG(:,:,frame) = im2double(TifLink.read());
    end
    TifLink.close();
clc;
% - - - - - - - - - - - 

% load('RENCellCalbryteDemo2.mat','imAll');
% 
% IMG = imAll; %IMG needs to be (dimX,dimY,frames)
% clear imAll
% meanIMG = imread('avgprojsample.jpg');
 meanIMG = mean(IMG,3); % preferably use this (whole image stack)
 

% added by SJG
% meansubIMG = double(IMG) - meanIMG;
% minIMG = min(IMG,3); 
% minsubIMG = (double(IMG) - double(minIMG));
% stdminsubIMG = std(minsubIMG,0,3); 

stdIMG = std(double(IMG),0,3);
pause(0.001); 
clc
disp('Loading Image Complete'); 
%% (SJG added) measure varriance of image
%https://en.wikipedia.org/wiki/Coefficient_of_variation#Similar_ratios

cvIMG = ((stdIMG ./ meanIMG).^1); %SJG Added  = Coefficient of variation or relative standard deviation
v2meanIMG = (stdIMG.^2)./ meanIMG; %variance to mean ratio
recipIMG = meanIMG ./ stdIMG;  % recipricol ratio (SNR)
recip_sqIMG = (meanIMG ./ stdIMG).^2; % recipricol ratio (SNR) squared
effIMG = cvIMG.^2 ; %efficicency 

figure(467582);
subplot(2,3,1); imagesc(meanIMG); axis('image'); title(["Mean Image";"\mu"]); 
subplot(2,3,2); imagesc(cvIMG); axis('image'); title(["Coefficient of variation (CV)";"relative standard deviation(RSD)";"\sigma/\mu"]); caxis([0 1]);
subplot(2,3,3); imagesc(effIMG); axis('image'); title(["SJG found ";"Efficiency";"(\sigma/\mu)^2"]);  caxis([0 1]);
subplot(2,3,4); imagesc(v2meanIMG); axis('image'); title(["Variance-to-mean ratio";"(relative variance)";"\sigma^2/\mu"]); 
subplot(2,3,5); imagesc(recipIMG); axis('image'); title(["reciprocal ratio";"\mu/\sigma"]);
subplot(2,3,6); imagesc(recip_sqIMG); axis('image'); title(["reciprocal ratio squared";"\mu/\sigma^2"]);
colormap(turbo);
%% image for background discrimination (SJG) 
% IMG_for_auto = meanIMG;
% 
% IMG_for_auto = imgaussfilt(IMG_for_auto,2);
% 
% se = strel('disk',32);
% bg_detect = imopen(IMG_for_auto,se);
% IMG_for_auto = IMG_for_auto - bg_detect;
% 
% vals = prctile(IMG_for_auto(:),[10 99.9]);
% IMG_for_auto = imadjust(IMG_for_auto, [vals(1) vals(2)],[0 1]);
% IMG_for_auto = meanIMG;
% 
% figure(1352); 
% subplot(1,2,1); imagesc(meanIMG); axis image; colormap( [0,0.6, 0;hot] );
% subplot(1,2,2); imagesc(IMG_for_auto); axis image;
%% identify cell centers (manual, load, automated[added by SJG using Diff of Gausian])

xc = [];
yc = [];
igg = v2meanIMG ;%meanIMG; % try meanIMG, cvIMG, effIMG, or stdIMG, v2meanIMG

answer = questdlg('Cell selection:', ... % Prompt
	'Get cell coordinates', ... % name figure window
	'Manually click','Load from file','Automated (DoG)',... % define options
    'Manually click'); % set default option
% Handle response
switch answer
    case 'Manually click'
%        figure; imagesc(meanIMG); colormap('gray'); axis('square')%orginal
        figure; imagesc(igg); colormap('turbo'); axis('square') %SJG changed
        
  
        fprintf('Click on Neuron Centers, press enter when done\n')
        [xc, yc] = getpts; %  manually select centers of the neurons

%SJG ADDED, does what getpts does
%         userStopped = false; 
%         pointhandles = gobjects(); 
%         while ~userStopped
%             a = drawpoint('MarkerSize',10,'Deletable',true); 
%             if ~isvalid(a) || isempty(a.Position)
%                 % End the loop
%                 userStopped = true;
%             else
%                 % store point object handle
%                 pointhandles(end+1) = a;
%                 xc(length(pointhandles)-1,1) = a.Position(1);
%                 yc(length(pointhandles)-1,1) = a.Position(1);
%             end
%         end
%         clear pointhandles
% --------
        
%         save(['CellDefFile_' date],'xc','yc'); %original 
        save([pathX 'CellDefFile_' date],'xc','yc'); %SJG changed
        close
        fprintf('...cell definition file saved!\n')
    case 'Load from file'
        uiopen
        
    case 'Automated (Otsu)' %Added by SJG
        level = graythresh(igg);
        BW_ = imbinarize(igg,level*0.5);
        BW_ = bwareafilt(BW_,[15 1000]);
        BW_l = bwlabel(BW_);
        imshowpair(igg,BW_,'montage');
        colormap(turbo);

        BW_props = regionprops(BW_,'centroid','circularity');
        for i = 1:length(BW_props)
            if BW_props(i).Circularity >=0.5
                xc(i,1) = BW_props(i).Centroid(1);
                yc(i,1) = BW_props(i).Centroid(2);
            else
                xc(i,1) = NaN;
                yc(i,1) = NaN;
            end
        end
        xc = rmmissing(xc);
        yc = rmmissing(yc);
        save([pathX 'CellDefFile_AutoOtsu' date],'xc','yc'); %SJG changed
        
    case 'Automated (DoG)' %difference of gausian 
        IMG_start = meanIMG;
        IMG_g = imgaussfilt(IMG_start,1); 
        bg = imopen(IMG_g,strel('disk',15));
        
        IMG_use = imsubtract(IMG_g,bg);

        sigma1 = 7; %This value should be ~ the radius of a cell/object; should be smaller than sigma2
        sigma2 = 10; %This value should be < sigma1
        G1_Img = imgaussfilt(IMG_use,sigma1);
        G2_Img = imgaussfilt(IMG_use,sigma2);
       
        D = mat2gray(G1_Img - G2_Img);
        BN = D > graythresh(D) * 0.570;
        
        SE = strel('disk',3);
        BN_rm = imregionalmax(D,8) .* BN > 0;
        BN_rm_d = imdilate(BN_rm,SE);
        BN_lab = bwlabel(BN_rm_d);
        
        figure(1); 
        subplot(2,2,1); imshow(imfuse(meanIMG,BN_rm_d));
%         subplot(2,2,2); imshow(imfuse(imadjust(IMG_for_auto, [0 0.2]),BN_rm_d));
        subplot(2,2,2); imshow(imfuse(imadjust(meanIMG, [prctile(meanIMG(:),0.5) prctile(meanIMG(:),70)]),BN_rm_d));

        subplot(2,2,3); imshow(imadjust(D));
        subplot(2,2,4); imshow(BN);

        
        BW_props = regionprops(BN_rm_d,'centroid','circularity');
        for i = 1:length(BW_props)
            if BW_props(i).Circularity >=0.5
                xc(i,1) = BW_props(i).Centroid(1);
                yc(i,1) = BW_props(i).Centroid(2);
            else
                xc(i,1) = NaN;
                yc(i,1) = NaN;
            end
        end
        xc = rmmissing(xc);
        yc = rmmissing(yc);
        save([pathX 'CellDefFile_AutoDoG' date],'xc','yc'); %SJG changed
end
%% Extract fluo over time from a circle around selected cell coords
% [rawFluo,roiBW] = extractFluo(IMG,xc,yc,expectedNeuronRadiusPix);
[rawFluo,roiBW] = extractFluo(IMG,xc,yc,expectedNeuronRadiusPix);

ROImap = sum(cat(3,roiBW{:}),3); %in case you want to see all ROIs together
%% dF/F calculation

% Normalizing raw F to get rid of negative vals before dF/F calc % Jul2020
 normF = rawFluo ./ max(rawFluo(:));

% Performs moving average baseline subtraction to compute dF/F
% dffBaselineSize = round(exptVars.frameRate * winSizeSeconds);
tic;
DFF = slideWinSub(normF',dffBaselineSize,dffBaselinePercent);

% ---- SJG 2-1-2023
DFF_active = zeros(height(DFF),width(DFF));
for i = 1:height(DFF)
%     [pks,locs] = findpeaks(DFF(i,:),'MinPeakHeight',minpkheight);

    prom = std(DFF(i,:))*3; 
    [pks,locs] = findpeaks((DFF(i,:)'),'MinPeakProminence',prom,'MinPeakHeight',0.2,'Annotate','extents');

    if length(locs) < 3 %== 0
         DFF_active(i,:) = NaN;
    else
        DFF_active(i,:) = DFF(i,:);
    end
    
    clear pks
    clear locs
end
DFF_active(any(isnan(DFF_active), 2), :) = [];
%----

c_max = size(xc,1);
toc;
%%  TEST look for peaks
fps = 0.5;
spike_width = 5;

figure(45432); 
for nom = 1:size(DFF,1)%[1:10:size(DFF,1)]
    prom = std(DFF(nom,:))*3; 
    findpeaks((DFF(nom,:)'),'MinPeakProminence',prom,'MinPeakHeight',0.2,'Annotate','extents');
    title("Object:"+nom+" at prom:"+prom);
   
    ylim([- max(DFF(:)) ,  max(DFF(:))]);
%     hold on
pause(0.4);
end 
%%      Correlation Coeff (SJG test)
%using Pearson Correlation Coefficient
cor_Mat = corrcoef (DFF');
%%      Cross Corr (SJG test)
% xcor_Mat = xcorr(DFF);

% cross correlation
for c1 = 1:c_max  
    for c2 = 1:c_max 
       [corrs,lags] = xcorr(DFF(c1,:),DFF(c2,:),750,'normalized');
%         indx = length(1:size(IMG,3));
        indx = abs(corrs) == max(abs(corrs)); %will get crosscor for max value
        XCorMat(c1,c2) = corrs(indx);
        XCorLagMat(c1,c2) = lags(indx);
    end
end 
%%      Distance Matrix (SJG test)
for uno  = [1:c_max]
    for dos = [1:c_max]
        xdistmat(uno,dos) = pdist([xc(uno);xc(dos)],'euclidean');
        ydistmat(uno,dos) = pdist([yc(uno);yc(dos)],'euclidean');
    end
end
distMat = hypot(xdistmat,ydistmat);
distMatUp = distMat; 
ii=ones(size(distMatUp));
idx=triu(ii,1);
distMatUp(~idx)=nan;
%% Plot DFF and Correlation Matrix
yes = 1; %makes color map for kymo
if yes == 1  
    resol1 = 100;
    s_f = 1;
    custom0 = [linspace(0,0,resol1)',linspace(1,1,resol1)',linspace(1,0,resol1)'];
    custom1 = [linspace(0,0,resol1)',linspace(0,1,resol1)',linspace(1,1,resol1)'];
    custom2 = [linspace(0.75,0,resol1*s_f)',linspace(0.75,0,resol1*s_f)',linspace(0.75,1,resol1*s_f)'];
    custom3 = [linspace(1,0.75,resol1*s_f)',linspace(0,0.75,resol1*s_f)',linspace(0,0.75,resol1*s_f)'];
    custom4 = [linspace(1,1,resol1)',linspace(1,0,resol1)',linspace(0,0,resol1)'];
    custom5 = [linspace(1,1,resol1)',linspace(0,1,resol1)',linspace(1,0,resol1)'];
    cmaprange = flipud([custom3;custom2]);
end

fig_n = 126;
figure(fig_n);
    subplot(2,3,1:3); plot(DFF','LineWidth',0.01); %colormap(fliplr(cmaprange));
    colororder(lines(height(DFF)));
    xlabel('Time (frames)'); ylabel('\DeltaF/F');


    subplot(2,3,4:6); imagesc(DFF); colormap(fliplr(cmaprange));
    crng = (max(DFF(:))); %axis('image');
    caxis([-crng crng]);
    xlabel('Time (frames)'); ylabel('Cell Index#'); colorbar;


figure(fig_n+1);
    subplot(2,2,1); imagesc(cor_Mat); colormap(fliplr(cmaprange));
    title(["Pearson Correlation Matrix"]);
    axis ('square'); caxis([-1,1]);
    xlabel('Cell Index#'); ylabel('Cell Index#'); colorbar;

    subplot(2,2,2); imagesc(XCorMat); colormap(fliplr(cmaprange));
    title("Cross Corr Matrix @ Max");
    axis ('square'); caxis([-1,1]);
    xlabel('Cell Index#'); ylabel('Cell Index#'); colorbar;

    subplot(2,2,3); imagesc(XCorLagMat); colormap(fliplr(cmaprange));
    title("Cross Corr Matrix Lag Values");
    axis ('square');
    xlabel('Cell Index#'); ylabel('Cell Index#'); colorbar;

    subplot(2,2,4); histogram(XCorLagMat,25); colormap(fliplr(cmaprange));
    title("Cross Corr Matrix Lag Values");
    axis ('square');
    xlabel('Count #'); ylabel('Lag Value#'); colorbar;


figure(fig_n+2);
    subplot(1,2,1); imagesc(meanIMG); axis image;

    subplot(1,2,2); imagesc(BN_lab); axis image; colormap ([0.15,0.15,0.15;turbo]);
%%      get triangles
tri_up_NaN = triu(ones(size(xc,1)),1);
tri_up_NaN(tri_up_NaN==0) = NaN;

cor_Mat_v = cor_Mat .* tri_up_NaN;
cor_Mat_v = reshape(cor_Mat_v,[],1);
cor_Mat_v = cor_Mat_v(~isnan(cor_Mat_v))';

distMat_v = distMat .* tri_up_NaN;
distMat_v = reshape(distMat_v,[],1);
distMat_v = distMat_v(~isnan(distMat_v))';

XCorMat_v = XCorMat .* tri_up_NaN;
XCorMat_v = reshape(XCorMat_v,[],1);
XCorMat_v = XCorMat_v(~isnan(XCorMat_v))';

XCorLagMat_v = XCorLagMat .* tri_up_NaN;
XCorLagMat_v = reshape(XCorLagMat_v,[],1);
XCorLagMat_v = XCorLagMat_v(~isnan(XCorLagMat_v))';
%% Plot Distance v correlation scatters (TEST) (using online found scatter kernal density estimate code)
figure(475); 
%distance vs correlation
colormap(jet);
subplot(1,3,1); scatter_kde (distMat_v',cor_Mat_v','filled','MarkerFaceAlpha',0.5);
title('Distance vs Correlation');
xlabel('Distance(px)'); ylabel('Correlation Value');
xlim([0,max(distMat_v)]); ylim([-1,1]);

subplot(1,3,2); scatter_kde (distMat_v',XCorMat_v','filled','MarkerFaceAlpha',0.5);
title('Distance vs CrossCorr');
xlabel('Distance(px)'); ylabel('Correlation Value');
xlim([0,max(distMat_v)]); ylim([-1,1]);

subplot(1,3,3); scatter_kde (distMat_v',XCorLagMat_v','filled','MarkerFaceAlpha',0.5);
title('Distance vs CrossCorr');
xlabel('Distance(px)'); ylabel('Lag Value');
xlim([0,max(distMat_v)]); ylim([min(XCorLagMat_v),max(XCorLagMat_v)]);


% scatter_kde(data1,data2,'.');
%% Plot Distance v correlation scatters (TEST)(using desu density scattercloud)
figure(27592)

num_xbins = 10;
num_contour_regions = 5;
pv = 1;

subplot(1,3,1);
[h,C,C_n,xEdges,yEdges] = scattercloud(distMat_v',cor_Mat_v',100,1,10,10,jet(256),0,max(distMat_v),-1,1);
colorbar; xlim([0,max(distMat_v)]); ylim([-1,1]);
hold on
contour(xEdges,yEdges,conv2(C_n,1/81*ones(9,9),'same'),'LineWidth',2,'LevelList',[0:1/num_contour_regions:1],'LineColor','k');
mean_scatter = binned_mean(distMat_v',cor_Mat_v',num_xbins,0,max(distMat_v),pv);
plot(mean_scatter(:,1),mean_scatter(:,2),'m','LineWidth',2);
plot(mean_scatter(:,1),mean_scatter(:,3),'m--','LineWidth',2);
plot(mean_scatter(:,1),mean_scatter(:,4),'m--','LineWidth',2);
xlabel('Pairwise Distance (pixel)')
ylabel('Correlation')
set(gcf,'color','w');
title('Pairwise Correlation vs. Pairwise Distance')
hold off; 

subplot(1,3,2);
[h,C,C_n,xEdges,yEdges] = scattercloud(distMat_v',XCorMat_v',100,1,10,10,jet(256),0,max(distMat_v),-1,1);
colorbar; xlim([0,max(distMat_v)]); ylim([-1,1]);
hold on
contour(xEdges,yEdges,conv2(C_n,1/81*ones(9,9),'same'),'LineWidth',2,'LevelList',[0:1/num_contour_regions:1],'LineColor','k');
mean_scatter = binned_mean(distMat_v',XCorMat_v',num_xbins,0,max(distMat_v),pv);
plot(mean_scatter(:,1),mean_scatter(:,2),'m','LineWidth',2);
plot(mean_scatter(:,1),mean_scatter(:,3),'m--','LineWidth',2);
plot(mean_scatter(:,1),mean_scatter(:,4),'m--','LineWidth',2);
xlabel('Pairwise Distance (pixel)')
ylabel('Max Cross Correlation')
set(gcf,'color','w');
title('Pairwise Cross Correlation vs. Pairwise Distance')
hold off; 

subplot(1,3,3);
[h,C,C_n,xEdges,yEdges] = scattercloud(distMat_v',XCorLagMat_v',100,1,10,10,jet(256),0,max(distMat_v),min(XCorLagMat_v),max(XCorLagMat_v));
colorbar; xlim([0,max(distMat_v)]); ylim([min(XCorLagMat_v),max(XCorLagMat_v)]);
hold on
contour(xEdges,yEdges,conv2(C_n,1/81*ones(9,9),'same'),'LineWidth',2,'LevelList',[0:1/num_contour_regions:1],'LineColor','k');
mean_scatter = binned_mean(distMat_v',XCorLagMat_v',num_xbins,0,max(distMat_v),pv);
plot(mean_scatter(:,1),mean_scatter(:,2),'m','LineWidth',2);
plot(mean_scatter(:,1),mean_scatter(:,3),'m--','LineWidth',2);
plot(mean_scatter(:,1),mean_scatter(:,4),'m--','LineWidth',2);
xlabel('Pairwise Distance (pixel)')
ylabel('Max Cross Correlation Lag (frame)')
set(gcf,'color','w');
title('Max Cross Correlation Lag vs. Pairwise Distance')
hold off; 
%% area under curve  (SJG test)

AUC = trapz(DFF,2); 

figure;
subplot(1,2,1); imagesc(AUC); colormap(fliplr(cmaprange)); caxis([-max(AUC) max(AUC)]);
subplot(1,2,2); imagesc(DFF); colormap(fliplr(cmaprange));caxis([-crng crng]);

%% TEST histogram of correlation coeff
figure; histogram(cor_Mat_v,'BinWidth',0.05); 
 
%percentage of active cells ( = active cells/ total cells found)
size(DFF_active,1)/size(DFF,1)



%% FUNCTION: Fluo extraction function -- I left in a lot of commented out code in case
% we go back and want to define cell-by-cell custom ROIs based on fluo.
% Right now all ROIs are the same size and shape.
 function [rawFluo,roiBW] = extractFluo(IMG,xc,yc,expectedNeuronRadiusPix)


% narginchk(4,7)
% if ~exist('winSizeSeconds','var'); winSizeSeconds = 10; end % default to 10 second window
% if ~exist('percentBaselineSub','var'); percentBaselineSub = 50; end % default to 50% baseline subtraction

numNeurons = length(xc);

%% PREALLOCATE
traceExtract_start = tic;

% pcimg = cell (numNeurons , 13 , 360 );
% imgCrop = cell ( numNeurons , 27, 27 );
% imgCropNorm = cell (numNeurons, 28 , 28 );
% roiBoundaries = cell ( numNeurons , 360 , 3 );
% smRoiBoundaries = cell ( numNeurons , 360 , 3 );
% ROIOut = zeros ( 360 , 2 );

% 2-1-2023 SJG commented out and it seems to allow code to run tell Anna
ignor = 1;
if ignor == 0; 
    ROIxvOut = cell ( numNeurons , 360 , 1 );
    ROIyvOut = cell ( numNeurons , 360 , 1 );
    roiBW = cell ( numNeurons , size(IMG,1) , size(IMG,2));
end 
%% PREALLOCATE FLUO AND NPFLUO MATRICES
rawFluo = zeros( size(IMG,3) , numNeurons );

%% FIND THE BOUNDARIES OF CLICKED NEURONS
roiBounds = [deg2rad(1:360)' repmat(expectedNeuronRadiusPix,360,1)];
for pp = 1:numNeurons
    
%     % Find fluorescent ring using local peak finding on each neuron
%     xpt=xc(pp);
%     ypt=yc(pp);
%     imgCrop{pp} = imcrop(meanIMG,[xpt-15 ypt-15 31 31]);
%     imgCropNorm{pp} = (imgCrop{pp} - min(imgCrop{pp}(:)))  ./ (max(imgCrop{pp}(:)) - min(imgCrop{pp}(:)));
%     pcimg{pp} = imgpolarcoord (imgCropNorm{pp} );  % this comes from Matlab Central Download
%     
%     RingPks = zeros(size(pcimg{pp},2),1); %reset vals to zero
%     
%     tmpNeuron = pcimg{pp}; %iterate this for each selected neuron
%     
%     for cc = 1:size(tmpNeuron,2) % for every direction - find the inner part of the ring
%         
%         pkTmp = find(diff(tmpNeuron(:,cc)) == min(diff(tmpNeuron(:,cc)))); % DW07122015_changed to make this more robust - seems to be working right now - continue testing
%         
%         if ~isempty(pkTmp)
%             if length(pkTmp) > 1 %more than one pixel identified - grab the first one
%                 if pkTmp(1) < expectedNeuronRadiusPix && pkTmp(1) > 2
%                     RingPks(cc) = pkTmp(1);
%                 else
%                     RingPks(cc) = expectedNeuronRadiusPix;
%                 end
%             else
%                 if pkTmp < expectedNeuronRadiusPix && pkTmp(1) > 2
%                     RingPks(cc) = pkTmp;
%                 else
%                     RingPks(cc) = expectedNeuronRadiusPix;
%                 end
%             end
%         elseif cc == 1 % if it's the first direction and no peaks are found
%             RingPks(cc) = expectedNeuronRadiusPix;  %made this dependent on mag factor DW_02022015
%         else
%             RingPks(cc) = RingPks(cc-1);
%         end
%         ROIOut(cc,:) = [ deg2rad(cc)  RingPks(cc) ];
% 
%     end
%     roiBoundaries{pp} = [ ROIOut(:,1) ROIOut(:,2)]; % [PolarCoords (0-2Pi)     OuterRing]
%     smRoiBoundaries{pp} = [ ROIOut(:,1) smooth(ROIOut(:,2),10)]; % [PolarCoords (0-2pi)     OuterRing]
    
    % CREATE MASKS FOR ALL CLICKED ROIS, THEN SHOW THEM -- DW 11232015
    % renamed variable for consistency
%     ROIxvOut{pp} =  xpt + smRoiBoundaries{pp}(:,2) .* (cos(smRoiBoundaries{pp}(:,1))) ;
%     ROIyvOut{pp} =  ypt + smRoiBoundaries{pp}(:,2) .* (sin(smRoiBoundaries{pp}(:,1))) ;
    ROIxvOut{pp} =  xc(pp) + roiBounds(:,2) .* (cos(roiBounds(:,1))) ;
    ROIyvOut{pp} =  yc(pp) + roiBounds(:,2) .* (sin(roiBounds(:,1))) ;
    roiBW{pp} = poly2mask( ROIxvOut{pp} , ROIyvOut{pp} , size(IMG,1) , size(IMG,2));
end

%DW 11232015 - adjusted for inclusion of neuropil correction
% correct for overlapping ROIs (exclude from both)
disp('Adjusting ROI masks for overlap....');
tStartROICorr = tic;
AllMasksTMP =  sum ( cat ( 3 , roiBW{:} ) , 3 ); % first term of cat (i.e., '3') points to element-wise alignement/stacking of arrays
[oLapRoiY, oLapRoiX] = find( AllMasksTMP > 1 );
for ii = 1:numNeurons
    for yy = 1:length(oLapRoiX)
        roiBW{ii}(oLapRoiY(yy),oLapRoiX(yy)) = 0;
    end
end

 ROImap =  sum ( cat ( 3 , roiBW{:} ) , 3 ); % first term of cat (i.e., '3') points to element-wise alignement/stacking of arrays
tElapsedROICorr = toc(tStartROICorr);
disp(['    Time elapsed for ROI mask Correction was ', num2str(tElapsedROICorr/60),' minutes']);

%% Extract fluorescence traces for each neuron
% Loop through each neuron to get somatic fluo and neuropil fluo
for nn = 1:numNeurons
    
    [r,c]=find(roiBW{nn}~=0);
    tmpPixels = NaN(length(r),size(IMG,3));
    for i = 1:length(r)
        tmpPixels(i,:) = IMG(r(i),c(i),:);
    end
    rawFluo(:,nn) = nanmean(tmpPixels,1);
    
end

traceExtract_finish = toc(traceExtract_start);
fprintf('Trace extraction took %.1f minutes\n',traceExtract_finish/60)


end
%% FUNCTION: Sliding window dF/F calculation
% Inputs: F       - raw fluorecense
%         winsize - baseline window size
%         percent - lower percent of baseline values to average
% Output: dFF     - relative change in fluorescense in percent
function dFF = slideWinSub(F,winsize,percent)
nRois = size(F,1);
nFrames = size(F,2);
dFF = zeros(size(F));
%M = zeros(1, nrois);
%SD = zeros(1, nrois);
for j = 1 : nRois
  
    for k = 1 : nFrames
            lWin = max(1, k-winsize);
            rWin = min(k+winsize, nFrames);
            percentWin = floor(percent/100*(rWin-lWin));
            tWin = sort(F(j,(lWin:rWin)),'ascend');
            F0 = mean(tWin(1:percentWin));
%             nDFF(j,k) = DFF(j,k)- twinm;   
            dFF(j,k) = 1 * (F(j,k)- F0) / F0; %SJG looks like ZB is multiplying 100, which is used in peak finding in other section of code, to correct can divide values by 100 as well to get true normalized DFF   
    end
    
    %twin = sort(nDFF(j,:));
    %pcentwin = floor(percent/100*(nframes));
    %M(j) = mean(twin(1:pcentwin));
    %SD(j) = std(twin(1:pcentwin));
    
end

end
%% FUNCTION:  ScatterCloud (SJG test, from Desu Chen) [NOT WORKING YET]
function [h,C,C_n,xEdges,yEdges] = scattercloud(x,y,n,l,n_interp,marker_size,cmap,minX,maxX,minY,maxY)
    %SCATTERCLOUD display density of scatter data
    %   SCATTERCLOUD(X,Y) creates a scatterplot of X and Y, displayed over a
    %   surface representing the smoothed density of the points.  The density is
    %   determined with a 2D histogram, using 25 equally spaced bins in both
    %   directions.
    %   SCATTERCLOUD(X,Y,N) uses N equally spaced bins.
    %   SCATTERCLOUD(X,Y,N,L) uses L as a parameter to the smoothing algorithm.
    %    Defaults to 1.  Larger values of L lead to a smoother density, but a
    %    worse fit to the original data.
    %   SCATTERCLOUD(X,Y,N,L,CLM) uses CLM as the color/linestyle/marker for
    %    the scatter plot.  Defaults to 'k+'.
    %   SCATTERCLOUD(X,Y,N,L,CLM,CMAP) uses CMAP as the figure's colormap.  The
    %    default is 'flipud(gray(256))'.
    %   H = SCATTERCLOUD(...) returns the handles for the surface and line
    %    objects created.
    %
    %   Example:
    %
    %     scattercloud(1:100 + randn(1,100), sin(1:100) + randn(1,100),...
    %                  50,.5,'rx',jet(256))
    % 
    %   References: 
    %     Eilers, Paul H. C. & Goeman, Jelle J. (2004). Enhancing scatterplots 
    %   with smoothed densities. Bioinformatics 20(5), 623-628.


    error(nargchk(2,11,nargin),'struct');

    x = x(:);
    y = y(:);

    idx_x = (x>=minX)&(x<=maxX);
    idx_y = (y>=minY)&(y<=maxY);

    x = x(idx_x&idx_y);
    y = y(idx_x&idx_y);


    if length(x) ~= length(y)
        error('SCATTERCLOUDDataVectorSizesDoNotMatch','The number of elements in x and y do not match')
    end

    if nargin < 6
        cmap = flipud(gray(256));
    end


    if nargin < 5
        clm = 'k+';
    end

    if nargin < 4
        l = 1;
    end    

    if nargin < 3
        n = 25;
    end

    % min/max of x and y
    % minX = min(x);
    % maxX = max(x);
    % minY = min(y);
    % maxY = max(y);

    % edge locations
    xEdges = linspace(minX,maxX,n);
    yEdges = linspace(minY,maxY,n);

    % shift edges
    xDiff = xEdges(2) - xEdges(1);
    yDiff = yEdges(2) - yEdges(1);
    xEdges = [-Inf, xEdges(2:end) - xDiff/2, Inf];
    yEdges = [-Inf, yEdges(2:end) - yDiff/2, Inf];

    % number of edges
    numX = numel(xEdges);
    numY = numel(yEdges);

    % hold counts
    C = zeros(numY,numX);

    % do counts
    for i = 1:numY-1
        for j = 1:numX-1
            C(i,j) = log(length(find(x >= xEdges(j) & x < xEdges(j+1) &...
                                 y >= yEdges(i) & y < yEdges(i+1)))+1);
        end
    end

    % get rid of Infs from the edges
    xEdges = [xEdges(2) - xDiff,xEdges(2:end-1), xEdges(end-1) + xDiff];
    yEdges = [yEdges(2) - yDiff,yEdges(2:end-1), yEdges(end-1) + yDiff];

    % smooth the density data, in both directions.
    C = localSmooth(localSmooth(C,l)',l)';
    % [idx_X,idx_Y] = ndgrid(1:numX-1,1:numY-1);
    [idx_Xq,idx_Yq] = meshgrid(1:1/n_interp:numX-1,1:1/n_interp:numY-1);
    Cq = interp2(C,idx_Xq,idx_Yq,'linear');

    minC = min(min(C));
    maxC = max(max(C));
    C_n = (C-minC)/(maxC-minC);

    minCq = min(min(Cq));
    maxCq = max(max(Cq));

    % cdata = zeros(length(x),1);
    scatter_color = ones(length(x),3);
    % scatter_color_idx = ones(length(x),1);

    for i = 1:1:length(x)

        cdata = Cq(floor((y(i)-minY)/yDiff*n_interp)+1,floor((x(i)-minX)/xDiff*n_interp)+1);
        scatter_color(i,:) = cmap(round((cdata-minCq)/(maxCq-minCq)*255)+1,:);
    %     scatter_color_idx(i) = round((cdata-minC)/(maxC-minC))*255+1;
    end


    % create the graphics
    ax = newplot;
    grid(ax,'on');
    colormap(ax,cmap);
    p = scatter(x,y,marker_size,scatter_color,'filled');
    % set(p,'CDataSource','Cq');
    axis(ax,'tight');
    % % holdstate = get(ax,'NextPlot');
    % % set(ax,'NextPlot','add');
    % % 
    % % 
    % % % alphamap(ax,[0:1/256:1]);
    % % s = surf(xEdges,yEdges,zeros(numY,numX),C,...
    % %          'EdgeColor','none',...
    % %          'FaceColor','interp',...
    % %          'AlphaData',C,...
    % %          'FaceAlpha','interp');
    % % colormap(ax,cmap);
    % % alphamap(ax,[0:1/256:1]);
    % % view(ax,2);
    % % 
    % % set(ax,'NextPlot',holdstate)
    % grid(ax,'off');
    % holdstate = get(ax,'NextPlot');
    % set(ax,'NextPlot','add');
    % p = plot(x,y,clm);
    % axis(ax,'tight');
    % set(ax,'NextPlot',holdstate)

    % outputs
    if nargout
    %     h = [s;p];
        h = p;
    end


        function B = localSmooth(A,L)
            r = size(A,1);
            I = eye(r);
            D1 = diff(I);
            D2 = diff(I,2);
            B = (I + L ^ 2 * D2' * D2 + 2 * L * D1' * D1) \ A;
        end 
end
%% FUNCTION: BINNED MEAN Function
function [ mean_data2 ] = binned_mean(data1,data2,bin_num,minb,maxb,p)
%BINNED_MEAN Summary of this function goes here
%   Detailed explanation goes here
mean_data2 = zeros(bin_num,6);
% bin_size = (max(data1)-min(data1))/bin_num;
bin_size = (maxb-minb)/bin_num;
    for i = 1:1:bin_num
    %     lb = min(data1)+(i-1)*bin_size;
    %     rb = min(data1)+i*bin_size;
        lb = minb+(i-1)*bin_size;
        rb = minb+i*bin_size;
        temp_idx = (data1>=lb)&(data1<rb);
        temp_data2 = data2(temp_idx,1);
        mean_data2(i,1) =  (lb+rb)/2;
        mean_data2(i,2) = mean(temp_data2);
        SEM = std(temp_data2)/sqrt(length(temp_data2));
        ts = tinv([0.025  0.975],length(temp_data2)-1); 
        mean_data2(i,3) = mean_data2(i,2)+ts(1)*SEM;
        mean_data2(i,4) = mean_data2(i,2)+ts(2)*SEM;

        p1 = prctile(temp_data2,p);
        p2 = prctile(temp_data2,100-p);
        mean_data2(i,5) = p1;
        mean_data2(i,6) = p2;
        mean_data2(i,7) = median(temp_data2);

    end
end
%% FUNCTION: Satter kernal density estimate
function h = scatter_kde(x, y, varargin)
% Scatter plot where each point is colored by the spatial density of nearby
% points. The function use the kernel smoothing function to compute the
% probability density estimate (PDE) for each point. It uses the PDE has
% color for each point.
%
% Input
%     x <Nx1 double> position of markers on X axis
%     y <Nx1 double> posiiton of markers on Y axis
%     varargin can be used to send a set of instructions to the scatter function
%           Supports the MarkerSize parameter
%           Does not support the MarkerColor parameter
%
% Output:
%     h returns handles to the scatter objects created
%
% Example
%     % Generate data
%     x = normrnd(10,1,1000,1);
%     y = x*3 + normrnd(10,1,1000,1);
%     % Plot data using probability density estimate as function
%     figure(1); 
%     scatter_kde(x, y, 'filled', 'MarkerSize', 100);
%     % Add Color bar
%     cb = colorbar();
%     cb.Label.String = 'Probability density estimate';
%
% author: Nils Haentjens
% created: Jan 15, 2018
% Use Kernel smoothing function to get the probability density estimate (c)
c = ksdensity([x,y], [x,y]);
if nargin > 2
  % Set Marker Size
  i = find(strcmp(varargin, 'MarkerSize'),1);
  if ~isempty(i); MarkerSize = varargin{i+1}; varargin(i:i+1) = [];
  else MarkerSize = []; end
  % Plot scatter plot
  h = scatter(x, y, MarkerSize, c, varargin{:});
else
  h = scatter(x, y, [], c);
end
end


