%#######################################################################
%
%                  * MUSCle THRESHold PTOA Program *
%
%          M-file which reads a T1 FFE image file, crops the image, 
%     and thresholds the image. The user then selects subcutaneous fat
%     and femur in the left and right thighs. The program finds all
%     connected muscle and subcutaneous fat tissues. The femur is
%     femur is filled and used to exclude the femur and marrow from the
%     noncontractile elements within the muscles. The program prompts
%     the user to create a polygon region of interest around the flexor
%     muscles, and then the extensor muscles. This is used to divide the 
%     muscles into two main muscle group areas. The cross-sectional areas 
%     for the muscles, subcutaneous fat and noncontractile elements are 
%     displayed in the command window and written to an MS-Excel 
%     spreadsheet. Plots are created of the raw image; threshold 
%     histogram; and left and right muscles, subcutaneous fat and 
%     noncontractile elements.  
%
%     NOTES:  1.  DICOM images are not scaled. Only raw pixel values
%             are used to segment the muscle tissue. Designed for
%             Philips 3T T1 FFE MRI images. The background is set to
%             zero by Philips.
%
%             2.  Program assumes images are gray scale images.
%
%             3.  Otsu's method is used to pick the two thresholds.
%             The program colors pixels below the lower threshold
%             (bone) red and pixels above the upper threshold (fat)
%             green. The thresholds are shown in a plot of the signal
%             intensity histogram.
%
%             4.  See Polygon_ROI_Guide.pdf for tips on creating the
%             polygon ROI. See musc_thresh_Guide.pdf for a guide to
%             using this program.
%
%             5.  Plots are written to Postscript or PDF file
%             mthresh_ptoa_*_**.ps/.pdf, where "*" is the the visit
%             name (Baseline, Y1, or Y2) and ** is the subject number.
%
%             6.  M-file function roi_mov.m must be in the current path
%             or directory.
%
%             7.  Results are written to the MS-Excel spreadsheet
%             mthresh_ptoa_*.xlsx in the PTOA_Muscle_CSA_* folder in the 
%             current path or directory, where * is the visit part of
%             the image files directory (Baseline, Y1 or Y2).
%                 If the file does not exist, this program creates the
%             file. If the file exists, the results are appended in a
%             row at the bottom of the file.
%                 The output MS-Excel spreadsheet, mthresh_ptoa_*.xlsx,
%             can NOT be open in another program (e.g. MS-Excel, text
%             editor, etc.) while using this program.
%
%       15-Oct-2021 * Mack Gardner-Morse
%       23-June-2022 * Kate French
%       07-August-2024 * Kate French
%       04-Feb-2025 * Mack Gardner-Morse
%

%#######################################################################
%
% Clear Workspace
clc;
clear;
close all;
fclose all;
%
% Get Matlab Version
%
v = version;
idot = strfind(v,'.');
idot = idot(1)-1;
v = str2double(v(1:idot));
v = v<10;
%
% Get PTOA Visit Names
%
vnams = {'Baseline';'Y1';'Y2'};        % PTOA visit names
%
% Get T1 FFE Image File Name
%
[fnam,pnam] = uigetfile( ...
   {'*.dcm*;*.dicom*','DICOM files (*.dcm*, *.dicom*)'; ...
   '*.tif*;*.png*;*.dcm*;*.dicom*;*.jpe*;*.jpg*', ...
   'Image files (*.tif*, *.png*, *.dcm*, *.dicom*, *.jpe*, *.jpg*)'; ...
   '*.*','All files (*.*)'},'Please Select Image File', ...
   'MultiSelect', 'off');
ffnam = fullfile(pnam,fnam);
%
% Get Results Directory
%
pnams = strsplit(pnam);
vnam = pnams{2};
%
if ~any((startsWith(vnams,vnam)))
  error([' *** ERROR in musc_thresh_ptoa:  Image directory does', ...
         ' not have the correct PTOA visit name!']);
end
% 
rdir = ['PTOA_Muscle_CSA_' vnam];      % Results folder
if ~exist(rdir,'dir')   % Create results directory if doesn't exist
  mkdir(rdir);
end 
%
% Output MS_Excel Spreadsheet File Name, Sheet Name and Headers/Units
%
xlsnam = fullfile(rdir,['mthresh_ptoa_' vnam '.xlsx']);    % Put in results directory
shtnam = 'PTOA Study';       % Sheet name
hdr = {'Subj ID','Subj #','MRI Date','Analysis Date', ...
       'L Mus CSA Ext','R Mus CSA Ext','L Mus CSA Flex', ...
       'R Mus CSA Flex','L Mus CSA total','R Mus CSA total', ...
       'L SubFat CSA','R SubFat CSA','L Non Con CSA', ...
       'R Non Con CSA'};     % Column headers
ulbls = [{'','','',''} cellstr(repmat('(cm^2)',10,1))'];    % Units
%
% Postscript Output File for Plots
%
pfile = fullfile(rdir,['mthresh_ptoa_' vnam]);   % Postscript output file
%
idot = strfind(fnam,'.');
idot = idot(end);            % Get last dot
ImageName = fnam(1:idot);    % Remove file extension
if v
  pfile = [pfile '_' ImageName 'ps'];  % Postscript output file
else
  pfile = [pfile '_' ImageName 'pdf']; % PDF output file
end
ImageName = ImageName(1:idot-1);       % Remove dot
%
% Read T1 FFE Image from File and Get Range
%
id = isdicom(ffnam);
if id                   % DICOM image format
    info = dicominfo(ffnam);
    img = dicomread(info);
else                    % Other image formats
    img = imread(ffnam);
end
%
image_sz = size(img);
if length(image_sz)>2
  img = img(:,:,1);     % Assume all channels are equal (gray scale)
end
%
rmin = min(img(:));
rmax = max(img(:));
%
% Setup T1 FFE Figure Window
%
hf1 = figure;
orient landscape;
set(hf1,'WindowState','maximized');
pause(0.1);
drawnow;
%
% Set Up Color Map for Thresholding
%
cmap = gray(128);
cmap(1,:) = [1 0 0];            % Red for lower threshold
cmap(128,:) = [0 0.7 0];        % Green for upper threshold
%
% Display T1 FFE Image
%
imagesc(img,[rmin rmax]);
colormap gray;
axis image;
axis off;
title({ImageName; 'Original T1 FFE Image'},'FontSize',16, ...
      'FontWeight','bold','Interpreter','none');
if v
  print('-dpsc2','-r600','-fillpage',pfile);
else
  exportgraphics(gcf,pfile,'Resolution',600,'Padding','figure');
end
%
% Crop Background From Image
%
[i,j] = find(img>rmin+1);   % Background
idi = min(i):max(i);
idj = min(j):max(j);
imgc = img(idi,idj);        % Cropped image

% Plot Cropped Image
hf2 = figure;
orient landscape;
set(hf2,'WindowState','maximized');
pause(0.1);
drawnow;
%
imagesc(imgc,[rmin rmax]);
colormap gray;
axis image;
axis off;
title({ImageName; 'Cropped T1 FFE Image'},'FontSize',16, ...
      'FontWeight','bold','Interpreter','none');
if v
  print('-dpsc2','-r600','-fillpage','-append',pfile);
else
  exportgraphics(gcf,pfile,'Resolution',600,'Padding','figure', ...
                 'Append',true);
end
%
% Get Thresholds and Plot Image Histogram
%
lvls = multithresh(imgc,2);     % Otsu's method
%
hf3 = figure;
orient landscape;
set(hf3,'WindowState','maximized');
pause(0.1);
drawnow;
%
histogram(imgc,'FaceAlpha',1,'FaceColor',[0 0 0.8]);
hold on;
axlim = axis;
plot(repmat(lvls(1),2,1),axlim(3:4)','r-','LineWidth',1.0);
plot(repmat(lvls(2),2,1),axlim(3:4)','g-','Color',[0 0.7 0], ...
      'LineWidth',1.0);
%
yt = get(gca,'YTick');
yt = yt(end-1);
text(double(lvls(1))+2,yt,int2str(lvls(1)),'FontSize',12,'Color','r');
text(double(lvls(2))+2,yt,int2str(lvls(2)),'FontSize',12, ...
     'Color',[0 0.7 0]);
%
xlabel('Signal Intensity','FontSize',12,'FontWeight','bold');
ylabel('Frequency','FontSize',12,'FontWeight','bold');
title('T1 FFE Image Histogram','FontSize',16,'FontWeight','bold');
if v
  print('-dpsc2','-r600','-fillpage','-append',pfile);
else
  exportgraphics(gcf,pfile,'Resolution',600,'Padding','figure', ...
                 'Append',true);
end
%
% Apply Thresholds and Plot
%
imgt = imgc;                % Image with thresholding
imgt(imgt<lvls(1)) = rmin;
imgt(imgt>lvls(2)) = rmax;
%
hf4 = figure;
orient landscape;
set(hf4,'WindowState','maximized');
pause(0.1);
drawnow;
%
imagesc(imgt,[rmin rmax]);
colormap(cmap);
axis image;
axis off;
title({ImageName; 'T1 FFE Image with Thresholds'},'FontSize',16, ...
      'FontWeight','bold','Interpreter','none');
if v
  print('-dpsc2','-r600','-fillpage','-append',pfile); 
else
  exportgraphics(gcf,pfile,'Resolution',600,'Padding','figure', ...
                 'Append',true);
end
%
% Identify Left Thigh Subcutaneous Fat
%
uiwait(msgbox(['Pick a point within the left thigh subcutaneous ', ...
       'fat.'],'Input','modal'));
[xl,yl] = ginput(1);
xl = round(xl);
yl = round(yl);
%
bwl = grayconnected(imgt,yl,xl,1);     % Get subcutaneous fat
imgl = imgc;
imgl(~bwl) = 0;
%
[i,j] = find(imgl>rmin+1);             % Crop background
idil = min(i):max(i);
idjl = min(j):max(j);
imgl = imgl(idil,idjl);                % Cropped image
%
rminl = min(imgl(:));
rmaxl = max(imgl(:));
%
% Identify Left Thigh Femur
%
uiwait(msgbox('Pick a point within the left thigh femur.', ...
       'Input','modal'));
[xl,yl] = ginput(1);
xl = round(xl);
yl = round(yl);
%
bwfl = grayconnected(imgt,yl,xl,1);    % Get femur
bwfl = imfill(bwfl,'holes');           % Femur and marrow
%
% Identify Left Thigh Muscle
%
% uiwait(msgbox('Pick a point within the left "bulk" muscle.', ... 
%     'Input','modal'));
% [xl,yl] = ginput(1);
% xl = round(xl);
% yl = round(yl);
%
imgbwl = imgt;
imgbwl(imgbwl==rmax) = 0;
imgbwl(imgbwl>=lvls(1)&imgbwl<=lvls(2)) = rmax;
imgbwl(bwfl) = 0;       % Remove femur and marrow
imgbwl = imgbwl(idil,idjl);            % Cropped image
%
cc = bwconncomp(imgbwl);
pidx = cc.PixelIdxList';
n = cell2mat(cellfun(@length,pidx,'UniformOutput',false));
nidx = n>200;           % Typical small muscles > 2000 pixels (10x)
pidx = pidx(nidx);
bwml = false(size(imgbwl));
bwml(cell2mat(pidx)) = true;
%
imgml = imgt;
imgml = imgml(idil,idjl);               % Cropped image
imgml(~bwml) = 0;
%
% Identify Right Thigh Subcutaneous Fat
%
uiwait(msgbox(['Pick a point within the right thigh subcutaneous ', ...
       'fat.'],'Input','modal'));
[xr,yr] = ginput(1);
xr = round(xr);
yr = round(yr);
%
bwr = grayconnected(imgt,yr,xr,1);     % Get subcutaneous fat
imgr = imgc;
imgr(~bwr) = 0;
%
[i,j] = find(imgr>rmin+1);             % Crop background
idir = min(i):max(i);
idjr = min(j):max(j);
imgr = imgr(idir,idjr);                % Cropped image
%
rminr = min(imgr(:));
rmaxr = max(imgr(:));
%
% Identify Right Thigh Femur
%
uiwait(msgbox('Pick a point within the right thigh femur.', ...
       'Input','modal'));
[xr,yr] = ginput(1);
xr = round(xr);
yr = round(yr);
%
bwfr = grayconnected(imgt,yr,xr,1);    % Get femur
bwfr = imfill(bwfr,'holes');           % Femur and marrow
%
% Identify Right Thigh Muscle
%
% uiwait(msgbox('Pick a point within the right "bulk" muscle.', ...
%     'Input','modal'));
% [xr,yr] = ginput(1);
% xr = round(xr);
% yr = round(yr);
%   
imgbwr = imgt;
imgbwr(imgbwr==rmax) = 0;
imgbwr(bwfr) = 0;       % Remove femur and marrow
imgbwr(imgbwr>=lvls(1)&imgbwr<=lvls(2)) = rmax;
imgbwr = imgbwr(idir,idjr);            % Cropped image
%
cc = bwconncomp(imgbwr);
pidx = cc.PixelIdxList';
n = cell2mat(cellfun(@length,pidx,'UniformOutput',false));
nidx = n>200;           % Typical small muscles > 2000 pixels (10x)
pidx = pidx(nidx);
bwmr = false(size(imgbwr));
bwmr(cell2mat(pidx)) = true;
%
imgmr = imgt;
imgmr = imgmr(idir,idjr);               % Cropped image
imgmr(~bwmr) = 0;
%
% Get ROI for Left Flexor Muscles
% Includes the hamstrings (biceps femoris short & long head, 
% semitendinosus, semimembranosus), hip adductors (adductor longus, 
% adductor magnus), gracilis and sartorius
%
imgfl = imgc(idil,idjl);
hf5 = figure;
orient landscape;
set(hf5,'WindowState','maximized');
pause(0.1);
drawnow;
%
imagesc(imgfl);
colormap gray;
axis image;
ha5 = gca;
axis off;
title({ImageName; 'Left Thigh Muscle T1 FFE Image'},'FontSize',16, ...
      'FontWeight','bold','Interpreter','none');
%
uiwait(msgbox({'Please digitize the left flexor muscles.'; ' '; ...
       'Press <Enter> when finished.'},'Input','modal'));
%
hlfroi = images.roi.Polygon(ha5,'LineWidth',1);
addlistener(hlfroi,'MovingROI',@roi_mov);
hlfroi.draw; 
%
kchk = true;
%
while kchk
    hlfroi.wait;
    lfm = hlfroi.createMask;           % Left flexor mask
    Fpts = hlfroi.Position;            % Flexor points
    Fpts = [Fpts; Fpts(1,:)];
%
% Get Flexors
%
    imgmlf = imgml;
    imgmlf(~lfm) = 0;
%
% Set Up Figure and Figure Menus
%
    hlf = figure;
    orient landscape;
    set(hlf,'WindowState','maximized');
    pause(0.1);
    drawnow;
%
    hw = imagesc(imgml);
    colormap gray;
    axis image;
    axis off;
    title({ImageName; 'Left Flexor Muscle T1 FFE Image'}, ...
           'FontSize',16,'FontWeight','bold','Interpreter','none');
    hold on;
%
    imgmle_t = imgml;                  % Temporary extensor muscle mask
    imgmle_t(lfm) = 0;
    he = imagesc(imgmle_t);
    set(he,'Visible','off');
    hf = imagesc(imgmlf);
    set(hf,'Visible','off');
    hm = plot(Fpts(:,1),Fpts(:,2),'r-','LineWidth',1);
%
    kchk = logical(menu('Flexor Muscles OK?','Yes','No')-1);
%
    close(hlf);
%
    if kchk
        figure(hf5);
    end
end
%
close(hf5);
%
% Get ROI for Left Extensor Muscles
% Include the quadriceps (rectus femoris, vastus medialis, vastus 
% intermedius, vastus lateralis)
%
imgel = imgc(idil,idjl);
hf6 = figure;
orient landscape;
set(hf6,'WindowState','maximized');
pause(0.1);
drawnow;
%
imagesc(imgel);
colormap gray;
axis image;
ha6 = gca;
axis off;
title({ImageName; 'Left Thigh Muscle T1 FFE Image'},'FontSize',16, ...
      'FontWeight','bold','Interpreter','none');
%
uiwait(msgbox({'Please digitize the left extensor muscles.'; ' '; ...
       'Press <Enter> when finished.'},'Input','modal'));
%
hleroi = images.roi.Polygon(ha6,'LineWidth',1);
addlistener(hleroi,'MovingROI',@roi_mov);
hleroi.draw;
%
kchk = true;
while kchk 
    hleroi.wait;
    lem = hleroi.createMask;           % Left extensor mask
    Epts = hleroi.Position;
    Epts = [Epts; Epts(1,:)];
%
% Get Extensors
%
    imgmle = imgml;
    imgmle(~lem) = 0;
%
% Set Up Figure and Figure Menus
%
    hle = figure;
    orient landscape;
    set(hle,'WindowState','maximized');
    pause(0.1);
    drawnow;
%
    hw = imagesc(imgml);
    colormap gray;
    axis image;
    axis off;
    title({ImageName; 'Left Extensor Muscle T1 FFE Image'}, ...
           'FontSize',16,'FontWeight','bold','Interpreter','none');
    hold on;
%
    he = imagesc(imgmle);
    set(he,'Visible','off');
    hf = imagesc(imgmlf);
    set(hf,'Visible','off');
    hm = plot(Epts(:,1),Epts(:,2),'r-','LineWidth',1);
%
    kchk = logical(menu('Extensor Muscles OK?','Yes','No')-1);
%
    close(hle);
%
    if kchk
        figure(hf6);
    end
end
%
close(hf6);
%
% Get ROI for Right Flexor Muscles
% Includes the hamstrings (biceps femoris short & long head, 
% semitendinosus, semimembranosus), hip adductors (adductor longus, 
% adductor magnus), gracilis and sartorius
%
imgfr = imgc(idir,idjr);
hf7 = figure;
orient landscape;
set(hf7,'WindowState','maximized');
pause(0.1);
drawnow;
%
imagesc(imgfr);
colormap gray;
axis image;
ha7 = gca;
axis off;
title({ImageName; 'Right Thigh Muscle T1 FFE Image'},'FontSize',16, ...
      'FontWeight','bold','Interpreter','none');
%
uiwait(msgbox({'Please digitize the right flexor muscles.'; ' '; ...
       'Press <Enter> when finished.'},'Input','modal'));
%
hrfroi = images.roi.Polygon(ha7,'LineWidth',1);
addlistener(hrfroi,'MovingROI',@roi_mov);
hrfroi.draw;
%
kchk = true;
%
while kchk
    hrfroi.wait;
    rfm = hrfroi.createMask;           % Right flexor mask
    Fpts = hrfroi.Position;            % Flexor points
    Fpts = [Fpts; Fpts(1,:)];
%
% Get Flexors
%
    imgmrf = imgmr;
    imgmrf(~rfm) = 0;
%
% Set Up Figure and Figure Menus
%
    hrf = figure;
    orient landscape;
    set(hrf,'WindowState','maximized');
    pause(0.1);
    drawnow;
%
    hw = imagesc(imgmr);
    colormap gray;
    axis image;
    axis off;
    title({ImageName; 'Right Flexor Muscle T1 FFE Image'}, ...
           'FontSize',16,'FontWeight','bold','Interpreter','none');
    hold on;
%
    imgmre_t = imgmr;                  % Temporary extensor muscle mask
    imgmre_t(rfm) = 0;
    he = imagesc(imgmre_t);
    set(he,'Visible','off');
    hf = imagesc(imgmrf);
    set(hf,'Visible','off');
    hm = plot(Fpts(:,1),Fpts(:,2),'r-','LineWidth',1);
%
    kchk = logical(menu('Flexor Muscles OK?','Yes','No')-1);
%
    close(hrf);
%
    if kchk
        figure(hf7);
    end
end
%
close(hf7);
%
% Get ROI for Right Extensor Muscles
% Include the quadriceps (rectus femoris, vastus medialis, vastus 
% intermedius, vastus lateralis)
%
imger = imgc(idir,idjr);
hf8 = figure;
orient landscape;
set(hf8,'WindowState','maximized');
pause(0.1);
drawnow;
%
imagesc(imger);
colormap gray;
axis image;
ha8 = gca;
axis off;
title({ImageName; 'Right Thigh Muscle T1 FFE Image'},'FontSize',16, ...
      'FontWeight','bold','Interpreter','none');
%
uiwait(msgbox({'Please digitize the right extensor muscles.'; ' '; ...
       'Press <Enter> when finished.'},'Input','modal'));
%
hreroi = images.roi.Polygon(ha8,'LineWidth',1);
addlistener(hreroi,'MovingROI',@roi_mov);
hreroi.draw;
%
kchk = true;
while kchk 
    hreroi.wait;
    rem = hreroi.createMask;           % Right extensor mask
    Epts = hreroi.Position;
    Epts = [Epts; Epts(1,:)];
%
% Get Extensors
%
    imgmre = imgmr;
    imgmre(~rem) = 0;
%
% Set Up Figure and Figure Menus
%
    hre = figure;
    orient landscape;
    set(hre,'WindowState','maximized');
    pause(0.1);
    drawnow;
%
    hw = imagesc(imgmr);
    colormap gray;
    axis image;
    axis off;
    title({ImageName; 'Right Extensor Muscle T1 FFE Image'},'FontSize',16, ...
           'FontWeight','bold','Interpreter','none');
    hold on;
%
    he = imagesc(imgmre);
    set(he,'Visible','off');
    hf = imagesc(imgmrf);
    set(hf,'Visible','off');
    hm = plot(Epts(:,1),Epts(:,2),'r-','LineWidth',1);
%
    kchk = logical(menu('Extensor Muscles OK?','Yes','No')-1);
%
    close(hre);
%
    if kchk
        figure(hf8);
    end
end
%
close(hf8);
%
% Get Left Thigh Noncontractile Elements
%
bwlf = imfill(bwml,'holes');           % Get whole muscle area
bwlf = bwlf&~bwfl(idil,idjl);          % Remove femur and marrow
imgncl = imgc(idil,idjl);
imgncl(~bwlf) = 0;
bwncl = imgncl>lvls(2);
imgncl(~bwncl) = 0;
%
% Get Right Thigh Noncontractile Elements
%
bwrf = imfill(bwmr,'holes');           % Get whole muscle area
bwrf = bwrf&~bwfr(idir,idjr);          % Remove femur and marrow
imgncr = imgc(idir,idjr);
imgncr(~bwrf) = 0;
bwncr = imgncr>lvls(2);
imgncr(~bwncr) = 0;
%
% Get Pixel Size
%
if isfield(info,'PixelSpacing')
    pix2cm = prod(info.PixelSpacing)/100;   % Conversion from pixel to cm
    units = 'cm^2';
else
    pix2cm = 1;
    units = 'pixel^2';
    warning([' *** Warning in musc_thresh:  Pixel size not found! ', ...
           ' Cross-sectional areas will be in pixel^2 and not cm^2!']);
end
%
% Calculate and Print Out Muscle Cross-Sectional Areas
%
bwmle = bwml;
bwmle(~lem) = 0;
areamle = sum(bwmle(:))*pix2cm;        % Left extensors
bwmlf = bwml;
bwmlf(~lfm) = 0;
areamlf = sum(bwmlf(:))*pix2cm;        % Left flexors
areaml = areamle+areamlf;              % Total left muscle
%
bwmre = bwmr;
bwmre(~rem) = 0;
areamre = sum(bwmre(:))*pix2cm;        % Right extensors
bwmrf = bwmr;
bwmrf(~rfm) = 0;
areamrf = sum(bwmrf(:))*pix2cm;        % Right flexors
areamr = areamre+areamrf;              % Total right muscle
%
fprintf(1,'\n\nMUSCLE CROSS-SECTIONAL AREAS FOR %s\n',ImageName);
fprintf(1,'Left Thigh Cross-Sectional Area = %.1f %s\n',areaml,units);
fprintf(1,'  Extensors Cross-Sectional Area = %.1f %s\n',areamle, ...
        units);
fprintf(1,'  Flexors Cross-Sectional Area = %.1f %s\n',areamlf,units);
%
fprintf(1,'Right Thigh Cross-Sectional Area = %.1f %s\n',areamr, ...
        units);
fprintf(1,'  Extensors Cross-Sectional Area = %.1f %s\n',areamre, ...
        units);
fprintf(1,'  Flexors Cross-Sectional Area = %.1f %s\n',areamrf,units);

% Calculate and Print Out Subcutaneous Fat Cross-Sectional Areas
areal = sum(bwl(:))*pix2cm;
arear = sum(bwr(:))*pix2cm;
%
fprintf(1,'\nSUBCUTANEOUS FAT CROSS-SECTIONAL AREAS FOR %s\n',ImageName);
fprintf(1,'Left Thigh Cross-Sectional Area = %.1f %s\n',areal,units);
fprintf(1,'Right Thigh Cross-Sectional Area = %.1f %s\n',arear,units);

% Calculate and Print Out Noncontractile Elements Cross-Sectional Areas
areancl = sum(bwncl(:))*pix2cm;
areancr = sum(bwncr(:))*pix2cm;
%
fprintf(1,['\nNONCONTRACTILE ELEMENTS CROSS-SECTIONAL AREAS ', ...
           'FOR %s\n'],ImageName);
fprintf(1,'Left Thigh Cross-Sectional Area = %.1f %s\n',areancl,units);
fprintf(1,'Right Thigh Cross-Sectional Area = %.1f %s\n\n\n', ...
        areancr,units);
%
% Check for Output MS-Excel Spreadsheet and Check for Headers/Units
%
if ~exist(xlsnam,'file')
  writecell(hdr,xlsnam,'Sheet',shtnam,'Range','A1');
  writecell(ulbls,xlsnam,'Sheet',shtnam,'Range','A2');
else
  [~,fshtnams] = xlsfinfo(xlsnam);     % Get sheet names in file
  idl = strcmp(shtnam,fshtnams);       % Sheet already exists in file?
  if all(~idl)        % Sheet name not found in file
    writecell(hdr,xlsnam,'Sheet',shtnam,'Range','A1');
    writecell(ulbls,xlsnam,'Sheet',shtnam,'Range','A2');
  else                  % Sheet in the file
    [~,txt] = xlsread(xlsnam,shtnam);
    if size(txt,1)<2    % No headers
      writecell(hdr,xlsnam,'Sheet',shtnam,'Range','A1');
      writecell(ulbls,xlsnam,'Sheet',shtnam,'Range','A2');
    end
  end
end
%
% Get MRI and Today's Dates
%
if isfield(info,'SeriesDate')
  sdate = info.SeriesDate;             % Date of MRI
  sdate = datenum(sdate,'yyyymmdd');
  sdate = datestr(sdate,1);            % Formatted date of MRI
else
  sdate = '';
end
%
adate = date;           % Analysis date (today)
%
% Put Results into a Table and Write to Output Spreadsheet File
%
sID = string(ImageName);
sName = string(fnam(1:3));
t = table(sID,sName,{sdate},{adate},areamle,areamre,areamlf, ...
          areamrf,areaml,areamr,areal,arear,areancl,areancr);
writetable(t,xlsnam,'WriteMode','append','Sheet',shtnam, ...
           'WriteVariableNames',false);
%
% Plot Final Left Thigh Muscle and Left Extensors and Flexors
%
hf5 = figure;
orient landscape;
set(hf5,'WindowState','maximized','InvertHardcopy','off');
pause(0.1);
drawnow;
%
imgmlef = imgmle;
idf = imgmle<1;
imgmlef(idf) = imgmlf(idf);
%
imagesc(imgmlef);
colormap gray;
axis image;
axis off;
hold on;
%
text(2,2,sprintf('Cross-sectional area = %.1f %s\n',areaml, ...
     units),'HorizontalAlign','left','VerticalAlignment','top', ...
     'Color','w','FontSize',11,'FontWeight','bold','Color','w');
%
title({ImageName; 'Left Thigh Muscle T1 FFE Image'},'FontSize',16, ...
      'FontWeight','bold','Interpreter','none');
if v
  print('-dpsc2','-r600','-fillpage','-append',pfile);
else
  exportgraphics(gcf,pfile,'Resolution',600,'Padding','figure', ...
                 'Append',true);
end
%
hf6 = figure;
orient tall;
set(hf6,'WindowState','maximized','InvertHardcopy','off');
pause(0.1);
drawnow;
%
subplot(2,1,1);
imagesc(imgmle);
colormap gray;
axis image;
axis off;
hold on;
%
dims = size(imgmle);
xt = dims(2)/2;
yt = 0.85*dims(1);
text(xt,yt,sprintf('Cross-sectional area = %.1f %s\n',areamle, ...
     units),'HorizontalAlign','center','VerticalAlignment','top', ...
     'Color','w','FontSize',11,'FontWeight','bold','Color','w');
%
title({'Left Thigh Muscle T1 FFE Image'; 'Extensor Muscles'}, ...
      'FontSize',16,'FontWeight','bold','Interpreter','none');
%
subplot(2,1,2);
imagesc(imgmlf);
colormap gray;
axis image;
axis off;
hold on;
%
dims = size(imgmlf);
xt = dims(2)/2;
yt = 0.02*dims(1);
text(xt,yt,sprintf('Cross-sectional area = %.1f %s\n',areamlf, ...
     units),'HorizontalAlign','center','VerticalAlignment','top', ...
     'Color','w','FontSize',11,'FontWeight','bold','Color','w');
%
title('Flexor Muscles','FontSize',16,'FontWeight','bold');
if v
  print('-dpsc2','-r600','-fillpage','-append',pfile);
else
  exportgraphics(gcf,pfile,'Resolution',600,'Padding','figure', ...
                 'Append',true);
end
%
% Plot Final Right Thigh Muscle and Right Extensors and Flexors
%
hf7 = figure;
orient landscape;
set(hf7,'WindowState','maximized','InvertHardcopy','off');
pause(0.1);
drawnow;
%
%
imgmref = imgmre;
idf = imgmre<1;
imgmref(idf) = imgmrf(idf);
%
imagesc(imgmref);
colormap gray;
axis image;
axis off;
hold on;
%
dims = size(imgmr,2);
text(dims-1,2,sprintf('Cross-sectional area = %.1f %s\n',areamr, ...
     units),'HorizontalAlign','right','VerticalAlignment','top', ...
     'Color','w','FontSize',11,'FontWeight','bold','Color','w');
%
title({ImageName; 'Right Thigh Muscle T1 FFE Image'},'FontSize',16, ...
      'FontWeight','bold','Interpreter','none');
if v
  print('-dpsc2','-r600','-fillpage','-append',pfile);
else
  exportgraphics(gcf,pfile,'Resolution',600,'Padding','figure', ...
                 'Append',true);
end
%
hf8 = figure;
orient tall;
set(hf8,'WindowState','maximized','InvertHardcopy','off');
pause(0.1);
drawnow;
%
subplot(2,1,1);
imagesc(imgmre);
colormap gray;
axis image;
axis off;
hold on;
%
dims = size(imgmre);
xt = dims(2)/2;
yt = 0.85*dims(1);
text(xt,yt,sprintf('Cross-sectional area = %.1f %s\n',areamre, ...
     units),'HorizontalAlign','center','VerticalAlignment','top', ...
     'Color','w','FontSize',11,'FontWeight','bold','Color','w');
%
title({'Right Thigh Muscle T1 FFE Image'; 'Extensor Muscles'}, ...
      'FontSize',16,'FontWeight','bold','Interpreter','none');
%
subplot(2,1,2);
imagesc(imgmrf);
colormap gray;
axis image;
axis off;
hold on;
%
dims = size(imgmrf);
xt = dims(2)/2;
yt = 0.02*dims(1);
text(xt,yt,sprintf('Cross-sectional area = %.1f %s\n',areamrf, ...
     units),'HorizontalAlign','center','VerticalAlignment','top', ...
     'Color','w','FontSize',11,'FontWeight','bold','Color','w');
%
title('Flexor Muscles','FontSize',16,'FontWeight','bold');
if v
  print('-dpsc2','-r600','-fillpage','-append',pfile);
else
  exportgraphics(gcf,pfile,'Resolution',600,'Padding','figure', ...
                 'Append',true);
end
%
% Display Left Thigh Subcutaneous Fat T1 FFE Image
%
hf9 = figure;
orient landscape;
set(hf9,'WindowState','maximized','InvertHardcopy','off');
pause(0.1);
drawnow;
%
imagesc(imgl,[rminl rmaxl]);
colormap gray;
axis image;
axis off;
hold on;
%
dims = size(imgl,2);
text(dims-1,2,sprintf('Cross-sectional area = %.1f %s\n',areal, ...
     units),'HorizontalAlign','right','VerticalAlignment','top', ...
     'Color','w','FontSize',11,'FontWeight','bold','Color','w');
%
title({ImageName; 'Left Thigh Subcutaneous Fat T1 FFE Image'}, ...
      'FontSize',16,'FontWeight','bold','Interpreter','none');
if v
  print('-dpsc2','-r600','-fillpage','-append',pfile);
else
  exportgraphics(gcf,pfile,'Resolution',600,'Padding','figure', ...
                 'Append',true);
end
%
% Display Right Thigh Subcutaneous Fat T1 FFE Image
%
hf10 = figure;
orient landscape;
set(hf10,'WindowState','maximized','InvertHardcopy','off');
pause(0.1);
drawnow;
%
imagesc(imgr,[rminr rmaxr]);
colormap gray;
axis image;
axis off;
hold on;
%
text(2,2,sprintf('Cross-sectional area = %.1f %s\n',arear, ...
     units),'HorizontalAlign','left','VerticalAlignment','top', ...
     'Color','w','FontSize',11,'FontWeight','bold','Color','w');
%
title({ImageName; 'Right Thigh Subcutaneous Fat T1 FFE Image'}, ...
      'FontSize',16,'FontWeight','bold','Interpreter','none');
if v
  print('-dpsc2','-r600','-fillpage','-append',pfile);
else
  exportgraphics(gcf,pfile,'Resolution',600,'Padding','figure', ...
                 'Append',true);
end
%
% Display Left Thigh Noncontractile Elements T1 FFE Image
%
hf11 = figure;
orient landscape;
set(hf11,'WindowState','maximized','InvertHardcopy','off');
pause(0.1);
drawnow;
%
imagesc(imgncl,[rminl rmaxl]);
colormap gray;
axis image;
axis off;
hold on;
%
text(2,2,sprintf('Cross-sectional area = %.1f %s\n',areancl, ...
     units),'HorizontalAlign','left','VerticalAlignment','top', ...
     'Color','w','FontSize',11,'FontWeight','bold','Color','w');
%
title({ImageName; 'Left Thigh Noncontractile Elements T1 FFE Image'}, ...
      'FontSize',16,'FontWeight','bold','Interpreter','none');
if v
  print('-dpsc2','-r600','-fillpage','-append',pfile);
else
  exportgraphics(gcf,pfile,'Resolution',600,'Padding','figure', ...
                 'Append',true);
end
%
% Display Right Thigh Noncontractile Elements T1 FFE Image
%
hf12 = figure;
orient landscape;
set(hf12,'WindowState','maximized','InvertHardcopy','off');
pause(0.1);
drawnow;
%
imagesc(imgncr,[rminr rmaxr]);
colormap gray;
axis image;
axis off;
hold on;
%
dims = size(imgncr,2);
text(dims-1,2,sprintf('Cross-sectional area = %.1f %s\n',areancr, ...
     units),'HorizontalAlign','right','VerticalAlignment','top', ...
     'Color','w','FontSize',11,'FontWeight','bold','Color','w');
%
title({ImageName; 'Right Thigh Noncontractile Elements T1 FFE Image'}, ...
      'FontSize',16,'FontWeight','bold','Interpreter','none');
if v
  print('-dpsc2','-r600','-fillpage','-append',pfile);
else
  exportgraphics(gcf,pfile,'Resolution',600,'Padding','figure', ...
                 'Append',true);
end
%
return