function varargout = cellFluorescence(varargin)

%--------------------------------------------------------------------------
%
% cellFluorescence:
%   The GUI calculates the total corrected fluorescence.
%
% Description:
%   The total corrected fluorescence intensity per cell (TCF) can be 
%   determined from fluorescence images. Using the GUI cellFluorescence, the 
%   boundary region of a single cell needs to be outlined with a freehand 
%   area selection tool. The boundary region is then used to calculate the 
%   cell area and the integrated density as the sum of all pixel intensities
%   within the boundary region. Then any region near the cell boundary, but 
%   without cellular fluorescence, needs to be selected to measure the 
%   average background fluorescence. The TCF is then calculated using the
%   formula:
%   TCF = Integrated Density – (cell Area * Mean Background Fluorescence Intensity) 
%   The TCF value can be used to quantify dye labeling.    
% 
% Syntax:  
%   cellFluorescence
% 
% Inputs:
%   The file explorer prompts the user to select Tiff image files. 
%   Multiple image files can be selected. 
% 
% Outputs:
%   TCF Value 
% 
% Author: 
%   P.G. Nov 2016
%
% Copyright:
%   2011-2020 Javitch Lab - Columbia Univiversity / RFMH
% 
%--------------------------------------------------------------------------

% CELLFLUORESCENCE MATLAB code for cellFluorescence.fig
%      CELLFLUORESCENCE, by itself, creates a new CELLFLUORESCENCE or raises the existing
%      singleton*.
%
%      H = CELLFLUORESCENCE returns the handle to a new CELLFLUORESCENCE or the handle to
%      the existing singleton*.
%
%      CELLFLUORESCENCE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLFLUORESCENCE.M with the given input arguments.
%
%      CELLFLUORESCENCE('Property','Value',...) creates a new CELLFLUORESCENCE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cellFluorescence_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cellFluorescence_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help cellFluorescence
% Last Modified by GUIDE v2.5 28-Nov-2016 12:13:12
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cellFluorescence_OpeningFcn, ...
                   'gui_OutputFcn',  @cellFluorescence_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
% -------------------------------------------------------------------------
%                           OPENING FUNCTIONS
% -------------------------------------------------------------------------
% --- Executes just before cellFluorescence is made visible.
function cellFluorescence_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cellFluorescence (see VARARGIN)
% Choose default command line output for cellFluorescence
handles.output = hObject;
cla(handles.axesIm,'reset');
cla(handles.axesInt,'reset');
cla(handles.axesBgd,'reset');
% Load Constants
constants        = smCellConstants('');
handles.photonConversion = constants.photonConvFactorDelta;
handles.QE               = 0.95; % set Camera Quantum efficiency if not specified
% Update handles structure
guidata(hObject, handles);
% --- Outputs from this function are returned to the command line.
function varargout = cellFluorescence_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;
% -------------------------------------------------------------------------
%                              OPEN FILE
% -------------------------------------------------------------------------
function openFile_Callback(hObject, eventdata, handles)
% hObject    handle to openFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Clear Axes Objects
cla(handles.axesIm,'reset')
cla(handles.axesInt,'reset')
cla(handles.axesBgd,'reset')
output = [ 0 0 0 0 0 ]; 
set( handles.editOutput,'String',num2str(output,...
       '%d   %d   %d   %.1f   %.1f'));
% Initialize intensity slider for Ch2
newIntRange = [];
sldRange    = [1 10];
stepSize1   = 0.05;
stepSize2   = 2*stepSize1;
set(handles.sliderIntensity,'Min',sldRange(1),'Max',sldRange(2),...
    'Value',max(sldRange),'SliderStep',[stepSize1,stepSize2]);
% Initialize toleranz value in 'stretchlim.m' 
handles.Tol = 0.0005;
set( handles.editToleranz,'String',num2str(handles.Tol));
% Initialize gui buttons' 
set(handles.btnPrevious,'Enable','on');
set(handles.btnNext,'Enable','on');
% Get traces filename by menu driven input
files=getFiles('*.tif', 'Open Files');
if isempty(files)
    return
else
    [handles.basepath,name,ext]=fileparts(files{1});
    handles.filename=[name ext];
    handles.nFiles=numel(files);
end
if handles.nFiles==0
    disp('NO FILES AVAILABLE IN FOLDER');
    return;
else
    cd(handles.basepath);
end
% Load first image
handles.fileNo = 1;
imInfo = imfinfo(files{1});
% Update directory and dataset
set(handles.editDir,'String',handles.basepath);
set(handles.editFilename,'String',handles.filename);
set(handles.editId,'String',num2str(handles.fileNo));
% Initialize coordinates for rectROI
rectRoiCoord = [10 10 50 50];
set(handles.editRectRoiStartPos,'String',num2str(rectRoiCoord));
% Initialize settings for local fluorescence
handles.widthPsf = 7'; % 7x7 pixel array
handles.stepPsf  = 4;  % for 7x7 array
handles.nSelMol  = 1; 
set(handles.editNCoordinates,'String',num2str(handles.nSelMol));
% Load image
currentImage=imread(files{1});
% Display the original image and resized image 
axes(handles.axesIm);
intRange     = stretchlim(currentImage,handles.Tol);
currentImage = imadjust(currentImage,[intRange(1) intRange(2)]);
imagesc(currentImage);colormap('gray');
axis([0 imInfo.Width 0 imInfo.Height]);
xlabel('x [Pixels]');
ylabel('y [Pixels]');
title(['Cell ' num2str(handles.fileNo) ' of ' num2str(handles.nFiles)]); 
% Initialze Axes Objects
title(handles.axesInt,'Cell');
xlabel(handles.axesInt,'Intensity');
ylabel(handles.axesInt,'Counts');
grid(handles.axesInt,'on');
title(handles.axesBgd,'Background');
xlabel(handles.axesBgd,'Intensity');
ylabel(handles.axesBgd,'Counts');
grid(handles.axesBgd,'on');
% initialize ROI 
handles.roiDataset = [0 0 imInfo.Width imInfo.Height];
% Set bitdepth
handles.bitDepth=imInfo.BitDepth;
% Set zoom reset value for images
infoAxesIm = getappdata(handles.axesIm,'matlab_graphics_resetplotview');
infoAxesIm.XLim =[0 imInfo.Width];
infoAxesIm.YLim =[0 imInfo.Height];
setappdata(handles.axesIm,'matlab_graphics_resetplotview',infoAxesIm);
zoom(handles.axesIm,'reset');
% Save all cell data in a structure stkData
%--- Save filenames
stkData.files     = files';
%--- Save images
imMat = zeros(imInfo.Width, imInfo.Height, numel(files));
for i=1:numel(files)
    tmpIm = imread(files{i});
    imMat(:,:,i) = tmpIm;
end
stkData.im  = imMat;
%--- Save ROI's
stkData.roiCell.x = zeros(handles.nFiles,imInfo.Height);
stkData.roiCell.y = zeros(handles.nFiles,imInfo.Height);
stkData.roiBgd.x  = zeros(handles.nFiles,imInfo.Height);
stkData.roiBgd.y  = zeros(handles.nFiles,imInfo.Height);
%--- Save Results
stkData.results.Isum  = zeros(numel(files),1);%Total Intensity of the cell
stkData.results.Iarea = zeros(numel(files),1);%Area of the cell 
stkData.results.Ibgd  = zeros(numel(files),1);%Mean bgd Intensity
stkData.results.TCF   = zeros(numel(files),1);%Total Corrected Fluorescence
stkData.results.IntDen = zeros(numel(files),1);%Total Inensity per cell area
% Since the image data are stored in ApplicationData
setappdata(handles.figure1,'stkData', stkData);
% Update handles structure
handles.nBins        = 96;
handles.intRange     = intRange; 
handles.sldRange     = sldRange;
handles.newIntRange  = newIntRange;
handles.currentImage = currentImage;
handles.selStat      = 1;
handles.rectRoi      = rectRoiCoord;
guidata(hObject, handles);
% -------------------------------------------------------------------------
%                       EDIT DIRECTORY AND FILENAME
% -------------------------------------------------------------------------
function editDir_Callback(hObject, eventdata, handles)
% hObject    handle to editDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes during object creation, after setting all properties.
function editDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editFilename_Callback(hObject, eventdata, handles)
% hObject    handle to editFilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes during object creation, after setting all properties.
function editFilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% -------------------------------------------------------------------------
%                             NAVIGATE 
% -------------------------------------------------------------------------
% --- Executes on button press in btnNext.
function btnNext_Callback(hObject, eventdata, handles)
% hObject    handle to btnNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set( handles.editId,'String',num2str(handles.fileNo+1) );
editId_Callback( handles.editId, [], handles );
% --- Executes on button press in btnPrevious.
function btnPrevious_Callback(hObject, eventdata, handles)
% hObject    handle to btnPrevious (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set( handles.editId,'String',num2str(handles.fileNo-1) );
editId_Callback( handles.editId, [], handles );
function editId_Callback(hObject, eventdata, handles)
% hObject    handle to editId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Set ID Number
id=str2double(get(hObject,'String'));
% If trace ID is invalid, reset it to what it was before.
if isnan(id) || id > handles.nFiles || id < 1,
    disp('WARNING in viewtraces: Invalid trace number. Resetting.');
    set( hObject,'String',num2str(handles.fileNo) );
    return;
else
    handles.fileNo=id;
end
% Make sure that the selected file actually exists.
if handles.fileNo+1 > handles.nFiles
    set(handles.btnNext,'Enable','off');
else
    set(handles.btnNext,'Enable','on');
end
if handles.fileNo-1 >= 1
    set(handles.btnPrevious,'Enable','on');
else
    set(handles.btnPrevious,'Enable','off');
end
% Clear Axes Objects
cla(handles.axesIm,'reset')
cla(handles.axesInt,'reset')
cla(handles.axesBgd,'reset')
%Switch on grids
grid(handles.axesInt,'on');
grid(handles.axesBgd,'on');
% Load cell data
stkData = getappdata(handles.figure1,'stkData');
% Get image size
filename=stkData.files{handles.fileNo};
info=imfinfo(filename);
imsize=[info.Width info.Height];
% Load image
currentImage=imread(filename);
% Display the original image
intRange = stretchlim(currentImage,handles.Tol);
currentImage    = imadjust(currentImage,intRange);
imagesc(handles.axesIm,currentImage);colormap('gray');
% Apply newIntensity range if exit
if ~isempty(handles.newIntRange)
    currentImageNew    = imadjust(currentImage,handles.newIntRange);
    imagesc(handles.axesIm,currentImageNew);colormap('gray');
end
% Re-initialize axes objects 
title(handles.axesIm,['Cell ' num2str(handles.fileNo) ' of ' num2str(handles.nFiles)]); 
xlabel(handles.axesIm,'x [Pixels]');
ylabel(handles.axesIm,'y [Pixels]');
axis(handles.axesIm,[0 imsize(1,1) 0 imsize(1,2)]);
title(handles.axesInt,'Intensity Distribution in ROI ( Cell )');
xlabel(handles.axesInt,'Intensity');
ylabel(handles.axesInt,'Counts');
grid(handles.axesInt,'on');
title(handles.axesBgd,'Intensity Distribution in ROI ( Background )');
xlabel(handles.axesBgd,'Intensity');
ylabel(handles.axesBgd,'Counts');
grid(handles.axesBgd,'on');
% Update Filename 
[~,name,ext]=fileparts(filename);
set(handles.editFilename,'String',[name ext]);
% Update handles structure
handles.intRange=intRange;
handles.currentImage=currentImage;
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function editId_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% -------------------------------------------------------------------------
%                          ADJUST IMAGE INTENSITY
% -------------------------------------------------------------------------
% --- Executes on slider movement.
function sliderIntensity_Callback(hObject, eventdata, handles)
% hObject    handle to sliderIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer=get(hObject,'Value');
set(handles.sliderIntensity,'Value',answer);
newRange=handles.intRange(2)+answer*0.1;
if newRange>1;newRange=1;end
currentImage=imadjust(handles.currentImage,...
                     [handles.intRange(1) newRange]);
imagesc(currentImage,'Parent', handles.axesIm);
% Axis Labels
xlabel(handles.axesIm, 'x [Pixels]');
ylabel(handles.axesIm, 'y [Pixels]');
title(handles.axesIm, ['Cell ' num2str(handles.fileNo) ' of '...
      num2str(handles.nFiles)]); 
% Update handles structure
handles.IntSliderValue = answer; 
handles.newIntRange    = [0 newRange];
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function sliderIntensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function editToleranz_Callback(hObject, eventdata, handles)
% hObject    handle to editToleranz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load cell data
stkData = getappdata(handles.figure1,'stkData');
% Get user input
answer=str2double(get(hObject,'String'));
handles.Tol = answer;
% Get image size
file=stkData.files{handles.fileNo};
info=imfinfo(file);
imsize=[info.Width info.Height];
% Load image
currentImage=imread(file);
% Display the image
intRange     = stretchlim(currentImage,handles.Tol);
currentImage = imadjust(currentImage,intRange);
imagesc(handles.axesIm, currentImage);colormap('gray');
% Re-Initialze Axes Objects
title(handles.axesIm,['Cell ' num2str(handles.fileNo)...
                      ' of ' num2str(handles.nFiles)]); 
xlabel(handles.axesIm,'x [Pixels]');
ylabel(handles.axesIm,'y [Pixels]');
axis(handles.axesIm,[0 imsize(1,1) 0 imsize(1,2)]);
% Update handles structure
handles.intRange=intRange;
handles.currentImage=currentImage;
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function editToleranz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editToleranz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% -------------------------------------------------------------------------
%                             DRAW ROI's
% -------------------------------------------------------------------------
% --- Executes on button press in btnDrawRoiCell.
function btnDrawRoiCell_Callback(hObject, eventdata, handles)
% hObject    handle to btnDrawRoiCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load cell data
stkData = getappdata(handles.figure1,'stkData');
% Get roi data from user
freeHandregion=imfreehand(handles.axesIm);hold on
pos  = getPosition(freeHandregion);
x    = pos(:,1)';
y    = pos(:,2)';
nPts = size(pos,1);
stkData.roiCell.x(handles.fileNo,1:nPts) = x;
stkData.roiCell.y(handles.fileNo,1:nPts) = y;
% Apply ROI to raw image 
rawIm    = stkData.im(:,:,handles.fileNo);
[width,height]=size(rawIm);
mask     = poly2mask(x,y,width,height);
roiIm    = rawIm.*mask;
Iarea = sum(sum(mask)); 
Isum  = sum(sum(roiIm));
% Plot Intensity Distribution
histInt=histogram(handles.axesInt,roiIm(roiIm>0),handles.nBins);
legend(handles.axesInt,['Area: ' num2str(Iarea)]);
xBinEnd    = histInt.BinEdges(2:end);
xBinCenter = (xBinEnd)-xBinEnd(1)/2;
intDist    = [xBinCenter; histInt.Values]';
assignin('base', 'intDist', intDist);
% Initialze Axes Objects
title(handles.axesInt,'Intensity Distribution in ROI ( Cell )');
xlabel(handles.axesInt,'Intensity');
ylabel(handles.axesInt,'Counts');
grid(handles.axesInt,'on');
% Save results in stkData
stkData.results.Isum(handles.fileNo) = Isum; 
stkData.results.Iarea(handles.fileNo) = Iarea;
% Output results
no=(1:handles.nFiles)';
output = [ no stkData.results.Isum stkData.results.Iarea...
           stkData.results.Ibgd stkData.results.TCF ];        
set( handles.editOutput,'String',num2str(output,...
        '%d   %d   %d   %.1f   %.1f'));
    
% Save Application Data
setappdata(handles.figure1,'stkData', stkData);
function editRectRoiStartPos_Callback(hObject, eventdata, handles)
% hObject    handle to editRectRoiStartPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get user input
answer  = get(hObject,'String');
rectRoiTmp = textscan(answer,'%f %f %f %f');
rectRoi=[rectRoiTmp{:}];
if size(rectRoi,2)==4
    handles.rectRoi = rectRoi;
else
    warning('ROI COORDINATES HAVE WRONG DIMENSION');
    set(handles.editRectRoiStartPos,'String',num2str(handles.rectRoi)); 
    return
end
% Update handles structure
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function editRectRoiStartPos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRectRoiStartPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in btnRectRoiCell.
function btnRectRoiCell_Callback(hObject, eventdata, handles)
% hObject    handle to btnRectRoiCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load cell data
stkData = getappdata(handles.figure1,'stkData');
% Read rectRoi data from the handles struct. 
rectRoiX0     = handles.rectRoi(1);
rectRoiY0     = handles.rectRoi(2);
rectRoiWidth  = handles.rectRoi(3);
rectRoiHeight = handles.rectRoi(4);
% Get roi data from user
h = imrect(gca, [rectRoiX0 rectRoiY0 rectRoiWidth rectRoiHeight]);
addNewPositionCallback(h,@(p) title(mat2str(p,3)));
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(h,fcn); 
pos = wait(h);
if isempty(pos);return;end
x = [pos(1) pos(1)+pos(3) pos(1)+pos(3) pos(1) pos(1)];
y = [pos(2) pos(2) pos(2)+pos(4) pos(2)+pos(4) pos(2)];  
nPts = 5;
stkData.roiCell.x(handles.fileNo,1:nPts) = x;
stkData.roiCell.y(handles.fileNo,1:nPts) = y;
% Apply ROI to raw image 
rawIm    = stkData.im(:,:,handles.fileNo);
[width,height]=size(rawIm);
mask     = poly2mask(x,y,width,height);
roiIm    = rawIm.*mask;
Iarea = sum(sum(mask)); 
Isum  = sum(sum(roiIm));
% Plot Intensity Distribution
hstInt=histogram(handles.axesInt,roiIm(roiIm>0),handles.nBins);
legend(handles.axesInt,['Area: ' num2str(Iarea)]);
% Assign hstInt values to the variable 'hstInt' in the matlab base workspace
assignin('base', 'hstInt', hstInt);
% Initialze Axes Objects
title(handles.axesInt,'Intensity Distribution in ROI ( Cell )');
xlabel(handles.axesInt,'Intensity');
ylabel(handles.axesInt,'Counts');
grid(handles.axesInt,'on');
% Save results in stkData
stkData.results.Isum(handles.fileNo) = Isum; 
stkData.results.Iarea(handles.fileNo) = Iarea;
% Save Application Data
setappdata(handles.figure1,'stkData', stkData);
% Display Results
editOutput_Callback( handles.editOutput, [], handles );
% --- Executes on button press in btnDrawRoiBgd.
function btnDrawRoiBgd_Callback(hObject, eventdata, handles)
% hObject    handle to btnDrawRoiBgd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load cell data
stkData = getappdata(handles.figure1,'stkData');
% Get roi data from user
freeHandregion=imfreehand(handles.axesIm);
pos = getPosition(freeHandregion);
xPos    = pos(:,1)';
yPos    = pos(:,2)';
nPts=size(pos,1);
stkData.roiBgd.x(handles.fileNo,1:nPts)=pos(:,1)';
stkData.roiBgd.y(handles.fileNo,1:nPts)=pos(:,2)';
% Apply ROI to raw image 
rawIm    = stkData.im(:,:,handles.fileNo);
[width,height]=size(rawIm);
mask     = poly2mask(xPos,yPos,width,height);
roiIm    = rawIm.*mask;
bgdPixArea = sum(sum(mask)); 
% Plot Intensity Distribution
edges=0:2:256;
histBgd=histogram(handles.axesBgd,roiIm(roiIm>0),edges);
hold(handles.axesBgd,'on');
xBinEnd    = histBgd.BinEdges(2:end);
xBinCenter = (xBinEnd)-xBinEnd(1)/2;
x          = xBinCenter;
y          = histBgd.Values;
% Assign histBgd values to the variable 'hstBgd' in the matlab base workspace
assignin('base', 'hstBgd', histBgd);
% Fit background distribution 
try 
    gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
    startPoints = [max(y) mean(roiIm(roiIm>0)) std(roiIm(roiIm>0)) 1];  
    f1 = fit(x',y',gaussEqn,'Start', startPoints);
    line(x,f1(x),'Parent',handles.axesBgd,'Color','red');
    stkData.results.Ibgd(handles.fileNo) = f1.b;
    legend(handles.axesBgd,['mean: ' num2str(f1.b); 'width: ' num2str(2*sqrt(f1.c))]);
catch
     disp('fit does not converge: saving approx. bgd mean');
     stkData.results.Ibgd(handles.fileNo) =  mean(roiIm(roiIm>0));
end
% Initialze Axes Objects
title(handles.axesBgd,'Intensity Distribution in ROI ( Background )');
xlabel(handles.axesBgd,'Intensity');
ylabel(handles.axesBgd,'Counts');
grid(handles.axesBgd,'on');
% Save Application Data
setappdata(handles.figure1,'stkData', stkData);
% Display Results
editOutput_Callback( handles.editOutput, [], handles );
% -------------------------------------------------------------------------
%                     CALCULATE LOCAL INTENSITY 
% -------------------------------------------------------------------------
% --- Executes on selection change in popUpSizePsf.
function popUpSizePsf_Callback(hObject, eventdata, handles)
% hObject    handle to popUpSizePsf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get User input
selPsf=get(handles.popUpSizePsf,'Value');
% Settings for the integration over a pixel area
switch selPsf
    case 1
        handles.widthPsf=7; handles.stepPsf=4; %for 7x7 area 
    case 2
        handles.widthPsf=5; handles.stepPsf=3; %for 5x5 area 
    case 3
        handles.widthPsf=3; handles.stepPsf=2; %for 3x3 area 
end
% Update handles structure
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function popUpSizePsf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popUpSizePsf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editNCoordinates_Callback(hObject, eventdata, handles)
% hObject    handle to editNCoordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get user input
answer  = get(hObject,'String');
nMol    = textscan(answer,'%u');
nSelMol = nMol{1};
if nSelMol==0
    warning('WRONG NUMBER');
    set(handles.editNCoordinates,'String',num2str(handles.nSelMol)); 
    return
end
% Update handles structure
handles.nSelMol=nSelMol;
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function editNCoordinates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNCoordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in btnGetCoordinates.
function btnGetCoordinates_Callback(hObject, eventdata, handles)
% hObject    handle to btnGetCoordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Clear Axes Objects
cla(handles.axesInt,'reset')
cla(handles.axesBgd,'reset')
%Switch on grids
grid(handles.axesInt,'on');
grid(handles.axesBgd,'on');
% Load cell data
stkData = getappdata(handles.figure1,'stkData');
% Initialize Parameters
widthPsf         = handles.widthPsf;
step             = handles.stepPsf;    
intPsf           = zeros(widthPsf,widthPsf)*NaN;
im               = handles.currentImage;
imRaw            = stkData.im(:,:,handles.fileNo);
[imHeigth, imWidth] = size(im);
QE               = handles.QE;
photonConversion = handles.photonConversion;
locInt           = [];
locBgd           = [];
allPht           = zeros(handles.nSelMol,1);
locSnr           = zeros(handles.nSelMol,1);
locSumPht        = zeros(handles.nSelMol,1);
locSumBgdPht     = zeros(handles.nSelMol,1);
% Create an imref2d object and associate it with the image.
R      = imref2d(size(im));
% Plot Image3
ax=handles.axesIm;
cla(ax);hold(ax,'on');
if ~isempty(handles.newIntRange)
    currentImageNew    = imadjust(im,handles.newIntRange);
    imagesc(ax,currentImageNew);colormap('gray');
else
    %currentImage   = imadjust(im,handles.intRange);
    imagesc(ax,im);colormap('gray');
end
for h = 1:handles.nSelMol
    
    % User input: Select coordinates from image
    [xCor,yCor] = ginput(1);
    
    % Convert world coordinates to subscripts
    [rowIdx, colIdx] = worldToSubscript(R,xCor,yCor);
    
    %--- Measure the intensity over the pixel area: for test purposes
    %--- use the function 'sumPixelArea.m'
    for i=1:widthPsf
        for j=1:widthPsf
            col=colIdx+i-step;
            row=rowIdx+j-step;
            plot(ax,col,row,'.','Color','r'); hold on
            if col>0 && row >0
                if col<=imHeigth && row <=imWidth
                    intPsf(i,j) = imRaw(row,col);
                end
            end
        end
    end
    allPht(h)=sum(sum(intPsf));
    
    % Rotate and flip matrirx (to have the same orientation than the image)  
    intPsf       = fliplr(rot90(intPsf,3));
    locIntTmp    = intPsf(2:widthPsf-1,2:widthPsf-1);
    locBgdTmp    = [intPsf(1,1:widthPsf);...                % upper row
                   intPsf(1:widthPsf,widthPsf)'; ...        % rigth cloumn
                   fliplr(intPsf(widthPsf,1:widthPsf));...  % bottom row
                   fliplr(intPsf(1:widthPsf,1)')];          % left cloumn
    locBgdTmp(:,1) = [];
    
    %Claculate the local SNR    
    %--- Integrate Intensities in the central region and in the
    %--- boundary region
    if sum(sum(isnan(intPsf)))>0
        disp('PSF Error: Found NAN Values');
    else
        %--- Intgrate intensity values over the central pixel area of 
        %--- the psf 
       
        photons    = sum(sum(locIntTmp))/QE;
        photons    = photons./photonConversion;
        bgdPsfTmp  = reshape(locBgdTmp,1,numel(locBgdTmp));
        %--- Elements of the bgd matrix need to match the number of
        %--- matrix elements of the intPsf. 
        %   size        size of                   number
        % widthPSF    intPsf kernel         surrounding pixels
        %    3x3     1x1 (1 element)   8  elem (choose 1 elem randm)
        %    5x5     3x3 (9 elements)  16 elem (choose 9 elem randm)
        %    7x7     5x5 (25 elements) 24 elem (choose all elem +1 randm)
        if widthPsf == 3
            bgdPsfTmp=bgdPsfTmp(5:8);
            bgdPsf=mean(bgdPsfTmp);
        elseif widthPsf == 5
            bgdPsfTmp=[bgdPsfTmp(1:4), bgdPsfTmp(9:12)];
            bgdPsf=[bgdPsfTmp,mean(bgdPsfTmp)];
        elseif widthPsf == 7
            bgdPsf=[bgdPsfTmp, mean(bgdPsfTmp)];
        end    
        photonsBgd  = sum(bgdPsf);
        photonsBgd  = photonsBgd./photonConversion;
        %--- test if the number of matrix elements are equal for the  
        %--- intPsf and bgdPsf before calculating the S/N ratio
        nIntPsf = numel(intPsf(2:widthPsf-1,2:widthPsf-1));
        nBgdPsf = numel(bgdPsf);
        if nIntPsf == nBgdPsf
            locSumPht(h) = photons;
            locSumBgdPht(h) = photonsBgd;
            locSnr(h) = photons/photonsBgd;
        end
    end
    
    % Concatenate arrays 
    locInt=cat(2,locInt,locIntTmp);
    locBgd=cat(2,locBgd,locBgdTmp);
end
% Distribution of the local background
edges=0:1:256;
hstLocBgd=histogram(handles.axesBgd,locBgd,edges);
xBinEnd    = hstLocBgd.BinEdges(2:end);
xBinCenter = (xBinEnd)-xBinEnd(1)/2;
% Label Axes local Background
title(handles.axesBgd,'Local Background Distribution');
xlabel(handles.axesBgd,'Intensity');
ylabel(handles.axesBgd,'Counts');
grid(handles.axesBgd,'on');
% Distribution of the local intensity
hstLocInt = histogram(handles.axesInt,locInt,edges);
% Label Axes local intensity
title(handles.axesInt,'Local Intensity Distribution');
xlabel(handles.axesInt,'Intensity');
ylabel(handles.axesInt,'Counts');
grid(handles.axesInt,'on');
% Assign variables to the matlab base workspace
locDist = [xBinCenter; hstLocBgd.Values; hstLocInt.Values]';
assignin('base', 'locDist', locDist);
assignin('base', 'locSnr', locSnr);
% Save Application Data
stkData.results.allPht       = allPht; 
stkData.results.locSumPht    = locSumPht;
stkData.results.locSumBgdPht = locSumBgdPht;
stkData.results.locSnr       = locSnr;
stkData.results.locDist      = locDist;
setappdata(handles.figure1,'stkData', stkData);
% -------------------------------------------------------------------------
%                          OUTPUT RESULTS
% -------------------------------------------------------------------------
function editOutput_Callback(hObject, eventdata, handles)
% hObject    handle to editOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load cell data
stkData = getappdata(handles.figure1,'stkData');
    
% Output results
Isum  = stkData.results.Isum(handles.fileNo); %Total Intensity of the cell
Iarea = stkData.results.Iarea(handles.fileNo);%Area of the cell 
Ibgd  = stkData.results.Ibgd(handles.fileNo); %Mean bgd Intensity
if Isum > 0 && Ibgd > 0 
    stkData.results.TCF(handles.fileNo) = ...
                           Isum - Iarea*Ibgd; %Total Corrected Fluorescence
end
if Isum > 0 && Ibgd > 0
    stkData.results.IntDen(handles.fileNo) = ...
                           (Isum - Iarea*Ibgd) / Iarea; % TCF Density
end
% Save Application Data
setappdata(handles.figure1,'stkData', stkData);
% Display Output
no=(1:handles.nFiles)';
output = [ no stkData.results.Isum stkData.results.Iarea...
           stkData.results.Ibgd stkData.results.TCF ]; 
    set( handles.editOutput,'String',num2str(output,...
       '%d   %d   %d   %.1f   %.1f'));
   
% Save Application Data
setappdata(handles.figure1,'stkData', stkData);
% --- Executes during object creation, after setting all properties.
function editOutput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% -------------------------------------------------------------------------
%                               HISTOGRAM 
% -------------------------------------------------------------------------
% --- Executes on selection change in listStat.
function listStat_Callback(hObject, eventdata, handles)
% hObject    handle to listStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
handles.selStat = contents{get(hObject,'Value')};
% Update handles structure
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function listStat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in btnHist.
function btnHist_Callback(hObject, eventdata, handles)
% hObject    handle to btnHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load cell data
stkData = getappdata(handles.figure1,'stkData');
figure;
switch handles.selStat
    
    case 'TCF / Area'
        bar(stkData.results.IntDen)
        title(gca,'Total Corrected Fluorescence ( TCF ) / Cell Area ');
        ylabel('TCF / Area');
        grid on
        
    case 'Sum I'
        bar(stkData.results.Isum)
        title(gca,'Integrated Intensity over Cell Area');
        ylabel(' Total Intensity ');
        grid on
        
    case 'Area'
        bar(stkData.results.Iarea)
        title(gca,'Cell Area');
        ylabel('Area (Pixel)^2');
        grid on
        
    case '<Bgd>'
        bar(stkData.results.Ibgd)
        title(gca,'Mean Background Intensity');
        ylabel('< Bgd >');
        grid on
    
    case 'TCF'
        bar(stkData.results.TCF)
        title(gca,'Total Corrected Fluorescence ( TCF )');
        ylabel('TCF');
        grid on
        
end
xlabel('Cell');
% -------------------------------------------------------------------------
%                           SAVE DATA TO FILE 
% -------------------------------------------------------------------------
% --- Executes on button press in btnSaveData.
function btnSaveData_Callback(hObject, eventdata, handles)
% hObject    handle to btnSaveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load cell data
stkData = getappdata(handles.figure1,'stkData'); %#ok<NASGU>
% Save data in a mat file
newFilename = 'cellFluor.mat';
outfile     = [handles.basepath filesep newFilename];
[filename,path] = uiputfile(outfile,'Save As...');
save([path filename],'-mat','stkData');