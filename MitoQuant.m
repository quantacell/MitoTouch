% ----------------------------------------
% Copyright Clariant, CNRS, Université de Montpellier, IRD, EPHE
% Contributors: Abdel Aouacheria, Victor Racine
% Contacts : abdel.aouacheria@umontpellier.fr, victor.racine@quantacell.com 

%  
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License v3.0 only.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
%  
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see https://spdx.org/licenses/GPL-3.0-only.html
% ----------------------------------------

function varargout = MitoQuant(varargin)
% MITOQUANT MATLAB code for MitoQuant.fig
%      MITOQUANT, by itself, creates a new MITOQUANT or raises the existing
%      singleton*.
%
%      H = MITOQUANT returns the handle to a new MITOQUANT or the handle to
%      the existing singleton*.
%
%      MITOQUANT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MITOQUANT.M with the given input arguments.
%
%      MITOQUANT('Property','Value',...) creates a new MITOQUANT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MitoQuant_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MitoQuant_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MitoQuant

% Last Modified by GUIDE v2.5 25-Mar-2019 15:45:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MitoQuant_OpeningFcn, ...
                   'gui_OutputFcn',  @MitoQuant_OutputFcn, ...
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


% --- Executes just before MitoQuant is made visible.
function MitoQuant_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MitoQuant (see VARARGIN)

% Choose default command line output for MitoQuant
handles.output = hObject;
handles.userdata=[];
handles.userdata.imageFilename='';
% Update handles structure
handles.text2.String=['Nuclei threshold ', num2str(round(handles.slider1.Value))];
handles.text3.String=['Cell threshold ', num2str(round(handles.slider2.Value))];
handles.text4.String=['Cell erode ', num2str(round(handles.slider3.Value))];
handles.text5.String=['Filter border cells ', num2str(handles.slider4.Value)];

handles.text6.String=['Mito threshold ', num2str(handles.slider5.Value/10)];
axis off;
guidata(hObject, handles);

% UIWAIT makes MitoQuant wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MitoQuant_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.text2.String=['Nuclei threshold ', num2str(round(get(hObject,'Value')))];

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.text3.String=['Cell threshold ', num2str(round(get(hObject,'Value')))];

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.text4.String=['Cell erode ', num2str(round(get(hObject,'Value')))];

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=round(get(hObject,'Value')*100); 
handles.text5.String=['Filter border cells ', num2str(val), '%'];

% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1, 'pointer', 'watch')
drawnow;


handles=guidata(hObject);
mergeFilename=handles.userdata.imageFilename;
if numel(mergeFilename)==0
    set(handles.figure1, 'pointer', 'arrow')
    drawnow;
    msgbox('Please Load a test image first', 'modal');
    return;
end
outputFilename=mergeFilename;
[path, file, ext]=fileparts(outputFilename);
file=strrep(file, '-MERGE', '-res');
outputFilename=[]; fullfile(path, file);
nucleiThreshold=round(handles.slider1.Value);
cellThreshold=round(handles.slider2.Value);
erodeDist=round(handles.slider3.Value);
filterBorder=handles.slider4.Value;
isCellBased=1-handles.checkbox1.Value;

mitoThreshold=handles.slider5.Value/10;

[data, dataHeaders, T, mask, maskName]=processImage(mergeFilename, outputFilename, nucleiThreshold, cellThreshold, mitoThreshold, erodeDist, filterBorder, isCellBased);
setResultsToUserdata(hObject, mask, maskName);

handles.radiobutton6.Value=1;

set(handles.figure1, 'pointer', 'arrow')
drawnow;

updateDisplay(hObject);

function setResultsToUserdata(hObject, mask, maskName)
handles=guidata(hObject);
handles.userdata.maskNuclei=getImageInList(mask, maskName, 'nuclei');
handles.userdata.maskCells=getImageInList(mask, maskName, 'CellsMask');
handles.userdata.maskMito=getImageInList(mask, maskName, 'MitoClusterMask');

handles.userdata.maskSkelMito=getImageInList(mask, maskName, 'MitoSkelMask');
guidata(hObject, handles);


function im=getImageInList(imList, imName, key)

im=[];
for i=1:numel(imList),
    if strcmpi(imName{i}, key)
        im=imList{i};
        return;
    end
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.figure1.Name='MitoTouch Processing ...';
handles.pushbutton2.String='Processing';
set(handles.figure1, 'pointer', 'watch')
drawnow;

nucleiThreshold=round(handles.slider1.Value);
cellThreshold=round(handles.slider2.Value);
erodeDist=round(handles.slider3.Value);
filterBorder=handles.slider4.Value;
isCellBased=1-handles.checkbox1.Value;
isParallel=handles.checkbox3.Value;
mitoThreshold=handles.slider5.Value/10;

[~, cptImage]=processAllImages([], [],[],[], nucleiThreshold, cellThreshold, mitoThreshold, erodeDist, filterBorder, isCellBased, isParallel);
if numel(cptImage)>0
    msgbox([num2str(cptImage), ' images were processed'],'modal');
end

handles.figure1.Name='MitoTouch';
handles.pushbutton2.String='Run';
set(handles.figure1, 'pointer', 'arrow')
drawnow;


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( {'*.tif', 'TIF image'; '*.tiff', 'TIFF image'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Pick a file', handles.userdata.imageFilename);
if isequal(filename,0)
    return;
end

set(handles.figure1, 'pointer', 'watch')
loadImage(hObject, handles, fullfile(pathname, filename));
set(handles.figure1, 'pointer', 'arrow')

updateDisplay(hObject);

function loadImage(hObject, handles, filename)
handles.userdata.imageFilename=filename;
if strcmpi(filename(end-4:end), '.tiff')
    image1=imread(filename);
    filenameChan2=strrep(filename, 'ch1', 'ch2');
    image2=imread(filenameChan2);
    filenameChan3=strrep(filename, 'ch1', 'ch3');
    image3=imread(filenameChan3);
    handles.userdata.image=cat(3, image1, image3, image2);
else
    
    handles.userdata.image=imread3D(handles.userdata.imageFilename);
end

handles.userdata.maskNuclei=[];
handles.userdata.maskCells=[];
handles.userdata.maskMito=[];
handles.userdata.maskSkelMito=[];
guidata(hObject, handles);

function rgb=autolevel(image)
image=double(image);
rgb=uint8(zeros(size(image, 1), size(image, 2), 3));
for i=1:size(image, 3)
    im0=image(:,:,i);
    mini=prctile(im0(:),30);
    maxi=prctile(im0(:),99.9);
    rgb(:,:,i)=uint8(255*(im0-mini)/(maxi-mini));
    
end
function updateDisplay(hObject)
handles=guidata(hObject);

set(handles.figure1, 'pointer', 'watch')
drawnow;

if handles.radiobutton1.Value==1 % original 
    imagesc(autolevel(handles.userdata.image));
    daspect([1 1 1]);
end

if handles.radiobutton2.Value==1 % Nuclei
    imagesc(autolevel(handles.userdata.image));
    drawRegions(handles.userdata.maskNuclei, [], false, [1,1,0]);
    daspect([1 1 1]);
end

if handles.radiobutton3.Value==1 % Cell
    imagesc(autolevel(handles.userdata.image));
    drawRegions(handles.userdata.maskCells, [], false, [1,0,1]);
    daspect([1 1 1]);
end

if handles.radiobutton4.Value==1 % Mito
    imagesc(autolevel(handles.userdata.image));
    drawRegions(handles.userdata.maskMito, [], false);
    daspect([1 1 1]);
    
end

if handles.radiobutton5.Value==1 % Skeleton
    imagesc(autolevel(handles.userdata.image));
    drawRegions(imdilate(handles.userdata.maskSkelMito, ones(3)), [], false, [1,0,1]);
    %drawRegions(handles.userdata.maskSkelMito, [], false);
    daspect([1 1 1]);
end

if handles.radiobutton6.Value==1 % Final
    imagesc(autolevel(handles.userdata.image));
    drawRegions(handles.userdata.maskCells, [], false, [1,0,1]);
    drawRegions(handles.userdata.maskMito, [], false);
    drawRegions(handles.userdata.maskNuclei, [], false, [1,1,0]);
    daspect([1 1 1]);
end
axis off

daspect([1 1 1]);
axis off

set(handles.figure1, 'pointer', 'arrow')
drawnow;

% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
updateDisplay(hObject);


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
updateDisplay(hObject);


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
updateDisplay(hObject);


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
updateDisplay(hObject);


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5
updateDisplay(hObject);


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6
updateDisplay(hObject);


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

posF=get(handles.figure1, 'Position');
pos=[37.3750    1.1392 , posF(3)-40, posF(4)-3];
set(handles.axes1, 'Position', pos);


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes during object creation, after setting all properties.
function pushbutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val=round(get(hObject,'Value')*10)/100; 
handles.text6.String=['Mito threshold ', num2str(val)];


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
