function varargout = Screen_Data_Interface(varargin)
% SCREEN_DATA_INTERFACE MATLAB code for Screen_Data_Interface.fig
%      SCREEN_DATA_INTERFACE, by itself, creates a new SCREEN_DATA_INTERFACE or raises the existing
%      singleton*.
%
%      H = SCREEN_DATA_INTERFACE returns the handle to a new SCREEN_DATA_INTERFACE or the handle to
%      the existing singleton*.
%
%      SCREEN_DATA_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCREEN_DATA_INTERFACE.M with the given input arguments.
%
%      SCREEN_DATA_INTERFACE('Property','Value',...) creates a new SCREEN_DATA_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Screen_Data_Interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Screen_Data_Interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Screen_Data_Interface

% Last Modified by GUIDE v2.5 04-Jan-2018 23:07:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Screen_Data_Interface_OpeningFcn, ...
                   'gui_OutputFcn',  @Screen_Data_Interface_OutputFcn, ...
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


% --- Executes just before Screen_Data_Interface is made visible.
function Screen_Data_Interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Screen_Data_Interface (see VARARGIN)

load('allBehaveMat.mat');
handles.dataCell = allDataCell;
load('allBehaveColNames.mat');
handles.cellNames = cellColNames;



% Choose default command line output for Screen_Data_Interface
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Screen_Data_Interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Screen_Data_Interface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in TypeBox.
function TypeBox_Callback(hObject, eventdata, handles)
% hObject    handle to TypeBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TypeBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TypeBox

str = cellstr(get(hObject,'String'));
handles.types = str(get(hObject,'Value'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function TypeBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TypeBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

annoFilepath = '/Users/kyobikakaria/Desktop/Data/split_screen/rubinAnnotation.csv';
importedAnno = importdata(annoFilepath);
fill = [importedAnno.textdata(2:end,1);importedAnno.textdata(1,2:end)'];

set(hObject,'String', fill, 'Max', length(fill),'Min', 1)

handles.importedAnno = importedAnno;
guidata(hObject,handles);




% --- Executes on selection change in AnalysisTypePop.
function AnalysisTypePop_Callback(hObject, eventdata, handles)
% hObject    handle to AnalysisTypePop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AnalysisTypePop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AnalysisTypePop

% Determine the selected data set.
handles.AnalysisType = get(hObject,'Value');
% Save the handles structure.
guidata(hObject,handles)




% --- Executes during object creation, after setting all properties.
function AnalysisTypePop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnalysisTypePop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in AnalysisGrpPop.
function AnalysisGrpPop_Callback(hObject, eventdata, handles)
% hObject    handle to AnalysisGrpPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AnalysisGrpPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AnalysisGrpPop
% Determine the selected data set.
handles.AnalysisGrp = get(hObject,'Value');
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function AnalysisGrpPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnalysisGrpPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2




% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in SubTurnsPop.
function SubTurnsPop_Callback(hObject, eventdata, handles)
% hObject    handle to SubTurnsPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SubTurnsPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SubTurnsPop
% Determine the selected data set.
handles.SubTurns = get(hObject,'Value');
% Save the handles structure.
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function SubTurnsPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SubTurnsPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StatisticChoice.
function StatisticChoice_Callback(hObject, eventdata, handles)
% hObject    handle to StatisticChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StatisticChoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StatisticChoice
handles.statistics = get(hObject,'Value');
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function StatisticChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StatisticChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Max', 5,'Min', 1)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


switch handles.AnalysisType
    case 1 
        handles.current_data = guiAnalysisFunction_1(handles);
    case 2
        handles.current_data = guiAnalysisFunction_2(handles);
    case 3
        handles.current_data = guiAnalysisFunction_3(handles);
end

guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on selection change in ZscoreOptionPop.
function ZscoreOptionPop_Callback(hObject, eventdata, handles)
% hObject    handle to ZscoreOptionPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ZscoreOptionPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ZscoreOptionPop

% Determine the selected data set.
handles.ZscoreOption = get(hObject,'Value');
% Save the handles structure.
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ZscoreOptionPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZscoreOptionPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TempSelectPop.
function TempSelectPop_Callback(hObject, eventdata, handles)
% hObject    handle to TempSelectPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TempSelectPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TempSelectPop
% Determine the selected data set.
handles.TempSelect = get(hObject,'Value');
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function TempSelectPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TempSelectPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
