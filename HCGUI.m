function varargout = HCGUI(varargin)
% HCGUI M-file for HCGUI.fig
%      HCGUI, by itself, creates a new HCGUI or raises the existing
%      singleton*.
%
%      H = HCGUI returns the handle to a new HCGUI or the handle to
%      the existing singleton*.
%
%      HCGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HCGUI.M with the given input arguments.
%
%      HCGUI('Property','Value',...) creates a new HCGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HCGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HCGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help HCGUI

% Last Modified by GUIDE v2.5 26-Jun-2022 18:52:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HCGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @HCGUI_OutputFcn, ...
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


% --- Executes just before HCGUI is made visible.
function HCGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HCGUI (see VARARGIN)

% Choose default command line output for HCGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS

cd (curdir)
load Auditory_handles
set(handles.play_3_tokenHCT,'Value',HCPARAMS.play_3_tokenHCT);

eval(sprintf('set(handles.frame_,''String'',''%5.3f'');',HCPARAMS.frame_dur))
set(handles.frame_,'Value',HCPARAMS.frame_dur);
%eval(sprintf('set(handles.exp_dur,''String'',''%5.3f'');',HCPARAMS.exp_dur))
%set(handles.exp_dur,'Value',HCPARAMS.exp_dur);
eval(sprintf('set(handles.gap_dur,''String'',''%5.3f'');',HCPARAMS.gap_dur))
set(handles.gap_dur,'Value',HCPARAMS.gap_dur);
eval(sprintf('set(handles.pip_dur,''String'',''%5.3f'');',HCPARAMS.pip_dur))
set(handles.pip_dur,'Value',HCPARAMS.pip_dur);
eval(sprintf('set(handles.total_frames_lines,''String'',''%i'');',HCPARAMS.total_frames_lines))
set(handles.total_frames_lines,'Value',HCPARAMS.total_frames_lines);
eval(sprintf('set(handles.stim_start,''String'',''%i'');',HCPARAMS.stim_start))
set(handles.stim_start,'Value',HCPARAMS.stim_start);
eval(sprintf('set(handles.ipi_val,''String'',''%8.3f'');',HCPARAMS.ipi_val))
set(handles.ipi_val,'Value',HCPARAMS.ipi_val);
eval(sprintf('set(handles.freq_lo,''String'',''%8.1f'');',HCPARAMS.freq_lo))
set(handles.freq_lo,'Value',HCPARAMS.freq_lo);
eval(sprintf('set(handles.freq_hi,''String'',''%8.1f'');',HCPARAMS.freq_hi))
set(handles.freq_hi,'Value',HCPARAMS.freq_hi);
eval(sprintf('set(handles.num_freq_steps,''String'',''%i'');',HCPARAMS.num_freq_steps))
set(handles.num_freq_steps,'Value',HCPARAMS.num_freq_steps);
eval(sprintf('set(handles.level_lo,''String'',''%3.1f'');',HCPARAMS.level_lo))
set(handles.level_lo,'Value',HCPARAMS.level_lo);
eval(sprintf('set(handles.level_hi,''String'',''%3.1f'');',HCPARAMS.level_hi))
set(handles.level_hi,'Value',HCPARAMS.level_hi);
eval(sprintf('set(handles.num_att_steps,''String'',''%i'');',HCPARAMS.num_att_steps))
set(handles.num_att_steps,'Value',HCPARAMS.num_att_steps);
eval(sprintf('set(handles.npips_per_train,''String'',''%i'');',HCPARAMS.npips_per_train))
set(handles.npips_per_train,'Value',HCPARAMS.npips_per_train);
eval(sprintf('set(handles.total_reps,''String'',''%i'');',HCPARAMS.total_reps))
set(handles.total_reps,'Value',HCPARAMS.total_reps);

%% FOR RANDOM HC
eval(sprintf('set(handles.numhc,''String'',''%i'');',HCPARAMS.numhc))
set(handles.numhc,'Value',HCPARAMS.numhc);
eval(sprintf('set(handles.freq_f0,''String'',''%i'');',HCPARAMS.freq_f0))
set(handles.freq_f0,'Value',HCPARAMS.freq_f0);
%%edited by ann 

%handles.text=str2double(HCPARAMS.str_inter_tok_int_ms);
% %%
%  eval(sprintf('set(handles.text,''String'',''%s'');',HCPARAMS.text))
%  set(handles.text,'Value',HCPARAMS.text);
%
% HCPARAMS.freq_f0=get(handles.freq_f0,'Value');freq_f0=HCPARAMS.freq_f0;
% HCPARAMS.text=get(handles.text,'Value');str_inter_tok_int_ms=HCPARAMS.text;
% text=str2num(str_inter_tok_int_ms);

% UIWAIT makes HCGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HCGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;



function frame__Callback(hObject, eventdata, handles)
% hObject    handle to frame_ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_ as text
%        str2double(get(hObject,'String')) returns contents of frame_ as a double
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.frame_dur=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.frame_,'Value',HCPARAMS.frame_dur);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function frame__CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pip_dur_Callback(hObject, eventdata, handles)
% hObject    handle to pip_dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pip_dur as text
%        str2double(get(hObject,'String')) returns contents of pip_dur as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.pip_dur=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.pip_dur,'Value',HCPARAMS.pip_dur);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function pip_dur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pip_dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function total_frames_lines_Callback(hObject, eventdata, handles)
% hObject    handle to total_frames_lines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of total_frames_lines as text
%        str2double(get(hObject,'String')) returns contents of total_frames_lines as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)

HCPARAMS.total_frames_lines=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.total_frames_lines,'Value',HCPARAMS.total_frames_lines);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function total_frames_lines_CreateFcn(hObject, eventdata, handles)
% hObject    handle to total_frames_lines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stim_start_Callback(hObject, eventdata, handles)
% hObject    handle to stim_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim_start as text
%        str2double(get(hObject,'String')) returns contents of stim_start as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.stim_start=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.stim_start,'Value',HCPARAMS.stim_start);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function stim_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ipi_val_Callback(hObject, eventdata, handles)
% hObject    handle to ipi_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ipi_val as text
%        str2double(get(hObject,'String')) returns contents of ipi_val as a double


curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.ipi_val=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.ipi_val,'Value',HCPARAMS.ipi_val);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function ipi_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ipi_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function freq_lo_Callback(hObject, eventdata, handles)
% hObject    handle to freq_lo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq_lo as text
%        str2double(get(hObject,'String')) returns contents of freq_lo as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.freq_lo=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.freq_lo,'Value',HCPARAMS.freq_lo);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function freq_lo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_lo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function freq_hi_Callback(hObject, eventdata, handles)
% hObject    handle to freq_hi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq_hi as text
%        str2double(get(hObject,'String')) returns contents of freq_hi as a double

load HCPARAMScurdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.freq_hi=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.freq_hi,'Value',HCPARAMS.freq_hi);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function freq_hi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_hi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_freq_steps_Callback(hObject, eventdata, handles)
% hObject    handle to num_freq_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_freq_steps as text
%        str2double(get(hObject,'String')) returns contents of num_freq_steps as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.num_freq_steps=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.num_freq_steps,'Value',HCPARAMS.num_freq_steps);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function num_freq_steps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_freq_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function level_lo_Callback(hObject, eventdata, handles)
% hObject    handle to level_lo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of level_lo as text
%        str2double(get(hObject,'String')) returns contents of level_lo as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.level_lo=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.level_lo,'Value',HCPARAMS.level_lo);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function level_lo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to level_lo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function level_hi_Callback(hObject, eventdata, handles)
% hObject    handle to level_hi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of level_hi as text
%        str2double(get(hObject,'String')) returns contents of level_hi as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.level_hi=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.level_hi,'Value',HCPARAMS.level_hi);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function level_hi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to level_hi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_att_steps_Callback(hObject, eventdata, handles)
% hObject    handle to num_att_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_att_steps as text
%        str2double(get(hObject,'String')) returns contents of num_att_steps as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.num_att_steps=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.num_att_steps,'Value',HCPARAMS.num_att_steps);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function num_att_steps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_att_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function npips_per_train_Callback(hObject, eventdata, handles)
% hObject    handle to npips_per_train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of npips_per_train as text
%        str2double(get(hObject,'String')) returns contents of npips_per_train as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.npips_per_train=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.npips_per_train,'Value',HCPARAMS.npips_per_train);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function npips_per_train_CreateFcn(hObject, eventdata, handles)
% hObject    handle to npips_per_train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function total_reps_Callback(hObject, eventdata, handles)
% hObject    handle to total_reps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of total_reps as text
%        str2double(get(hObject,'String')) returns contents of total_reps as a double
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
load HCPARAMS
HCPARAMS.total_reps=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.total_reps,'Value',HCPARAMS.total_reps);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function total_reps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to total_reps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in done_button.
function done_button_Callback(hObject, eventdata, handles)
% hObject    handle to done_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)


HCPARAMS.frame_dur=get(handles.frame_,'Value');fdur=HCPARAMS.frame_dur;
%HCPARAMS.exp_dur=get(handles.exp_dur,'Value');edur=HCPARAMS.exp_dur;
HCPARAMS.gap_dur=get(handles.gap_dur,'Value');gapdur=HCPARAMS.gap_dur;
HCPARAMS.pip_dur=get(handles.pip_dur,'Value');
HCPARAMS.total_frames_lines=get(handles.total_frames_lines,'Value');tfl=HCPARAMS.total_frames_lines;
HCPARAMS.stim_start=get(handles.stim_start,'Value');sf=HCPARAMS.stim_start;
HCPARAMS.ipi_val=get(handles.ipi_val,'Value');ipi=HCPARAMS.ipi_val;
HCPARAMS.freq_lo=get(handles.freq_lo,'Value');fl=HCPARAMS.freq_lo;
HCPARAMS.freq_hi=get(handles.freq_hi,'Value');fh=HCPARAMS.freq_hi;
HCPARAMS.num_freq_steps=get(handles.num_freq_steps,'Value');nf=HCPARAMS.num_freq_steps;
HCPARAMS.level_lo=get(handles.level_lo,'Value');ll=HCPARAMS.level_lo;
HCPARAMS.level_hi=get(handles.level_hi,'Value');lh=HCPARAMS.level_hi;
HCPARAMS.num_att_steps=get(handles.num_att_steps,'Value');na=HCPARAMS.num_att_steps;
HCPARAMS.npips_per_train=get(handles.npips_per_train,'Value');nppt=HCPARAMS.npips_per_train;
HCPARAMS.total_reps=get(handles.total_reps,'Value');treps=HCPARAMS.total_reps;

HCPARAMS.numhc=get(handles.numhc,'Value');numhc=HCPARAMS.numhc;
HCPARAMS.freq_f0=get(handles.freq_f0,'Value');freq_f0=HCPARAMS.freq_f0;

%HCPARAMS.text=get(handles.text,'Value');
% str_inter_tok_int_ms=HCPARAMS.text;
%  HCPARAMS.text=[60 90 150 280];
%text=str2num(HCPARAMS.str_inter_tok_int_ms);

HCPARAMS.play_3_tokenHCT=get(handles.play_3_tokenHCT,'Value');
play_3_tokenHCT=HCPARAMS.play_3_tokenHCT;

if play_3_tokenHCT==1
    set(handles.Cover3,'Visible','Off')
    set(handles.Cover4,'Visible','Off')
    set(handles.Cover1,'Visible','Off')
    set(handles.Cover2,'Visible','Off')
elseif play_3_tokenHCT==0
    set(handles.Cover3,'Visible','On')
    set(handles.Cover4,'Visible','On')
    set(handles.Cover1,'Visible','Off')
    set(handles.Cover2,'Visible','Off')
end

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
%cd (curdir)

load AudStimPARAMS
AudStimPARAMS.type='HC';
AudStimPARAMS.stim_protocol=HCPARAMS;
cd(curdir);
close(handles.figure1)
total_prairie_iterations=nf*na*treps;
each_stim_dur=ipi*nppt;
each_prairie_dur=tfl*fdur;
lines_per_set=nf*na;

nframes_per_iter=sf+(fix((each_stim_dur+gapdur)/fdur)+1)*total_prairie_iterations;
if tfl==nframes_per_iter
    durcheck_str='MicroManager: Number of frames OKAY';
    disp_bgcolor=[0 .8 1];
else
    eval(sprintf('durcheck_str=''MicroManager: Set frames to %i'';',nframes_per_iter))
    disp_bgcolor=[1 0 0];
end
eval(sprintf('add_dur_check='' ...[Must be %i : currently %i]'';',nframes_per_iter,tfl))
eval(sprintf('itercheck_str=''No Iteration check for MicroManager'';'))
display_str=strcat(durcheck_str,add_dur_check,'*\*/\*/*',itercheck_str);
load HC_handles
set(AA_handles.check_disp,'String',display_str)
set(AA_handles.check_disp,'BackgroundColor',disp_bgcolor)
set(AA_handles.check_disp,'FontSize',9)
eval(sprintf('prot_str=''HC: %i freqs from %8.1f to %8.1f Hz, %i atten steps from %4.1f to %4.1f dB, %i Tone pips per train, Starts at frame %i'';',nf,fl,fh,na,ll,lh,nppt,sf))
set(AA_handles.disp_protocol,'String',prot_str)

AudStimPARAMS.protocol_str=prot_str;
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save AudStimPARAMS AudStimPARAMS
cd (curdir)

%save AudStimPARAMS AudStimPARAMS

function exp_dur_Callback(hObject, eventdata, handles)
% hObject    handle to exp_dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of exp_dur as text
%        str2double(get(hObject,'String')) returns contents of exp_dur as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.exp_dur=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.exp_dur,'Value',HCPARAMS.exp_dur);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

%--- Executes during object creation, after setting all properties.
function exp_dur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exp_dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gap_dur_Callback(hObject, eventdata, handles)
% hObject    handle to gap_dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gap_dur as text
%        str2double(get(hObject,'String')) returns contents of gap_dur as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.gap_dur=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.gap_dur,'Value',HCPARAMS.gap_dur);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

% --- Executes during object creation, after setting all properties.
function gap_dur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gap_dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function numhc_Callback(hObject, eventdata, handles)
% hObject    handle to numhc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numhc as text
%        str2double(get(hObject,'String')) returns contents of numhc as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)
HCPARAMS.numhc=str2double(get(hObject,'String')); %returns contents of num_pt_val as a double
set(handles.numhc,'Value',HCPARAMS.numhc);
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)


% --- Executes during object creation, after setting all properties.
function numhc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numhc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in play_3_tokenHCT.
function play_3_tokenHCT_Callback(hObject, eventdata, handles)
% hObject    handle to play_3_tokenHCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of play_3_tokenHCT

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
%HCPARAMS.play_3_tokenHCT=0;
%HCPARAMS.play_3_tokenHCT=get(hObject,'Value');
HCPARAMS.play_3_tokenHCT=get(handles.play_3_tokenHCT,'Value');
if HCPARAMS.play_3_tokenHCT==1
    set(handles.Cover3,'Visible','Off')
    set(handles.Cover4,'Visible','Off')
    set(handles.Cover1,'Visible','Off')
    set(handles.Cover2,'Visible','Off')
elseif HCPARAMS.play_3_tokenHCT==0
    set(handles.Cover3,'Visible','On')
    set(handles.Cover4,'Visible','On')
    set(handles.Cover1,'Visible','Off')
    set(handles.Cover2,'Visible','Off')
end
set(handles.play_3_tokenHCT,'Value',HCPARAMS.play_3_tokenHCT);
%play_3_tokenHCT=HCPARAMS.play_3_tokenHCT;
curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)

function play_3_tokenHCT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numhc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% 
% %     
% function inter_tok_int_ms_Callback(hObject, eventdata, handles)
% % hObject    handle to text (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of text as text
% %        str2double(get(hObject,'String')) returns contents of text as a double
% 
% curdir=pwd;
% cd C:\EXPERIMENTS\CODES\PARAMS
% load HCPARAMS
% cd (curdir)
% % str_inter_tok_int_ms=get(hObject,'Value');HCPARAMS.text=str2num(str_inter_tok_int_ms);
% % 
%  %HCPARAMS.text=get(handles.text,'Value');str_inter_tok_int_ms=HCPARAMS.text;
% HCPARAMS.str_inter_tok_int_ms=str2num(get(hObject,'String')); %returns contents of num_tones as a double
% set(handles.text,'Value',HCPARAMS.text);
% 
% 
% curdir=pwd;
% cd C:\EXPERIMENTS\CODES\PARAMS
% save HCPARAMS HCPARAMS
% cd (curdir)
% % 
% % % --- Executes during object creation, after setting all properties.
% function inter_tok_int_ms_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to text (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


function freq_f0_Callback(hObject, eventdata, handles)
% hObject    handle to freq_f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq_f0 as text
%        str2double(get(hObject,'String')) returns contents of freq_f0 as a double

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
load HCPARAMS
cd (curdir)

%str_inter_tok_int_ms=get(hObject,'Value');HCPARAMS.text=play_3_tokenHCT;
HCPARAMS.freq_f0=get(handles.freq_f0,'Value');freq_f0=HCPARAMS.freq_f0;
%HCPARAMS.text=get(handles.text,'Value');str_inter_tok_int_ms=HCPARAMS.text;
%text=str2num(str_inter_tok_int_ms);

curdir=pwd;
cd C:\EXPERIMENTS\CODES\PARAMS
save HCPARAMS HCPARAMS
cd (curdir)


% --- Executes during object creation, after setting all properties.
function freq_f0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

