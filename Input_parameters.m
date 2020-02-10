function varargout = Input_parameters(varargin)
% INPUT_PARAMETERS MATLAB code for Input_parameters.fig
%      INPUT_PARAMETERS, by itself, creates a new INPUT_PARAMETERS or raises the existing
%      singleton*.
%
%      H = INPUT_PARAMETERS returns the handle to a new INPUT_PARAMETERS or the handle to
%      the existing singleton*.
%
%      INPUT_PARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INPUT_PARAMETERS.M with the given input arguments.
%
%      INPUT_PARAMETERS('Property','Value',...) creates a new INPUT_PARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Input_parameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Input_parameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Input_parameters

% Last Modified by GUIDE v2.5 16-Oct-2017 09:20:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Input_parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @Input_parameters_OutputFcn, ...
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


% --- Executes just before Input_parameters is made visible.
function Input_parameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Input_parameters (see VARARGIN)

% Choose default command line output for Input_parameters
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Input_parameters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Input_parameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
