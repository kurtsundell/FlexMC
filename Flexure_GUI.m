function varargout = Flexure_GUI(varargin)
% FLEXURE_GUI MATLAB code for Flexure_GUI.fig
%      FLEXURE_GUI, by itself, creates a new FLEXURE_GUI or raises the existing
%      singleton*.
%
%      H = FLEXURE_GUI returns the handle to a new FLEXURE_GUI or the handle to
%      the existing singleton*.
%
%      FLEXURE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLEXURE_GUI.M with the given input arguments.
%
%      FLEXURE_GUI('Property','Value',...) creates a new FLEXURE_GUI or raises the
%      existing singleton*.  Starting from the topo_rect, property value pairs are
%      applied to the GUI before Flexure_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Flexure_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Flexure_GUI

% Last Modified by GUIDE v2.5 12-Dec-2019 10:54:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Flexure_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Flexure_GUI_OutputFcn, ...
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

% --- Executes just before Flexure_GUI is made visible.
function Flexure_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Flexure_GUI (see VARARGIN)

% Choose default command line output for DZmix
handles.output = hObject;
% Choose default command line output for Flexure_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Flexure_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Flexure_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in run_model.
function run_model_Callback(hObject, eventdata, handles)
setting_up=msgbox('Setting up calculation');
%set variable "pin_m" to 0 to pin the back of the load, alternatively the
%front of the load will be pinned at "pin"
rad_on=get(handles.ui_pin_location,'selectedobject');
switch rad_on
	case handles.front
Pin_location = 'Front';
	pin_front=1;
	case handles.back
Pin_location = 'Back';
	pin_front=0;
end
handles
% set variable limits
trials = str2num(get(handles.model_trials,'String')); %number of model simulations

numblocks = str2num(get(handles.num_blocks,'String'));
rho_in = str2num(get(handles.infill_p,'String')); %define infill density (kg/m^3)
rho_m = str2num(get(handles.mantle_p,'String')); %define mantle density (kg/m^3)
g = str2num(get(handles.gravity,'String')); %gravity (m/s^2)
E = str2num(get(handles.youngs_mod,'String')); %Young's Modulus (Pa)
v = str2num(get(handles.Poissons,'String')); %Poisson's ratio
pin = str2num(get(handles.pin_m,'String')); %pin location



h_min = str2num(get(handles.h_min,'String')); %minimum height of central load block (m)
h_max = str2num(get(handles.h_max,'String')); %maximum height of central load block (m)
width_min = str2num(get(handles.width_min,'String')); %minimum width of central load block (m)
width_max = str2num(get(handles.width_max,'String')); %maximum height of central load block (m)
rho_c_min = str2num(get(handles.density_min,'String')); %minimum load density (kg/m^3)
rho_c_max = str2num(get(handles.density_max,'String')); %maximum load density (kg/m^3)
Te_min = str2num(get(handles.Te_min,'String')); %minimum effecive elastic thickness (m) 
Te_max = str2num(get(handles.Te_max,'String')); %maximum effecive elastic thickness (m) 
model_space=str2num(get(handles.model_space,'String')); %model space (m)
exp_max=str2num(get(handles.exp_max,'String'));
exp_min=str2num(get(handles.exp_min,'String'));
exponent_all=(exp_max-exp_min).*rand(trials,1)+exp_min;
max_angle = str2num(get(handles.max_angle,'String')); %max slope angle



%compare number of load blocks to minimum width
while width_min/numblocks<1000
    delete(setting_up);%delete setting up message box
    err_dlg=errordlg('Minimum width/number of load bloacks must be >= 1 km','Input Error');
    set(handles.num_blocks, 'BackgroundColor','Red');
    set(handles.width_min, 'BackgroundColor', 'Red');
    waitfor(err_dlg);
    numblocks = str2num(get(handles.num_blocks,'String'));
    width_min = str2num(get(handles.width_min,'String')); %minimum width of central load block (m)
end
set(handles.num_blocks, 'BackgroundColor','White');
set(handles.width_min, 'BackgroundColor', 'White');

%compare number of load blocks to minimum width
while numblocks>10
    delete(setting_up);%delete setting up message box
    err_dlg=errordlg('Maximum number of blocks is 10','Input Error');
    set(handles.num_blocks, 'BackgroundColor','Red');
    waitfor(err_dlg);
    numblocks = str2num(get(handles.num_blocks,'String'));
end
set(handles.num_blocks, 'BackgroundColor','White');


x_dim = transpose(0:1000:model_space); %define model space (m)



%%%%%%%%%%%%%%%% Get filter information %%%%%%%%%%%%%%%%%%%%%%%%%%
%Filter 1 Subsidence
if isempty(str2num(get(handles.filter1_amount,'String'))) == 1
    subsidence = NaN;
    set(handles.filter1_amount, 'String', NaN);
    else
    subsidence = str2num(get(handles.filter1_amount,'String')); %magnitude of tectonic subsidence used in filter 1 (m)
end

%Filter 2 Subsidence
if isempty(str2num(get(handles.filter2_amount,'String'))) == 1
    subsidence2 = NaN;
    set(handles.filter2_amount, 'String', NaN);
    else
    subsidence2 = str2num(get(handles.filter2_amount,'String')); %magnitude of tectonic subsidence used in filter 2 (m)
end

%Filter 3 Subsidence
if isempty(str2num(get(handles.filter3_amount,'String'))) == 1
    subsidence3 = NaN;
    set(handles.filter3_amount, 'String', NaN);
    else
    subsidence3 = str2num(get(handles.filter3_amount,'String')); %magnitude of tectonic subsidence used in filter 3 (m)
end

%Filter 4 Subsidence
if isempty(str2num(get(handles.filter4_amount,'String'))) == 1
    subsidence4 = NaN;
    set(handles.filter4_amount, 'String', NaN);
    else
    subsidence4 = str2num(get(handles.filter4_amount,'String')); %magnitude of tectonic subsidence used in filter 4 (m)
end

%Filter 5 Subsidence
if isempty(str2num(get(handles.filter5_amount,'String'))) == 1
    subsidence5 = NaN;
    set(handles.filter5_amount, 'String', NaN);
    else
    subsidence5 = str2num(get(handles.filter5_amount,'String')); %magnitude of tectonic subsidence used in filter 5 (m)
end

%Filter 1 Location
if isempty(str2num(get(handles.filter1_loc,'String'))) == 1
    subsidence_loc = NaN;
    set(handles.filter1_loc, 'String', NaN);
    else
    subsidence_loc = str2num(get(handles.filter1_loc,'String')); %subsidence at location 1 (m)
end

%Filter 2 Location
if isempty(str2num(get(handles.filter2_loc,'String'))) == 1
    subsidence_loc2 = NaN;
    set(handles.filter2_loc, 'String', NaN);
    else
    subsidence_loc2 = str2num(get(handles.filter2_loc,'String')); %subsidence at location 2 (m)
end

%Filter 3 Location
if isempty(str2num(get(handles.filter3_loc,'String'))) == 1
    subsidence_loc3 = NaN;
    set(handles.filter3_loc, 'String', NaN);
    else
    subsidence_loc3 = str2num(get(handles.filter3_loc,'String')); %subsidence at location 3 (m)
end

%Filter 4 Location
if isempty(str2num(get(handles.filter4_loc,'String'))) == 1
    subsidence_loc4 = NaN;
    set(handles.filter4_loc, 'String', NaN);
    else
    subsidence_loc4 = str2num(get(handles.filter4_loc,'String')); %subsidence at location 4 (m)
end

%Filter 5 Location
if isempty(str2num(get(handles.filter5_loc,'String'))) == 1
    subsidence_loc5 = NaN;
    set(handles.filter5_loc, 'String', NaN);
    else
    subsidence_loc5 = str2num(get(handles.filter5_loc,'String')); %subsidence at location 5 (m)
end

%Filter 1 Uncertainty
if isempty(str2num(get(handles.filter1_loc_unc,'String'))) == 1
    subuncertainty1 = NaN;
    set(handles.filter1_loc_unc, 'String', NaN);
    else
    subuncertainty1 = str2num(get(handles.filter1_loc_unc,'String')); %uncertainty in position of subsidence location 1 (m)
end

%Filter 2 Uncertainty
if isempty(str2num(get(handles.filter2_loc_unc,'String'))) == 1
    subuncertainty2 = NaN;
    set(handles.filter2_loc_unc, 'String', NaN);
    else
    subuncertainty2 = str2num(get(handles.filter2_loc_unc,'String')); %uncertainty in position of subsidence location 2 (m)
end

%Filter 3 Uncertainty
if isempty(str2num(get(handles.filter3_loc_unc,'String'))) == 1
    subuncertainty3 = NaN;
    set(handles.filter3_loc_unc, 'String', NaN);
    else
    subuncertainty3 = str2num(get(handles.filter3_loc_unc,'String')); %uncertainty in position of subsidence location 3 (m)
end

%Filter 4 Uncertainty
if isempty(str2num(get(handles.filter4_loc_unc,'String'))) == 1
    subuncertainty4 = NaN;
    set(handles.filter4_loc_unc, 'String', NaN);
    else
    subuncertainty4 = str2num(get(handles.filter4_loc_unc,'String')); %uncertainty in position of subsidence location 4 (m)
end

%Filter 5 Uncertainty
if isempty(str2num(get(handles.filter5_loc_unc,'String'))) == 1
    subuncertainty5 = NaN;
    set(handles.filter5_loc_unc, 'String', NaN);
    else
    subuncertainty5 = str2num(get(handles.filter5_loc_unc,'String')); %uncertainty in position of subsidence location 5 (m)
end
%%%END SET FILTER PARAMETERS

min_filter = subsidence_loc - subuncertainty1 + pin; %minimum distance used in filter 1 (m)
max_filter = subsidence_loc + subuncertainty1 + pin; %maximum distance used in filter 2 (m)
min_filter2 = subsidence_loc2 - subuncertainty2 + pin; %minimum distance used in filter 3 (m)
max_filter2 = subsidence_loc2 + subuncertainty2 + pin; %maximum distance used in filter 4 (m)
min_filter3 = subsidence_loc3 - subuncertainty3 + pin; %minimum distance used in filter 5 (m)
max_filter3 = subsidence_loc3 + subuncertainty3 + pin; %maximum distance used in filter 1 (m)
min_filter4 = subsidence_loc4 - subuncertainty4 + pin; %minimum distance used in filter 2 (m)
max_filter4 = subsidence_loc4 + subuncertainty4 + pin; %maximum distance used in filter 3 (m)
min_filter5 = subsidence_loc5 - subuncertainty5 + pin; %minimum distance used in filter 4 (m)
max_filter5 = subsidence_loc5 + subuncertainty5 + pin; %maximum distance used in filter 5 (m)

rho_c_rand = rand(trials,1); %random number array for load density (between 0 and 1)
rho_c = rho_c_min + (rho_c_max-rho_c_min).*rho_c_rand; %generate random values for load density (kg/m^3)

delta_rho = rho_m - rho_in; %load to mantle density contrast (kg/m^3)



%Calculate Te gradient for each trial
complete=0
while complete==0
rad_EET=get(handles.EET_gradient,'selectedobject');
switch rad_EET
    case handles.EET_constant %Constant Te=Te_low
        Alpha=zeros(size(x_dim,1),trials);
        Te_rand = rand(trials,1); %random number array for effective elastic thickness  (between 0 and 1)
        Te_low = Te_min + (Te_max-Te_min).*Te_rand; %generate random values for effective elastic thickness (m)
        D = (E*Te_low.*Te_low.*Te_low)/(12*(1-v.*v)); %calculate flexural rigidity (D) from random values above (Nm) (Turcotte and Schubert, 2014)
        Alpha_left=transpose((4*D./(delta_rho.*g)).^.25); %calculate flexural parameter (alpha) from random values above (m) (Turcotte and Schubert, 2014)
        Alpha=repmat(Alpha_left, length(x_dim), 1);
        complete=1
        EET_grad='Constant';
        exp_max=NaN([1 1]);
        exp_min=NaN([1 1]);
        
        clear Te_rand Te_low;
    case handles.EET_linear
        Te_min_rand = 0.5*rand(trials,1); %random number array for MINIMUM effective elastic thickness  (between 0 and 0.5)
        Te_max_rand = (ones(size(Te_min_rand))-Te_min_rand); %random number array for MAXIMUM effective elastic thickness  (complement of Te_min_rand)
        Te_low = Te_min + (Te_max-Te_min).*Te_min_rand; %generate random values for MINIMUM effective elastic thickness (m)
        Te_high= Te_min + (Te_max-Te_min).*Te_max_rand; %generate random values for MAXIMUM effective elastic thickness (m)

        %Calculate min and max D and alpha
        D_min = (E*Te_low.*Te_low.*Te_low)/(12*(1-v.*v)); %calculate flexural rigidity (D) from random values above (Nm) (Turcotte and Schubert, 2014)
        Alpha_left =transpose((4*D_min./(delta_rho.*g)).^.25); %calculate flexural parameter (alpha) from random values above (m) (Turcotte and Schubert, 2014)
        D_max = (E*Te_high.*Te_high.*Te_high)/(12*(1-v.*v)); %calculate flexural rigidity (D) from random values above (Nm) (Turcotte and Schubert, 2014)
        Alpha_right =transpose((4*D_max./(delta_rho.*g)).^.25); %calculate flexural parameter (alpha) from random values above (m) (Turcotte and Schubert, 2014)
        %randomize a right or left thinning EET
        randomize=rand(trials,1);
        for j=1:trials
            if randomize(j)>0.5
                alpha_l(j)=Alpha_left(j);
                alpha_r(j)=Alpha_right(j);
            else
                alpha_r(j)=Alpha_left(j);
                alpha_l(j)=Alpha_right(j);
            end
        end
        %end randomize
        
        alpha_step=(alpha_l-alpha_r)/(size(x_dim,1)-1);
        Alpha=zeros(size(x_dim,1),trials);
        for i=1:trials
           if alpha_l(i)>alpha_r(i)
              alpha_step(i)=-abs(alpha_step(i));
           else
               alpha_step(i)=abs(alpha_step(i));
           end
        end
        %Define Alpha for each x step.  
        for i=1:trials
                if alpha_step(i)==0
                   Alpha(1:length(x_dim),i)=alpha_l(i);
                else
                   Alpha(:,i)=transpose(alpha_l(i):alpha_step(i):alpha_r(i));
                end
        end
        Alpha_left=alpha_l;
        Alpha_right=alpha_r;
    
        complete=1;
        EET_grad='Linear';
        exp_max=NaN([1 1]);
        exp_min=NaN([1 1]);
        clear Te_rand_min Te_max_rand rho_c_rand Te_high Te_low Alpha_left Alpha_right alpha_l alpha_r; 
    case handles.EET_exp
        Te_min_rand = 0.5*rand(trials,1); %random number array for MINIMUM effective elastic thickness  (between 0 and 0.5)
        Te_max_rand = (ones(size(Te_min_rand))-Te_min_rand); %random number array for MAXIMUM effective elastic thickness  (complement of Te_min_rand)
        Te_low = Te_min + (Te_max-Te_min).*Te_min_rand; %generate random values for MINIMUM effective elastic thickness (m)
        Te_high= Te_min + (Te_max-Te_min).*Te_max_rand; %generate random values for MAXIMUM effective elastic thickness (m)

        %Calculate min and max D and alpha
        D_min = (E*Te_low.*Te_low.*Te_low)/(12*(1-v.*v)); %calculate flexural rigidity (D) from random values above (Nm) (Turcotte and Schubert, 2014)
        Alpha_left =transpose((4*D_min./(delta_rho.*g)).^.25); %calculate flexural parameter (alpha) from random values above (m) (Turcotte and Schubert, 2014)
        D_max = (E*Te_high.*Te_high.*Te_high)/(12*(1-v.*v)); %calculate flexural rigidity (D) from random values above (Nm) (Turcotte and Schubert, 2014)
        Alpha_right =transpose((4*D_max./(delta_rho.*g)).^.25); %calculate flexural parameter (alpha) from random values above (m) (Turcotte and Schubert, 2014)
        
        for i=1:trials
            exponent_i=exponent_all(i);
            if exponent_i>=0
            for j=1:length(x_dim)
                Alpha(j,i)=Alpha_left(i)+(Alpha_right(i)-Alpha_left(i))*j^exponent_i/max(length(x_dim))^exponent_i;
            end
            else
                for j=1:length(x_dim)
                Alpha(j,i)=Alpha_right(i)*((Alpha_right(i)-(Alpha_right(i)-Alpha_left(i))*j^exponent_i)/Alpha_right(i))/...
                    ((Alpha_right(i)-(Alpha_right(i)-Alpha_left(i))*max(length(x_dim))^exponent_i)/Alpha_right(i));
                end
            end
        end
        complete=1;
        EET_grad='Exponential';
        clear Te_rand_min Te_max_rand rho_c_rand Te_high Te_low Alpha_left Alpha_right alpha_l alpha_r;
    case handles.EET_stepped
        Te_min_rand = 0.5*rand(trials,1); %random number array for MINIMUM effective elastic thickness  (between 0 and 0.5)
        Te_max_rand = (ones(size(Te_min_rand))-Te_min_rand); %random number array for MAXIMUM effective elastic thickness  (complement of Te_min_rand)
        Te_low = Te_min + (Te_max-Te_min).*Te_min_rand; %generate random values for MINIMUM effective elastic thickness (m)
        Te_high= Te_min + (Te_max-Te_min).*Te_max_rand; %generate random values for MAXIMUM effective elastic thickness (m)

        %Calculate min and max D and alpha
        D_min = (E*Te_low.*Te_low.*Te_low)/(12*(1-v.*v)); %calculate flexural rigidity (D) from random values above (Nm) (Turcotte and Schubert, 2014)
        Alpha_left =transpose((4*D_min./(delta_rho.*g)).^.25); %calculate flexural parameter (alpha) from random values above (m) (Turcotte and Schubert, 2014)
        D_max = (E*Te_high.*Te_high.*Te_high)/(12*(1-v.*v)); %calculate flexural rigidity (D) from random values above (Nm) (Turcotte and Schubert, 2014)
        Alpha_right =transpose((4*D_max./(delta_rho.*g)).^.25); %calculate flexural parameter (alpha) from random values above (m) (Turcotte and Schubert, 2014)
        
        step_number=str2num(get(handles.step_number,'String'))%get number of steps
        while step_number>3 && get(handles.EET_stepped,'Value')==1;
            err_dlg=errordlg('Maximum number of steps is 3','Input Error');
            set(handles.step_number, 'BackgroundColor','Red');
            waitfor(err_dlg);
            step_number=str2num(get(handles.step_number,'String'))%check number of steps
            rad_EET=get(handles.EET_gradient,'selectedobject');
        end
        set(handles.step_number, 'BackgroundColor','White');
        if step_number == 1
            for i=1:trials
            Alpha(1:round(length(x_dim)/2),i)=Alpha_left(i);
            Alpha(1+round(length(x_dim)/2):length(x_dim),i)=Alpha_right(i);
            end
            complete=1;
        elseif step_number == 2
            for i=1:trials
            Alpha(1:round(length(x_dim)/3),i)=Alpha_left(i);
            Alpha(round(length(x_dim)/3)+1:2*round(length(x_dim)/3),i)=(Alpha_right(i)+Alpha_left(i))/2;
            Alpha(2*round(length(x_dim)/3)+1:length(x_dim),i)=Alpha_right(i);
            end
            complete=1;
        elseif step_number == 3
            for i=1:trials
            Alpha(1:round(length(x_dim)/4),i)=Alpha_left(i);
            Alpha(round(length(x_dim)/4)+1:2*round(length(x_dim)/4),i)=Alpha_left(i)+(Alpha_right(i)-Alpha_left(i))/3;
            Alpha(2*round(length(x_dim)/4)+1:3*round(length(x_dim)/4),i)=Alpha_left(i)+2*(Alpha_right(i)-Alpha_left(i))/3;
            Alpha(3*round(length(x_dim)/4)+1:length(x_dim),i)=Alpha_right(i);
            end
            complete=1;
        else
            complete=0;
        end
        EET_grad='Stepped';
        exp_max=step_number;
        exp_min=NaN([1 1]);
        clear Te_rand_min Te_max_rand rho_c_rand Te_high Te_low Alpha_left Alpha_right alpha_l alpha_r;
    otherwise
        %
    end
end
%END build Alpha

h_rand = rand(trials,1); %random number array for load height (between 0 and 1)
h= h_min + (h_max - h_min)*h_rand; %generate random values for load height (m)

test_plateau_width = 0;
while test_plateau_width < 1
    width_rand = rand(trials,1); %random number array for central load block width (between 0 and 1)
    width_total = (width_min + (width_max - width_min)*width_rand); %generate random values for central load block width (m)
    
    rad_on=get(handles.ui_topography,'selectedobject');
    switch rad_on
        case handles.topo_rect
            Topo = 1;
        case handles.topo_right
            Topo = 2;
        case handles.topo_left
            Topo = 3;
    end
    if Topo >1;
            
            low_angle=0;
            for i=1:size(width_total)
                count=1;
                while max_angle<atand(h(i)/width_total(i))
                    ma(i)=atand(h(i)/width_total(i));
                width_total(i)=(width_min + (width_max - width_min)*rand());
                h(i)= h_min + (h_max - h_min)*rand();
                angle_catch=i;
                count=count+1
                if count>100
                    break
                    low_angle=1;
                end
                if low_angle==1
                break
                end
                end
            end
    end
    
    if (pin_front==1) && (max(width_total) > pin)
        test_plateau_width = 0;
        width_max = pin;
      else
        test_plateau_width = 2;
    end
end
block_width = round(width_total/numblocks);

for i=1:trials
    x_step=0;
    
    for j=1:numblocks
        if pin_front==0 %pin load in back
            xc(i,j) = x_step+(block_width(i,1)/2+pin); %determine xc of central load block for each simulation
            x_step=x_step+block_width(i,1);
        else %pin load in front
            
            xc(i,j) = pin + x_step-(block_width(i,1)/2);
            x_step=x_step-block_width(i,1);
        end
        
    end
end
loadwidth_all = repmat(block_width,1,numblocks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Build Topography                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%     Rectangular load geometry       %%%%%%%%%%%%%%%%%%%
rad_on=get(handles.ui_topography,'selectedobject');
switch rad_on
case handles.topo_rect
Topo_type = 'Rectangular';
    % build matrix with all xc locations
    xc_all = zeros(trials,numblocks); %preallocate
    for i = 1:trials
        angle(i)=nan;
        xc_all(i,:) = [xc(i,:)];
    end

    % build matrix with all load heights
    h_all = zeros(numblocks,trials); %preallocate
    for i = 1:trials
        h_all(:,i) = [h(i,1)];
    end
end
 
%%%%%%%%%%%%%%%     if pinned in back = Thrust belt load geometry       %%%%%%%%%%%%%%%%%%%
rad_on=get(handles.ui_topography,'selectedobject');
switch rad_on
case handles.topo_right
Topo_type = 'Right';
    % build matrix with all xc locations
    xc_all = zeros(trials,numblocks); %preallocate
    for i = 1:trials
        xc_all(i,:) = [xc(i,:)];
    end

    % build matrix with all load heights
    h_all = zeros(numblocks,trials); %preallocate
    for i = 1:trials
        angle(i)=atand(h(i,1)/width_total(i));
        for j=1:numblocks
        h_all(j,i) = tand(angle(i))*(width_total(i)-(j*block_width(i)-block_width(i)));
        end
    end
end

%%%%%%%%%%%%%%%     if pinned in back = Laramide load geometry       %%%%%%%%%%%%%%%%%%%
rad_on=get(handles.ui_topography,'selectedobject');
switch rad_on
case handles.topo_left
Topo_type = 'Left';
    % build matrix with all xc locations
    xc_all = zeros(trials,numblocks); %preallocate
    for i = 1:trials
        xc_all(i,:) = [xc(i,:)];
    end
    
    % build matrix with all load heights
    h_all = zeros(numblocks,trials); %preallocate
    for i = 1:trials
        angle(i)=atand(h(i,1)/width_total(i));
        for j=1:numblocks
        h_all(j,i) = tand(angle(i))*(width_total(i)-((numblocks-j)*block_width(i)));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate individual flexural profiles (w) and dimensionalize       %
%           (w_dim) for each load block                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_dim = zeros(length(x_dim), numblocks,  trials); %preallocate

perc_done = 0;
%Broken beam model from Hetenyi (1979)
rad_on=get(handles.ui_model_type,'selectedobject');

x = zeros(length(x_dim), numblocks,  trials); %preallocate
w = zeros(length(x_dim), numblocks,  trials); %preallocate
w_dim = zeros(length(x_dim), numblocks,  trials); %preallocate
delete(setting_up);%delete setting up message box
f=waitbar(j/trials, "Calculating",'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');%calculation waitbar
setappdata(f,'canceling',0);

switch rad_on
	case handles.broken_plate
Plate_type = 'Broken';
    for i = 1:numblocks
        for j = 1:trials
            if getappdata(f,'canceling')
                break
            end
            waitbar(((i-1)*trials+j)/(trials*numblocks));
            [x(:,i,j),w(:,i,j),w_dim(:,i,j)] = Hetenyi_variable_EET(j,h_all(i,j), xc_all(j,i), ...
                loadwidth_all(j,i), Alpha(:,j), rho_c(j,1), delta_rho, x_dim, g);
%			[x(:,i,j), w(:,i,j), w_dim(:,i,j)] = Hetenyi2(j,h_all(i,j), xc_all(j,i), loadwidth_all(j,i), alpha(j,1), rho_c(j,1), delta_rho);
        end
        if getappdata(f,'canceling')
        break
        end
    end
    delete(f);
end

%Infinite beam model from Wangen (2010)
rad_on=get(handles.ui_model_type,'selectedobject');
switch rad_on
	case handles.infinite_plate
Plate_type = 'Infinite';
    for i = 1:numblocks
        for j = 1:trials
            if getappdata(f,'canceling')
                break
            end
            waitbar(((i-1)*trials+j)/(trials*numblocks));
            perc_done = perc_done+1/(numblocks)/(trials/100);
            [w_dim(:,i,j)] = Wangen2010_variable_EET(h_all(i,j), xc_all(j,i), loadwidth_all(j,i), Alpha(:,j), rho_c(j,1), delta_rho, x_dim, g);
        end
        if getappdata(f,'canceling')
        break
        end
    end
    delete(f);
end

w_dim_sum = sum(w_dim,2); %sum all individual flexural profiles

% final dimensionalized flexural profile output
w_dim_out = zeros(length(x_dim), trials);

for j = 1:trials
    w_dim_out(:,j) = w_dim_sum(:,:,j);
end

%{
clear width_total
clear delta_rho
clear h
clear h_all
clear halfwidth_total
clear halfwidth_total
clear loadwidth_all
clear w
clear w_dim_sum
clear x
clear xc_all
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       filter model results                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xc = round(xc); %round each xc to nearest integer

w_dim_out_trim = w_dim_out; %preallocate

% find location along each flexural profile that matches subsidence filter

        
    xfilt_loc = zeros(trials,1); %preallocate 
    for i=1:trials
        if(isnan(subsidence))==0
            
            wfilt_idx = min(find(w_dim_out_trim(:,i)>-subsidence));
            try
            xfilt_loc(i,1) = x_dim(wfilt_idx,1);
            catch
                err_dlg=errordlg('Input subsidence not found within model space.  Try changing Subsidence Filter 1 or increasing model space','Input Error');
                waitfor(err_dlg);
                error;
            end
        else
            xfilt_loc(i,1)=nan;
        end
    end

% find location along each flexural profile that matches subsidence filter
% 2

xfilt_loc2 = zeros(trials,1); %preallocate 
for i=1:trials
    if(isnan(subsidence2))==0
        wfilt_idx2 = min(find(w_dim_out_trim(:,i)>-subsidence2));
           try 
            xfilt_loc2(i,1) = x_dim(wfilt_idx2,1);
            catch
                err_dlg=errordlg('Input subsidence not found within model space.  Try changing Subsidence Filter 2 or increasing model space','Input Error');
                waitfor(err_dlg);
                error;
            end
    else
        xfilt_loc2(i,1)=nan;
end
end

% find location along each flexural profile that matches subsidence filter
% 3

xfilt_loc3 = zeros(trials,1); %preallocate 
for i=1:trials
    if(isnan(subsidence3))==0
		wfilt_idx3 = min(find(w_dim_out_trim(:,i)>-subsidence3));
        try
		xfilt_loc3(i,1) = x_dim(wfilt_idx3,1);
        catch
                err_dlg=errordlg('Input subsidence not found within model space.  Try changing Subsidence Filter 3 or increasing model space','Input Error');
                waitfor(err_dlg);
                error;
            end
    else
        xfilt_loc3(i,1)=nan;
end
end

% find location along each flexural profile that matches subsidence filter
% 4

xfilt_loc4 = zeros(trials,1); %preallocate 
for i=1:trials
    if(isnan(subsidence4))==0
		wfilt_idx4 = min(find(w_dim_out_trim(:,i)>-subsidence4));
        try
		xfilt_loc4(i,1) = x_dim(wfilt_idx4,1);
        catch
                err_dlg=errordlg('Input subsidence not found within model space.  Try changing Subsidence Filter 4 or increasing model space','Input Error');
                waitfor(err_dlg);
                error;
            end
    else
        xfilt_loc4(i,1)=nan;
end
end

% find location along each flexural profile that matches subsidence filter
% 5

xfilt_loc5 = zeros(trials,1); %preallocate 
for i=1:trials
    if(isnan(subsidence5))==0
		wfilt_idx5 = min(find(w_dim_out_trim(:,i)>-subsidence5));
        try
		xfilt_loc5(i,1) = x_dim(wfilt_idx5,1);
        catch
                err_dlg=errordlg('Input subsidence not found within model space.  Try changing Subsidence Filter 5 or increasing model space','Input Error');
                waitfor(err_dlg);
                error;
        end
    else
        xfilt_loc5(i,1)=nan;
end
end


Te_low=transpose((3*Alpha(1,:).^4*delta_rho*g*(1-v^2)/E).^(1/3));
Te_high=transpose((3*Alpha(size(w_dim_out,1),:).^4*delta_rho*g*(1-v^2)/E).^(1/3));
D_min = (E*Te_low.*Te_low.*Te_low)/(12*(1-v.*v));
D_max = (E*Te_high.*Te_high.*Te_high)/(12*(1-v.*v));
Alpha_left = Alpha(1,:);
Alpha_right = Alpha(size(w_dim_out,1),:);

% concatenate and sort results to filter based on minimum and maximum filter distances
x_filter = horzcat(...
    xfilt_loc, ...
    xfilt_loc2, ...
    xfilt_loc3, ...
    xfilt_loc4, ...
    xfilt_loc5, ...
    transpose(w_dim_out), ...
    h, ...
    width_total, ...
    rho_c, ...
    Te_low, ...
    Te_high, ...
    D_min, ...
    D_max,...
    exponent_all,...
    transpose(Alpha_left), ...
    transpose(Alpha_right), ...
    transpose(mean(Alpha)),...
    xc, ...
    transpose(h_all), ...
    transpose(angle), ...
    loadwidth_all,...
    transpose(Alpha));
x_filter_sort=sortrows(x_filter,1);

% filter results based on minimum and maximum filter distances and replace unsuccessful model fits with zeros
uncert_orientation=get(handles.uncert_orientation,'selectedobject');
switch uncert_orientation
	case handles.Uncert_horz
        Uncertainty_orientation = 'Horizontal';
if(isnan(subsidence))==0
for i = 1:trials
if x_filter_sort(i,1) >= min_filter && x_filter_sort(i,1) <= max_filter
    x_filter_sort(i,:) = x_filter_sort(i,:);
else
    x_filter_sort(i,:) = 0;
end
end
x_filter_sort( ~any(x_filter_sort,2), : ) = []; %remove all zeros (unsuccessful model fits)
end

if(isnan(subsidence2))==0
for i = 1:length(x_filter_sort(:,1))
if x_filter_sort(i,2) >= min_filter2 && x_filter_sort(i,2) <= max_filter2
    x_filter_sort(i,:) = x_filter_sort(i,:);
else
    x_filter_sort(i,:) = 0;
end
end
x_filter_sort( ~any(x_filter_sort,2), : ) = []; %remove all zeros (unsuccessful model fits)
end

if(isnan(subsidence3))==0
for i = 1:length(x_filter_sort(:,1))
if x_filter_sort(i,3) >= min_filter3 && x_filter_sort(i,3) <= max_filter3
    x_filter_sort(i,:) = x_filter_sort(i,:);
else
    x_filter_sort(i,:) = 0;
end
end
x_filter_sort( ~any(x_filter_sort,2), : ) = []; %remove all zeros (unsuccessful model fits)
end

if(isnan(subsidence4))==0
for i = 1:length(x_filter_sort(:,1))
if x_filter_sort(i,4) >= min_filter4 && x_filter_sort(i,4) <= max_filter4
    x_filter_sort(i,:) = x_filter_sort(i,:);
else
    x_filter_sort(i,:) = 0;
end
end
x_filter_presort3=x_filter_sort;
x_filter_sort( ~any(x_filter_sort,2), : ) = []; %remove all zeros (unsuccessful model fits)
end

if(isnan(subsidence5))==0
for i = 1:length(x_filter_sort(:,1))
if x_filter_sort(i,5) >= min_filter5 && x_filter_sort(i,5) <= max_filter5
    x_filter_sort(i,:) = x_filter_sort(i,:);
else
    x_filter_sort(i,:) = 0;
end
end
x_filter_sort( ~any(x_filter_sort,2), : ) = []; %remove all zeros (unsuccessful model fits)
end
end
%END filter based on horizontal uncertainties

%OR filter results based on vertical uncertainties
switch uncert_orientation
	case handles.Uncert_vert
    Uncertainty_orientation = 'Vertical';
    %Vertical Filter 1    
    %xfilt_loc = zeros(trials,1); %preallocate 
    if(isnan(subsidence))==0    
    xfilt_subsidence_loc=round(subsidence_loc + pin,-3); %round subsidence to nearest km
        for i=1:length(x_filter_sort(:,1))
            [r_ind, c_ind]=find(x_dim==xfilt_subsidence_loc);
            if -x_filter_sort(i, r_ind+5)<=subsidence+subuncertainty1 && -x_filter_sort(i, r_ind+5)>=subsidence-subuncertainty1
                x_filter_sort(i,:) = x_filter_sort(i,:);
            else
                x_filter_sort(i,:) = 0;
            end
        end
    else
        xfilt_loc(i,1)=nan;
       
    end
    x_filter_sort( ~any(x_filter_sort,2), : ) = []; %remove all zeros (unsuccessful model fits)
    
    %Vertical Filter 2    
    %xfilt_loc = zeros(trials,1); %preallocate 
    if(isnan(subsidence2))==0
    xfilt_subsidence_loc2=round(subsidence_loc2 + pin,-3); %round subsidence to nearest km
       for i=1:length(x_filter_sort(:,1))
            [r_ind, c_ind]=find(x_dim==xfilt_subsidence_loc2);
            if -x_filter_sort(i, r_ind+5)<=subsidence2+subuncertainty2 && -x_filter_sort(i, r_ind+5)>=subsidence2-subuncertainty2
                x_filter_sort(i,:) = x_filter_sort(i,:);
            else
                x_filter_sort(i,:) = 0;
            end
       end
    else
        xfilt_loc(i,1)=nan;
        xfilt_subsidence_loc2=NaN;
    end
    x_filter_sort( ~any(x_filter_sort,2), : ) = []; %remove all zeros (unsuccessful model fits)
    
    %Vertical Filter 3    
    %xfilt_loc = zeros(trials,1); %preallocate 
    if(isnan(subsidence3))==0
    xfilt_subsidence_loc3=round(subsidence_loc3 + pin,-3); %round subsidence to nearest km
       for i=1:length(x_filter_sort(:,1))
            [r_ind, c_ind]=find(x_dim==xfilt_subsidence_loc3);
            if -x_filter_sort(i, r_ind+5)<=subsidence3+subuncertainty3 && -x_filter_sort(i, r_ind+5)>=subsidence3-subuncertainty3
                x_filter_sort(i,:) = x_filter_sort(i,:);
            else
                x_filter_sort(i,:) = 0;
            end
       end
   else
       xfilt_loc(i,1)=nan;
        xfilt_subsidence_loc3=NaN;
    end
    x_filter_sort( ~any(x_filter_sort,2), : ) = []; %remove all zeros (unsuccessful model fits)
    
    %Vertical Filter 4    
    %xfilt_loc = zeros(trials,1); %preallocate 
    if(isnan(subsidence4))==0
    xfilt_subsidence_loc4=round(subsidence_loc4 + pin,-3); %round subsidence to nearest km
       for i=1:length(x_filter_sort(:,1))
            [r_ind, c_ind]=find(x_dim==xfilt_subsidence_loc4);
            if -x_filter_sort(i, r_ind+5)<=subsidence4+subuncertainty4 && -x_filter_sort(i, r_ind+5)>=subsidence4-subuncertainty4
                x_filter_sort(i,:) = x_filter_sort(i,:);
            else
                x_filter_sort(i,:) = 0;
            end
       end
    else
       xfilt_loc(i,1)=nan;
        xfilt_subsidence_loc4=NaN;
    end
    x_filter_sort( ~any(x_filter_sort,2), : ) = []; %remove all zeros (unsuccessful model fits)
    
    %Vertical Filter 5    
    %xfilt_loc = zeros(trials,1); %preallocate 
    if(isnan(subsidence5))==0
        xfilt_subsidence_loc5=round(subsidence_loc5 + pin,-3); %round subsidence to nearest km
        for i=1:length(x_filter_sort(:,1))
            [r_ind, c_ind]=find(x_dim==xfilt_subsidence_loc5);
            if -x_filter_sort(i, r_ind+5)<=subsidence5+subuncertainty5 && -x_filter_sort(i, r_ind+5)>=subsidence5-subuncertainty5
                x_filter_sort(i,:) = x_filter_sort(i,:);
            else
                x_filter_sort(i,:) = 0;
            end
        end
    else
        xfilt_loc(i,1)=nan;
        xfilt_subsidence_loc5=NaN;
    end
    x_filter_sort( ~any(x_filter_sort,2), : ) = []; %remove all zeros (unsuccessful model fits)
end
%END filter based on vertical uncertainties

if size(x_filter_sort,1) == 0
	err_dlg=errordlg('No model fits... Try again.','Nope.');
    waitfor(err_dlg);
    else
        %
end


% create individual arrays of filtered model results
i=1;
x_filter_sort_header(1,i)={'Subsidence filter locations'};
filtered_dist_from_origin = x_filter_sort(:,i:i+numblocks-1);%populate locations of subsidence filters 

i=i+5;
x_filter_sort_header(1,i)={'Subsidence Profile'};
filtered_w_dim_out = x_filter_sort(:,i:length(x_dim)+i-1);

i=i+length(x_dim);
filtered_h = x_filter_sort(:,i);
x_filter_sort_header(1,i)={'Load height'};

i=i+1;
filtered_width_total = x_filter_sort(:,i);
x_filter_sort_header(1,i)={'Total Load Width'};

i=i+1;
filtered_rho_c = x_filter_sort(:,i);
x_filter_sort_header(1,i)={'Crustal density'};

i=i+1;
filtered_Te = x_filter_sort(:,i);
x_filter_sort_header(1,i)={'EET'};


i=i+1;
filtered_Te_max=x_filter_sort(:,i);  % holding place for Te_max
x_filter_sort_header(1,i)={'EET max'};


i=i+1;
filtered_D = x_filter_sort(:,i);
x_filter_sort_header(1,i)={'D'};


i=i+1;
filtered_D_max= x_filter_sort(:,i); % holding place for D_max
x_filter_sort_header(1,i)={'D max'};


i=i+1;
filtered_exponent = x_filter_sort(:,i);
x_filter_sort_header(1,i)={'Exponent'};

i=i+1;
filtered_alpha_left = x_filter_sort(:,i);
x_filter_sort_header(1,i)={'Alpha left'};

i=i+1;
filtered_alpha_right= x_filter_sort(:,i);
x_filter_sort_header(1,i)={'Alpha right'};

i=i+1;
filtered_mean_alpha=x_filter_sort(:,i);
x_filter_sort_header(1,i)={'Alpha mean'};

i=i+1;
filtered_xc = x_filter_sort(:,i:i+numblocks-1);
x_filter_sort_header(1,i)={'Load centers'};

i=i+numblocks;
filtered_h_all = x_filter_sort(:,i:i+numblocks-1);
x_filter_sort_header(1,i)={'Load heights'};

i=i+numblocks;
filtered_angle = x_filter_sort(:,i);
x_filter_sort_header(1,i)={'Angle'};

i=i+1
filtered_loadwidth_all = x_filter_sort(:,i:i+(numblocks-1));
x_filter_sort_header(1,i)={'Load widths'};

i=i+numblocks;
filtered_Alpha_x=x_filter_sort(:,i:end);
x_filter_sort_header(1,i)={'Alpha'};

filtered_block_width = filtered_loadwidth_all(:,1);

x_filter_sort_with_headers=cell(size(x_filter_sort,1)+1, size(x_filter_sort, 2));


x_filter_sort_with_headers(2:end, :)=num2cell(x_filter_sort);
%{
switch rad_EET
    case handles.EET_constant
        %use filtered_Te from above
        filtered_Te_left=[];
        filtered_Te_right=[];
    otherwise
        for i=1:length(filtered_Te)
            
                filtered_Te(i)=(3*filtered_alpha_mean(i)^4*delta_rho*g*(1-v^2)/E)^(1/3);
                filtered_Te_left(i)=(3*filtered_alpha_left(i)^4*delta_rho*g*(1-v^2)/E)^(1/3);
                filtered_Te_right(i)=(3*filtered_alpha_right(i)^4*delta_rho*g*(1-v^2)/E)^(1/3);
        end
end
%}

%calculate means and standard deviations for individual arrays of filtered model results
filtered_h_mean = mean(filtered_h);
filtered_h_std = std(filtered_h);
filtered_h_all_mean = mean(filtered_h_all,2);
filtered_h_all_std = std(filtered_h_all,0,2);
filtered_h_all_max = max(filtered_h_all,[],2);
filtered_h_all_min = min(filtered_h_all,[],2);
filtered_width_total_mean = mean(filtered_width_total);
filtered_width_total_std = std(filtered_width_total);
filtered_block_width_mean = mean(filtered_loadwidth_all(:,1));
filtered_block_width_std = std(filtered_loadwidth_all(:,1));
filtered_rho_c_mean = mean(filtered_rho_c);
filtered_rho_c_std = std(filtered_rho_c);
filtered_Te_mean = (3*mean(filtered_mean_alpha).^4*delta_rho*g*(1-v^2)/E).^(1/3);
    filtered_Te_variation_between_profiles=std((3*(filtered_mean_alpha).^4*delta_rho*g*(1-v^2)/E).^(1/3));
    filtered_Te_variation_within_profiles=rssq((3*(std(filtered_Alpha_x,0,2)).^4*delta_rho*g*(1-v^2)/E).^(1/3));
    total_variation=[filtered_Te_variation_between_profiles filtered_Te_variation_within_profiles];
filtered_Te_std = rssq(total_variation);
filtered_D_mean = mean(filtered_D);
filtered_D_std = std(filtered_D);
filtered_alpha_mean = mean(filtered_mean_alpha);
filtered_alpha_std = std(filtered_mean_alpha);
filtered_angle_mean = mean(filtered_angle);
filtered_angle_std = std(filtered_angle);

filtered_alpha_mean2=mean(mean(Alpha));

if contains(EET_grad,'Exp') || contains(EET_grad,'exp') 
    filtered_exponent_mean = mean(filtered_exponent);
    filtered_exponent_std = std(filtered_exponent);
else
    filtered_exponent_mean = 'NaN';
    filtered_exponent_std = 'NaN';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       sum filtered model result flexural profiles and loads          %
%                   to determine topographic height                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sum_plateau = filtered_w_dim_out; %preallocate
%calculate depth on flexed profile (w) beneath the center of...
%each block (xc) 
for j=1:size(x_filter_sort,1)
    x_step=0;
    for i=1:numblocks
        if pin_front==0
            w_xc_blocks(j,i) = filtered_w_dim_out(j,round(filtered_xc(j,i)/1000)); %this line will throw an error for more than 10 blocks
        else
            xf=filtered_xc(j,i);
            w_xc_blocks(j,i) = filtered_w_dim_out(j,round(filtered_xc(j,i)/1000));
        end
    end
end

for j=1:size(x_filter_sort,1)
    for i=1:numblocks
        block_base(j,i) = w_xc_blocks(j,i)*(1-(rho_in/filtered_rho_c(j)));
    end
end 

for j=1:size(x_filter_sort,1)
    for i=1:numblocks
        topo(j,i) = block_base(j,i)+filtered_h_all(j,i);
    end
end 


filtered_topo_mean = mean(mean(topo, 2));
filtered_topo_std = std(mean(topo, 2));

% tabulate results
percent_fits = 100*size(x_filter_sort,1)/trials;

if get(handles.meters,'Value') == 1
Units = 'm';
else
Units = 'km';
end

if get(handles.comp_corr_yes,'Value') == 1
Comp_corr = 'yes';
else
Comp_corr = 'no';
end

x_min = get(handles.plot_opt_xmin,'String');
x_max = get(handles.plot_opt_xmax,'String');
y_min = get(handles.plot_opt_ymin,'String');
y_max = get(handles.plot_opt_ymax,'String');

Results = {'Model Results','Model fits (%)', 'Maximum load height mean (m)', 'Maximum load height standard deviation (m)', 'Topography height mean (m)', 'Topography height standard deviation (m)', 	'Topography slope mean (degrees)', 'Topography slope standard deviation (degrees)', 'Load width mean (m)', 'Load width standard deviation (m)', 	'Block width mean (m)', 'Block width standard deviation (m)', 'Load density mean (kg/m^3)', 'Load density standard deviation (kg/m^3)',	'Effective elastic thickness mean (m)', 'Effective elastic thickness standard deviation (m)','mean coefficient (exponential EET gradient only)', 'Coefficient standard deviation','Flexural rigidity mean (Nm)', 'Flexural rigidity standard deviation (Nm)', 'Flexural parameter mean (Nm)', 'Flexural parameter standard deviation (Nm)', [],'Model Parameters','Number of trials', 'Model Space', 'Plate type', 'Load pin front or back', 'Pin Location (m)', 'Topographic slope type', 'Max Angle','Number of blocks', 'Infill density (kg/m^3)', 'Mantle density (kg/m^3)', 'Gravity', 'Young''s Modulus (Pa)', 'Poisson''s ratio', 'Load height minimum (m)', 'Load height maximum (m)', 'Load width minimum (m)', 'Load width maximum (m)', 'Load density minimum (kg/m^3)', 'Load density maximum (kg/m^3)', 'Effective elastic thickness minimum (m)', 	'Effective elastic thickness maximum (m)', 'Effective Elastic Thickness gradient', 'Coefficient Min (exponential gradient only)','Coefficient Max','Subsidence uncertainty orientation','Filter 1 subsidence amount (m)', 'Filter 1 subsidence location (m)', 'Filter 1 subsidence uncertainty (m)', 'Filter 2 subsidence amount (m)', 'Filter 2 subsidence location (m)', 'Filter 2 subsidence uncertainty (m)', 	'Filter 3 subsidence amount (m)', 'Filter 3 subsidence location (m)', 'Filter 3 subsidence uncertainty (m)', 'Filter 4 subsidence amount (m)', 'Filter 4 subsidence location (m)', 'Filter 4 subsidence uncertainty (m)', 	'Filter 5 subsidence amount (m)', 'Filter 5 subsidence location (m)', 'Filter 5 subsidence uncertainty (m)', 'Units', 'Plot x-axis minimum (m)', 'Plot x-axis maximum (m)', 'Plot y-axis minimum (m)', 'Plot y-axis maximum (m)', 'Compaction-corrected'; 
            [],             percent_fits,   filtered_h_mean,                    filtered_h_std,                             filtered_topo_mean,                 filtered_topo_std,                                filtered_angle_mean,                  filtered_angle_std,                filtered_width_total_mean, filtered_width_total_std,         filtered_block_width_mean,  filtered_block_width_std,           	filtered_rho_c_mean,            filtered_rho_c_std,                         filtered_Te_mean,                       filtered_Te_std,                                          filtered_exponent_mean,                     filtered_exponent_std,               filtered_D_mean,                  filtered_D_std,                             	filtered_alpha_mean,         filtered_alpha_std,                  [],        [],             	trials,       model_space,    Plate_type,      Pin_location,           pin,                 Topo_type,            max_angle,       numblocks,              rho_in,                     rho_m,              	g,          E,                      v,                  h_min,                      h_max,                      width_min,          	width_max,                  rho_c_min,                      rho_c_max,                              Te_min,                                 Te_max,                                 EET_grad,                              exp_min,                                     exp_max,          Uncertainty_orientation ,        subsidence,                     subsidence_loc,                                 subuncertainty1,            	subsidence2,                        subsidence_loc2,                    subuncertainty2,                    	subsidence3,                        subsidence_loc3,                    subuncertainty3,                	subsidence4,                        subsidence_loc4,                    subuncertainty4,                    subsidence5,                        subsidence_loc5,                        subuncertainty5,                	Units,          x_min,                  x_max,                  y_min,                      y_max,                      Comp_corr};

Results = transpose(cell(Results));

set(handles.uitable1,'data', Results)

Model_fits = size(x_filter_sort,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     plot filtered model results                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x and y variables for plotting zero elevation
xx2 = 0:100:max(x_dim);
yy2(1,1:length(xx2)) = 0;

%plot figure with load geometry
    cla(handles.plot_results,'reset');
    axes(handles.plot_results);
    hold on
    %plot flexural profiles
    for i = 1:Model_fits
       plot(x_dim, filtered_w_dim_out(i,:), 'Color','k')
    end
    %plot control points for horizontal uncertainties
switch uncert_orientation
	case handles.Uncert_horz
        plot(xx2,yy2, 'LineWidth', 1,'Color','k')
        if(isnan(subsidence))==0
            scatter(((max_filter-min_filter)/2+min_filter),-subsidence, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter,-subsidence, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter,-subsidence, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence2))==0
            scatter(((max_filter2-min_filter2)/2+min_filter2),-subsidence2, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter2,-subsidence2, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter2,-subsidence2, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence3))==0
            scatter(((max_filter3-min_filter3)/2+min_filter3),-subsidence3, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter3,-subsidence3, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter3,-subsidence3, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence4))==0
            scatter(((max_filter4-min_filter4)/2+min_filter4),-subsidence4, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter4,-subsidence4, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter4,-subsidence4, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence5))==0
            scatter(((max_filter5-min_filter5)/2+min_filter5),-subsidence5, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter5,-subsidence5, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter5,-subsidence5, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
end
%plot control points for horizontal uncertainties
switch uncert_orientation
	case handles.Uncert_vert
        plot(xx2,yy2, 'LineWidth', 1,'Color','k')
        if(isnan(subsidence))==0
            scatter(xfilt_subsidence_loc,-subsidence, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc,-subsidence-subuncertainty1, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc,-subsidence+subuncertainty1, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence2))==0
            scatter(xfilt_subsidence_loc2,-subsidence2, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc2,-subsidence2-subuncertainty2, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc2,-subsidence2+subuncertainty2, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence3))==0
            scatter(xfilt_subsidence_loc3,-subsidence3, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc3,-subsidence3-subuncertainty3, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc3,-subsidence3+subuncertainty3, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence4))==0
            scatter(xfilt_subsidence_loc4,-subsidence4, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc4,-subsidence4-subuncertainty4, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc4,-subsidence4+subuncertainty4, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence5))==0
            scatter(xfilt_subsidence_loc5,-subsidence5, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc5,-subsidence5-subuncertainty5, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc5,-subsidence5+subuncertainty5, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
end
%END plot control points


%plot load blocks
rad_on=get(handles.ui_plot_options,'selectedobject');
switch rad_on
	case handles.comp_corr_no
    Comp_corr = 'no';
    
    for j=1:size(x_filter_sort,1)
        for i=1:numblocks
        %rectangle('Position',[(i-1)*loadwidth_all(j,i),block_base(j,i),loadwidth_all(j,i),h_all(i,j)])
        rectangle('Position',[filtered_xc(j,i)-filtered_block_width(j)/2 0 filtered_loadwidth_all(j,i) filtered_h_all(j,i)])
        end
    end
case handles.comp_corr_yes
    Comp_corr = 'yes';
    
    for j=1:size(x_filter_sort,1)
        for i=1:numblocks
        %rectangle('Position',[(i-1)*loadwidth_all(j,i),block_base(j,i),loadwidth_all(j,i),h_all(i,j)])
        rectangle('Position',[filtered_xc(j,i)-filtered_block_width(j)/2,block_base(j,i),...
            filtered_loadwidth_all(j,i),filtered_h_all(j,i)])
        end
    end
end

%plot figure with EET_profiles
cla(handles.EET_profiles,'reset');
axes(handles.EET_profiles);
hold on
%plot EET_profiles
for i=1:size(filtered_Alpha_x,1)
    for j=1:size(filtered_Alpha_x,2)
        Te_plotting(i,j)=(3*filtered_Alpha_x(i,j)^4*delta_rho*g*(1-v^2)/E)^(1/3);
    end
end
for i = 1:Model_fits
    plot(x_dim, Te_plotting(i,:), 'Color','k')
end

if min(Te_plotting)==max(Te_plotting)
    axis([str2num(get(handles.plot_opt_xmin,'String')), str2num(get(handles.plot_opt_xmax,'String')),...
        floor(min(min(Te_plotting-1000))), ceil(max(max(Te_plotting+1000)))])
else
    axis([str2num(get(handles.plot_opt_xmin,'String')), str2num(get(handles.plot_opt_xmax,'String')),...
        floor(min(min(Te_plotting))), ceil(max(max(Te_plotting)))])
end
e=title('Effective Elastic Thickness Profiles');
set(e, 'FontSize', 10);
x_label=xlabel('Distance (m)','FontSize', 10)
set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.100, 0]);
ylabel('EET (m)','FontSize', 10)
hold off
%END plot EET profilles

%Set parameters for Flexure Profiles axes
axes(handles.plot_results);
axis([str2num(get(handles.plot_opt_xmin,'String')), str2num(get(handles.plot_opt_xmax,'String')), ...
str2num(get(handles.plot_opt_ymin,'String')), str2num(get(handles.plot_opt_ymax,'String'))])
xt = get(gca,'XTick');
set(gca,'XTickLabel', xt)
yt = get(gca,'YTick');
set(gca,'YTickLabel', yt)
 
t=title('Flexure Model Results');
set(t, 'FontSize', 12);
x_label=xlabel('Distance (m)','FontSize', 10)
set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.050, 0]);
ylabel('Height (m)','FontSize', 10)
hold off
set(gca,'FontSize',10)
 
%Set axes for EET profiles
axes(handles.EET_profiles)

xt = get(gca,'XTick');
set(gca,'XTickLabel', xt)
 
yt = get(gca,'YTick');
set(gca,'YTickLabel', yt)
 
x_label=xlabel('Distance (m)','FontSize', 10)
set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.1, 0]);
ylabel('Effective Elastic Thickness (m)','FontSize', 10)
 
%Modify if km selected
rad_on=get(handles.ui_m_km,'selectedobject');

switch rad_on
case handles.km
% Modify Flexure Results axes
axes(handles.plot_results);
xt = get(gca,'XTick');
xt = num2cell(xt/1000);
set(gca,'XTickLabel', xt)
 
yt = get(gca,'YTick');
yt = num2cell(yt/1000)
set(gca,'YTickLabel', yt)
 
x_label=xlabel('Distance (km)','FontSize', 10)
set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.07, 0]);
ylabel('Height (km)','FontSize', 10)
 
%Modify EET profiles axes
axes(handles.EET_profiles)
xt = get(gca,'XTick');
xt = num2cell(xt/1000);
set(gca,'XTickLabel', xt)
 
yt = get(gca,'YTick');
yt = num2cell(yt/1000)
set(gca,'YTickLabel', yt)
 
 
x_label=xlabel('Distance (km)','FontSize', 10)
set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.1, 0]);
ylabel('Effective Elastic Thickness (km)','FontSize', 10)
end


%{
%set axis parameters for Flexure Model Results
axis([str2num(get(handles.plot_opt_xmin,'String')), str2num(get(handles.plot_opt_xmax,'String')), ...
str2num(get(handles.plot_opt_ymin,'String')), str2num(get(handles.plot_opt_ymax,'String'))])
t=title('Flexure Model Results');
set(t, 'FontSize', 10);
x_label=xlabel('Distance (m)','FontSize', 10)
set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.050, 0]);
ylabel('Elevation (m)','FontSize', 10)
hold off
set(gca,'FontSize',10)

rad_on=get(handles.ui_m_km,'selectedobject');
switch rad_on
    case handles.km
    axes(handles.plot_results);
    xt = get(gca,'XTick');
    xt = num2cell(xt/1000);
    set(gca,'XTickLabel', xt)

    yt = get(gca,'YTick');
    yt = num2cell(yt/1000)
    set(gca,'YTickLabel', yt)

    xlabel('Distance (km)','FontSize', 10)
    set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.050, 0]);
    ylabel('Elevation (km)','FontSize', 10)

end
%}
handles.x_dim = x_dim;
handles.Model_fits = Model_fits;
handles.filtered_w_dim_out = filtered_w_dim_out;
handles.x_filter_sort = x_filter_sort;
handles.numblocks = numblocks;
handles.xc = xc;
handles.block_width = block_width;
handles.block_base = block_base;
handles.loadwidth_all = loadwidth_all;
handles.filtered_xc = filtered_xc;
handles.filtered_block_width = filtered_block_width;
handles.block_base = block_base;
handles.filtered_loadwidth_all = filtered_loadwidth_all;
handles.h_all = h_all;
handles.filtered_h_all = filtered_h_all;
handles.subsidence = subsidence;
handles.subsidence2 = subsidence2;
handles.subsidence3 = subsidence3;
handles.subsidence4 = subsidence4;
handles.subsidence5 = subsidence5;
handles.min_filter = min_filter;
handles.min_filter2 = min_filter2;
handles.min_filter3 = min_filter3;
handles.min_filter4 = min_filter4;
handles.min_filter5 = min_filter5;
handles.max_filter = max_filter;
handles.max_filter2 = max_filter2;
handles.max_filter3 = max_filter3;
handles.max_filter4 = max_filter4;
handles.max_filter5 = max_filter5;
handles.pin=pin;
handles.filtered_Alpha_x=filtered_Alpha_x;
handles.delta_rho=delta_rho;
handles.g=g;
handles.v=v;
handles.E=E;
handles.filtered_Te = filtered_Te;
handles.filtered_Te_left = filtered_Te;
handles.filtered_Te_right = filtered_Te_max;
%handles.xfilt_subsidence_loc= xfilt_subsidence_loc;
handles.subuncertainty1=subuncertainty1;
handles.subuncertainty2=subuncertainty2;
handles.subuncertainty3=subuncertainty3;
handles.subuncertainty4=subuncertainty4;
handles.subuncertainty5=subuncertainty5;
handles.exp_minval=exp_min;
handles.exp_maxval=exp_max;
handles.filtered_exponent=filtered_exponent;
switch uncert_orientation
	case handles.Uncert_vert
    handles.xfilt_subsidence_loc=xfilt_subsidence_loc;
    handles.xfilt_subsidence_loc2=xfilt_subsidence_loc2;
    handles.xfilt_subsidence_loc3=xfilt_subsidence_loc3;
    handles.xfilt_subsidence_loc4=xfilt_subsidence_loc4;
    handles.xfilt_subsidence_loc5=xfilt_subsidence_loc5;
end
%{
if(isnan(subsidence2))==0
    handles.xfilt_subsidence_loc2=xfilt_subsidence_loc2;
    handles.subuncertainty2=subuncertainty2;
end
if(isnan(subsidence3))==0
    handles.xfilt_subsidence_loc3=xfilt_subsidence_loc3;
    handles.subuncertainty3=subuncertainty3;
end
if(isnan(subsidence4))==0
    handles.xfilt_subsidence_loc4=xfilt_subsidence_loc4;
    handles.subuncertainty4=subuncertainty4;
end
if(isnan(subsidence5))==0
    handles.xfilt_subsidence_loc5=xfilt_subsidence_loc5;
    handles.subuncertainty5=subuncertainty5;
end
%}

handles.Results = Results;

guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         LOAD PARAMETERS BUTTON        %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in load_params.
function load_params_Callback(hObject, eventdata, handles)

try
    parameters_load_save=handles.parameters_load_save;
catch
    parameters_load_save=pwd;
end
fullpath=horzcat(parameters_load_save,'\*.xls');
[filename parameters_load_save] = uigetfile({'*'},'File Selector', fullpath);

[numbers text, params] = xlsread(strcat(parameters_load_save, filename));

%set trials
i=2;
set(handles.model_trials, 'String', cell2mat(params(i,2)));
i=i+1;

%set model space
set(handles.model_space, 'String', cell2mat(params(i,2)));
i=i+1;

%set plate configuration
if contains(params(i,2), 'Inf') || contains(params(i,2), 'inf') == 1
set(handles.ui_model_type,'SelectedObject',handles.infinite_plate);
end
if contains(params(i,2), 'Br') || contains(params(i,2), 'br') == 1
set(handles.ui_model_type,'SelectedObject',handles.broken_plate);
end
i=i+1;

%set pin type
if contains(params(i,2), 'Fr') || contains(params(i,2), 'fr')== 1
set(handles.ui_pin_location,'SelectedObject',handles.front);
end
if contains(params(i,2), 'B') || contains(params(i,2), 'b') == 1
set(handles.ui_pin_location,'SelectedObject',handles.back);
end
i=i+1;

%set pin location
set(handles.pin_m, 'String', cell2mat(params(i,2)));
i=i+1;

%set load geometry
if contains(params(i,2), 'Rect') ||  contains(params(i,2), 'rect')== 1
set(handles.ui_topography,'SelectedObject',handles.topo_rect);
end
if contains(params(i,2), 'Ri') == 1||contains(params(i,2), 'ri')
set(handles.ui_topography,'SelectedObject',handles.topo_right);
end
if contains(params(i,2), 'Lef') == 1||contains(params(i,2), 'lef')
set(handles.ui_topography,'SelectedObject',handles.topo_left);
end
i=i+1;

set(handles.max_angle, 'String', cell2mat(params(i,2)));%max angle
i=i+1;
set(handles.num_blocks, 'String', cell2mat(params(i,2)));%number of blocks
i=i+1;
set(handles.infill_p, 'String', cell2mat(params(i,2)));%infill density
i=i+1;
set(handles.mantle_p, 'String', cell2mat(params(i,2)));%mantle density
i=i+1;
set(handles.gravity, 'String', cell2mat(params(i,2)));%gravity
i=i+1;
set(handles.youngs_mod, 'String', cell2mat(params(i,2)));%youngs modulus
i=i+1;
set(handles.Poissons, 'String', cell2mat(params(i,2)));%poissons ratio
i=i+1;
set(handles.h_min, 'String', cell2mat(params(i,2)));%min load height
i=i+1;
set(handles.h_max, 'String', cell2mat(params(i,2)));%max load height
i=i+1;
set(handles.width_min, 'String', cell2mat(params(i,2)));%min load width
i=i+1;
set(handles.width_max, 'String', cell2mat(params(i,2)));%max load width
i=i+1;
set(handles.density_min, 'String', cell2mat(params(i,2)));%min load density
i=i+1;
set(handles.density_max, 'String', cell2mat(params(i,2)));%max load density
i=i+1;
set(handles.Te_min, 'String', cell2mat(params(i,2)));%min EET 
i=i+1;
set(handles.Te_max, 'String', cell2mat(params(i,2)));%max EET
i=i+1;

if contains(params(i,2),'Const') || contains(params(i,2),'const') == 1
    set(handles.EET_gradient,'SelectedObject',handles.EET_constant);
elseif contains(params(i,2),'Lin') || contains(params(i,2),'lin') == 1
    set(handles.EET_gradient,'SelectedObject',handles.EET_linear);
elseif contains(params(i,2),'Exp') || contains(params(i,2),'exp') == 1
    set(handles.EET_gradient,'SelectedObject',handles.EET_exp);
    set(handles.exp_min,'String',cell2mat(params(i+1,2)));
    set(handles.exp_max,'String',cell2mat(params(i+2,2)));
elseif contains(params(i,2),'Step') || contains(params(i,2),'step') == 1
    set(handles.EET_gradient,'SelectedObject',handles.EET_stepped);
    set(handles.step_number,'String',cell2mat(params(i+1,2)));
else
    f = errordlg('Invalid parameter entry','Invalid Entry','modal');
end
i=i+3;

%uncertainty orientation
if contains(params(i,2),'V')||contains(params(i,2),'v') == 1
    set(handles.uncert_orientation,'SelectedObject',handles.Uncert_vert);
elseif contains(params(i,2), 'H')||contains(params(i,2),'h') == 1
    set(handles.uncert_orientation,'SelectedObject',handles.Uncert_horz);
else
    f = errordlg('Invalid uncertainty orientation','Invalid Entry','modal');
end
i=i+1;

set(handles.filter1_amount, 'String', cell2mat(params(i,2)));%Filter 1 amount
i=i+1;
set(handles.filter1_loc, 'String', cell2mat(params(i,2)));%Filter 1 location
i=i+1;
set(handles.filter1_loc_unc, 'String', cell2mat(params(i,2)));%Filter 1 uncertainty
i=i+1;
set(handles.filter2_amount, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.filter2_loc, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.filter2_loc_unc, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.filter3_amount, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.filter3_loc, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.filter3_loc_unc, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.filter4_amount, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.filter4_loc, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.filter4_loc_unc, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.filter5_amount, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.filter5_loc, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.filter5_loc_unc, 'String', cell2mat(params(i,2)));
i=i+1;

%set units
if strcmp(params(i,2), 'm') || strcmp(params(i,2), 'meters')== 1
set(handles.ui_m_km,'SelectedObject',handles.meters);
end
if strcmp(params(i,2), 'km') || strcmp(params(i,2), 'kilometers')== 1
set(handles.ui_m_km,'SelectedObject',handles.km);
end
i=i+1;

set(handles.plot_opt_xmin, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.plot_opt_xmax, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.plot_opt_ymin, 'String', cell2mat(params(i,2)));
i=i+1;
set(handles.plot_opt_ymax, 'String', cell2mat(params(i,2)));
i=i+1;

if contains(params(i,2), 'Y') || contains(params(i,2), 'y')== 1
    set(handles.ui_plot_options,'SelectedObject',handles.comp_corr_yes);
elseif contains(params(i,2), 'N') || contains(params(i,2), 'n')== 1
    set(handles.ui_plot_options,'SelectedObject',handles.comp_corr_no);
else
    f = errordlg('Invalid parameter entry','Invalid Entry','modal');

end

handles.parameters_load_save=parameters_load_save;
guidata(hObject,handles);









% --- Executes on button press in example_params.
function example_params_Callback(hObject, eventdata, handles)

Input_parameters;





% --- Executes on button press in export_results.
function export_results_Callback(hObject, eventdata, handles)

Results = handles.Results;
try
    path=handles.path;
catch
    path=pwd;
end
fullpath=horzcat(path,'\*.xls');
[file,path] = uiputfile('*.xls','Save file', path);

handles.path=path;
    
xlswrite([path file], Results);
guidata(hObject,handles);










% --- Executes on button press in clear_and_plot.
function clear_and_plot_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     plot filtered model results                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
block_base = handles.block_base;
delta_rho=handles.delta_rho;
E=handles.E;
filtered_Alpha_x=handles.filtered_Alpha_x;
filtered_w_dim_out = handles.filtered_w_dim_out;
filtered_loadwidth_all = handles.filtered_loadwidth_all;
filtered_xc = handles.filtered_xc;
filtered_block_width = handles.filtered_block_width;
filtered_h_all = handles.filtered_h_all;
g=handles.g;
Model_fits = handles.Model_fits;
min_filter = handles.min_filter;
min_filter2 = handles.min_filter2;
min_filter3 = handles.min_filter3;
min_filter4 = handles.min_filter4;
min_filter5 = handles.min_filter5;
max_filter = handles.max_filter;
max_filter2 = handles.max_filter2;
max_filter3 = handles.max_filter3;
max_filter4 = handles.max_filter4;
max_filter5 = handles.max_filter5;
numblocks = handles.numblocks;
pin=handles.pin;
subsidence = handles.subsidence;
subsidence2 = handles.subsidence2;
subsidence3 = handles.subsidence3;
subsidence4 = handles.subsidence4;
subsidence5 = handles.subsidence5;
v=handles.v;
x_filter_sort = handles.x_filter_sort;

%%%%%%%%%%%%%%%% Get filter information %%%%%%%%%%%%%%%%%%%%%%%%%%
%Filter 1 Subsidence
if isempty(str2num(get(handles.filter1_amount,'String'))) == 1
    subsidence = NaN;
    set(handles.filter1_amount, 'String', NaN);
    else
    subsidence = str2num(get(handles.filter1_amount,'String')); %magnitude of tectonic subsidence used in filter 1 (m)
end

%Filter 2 Subsidence
if isempty(str2num(get(handles.filter2_amount,'String'))) == 1
    subsidence2 = NaN;
    set(handles.filter2_amount, 'String', NaN);
    else
    subsidence2 = str2num(get(handles.filter2_amount,'String')); %magnitude of tectonic subsidence used in filter 2 (m)
end

%Filter 3 Subsidence
if isempty(str2num(get(handles.filter3_amount,'String'))) == 1
    subsidence3 = NaN;
    set(handles.filter3_amount, 'String', NaN);
    else
    subsidence3 = str2num(get(handles.filter3_amount,'String')); %magnitude of tectonic subsidence used in filter 3 (m)
end

%Filter 4 Subsidence
if isempty(str2num(get(handles.filter4_amount,'String'))) == 1
    subsidence4 = NaN;
    set(handles.filter4_amount, 'String', NaN);
    else
    subsidence4 = str2num(get(handles.filter4_amount,'String')); %magnitude of tectonic subsidence used in filter 4 (m)
end

%Filter 5 Subsidence
if isempty(str2num(get(handles.filter5_amount,'String'))) == 1
    subsidence5 = NaN;
    set(handles.filter5_amount, 'String', NaN);
    else
    subsidence5 = str2num(get(handles.filter5_amount,'String')); %magnitude of tectonic subsidence used in filter 5 (m)
end

%Filter 1 Location
if isempty(str2num(get(handles.filter1_loc,'String'))) == 1
    subsidence_loc = NaN;
    set(handles.filter1_loc, 'String', NaN);
    else
    subsidence_loc = str2num(get(handles.filter1_loc,'String')); %subsidence at location 1 (m)
end

%Filter 2 Location
if isempty(str2num(get(handles.filter2_loc,'String'))) == 1
    subsidence_loc2 = NaN;
    set(handles.filter2_loc, 'String', NaN);
    else
    subsidence_loc2 = str2num(get(handles.filter2_loc,'String')); %subsidence at location 2 (m)
end

%Filter 3 Location
if isempty(str2num(get(handles.filter3_loc,'String'))) == 1
    subsidence_loc3 = NaN;
    set(handles.filter3_loc, 'String', NaN);
    else
    subsidence_loc3 = str2num(get(handles.filter3_loc,'String')); %subsidence at location 3 (m)
end

%Filter 4 Location
if isempty(str2num(get(handles.filter4_loc,'String'))) == 1
    subsidence_loc4 = NaN;
    set(handles.filter4_loc, 'String', NaN);
    else
    subsidence_loc4 = str2num(get(handles.filter4_loc,'String')); %subsidence at location 4 (m)
end

%Filter 5 Location
if isempty(str2num(get(handles.filter5_loc,'String'))) == 1
    subsidence_loc5 = NaN;
    set(handles.filter5_loc, 'String', NaN);
    else
    subsidence_loc5 = str2num(get(handles.filter5_loc,'String')); %subsidence at location 5 (m)
end

%Filter 1 Uncertainty
if isempty(str2num(get(handles.filter1_loc_unc,'String'))) == 1
    subuncertainty1 = NaN;
    set(handles.filter1_loc_unc, 'String', NaN);
    else
    subuncertainty1 = str2num(get(handles.filter1_loc_unc,'String')); %uncertainty in position of subsidence location 1 (m)
end

%Filter 2 Uncertainty
if isempty(str2num(get(handles.filter2_loc_unc,'String'))) == 1
    subuncertainty2 = NaN;
    set(handles.filter2_loc_unc, 'String', NaN);
    else
    subuncertainty2 = str2num(get(handles.filter2_loc_unc,'String')); %uncertainty in position of subsidence location 2 (m)
end

%Filter 3 Uncertainty
if isempty(str2num(get(handles.filter3_loc_unc,'String'))) == 1
    subuncertainty3 = NaN;
    set(handles.filter3_loc_unc, 'String', NaN);
    else
    subuncertainty3 = str2num(get(handles.filter3_loc_unc,'String')); %uncertainty in position of subsidence location 3 (m)
end

%Filter 4 Uncertainty
if isempty(str2num(get(handles.filter4_loc_unc,'String'))) == 1
    subuncertainty4 = NaN;
    set(handles.filter4_loc_unc, 'String', NaN);
    else
    subuncertainty4 = str2num(get(handles.filter4_loc_unc,'String')); %uncertainty in position of subsidence location 4 (m)
end

%Filter 5 Uncertainty
if isempty(str2num(get(handles.filter5_loc_unc,'String'))) == 1
    subuncertainty5 = NaN;
    set(handles.filter5_loc_unc, 'String', NaN);
    else
    subuncertainty5 = str2num(get(handles.filter5_loc_unc,'String')); %uncertainty in position of subsidence location 5 (m)
end
%%%END SET FILTER PARAMETERS

xfilt_subsidence_loc=round(subsidence_loc + pin,-3); %round subsidence to nearest km
xfilt_subsidence_loc2=round(subsidence_loc2 + pin,-3); %round subsidence to nearest km
xfilt_subsidence_loc3=round(subsidence_loc3 + pin,-3); %round subsidence to nearest km
xfilt_subsidence_loc4=round(subsidence_loc4 + pin,-3); %round subsidence to nearest km
xfilt_subsidence_loc5=round(subsidence_loc5 + pin,-3); %round subsidence to nearest km

model_space=str2num(get(handles.model_space,'String')); %model space (m)
x_dim = transpose(0:1000:model_space); %define model space (m) %define model space (m)

% x and y variables for plotting zero elevation
xx2 = 0:100:max(x_dim);
yy2(1,1:length(xx2)) = 0;

%plot figure with load geometry
    cla(handles.plot_results,'reset');
    axes(handles.plot_results);
    hold on
    %plot flexural profiles
    for i = 1:Model_fits
       plot(x_dim, filtered_w_dim_out(i,:), 'Color','k')
    end
    %plot control points for horizontal uncertainties
uncert_orientation=get(handles.uncert_orientation,'selectedobject');
switch uncert_orientation
	case handles.Uncert_horz
        plot(xx2,yy2, 'LineWidth', 1,'Color','k')
        if(isnan(subsidence))==0
            scatter(((max_filter-min_filter)/2+min_filter),-subsidence, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter,-subsidence, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter,-subsidence, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence2))==0
            scatter(((max_filter2-min_filter2)/2+min_filter2),-subsidence2, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter2,-subsidence2, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter2,-subsidence2, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence3))==0
            scatter(((max_filter3-min_filter3)/2+min_filter3),-subsidence3, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter3,-subsidence3, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter3,-subsidence3, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence4))==0
            scatter(((max_filter4-min_filter4)/2+min_filter4),-subsidence4, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter4,-subsidence4, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter4,-subsidence4, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence5))==0
            scatter(((max_filter5-min_filter5)/2+min_filter5),-subsidence5, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter5,-subsidence5, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter5,-subsidence5, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
end
%plot control points for horizontal uncertainties
switch uncert_orientation
	case handles.Uncert_vert
        plot(xx2,yy2, 'LineWidth', 1,'Color','k')
        if(isnan(subsidence))==0
            scatter(xfilt_subsidence_loc,-subsidence, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc,-subsidence-subuncertainty1, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc,-subsidence+subuncertainty1, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence2))==0
            scatter(xfilt_subsidence_loc2,-subsidence2, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc2,-subsidence2-subuncertainty2, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc2,-subsidence2+subuncertainty2, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence3))==0
            scatter(xfilt_subsidence_loc3,-subsidence3, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc3,-subsidence3-subuncertainty3, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc3,-subsidence3+subuncertainty3, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence4))==0
            scatter(xfilt_subsidence_loc4,-subsidence4, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc4,-subsidence4-subuncertainty4, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc4,-subsidence4+subuncertainty4, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence5))==0
            scatter(xfilt_subsidence_loc5,-subsidence5, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc5,-subsidence5-subuncertainty5, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc5,-subsidence5+subuncertainty5, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
end
%END plot control points

%Plot load blocks
rad_on=get(handles.ui_plot_options,'selectedobject');
switch rad_on
	case handles.comp_corr_no
    Comp_corr = 'no';
    for j=1:size(x_filter_sort,1)
        for i=1:numblocks
        %rectangle('Position',[(i-1)*loadwidth_all(j,i),block_base(j,i),loadwidth_all(j,i),h_all(i,j)])
        rectangle('Position',[filtered_xc(j,i)-filtered_block_width(j)/2 0 filtered_loadwidth_all(j,i) filtered_h_all(j,i)])
        end
    end
end

rad_on=get(handles.ui_plot_options,'selectedobject');
switch rad_on
	case handles.comp_corr_yes
    Comp_corr = 'yes';
    axes(handles.plot_results);
    hold on
    for j=1:size(x_filter_sort,1)
        for i=1:numblocks
        %rectangle('Position',[(i-1)*loadwidth_all(j,i),block_base(j,i),loadwidth_all(j,i),h_all(i,j)])
        rectangle('Position',[filtered_xc(j,i)-filtered_block_width(j)/2,block_base(j,i),filtered_loadwidth_all(j,i),filtered_h_all(j,i)])
        end
    end
end


%plot figure with EET_profiles
cla(handles.EET_profiles,'reset');
axes(handles.EET_profiles);
hold on
%plot EET_profiles
for i=1:size(filtered_Alpha_x,1)
    for j=1:size(filtered_Alpha_x,2)
        Te_plotting(i,j)=(3*filtered_Alpha_x(i,j)^4*delta_rho*g*(1-v^2)/E)^(1/3);
    end
end
for i = 1:Model_fits
    plot(x_dim, Te_plotting(i,:), 'Color','k')
end
hold off
if min(Te_plotting)==max(Te_plotting)
    axis([str2num(get(handles.plot_opt_xmin,'String')), str2num(get(handles.plot_opt_xmax,'String')),...
        floor(min(min(Te_plotting-1000))), ceil(max(max(Te_plotting+1000)))])
else
    axis([str2num(get(handles.plot_opt_xmin,'String')), str2num(get(handles.plot_opt_xmax,'String')),...
        floor(min(min(Te_plotting))), ceil(max(max(Te_plotting)))])
end
e=title('Effective Elastic Thickness Profiles');
set(e, 'FontSize', 10);
x_label=xlabel('Distance (m)','FontSize', 10)
set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.100, 0]);
ylabel('Effective Elastic Thickness (m)','FontSize', 10)
%END plot EET profilles

rad_on=get(handles.ui_m_km,'selectedobject');
switch rad_on
case handles.meters
%Set parameters for Flexure Profiles axes
axes(handles.plot_results);
axis([str2num(get(handles.plot_opt_xmin,'String')), str2num(get(handles.plot_opt_xmax,'String')), ...
str2num(get(handles.plot_opt_ymin,'String')), str2num(get(handles.plot_opt_ymax,'String'))])
xt = get(gca,'XTick');
set(gca,'XTickLabel', xt)
yt = get(gca,'YTick');
set(gca,'YTickLabel', yt)
 
t=title('Flexure Model Results');
set(t, 'FontSize', 12);
x_label=xlabel('Distance (m)','FontSize', 10)
set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.050, 0]);
ylabel('Height (m)','FontSize', 10)
hold off
set(gca,'FontSize',10)

%Modify axes for EET profiles
axes(handles.EET_profiles)
xt = get(gca,'XTick');
set(gca,'XTickLabel', xt)

yt = get(gca,'YTick');
set(gca,'YTickLabel', yt)

x_label=xlabel('Distance (m)','FontSize', 10)
set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.1, 0]);
ylabel('Effective Elastic Thickness (m)','FontSize', 10)
end

rad_on=get(handles.ui_m_km,'selectedobject');
switch rad_on
case handles.km
    %Set parameters for Flexure Profiles axes
    axes(handles.plot_results);
    axis([str2num(get(handles.plot_opt_xmin,'String')), str2num(get(handles.plot_opt_xmax,'String')), ...
    str2num(get(handles.plot_opt_ymin,'String')), str2num(get(handles.plot_opt_ymax,'String'))])
    xt = get(gca,'XTick');
    xt = num2cell(xt/1000);
    set(gca,'XTickLabel', xt)
    yt = get(gca,'YTick');
    yt = num2cell(yt/1000);
    set(gca,'YTickLabel', yt)
 
    t=title('Flexure Model Results');
    set(t, 'FontSize', 12);
    x_label=xlabel('Distance (m)','FontSize', 10)
    set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.050, 0]);
    ylabel('Height (m)','FontSize', 10)
    hold off
    set(gca,'FontSize',10)
    
    %Modify EET profiles axes
    axes(handles.EET_profiles)
    xt = get(gca,'XTick');
    xt = num2cell(xt/1000);
    set(gca,'XTickLabel', xt)

    yt = get(gca,'YTick');
    yt = num2cell(yt/1000)
    set(gca,'YTickLabel', yt)

    x_label=xlabel('Distance (km)','FontSize', 10)
    set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.1, 0]);
    ylabel('Effective Elastic Thickness (km)','FontSize', 10)

end











% --- Executes on button press in save_params.
function save_params_Callback(hObject, eventdata, handles)

Results = handles.Results;
params_out = Results(24:end,:);
i=2;
%set number of trials
params_out(i,2) = cellstr(get(handles.model_trials, 'String'));
i=i+1;

%set model space
params_out(i,2)=cellstr(get(handles.model_space, 'String'));
i=i+1;

%set plate configuration
if get(handles.infinite_plate, 'Value') == 1
params_out(i,2) = cellstr('Infinite');
else
params_out(i,2) = cellstr('Broken');
end
i=i+1;

%set pin type
if get(handles.front, 'Value') == 1
params_out(i,2) = cellstr('Front');
else
params_out(i,2) = cellstr('Back');
end
i=i+1;

%set pin location
params_out(i,2)=cellstr(get(handles.pin_m, 'String'));
i=i+1;

%set topography type
if get(handles.topo_rect,'Value') == 1
params_out(i,2) = cellstr('Rectangular');
end
if get(handles.topo_right,'Value') == 1
params_out(i,2) = cellstr('Right');
end
if get(handles.topo_left,'Value') == 1
params_out(i,2) = cellstr('Left');
end
i=i+1;

params_out(i,2) = cellstr(get(handles.max_angle,'String'));%%max angle
i=i+1;

params_out(i,2) = cellstr(get(handles.num_blocks, 'String'));%number of blocks
i=i+1;

params_out(i,2) = cellstr(get(handles.infill_p, 'String'));%infill density
i=i+1;

params_out(i,2) = cellstr(get(handles.mantle_p, 'String'));%mantle density 
i=i+1;

params_out(i,2) = cellstr(get(handles.gravity, 'String'));%gravity
i=i+1;

params_out(i,2) = cellstr(get(handles.youngs_mod, 'String'));%youngs modulus
i=i+1;

params_out(i,2) = cellstr(get(handles.Poissons, 'String'));%poisson's ratio
i=i+1;

params_out(i,2) = cellstr(get(handles.h_min, 'String'));%minimum load height (m)
i=i+1;

params_out(i,2) = cellstr(get(handles.h_max, 'String'));%max load height (m)
i=i+1;

params_out(i,2) = cellstr(get(handles.width_min, 'String'));%minimum load width (m)
i=i+1;

params_out(i,2) = cellstr(get(handles.width_max, 'String'));%max load width (m)
i=i+1;

params_out(i,2) = cellstr(get(handles.density_min, 'String'));%minimum load density (kgm^-3)
i=i+1;

params_out(i,2) = cellstr(get(handles.density_max, 'String'));%max load density (kgm^-3)
i=i+1;

params_out(i,2) = cellstr(get(handles.Te_min, 'String'));%minimum Effective Elastic Thickness (m)
i=i+1;

params_out(i,2) = cellstr(get(handles.Te_max, 'String'));%max Effective elastic thickness (m)
i=i+1;

%set EET gradient
if get(handles.EET_constant,'Value') == 1
params_out(i,2) = cellstr('Constant');
end
if get(handles.EET_exp,'Value') == 1
params_out(i,2) = cellstr('Exponential');
end
if get(handles.EET_linear,'Value') == 1
params_out(i,2) = cellstr('Linear');
end
if get(handles.EET_stepped,'Value') == 1
params_out(i,2) = cellstr('Stepped');
end
i=i+1;

%set EET parameters
if get(handles.EET_exp,'Value') == 1
params_out(i,2) = cellstr(get(handles.exp_min, 'String'));
params_out(i+1,2) = cellstr(get(handles.exp_max,'String'));
elseif get(handles.EET_stepped,'Value') == 1
params_out(i,2) = cellstr(get(handles.step_number, 'String'));
params_out(i+1,2) = cellstr('NaN');
else
params_out(i,2) = cellstr('NaN');
params_out(i+1,2) = cellstr('NaN');    
end
i=i+2;

%set uncertainty orientation
if get(handles.Uncert_vert,'Value')==1
    params_out(i,2) = cellstr('Vertical');
else
    params_out(i,2) = cellstr('Horizontal');
end
i=i+1;

params_out(i,2) = cellstr(get(handles.filter1_amount, 'String'));%Filter 1 amount
i=i+1;

params_out(i,2) = cellstr(get(handles.filter1_loc, 'String'));%Filter 1 Location
i=i+1;

params_out(i,2) = cellstr(get(handles.filter1_loc_unc, 'String'));%Filter 1 uncertainty
i=i+1;

params_out(i,2) = cellstr(get(handles.filter2_amount, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.filter2_loc, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.filter2_loc_unc, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.filter3_amount, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.filter3_loc, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.filter3_loc_unc, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.filter4_amount, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.filter4_loc, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.filter4_loc_unc, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.filter5_amount, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.filter5_loc, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.filter5_loc_unc, 'String'));
i=i+1;

%set units
if get(handles.meters,'Value') == 1
params_out(i,2) = cellstr('m');
else
params_out(i,2) = cellstr('km');
end
i=i+1;

params_out(i,2) = cellstr(get(handles.plot_opt_xmin, 'String'));%plotting extent
i=i+1;

params_out(i,2) = cellstr(get(handles.plot_opt_xmax, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.plot_opt_ymin, 'String'));
i=i+1;

params_out(i,2) = cellstr(get(handles.plot_opt_ymax, 'String'));
i=i+1;

%set compensation correction
if get(handles.comp_corr_yes,'Value') == 1
params_out(i,2) = cellstr('Yes');
else
params_out(i, 2) = cellstr('No');
end
i=i+1;

try
    parameters_load_save=handles.parameters_load_save;
catch
    parameters_load_save=pwd;
end
fullpath=horzcat(parameters_load_save,'\*.xls');


%[numbers text, params] = xlsread(strcat(load_path, filename));
[file,parameters_load_save] = uiputfile('*.xls','Save file', fullpath);
xlswrite([parameters_load_save file], params_out);

handles.parameters_load_save=parameters_load_save;
guidata(hObject,handles);










% --- Executes on button press in flexure_profiles.
function flexure_profiles_Callback(hObject, eventdata, handles)
delta_rho=handles.delta_rho;
g=handles.g;
v=handles.v;
E=handles.E;
x_dim = handles.x_dim;
filtered_w_dim_out = handles.filtered_w_dim_out;
filtered_Alpha_x=handles.filtered_Alpha_x;
filtered_Te_left = handles.filtered_Te_left;
filtered_Te_right = handles.filtered_Te_right;
filtered_Te_mean = handles.filtered_Te;

for i=1:size(filtered_Alpha_x,1)
    for j=1:size(filtered_Alpha_x,2)
    Te(i,j)=(3*filtered_Alpha_x(i,j)^4*delta_rho*g*(1-v^2)/E)^(1/3);
    end
end

mean_profile=transpose(mean(filtered_w_dim_out,1));
max_profile=transpose(min(filtered_w_dim_out,[],1));
min_profile=transpose(max(filtered_w_dim_out,[],1));
profiles = horzcat(x_dim,mean_profile, max_profile, min_profile,transpose(filtered_w_dim_out));
profiles_out=cell(size(profiles,1)+1, size(profiles,2));
headers1 = {'X-axis (m)','Mean Deflection','Max Deflection', 'Min Deflection','All Flexural Deflections (m)'};
profiles_cell = num2cell(profiles);
profiles_out(2:end,:)=profiles_cell;
profiles_out(1,1:length(headers1))=headers1;


EET=num2cell(transpose(Te));
EET_out=cell(size(EET,1)+1, size(EET,2)+4);
headers2 = {'X-axis (m)', 'Mean Effective Elastic Thickness (m)','Min Effective Elastic Thickness (m)','Max Effective Elastic Thickness (m)','Effective Elastic Thickness (m)',};
EET_mean=num2cell(transpose(mean(Te)));
ind_min=Te(:,1)-min(Te(:,1));
EET_min=num2cell(transpose(Te(find(~ind_min),:)));
ind_max=Te(:,1)-max(Te(:,1));
EET_max=num2cell(transpose(Te(find(~ind_max),:)));
x_dim=num2cell(x_dim);

%Trim min and max if longer than (x,1)
if size(EET_min,2)>1
    EET_min= EET_min(:,1);
end
if size(EET_max,2)>1
    EET_max= EET_max(:,1);
end
EET_out(2:end,1)=x_dim;
EET_out(2:end,2)=EET_mean;
EET_out(2:end,3)=EET_min;
EET_out(2:end, 4)=EET_max;
EET_out(2:end,5:end) = EET;
EET_out(1,1:5) = headers2;

%{
EET_trials = horzcat(transpose(filtered_Te_left),transpose(filtered_Te_right),filtered_Te_mean);
headers3 = {'Eff. El. Thickness adjacent to load (m)','Eff. El. Thickness adjacent to foreland (m)','Mean Eff. El. Thickness(m)'}
EET_endpoints_out=cell(size(EET_trials,1)+1, size(EET_trials,2));
EET_endpoints_cell=num2cell(EET_trials);
EET_endpoints_out(2:end, :)=EET_endpoints_cell;
EET_endpoints_out(1,1:3)=headers3;
%}

[file,path] = uiputfile('*.xls','Save file');
writecell(profiles_out, horzcat(path, file), 'Sheet','Flexural Profiles')
writecell(EET_out, horzcat(path, file), 'Sheet','Effective Elastic Thickness')




% --- Executes on button press in export_loads.
function export_loads_Callback(hObject, eventdata, handles)

Model_fits = handles.Model_fits;
filtered_xc = handles.filtered_xc;
filtered_block_width = handles.filtered_block_width;
filtered_h_all = handles.filtered_h_all;
numblocks = handles.numblocks;

loads_out = zeros(numblocks*2+6,Model_fits);

loads_out(2:numblocks+1,:) = transpose(filtered_xc);
loads_out(numblocks+4,:) = transpose(filtered_block_width);
loads_out(numblocks+7:end,:) = transpose(filtered_h_all)

loads_out = cellfun(@(x)x(logical(x)),num2cell(loads_out),'uni',false) %cell with empties where zeros were
loads_out(1,1) = cellstr('Load centers (m)')
loads_out(numblocks+3,1) = cellstr('Load widths (m)')
loads_out(numblocks+6,1) = cellstr('Load heights (m)')


[file,path] = uiputfile('*.xlsx','Save file');
xlswrite([path file], loads_out);



%loads_out(1,1) = cellstr('Load center (m)')




% --- Executes on button press in export_plot.
function export_plot_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     plot filtered model results                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_dim = handles.x_dim;
Model_fits = handles.Model_fits;
filtered_w_dim_out = handles.filtered_w_dim_out;
x_filter_sort = handles.x_filter_sort;
numblocks = handles.numblocks;
filtered_xc = handles.filtered_xc;
filtered_block_width = handles.filtered_block_width;
block_base = handles.block_base;
filtered_loadwidth_all = handles.filtered_loadwidth_all;
filtered_h_all = handles.filtered_h_all;
subsidence = handles.subsidence;
subsidence2 = handles.subsidence2;
subsidence3 = handles.subsidence3;
subsidence4 = handles.subsidence4;
subsidence5 = handles.subsidence5;
min_filter = handles.min_filter;
min_filter2 = handles.min_filter2;
min_filter3 = handles.min_filter3;
min_filter4 = handles.min_filter4;
min_filter5 = handles.min_filter5;
max_filter = handles.max_filter;
max_filter2 = handles.max_filter2;
max_filter3 = handles.max_filter3;
max_filter4 = handles.max_filter4;
max_filter5 = handles.max_filter5;
filtered_Alpha_x=handles.filtered_Alpha_x;
delta_rho=handles.delta_rho;
g=handles.g;
v=handles.v;
E=handles.E;

uncert_orientation=get(handles.uncert_orientation,'selectedobject');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     plot filtered model results                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x and y variables for plotting zero elevation
xx2 = 0:100:max(x_dim);
yy2(1,1:length(xx2)) = 0;

%plot figure with load geometry
figure
hold on
    %plot flexural profiles
    for i = 1:Model_fits
       plot(x_dim, filtered_w_dim_out(i,:), 'Color','k')
    end
    %plot control points for horizontal uncertainties
switch uncert_orientation
	case handles.Uncert_horz
        plot(xx2,yy2, 'LineWidth', 1,'Color','k')
        if(isnan(subsidence))==0
            scatter(((max_filter-min_filter)/2+min_filter),-subsidence, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter,-subsidence, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter,-subsidence, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence2))==0
            scatter(((max_filter2-min_filter2)/2+min_filter2),-subsidence2, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter2,-subsidence2, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter2,-subsidence2, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence3))==0
            scatter(((max_filter3-min_filter3)/2+min_filter3),-subsidence3, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter3,-subsidence3, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter3,-subsidence3, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence4))==0
            scatter(((max_filter4-min_filter4)/2+min_filter4),-subsidence4, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter4,-subsidence4, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter4,-subsidence4, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence5))==0
            scatter(((max_filter5-min_filter5)/2+min_filter5),-subsidence5, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(min_filter5,-subsidence5, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(max_filter5,-subsidence5, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
end
%plot control points for horizontal uncertainties
switch uncert_orientation
	case handles.Uncert_vert
        plot(xx2,yy2, 'LineWidth', 1,'Color','k')
        if(isnan(subsidence))==0
            xfilt_subsidence_loc=handles.xfilt_subsidence_loc;
            subuncertainty1=handles.subuncertainty1;

            scatter(xfilt_subsidence_loc,-subsidence, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc,-subsidence-subuncertainty1, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc,-subsidence+subuncertainty1, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence2))==0
            xfilt_subsidence_loc2=handles.xfilt_subsidence_loc2;
            subuncertainty2=handles.subuncertainty2;

            scatter(xfilt_subsidence_loc2,-subsidence2, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc2,-subsidence2-subuncertainty2, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc2,-subsidence2+subuncertainty2, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence3))==0
            xfilt_subsidence_loc3=handles.xfilt_subsidence_loc3;
            subuncertainty3=handles.subuncertainty3;

            scatter(xfilt_subsidence_loc3,-subsidence3, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc3,-subsidence3-subuncertainty3, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc3,-subsidence3+subuncertainty3, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence4))==0
            xfilt_subsidence_loc4=handles.xfilt_subsidence_loc4;
            subuncertainty4=handles.subuncertainty4;

            scatter(xfilt_subsidence_loc4,-subsidence4, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc4,-subsidence4-subuncertainty4, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc4,-subsidence4+subuncertainty4, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
        if(isnan(subsidence5))==0
            xfilt_subsidence_loc5=handles.xfilt_subsidence_loc5;
            subuncertainty5=handles.subuncertainty5;
            
            scatter(xfilt_subsidence_loc5,-subsidence5, 75, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc5,-subsidence5-subuncertainty5, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
            scatter(xfilt_subsidence_loc5,-subsidence5+subuncertainty5, 10, 'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1.5)
        end
end
%END plot control points


rad_on=get(handles.ui_plot_options,'selectedobject');
switch rad_on
	case handles.comp_corr_no
    Comp_corr = 'no';
    
    
    
    %plot load blocks
    for j=1:size(x_filter_sort,1)
        for i=1:numblocks
        %rectangle('Position',[(i-1)*loadwidth_all(j,i),block_base(j,i),loadwidth_all(j,i),h_all(i,j)])
        rectangle('Position',[filtered_xc(j,i)-filtered_block_width(j)/2 0 filtered_loadwidth_all(j,i) filtered_h_all(j,i)])
        end
    end


	case handles.comp_corr_yes
    Comp_corr = 'yes';
    
    %plot load blocks
    for j=1:size(x_filter_sort,1)
        for i=1:numblocks
        %rectangle('Position',[(i-1)*loadwidth_all(j,i),block_base(j,i),loadwidth_all(j,i),h_all(i,j)])
        rectangle('Position',[filtered_xc(j,i)-filtered_block_width(j)/2,block_base(j,i),filtered_loadwidth_all(j,i),filtered_h_all(j,i)])
        end
    end
end



axis([str2num(get(handles.plot_opt_xmin,'String')), str2num(get(handles.plot_opt_xmax,'String')), ...
str2num(get(handles.plot_opt_ymin,'String')), str2num(get(handles.plot_opt_ymax,'String'))])
t=title('Flexure Model Results');
set(t, 'FontSize', 10);
x_label=xlabel('Distance (m)','FontSize', 10)
set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.050, 0]);
ylabel('Elevation (m)','FontSize', 10)
hold off
set(gca,'FontSize',10)

%plot figure with EET_profiles
figure
hold on
%plot EET_profiles
for i=1:size(filtered_Alpha_x,1)
    for j=1:size(filtered_Alpha_x,2)
        Te_plotting(i,j)=(3*filtered_Alpha_x(i,j)^4*delta_rho*g*(1-v^2)/E)^(1/3);
    end
end
for i = 1:Model_fits
    plot(x_dim, Te_plotting(i,:), 'Color','k')
end

axis([str2num(get(handles.plot_opt_xmin,'String')), str2num(get(handles.plot_opt_xmax,'String')), ...
min(min(Te_plotting)), max(max(Te_plotting))])
e=title('Effective Elastic Thickness Profiles');
set(e, 'FontSize', 10);
x_label=xlabel('Distance (m)','FontSize', 10)
set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.100, 0]);
ylabel('EET (m)','FontSize', 10)
hold off
%END plot EET profilles

rad_on=get(handles.ui_m_km,'selectedobject');
switch rad_on
    case handles.km
    axes(handles.plot_results);
    xt = get(gca,'XTick');
    xt = num2cell(xt/1000);
    set(gca,'XTickLabel', xt)

    yt = get(gca,'YTick');
    yt = num2cell(yt/1000)
    set(gca,'YTickLabel', yt)

    xlabel('Distance (km)','FontSize', 10)
    set(x_label, 'Units', 'Normalized', 'Position', [0.5, -0.050, 0]);
    ylabel('Elevation (km)','FontSize', 10)

end








function model_trials_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function model_trials_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in results_listbox.
function results_listbox_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function results_listbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in infinite_plate.
function infinite_plate_Callback(hObject, eventdata, handles)

% --- Executes on button press in broken_plate.
function broken_plate_Callback(hObject, eventdata, handles)

function h_min_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function h_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in front.
function front_Callback(hObject, eventdata, handles)

% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)

% --- Executes on button press in topo_rect.
function topo_rect_Callback(hObject, eventdata, handles)

% --- Executes on button press in topo_left.
function topo_left_Callback(hObject, eventdata, handles)

% --- Executes on button press in topo_right.
function topo_right_Callback(hObject, eventdata, handles)

function h_max_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function h_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function width_min_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function width_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function width_max_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function width_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function density_min_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function density_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function density_max_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function density_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Te_min_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Te_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Te_max_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Te_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter1_amount_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter1_amount_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter1_loc_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter1_loc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter1_loc_unc_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter1_loc_unc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter2_amount_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter2_amount_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter2_loc_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter2_loc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter2_loc_unc_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter2_loc_unc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter3_amount_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter3_amount_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter3_loc_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter3_loc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter3_loc_unc_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter3_loc_unc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter4_amount_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter4_amount_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter4_loc_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter4_loc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter4_loc_unc_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter4_loc_unc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter5_amount_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter5_amount_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter5_loc_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter5_loc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filter5_loc_unc_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function filter5_loc_unc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function num_blocks_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function num_blocks_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function infill_p_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function infill_p_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mantle_p_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function mantle_p_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gravity_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function gravity_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function youngs_mod_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function youngs_mod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Poissons_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Poissons_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_opt_xmin_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function plot_opt_xmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_opt_xmax_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function plot_opt_xmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_opt_ymin_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function plot_opt_ymin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_opt_ymax_Callback(hObject, eventdata, handles)

function plot_opt_ymax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in comp_corr_yes.
function comp_corr_yes_Callback(hObject, eventdata, handles)

% --- Executes on button press in comp_corr_no.
function comp_corr_no_Callback(hObject, eventdata, handles)

function pin_m_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function pin_m_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pin_back_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function pin_back_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function model_space_Callback(hObject, eventdata, handles)
% hObject    handle to model_space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of model_space as text
%        str2double(get(hObject,'String')) returns contents of model_space as a double


% --- Executes during object creation, after setting all properties.
function model_space_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Uncert_vert.
function Uncert_vert_Callback(hObject, eventdata, handles)
% hObject    handle to Uncert_vert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Uncert_vert


% --- Executes on button press in Uncert_horz.
function Uncert_horz_Callback(hObject, eventdata, handles)
% hObject    handle to Uncert_horz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Uncert_horz





function step_number_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function step_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function exp_max_Callback(hObject, eventdata, handles)
% hObject    handle to exp_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of exp_max as text
%        str2double(get(hObject,'String')) returns contents of exp_max as a double


% --- Executes during object creation, after setting all properties.
function exp_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exp_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function exp_min_Callback(hObject, eventdata, handles)
% hObject    handle to exp_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of exp_min as text
%        str2double(get(hObject,'String')) returns contents of exp_min as a double


% --- Executes during object creation, after setting all properties.
function exp_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exp_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%        RESTORE DEFAULTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in defaults.
function defaults_Callback(hObject, eventdata, handles)
% hObject    handle to defaults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.plot_results,'reset');
cla(handles.EET_profiles,'reset');
emptydata={};
set(handles.uitable1,'Data',emptydata)

%set trials
i=2;
set(handles.model_space,'String', 300000);
set(handles.model_trials, 'String', 1000);
i=i+1;

%set plate configuration
set(handles.ui_model_type,'SelectedObject',handles.infinite_plate);
i=i+1;

%set pin point
set(handles.ui_pin_location,'SelectedObject',handles.back);
i=i+1;

%set pin location
set(handles.pin_m,'String', 0);
i=i+1;

%set load geometry
set(handles.ui_topography,'SelectedObject',handles.topo_right);
i=i+1;

set(handles.max_angle, 'String', 10);%max angle
i=i+1;
set(handles.num_blocks, 'String', 5);%number of blocks
i=i+1;
set(handles.infill_p, 'String', 2400);%infill density
i=i+1;
set(handles.mantle_p, 'String', 3300);%mantle density
i=i+1;
set(handles.gravity, 'String', 9.81);%gravity
i=i+1;
set(handles.youngs_mod, 'String', 7E10);%youngs modulus
i=i+1;
set(handles.Poissons, 'String', 0.25);%poissons ratio
i=i+1;
set(handles.h_min, 'String', 3000);%min load height
i=i+1;
set(handles.h_max, 'String', 5000);%max load height
i=i+1;
set(handles.width_min, 'String', 75000);%min load width
i=i+1;
set(handles.width_max, 'String', 125000);%max load width
i=i+1;
set(handles.density_min, 'String', 2500);%min load density
i=i+1;
set(handles.density_max, 'String', 2700);%max load density
i=i+1;
set(handles.Te_min, 'String', 10000);%min EET 
i=i+1;
set(handles.Te_max, 'String', 30000);%max EET
i=i+1;
set(handles.EET_gradient,'SelectedObject',handles.EET_constant);
set(handles.exp_min,'String',-0.1);
set(handles.exp_max,'String',5);
i=i+3;

%uncertainty orientation
set(handles.uncert_orientation,'SelectedObject',handles.Uncert_vert);
i=i+1;

set(handles.filter1_amount, 'String', 2500);%Filter 1 amount
i=i+1;
set(handles.filter1_loc, 'String', 125000);%Filter 1 location
i=i+1;
set(handles.filter1_loc_unc, 'String', 200);%Filter 1 uncertainty
i=i+1;
set(handles.filter2_amount, 'String', 'NaN');
i=i+1;
set(handles.filter2_loc, 'String', 'NaN');
i=i+1;
set(handles.filter2_loc_unc, 'String', 'NaN');
i=i+1;
set(handles.filter3_amount, 'String', 'NaN');
i=i+1;
set(handles.filter3_loc, 'String', 'NaN');
i=i+1;
set(handles.filter3_loc_unc, 'String', 'NaN');
i=i+1;
set(handles.filter4_amount, 'String', 'NaN');
i=i+1;
set(handles.filter4_loc, 'String', 'NaN');
i=i+1;
set(handles.filter4_loc_unc, 'String', 'NaN');
i=i+1;
set(handles.filter5_amount, 'String', 'NaN');
i=i+1;
set(handles.filter5_loc, 'String', 'NaN');
i=i+1;
set(handles.filter5_loc_unc, 'String', 'NaN');
i=i+1;

%set units

set(handles.ui_m_km,'SelectedObject',handles.km);
i=i+1;

set(handles.plot_opt_xmin, 'String', 0);
i=i+1;
set(handles.plot_opt_xmax, 'String', 300000);
i=i+1;
set(handles.plot_opt_ymin, 'String', -10000);
i=i+1;
set(handles.plot_opt_ymax, 'String', 10000);
i=i+1;


set(handles.ui_plot_options,'SelectedObject',handles.comp_corr_yes);



function max_angle_Callback(hObject, eventdata, handles)
% hObject    handle to max_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_angle as text
%        str2double(get(hObject,'String')) returns contents of max_angle as a double


% --- Executes during object creation, after setting all properties.
function max_angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
