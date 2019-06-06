function varargout = AudioCompression(varargin)
% AUDIOCOMPRESSION MATLAB code for AudioCompression.fig
%      AUDIOCOMPRESSION, by itself, creates a new AUDIOCOMPRESSION or raises the existing
%      singleton*.
%
%      H = AUDIOCOMPRESSION returns the handle to a new AUDIOCOMPRESSION or the handle to
%      the existing singleton*.
%
%      AUDIOCOMPRESSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUDIOCOMPRESSION.M with the given input arguments.
%
%      AUDIOCOMPRESSION('Property','Value',...) creates a new AUDIOCOMPRESSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AudioCompression_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AudioCompression_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AudioCompression

% Last Modified by GUIDE v2.5 21-Nov-2014 17:35:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AudioCompression_OpeningFcn, ...
                   'gui_OutputFcn',  @AudioCompression_OutputFcn, ...
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


% --- Executes just before AudioCompression is made visible.
function AudioCompression_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AudioCompression (see VARARGIN)

% Choose default command line output for AudioCompression
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AudioCompression wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AudioCompression_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global file_name;
%guidata(hObject,handles)
file_name=uigetfile({'*.wav'},'Select an Audio File');
fileinfo = dir(file_name);
SIZE = fileinfo.bytes;
Size = SIZE/1024;
[x,Fs,bits] = wavread(file_name);
xlen=length(x);
t=0:1/Fs:(length(x)-1)/Fs;
set(handles.text2,'string',Size);
%plot(t,x);
axes(handles.axes3) % Select the proper axes
plot(t,x)
set(handles.axes3,'XMinorTick','on')
grid on

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global file_name;
if(~ischar(file_name))
   errordlg('Please select Audio first');
else
[x,Fs,bits] = wavread(file_name);
xlen=length(x);
t=0:1/Fs:(length(x)-1)/Fs;
wavelet='haar';
level=5;
frame_size=2048;
psychoacoustic='on '; %if it is off it uses 8 bits/frame as default
wavelet_compression = 'on ';
heavy_compression='off';
compander='on ';
quantization ='on ';

% ENCODER 

step=frame_size;
N=ceil(xlen/step);
%computational variables
Cchunks=0;
Lchunks=0;
Csize=0;
PERF0mean=0;
PERFL2mean=0;
n_avg=0;
n_max=0;
n_0=0;
n_vector=[];
for i=1:1:N
if (i==N);
frame=x([(step*(i-1)+1):length(x)]);
else
frame=x([(step*(i-1)+1):step*i]);
end
%wavelet decomposition of the frame
[C,L] = wavedec(frame,level,wavelet);
%wavelet compression scheme
if wavelet_compression=='on '
[thr,sorh,keepapp] = ddencmp('cmp','wv',frame);
if heavy_compression == 'on '
thr=thr*10^6;
end
[XC,CXC,LXC,PERF0,PERFL2] = wdencmp('gbl',C, L, wavelet,level,thr,sorh,keepapp);
C=CXC;
L=LXC;
PERF0mean=PERF0mean + PERF0;
PERFL2mean=PERFL2mean+PERFL2;
end
%Psychoacoustic model
if psychoacoustic=='on '
P=10.*log10((abs(fft(frame,length(frame)))).^2);
Ptm=zeros(1,length(P));
%Inspect spectrum and find tones maskers
for k=1:1:length(P)
if ((k<=1) | (k>=250))
bool = 0;
elseif ((P(k)<P(k-1)) | (P(k)<P(k+1))),
bool = 0;
elseif ((k>2) & (k<63)),
bool = ((P(k)>(P(k-2)+7)) & (P(k)>(P(k+2)+7)));
elseif ((k>=63) & (k<127)),
bool = ((P(k)>(P(k-2)+7)) & (P(k)>(P(k+2)+7)) & (P(k)>(P(k-3)+7)) & (P(k)>(P(k+3)+7)));
elseif ((k>=127) & (k<=256)),
bool = ((P(k)>(P(k-2)+7)) & (P(k)>(P(k+2)+7)) & (P(k)>(P(k-3)+7)) & (P(k)>(P(k+3)+7)) & (P(k)>(P(k-4)+7)) & (P(k)>(P(k+4)+7)) &(P(k)>(P(k-5)+7)) & (P(k)>(P(k+5)+7)) & (P(k)>(P(k-6)+7)) &(P(k)>(P(k+6)+7)));
else
bool = 0;
end
if bool==1
Ptm(k)=10*log10(10.^(0.1.*(P(k-1)))+10.^(0.1.*(P(k)))+10.^(0.1.*P(k+1)));
end
end
sum_energy=0;
for k=1:1:length(Ptm)
sum_energy=10.^(0.1.*(Ptm(k)))+sum_energy;
end
E=10*log10(sum_energy/(length(Ptm)));
SNR=max(P)-E;
n=ceil(SNR/6.02);
if n<=3
n=4;
n_0=n_0+1;
end
if n>=n_max
n_max=n;
end
n_avg=n+n_avg;
n_vector=[n_vector n];
end
%Compander(compressor)
if compander=='on '
Mu=255;
C = compand(C,Mu,max(C),'mu/compressor');
end
%Quantization
if quantization=='on '
if psychoacoustic=='off'
n=8;
end
partition = [min(C):((max(C)-min(C))/2^n):max(C)];
codebook = [1 min(C):((max(C)-min(C))/2^n):max(C)];
[index,quant,distor] = quantiz(C,partition,codebook);
%find and correct offset
offset=0;
for j=1:1:N
if C(j)==0
offset=-quant(j);
break;
end
end
quant=quant+offset;
C=quant;
end
%Put together all the chunks
Cchunks=[Cchunks C]; 
Lchunks=[Lchunks L];
Csize=[Csize length(C)];
Encoder = round((i/N)*100); %indicator of progess
end
Cchunks=Cchunks(2:length(Cchunks));
%wavwrite(Cchunks,Fs,bits,'output1.wav')
Csize=[Csize(2) Csize(N+1)];
Lsize=length(L);
Lchunks=[Lchunks(2:Lsize+1) Lchunks((N-1)*Lsize+1:length(Lchunks))];
PERF0mean=PERF0mean/N; %indicator
PERFL2mean=PERFL2mean/N;%indicator
n_avg=n_avg/N;%indicator
n_max;%indicator
end_of_encoder='done';
xdchunks=0;
for i=1:1:N;
if i==N;
Cframe=Cchunks([((Csize(1)*(i-1))+1):Csize(2)+(Csize(1)*(i-1))]);
%Compander (expander)
if compander=='on '
if max(Cframe)==0
else
Cframe = compand(Cframe,Mu,max(Cframe),'mu/expander');
end
end
xd = waverec(Cframe,Lchunks(Lsize+2:length(Lchunks)),wavelet);
else
Cframe=Cchunks([((Csize(1)*(i-1))+1):Csize(1)*i]);
%Compander (expander)
if compander=='on '
if max(Cframe)==0
else
Cframe = compand(Cframe,Mu,max(Cframe),'mu/expander');
end
end
xd = waverec(Cframe,Lchunks(1:Lsize),wavelet);
end
xdchunks=[xdchunks xd];
Decoder = round((i/N)*100); %indicator of progess
end
xdchunks=xdchunks(2:length(xdchunks));
%distorsion = sum((xdchunks-x').^2)/length(x)
end_of_decoder='done';
%creating audio files with compressed schemes
wavwrite(xdchunks,Fs,bits,'output1.wav');
end_of_writing_file='done';%indicator of progess;
[x,Fs,bits] = wavread('output1.wav');
fileinfo = dir('output1.wav');
SIZE = fileinfo.bytes;
Size = SIZE/1024;
set(handles.text3,'string',Size)
xlen=length(x);
t=0:1/Fs:(length(x)-1)/Fs;
axes(handles.axes4) % Select the proper axes
plot(t,xdchunks)
set(handles.axes4,'XMinorTick','on')
grid on

[y1,fs1, nbits1,opts1]=wavread(file_name);
[y2,fs2, nbits2,opts2]=wavread('output1.wav');
[c1x,c1y]=size(y1);
[c2x,c2y]=size(y1);
if c1x ~= c2x
    disp('dimeonsions do not agree');
 else
 R=c1x;
 C=c1y;
  err = (sum(y1(2)-y2).^2)/(R*C);
 MSE=sqrt(err);
 MAXVAL=255;
  PSNR = 20*log10(MAXVAL/MSE); 
  MSE= num2str(MSE);
  if(MSE > 0)
  PSNR= num2str(PSNR);
  else
PSNR = 99;
end
fileinfo = dir(file_name);
SIZE = fileinfo.bytes;
Size = SIZE/1024;
fileinfo1 = dir('output1.wav');
SIZE1 = fileinfo1.bytes;
Size1 = SIZE1/1024;

CompressionRatio = Size/Size1;

  set(handles.text14,'string',PSNR)
  set(handles.text16,'string',MSE)
  set(handles.text17,'string',CompressionRatio)
  
end


end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
