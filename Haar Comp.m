function varargout = AudioCompression(varargin)
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
function AudioCompression_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

guidata(hObject, handles);

function varargout = AudioCompression_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;
function pushbutton1_Callback(hObject, eventdata, handles)
global file_name;
file_name=uigetfile({'*.wav'},'Select an Audio File');
fileinfo = dir(file_name);
SIZE = fileinfo.bytes;
Size = SIZE/1024;
[x,Fs,bits] = wavread(file_name);
xlen=length(x);
t=0:1/Fs:(length(x)-1)/Fs;
set(handles.text2,'string',Size);
axes(handles.axes3) 
plot(t,x)
set(handles.axes3,'XMinorTick','on')
grid on
function pushbutton2_Callback(hObject, eventdata, handles)
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
psychoacoustic='on ';
wavelet_compression = 'on ';
heavy_compression='off';
compander='on ';
quantization ='on ';

step=frame_size;
N=ceil(xlen/step);

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
[C,L] = wavedec(frame,level,wavelet);

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

if psychoacoustic=='on '
P=10.*log10((abs(fft(frame,length(frame)))).^2);
Ptm=zeros(1,length(P));
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
if compander=='on '
Mu=255;
C = compand(C,Mu,max(C),'mu/compressor');
end
if quantization=='on '
if psychoacoustic=='off'
n=8;
end
partition = [min(C):((max(C)-min(C))/2^n):max(C)];
codebook = [1 min(C):((max(C)-min(C))/2^n):max(C)];
[index,quant,distor] = quantiz(C,partition,codebook);
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
Cchunks=[Cchunks C]; 
Lchunks=[Lchunks L];
Csize=[Csize length(C)];
Encoder = round((i/N)*100); 
end
Cchunks=Cchunks(2:length(Cchunks));
Csize=[Csize(2) Csize(N+1)];
Lsize=length(L);
Lchunks=[Lchunks(2:Lsize+1) Lchunks((N-1)*Lsize+1:length(Lchunks))];
PERF0mean=PERF0mean/N; 
PERFL2mean=PERFL2mean/N;
n_avg=n_avg/N;
n_max;
end_of_encoder='done';
xdchunks=0;
for i=1:1:N;
if i==N;
Cframe=Cchunks([((Csize(1)*(i-1))+1):Csize(2)+(Csize(1)*(i-1))]);

if compander=='on '
if max(Cframe)==0
else
Cframe = compand(Cframe,Mu,max(Cframe),'mu/expander');
end
end
xd = waverec(Cframe,Lchunks(Lsize+2:length(Lchunks)),wavelet);
else
Cframe=Cchunks([((Csize(1)*(i-1))+1):Csize(1)*i]);

if compander=='on '
if max(Cframe)==0
else
Cframe = compand(Cframe,Mu,max(Cframe),'mu/expander');
end
end
xd = waverec(Cframe,Lchunks(1:Lsize),wavelet);
end
xdchunks=[xdchunks xd];
Decoder = round((i/N)*100); 
end
xdchunks=xdchunks(2:length(xdchunks));

end_of_decoder='done';
wavwrite(xdchunks,Fs,bits,'output1.wav');
end_of_writing_file='done';
[x,Fs,bits] = wavread('output1.wav');
fileinfo = dir('output1.wav');
SIZE = fileinfo.bytes;
Size = SIZE/1024;
set(handles.text3,'string',Size)
xlen=length(x);
t=0:1/Fs:(length(x)-1)/Fs;
axes(handles.axes4) 
plot(t,xdchunks)
set(handles.axes4,'XMinorTick','on')
grid on

[y1,fs1, nbits1,opts1]=wasvread(file_name);
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

function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
