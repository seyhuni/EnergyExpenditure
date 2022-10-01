clearvars; clc; close all;
%% 
% Eldeki veriler:
% 0 2 ve 4 derecede RT3, GARMİN, HR, AEE, REE, HY3, 3D Gait System verileri

% RT3: Kalori, Garmin: Kalori, HR: Kalp Ritmi, AEE: IC gerçek EE veri, HY3:
% Accelerometre verisi ve algoritma, 3D Gait: Yürüme Karakteristiği Çıkarım
% verisi

% Yapılacaklar

% Eldeki gerçek kalori değeri ile eşleştir ve algoritmanı HY3 ve Gait
% üzerinden adapte et! Hepsini AEE verisi ile kıyasla.
% Sonuçlar incelendiğinde eğimin etkisi, yürüme karakteristiğinin etkisi,
% EKG ve HR etkisi nasıl bir fark yaratıyor sun!

%% Gait Features 
% RT2  RT5	Rankle	Rheel	Rknee	Rhip  
% LT2  LT5	Lankle	Lheel	Lknee	Lhip
addpath(genpath('C:\Users\asus\Desktop\energy_expenditure_paper'));
cd('C:\Users\asus\Desktop\energy_expenditure_paper');
[GaitParameters,path] = uigetfile('*.*','.mat');
cd(path)
GaitFile=load('-mat',GaitParameters);
gaitdata=GaitFile.gait; % 19.98, 21,55, 24.28 Erdem Uysal Örneği
[s1,s2]=size(gaitdata);
for i=2:s1*s2
     if(isnan(gaitdata(i))==1) 
         gaitdata(i)=gaitdata(i-1);
     end
end
SR=60;
YP=gaitdata;
N=size(gaitdata);
xCoordinates=1:3:N(2);
yCoordinates=2:3:N(2);
zCoordinates=3:3:N(2);
%% Right Foot - Walking Characteristics
RtoeXs=YP(:,xCoordinates(1));
RtoeYs=YP(:,yCoordinates(1));
RtoeZs=YP(:,zCoordinates(1));
RheeXs=YP(:,xCoordinates(4));
RheeYs=YP(:,yCoordinates(4));
RheeZs=YP(:,zCoordinates(4));

RightWalkingAngle_Xcoord=RtoeXs-RheeXs;
RightWalkingAngle_Ycoord=RtoeYs-RheeYs;
RightWalkingAngle_Zcoord=RtoeZs-RheeZs;
RightWalkingAngle=atand(RightWalkingAngle_Zcoord./RightWalkingAngle_Xcoord);
RightWalkingAngle_last=mean(RightWalkingAngle);
%Heel Strike
threshol = max(RheeXs);
thresholdr=0.5*threshol;

y=1;
for i=1:60
    gecici=RheeXs(i-y+1:i+y);
    if(gecici(ceil(length(gecici)/2))==max(gecici) & gecici(ceil(length(gecici)/2))>thresholdr)
        stancefirstr(i) = RheeXs(i);
    end
end

stancefirstr1=find(stancefirstr==max(stancefirstr));

x=17;
for k=x:(length(RheeXs)-x)
    gecici=RheeXs(k-x+1:k+x);
    if(gecici(ceil(length(gecici)/2))==max(gecici) & gecici(ceil(length(gecici)/2))>thresholdr)
        stancesr(k) = RheeXs(k);
    end
end

figure
stancestartr=find(stancesr>0);
stancestartr1=[stancefirstr1 stancestartr];
plot(RheeXs)
hold on
plot(stancestartr1,RheeXs(stancestartr1),'r+','MarkerFaceColor','r')
str = 'Heel Strike';
text(stancestartr1,RheeXs(stancestartr1),str)

HeelStrikePointsR=stancestartr1;
% Toe Off
threshol1 = max(RtoeXs);
thresholdr1=0.5*threshol1;

y=1;
for i=1:60
    gecici=RtoeXs(i-y+1:i+y);
    if(gecici(ceil(length(gecici)/2))==min(gecici) & gecici(ceil(length(gecici)/2))<thresholdr1)
        swingfirstr(i) = RtoeXs(i);
    end
end

swingfirstr1=find(swingfirstr==max(swingfirstr));

x=15;
for k=x:(length(RtoeXs)-x)
    gecici=RtoeXs(k-x+1:k+x);
    if(gecici(ceil(length(gecici)/2))==min(gecici) & gecici(ceil(length(gecici)/2))<thresholdr1)
        swingsr(k) = RtoeXs(k);
    end
end

figure
swingstartr=find(swingsr>0);
swingfirstr1=[stancefirstr1 swingstartr];

hold on
plot(RtoeXs,'r')
hold on
plot(swingstartr,RtoeXs(swingstartr),'k+','MarkerFaceColor','k')
str = 'Toe Off';
text(swingstartr,RtoeXs(swingstartr),str)

ToeOffPointsR=swingstartr;
%% Left Foot - Walking Characteristics

LtoeXs=YP(:,xCoordinates(7));
LtoeYs=YP(:,yCoordinates(7));
LtoeZs=YP(:,zCoordinates(7));
LheeXs=YP(:,xCoordinates(10));
LheeYs=YP(:,yCoordinates(10));
LheeZs=YP(:,zCoordinates(10));

LeftWalkingAngle_Xcoord=LtoeXs-LheeXs;
LeftWalkingAngle_Ycoord=LtoeYs-LheeYs;
LeftWalkingAngle_Zcoord=LtoeZs-LheeZs;
LeftWalkingAngle=atand(LeftWalkingAngle_Zcoord./LeftWalkingAngle_Xcoord);
LeftWalkingAngle_last=mean(LeftWalkingAngle);
%Heel Strike
threshol = max(LheeXs);
thresholdl=0.5*threshol;

y=1;
for i=1:60
    gecici=LheeXs(i-y+1:i+y);
    if(gecici(ceil(length(gecici)/2))==max(gecici) & gecici(ceil(length(gecici)/2))>thresholdl)
        stancefirstl(i) = LheeXs(i);
    end
end

stancefirstl1=find(stancefirstl==max(stancefirstl));

x=17;
for k=x:(length(LheeXs)-x)
    gecici=LheeXs(k-x+1:k+x);
    if(gecici(ceil(length(gecici)/2))==max(gecici) & gecici(ceil(length(gecici)/2))>thresholdl)
        stancesl(k) = LheeXs(k);
    end
end

figure
stancestartl=find(stancesl>0);
stancestartl1=[stancefirstl1 stancestartl];
plot(LheeXs)
hold on
plot(stancestartl1,LheeXs(stancestartl1),'r+','MarkerFaceColor','r')
str = 'Heel Strike';
text(stancestartl1,LheeXs(stancestartl1),str)

HeelStrikePointsL=stancestartl1;

threshol1 = max(LtoeXs);
thresholdl1=0.5*threshol1;

y=1;
for i=1:60
    gecici=LtoeXs(i-y+1:i+y);
    if(gecici(ceil(length(gecici)/2))==min(gecici) & gecici(ceil(length(gecici)/2))<thresholdl1)
        swingfirstl(i) = LtoeXs(i);
    end
end

swingfirstl1=find(swingfirstl==max(swingfirstl));

x=17;
for k=x:(length(LtoeXs)-x)
    gecici=LtoeXs(k-x+1:k+x);
    if(gecici(ceil(length(gecici)/2))==min(gecici) & gecici(ceil(length(gecici)/2))<thresholdl1)
        swingsl(k) = LtoeXs(k);
    end
end

figure
swingstartl=find(swingsl>0);
swingfirstl1=[stancefirstl1 swingstartl];

hold on
plot(LtoeXs,'r')
hold on
plot(swingstartl,LtoeXs(swingstartl),'k+','MarkerFaceColor','k')
str = 'Toe Off';
text(swingstartl,LtoeXs(swingstartl),str)

ToeOffPointsL=swingstartl;
%% StanceveSwing faz belirleme

%HeelStrikeNoktalariSol ve ToeOffnoktalariSol ve sağ da mevcut
[b,m1,n1] = unique(HeelStrikePointsL,'first');
[c1,d1] =sort(m1);
b = b(d1);
HeelStrikePointsL=b;

if length(HeelStrikePointsL)>length(ToeOffPointsL)
    HeelStrikePointsL=HeelStrikePointsL(1:length(ToeOffPointsL));
end

StancePhaseLeft=HeelStrikePointsL-ToeOffPointsL(2:end);
% Burada ayağın yere bastıığı aralık. o yüzden önce heel strike gelirken
% sonra toe of gelmeli
SwingPhaseLeft=ToeOffPointsL(1:end-1)-HeelStrikePointsL(2:end);
% burada ise ayak havada kalacak. hava kalma işlemi toe dan başladığı için
% önce toe , sonra heel gelecek burda da. burada yapılan işler heel ve toe
% tespit edilen koordinatlara bağlı frame değerleri üzerinden yapıldığı
% için zaman tespiti yapıyoruz , FRAME cinsinden. Frame'ide 60'a bölersek
% sn cinsinden süreyi de hesaplayabiliriz.
StancePhaseRight=HeelStrikePointsR(2:end-1)-ToeOffPointsR(2:end);
%önce heel
SwingPhaseRight=ToeOffPointsR-HeelStrikePointsR(2:end);
%önce toe
%% Double Support Single Support Belirleme

%DS:Stance-(100-Stance)
% DSsol=abs(StancefaziSol)-(100-abs(StancefaziSol));
% DSsag=StancefaziSag-(100-StancefaziSag);
% %Single support:
% SingleSupportSol=ToeOffnoktalariSol-HeelStrikePointsSol(3:end);
% SingleSupportSag=ToeOffnoktalari-HeelStrikePoints(2:end);

firstDs=HeelStrikePointsR(2:end-1)-ToeOffPointsL(1:end);
secondDs=HeelStrikePointsL(1:end)-ToeOffPointsR(1:end-1);
singleSupRight=ToeOffPointsL(1:end-1)-HeelStrikePointsL(2:end);
singleSupLeft=ToeOffPointsR(1:end)-HeelStrikePointsR(2:end);


StanceRatioL=abs(StancePhaseLeft(1:end-1))/(abs(StancePhaseLeft(1:end-1)+SwingPhaseLeft));
StanceRatioR=abs(StancePhaseRight)/(abs(StancePhaseRight+SwingPhaseRight(1:end-1)));
% Paydayı stanceMeanl+swingMeanl+stanceMeanR+swingMeanR
StanceRatioL_FC=abs(StancePhaseLeft(1:end-1))/(abs(StancePhaseLeft(1:end-1)+SwingPhaseLeft+StancePhaseRight(1:end-1)+SwingPhaseRight(1:end-2)));
StanceRatioR_FC=abs(StancePhaseRight(1:end-1))/(abs(StancePhaseLeft(1:end-1)+SwingPhaseLeft+StancePhaseRight(1:end-1)+SwingPhaseRight(1:end-2)));

fDSRatio=abs(firstDs(1:end-1))/(abs(StancePhaseLeft(1:end-1)+SwingPhaseLeft+StancePhaseRight(1:end-1)+SwingPhaseRight(1:end-2)));
sDSRatio=abs(secondDs(1:end-1))/(abs(StancePhaseLeft(1:end-1)+SwingPhaseLeft+StancePhaseRight(1:end-1)+SwingPhaseRight(1:end-2)));
%% Spatio Temporal parameters

StepLengthL=RheeXs(HeelStrikePointsR(2:end-2))-RheeXs(HeelStrikePointsL(2:end));
MeanStepLengthL=mean(StepLengthL)/10;

StepLengthR=LheeXs(HeelStrikePointsL)-LheeXs(HeelStrikePointsR(2:end-1));
MeanStepLengthR=mean(StepLengthR)/10;

StepTimeL=(HeelStrikePointsR(2:end-2)-HeelStrikePointsL(2:end))/SR;
MeanStepTimeL=mean(abs(StepTimeL));

StepTimeR=(HeelStrikePointsL-HeelStrikePointsR(2:end-1))/SR;
MeanStepTimeR=mean(abs(StepTimeR));
StrideLength=StepLengthL+StepLengthR(1:end-1);
format longG
MeanStrideLength=mean(StrideLength)/1000;
StrideTimeL=0;
for i=1:length(HeelStrikePointsL)-1
    StrideTimeL(i)=(HeelStrikePointsL(i+1))-(HeelStrikePointsL(i));
end
StrideTimeR=0;
for i=1:length(HeelStrikePointsR)-1
    StrideTimeR(i)=(HeelStrikePointsR(i+1))-(HeelStrikePointsR(i));
end
StrideTimeL=mean(StrideTimeL/SR);
StrideTimeR=mean(StrideTimeR(2:end)/SR);
MeanStrideTime=(StrideTimeL+StrideTimeR)/2;

StepWidth=RheeYs(HeelStrikePointsR(2:end-2))-LheeYs(HeelStrikePointsL(2:end));
MeanStepWidth=mean(StepWidth)/10;
StepCountLeft=length(HeelStrikePointsR);
StepCountRight=length(HeelStrikePointsL);
TotalSteps=StepCountLeft+StepCountRight;
TotalSteps_Minute= 60*SR*TotalSteps/length(gaitdata);
% LFPA
LeftAciXs=LtoeXs-LheeXs;
LeftAciYs=LtoeYs-LheeYs;
LeftFootProgressionAngle=atand(LeftAciYs./LeftAciXs);
AverageLFPA=mean(LeftFootProgressionAngle);
% RFPA
RightAciXs=RtoeXs-RheeXs;
RightAciYs=RtoeYs-RheeYs;
RightFootProgressionAngle=atand(RightAciYs./RightAciXs);
AverageRFPA=mean(RightFootProgressionAngle);

GaitCell = cell(1,18);
GaitCell{1} = RightWalkingAngle_last;
GaitCell{2} = firstDs;
GaitCell{3} = secondDs;
GaitCell{4} = singleSupRight;
GaitCell{5} = singleSupLeft;
GaitCell{6} = fDSRatio;
GaitCell{7} = sDSRatio;
GaitCell{8} = MeanStepLengthL;
GaitCell{9} = MeanStepLengthR;
GaitCell{10} = MeanStepTimeL;
GaitCell{11} = MeanStepTimeR;
GaitCell{12} = MeanStrideLength;
GaitCell{13} = MeanStrideTime;
GaitCell{14} = MeanStepWidth;
GaitCell{15} = TotalSteps_Minute*7; 
GaitCell{16} = AverageLFPA;
GaitCell{17} = AverageRFPA;
GaitCell{18} = LeftWalkingAngle_last;

%% HY3 verileri

x = data(:,1);  %Ön Arka
y = data(:,2);   % sağ sol
z = data(:,3);   % yukarı aşağı
figure(1)

fftX = fft(x);
lengthOfFft = length(fftX);
fdiv = linspace(-fs/2,fs/2,lengthOfFft);
% plot(fdiv,abs(fftshift(fftX)))
stepSpeedRegionLow = 0.5; % 0.5 Hz step speed region low
stepSpeedRegionHigh = 2.5; % 2.5 Hz step speed region high
frequencyMarkLow1 = stepSpeedRegionLow/fs*lengthOfFft;
frequencyMarkLow2 = lengthOfFft - stepSpeedRegionLow/fs*lengthOfFft;
frequencyMarkHigh1 = stepSpeedRegionHigh/fs*lengthOfFft;
frequencyMarkHigh2 = lengthOfFft - stepSpeedRegionHigh/fs*lengthOfFft;
fftX(1:frequencyMarkLow1) = 0;
fftX(frequencyMarkHigh1:frequencyMarkHigh2) = 0;
fftX(frequencyMarkLow2:end) = 0;
filtered_step_x = real(ifft(fftX));
plot(timeDiv,filtered_step_x)
hold on

fftY = fft(y);
lengthOfFft = length(fftY);
fdiv = linspace(-fs/2,fs/2,lengthOfFft);
% plot(fdiv,abs(fftshift(fftY)))
stepSpeedRegionLow = 0.5; % 0.5 Hz step speed region low
stepSpeedRegionHigh = 2.5; % 2.5 Hz step speed region high
frequencyMarkLow1 = stepSpeedRegionLow/fs*lengthOfFft;
frequencyMarkLow2 = lengthOfFft - stepSpeedRegionLow/fs*lengthOfFft;
frequencyMarkHigh1 = stepSpeedRegionHigh/fs*lengthOfFft;
frequencyMarkHigh2 = lengthOfFft - stepSpeedRegionHigh/fs*lengthOfFft;
fftY(1:frequencyMarkLow1) = 0;
fftY(frequencyMarkHigh1:frequencyMarkHigh2) = 0;
fftY(frequencyMarkLow2:end) = 0;
filtered_step_y = real(ifft(fftY));
plot(timeDiv,filtered_step_y,'r')
hold on

fftZ = fft(z);
lengthOfFft = length(fftZ);
fdiv = linspace(-fs/2,fs/2,lengthOfFft);
% plot(fdiv,abs(fftshift(fftY)))
stepSpeedRegionLow = 0.5; % 0.5 Hz step speed region low
stepSpeedRegionHigh = 2.5; % 2.5 Hz step speed region high
frequencyMarkLow1 = stepSpeedRegionLow/fs*lengthOfFft;
frequencyMarkLow2 = lengthOfFft - stepSpeedRegionLow/fs*lengthOfFft;
frequencyMarkHigh1 = stepSpeedRegionHigh/fs*lengthOfFft;
frequencyMarkHigh2 = lengthOfFft - stepSpeedRegionHigh/fs*lengthOfFft;
fftZ(1:frequencyMarkLow1) = 0;
fftZ(frequencyMarkHigh1:frequencyMarkHigh2) = 0;
fftZ(frequencyMarkLow2:end) = 0;
filtered_step_z = real(ifft(fftZ));
plot(timeDiv,filtered_step_z,'g')
title("De-noise & Time decomposition")


figure(2)

fftX = fft(x);
fftY = fft(y);
fftZ = fft(z);
lengthOfFft = length(fftX);
fdiv = linspace(-fs/2,fs/2,lengthOfFft);
stepSpeedRegionLow = 0.01; % Only Constant companents
frequencyMarkLow1 = stepSpeedRegionLow/fs*lengthOfFft;
frequencyMarkLow2 = lengthOfFft - stepSpeedRegionLow/fs*lengthOfFft;
fftX(frequencyMarkLow1:end) = 0;
fftY(frequencyMarkLow1:end) = 0;
fftZ(frequencyMarkLow1:end) = 0;
filtered_magnitude_x = real(ifft(fftX));
filtered_magnitude_y = real(ifft(fftY));
filtered_magnitude_z = real(ifft(fftZ));
plot(timeDiv,filtered_magnitude_x)
hold on
plot(timeDiv,filtered_magnitude_y,'r')
hold on
plot(timeDiv,filtered_magnitude_z,'g')
title("De-noise & Magnitude decomposition")

figure(3)
alpha = atand(filtered_magnitude_x./sqrt(filtered_magnitude_y.^2+filtered_magnitude_z.^2));
beta = atand(filtered_magnitude_y./sqrt(filtered_magnitude_x.^2+filtered_magnitude_z.^2));
theta = atand(sqrt(filtered_magnitude_x.^2+filtered_magnitude_y.^2)./(filtered_magnitude_z.^2));
plot(timeDiv,alpha)
hold on
plot(timeDiv,beta,'r')
hold on
plot(timeDiv,theta,'g')
title("Filtered Tilt Angles decomposition (degrees)")

figure (4)

vectormag=sqrt ((filtered_step_x .* filtered_step_x)+(filtered_step_y .* filtered_step_y)+(filtered_step_z .* filtered_step_z));
plot(timeDiv, vectormag)

mean(vectormag)

vectormagnew=vectormag/max(vectormag); 
% initialize filtered signal
eogF = vectormag;

% TKEO basic  % Teager–Kaiser energy operator 
for i=2:length(eogF)-1
    eogF(i) = vectormagnew(i)^2 - vectormagnew(i-1)*vectormagnew(i+1);
end

[c,l] = wavedec(eogF,8,'db6'); % 8 level decomposition  

for t=1:8
    D(:,t)=wrcoef('d',c,l,'db6',t);
end
for t=1:8
    A(:,t)=wrcoef('a',c,l,'db6',t);
end

A8=A(:,8);  %


denemece=filtered_step_x;
denemece2=denemece.^2;

steps_negative = zeros(1,length(filtered_step_x));
for k=1:length(filtered_step_x)   
if filtered_step_x(k)<0
steps_negative(k)=filtered_step_x(k);
end
end

steps_negative2=steps_negative.^2;

steps_last=steps_negative2;
threshold=mean(steps_last);
% threshold=threshold*0.01;
stepresult=0;
x=5;

for k=x:(length(steps_last)-x)
        gecici=steps_last(k-x+1:k+x);
        if(gecici(ceil(length(gecici)/2))==max(gecici) & gecici(ceil(length(gecici)/2))>threshold)
        stepresult(k) = steps_last(k);
    end
end
stepresult_last=find(stepresult>0);

figure;
hold on;
plot(steps_last)
plot(stepresult_last,steps_last(stepresult_last),'r+','MarkerFaceColor','r')

HY3Cell= cell(1,8);
HY3Cell{1} = length(stepresult_last)*2;
HY3Cell{2} = 
