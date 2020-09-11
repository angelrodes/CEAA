function [] = start()
% Cosmogenic Exposure Age Averages (CEAA)
% Shows different methods to calculate an average (average ± deviation)
% from a set of gaussian data (typically Cosmogenic exposure ages).
% Input data in .csv file containing 4 columns: 
% sample names ; ages ; internal uncertainties ; external uncertainties
% First line of the csv file (header) is ignored.
%
% Angel Rodes, 2020
% www.angelrodes.com

%% clear previous data and plots
clear
close all hidden
clc

%% Version
scriptversion='1.2';


%% Import data
[file,path] = uigetfile({'*.csv','*.CSV'}, 'Select input file');
selectedfile = fullfile(path,file);
fid = fopen(selectedfile);
% Sample name , age , internal uncertainty , external uncertainty
mydata = textscan(fid, '%s %f %f %f',...
    'HeaderLines', 1,'Delimiter',',');
fclose(fid);
names=mydata{1};
ages=mydata{2};
errors=mydata{3};
exterrors=mydata{4};

%% Generate fake data (testing)

% for n=1:2+round(rand*15)
%     names{n}=['Sample-' num2str(n)];
%     ages(n)=15+rand*5;
%     errors(n)=ages(n)*(0.02+5*rand/100);
%     exterrors(n)=(errors(n)^2+(ages(n)*0.082)^2)^0.5;
% end

%% Calculate calibration error
caliberrors=((exterrors).^2-(errors).^2).^0.5;
calibpercents=caliberrors./ages;
calibpercent=median(calibpercents);
if std(calibpercents)/mean(calibpercents)>0.1
    warning (['Calibration error between ' num2str(min(calibpercents)*100) ' and ' num2str(max(calibpercents)*100) '%. Setting calibration error to ' num2str(calibpercent*100) '%']);
end

%% Calculate agerage
Av=mean(ages);
dAv=std(ages);
extdAv=(dAv^2+(Av*calibpercent)^2)^0.5;

%% Calculate weighted average
Wa=sum(ages./errors.^2)/sum(1./errors.^2);
variance=sum((ages-Wa).^2./errors.^2)/sum(1./errors.^2);
dWa=variance^0.5;
extdWa=(dWa^2+(Wa*calibpercent)^2)^0.5;

%% Calculate camelplot
% define space
nx=10000;
% minx=min(min(ages-errors*3),mean(ages)-3*std(ages));
% maxx=max(max(ages+errors*3),mean(ages)+3*std(ages));
minx=min(Av-dAv*2,Wa-dWa*2);
maxx=max(Av+dAv*2,Wa+dWa*2);
x=linspace(minx,maxx,nx);
y=0.*(x);

% define gussian pdf
normalprob=@( x, mu, s )...
    1/(s*(2*pi)^0.5)*...
    exp(-(x-mu).^2/(2*s^2));

% camelplot sum
for n=1:length(ages)
    mu = ages(n);
    sd = errors(n);
    iy = normalprob(x, mu, sd);
    y=y+iy;
end

%% Calculate one sigma limits
onesigmapercentage=sum(normalprob(-1000:1000,0,1000));

sortedy = sort(y,'descend');
prob=0;
n=0;
while prob<onesigmapercentage
    n=n+1;
    prob=sum(sortedy(1:n))/sum(y);
end
ylimit=sortedy(n);
xonesigma=x;
xonesigma(y<ylimit)=NaN;

%% Calculate best gaussian fit
[Bgf,dBgf] = BFG(x,y);
extdBgf=(dBgf^2+(Bgf*calibpercent)^2)^0.5;
aBgf=(y/normalprob( x, Bgf, dBgf));

%% Calculate possible outliers
outliers= ( (ages-Bgf).^2>2*(errors.^2+dBgf.^2) );

%% Display data and results
% Display credits
disp(['Cosmogenic Exposure Age Averages (CEAA) v.' scriptversion])
disp('Angel Rodes, 2020. www.angelrodes.com')
disp(' ')

% Display samples
disp('Samples:')
for n=1:length(ages)
    roundednumbers = significant_figures([ages(n),errors(n),exterrors(n)]);
    strnames{n}=[names{n} ': ' num2str(roundednumbers(1)) ' ± ' num2str(roundednumbers(2)) ' (' num2str(roundednumbers(3)) ')'];
    disp(['    ' strnames{n}])
end
disp('----------------------------------------------')

% Display average
roundednumbers = significant_figures([Av,dAv,extdAv]);
strAv=['Average: ' num2str(roundednumbers(1)) ' ± ' num2str(roundednumbers(2)) ' (' num2str(roundednumbers(3)) ')'];
disp(strAv)

% Display weighted average
roundednumbers = significant_figures([Wa,dWa,extdWa]);
strWa=['Waighted average: ' num2str(roundednumbers(1)) ' ± ' num2str(roundednumbers(2)) ' (' num2str(roundednumbers(3)) ')'];
disp(strWa)

% Display Best-Gaussian-Fit
roundednumbers = significant_figures([Bgf,dBgf,extdBgf]);
strBgf=['Best gaussin fit: ' num2str(roundednumbers(1)) ' ± ' num2str(roundednumbers(2)) ' (' num2str(roundednumbers(3)) ')'];
disp(strBgf)


% Display one sigma range
roundednumbers = significant_figures([min(xonesigma),max(xonesigma),dWa]);
strrange=['One-sigma range: [ ' num2str(roundednumbers(1)) ' - ' num2str(roundednumbers(2)) ' ]'];
disp(strrange)

% Display outliers
if sum(outliers)>0
    disp(' ')
    disp('Apparent outliers:')
    for n=1:length(ages)
        if outliers(n)==1
            disp(['    ' strnames{n}])
        end
    end
end
%% Plot data and results
figure; set(gcf, 'color', 'w')
subplot(5,5,[1 19]); hold on

% plot samples
top=length(ages)+1;
for n=1:length(ages)
    plot([ages(n)-errors(n),ages(n)+errors(n)],top-[n,n],'-k')
    caplength=0.2;
    plot([ages(n)-errors(n),ages(n)-errors(n)],top-[n-caplength/2,n+caplength/2],'-k')
    plot([ages(n)+errors(n),ages(n)+errors(n)],top-[n-caplength/2,n+caplength/2],'-k')
    plot(ages(n),top-n,'xk')
end

% plot one sigma limits
plot(xonesigma,xonesigma.*0-1,'-','Color',[0.5 0.5 0.5])

% plot averages
n=length(ages)+4;
meanvalue=Av; err=dAv;
plot([meanvalue-err,meanvalue+err],[n,n],'-g')
caplength=0.2;
plot([meanvalue-err,meanvalue-err],[n-caplength/2,n+caplength/2],'-g')
plot([meanvalue+err,meanvalue+err],[n-caplength/2,n+caplength/2],'-g')
plot(meanvalue,n,'xg')

n=length(ages)+3;
meanvalue=Wa; err=dWa;
plot([meanvalue-err,meanvalue+err],[n,n],'-r')
caplength=0.2;
plot([meanvalue-err,meanvalue-err],[n-caplength/2,n+caplength/2],'-r')
plot([meanvalue+err,meanvalue+err],[n-caplength/2,n+caplength/2],'-r')
plot(meanvalue,n,'xr')

n=length(ages)+2;
meanvalue=Bgf; err=dBgf;
plot([meanvalue-err,meanvalue+err],[n,n],'-b')
caplength=0.2;
plot([meanvalue-err,meanvalue-err],[n-caplength/2,n+caplength/2],'-b')
plot([meanvalue+err,meanvalue+err],[n-caplength/2,n+caplength/2],'-b')
plot(meanvalue,n,'xb')

xlim([minx maxx])
ylim([-2 length(ages)+5])
set(gca,'ytick',[])
set(gca,'XAxisLocation','Top');

% plot labels
subplot(5,5,[5 20]); hold on

top=length(ages)+1;
for n=1:length(ages)
    text(0,top-n,strnames{n},'Color','k')
end

n=length(ages)+4;
text(0,n,strAv,'Color','g')
n=length(ages)+3;
text(0,n,strWa,'Color','r')
n=length(ages)+2;
text(0,n,strBgf,'Color','b')
n=-1;
text(0,n,strrange,'Color',[0.5 0.5 0.5])

ylim([-2 length(ages)+5])
axis off

% plot ditributions
subplot(5,5,[21 24]); hold on
plot(x,y,'-k')
plot(x,aBgf*normalprob( x, Bgf, dBgf),'-b')
set(gca,'ytick',[])
xlim([minx maxx])

% plot credits
subplot(5,5,25); hold on
text(0,1,['CEAA v.' scriptversion])
text(0,0.5,'Angel Rodes, 2020')
text(0,0,'www.angelrodes.com')
ylim([0 max(1.3,(length(ages)+5)/4/3)])
axis off

%% Bye
disp(' ')
input ('Done. Press any key to exit.') % avoid closing everything if the script is run as "octave start"
end
