%Kai Kindred
%Glasgow Uni Aerospace MSc Project
%TLE File Reader and Plotter (using txt files from celestrak.com)
%Must run after solar and magnetic index parsers for plots at end

%Read File
[file,path]=uigetfile('.txt');
fileID=fopen(strcat(path,file)); %can open files outside current directory
text=fscanf(fileID,'%c');

epochYear_text=[];
epochDay_text=[];
inclination_text=[];
RAAN_text=[];
eccentricity_text=[];
argPerigee_text=[];
meanAnomaly_text=[];
meanMotion_text=[];
revNumber_text=[];

%parse file
i=1;
while i<length(text)
   epochYear_text=[epochYear_text; text(i+18:i+19)];
   epochDay_text=[epochDay_text; text(i+20:i+31)];
   inclination_text=[inclination_text; text(i+70+9:i+70+16)];
   RAAN_text=[RAAN_text; text(i+70+18:i+70+25)];
   eccentricity_text=[eccentricity_text; text(i+70+27:i+70+33)];
   argPerigee_text=[argPerigee_text; text(i+70+35:i+70+42)];
   meanAnomaly_text=[meanAnomaly_text; text(i+70+44:i+70+51)];
   meanMotion_text=[meanMotion_text; text(i+70+53:i+70+63)];
   revNumber_text=[revNumber_text; text(i+70+64:i+70+68)];
   
   i=i+71+70+1;
   
end

%convert yy to yyyy
epochYear=str2num(epochYear_text);
for j=1:length(epochYear)
    if epochYear(j)<50 %basically assuming all dates fall between 1950-2049
        epochYear(j)=epochYear(j)+2000;
    else
        epochYear(j)=epochYear(j)+1900;
    end
end
%break down day as a decimal in hours, minutes, and seconds
epochDay=str2num(epochDay_text); %grabs everything before decimal
epochTime=epochDay-floor(epochDay); %grabs everything after decimal
epochHour=floor(epochTime*24); %hour
epochTime2=epochTime*24-epochHour;
epochMinute=floor(epochTime2*60); %m
epochTime3=epochTime2*60-epochMinute;
epochSecond=floor(epochTime3*60); %s

%construct epochs
filler=zeros(length(epochYear),1);
epochin=datetime([epochYear filler filler filler filler filler])+calmonths(1)+days(1); %for some reason was setting as 30 Nov 1979 instead of 1/1/1980 so adjusted
epoch=epochin+days(floor(epochDay))+hours(epochHour)+minutes(epochMinute)+seconds(epochSecond);

%convert other values to numbers

%plugs holes in RAAN data, sometimes observed in TLE files if conversion
%fails
inclination=str2num(inclination_text); %deg
RAAN=str2num(RAAN_text);    %deg
if isempty(RAAN)==1
    for m=1:length(RAAN_text) 
        RAAN=[RAAN;str2num(RAAN_text(m,:))]; 
    end
end
%adds decimal point to eccentricity
eccentricity_text=strcat('.',eccentricity_text(:,:));
eccentricity=str2num(eccentricity_text);
argPerigee=str2num(argPerigee_text); %deg
meanAnomaly=deg2rad(str2num(meanAnomaly_text)); %rad
meanMotion=str2num(meanMotion_text);
revNumber=str2num(revNumber_text); %#

%Plot Orbit akin to model
mu=3.986E5;         %km^3/s^2
RE=6378.2;          %km
orbitDuration=seconds(epoch(end)-epoch(1));
error=0.001; %for Newton-Raphson Method
trueAnomaly=[];
semiMajorAxis=[];

figure
E=wgs84Ellipsoid("kilometer");
ellipsoid(0,0,0,E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis)
hold on

for k=1:length(epoch)
    
    %compute semi-major axis at each epoch
    Period=3600*24/meanMotion(k); %s
    a = (1/1000)*((Period/(2*pi))^2*3.986e14)^(1/3); %km
    semiMajorAxis=[semiMajorAxis;a]; %km
    
    %Newton-Raphson method, with initial guesses, to find eccentric anomaly
    if eccentricity(k)<=0.2 || meanAnomaly(k)==pi
        initial_guess=meanAnomaly(k);
    elseif meanAnomaly(k)<pi
        initial_guess=meanAnomaly(k)+eccentricity(k)/2;
    elseif meanAnomaly(k)>pi
        initial_guess=meanAnomaly(k)-eccentricity(k)/2;
    end
    eccentric_anomaly=newtonraphson(meanAnomaly(k),initial_guess,eccentricity(k),error); %rad
    theta=2*atand(tan(eccentric_anomaly/2)/sqrt((1-eccentricity(k))/(1+eccentricity(k)))); %deg
    trueAnomaly=[trueAnomaly;theta]; %deg
    
    %convert Keplerian elements to cartesian coordinates and plot
    kep=[semiMajorAxis(k), eccentricity(k), inclination(k), RAAN(k), argPerigee(k), trueAnomaly(k)];
    [r,v]=kep2car(kep,mu);
    plot3(r(1),r(2),r(3),'r*');
    hold on 
    
end

axis equal
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title(file)

fclose('all');

%Comparison Calculations

%altitude check to ensure not larger than 1000km & thus invalid for atm model
radius=semiMajorAxis.*(1-eccentricity.^2)./(1+eccentricity.*cosd(trueAnomaly)); %km
height=radius-RE;   %km
if height(:)<1000
    disp('Good Altitudes')
end

%Finds Initial Keplerian Elements to Input into Orbital Model
comparisonstartIDX=1; %7500 for NOAA-12, allows to only compare parts of orbit
fprintf('Epoch Start: %s\n', epoch(comparisonstartIDX))
fprintf('Epoch End: %s\n', epoch(end))
fprintf('Starting Semi-Major Axis: %d km\n',semiMajorAxis(comparisonstartIDX))
fprintf('Starting Eccentricity: %d\n',eccentricity(comparisonstartIDX))
fprintf('Starting True Anomaly: %d deg\n',trueAnomaly(comparisonstartIDX))
fprintf('Inclination: %d deg\n',inclination(comparisonstartIDX))
fprintf('Starting RAAN: %d deg\n',RAAN(comparisonstartIDX))
fprintf('Starting Argument of Perigee: %d deg\n',argPerigee(comparisonstartIDX))

%Finds Differences in Keplerian elements at each step in TLE data
epochdiff=[];
semimajaxdiff=[];
eccentricitydiff=[];
RAANdiff=[];
argperigeediff=[];
for n=comparisonstartIDX:length(epoch)-1
    epochdiff=[epochdiff;epoch(n+1)-epoch(n)];
    semimajaxdiff=[semimajaxdiff;semiMajorAxis(n+1)-semiMajorAxis(n)];
    eccentricitydiff=[eccentricitydiff;eccentricity(n+1)-eccentricity(n)];
    RAANdiff=[RAANdiff;RAAN(n+1)-RAAN(n)];
    argperigeediff=[argperigeediff;argPerigee(n+1)-argPerigee(n)];
end

%average time between TLE data points
avgEpochDiff=mean(seconds(epochdiff)); %s

%Finds average changes in RAAN and Argument of Perigee, if don't use
%oblateness in orbital model
avgRAANDiff=mean(RAANdiff); %deg
avgArgPerigeeDiff=mean(argperigeediff);
avgDeltaRAAN=avgRAANDiff/avgEpochDiff;
avgDeltaArgPerigee=avgArgPerigeeDiff/avgEpochDiff;
fprintf('Average RAAN Delta (deg/s): %d\n',avgDeltaRAAN)
fprintf('Average Argument of Perigee Delta (deg/s): %d\n',avgDeltaArgPerigee)

%Finds average changes in semi-major axis and eccentricity to compare with
%orbital model
avgSemiMajAxisDiff=mean(semimajaxdiff); %km
avgEccentricDiff=mean(eccentricitydiff);
avgDeltaSMA=avgSemiMajAxisDiff/avgEpochDiff; %km/s
avgDeltaEcc=avgEccentricDiff/avgEpochDiff; %/s
fprintf('Average Semi-Major Axis Delta (km/s): %d\n',avgDeltaSMA)
fprintf('Average Eccentricity Delta (/s): %d\n',avgDeltaEcc)

%Finds total differences in semi-major axis and eccentricity
epochduration=seconds(epoch(end)-epoch(comparisonstartIDX));
finalSMAdiff=semiMajorAxis(comparisonstartIDX)-semiMajorAxis(end);
finalECCdiff=eccentricity(comparisonstartIDX)-eccentricity(end);
finalSMAdelta=finalSMAdiff/epochduration;
finalECCdelta=finalECCdiff/epochduration;

%constructs epochs
filler3=zeros(length(solarYear),1);
epochin3=datetime([solarYear filler3 filler3 filler3 filler3 filler3])+calmonths(1)+days(1);
solarepoch=epochin3+days(solarDoY);

%Plots Semi-Major Axis and Eccentricity alongside solar flux
figure
subplot(2,1,1)
yyaxis left
plot(solarepoch(41:end-40),f107(41:end-40))
ylabel('Daily F107 Index')
yyaxis right
plot(epoch,semiMajorAxis)
ylabel('Semi-Major Axis (km)')
xlabel('Date')
title('Daily F107 Index vs. Semi-Major Axis')
subplot(2,1,2)
yyaxis left
plot(solarepoch(41:end-40),f107(41:end-40))
ylabel('Daily F107 Index')
yyaxis right
plot(epoch(1:end-1),semimajaxdiff)
ylabel('Semi-Major Axis Step Change (km)')
xlabel('Date')
title('Daily F107 Index vs. Semi-Major Axis Step Change')

figure
subplot(2,1,1)
yyaxis left
plot(solarepoch(41:end-40),f107(41:end-40))
ylabel('Daily F107 Index')
yyaxis right
plot(epoch,eccentricity)
ylabel('Eccentricity')
xlabel('Date')
title('Daily F107 Index vs. Eccentricity')
subplot(2,1,2)
yyaxis left
plot(solarepoch(41:end-40),f107(41:end-40))
ylabel('Daily F107 Index')
yyaxis right
plot(epoch(1:end-1),eccentricitydiff)
ylabel('Eccentricity Step Change')
xlabel('Date')
title('Daily F107 Index vs. Eccentricity Step Change')

rp=semiMajorAxis(end)*(1-eccentricity(end));
fprintf('Satellite still needs to drop by %d km to burn up\n',rp-100)
