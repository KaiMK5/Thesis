
%Kai Kindred
%Glasgow Uni Aerospace MSc Project
%Magnetic Index Data Reader (Ap index from
%https://omniweb.gsfc.nasa.gov/form/dx1.html)
%Needs to be run before orbital model for NRLMSISE-00 atm model input
%mag idx data (hourly) needs to be at least 59 hours before epoch

%opens files and loads text
[file,path]=uigetfile('.txt');
fileID=fopen(strcat(path,file)); %can open files outside current directory
text=fscanf(fileID,'%c');

year_text=[];
dayofyear_text=[];
hour_text=[];
ap_text=[];

%parses text
i=1;
while i<length(text)

   year_text=[year_text; text(i:i+3)];
   dayofyear_text=[dayofyear_text; text(i+5:i+7)];
   hour_text=[hour_text; text(i+9:i+10)];
   ap_text=[ap_text; text(i+13:i+14)];
   
   i=i+16;
   
end

%converts text to numbers
magYear=str2num(year_text);
magDoY=str2num(dayofyear_text);
magHour=str2num(hour_text);
apIDX=str2num(ap_text);

%finds daily Ap index values by averaging hourly valus
dailyAP=[];
j=1;
while j<length(apIDX)
    dailyAP=[dailyAP; round(mean(apIDX(j:j+23)))];
    j=j+24;
end

%find 3-hour values for current time by averaging past 3 hours
avg3hoursAPcurrent=[];
for k=4:length(apIDX)
    avg3hoursAPcurrent=[avg3hoursAPcurrent; mean(apIDX(k-3:k-1))];
end

%find 3-hour values for 3 hours before current time by averaging past 3-6 hours
avg3hoursAP3before=[];
for l=7:length(apIDX)
    avg3hoursAP3before=[avg3hoursAP3before; mean(apIDX(l-6:l-4))];
end

%find 3-hour values for 6 hours before current time by averaging past 6-9 hours
avg3hoursAP6before=[];
for m=10:length(apIDX)
    avg3hoursAP6before=[avg3hoursAP6before; mean(apIDX(m-9:m-7))];
end

%find 3-hour values for 9 hours before current time by averaging past 9-12 hours
avg3hoursAP9before=[];
for n=13:length(apIDX)
    avg3hoursAP9before=[avg3hoursAP9before; mean(apIDX(n-12:n-10))];
end

%find average of 8 3-hour values for 12+ hours before current time by 
%averaging past 3 hour increments
avg3hoursAP12before=[];
for o=37:length(apIDX)
    i1=mean(apIDX(o-15:o-13));
    i2=mean(apIDX(o-18:o-16));
    i3=mean(apIDX(o-21:o-19));
    i4=mean(apIDX(o-24:o-22));
    i5=mean(apIDX(o-27:o-25));
    i6=mean(apIDX(o-30:o-28));
    i7=mean(apIDX(o-33:o-31));
    i8=mean(apIDX(o-36:o-34));
    avg3hoursAP12before=[avg3hoursAP12before; mean([i1,i2,i3,i4,i5,i6,i7,i8])];
end

%find average of 8 3-hour values for 36+ hours before current time by 
%averaging past 3 hour increments
avg3hoursAP24before=[];
for p=60:length(apIDX)
    id1=mean(apIDX(p-38:p-36));
    id2=mean(apIDX(p-41:p-39));
    id3=mean(apIDX(p-44:p-42));
    id4=mean(apIDX(p-47:p-45));
    id5=mean(apIDX(p-50:p-48));
    id6=mean(apIDX(p-53:p-51));
    id7=mean(apIDX(p-56:p-54));
    id8=mean(apIDX(p-59:p-57));
    avg3hoursAP24before=[avg3hoursAP24before; mean([id1,id2,id3,id4,id5,id6,id7,id8])];
end

fclose('all');

clearvars -except dailyAP avg3hoursAPcurrent avg3hoursAP3before avg3hoursAP6before avg3hoursAP9before avg3hoursAP12before avg3hoursAP24before magYear magDoY magHour apIDX solarYear solarDoY f107 avg81daysF107