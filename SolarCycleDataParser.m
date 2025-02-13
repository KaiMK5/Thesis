%Kai Kindred
%Glasgow Uni Aerospace MSc Project
%Solar Index Data Reader (F10.7 index from
%https://omniweb.gsfc.nasa.gov/form/dx1.html)
%for NRLMSISE-00 atmospheric model input
%solar data (daily) needs to be 40 days before epoch start in orbit model
%...and 40 days after epoch end

%opens file and imports text
[file,path]=uigetfile('.txt');
fileID=fopen(strcat(path,file)); %can open files outside current directory
text=fscanf(fileID,'%c');

year_text=[];
dayofyear_text=[];
hour_text=[];
f107_text=[];

%parses text
i=1;
while i<length(text)

   year_text=[year_text; text(i:i+3)];
   dayofyear_text=[dayofyear_text; text(i+5:i+7)];
   hour=[hour_text; text(i+9:i+10)];
   f107_text=[f107_text; text(i+12:i+16)];
   
   i=i+18;
   
end

%convert text to numbers
solarYear=str2num(year_text);
solarDoY=str2num(dayofyear_text);
f107=str2num(f107_text);
%smooths out extraneously large F10.7 values
%could not interpolate readily as they'd appear in a series of points
for m=1:length(f107)
    if f107(m)==999.9
        f107(m)=f107(m-1);
    end
end

%Finds centered 81-day average
avg81daysF107=[];
for j=41:length(f107)-40
    avg81daysF107=[avg81daysF107; mean(f107(j-40:j+40))];
end

fclose('all');

date=[];
for k=1:length(dayofyear_text)
    date=[date;strcat(year_text(k,:),'-',dayofyear_text(k,:))];
end

clearvars -except solarYear solarDoY f107 avg81daysF107 dailyAP avg3hoursAPcurrent avg3hoursAP3before avg3hoursAP6before avg3hoursAP9before avg3hoursAP12before avg3hoursAP24before magYear magDoY magHour apIDX 
