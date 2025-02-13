%Kai Kindred
%Glasgow Uni Aerospace MSc Project
%Orbit Model
%Needs solar flux data (40 days before and after timespan) in workspace
%Needs magnetic flux data (3 days before start up until end) in workspace

global atm_density v_magnitude delta true_anomaly mu RE inclination

%Space Vehicle Inputs - currently for DMSP satellite
mass=830; %kg
Cd=2.2;
S=3.7*1.2+9.29; %m^2 SV surface area

%Attacting Body Inputs (currently Earth)
mu=3.986E5;         %km^3/s^2
RE=6378.2;          %km
w=2*pi/24/3600;     %/s rotational vel atm

%Orbital Parameters/Keplerian Elements
initial_eccentricity=0.8;
inclination=99;      %degrees
initial_radius_pericenter=7228; %RE+500; %km
initial_semi_major_axis=initial_radius_pericenter/(1-initial_eccentricity); %km
initial_true_anomaly=95;    %degrees
initial_RAAN=175;            %degrees
initial_anomaly_of_pericenter=260; %degrees

error=0.001; %for Newton-Raphson Method

%Epoch Definition - can specify duration or epoch end
epoch_start=datetime('1996-02-10 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');
%epoch_duration=180; %days
%epoch_end=epoch_start+epoch_duration;
epoch_end=datetime('2009-11-21 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');
stepsize=3*60*60; %seconds

epoch_duration=epoch_end-epoch_start; %seconds
endtime=seconds(epoch_duration); %seconds

%if epoch is not divisible by stepsize, adjust epoch to next smallest divisible
if mod(endtime,stepsize)~=0
    modified_endtime=endtime-mod(endtime,stepsize);
    epoch_end=epoch_start+seconds(modified_endtime);
    endtime=modified_endtime;
end

initial_period = 2*pi*sqrt(initial_semi_major_axis^3/mu); %seconds

%if not at pericenter to start, find timeshift for iteration
time_from_pericenter=0;           %s
if initial_true_anomaly ~= 0
    initial_eccentric_anomaly=2*atan(sqrt((1-initial_eccentricity)/(1+initial_eccentricity))*tand(initial_true_anomaly/2));
    initial_mean_anomaly=initial_eccentric_anomaly-initial_eccentricity*sin(initial_eccentric_anomaly);
    time_from_pericenter=initial_mean_anomaly/sqrt(mu/initial_semi_major_axis^3);
end

semi_major_axis=initial_semi_major_axis;
eccentricity=initial_eccentricity;
RAAN=initial_RAAN;
anomaly_of_pericenter=initial_anomaly_of_pericenter;

%Plot Center of Attracting Body
figure
E=wgs84Ellipsoid("kilometer");
ellipsoid(0,0,0,E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis)
hold on

%time iterate through orbit
derivatives=[];
plottingdata=[];
plottingepoch=[];
times=[];
position=[];
for i=0:stepsize:endtime

    %Determine ECI Position and Velocity @ Point in Time
    
    %Newton-Raphson method w/ intial guess for eccentric anomaly
    mean_anomaly=(i+time_from_pericenter)*sqrt(mu/semi_major_axis^3); %uses offset time to properly place satellite in its orbit
    if eccentricity<=0.2 || mean_anomaly==pi
        initial_guess=mean_anomaly;
    elseif mean_anomaly<pi
        initial_guess=mean_anomaly+eccentricity/2;
    elseif mean_anomaly>pi
        initial_guess=mean_anomaly-eccentricity/2;
    end
    eccentric_anomaly=newtonraphson(mean_anomaly,initial_guess,eccentricity,error);
    %Current Epoch's true anomaly
    true_anomaly=2*atand(tan(eccentric_anomaly/2)/sqrt((1-eccentricity)/(1+eccentricity)));
    %COnvert to cartesian coordinates
    kep=[semi_major_axis, eccentricity, inclination, RAAN, anomaly_of_pericenter, true_anomaly];
    radius=semi_major_axis*(1-eccentricity^2)/(1+eccentricity*cosd(true_anomaly)); %km
    height=abs(radius)-RE; %not used for anything, not really accurate either
    [r,v]=kep2car(kep,mu); %km
    position=[position;r']; %for orbit visualization
    
    %Plot Current Point (every 3 hours to make figure less cumbersome)
    if mod(i,3*3600)==0
        plot3(r(1),r(2),r(3),'r*');
        hold on
    end
    
    %Determine Current Epoch in UTC Time
    epoch=epoch_start+seconds(i);
    times=[times;epoch];
    epoch_year=year(epoch);
    epoch_month=month(epoch);
    epoch_day=day(epoch);
    epoch_hour=hour(epoch);
    epoch_minute=minute(epoch);
    epoch_second=second(epoch);
    utc=[epoch_year epoch_month epoch_day epoch_hour epoch_minute epoch_second];
    
    %Determine Geodetic Position of Space Vehicle
    lla = eci2lla(1000*r',utc); %take in position as meters
    altitude=abs(lla(3)); %m
    
    %Find Atmospheric Density 
    dayOfYear = day(epoch,'dayofyear'); %finds day of year (ex: #5/365) of current epoch for atm model fxn input
    UTseconds=epoch_hour*3600+epoch_minute*60+epoch_second; %finds number of seconds in day of current epoch for atm model fxn input
    if altitude<=100*1e3 %m %Karman line defined as reentry
        disp(epoch)
        final_duration=epoch-epoch_start;
        disp(final_duration)
        fprintf('Satellite below Karman Line\n')
        break
    elseif altitude<=1e6 %atm model only valid through 1,000,000 meters
        %incorporate solar and magnetic index data
        for k=1:length(avg81daysF107)
            if solarYear(k+40)==epoch_year && solarDoY(k+40)==dayOfYear
                dailysolar=f107(k+40);
                avgsolar=avg81daysF107(k);
            end
        end
        for j=60:length(magHour)
            if magYear(j)==epoch_year && magDoY(j)==dayOfYear && magHour(j)==epoch_hour
                mag=[dailyAP(floor(j/24)) avg3hoursAPcurrent(j-3) avg3hoursAP3before(j-6)...
                    avg3hoursAP6before(j-9) avg3hoursAP9before(j-12) ...
                    avg3hoursAP12before(j-36) avg3hoursAP24before(j-59)];
            end
        end
        [T rho]=atmosnrlmsise00(altitude,lla(1),lla(2),epoch_year,dayOfYear,UTseconds,avgsolar, dailysolar, mag); %takes in meters and degrees 
        atm_density=rho(6);
    else
        atm_density=0; %total atm mass density assumed zero above 1000km
    end
    
    %Calculate change in Keplerian elements
    v_magnitude=sqrt(v(1)^2+v(2)^2+v(3)^2); %km/s
    F=(1-radius*w*cosd(inclination)/v_magnitude)^2; %dimensionless
    delta=F*S*Cd/mass; %m^2/kg
    x=[semi_major_axis; eccentricity; RAAN; anomaly_of_pericenter]; % 
    xdot = model(x);
    derivatives=[derivatives;xdot(1:2,:)'];
    %comparison data only every day - for investigation, change when
    %comparing to TLEs
    if mod(i,86400)==0 || i==0
        plottingdata=[plottingdata; x(1),xdot(1),x(2),xdot(2)];
        plottingepoch=[plottingepoch;epoch];
    end
    x = rk4int('model', stepsize, x);
    semi_major_axis=x(1);
    if x(2)>=0 && initial_eccentricity>0 %prevents eccentricity from <0 and from changing when in circular orbit
        eccentricity=x(2);
    else
        eccentricity=0;
    end
    RAAN=x(3);
    anomaly_of_pericenter=x(4);
    
end
rp=semi_major_axis-eccentricity*semi_major_axis;

axis equal
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title(['Eccentricity: ',num2str(initial_eccentricity),', Inclination: ',num2str(inclination),', Semi-Major Axis: ',num2str(initial_semi_major_axis),', RAAN: ',num2str(initial_RAAN),', Anomaly of Pericenter: ', num2str(initial_anomaly_of_pericenter)])

%Comparison Calculations
avgDeltaSMA=mean(derivatives(:,1)); %km/s
avgDeltaEcc=mean(derivatives(:,2)); %/s
fprintf('Average Semi-Major Axis Delta (km/s): %d\n',avgDeltaSMA)
fprintf('Average Eccentricity Delta (/s): %d\n',avgDeltaEcc)

%Manually Paste TLE avg semi-major axis and eccentricty rates for errors
% deltaSMAerror=100*(avgDeltaSMA--2.331609e-08)/-2.331609e-08;
% deltaECCerror=100*(avgDeltaEcc--2.571204e-13)/-2.571204e-13;
% fprintf('Semi-Major Axis Error: %d\n',deltaSMAerror)
% fprintf('Eccentricity Error: %d\n',deltaECCerror)

finalSMAdiff=initial_semi_major_axis-semi_major_axis;
finalECCdiff=initial_eccentricity-eccentricity;
fprintf('Total Semi-Major Axis Change (km): %d\n',finalSMAdiff)
fprintf('Total Eccentricity Change: %d\n',finalECCdiff)
finalSMAdelta=finalSMAdiff/endtime;
finalECCdelta=finalECCdiff/endtime;

%Construct epoch for solar data
filler2=zeros(length(solarYear),1);
epochin2=datetime([solarYear filler2 filler2 filler2 filler2 filler2])+calmonths(1)+days(1);
solarepoch=epochin2+days(solarDoY);

%Plot semi-major axis and eccentricity with solar flux, & find correlations
figure
subplot(2,1,1)
yyaxis left
plot(solarepoch(41:end-40),f107(41:end-40))
ylabel('Daily F107 Index')
yyaxis right
plot(plottingepoch,plottingdata(:,1))
ylabel('Semi-Major Axis (km)')
xlabel('Date')
title('Daily F107 Index vs. Semi-Major Axis')
subplot(2,1,2)
yyaxis left
plot(solarepoch(41:end-40),f107(41:end-40))
ylabel('Daily F107 Index')
yyaxis right
plot(plottingepoch,plottingdata(:,2))
ylabel('Semi-Major Axis Derivative (km/s)')
xlabel('Date')
title('Daily F107 Index vs. Semi-Major Axis Derivative')

SMAcorr=corrcoef(f107(41:end-40),plottingdata(:,1));
adotcorr=corrcoef(f107(41:end-40),plottingdata(:,2));
fprintf('F107 and SMA Correlation:  %d\n',SMAcorr(2))
fprintf('F107 and SMA Derivative Correlation:  %d\n',adotcorr(2))

figure
subplot(2,1,1)
yyaxis left
plot(solarepoch(41:end-40),f107(41:end-40))
ylabel('Daily F107 Index')
yyaxis right
plot(plottingepoch,plottingdata(:,3))
ylabel('Eccentricity')
xlabel('Date')
title('Daily F107 Index vs. Eccentricity')
subplot(2,1,2)
yyaxis left
plot(solarepoch(41:end-40),f107(41:end-40))
ylabel('Daily F107 Index')
yyaxis right
plot(plottingepoch,plottingdata(:,4))
ylabel('Eccentricity Derivative (/s)')
xlabel('Date')
title('Daily F107 Index vs. Eccentricity Derivative')

ECCcorr=corrcoef(f107(41:end-40),plottingdata(:,3));
edotcorr=corrcoef(f107(41:end-40),plottingdata(:,4));
fprintf('F107 and Eccentricity Correlation:  %d\n',ECCcorr(2))
fprintf('F107 and Eccentricity Derivative Correlation:  %d\n',edotcorr(2))

fprintf('Satellite still needs to drop by %d km to burn up\n',rp-100)
