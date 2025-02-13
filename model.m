%Kai Kindred
%Glasgow Uni Aerospace MSc Project
%State Space Model Function - finds derivatives of Keplerian elements
%Paired with integration method

function xdot = model(x) %takes in state variables

global atm_density v_magnitude delta true_anomaly mu RE	inclination	% global parameters transferred from main program


xdot(1,1) = -atm_density*(1000*x(1))^2*(v_magnitude*1000)^3*delta/(1000^3*mu)/1000;  %semi-major axis derivative
xdot(2,1) = -atm_density*1000*v_magnitude*delta*(x(2)+cosd(true_anomaly));         %eccentricity derivative
xdot(3,1) = -9.964/3600/24*(RE/x(1))^3.5*(1-x(2)^2)^(-2)*cosd(inclination);   %RAAN derivative
xdot(4,1) = 4.982/3600/24*(RE/x(1))^3.5*(1-x(2)^2)^(-2)*(5*(cosd(inclination))^2-1); %argument of Perigee derivative
