%Kai Kindred
%Glasgow Uni Aerospace MSc Project
%Converts Keplerian to Cartesian Coordinates and Vice-Versa

%Function Converting from Keplerian Elements to Cartesian Coordinates
function [r,v]=kep2car(kep,mu)
    
    a=kep(1);e=kep(2);inc=kep(3);omega=kep(4);w=kep(5);theta=kep(6);
    
    %intermediate perifocal plane parameter calculations
    p=a*(1-e^2); %distance from focus to perigee
    r_peri=p/(1+e*cosd(theta)); %distance of orbiting body to focus
    rx=r_peri*cosd(theta); %x-component of orbiting distance
    ry=r_peri*sind(theta); %y-component of orbiting distance
    theta_dot=sqrt(mu*p)/r_peri^2; %angular velocity of orbiting body in plane
    r_dot=sqrt(mu/p)*e*sind(theta); %tangential velocity of orbiting body
    vx=r_dot*cosd(theta)-r_peri*sind(theta)*theta_dot; %x-component of velocity
    vy=r_dot*sind(theta)+r_peri*cosd(theta)*theta_dot; %y-component of velocity
    
    %rotation matrix calculation
    R11=cosd(w)*cosd(omega)-sind(w)*cosd(inc)*sind(omega);
    R21=cosd(w)*sind(omega)+sind(w)*cosd(inc)*cosd(omega);
    R31=sind(w)*sind(inc);
    R12=-sind(w)*cosd(omega)-cosd(w)*cosd(inc)*sind(omega);
    R22=-sind(w)*sind(omega)+cosd(w)*cosd(inc)*cosd(omega);
    R32=cosd(w)*sind(inc);
    R13=sind(inc)*sind(omega);
    R23=-sind(inc)*cosd(omega);
    R33=cosd(inc);
    
    R = [R11 R12 R13; R21 R22 R23; R31 R32 R33];
    
    %rotating perifocal plane values using above matrix to get position and
    %velocity in cartesian coordinates
    r=[R11*rx+R12*ry; R21*rx+R22*ry; R31*rx+R32*ry];
    v=[R11*vx+R12*vy; R21*vx+R22*vy; R31*vx+R32*vy];
    %r=R*[rx;ry;0];
    %v=R*[vx;vy;0];
    
end