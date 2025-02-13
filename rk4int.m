%Kai Kindred
%Glasgow Uni Aerospace MSc Project
%Fourth Order Runge-Kutta Numerical Integration Function

function x = rk4int(modelname, h, x) %takes in step size and state variables

k1 = h*feval(modelname, x);         % evaluate derivative k1
k2 = h*feval(modelname, x+k1/2);	% evaluate derivative k2
k3 = h*feval(modelname, x+k2/2);	% evaluate derivative k3
k4 = h*feval(modelname, x+k3);		% evaluate derivative k4

x = x + (k1 + 2*k2 + 2*k3 + k4)/6;  % averaged output
