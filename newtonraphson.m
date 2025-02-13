%Kai Kindred
%Glasgow Uni Aerospace MSc Project
%Newton-Raphson Method
%Used to solve for eccentric anomaly from mean anomaly and eccentrictiy
%Takes in an initial guess and desired error tolerance

function eccentric_anomaly=newtonraphson(M,initial_guess,e,error)

    E=initial_guess;
    deltaE=1;

    while abs(deltaE)>error
        func=E-e*sin(E)-M;
        funcprime=1-e*cos(E);
        deltaE=-func/funcprime;
        E=E+deltaE;
    end

    eccentric_anomaly=E;

end