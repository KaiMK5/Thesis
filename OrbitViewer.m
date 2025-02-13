%Kai Kindred
%Glasgow Uni Aerospace MSc Project
%Orbit Viewer

%run this after running the primary orbit model to visualize an orbit akin to STK or GMAT 
%ensure all variables are still in the workspace
%can add indexing to times and position if only want to visualize parts of
%orbit
%%%REQUIRES R2021a OR LATER%%%

timeliststr=datestr(times);
satdata=timeseries(1000*position,timeliststr); %uses meters as position
sc = satelliteScenario(epoch_start,epoch_end,stepsize);
sat = satellite(sc, satdata,'CoordinateFrame',"inertial");
satelliteScenarioViewer(sc)
