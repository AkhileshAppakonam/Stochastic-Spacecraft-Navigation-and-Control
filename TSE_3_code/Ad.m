% Ad function 
function Ad_ = Ad(g)
R = g(1:3,1:3); % rotation matrix
r = g(4,1:3); % translational displacement
Ad_ = [R        zeros(3,3);
	crossm(r)*R     R];
