function ad_ = adj_(X)
Theta = X(1:3);
p = X(4:6);
ad_ = [crossm(Theta) zeros(3,3);
    crossm(p) crossm(Theta)];
