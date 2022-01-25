function h = h_sys(g,V)

y1 = vedge_inv(logSE3(g));
y2 = V;

% h = [y1(1:3);g(1:3,4);y2];
h = [y1;y2];
end