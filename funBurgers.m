function yp = funBurgers(t,y)

global epsilon L1 L2
yp = epsilon*L1*y - (1/2)*L2*(y.^2);

end