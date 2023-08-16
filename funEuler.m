function yp = funEuler(t,y)

yp(1) = -2*y(2)*y(3);
yp(2) = 5/4*y(3)*y(1);
yp(3) = -1/2*y(1)*y(2);
yp = [yp(1);yp(2);yp(3)];

end