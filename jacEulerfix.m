function J = jacEulerfix()

global y0
J = [0 -2*y0(3) -2*y0(2);
    5/4*y0(3) 0 5/4*y0(1);
    -1/2*y0(2) -1/2*y0(1) 0];

end