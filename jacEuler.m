function J = jacEuler(t,y)

J = [0 -2*y(3) -2*y(2);
    5/4*y(3) 0 5/4*y(1);
    -1/2*y(2) -1/2*y(1) 0];

end