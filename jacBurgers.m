function J = jacBurgers(t,y)

global epsilon L1 L2 M
J = epsilon*L1-L2.*repmat(y,1,M)';

end