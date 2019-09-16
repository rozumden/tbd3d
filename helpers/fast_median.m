function B = fast_median(x,y,z)
xy = x > y;
xz = x > z;
yz = y > z;

y_ind = (xy & yz) | (~xy & ~yz);
x_ind = (~xy & yz & xz) | (xy & ~yz & ~xz);
B = x.*x_ind + y.*y_ind + z.*(~x_ind & ~y_ind);