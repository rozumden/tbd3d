function bb = bbextend(BoundingBox, ex, Size)
bb = round(BoundingBox) - [ex ex -2*ex -2*ex];
bb = round(bb);
bb(3:4) = bb(1:2) + bb(3:4)-1;
bb(bb<1) = 1;
if bb(3) > Size(2), bb(3) = Size(2); end
if bb(4) > Size(1), bb(4) = Size(1); end
