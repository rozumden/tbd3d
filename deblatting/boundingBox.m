function bbox = boundingBox(bw)
% returns the bounding box of nonzero pixels in bw as rectangle specified as [top bottom left right]

p1 = any(bw,2); p2 = any(bw,1); % projections
p11 = find(p1, 1, 'first'); p12 = find(p1, 1, 'last');
p21 = find(p2, 1, 'first'); p22 = find(p2, 1, 'last');
bbox = [p11 p12 p21 p22];
end