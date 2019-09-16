function bwe = bwedge(bw)
B = bwboundaries(bw);
bwe = zeros(size(bw)); bwe(sub2ind(size(bw),B{1}(:,1),B{1}(:,2))) = 1;