function [res] = reconstruct(matF,vis)
if ~exist('vis','var')
	vis = false;
end
if vis
	hold on;
end

mag_th = 0.1;
dist_th = 1.5;
grad_th = 0.00;

M = double(diskMask(size2(matF),0.99*size(matF,1)/2));

K = zeros(3); K(3,3) = 1;
K(3,1) = size(matF,1)/2; K(3,2) = K(3,1);
flen = K(3,1); 
K(1,1) = flen; K(2,2) = K(1,1);

KP = cameraParameters('IntrinsicMatrix', K);
% iv = [251:390];
iv = [1:(size(matF,4)-1)];
vSet = [];

pnts = cell(1, iv(end));
Rs = cell(1, iv(end));
Cs = cell(1, iv(end));
indexPairs = cell(1, iv(end));
res = [];
for kk = iv
	fprintf('Processing view %d\n', kk);
	I1 = matF(:,:,:,kk); I2 = matF(:,:,:,kk+1);
	flow = TD.estimate_flow(I1,I2,mag_th,vis);
	[pnts1, pnts2] = TD.estimate_corr(flow, M, I1, I2, mag_th, grad_th);
	
	fprintf(' - number of points %d\n', size(pnts2,1));

	R = [];
	try
		[R, C, inlierIdx] = helperEstimateRelativePose(pnts1, pnts2, KP);
	catch
		fprintf('[ERROR] View %d\n', kk);
	end

	if ~isempty(R)
		Rs{kk} = R;
		Cs{kk} = C;
		pnts1 = pnts1(inlierIdx,:);
		pnts2 = pnts2(inlierIdx,:);
	else
		Rs{kk} = Rs{kk-1};
		Cs{kk} = Cs{kk-1};
	end

	
	if ~isempty(pnts{kk})
		Vx = interp2(flow.Vx, pnts{kk}(:,1), pnts{kk}(:,2));
		Vy = interp2(flow.Vy, pnts{kk}(:,1), pnts{kk}(:,2));
		pntsTracked = pnts{kk} + [Vx Vy];
		mag = sqrt(Vx.^2 + Vy.^2);
		[~,ind]	= find(sum(~isnan(pntsTracked)') == 2 & mag' > mag_th);
		pntsTracked = pntsTracked(ind,:);
		indexPairs{kk} = [ind; [1:size(pntsTracked)]]';
		pnts{kk+1} = pntsTracked;
		[d,~] = pdist2(pnts{kk}, pnts1, 'euclidean' ,'smallest', 1);
		inl = d > dist_th;
		pnts1 = pnts1(inl, :);
		pnts2 = pnts2(inl, :);
	end

	s1 = size(pnts{kk},1);
	s2 = size(pnts{kk+1},1);

	indexPairs{kk} = [ indexPairs{kk}; [[(s1+1):(s1+size(pnts1,1))]; [(s2+1):(s2+size(pnts2,1))]]' ];
	pnts{kk} = [pnts{kk}; pnts1];
	pnts{kk+1} = [pnts{kk+1}; pnts2];

	p1 = pnts{kk}(indexPairs{kk}(:,1),:);
	p2 = pnts{kk+1}(indexPairs{kk}(:,2),:);
	[Rs{kk}, Cs{kk}, inlierIdx] = helperEstimateRelativePose(p1, p2, KP);
	
	res(kk).flow = flow;
	res(kk).img = I1;
	res(kk+1).img = I2;

	if vis, hold on; plot(pnts{kk}(:,1), pnts{kk}(:,2),'.r'); plot(pnts{kk+1}(:,1), pnts{kk+1}(:,2),'.g'); drawnow; end
end

vSet = viewSet;
vSet = addView(vSet, 1, 'Points', pnts{1}, 'Orientation', eye(3), 'Location', zeros(1, 3));
for kk = iv
	ik = kk + 1;
	vSet = addView(vSet, ik, 'Points', pnts{ik});
	vSet = addConnection(vSet, kk, ik, 'Matches', indexPairs{kk});
	res(kk).pnts = pnts{kk};
	res(kk).indexPairs = indexPairs{kk};

	prevPose = poses(vSet, kk);
    prevOrientation = prevPose.Orientation{1};
    prevLocation    = prevPose.Location{1};
    orientation = Rs{kk} * prevOrientation;
    location    = prevLocation + Cs{kk} * Rs{kk};

	vSet = updateView(vSet, ik, 'Orientation', orientation, 'Location', location);
	tracks = findTracks(vSet);
	camPoses = poses(vSet);
	xyzPoints = triangulateMultiview(tracks, camPoses, KP);
	[xyzPoints, camPoses, reprojectionErrors] = bundleAdjustment(xyzPoints, ...
		tracks, camPoses, KP, 'FixedViewId', 1, ...
		'PointsUndistorted', true);
	vSet = updateView(vSet, camPoses);
end

norms = sqrt(sum(xyzPoints'.^2,1))';
goodIdx = (reprojectionErrors < 2.5 & norms < size(matF,1)*1.3);

if vis
	figure;
	pcshow(xyzPoints(goodIdx, :), 'VerticalAxis', 'y', 'VerticalAxisDir', 'down','MarkerSize', 45);
	hold on;
	camPoses = poses(vSet);
	% plotCamera(camPoses, 'Size', 0.2);

	% model = []; [model.Center,model.Radius] = sphereFit(xyzPoints(goodIdx, :)); 
	pcl = pointCloud(xyzPoints(goodIdx, :)); model = pcfitsphere(pcl,15,'Confidence',99.999,'MaxNumTrials',1e5);	
	[x,y,z] = sphere; plot3(x*model.Radius+model.Center(1), y*model.Radius+model.Center(2), z*model.Radius+model.Center(3));

	v1 = camPoses.Location{1} - model.Center;
	v2 = camPoses.Location{end} - model.Center;
	ang = atan2(norm(cross(v1,v2)), dot(v1,v2));

	keyboard
end

