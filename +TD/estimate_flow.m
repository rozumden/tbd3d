function [flow] = estimate_flow(I1,I2,mag_th,vis)
if ~exist('vis','var')
	vis = false;
end

opticFlow = opticalFlowFarneback;
flow = estimateFlow(opticFlow,rgb2gray(I1));
flow = estimateFlow(opticFlow,rgb2gray(I2));

vx = mean(flow.Vx (flow.Magnitude > mag_th));
vy = mean(flow.Vy (flow.Magnitude > mag_th));

if vis
	clf;
	imshow(I1); hold on;
	plot(flow,'DecimationFactor',[1 1],'ScaleFactor',1);
	set(gca,'Unit','normalized','Position',[0 0 1 1]);
	quiver(14,14,vx,vy,1,'Color','r');
	drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opticalFlow = vision.OpticalFlow('ReferenceFrameSource', 'Input port');                                                          
% flow = step(opticalFlow, rgb2gray(I1), rgb2gray(I2));
% opticFlow = opticalFlowHS;
% opticFlow = opticalFlowLK('NoiseThreshold',0.009);
% opticFlow = opticalFlowFarneback;
% opticFlow = opticalFlowLKDoG;
