function [] = doubleslitSimulationWithHusimi()
%Distance between slits in meter
dSlit = 133e-6;
%Wavelength in m 
lambda = 830e-9;
% length of free space path in meter
D = 0.6;

% %Intensity ratio between the slits:
% ratio = 0:0.1:1;
% ratioMatrix = ones(1,1, length(ratio));
% ratioMatrix(1,1,:) =  ratio;
% %phase difference between the two light sources:
% DeltaTheta = 0:0.1:2*pi;
% DeltaThetaMatrix = ones(1,1,1, length(DeltaTheta));
% DeltaThetaMatrix(1,1,1,:) =  DeltaTheta;
%Intensity and theta of first slit
intens1=3;
theta1 = 0.3;
%%
%Simulation of a husimifunction with phase theta between 0 and 2pi and
%intensity with standard deviation of 20
theta = 0 + 2*pi*rand(10000,1);  
 intens = normrnd(0,20,[10000,1]);
 x=zeros(10000,1);
 y=zeros(10000,1);

for i=1:length(x)
    x(i) = intens(i)*cos(theta(i));
    y(i) = intens(i)*sin(theta(i));
end
 
%amplitude and phase selection
minTheta = 0;
maxTheta = 0.5;
intensR = 10;
intensW = 1;
minInt = intensR-intensW/2;
maxInt = intensR+intensW/2;

hustheta =atan2(x,y);
Sel = find(hustheta < maxTheta & hustheta > minTheta &  sqrt(x.^2+y.^2)>minInt & sqrt(x.^2+y.^2)<maxInt);

% xsel=x(thetaSel);
% ysel=y(thetaSel);
intensSel=intens(Sel);
intensSel = intensSel*(intens1/intensR);
thetaSel=theta(Sel);

intens1 = ones(length(Sel),1)*intens1;
theta1 = ones(length(Sel),1)*theta1;
%% 

%horizontal Position after the free space path (zero is middle) in meter:
position = -0.01:0.0001:0.01;

%Path lengths from doubleslit to measurement points for slit 1 (left from
%the middle):
s1 = sqrt((position + dSlit/2).^2 + D^2);
s2 = sqrt((position - dSlit/2).^2 + D^2); 

%Phases from the path:
phase1 = mod(2*pi/lambda * s1, 2*pi);
phase2 = mod(2*pi/lambda * s2, 2*pi);
phase1Matrix=repmat(phase1,length(Sel),1);
phase2Matrix=repmat(phase2,length(Sel),1);


%iterate through all phases
allPhases = 0:0.01:2*pi;
allPhasesMatrix = repmat(allPhases,length(Sel),1);
allPhasesMatrix=allPhasesMatrix.';
two=allPhases+phase1;

% Intensity pattern:
pattern = zeros(length(position),2);
for i=1:length(position)
    Int = zeros(length(Sel),1);
 for j=1:length(Sel)
        Int(j) = abs(sqrt(intens1(j))*cos(theta1(j) + phase1Matrix(j,i))  +  sqrt(intensSel(j))*cos(thetaSel(j) + phase2Matrix(j,i)) ).^2;
 end
 meanInt = mean(Int);
 pattern(i,1) = position(i);
 pattern(i,2) = meanInt;
end


%averaging over all phases of one light source 
%(corresponding to averaging over one cycle of the electrical field (with fixed phase difference
%between the light sources):

Int = mean(Int,1);

%plot on one position, all deltaThetas for intensity ratio 1:
IntReshaped = Int(:,:,end,:);
IntReshaped = reshape(IntReshaped,[size(IntReshaped,2) size(IntReshaped,4)]);
plot(DeltaTheta,IntReshaped(1,:), '-');
xlabel('Delta Theta');
ylabel('Intensity on one position');
ylim([0 2.2]);
graphicsSettings();
savefig(['PatternForDeltaThetasOnePosition-' num2str(min(DeltaTheta)) '-' num2str(max(DeltaTheta)) '-IntensityRatio-' num2str(ratio(end)) '.fig']);
clf();

%plot all intensity ratios for one DeltaTheta:
IntReshaped = Int(:,:,:,1);
IntReshaped = reshape(IntReshaped,[size(IntReshaped,2) size(IntReshaped,3)]);
plot(IntReshaped);
plot(position*1000,IntReshaped, '-');
xlabel('Position (mm)');
ylabel('Intensity');
ylim([0 2.2]);
graphicsSettings();
savefig(['PatternForIntensityRatios-' num2str(min(ratio)) '-' num2str(max(ratio)) '-DeltaTheta-' num2str(DeltaTheta(1)) '.fig']);
clf();

%averaging over the intensity ratios:
IntAveraged = mean(Int,3);
plot(position*1000,IntAveraged(1,:,1,1), '-');
xlabel('Position (mm)');
ylabel('Intensity');
ylim([0 2.2]);
graphicsSettings();
savefig(['PatternForAveragedIntensityRatios-' num2str(min(ratio)) '-' num2str(max(ratio)) '-DeltaTheta-' num2str(DeltaTheta(1)) '.fig']);
clf();

%plot all deltaThetas for intensity ratio 1:
IntReshaped = Int(:,:,end,:);
IntReshaped = reshape(IntReshaped,[size(IntReshaped,2) size(IntReshaped,4)]);
plot(IntReshaped);
plot(position*1000,IntReshaped, '-');
xlabel('Position (mm)');
ylabel('Intensity');
ylim([0 2.2]);
graphicsSettings();
savefig(['PatternForDeltaThetas-' num2str(min(DeltaTheta)) '-' num2str(max(DeltaTheta)) '-IntensityRatio-' num2str(ratio(end)) '.fig']);
clf();

%averaging over the deltaThetas for one intensity ratio:
IntAveraged = mean(Int,4);
plot(position*1000,IntAveraged(1,:,end,1), '-');
xlabel('Position (mm)');
ylabel('Intensity');
ylim([0 2.2]);
graphicsSettings();
savefig(['PatternForAveragedDeltaThetas-' num2str(min(DeltaTheta)) '-' num2str(max(DeltaTheta)) '-IntensityRatio-' num2str(ratio(end)) '.fig']);
clf();

%averaging over deltaThetas and intensity ratios:
IntAveraged = mean(Int,4);
IntAveraged = mean(IntAveraged,3);
plot(position*1000,IntAveraged(1,:,1,1), '-');
xlabel('Position (mm)');
ylabel('Intensity');
ylim([0 2.2]);
graphicsSettings();
savefig(['PatternForAveragedDeltaThetas-' num2str(min(DeltaTheta)) '-' num2str(max(DeltaTheta)) '-AveragedIntensityRatios-' num2str(min(ratio)) '-' num2str(max(ratio)) '.fig']);
clf();


end