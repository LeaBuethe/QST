%clear;
%Distance between slits in meter
dSlit = 133e-6;
%Wavelength in m 
lambda = 830e-9;
% length of free space path in meter
D = 0.062;
%StartPhase, wenn husimiSelection vorausgesetzt wird, sonst wäre
%Phasendifferenz nicht fest.
theta1 = 0;
theta2 = 0;
%StartAmplitude
A1=1.1;
A2=0.6;


%factor for calculating offsets
const=0;
%Offsets in doubleslitfiber
C1=const*rand(1,1);
C2=const*rand(1,1);
%other Offsets
C_LO=const*rand(1,1);
C_LOtg=const*rand(1,1);
Ctg=const*rand(1,1);
%Noise from doubleslitfiber
% D1=0*rand(40,1)+0.8;
% D2=0*rand(40,1)+0.8;
% D_LO=0.2*rand(1,1);
% D_LOtg=0.2*rand(1,1);
% Dtg=0.2*rand(1,1);

 
% for i = -0.002:0.00004:0.002
%     theta2(i)=pi*(rand(1,100);
% end
% theta = pi*rand(1,100);
% theta2 = rand(mittelung über Phasen für jede Position)!!!!!!!!!!!!!!!!
 
 
%horizontal Position after the free space path (zero is middle) in meter:
position = -0.003:0.000005:0.004;
%Noise for each position (for addition)
%D3=0.5*pi*rand(1,length(position));
%D4=0.7*pi*rand(1,length(position));
D3=0;
D4=0;
%amplitude noise because of husimiSelection of a endless amplitude range
%for selection:
%A1_range=0.8*rand(1,length(position))+0.9;
A2_range=0*rand(1,length(position));
%A1=A1*A1_range;
%A2=A2*A2_range;

%theta1 = ones(1,length(position))*theta1;
%Path lengths from doubleslit to measurement points for slit 1 (left from
%the middle):
s1 = sqrt((position + dSlit/2).^2 + D^2);
s2 = sqrt((position - dSlit/2).^2 + D^2); 
 
 
%Phases from the path:
% phase1 = mod(2*pi/lambda * s1,2*pi);
% phase2 = mod(2*pi/lambda * s2,2*pi);
 
 
%Phases from the path:
phase1 = 2*pi/lambda * s1;
phase2 = 2*pi/lambda * s2;
 
%Phasen einzeln am Sammelpunkt. Hier könnten noch Konstanten und Rauschen addiert werden.
thetaM = theta1 +  phase1+D3+C1;
thetaD = theta2 + phase2+D4+C2;
%phiLOtg=theta1+C_LO+D_LO+C_LOtg+D_LOtg;
phiLOtg=theta1+C_LO+C_LOtg+Ctg;

%für Korrektur mit ThetaMira wird die Phase1 wieder abgezogen. Anderes
%PhaseshiftMuster entsteht!
%phiLOtg=theta1+C_LO+C_LOtg+Ctg+phase1;


%phiSammel=atan((A1*sin(thetaD)+A2*sin(thetaM))./(A1*cos(thetaD)+A2*cos(thetaM)));
%for adding amplitude noise to the calculation:
% phiSammel=atan2((A1.*sin(thetaD)+A2.*sin(thetaM)),(A1.*cos(thetaD)+A2.*cos(thetaM)))+pi;
% amplSammel= sqrt(A1.^2+A2.^2+2*A1.*A2.*cos(thetaD-thetaM));
% % 
phiSammel=atan2((A1*sin(thetaD)+A2*sin(thetaM)),(A1*cos(thetaD)+A2*cos(thetaM)))+pi;
amplSammel= sqrt(A1^2+A2^2+2*A1*A2*cos(thetaD-thetaM));
phiDetection=mod(phiSammel-phiLOtg,2*pi); 

%Intensity proportional to amplitude^2 at collecting point
I=amplSammel.^2; 
 
 
thetaDiff=mod(thetaM-thetaD,2*pi);
thetaSum=mod(thetaM+thetaD, 2*pi);
onlyMira=mod(thetaM-phase1, 2*pi);
%figure(1);
%plot(position, onlyMira);
figure(1);
%plot(position,thetaM,'-o');
%plot(position,thetaD,'-o');
%plot(position, thetaDiff, '-o');
%plot(position, thetaSum, '-o');
plot(position, phiDetection);
ylabel('Phase');
xlabel('Position in m');
%title('Mit Rauschen D3=D4=0.4*pi const=0.2 Amplitudenrange 0.4');
% savefig('Phase_mit_Rauschen_D0,4_C0,2_Arange0,4');
% 
% figure(3);
% plot(position, I,'-o');
% ylabel('Intensität');
% xlabel('Position in m');
% title('Mit Rauschen D3=D4=0.4*pi const=0.2 Amplitudenrange 0.4');
% savefig('Intensitaet_mit_Rauschen_D0,4_C0,2_Arange0,4');
