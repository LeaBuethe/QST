folder='mat-data';

[filenames,~,positions]= getParametersFromFilenames('Folder',folder,'Parameter','position');
%fitCoeff = struct; 

photonnumbers=struct;

for j = 1:length(filenames)

filename = cell2mat(filenames(j));
load([folder '/' filename],'X2');


[~,ntg] = nPhotons(X2,X2,X2);
%photonnumbers.n1(j)=nps1;
%photonnumbers.n2(j)=nps2;
photonnumbers.ntg(j)=ntg;

end