
folder='mat-data';

[filenames,~,positions]= getParametersFromFilenames('Folder',folder,'Parameter','position');
%fitCoeff = struct; 
quantities = struct;

for j = 1:length(filenames)

filename = cell2mat(filenames(j));
load([folder '/' filename], 'X1','X2','X3','O1','O2','O3','oTheta','oThetaMira');
  quantities.position(j) = positions(j);
chAssign = [2,1,3];

if exist('X1','var')
        % make sure all have the same number of pulses. 
         if ~isequal(size(X1,1),size(X2,1),size(X3,1))
             X1 = X1(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
             X2 = X2(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
             X3 = X3(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
         end
 
        % set which channel ist the target channel etc
        quadratures = zeros([size(X1) 3]);
        quadratures(:,:,:,1) = X1;
        quadratures(:,:,:,2) = X2;
        quadratures(:,:,:,3) = X3;
        Xtg = quadratures(:,:,:,chAssign(1));  
        XpsFast = quadratures(:,:,:,chAssign(2));
        XpsSlow = quadratures(:,:,:,chAssign(3));
        clear('quadratures'); 
end
[nTg,nPsFast,nPsSlow] = nPhotons(Xtg,XpsFast,XpsSlow);            
%phiSimulation=phiDifference(j);
 
%%
% phaseselection = -pi:0.5:pi ;
p = pi ;

% radius and thickness for amplitude selection
r=10;
d=2;
w=0.5;

 [nPsFast,nPsSlow,ntg] = nPhotons(O1,O2,O3);
 At=sqrt(2)*10*sqrt((nPsFast+nPsSlow)*ntg)/(1+nPsFast+nPsSlow);

%nTg = nPhotons(X2,X2,X2);
%[nDs,nDsStd] = deal(length(phaseselection));

  
% 
       selParams = struct('Type','phase','Position',[p,w,r,d]);
%        %selParams = struct('Type','phaseAndAmplitude','Position',[0.5,0.01,r,d]);
 % [~,~,iSelect] = selectRegion(O1,O2,O3,oTheta,phiSimulation,selParams);%,'Plot','show','Filename',['testTimes2' num2str(phaseselection(i))]);
[~,~,iSelect] = selectRegion(O1,O2,O3,oTheta,selParams);
  % 
% % 
%       
        selO3=O3(iSelect);
        selO2=O2(iSelect);
        selO1=O1(iSelect);
        [n1sel,n2sel,n3sel] = nPhotons(selO1,selO2,selO3);
%         nDs1(i) = n1sel;
%         nDs2(i) = n2sel;
%         nDs3(i) = n3sel;
%        nOtg = mean(nDs3);
        Aps= At*(1+(n1sel+n2sel))/sqrt(2*((n1sel+n2sel)*ntg));


 

        r=Aps;
       selParams = struct('Type','phaseAndAmplitude','Position',[p,w,r,d]);
      % selParams = struct('Type','phase','Position',[phaseselection(i),0.5]);
%         [selX,selTheta,iSelect] = selectRegion(O1,O2,O3,oTheta,selParams); %,'Plot','show','Filename',['test' num2str(phaseselection(i)) '-assessTheta']);
%       
        %[selX,selTheta,iSelect] = selectRegionAroundZero(O1,O2,O3,oTheta,oThetaMira,phiSimulation,selParams);
        [selX,selTheta,iSelect] = selectRegionAroundZero(O1,O2,O3,oTheta,oThetaMira,selParams);
        thetaMiraSel = oThetaMira(iSelect);
        
     % compute photon number of postselected quadratures in the doubleslit
    
    [nValues] = deal(zeros(10,1));
    for iN=1:10
        try
            uniformX = seriesUniformSampling(selX,thetaMiraSel,'NBins',100); %thetaMira oder selTheta?
        catch
            disp([num2str(positions(j)) 'Problem with uniformSampling.', ...
                'Use X without uniform sampling.']);
            uniformX = selX;
        end  
        [nValues(iN),~,~] = nPhotons(uniformX,uniformX,uniformX);

    end 
     
    
    nDs = mean(nValues); % postselected photon numbers in the doubleslit
    quantities.nDs(j) = nDs;
    quantities.nDsStd(j) = std(nValues);

    if (exist('nTg','var'))
        quantities.nTg(j) = nTg;
        quantities.nPsFast(j) = nPsFast;
        quantities.nPsSlow(j) = nPsSlow;
    end




end
figure(20);
plot(quantities.position,quantities.nDs,'o-');
xlabel('Position (mm)');
ylabel('Postselected photon number');
graphicsSettings();
ylim([0.5 2.5]);
savefig('nDsVsPosition-0-2pi_Phasepi.fig');
print('nDsVsPosition-0-2pi_PhasePi.png','-dpng','-r300');

figure(22);
plot(quantities.position,quantities.nDs./quantities.nTg,'o-');
xlabel('Position (mm)');
ylabel('n_{postselected}/n_{Tg,not postselected} ');
graphicsSettings();
ylim([0 2]);
savefig('nDsNormalizedVsPosition-0-2pi_Phasepi.fig');
print('nDsNormalizedVsPosition-0-2pi_PhasePi.png','-dpng','-r300');



            
