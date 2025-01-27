% phaseselection = -pi:0.5:pi ;
phaseselection = 0:0.5:2*pi ;

nTg = nPhotons(X2,X2,X2);
[nDs,nDsStd] = deal(length(phaseselection));
for i  = 1:length(phaseselection)  

       selParams = struct('Type','phaseAndAmplitude','Position',[phaseselection(i),0.5,10,1]);
      % selParams = struct('Type','phase','Position',[phaseselection(i),0.5]);
%         [selX,selTheta,iSelect] = selectRegion(O1,O2,O3,oTheta,selParams); %,'Plot','show','Filename',['test' num2str(phaseselection(i)) '-assessTheta']);
%       
        [selX,selTheta,iSelect] = selectRegionOfTotalPhase(O1,O2,O3,oTheta,oThetaMira,selParams);
         thetaMiraSel = oThetaMira(iSelect);
     % compute photon number of postselected quadratures in the doubleslit
    
    [nValues] = deal(zeros(10,1));
    for iN=1:10
        try
            uniformX = seriesUniformSampling(selX,thetaMiraSel,'NBins',100); %thetaMira oder selTheta?
        catch
            warning(['Problem with uniformSampling.', ...
                'Use X without uniform sampling.']);
            uniformX = selX;
        end  
        [nValues(iN),~,~] = nPhotons(uniformX,uniformX,uniformX);

    end 
    
    nDs(i) = mean(nValues);
    nDsStd(i) = std(nValues);
end

figure(16);
hold on;
errorbar(phaseselection,nDs./nTg,nDsStd,'o-');
ylim([0.8 2]);
xlabel('postselected phase');
ylabel('postselected photon number');