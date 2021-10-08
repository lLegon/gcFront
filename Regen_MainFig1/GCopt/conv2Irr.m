function [modelIrrev,irrevFluxDistWT] = conv2Irr(model,fluxDist)
% Converts reversible model to an irreverseble structure and transforms the
% flux distribution "fluxDist" accordingly

% Convert reversible model to irreversible format
[modelIrrev,~,rev2irrev,irrev2rev] = convertToIrreversible(model);

% delivers a conversion matrix for irr 2 rev flux distributions
mapIrr2Rev  = zeros(length(rev2irrev),length(irrev2rev));

if nargin < 2
    irrevFluxDistWT     = [];
    flagFluxDist        = 0;
else
    % Convert wild type flux distribution accordingly
    irrevFluxDistWT     = zeros(length(irrev2rev),1);
    flagFluxDist        = 1;
end

for i=1:length(rev2irrev)
    irrRxnNum   = rev2irrev{i};
    if length(irrRxnNum)>1
        mapIrr2Rev(i,irrRxnNum(1))      = 1;
        mapIrr2Rev(i,irrRxnNum(2))      = -1;
        if flagFluxDist
            if fluxDist(i)<0
                irrevFluxDistWT(irrRxnNum(1))   = 0;
                irrevFluxDistWT(irrRxnNum(2))   = -fluxDist(i);
            else
                irrevFluxDistWT(irrRxnNum(1))   = fluxDist(i);
                irrevFluxDistWT(irrRxnNum(2))   = 0;
            end
        end
    else
        % account for irreversible back fluxes
        stoichCoeff         = find(model.S(:,i));
        if ~isempty(stoichCoeff)
            stoichCoeff         = stoichCoeff(1);
            if model.S(stoichCoeff,i)*modelIrrev.S(stoichCoeff,irrRxnNum(1)) > 0
                mapIrr2Rev(i,irrRxnNum(1))      = 1;
                if flagFluxDist
                    irrevFluxDistWT(irrRxnNum)      = fluxDist(i);
                end
            else
                mapIrr2Rev(i,irrRxnNum(1))      = -1;
                if flagFluxDist
                    irrevFluxDistWT(irrRxnNum)      = -fluxDist(i);
                end
            end
        else
            mapIrr2Rev(i,irrRxnNum(1))      = 1;
            if flagFluxDist
                irrevFluxDistWT(irrRxnNum)      = fluxDist(i);
            end
        end
    end
end
modelIrrev.mapIrr2Rev   = mapIrr2Rev;

end