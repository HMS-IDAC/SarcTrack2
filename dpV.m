function ijPath = dpV(C,r0,c0)

[nRows, nCols, nPlns] = size(C);

pathCosts = zeros(size(C));

maxDisp = 3;

pathCosts(:,:,1) = Inf;
pathCosts(r0,c0,1) = 0;
pathCosts(:,[1:maxDisp end-maxDisp+1:end],:) = Inf;
pathCosts([1:maxDisp end-maxDisp+1:end],:,:) = Inf;

predecessorsI = nan(size(C));
predecessorsJ = nan(size(C));

for k = 2:nPlns
    for i = maxDisp+1:nRows-maxDisp
        for j = maxDisp+1:nCols-maxDisp
            localCosts = pathCosts(i-maxDisp:i+maxDisp,j-maxDisp:j+maxDisp,k-1);
            if ~any(~isinf(localCosts(:))) % all local costs = Inf
                rim = maxDisp+1; cim = maxDisp+1; m = Inf;
            else
                m = min(localCosts(:));
                [rim,cim] = find(localCosts == m);
                rim = rim(1); cim = cim(1);
            end
            predRow = i+rim-(maxDisp+1);
            predCol = j+cim-(maxDisp+1);
            Cijk = C(i,j,k);
            updatedCost = m+Cijk;%+5*abs(C(predRow,predCol,k-1)-Cijk);%+0.1*sqrt((rim-maxDisp-1)^2+(cim-maxDisp-1)^2);
            pathCosts(i,j,k) = updatedCost;
            predecessorsI(i,j,k) = predRow;
            predecessorsJ(i,j,k) = predCol;
        end
    end
end

m = min(min(pathCosts(:,:,end)));
[rim,cim] = find(pathCosts(:,:,end) == m);
rim = rim(1); cim = cim(1);

ijPath = zeros(2,nPlns);
ijPath(:,nPlns) = [rim,cim];
for k = nPlns-1:-1:1
    ijPath(:,k) = [predecessorsI(ijPath(1,k+1),ijPath(2,k+1),k+1);...
                   predecessorsJ(ijPath(1,k+1),ijPath(2,k+1),k+1)];
end

end