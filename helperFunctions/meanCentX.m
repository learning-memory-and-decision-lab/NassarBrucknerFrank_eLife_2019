function xMat=meanCentX(xMat)

if std(xMat(:,1))==0
xMat(:, 2:end)=xMat(:, 2:end)-repmat(nanmean(xMat(:, 2:end), 1), length(xMat(:, 2:end)), 1);

else
 xMat=xMat-repmat(nanmean(xMat, 1), length(xMat), 1);
  
end

if any(any(~isfinite(xMat)))
    for i = 1:size(xMat,2)
        xMat(~isfinite(xMat(:,i)),i)=0;
    end
end
