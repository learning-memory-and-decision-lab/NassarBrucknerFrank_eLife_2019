function z=nanzscore(x)
    %zscore x while ignoring nans;
    z=nan(size(x));
    for j = 1:size(x, 2)
        sel=isfinite(x(:,j));
        z(sel,j)=zscore(x(sel,j))
    end
    

