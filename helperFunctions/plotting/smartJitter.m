function [xOffset] = smartJitter(yData,xMargin,yMargin)
% Helps plot individual points in categorical data.
%
% For a column of y values, determines appropriate offsets on the x axis
% so all points are visible. 
%
% After running this function, xOffset will need to be added to the 
% mean x value for the category. 
%
% Inputs:
%   yData is a column vector, matrix, or cell array, holding the y values
%       to be plotted. If a matrix, each column is treated as a group. If a
%       cell array, each entry is treated as a group (and must hold a
%       column vector). 
%   xMargin (scalar) determines how far points will be jittered on x if
%       they overlap on y. This determines the scaling of the results in
%       xOffset.
%   yMargin (scalar) determines the threshold for treating 2 points as
%       overlapping on y, in units of the y axis. 
%
% Output:
%   xOffset is formatted the same way as yData (column, matrix, or cell
%       array). The offset values in each group are zero-centered, so means
%       for each category subsequently need to be added back in. 
%
% by JTM, 2/13/2013 (adapted from previous code)



% process input
% check the two scalar args
assert(isscalar(xMargin),'xMargin must be a scalar.');
assert(isscalar(yMargin),'yMargin must be a scalar.');
% yData can originally be a column vector, a matrix, or a cell array
if ~iscell(yData) % if a column or matrix
    nGrps = size(yData,2);
    xOffset = nan(size(yData)); % initialize output
else
    nGrps = numel(yData);
    xOffset = cell(size(yData)); % initialize output
end

% loop over groups
for g = 1:nGrps

    % get data from the input variable
    if iscell(yData)
        y = yData{g};
    else
        y = yData(:,g);
    end
    assert(iscolumn(y),'y values in each group must be a column vector');
    nObs = length(y);

    % initialize results
    x = zeros(size(y)); 

    % adjust x coords so points do not overlap
    for j = 2:nObs
        % index those _previous_ points that vertically overlap the current one
        conflicting = abs(y(j) - y(1:(j-1))) < yMargin;
        if any(conflicting)
            incrSign = 1; % will add a pos or neg value
            incrN = 1;
            while ismember(x(j),x(conflicting)) % so long as this point still overlaps something
                x(j) = incrSign*incrN*xMargin;
                incrSign = incrSign*-1; % switch between leftward and rightward offset
                if incrSign==1, incrN = incrN+1; end % increase horizontal offset
            end
        end
    end

    % place results into the output variable, which is either a cell array
    % or a matrix depending on the format of the input
    if iscell(xOffset)
        xOffset{g} = x;
    else
        xOffset(:,g) = x;
    end
    
end % loop over groups





