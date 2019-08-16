function allBehav=catBehav(behav, allBehav, forceCat)

% if forceCat is equal to one, the concatanation will proceed even if
% not all fields are present in both structures.  Fields that are not
% in both structures will be discarded.

if nargin<3|isempty(forceCat)
    forceCat=0;
end


allNames=fieldnames(behav);
maxLn=eval(sprintf('length(behav.%s);', char(allNames(1))));


for i = 1:length(allNames)
    
    if forceCat==0||isfield(allBehav, (allNames(i)))
        
        % get data from a field
        eval(sprintf('a=allBehav.%s;', char(allNames(i))));
        eval(sprintf('b=behav.%s;', char(allNames(i))));
        
        ln=eval(sprintf('length(behav.%s);', char(allNames(i))));
        if ln~=maxLn
            %disp('variable length mismatch')
            %input('hit enter to coninue')
        end
        
        
        
        % flip it if its sideways
        if size(b,1) < size(b, 2)
            b=b';
        end
        
        if size(a,1) < size(a, 2)
            a=a';
        end
        
        try
        % concatanate
        eval(sprintf('allBehav.%s = cat(1, a, b);', char(allNames(i))));
        catch
            keyboard
        end
        
            
            
    end
    
end

      