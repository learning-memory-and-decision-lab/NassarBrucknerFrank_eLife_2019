function selBehav=selBehav(behav, sel, padNan)

%% Selects all elements specified in sel from each field in structure
%% behav.  If padNan=1, then array indices beyond the end of a given field
%% array will be returned as nan.  Otherwise a non-existent selection will
%% lead to crash...


if nargin<3|isempty(padNan)
    padNan=0;
end


if padNan==1
    behav=nanPadBehav(behav);
end

    
    




allNames=fieldnames(behav);



eval(sprintf('maxLn=length(behav.%s);', char(allNames(1))))


for i = 1:length(allNames)
    % get data from a field
        ln=eval(sprintf('length(behav.%s);', char(allNames(i))));
%         if ln~=maxLn
%             disp('variable length mismatch')
%             input('hit enter to coninue')
%         end
%     
        
        
            
    %     eval(sprintf('selBehav.%s(end+1:maxLn)=nan',  char(allNames(i))))
    %     end
    
    % not sure why i can't do this with cell arrays...
    try
        if ~eval(sprintf('iscell(behav.%s)', char(allNames(i))));
            eval(sprintf('selBehav.%s=behav.%s(sel);', char(allNames(i)), char(allNames(i))));
        else
            eval(sprintf('selBehav.%s=behav.%s(sel);', char(allNames(i)), char(allNames(i))));
        end
        
    catch
        keyboard
    end
    
    
end

