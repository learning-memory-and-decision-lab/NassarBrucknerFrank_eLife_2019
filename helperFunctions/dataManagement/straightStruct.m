function behav=straightStruct(behav)

    allNames=fieldnames(behav);
    
    for i = 1:length(allNames)
      eval(sprintf('a=behav.%s;', char(allNames(i))));
      if size(a, 2)>size(a,1)
          a=a';
      end
      
      eval(sprintf('behav.%s = a;', char(allNames(i))));
       
    end
    
      