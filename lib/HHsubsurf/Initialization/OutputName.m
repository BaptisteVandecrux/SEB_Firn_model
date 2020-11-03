function [RunName, c] = OutputName(c)
if c.ConductionModel == 1
    RunName = sprintf('%s_%i_ConductionOnly',...
        tag, c.year);
else
    RunName = c.station;
    RunName = [RunName, sprintf('_%i',c.year)];
    if c.calc_CLliq == 1
        text_Si = 'CL';
    else
        text_Si = sprintf('%0.2f',c.liqmax);
    end
    RunName = [RunName, '_IWC_', text_Si];
    RunName = [RunName, sprintf('_%i_layers',c.jpgrnd-1)];
    
    c.OutputFolder = sprintf('%s/%s',c.OutputRoot,RunName);
    [~,~,id] =  mkdir(c.OutputFolder);
    count = 1;
    while ~isempty(strfind(id,'DirectoryExists'))
       count =count+1;

       c.OutputFolder = sprintf('./Output/%s_%i',RunName,count);
       [~,~,id] =  mkdir(c.OutputFolder);
    end
    if count>1
           RunName = sprintf('%s_%i',RunName,count);
    end
end