function [RunName, c] = OutputName(c, tag)
if c.ConductionModel == 1
    RunName = sprintf('%s_%i_ConductionOnly',...
        tag, c.year);
else
    if c.do_no_darcy == 1
        s3 = '_no_darcy';
    else
        s3 = sprintf('_darcy_wh%0.2f',c.whwice);
    end

    if c.calc_CLliq == 1
        text_Si = 'CL';
    else
        text_Si = spritnf('%0.2f',c.liqmax);
    end

    if length(c.year)==1
        RunName = sprintf('%s_%i_Si%s_pr%0.3f_Ck%0.2f%s',...
            tag, c.year,  text_Si ,c.prec_rate, c.Ck, s3);
    else
        RunName = sprintf('%s_%i-%i_Si%s_pr%0.3f_Ck%0.2f%s',...
            tag, c.year(1),c.year(2), text_Si ,c.prec_rate, c.Ck, s3);
    end

    c.OutputFolder = sprintf('./Output/%s',RunName);
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