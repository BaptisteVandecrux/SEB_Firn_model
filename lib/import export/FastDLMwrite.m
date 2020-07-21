    
function FastDLMwrite(filename,data,del)
    ms  = permute( data, [ 2, 1 ] );
    fid = fopen(filename, 'a');
    format = [repmat(sprintf('%%f%s',del),1,size(data,2)) '\n'];
    fprintf( fid, format, ms );
    fclose( fid );
end