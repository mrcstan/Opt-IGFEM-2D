%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/26/2016
%%% Copyright 2016 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function append a row of values to a file
%%% If the file does not exist, the header will be written first
%%% If the file already exists, the header will be omitted
function output_values (outfile, header, counter, rowValues)
if counter == 1
    if ~ischar(header)
        error('header should be a string')
    end
    fileID = fopen(outfile, 'w');
    fprintf(fileID,'%s\n',header);
else
    fileID = fopen(outfile, 'a');
end

fprintf(fileID,'%i\t',counter);

if ~isempty(rowValues)
    fprintf(fileID,'%g\t',rowValues);
end
fprintf(fileID,'\n');
fclose(fileID);
end