function [] = write_fileLoc( AQpath,time )
    % writes directories where input and output files are located on the disk to a text file

file=['FileLocations_' num2str(time) '.txt'];
filename=[AQpath file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% Enter Location of Input files %%');
fprintf( fid, '%s\n',[AQpath 'Input']);
fprintf( fid, '%s\n','%% Enter Location of Output files %%');
fprintf( fid, '%s\n',[AQpath 'Output']);
fclose(fid);


end

