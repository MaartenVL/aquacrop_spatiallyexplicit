function [x,row,col,hres,vres] = Read_Idrisi(fname)

%this function reads an Idrisi raster and rdc file and returns indicated
%variables

rasfname=[fname '.rst'];
docfname=[fname '.rdc'];

fclose('all');

%read doc file into text array
mydocFile=fopen(docfname,'r');
%fseek(mydocFile,0,-1);
idrformat=textscan(mydocFile,'%14c %20s',1,'whitespace','\n');
% idrfiletitle=textscan(mydocFile,'%14c %128s',1,'whitespace','\n');
idrdatatype=textscan(mydocFile,'%14c %20s',1,'whitespace','\n');
idrfiletype=textscan(mydocFile,'%14c %20s',1,'whitespace','\n');
idrcol=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrrow=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrrefsyst=textscan(mydocFile,'%14c %25s',1,'whitespace','\n');
idrrefunits=textscan(mydocFile,'%14c %25s',1,'whitespace','\n');
idrunitdist=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrminx=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrmaxx=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrminy=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrmaxy=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrposerr=textscan(mydocFile,'%14c %25s',1,'whitespace','\n');
idrres=textscan(mydocFile,'%14c %25s',1,'whitespace','\n');
idrminval=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrmaxval=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrdispmin=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrdispmax=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrvalunits=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrvalerror=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrflagval=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrflagdef=textscan(mydocFile,'%14c %f',1,'whitespace','\n');
idrlegcats=textscan(mydocFile,'%14c %u',1,'whitespace','\n');

row=idrrow{2};
col=idrcol{2};

hres=round((idrmaxx{2}-idrminx{2})/col);
vres=round((idrmaxy{2}-idrminy{2})/row);

if hres ~= vres
    disp('fout:de horizontale resolutie is gelijk aan verticale resultie. STOP');
    fprintf('hres %u vres %u',hres,vres);
end

close all

myrasfile=fopen(rasfname,'r');
x=fread(myrasfile,[col,row],'float32');
x=x';

close all
