function varargout = bin2mat(file,precf,mff)
% data = BIN2MAT(file,precf,mff)
%
% Reads data stored in a binary file written by, for example, a FORTRAN
% program. The function assumes that data chunk(s), i.e. write steps, are
% accompanied by header(s) and tail(s) that give the size in bytes of that
% chunk, and it uses that information to only return the actual stored data
% of interest as a single vector. Headers/tails are assumed to be store 
% with 32-bit integer precision, which seems to be standard for FORTRAN.
%
% NOTE: For now, this function can only handle data stored as 32-bit signed
% integers or 64-bit double-precision floats.
%
% INPUT:
% 
% file    String giving the binary file to read in
% precf   Precision flag giving precision class of binary file data
%            OPTIONS:  'int'    = signed integer (4 bytes)
%                      'double' = double-precision floating point (8 bytes)
% mff     Machine format flag giving the encoding of the binary file data
%            OPTIONS:  'b' = Big-ending byte ordering 
%                      'l' = Little-ending byte ordering
%                      See FOPEN for more options.
%
% OUTPUT:
%
% data    Data stored in the binary file, excluding headers/tails
%
% SEE ALSO: FOPEN, FREAD, FSEEK, FTELL
%
% Last modified by gleggers-at-gatech.edu, 11/04/2018
% This function was derived from a similar one by charper4-at-uoregon.edu

% Precision class of the binary file headers and footers
hfp = 'int32';

% Based on the passed precision flag, get the conversion factor AKA the
% size in bytes of one entry of the passed precision class
switch precf
    case 'int'
        pcf = 4;
    case 'double'
        pcf = 8;
end

% Open the binary data file
fid = fopen(file,'rb',mff);

% Get the initial current file position indicator
cof = ftell(fid);

% Move to the end of the file, get the file position indicator, and return
fseek(fid,0,'eof');
eof = ftell(fid);
fseek(fid,0,'bof');

% Initialize the data array to store read-in binary data;
data = [];

% So long as file position indicator has not reached the end of file...
while cof < eof
    
    % Read the header giving the size for a data chunk
    header = fread(fid,1,hfp);
    
    % Number of lines in the data chunk
    num = header/pcf;
    
    % Read in the data chunk and append it to prior chunks
    chunk = fread(fid,num,precf);
    data = [data; chunk];
    
    % Read the tail    
    fread(fid,1,hfp);
    
    % Update the current file position indicator
    cof = ftell(fid);
end

% Close the binary data file
fclose(fid);

% Generate output as needed
vars = {data};
varargout = vars(1:nargout);