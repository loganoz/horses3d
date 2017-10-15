%
%   Start
clear all
clc
close all
format long

[fileName,pathName] = uigetfile();
file = [pathName,fileName];


%
%   Open file
%   ---------
    fid = fopen(file,'rb');
%
%   Create empty structure
%   ----------------------
    file = struct();
%
%   Read file name
%   --------------
    hr=fread(fid, 1, 'int32');            
    file.name=setstr(fread(fid, hr, 'uchar'))';
    hr=fread(fid, 1, 'int32');             
%
%   Get file type
%   -------------
    hr=fread(fid, 1, 'int32');
    file.type=fread(fid, hr/4, 'int32');
    hr=fread(fid, 1, 'int32');
%
%   Read the number of elements
%   ---------------------------
    hr = fread(fid, 1, 'int32');
    file.no_of_elements = fread(fid, hr/4, 'int32');
    hr = fread(fid, 1, 'int32');
%
%   Read iteration and time
%   -----------------------
    hr=fread(fid,1,'int32');
    file.iter = fread(fid, hr/4, 'int32');
    hr=fread(fid,1,'int32');
        
    hr=fread(fid,1,'int32');
    file.time = fread(fid, hr/8, 'float64');
    hr=fread(fid,1,'int32');
%
%   Read reference values
%   ---------------------
    hr=fread(fid,1,'int32');
    file.refs = fread(fid, hr/8, 'float64');
    hr=fread(fid,1,'int32');
%
%   Read the beginning data code
%   ----------------------------
    hr = fread(fid,1,'int32');
    beginning = fread(fid, hr/4, 'int32');
    if ( beginning ~= 88 )
        stop;
    end
    hr = fread(fid,1,'int32');
%
%   Read data
%   ---------
    file.data = cell(1,file.no_of_elements);
    
    for e = 1 : file.no_of_elements
%
%       Get number of dimensions
%       ------------------------
        hr = fread(fid,1,'int32')/4;
        nDim = fread(fid,hr,'int32');
        hr = fread(fid,1,'int32');
%
%       Get dimension vector
%       --------------------
        hr = fread(fid,1,'int32')/4;
        N  = fliplr(fread(fid,hr,'int32')');
        hr = fread(fid,1,'int32');
%
%       Get data
%       --------
        hr = fread(fid,1,'int32')/8;
        data = fread(fid,hr,'float64');
        hr = fread(fid,1,'int32');
        
        if ( nDim == 4 )
            file.data{e} = zeros(N);
            counter = 0;
            for k = 1:N(1)
                for j = 1:N(2)
                    for i = 1:N(3)
                        for d = 1:N(4)
                            counter = counter + 1;
                            file.data{e}(k,j,i,d) = data(counter);
                        end
                    end
                end
            end
            
        else
            fprintf('Number of dimensions %d not implemented\n',nDim);
        end
    end
    
%
%   Close file
%   ----------
    fclose(fid);