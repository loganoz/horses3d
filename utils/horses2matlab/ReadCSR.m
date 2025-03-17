function A = ReadCSR(FileName)

%% //////////////////////////////////////////////////////
%
%    @File:    ReadCSR.m
%    @Author:  Wojciech Laskowski (wj.laskowski@upm.es)
%    @Created: Mon Jul 24 08:11:10 2021
%
%    Function translates matrix in the CSR format
%    obtained from HORSES3D (CSRMatrix % Visualize) to
%    MATLAB ijk format.
%
%    'Filename' is an ASCII file with the CSR matrix.
%
%% //////////////////////////////////////////////////////

% ************ OPEN *************
fileID = fopen(FileName,'r');
% *******************************

% dimensions
dimy = fgetl(fileID);
dimy = str2num(dimy);
n = dimy(1);
nnz = dimy(2);

% allocate
Rows = zeros(n+1,1);
Cols = zeros(nnz,1);
Vals = zeros(nnz,1);

% rows
for i=1:n+1
    Rows(i) = str2num(fgetl(fileID));
end
disp('Rows read!');

% cols
for i=1:nnz
    Cols(i) = str2num(fgetl(fileID));
end
disp('Columns read!');


% vals
for i=1:nnz
    Vals(i) = str2double(fgetl(fileID));
end
disp('Values read!');

% *********** CLOSE *************
fclose(fileID);
% *******************************


%% Rows2
Rows2 = Cols*0;

for i = 1 : length(Rows) - 1
    tmp = Rows(i):Rows(i+1)-1;
    Rows2(tmp) = i;
end

%% transfer to sparse matrix
A = sparse(Rows2,Cols,Vals,max(Cols),max(Cols));

return
