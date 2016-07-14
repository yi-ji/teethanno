%
% This function reads an smf file and returns its vertices and faces
%
% [F X status] = readSmf(filename)
%
% Input:
%   - filename: name of smf file to load
%
% Output:
%   - F: <m x 3> face list. F(i, :) denotes the three vertices of face
%   'i', indexed in X.
%   - X: <n x 3> vertex list. X(i, :) denotes the three coordinates of
%   vertex 'i'
%   - status: is 0 if no error ocurred when loading the file, otherwise
%   it is 1
%
function [F X status] = readSmf(filename)

    fid = fopen(filename);
    status = 0;
    if fid == -1
        disp('ERROR: could not open file');
        F = []; X = []; status = 1;
        return;
    end
    vnum = 1;
    fnum = 1;
    while(feof(fid) ~= 1)
        line = '';
        line = fgetl(fid);
        if length(line) == 0
            continue;
        end
        if line(1) == 'v'
            line = line(3:length(line));
            X(vnum, :) = sscanf(line, '%f %f %f');
            vnum = vnum + 1;
        elseif line(1) == 'f'
            line = line(3:length(line));
            F(fnum, :) = sscanf(line, '%d %d %d');
            fnum = fnum + 1;
        end
    end
    %X1 = [inf inf inf];
    %for i = 1:vnum
    %    for j = 1:length(X1(:,1))
    %        if X(i,:) == X1(j,:)
    %            break;
    %        end
    %    end
    %    if j == length(X1(:,1));
    %        X1 = [X1;X(i,:)];
    %    end
    %end
    %X = X1(2:length(X1(:,1)),:);
    fclose(fid);
