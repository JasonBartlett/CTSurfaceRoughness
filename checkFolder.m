function [] = checkFolder(folder)
%checkFolder Checks to see if a folder exists and if it does not creates it
if exist(folder,'file')~=7
    mkdir(folder)
end
end