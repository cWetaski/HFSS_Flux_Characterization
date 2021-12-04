% Use this script to transform oblique plane to orthonormal view
% Using the method from John Loomis johnloomis.org/ece564/notes/tform/

clear all;
close all;
clc

% SET PARAMETERS
% for the transformation
square_size = 2000; % set side length (in pixels) of transformed square
border_size = 1000;  % set size of border (in pixels) around transformed square
% the transformed image will be a square of with a side length of
% (square_size +  2*border_size).
% The total pixels in the image will be (square_size + 2*border_size)^2

folder = pwd; % get current folder
folder = strcat(folder,'\images\original');
folder = uigetdir(folder); % the user will select the set of images they wish to use.
folder_split = regexp(folder,'\','split'); % Used to get folder name when saving

% Get the flux readings from the spreadsheet which is in the same folder as
% the images. The sensor readings are recorded in column B, with no header.
flux_readings = dir(fullfile(folder,'*.xlsx')); % get the list of spreadsheets from the folder
flux_readings = fullfile(folder,flux_readings(1).name); % get the spreadsheet file name (only 1 spreadsheet should be in the folder)
flux_readings = readtable(flux_readings); % get the actual data from the spreadsheet
flux_readings = flux_readings(:,2); % the readings are stored in the second column
flux_readings = flux_readings{:,:}; % convert from a table to an array of type double
flux_readings = cellstr(num2str(flux_readings,'%06.3f')); % convert from doubles to cell array of strings with 2 digits before decimal and 3 after

image_files = dir(fullfile(folder,'*.TIF')); % get the image files
image_files_overexp = dir(fullfile(folder,'overexposed','*.TIF')); % get the overexposed image files
for k = 1 : length(image_files) % get a list of images in the images folder
% 	fprintf('%d: %s\n', k, files(k).name);
	[~, base_file_name, extension] = fileparts(image_files(k).name);
	image_names{k} = base_file_name;
    [~, base_file_name_overexp, extension] = fileparts(image_files_overexp(k).name);
	image_names_overexp{k} = base_file_name_overexp;
end

% Pick one image from a set that will be transformed to be orthonomal.
% All images must be shot from the same location and have the flux in the
% same absolute location within the image.
button = menu('Use which gray scale demo image?', image_names_overexp); % Display all image file names in a popup menu.
% Get the base filename.
base_file_name = image_files_overexp(button).name; % Assign the one on the button that they clicked on.
% Get the full filename
bright_image_file = dir(fullfile(fullfile(folder,'overexposed','*.TIF')));
bright_image_name = fullfile(folder,'overexposed',bright_image_file(button).name);
bright_image = im2double(imread(bright_image_name));

% apply gamma adjustment to see contrast in shadows.
bright_img_gamma = imadjust(bright_image,[],[],0.15);
f1 = figure;
imshow(bright_img_gamma)
axis on
hold on

% click to identify corners of alumina sheet to get transform
% note that you must select corners in clockwise order starting from the
% top-left. Otherwise it will not warp correctly. 
go_to_tform = false;

while ~go_to_tform
    waitfor(msgbox('select corners of square plane (clockwise starting from top-left corner)'));
    [x y] = ginput(4);
    corners_plot = plot([x;x(1)],[y;y(1)],'Marker','x','MarkerSize',10,'Color','r','MarkerEdgeColor','Magenta','Linewidth',0.5); % burn points into image
    base = [0 0; 1 0; 1 1; 0 1]; % transform to square, since alumina sheet is a square
    answer = questdlg('Continue to transform or re-pick corners?','Continue?','Continue','Re-pick corners','Continue');
    switch answer
        case 'Continue'
            waitfor(msgbox('continue'));
            go_to_tform = true;
        case 'Re-pick corners'
            delete(corners_plot);
            go_to_tform = false;
    end
end

tf = fitgeotrans([x y],base*square_size,'projective'); % choose to make alumina square have 2000 'pixel' width in warped view

% create destination folder if it does not exist.
destinationFolder = pwd;
destinationFolder = char(strcat(destinationFolder,'/images/transformed/',folder_split(length(folder_split))));
if ~exist(destinationFolder, 'dir')
    mkdir(destinationFolder);
end

% create destination folder for overexposed images if it does not exist.
destinationFolder_overexp = char(strcat(destinationFolder,'/overexposed'));
if ~exist(destinationFolder_overexp, 'dir')
    mkdir(destinationFolder_overexp);
end

% Transform each image in the regular set and save
for k = 1 : length(image_files) % For each image in the set
    
    base_file_name = char(image_names(k));
    fullFileName = fullfile(folder,strcat(base_file_name,'.tif')); % get the full file name
    img = im2double(imread(fullFileName)); % get the image and convert to double
    if size(img,3) > 1 % if there is more than 1 colour channel (this is the case if image is from rawproc)
        img = img(:,:,1); % just take the first channel, since all three channels should be identical 
    end
    [img2, img2_ref] = imwarp(img,tf); % make image orthonormal with the transform we determined
    
    % change world limits so that there is 500 pixel border around corners
    % of alumina sheet. 
    
    img2_ref.XWorldLimits = [-border_size, square_size + border_size];
    img2_ref.YWorldLimits = [-border_size, square_size + border_size]; 
    img2_ref.ImageSize = [square_size + 2*border_size, square_size + 2*border_size];
    [img3, img3_ref] = imwarp(img,tf,'OutputView',img2_ref); % reapply transform using the world limits we altered.
    
    
    img3 = im2uint16(img3); % convert img back to 16bit img
    outputName = [base_file_name '_' flux_readings{k} '_tform.tif']; % make new file name
    fullDestinationName = fullfile(destinationFolder,outputName); % get full file name
    imwrite(img3,fullDestinationName); % save the file
end

% Transform each image in the overexposed set and save
for k = 1 : length(image_files_overexp) % For each image in the set
    
    base_file_name = char(image_names_overexp(k));
    fullFileName = fullfile(folder,'overexposed',strcat(base_file_name,'.tif')); % get the full file name
    img = im2double(imread(fullFileName)); % get the image and convert to double
    if size(img,3) > 1 % if there is more than 1 colour channel (this is the case if image is from rawproc)
        img = img(:,:,1); % just take the first channel, since all three channels should be identical 
    end
    [img2, img2_ref] = imwarp(img,tf); % make image orthonormal with the transform we determined
    
    % change world limits so that there is 500 pixel border around corners
    % of alumina sheet. 
    
    img2_ref.XWorldLimits = [-border_size, square_size + border_size];
    img2_ref.YWorldLimits = [-border_size, square_size + border_size]; 
    img2_ref.ImageSize = [square_size + 2*border_size, square_size + 2*border_size];
    [img3, img3_ref] = imwarp(img,tf,'OutputView',img2_ref); % reapply transform using the world limits we altered.
    
    
    img3 = im2uint16(img3); % convert img back to 16bit img
    outputName = [base_file_name '_' flux_readings{k} '_tform.tif']; % make new file name
    fullDestinationName = fullfile(destinationFolder_overexp,outputName); % get full file name
    imwrite(img3,fullDestinationName); % save the file
    
    if k == button
        transImg = img3;
    end
end

waitfor(msgbox('transforms complete!'));

% we will now display the original and transformed images side by side.
% close(f1); % close the original figure

f2 = figure; % create a new figure
sp(1) = subplot(1,2,1); % create the left subplot
imshow(bright_img_gamma) % display the original image
axis on
hold on
% plot the points the user selected
plot([x;x(1)],[y;y(1)],'Marker','x','MarkerSize',10,'Color','r','MarkerEdgeColor','Magenta','Linewidth',0.5); % burn points into image

title('Original Image')

sp(2) = subplot(1,2,2); % create the right subplot
transImg_gamma = imadjust(transImg,[],[],0.15); % adjust the gamma on the transformed image
imshow(transImg_gamma); % display the transformed image
axis square;
axis on

% create the transformed points
x_tf = [0; 1; 1; 0]*square_size + ones(4,1)*border_size;
y_tf = [0; 0; 1; 1]*square_size + ones(4,1)*border_size;

% plot the new points
hold on
points_plot = plot([x_tf;x_tf(1)],[y_tf;y_tf(1)],'Marker','x','MarkerSize',10,'Color','r','MarkerEdgeColor','Magenta','Linewidth',0.5); % burn points into image

title('Transformed Image')

f3 = figure;
imshow(transImg_gamma);
axis off
hold on
points_plot = plot([x_tf;x_tf(1)],[y_tf;y_tf(1)],'Marker','x','MarkerSize',10,'Color','r','MarkerEdgeColor','Magenta','Linewidth',0.5); % burn points into image
