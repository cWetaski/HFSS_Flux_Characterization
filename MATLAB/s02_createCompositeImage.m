clc;  % Clear command window.
clear;  % Delete all variables.
close all;  % Close all figure windows except those created by imtool.
imtool close all;  % Close all figure windows created by imtool.
workspace;  % Make sure the workspace panel is showing.
fontSize = 16;

% Iterate through a set of images.
% For each image, the user draws a polygon around areas to omit from
% average of flux distribution(flux sensor or visible bolt heads).
% A weighted average of pixel values is produced based on omitting
% the specified areas from the average.

% Get some dimensions
radius_alumina_hole = 1.125*25.4/2 *1.3; % mm. hole is actually 1.125" diameter, but we multiply by 1.3 to create a slightly larger roi
radius_alumina_hole_pixels = radius_alumina_hole*2000/(12*25.4); % 2000 pixels per 12 inchs
radius_bolt = 7/16*25.4/2 * 1.5; %mm. bolt head has 7/16" diameter. Multiply by 1.5 since it's harder to get the center accurate
radius_bolt_pixels = radius_bolt*2000/(12*25.4);

% select folder which contains set of N images, which we store in array
% called imgs. We will convert each image to type double. 
folder = strcat(pwd,'/images/transformed');
folder = uigetdir(folder);

% Create filename for saving composite image
folderSplit = regexp(folder,'\','split');
compositeImageName = strcat(char(folderSplit(length(folderSplit))),'_composite.tif'); 

imageFiles = dir(fullfile(folder,'*.TIF')); % get the set of images
folder_overexp = strcat(folder,'\overexposed');
imageFiles_overexp = dir(fullfile(folder_overexp,'*.TIF')); % get the set of overexposed images

image_size = size(imread(fullfile(folder,imageFiles(1).name)));
imgs = zeros(image_size(1),image_size(2),length(imageFiles)); % create array to store images
imgs_overexp = imgs; % create array to store overexposed images (it will be the same size as the array "imgs")
for k = 1 : length(imageFiles)
    imgs(:,:,k) = im2double(imread(fullfile(folder,imageFiles(k).name)));
    imgs_overexp(:,:,k) = im2double(imread(fullfile(folder_overexp,imageFiles_overexp(k).name)));
end

flux_img_sum = zeros(image_size);
flux_img_sum_overexp = zeros(image_size);
img_weights = uint8(zeros(image_size));

for n = 1:length(imageFiles)
    % display imgs(n) user draws regions to be omitted, 
    % which are returned as masks.
    % First Image
    img = imgs(:,:,n); % get the current img
    img_overexp = imgs_overexp(:,:,n); % get the overexposed version of the current image
    gamma_img = imadjust(img_overexp,[],[],0.2); % gamma adjusted img to better see dark regions
    imshow(gamma_img);
    
    
    waitfor(msgbox('select just inside the four corners of the alumina plate (counterclockwise starting from top left)'));
    [corners_x, corners_y] = ginput(4);
    
    % omit_mask will be a mask of all locations in the image which will be
    % omitted in the averaging. It will include: the area outside the 
    % alumina plate, the area near each of the 4 bolts, and the area near
    % the hole in the alumina plate.
    
    omit_mask = (~poly2mask(corners_x,corners_y,image_size(1),image_size(2))); 
    
    % show the masked image
    gamma_img(omit_mask) = 0; 
    imshow(gamma_img);
    
    
    waitfor(msgbox('select the center of the four bolts'));
    [bolt_x, bolt_y] = ginput(4);
    
    % we will make a circular mask around each bolt location selected
    num_points = 50; % circle will be 50 point polygon
    theta = (0:num_points-1)*(2*pi/num_points); % generate vector of theta values
    for k = 1:4
        bolt_circle_x = bolt_x(k) + radius_bolt_pixels*cos(theta);
        bolt_circle_y = bolt_y(k) + radius_bolt_pixels*sin(theta);
        mask = poly2mask(bolt_circle_x,bolt_circle_y,image_size(1),image_size(2)); % create the circular mask
        omit_mask = or(omit_mask,mask); % use logical OR to add the mask into omit_mask
        gamma_img(omit_mask)=0;
    end
    
    imshow(gamma_img)
    
    waitfor(msgbox('select the center of the hole in the alumina sheet'));
    [hole_x, hole_y] = ginput(1);
    
    % make a circular mask around the center of the hole in the alumina.
    num_points = 50; 
    theta = (0:num_points-1)*(2*pi/num_points);
    hole_circle_x = hole_x + radius_alumina_hole_pixels*cos(theta);
    hole_circle_y = hole_y + radius_alumina_hole_pixels*sin(theta);
    mask = poly2mask(hole_circle_x,hole_circle_y,image_size(1),image_size(2));
    omit_mask = or(omit_mask,mask);
    gamma_img(omit_mask) = 0;
    imshow(gamma_img);
    
    masked_img = img; % mask the actual img
    masked_img(omit_mask) = 0;
    
    masked_img_overexp = img_overexp; % mask the overexposed image
    masked_img_overexp(omit_mask) = 0;
    
    incl_mask = uint8(~omit_mask); % convert the mask to uint8 so we can add masks and then use the values as weights;
    flux_img_sum = flux_img_sum + masked_img; % add the masked images together
    flux_img_sum_overexp = flux_img_sum_overexp + masked_img_overexp; % add the masked overexposed images together
    img_weights = img_weights + incl_mask; % add the masks together to get weights
    
    waitfor(msgbox('continue'));
end

% We now have fluxImgSum which carries the sum total of all images after
% omitting certain areas. We also have the matrix imgWeights which carries
% the number of images which were included in each position of fluxImg. 
% e.g. If you have 3 images total, then 3s in imgWeights represent regions
% which were omitted from none of the images, 2s represent regions which
% were omitted from only one 1 image, and so on. 

% We will now build the weighted average fluxImgSum
composite_flux_img = zeros(image_size);
composite_flux_img_overexp = zeros(image_size);
for n = 1:length(imageFiles)
    % Let 'mask' equal roi where imgWeights == n
    mask = (img_weights == n);
    masked_flux_sum = flux_img_sum;
    masked_flux_sum(~mask) = 0; % only take the values of masked_flux_sum which coincide with weights which are equal to n
    composite_flux_img = composite_flux_img + masked_flux_sum/n;
    
    % and for the overexposed set
    masked_flux_sum_overexp = flux_img_sum_overexp;
    masked_flux_sum_overexp(~mask) = 0;
    composite_flux_img_overexp = composite_flux_img_overexp + masked_flux_sum_overexp/n;
end
% Save the resulting image
destinationFolder = strcat(pwd,'\images\composite\');
full_destination_file_path = strcat(destinationFolder,compositeImageName);
composite_flux_img_16 = im2uint16(composite_flux_img); % convert to uint16
imwrite(composite_flux_img_16,full_destination_file_path); % save the composite image

full_destination_file_path_overexp = strcat(destinationFolder,'overexposed\',compositeImageName);
composite_flux_img_16_overexp = im2uint16(composite_flux_img_overexp);
imwrite(composite_flux_img_16_overexp,full_destination_file_path_overexp);

average_flux_img_overexp_gamma = imadjust(composite_flux_img_overexp,[],[],0.2);
imshow(average_flux_img_overexp_gamma); % show the composite image for the overexposed set with gamma adjustment

    

        
        



