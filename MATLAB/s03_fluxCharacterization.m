clc;  % Clear command window.
clear;  % Delete all variables.
close all;  % Close all figure windows except those created by imtool.
imtool close all;  % Close all figure windows created by imtool.
workspace;  % Make sure the workspace panel is showing.
font_size = 16;

% Using a composite image, and the flux measurement values from each image
% that make up the composite image, characterize the flux.

define_sensors = false; % set to false if the sensor masks are already defined and you just want to run the code using existing saved masks
L = 2000; % This variable must be equal to the side length of the square defined in transformation in the oblique_transform script
sensor_sensitivity = 0.112; % mV/(kW/m^2) or 0.112 x 10^(-6) V/(W/m^2)


% select folder which contains set of images, which we store in 3d array
% called imgs.
folder = uigetdir(strcat(pwd,'/images/transformed')); % choose the folder with the set of images
image_files = dir(fullfile(folder,'*.TIF')); % get the set of the images from the selected folder
folder_overexp = strcat(folder,'\overexposed');
image_files_overexp = dir(fullfile(folder_overexp,'*.TIF')); % get the set of overexposed images

N = length(image_files); % get the number of images
image_size = size(imread(fullfile(folder,image_files(1).name))); % get the dimensions of the images

imgs = zeros(image_size(1),image_size(2),N); % initialize 3D matrix which stores each image
imgs_overexp = imgs; % create array to store overexposed images (it will be the same size as the array "imgs")
img_names = cell(N,1); % initialize cell array of the image file names
sensor_readings = zeros(N,1); % initialize vector which stores the sensor measurement associated with each image

for n = 1 : N % iterate through N images
    imgs(:,:,n) = imread(fullfile(folder,image_files(n).name)); % get image data
    imgs_overexp(:,:,n) = imread(fullfile(folder_overexp,image_files_overexp(n).name)); % get overexposed image data
    
     % get the filename of each img, excluding file extension
    [lin_fit,f,e] = fileparts(image_files(n).name);
    img_names{n} = fullfile(lin_fit,f);
    
    % extract the sensor measurement from the filename assuming filename includes measurement as "01.234" -> 1.234 mV
    sensor_readings(n) = str2double(regexp(img_names{n},'\d+\.+\d*','match'));
end

% get associated composite image
folderSplit = regexp(folder,'\','split');
composite_image_name = strcat(char(folderSplit(length(folderSplit))),'_composite.tif'); 
composite_folder = strcat(pwd,'\images\composite\');
composite_file_path = strcat(composite_folder,composite_image_name);
composite_img = imread(composite_file_path);

composite_folder_overexp = strcat(pwd,'\images\composite\overexposed\');
composite_file_path = strcat(composite_folder_overexp,composite_image_name);
composite_img_overexp = imread(composite_file_path);

% compute circle radius
alumina_area = 304.8^2; % mm^2 area of the alumina sheet
sensor_area = pi/4*10^2; % mm^2 area of the sensor (diameter 10 mm)

area_per_pixel = alumina_area/L^2*1/1000^2; % area of each pixel m^2 / pixel
pix_per_mm = L/(25.4*12); % pixels/mm
sensor_pixels = sensor_area*pix_per_mm^2; % number of pixels in the sensor area
sensor_radius = sqrt(sensor_pixels/pi); % calculated sensor radius (in number of pixels)

% initialize some variables used when the user picks the sensor locations
c = zeros(N,2); % stores the center of the flux sensor as defined by the user (or existing .mat files)
r = zeros(N,1); % stores radii of circles defined by the user (should be roughly equal to sensor_radius)
ravg = 0; % Will calculate the average of r, to check that it is close to sensor_radius
pixel_vals = zeros(N,1); % Will store the average pixel value of the composite image in the region of each sensor measurement.

flux_measurements = sensor_readings/sensor_sensitivity*2.446; % 2.446 is conversion factor for multimeter

if define_sensors == true % if we are defining the sensor locations, we want to visualize them 
    f1 = figure;
    imshow(imadjust(composite_img_overexp,[],[],0.4));
    zoom(2)
    hold on;
    f2 = figure;
end

for n = 1 : N % obtain sensor masks
    mask_file_name = strcat(folder, '/', img_names{n}, '_sensorMask.mat');
    if define_sensors == false % if we want to use already defined sensor masks
        try
            sensor_mask = load(mask_file_name).sensor_mask; % try loading existing sensor mask
            mask_area = sum(sensor_mask(:));
        catch e
            waitfor(msgbox('No pre-defined sensor masks found, you will have to define new sensor locations'));
            define_sensors = true; % if sensor masks don't already exist, then we must create them.
            f1 = figure;
            imshow(imadjust(composite_img_overexp,[],[],0.4));
            zoom(2)
            hold on;
            f2 = figure;
        end
    end
    if define_sensors == true
        set(0, 'currentfigure', f2)
        clf(f2);
        img = uint16(imgs_overexp(:,:,n));
        imgGamma = imadjust(img,[],[],0.4);
        imshow(imgGamma);
        zoom(2);
        hold on
        % draw 3 points on edge of flux sensor to define it
        pts(1) = drawpoint;
        pts(2) = drawpoint;
        pts(3) = drawpoint;
        p1 = pts(1).Position;
        p2 = pts(2).Position;
        p3 = pts(3).Position;

        [center radius] = def3ptCircle(p1,p2,p3); % custom function which defines circle from three points
        c(n,:) = center;
        r(n) = radius;
        ravg = ravg + r(n)/N;

        circle = drawcircle('Center',c(n,:),'Radius',sensor_radius);
        sensor_mask = createMask(circle);
        save(mask_file_name,'sensor_mask');
        num_points = 100;   
        theta = (0:num_points-1)*(2*pi/num_points);
        circle_x = c(n,2) + sensor_radius*cos(theta);
        circle_y = c(n,1) + sensor_radius*sin(theta);
        circle_poly = polyshape(circle_y,circle_x);
        plot(circle_poly,'FaceColor','[1 0 0]','EdgeColor',[0 0 0]);
        delete(circle);
        set(0, 'currentfigure', f1)
        plot(circle_poly,'FaceColor','[1 0 0]','EdgeColor',[0 0 0]);
        waitfor(msgbox('continue'));
            
    end

    flux_mask = composite_img;
    flux_mask(~sensor_mask) = 0;
    % Get the average pixel value in the region defined by the sensor mask.
    % This is obtained by dividing the sum of the pixel values in that
    % region by the number of pixels in that region. 
    pixel_vals(n) = sum(flux_mask(:))/sum(sensor_mask(:)); 
end

% get max pixel value
max_pixel = double(max(composite_img(:)));

% calculate measurement errors
flux_error = (flux_measurements.^2*0.00526 + 16.13).^(1/2); % calculated in uncertainty section
weights = 1./flux_error.^2; % weights for weighted linear regression are inverse of variance

% weighted linear regression through origin
fit_WLS_origin = fit(pixel_vals,flux_measurements,'poly1','Lower', [-Inf 0], 'Upper', [Inf, 0],'Weights',weights);
p = fit_WLS_origin.p1;
poly_bounds = confint(fit_WLS_origin); % get 95% confidence interval bounds for polynomial
p_lower = poly_bounds(1,1); % get smallest bounding polynomial
p_upper = poly_bounds(2,1); % get largest bounding polynomial

% calculate r_squared
fx = p*pixel_vals;
residuals = fx - flux_measurements; % vector of residuals
yavg = mean(flux_measurements); % mean of y values
ssr = sum(residuals.^2); % sum of the square of the resididuals
sst = sum((flux_measurements - yavg).^2); % total sum of squares
r_squared = 1 - ssr/sst; % r squared value

f1 = figure;
grid on;
hold on
errorbar(pixel_vals,flux_measurements,flux_error,'o');
axis([0 max_pixel 0 max_pixel*p_upper]);
xlabel('Pixel Value')
ylabel('Flux (kW/m^2)');

% get line of best fit
x = linspace(0,max_pixel);
y_WLS_origin = p*x;
y_upper = p_upper*x;
y_lower = p_lower*x;


% plot the correlation, in addition to the upper and lower bounds
hold on;
% plot(x,y_WLS,'k');
% plot(x,y_WLS_origin,'--','Color','r');
% plot(x,y_OLS,'--','Color','b');
plot(x,y_WLS_origin,'Color','r');
plot(x,y_upper,'--','Color','k');
plot(x,y_lower,'--','Color','k');
title('Flux vs Pixel Value')
leg = legend('Flux Measurements','Best fit','95% CI')
leg.Location = 'southeast'

% create flux map
f2 = figure;
flux_img = (im2double(composite_img))*(2^16-1)*p;
flux_img_upper = (im2double(composite_img))*(2^16-1)*p_upper;
flux_img_lower = (im2double(composite_img))*(2^16-1)*p_lower;
peak_val = max(flux_img(:));
for row=1:image_size(1)
    for col=1:image_size(2)
        if(flux_img(row,col)== peak_val)
        peak_x = col;
        peak_y = row;
        end
    end
end

% we will center the flux map on the peak flux
origin_curtain = [peak_x, peak_y]; 

% offset x and y data by location of peak flux
xdata = -origin_curtain(1) : image_size(2) - origin_curtain(1);
ydata = -origin_curtain(2) : image_size(1) - origin_curtain(2);

ticks_pixels = [-200 -100 0 100 200]*pix_per_mm;

total_flux_sc = imagesc(flux_img, 'XData', xdata, 'YData', ydata);
hold on;
cbar_curtain = colorbar;
title(cbar_curtain,'Flux (kW/m^2)')
colormap('gray')
axis equal
axis on;
xlim([-200*pix_per_mm,200*pix_per_mm]);
ylim([-200*pix_per_mm,200*pix_per_mm]);
xticks(ticks_pixels);
yticks(ticks_pixels);
% conversion of ticks to units of mm from https://www.mathworks.com/matlabcentral/answers/364462-displaying-an-axis-scale-in-mm-rather-than-pixels
addMM=@(x) sprintf('%.0f',x/pix_per_mm);
xlabels = cellfun(addMM,num2cell(xticks'),'UniformOutput',false);
ylabels = cellfun(addMM,num2cell(-yticks'),'UniformOutput',false);
for i = 1:length(ylabels)
    if strcmp(ylabels{i},'-0') % flipping the y axis results in a '-0' label, need to fix this
        ylabels{i} = '0';
    end
end
xticklabels(xlabels);
yticklabels(ylabels);
xlabel('x (mm)');
ylabel('y (mm)');
grid on;
title('Flux Map');


% obtaining a radial flux profile, centered on the peak flux value - credit
% to ImageAnalyst again https://www.mathworks.com/matlabcentral/answers/276298-how-to-plot-the-radial-profile-of-a-2d-image

% Find out what the max distance will be by computing the distance to each corner.
distanceToUL = sqrt((1-peak_y)^2 + (1-peak_x)^2);
distanceToUR = sqrt((1-peak_y)^2 + (image_size(2)-peak_x)^2);
distanceToLL = sqrt((image_size(1)-peak_y)^2 + (1-peak_x)^2);
distanceToLR= sqrt((image_size(1)-peak_y)^2 + (image_size(2)-peak_x)^2);
maxDistance = ceil(max([distanceToUL, distanceToUR, distanceToLL, distanceToLR]));

profileSums = zeros(1, maxDistance);
profileCounts = zeros(1, maxDistance);
zero_dist = 0;
count = 0;

for col = 1 : image_size(2)
	for row = 1 : image_size(1)
		thisDistance = round(sqrt((col-peak_x)^2 + (row-peak_y)^2));
		if thisDistance <= 0
			zero_dist = zero_dist + flux_img(row,col);
            count = count + 1;
            continue;
		end
		profileSums(thisDistance) = profileSums(thisDistance) + flux_img(row, col);
		profileCounts(thisDistance) = profileCounts(thisDistance) + 1;
	end
end
% Divide the sums by the counts at each distance to get the average profile
averageRadialProfile = profileSums ./ profileCounts;
zero_dist = zero_dist/count;
% Plot it.
rad_profile_fig = figure;
plot((-1*length(averageRadialProfile)):length(averageRadialProfile), [flip(averageRadialProfile),zero_dist,averageRadialProfile], 'b-', 'LineWidth', 3);
xlim([-1000 1000]);
title('Radial Flux Profile');
xlabel('Distance from peak flux (pixels)');
ylabel('Average flux (kW/m^2)');


% calculate some stats about the flux map
total_flux = sum(flux_img(:))*area_per_pixel
flux_upper = sum(flux_img_upper(:))*area_per_pixel
flux_lower = sum(flux_img_lower(:))*area_per_pixel

p_error_rel = (p_upper - p)/p;
total_flux_uncertainty = total_flux*sqrt(0.000325 + p_error_rel^2) % Determined in uncertainty section
