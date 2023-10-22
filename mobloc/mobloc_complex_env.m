% This script is designed to run MobLoc tests on DLoc-designed test cases.

clear; close all;

addpath('../common/data');
addpath('../common/localization');
addpath('../common/music');
addpath('./helpers');

% Grid configurations for complex environment
multiploodier = 3;
area_height_m = 10;
area_width_m = 14;
grid_height = area_height_m * multiplier;
grid_width = area_width_m * multiplier;
tile_height = area_height_m / grid_height;
tile_width = area_width_m / grid_width;

% Available Datasets
DATASETS = [
    "../data/channels/channels_July18.mat";
    "../data/channels/channels_jacobs_July28.mat";
    "../data/channels/channels_jacobs_July28_2.mat";
    "../data/channels/channels_jacobs_Aug16_3.mat";
    "../data/channels/channels_jacobs_Aug16_4_ref.mat";
    "../data/channels/channels_jacobs_Aug16_1.mat"];

% Available split indexes provided by DLoc team
SPLIT_IDX = [
    "../data/split_idx/data_split_ids_July18.mat";
    "../data/split_idx/data_split_ids_jacobs_July28.mat";
    "../data/split_idx/data_split_ids_jacobs_July28_2.mat";
    "../data/split_idx/data_split_ids_jacobs_Aug16_3";
    "../data/split_idx/data_split_ids_jacobs_Aug16_4_ref";
    "../data/split_idx/data_split_ids_jacobs_Aug16_1"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('1. Prepare dataset for evaluation.');
[csi_train, csi_test, labels_train, labels_test, ap_coords, freq_center, freq, ant_d, ap_n, cdf_title] = dloc_split.prepare_data_test_figure10b(DATASETS, SPLIT_IDX);

disp('2. Plotting labels for split datasets.');
figure
data.plot_dataset(ap_coords, labels_train, false, 'black.');
data.plot_dataset(ap_coords, labels_test, false, 'black.');

disp('3. Group data into physical "sectors"');
data_grid_train = data.find_sectors(labels_train, area_height_m, area_width_m, grid_height, grid_width);
data_grid_test = data.find_sectors(labels_test, area_height_m, area_width_m, grid_height, grid_width);

disp('4. Build a supporting map to show which sectors have data');
MIN_PACKETS_OFFLINE = 20;
grid_availability_train = data.build_availability_map_thresh(data_grid_train, MIN_PACKETS_OFFLINE);
grid_availability_test = data.build_availability_map_thresh(data_grid_test, MIN_PACKETS_OFFLINE);

disp('5. Adapt unfiltered grid for further computation (to compare with/without IAO performance)');
data_grid_train_ready = data.adapt_grid(data_grid_train, ap_n);

disp('6. Build fingerprints (MUSIC-based)');
fpdb_offline = mobloc_v1.build_fp_database(grid_availability_train, csi_train, data_grid_train_ready, ap_n, ant_d, freq_center);

disp('7. Build test fingerprints');
MIN_PACKETS_ONLINE = 2;
[fpdb_online_fps, fpdb_online_labels, ~] = mobloc_v1.build_test_fingerprints(csi_test, labels_test, data_grid_test, grid_availability_test, ap_n, ant_d, freq_center, MIN_PACKETS_ONLINE);

disp('8. Find most optimal kernel coefficient');
% rho = mobloc_v1.find_kernel_coefficient(fpdb_offline, grid_availability_train, tile_height, tile_width, 5, 18, 1);
rho = 11;

disp('9. Perform evaluation');
[errors, cdf_x, cdf_y] = mobloc_v1.test_localization(fpdb_offline, fpdb_online_fps, fpdb_online_labels, tile_height, tile_width, grid_availability_test, cdf_title, rho);
