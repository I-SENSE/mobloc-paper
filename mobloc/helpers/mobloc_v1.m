% This script contains functions for implementing MobLoc v1 fingerprinting
% methods based on 1D pseudospectrum of MUSIC.
classdef mobloc_v1
    methods (Static)
        function [rho] = find_kernel_coefficient(fpdb_offline, grid_availability, tile_height, tile_width, rho_min, rho_max, rho_step)
            % Runs a leave-one-out cross validation algorithm to find the most
            % optimal value of kernel regression coefficient 'rho'.
        
            % Generate possible parameter values
            rho_values = rho_min:rho_step:rho_max;
            rho_value_count = length(rho_values);
        
            % Prepare array for storing errors for each rho value
            rho_errors = zeros(length(rho_values), 1);
        
            % Run leave-one-out cross validation for each rho value
            for rho_idx = 1:length(rho_values)
                disp(strcat('Running tests for rho: [', int2str(rho_idx), '/', int2str(rho_value_count), ']'));
                rho = rho_values(rho_idx);
                rho_error_values = zeros(0, 1);
                for i = 1:size(fpdb_offline, 1)
                    for j = 1:size(fpdb_offline, 2)
                        if grid_availability(i, j) == 1
                            % Calculate ground truth location of the tile
                            grid_cell_x = j * tile_width - tile_width/2;
                            grid_cell_y = i * tile_height - tile_height/2;
        
                            % Retrieve online FP
                            fp_music = squeeze(fpdb_offline(i, j, :, :, :));
            
                            % Mark tile as unavailable
                            grid_avail_local = grid_availability;
                            grid_avail_local(i, j) = 0;
            
                            % Run localization
                            map_d_mobloc = mobloc_v1.run_fp_matching(grid_avail_local, fpdb_offline, fp_music, false, '', false);
            
                            % Calculate kernel coefficient & use it for weighted centroid
                            k_grid = localization.calculate_kernel_function(map_d_mobloc, grid_avail_local, rho, false);
                            centroid_coord = localization.calculate_centroid(k_grid, tile_height, tile_width);
                    
                            % Calculate Gaussian distance between predicted & real coordinates
                            err = sqrt((centroid_coord(1) - grid_cell_x)^2 + (centroid_coord(2) - grid_cell_y)^2);
            
                            % Save error
                            rho_error_values(end+1, 1) = err;
                        end
                    end
                end
        
                rho_errors(rho_idx, 1) = mean(rho_error_values, "omitnan");
            end
        
            % Plot errors
            figure;
            hold on;
            plot(rho_values, rho_errors, 'blue');
            scatter(rho_values, rho_errors, 'blue');
            hold off;
            title('Localization error vs. Rho values');
        
            % Find rho associated with the smallest mean localization error
            [~, err_idx] = min(rho_errors, [], 'omitnan');
        
            rho = rho_values(err_idx);
        end
        
        % Perform localization testing based on offline FPDB and a list of
        % prepared online fingerprints & their labels
        function [errors, cdf_x, cdf_y] = test_localization(fpdb_offline, fpdb_online_fps, fpdb_online_labels, tile_height, tile_width, grid_availability, label, rho)
            errors = zeros(size(fpdb_online_fps, 1), 1);
            for test_idx = 1:size(fpdb_online_labels, 1)
                % Retrieve fingerprit, label & grid index
                fp_online = squeeze(fpdb_online_fps(test_idx, :, :, :));
                fp_online_label = squeeze(fpdb_online_labels(test_idx, :));
                % fp_online_label_grid = squeeze(fpdb_online_labels_grid(test_idx, :));
            
                % Localize fingerprint
                sector_result = mobloc_v1.run_fp_matching(grid_availability, fpdb_offline, fp_online, false, '', false);

                k_grid = localization.calculate_kernel_function(sector_result, grid_availability, rho, false);
                centroid_coord = localization.calculate_centroid(k_grid, tile_height, tile_width);

                % Calculate Gaussian distance between predicted & real coordinates
                dist = sqrt((centroid_coord(1) - fp_online_label(1))^2 + (centroid_coord(2) - fp_online_label(2))^2);

                % Save into errors array
                errors(test_idx, 1) = dist;
            end
            
            % Plot CDF
            figure;
            h = cdfplot(errors);
            hold on;
            title(label);
            xlim([0, 5]);
            ylim([0, 1]);
            hold off;

            % Extract XY values of the plot
            cdf_x = get(h, 'Xdata');
            cdf_y = get(h, 'Ydata');
        end

        % Performs localization of a specified sector using online &
        % offline FPDBs (for a merged fingerprint)
        function [localization_map] = run_fp_matching(grid_availability, fpdb_offline, fp_online, apply_blur, label, plot_map)
            test_sector_online_fp = mobloc_v1.merge_fps(fp_online);
            
            % Localization map dimensions: (map y, map x, aps (comparison for each ap))
            localization_map = zeros(size(fpdb_offline, 1), size(fpdb_offline, 2));
            
            % Go through localization map and record comparisons of FPs for all APs
            for y = 1:size(localization_map, 1)
                for x = 1:size(localization_map, 2)
                    if grid_availability(y, x) == 1
                        fp_offline = mobloc_v1.merge_fps(squeeze(fpdb_offline(y, x, :, :, :)));
                        res = mobloc_v1.euclid_weighted(fp_offline, test_sector_online_fp);
                        localization_map(y, x) = res;
                    end
                end
            end

            localization_map(~grid_availability) = NaN;
            
            % Apply Gaussian blurring
            if apply_blur
                blur_factor = 0.8;
                localization_map = imgaussfilt(localization_map, blur_factor);
            end
            
            % Plot heatmap
            if plot_map
                figure
                h = heatmap(localization_map);
                colormap(autumn)
                h.FontSize = 10;
                h.NodeChildren(3).YDir='normal';
                title(label);
            end
        end

        % This function performs the following operations:
        % 1. Runs MUSIC on samples within the specified range, compiles into a single array
        % 2. For each subcarrier, estimates PDF (probability density function)
        % 3. For each PDF, takes out the variation, and a maximum, stores in a fingerprint array
        function [fingerprint] = build_fp(csi, ant_d, freq)
            if size(csi, 1) == 0
                return;
            end

            ant_n = size(csi, 3);

            % Sector pseudo dims: (packets, angles, aps)
            sector_pseudo = zeros(size(csi, 1), 181, size(csi, 4));

            % For each AP, run MUSIC for each available packet
            for ap_idx = 1:size(csi, 4)
                % Initial CSI shape: [packets, subcarriers, antennas]
                csi_ap = permute(squeeze(csi(:, :, :, ap_idx)), [1, 3, 2]);
            
                % Run 1D MUSIC to obtain pseudospectrums for each packet
                for packet_idx = 1:size(csi_ap, 1)
                    X = squeeze(csi_ap(packet_idx, :, :));
                    [~, music_pseudo] = music.run_music(X, ant_d, freq, ant_n);
                    sector_pseudo(packet_idx, :, ap_idx) = music_pseudo;
                end

                sector_pseudo(:, :, ap_idx) = rescale(sector_pseudo(:, :, ap_idx), 0, 1);
            end

            % Prepare fingerprint
            fingerprint = zeros(2, 181, size(csi, 4));
            
            % For each AP, estimate PDFs for each subcarrier
            for ap_idx = 1:size(csi, 4)
                for subcarr_idx = 1:181
                    % Run KDE on P for subcarrier, for a given range of possible P values
                    p_subcarr = sector_pseudo(:, subcarr_idx, ap_idx);
                    p_range = 0:0.01:1;
                    f = ksdensity(p_subcarr, p_range)/1000;
            
                    % Find value of the peak
                    [~, f_max_idx] = max(f);
            
                    % Find peak for the fingerprint
                    fingerprint(1, subcarr_idx, ap_idx) = p_range(f_max_idx);
                end

                sector_ap_std = var(sector_pseudo(:, :, ap_idx), 0, 1) + 0.001;
                sector_ap_weights = 1 ./ sqrt(sector_ap_std);
                fingerprint(2, :, ap_idx) = sector_ap_weights;

                % Note: re-scaling is performed across all sectors outside
                % of this function.
            end
        end

        % Performs comparison of two fingerprints using Euclidean distance
        % Note: the algorithm uses weights from first FP for calculations
        % Note: fingerprint must have dimensions 2x181
        % Ref: https://math.stackexchange.com/a/917133
        function [distance] = euclid_weighted(fp1, fp2)
            euclid_total = 0;

            for i = 1:size(fp1, 2)
                fp1_p = fp1(1, i);
                fp2_p = fp2(1, i);
                weight = mean([fp1(2, i), fp2(2, i)]); % THIS IS TEMPORARY
                euclid_total = euclid_total + (weight * (fp1_p - fp2_p)) ^ 2;
            end

            if euclid_total == 0
                distance = 0;
            else
                distance = sqrt(euclid_total);
            end
        end

        % Builds database of fingerprints based on offline dataset
        % Dimensions: (grid y, grid x, 2, 181, aps)
        function [fp_database] = build_fp_database(grid_availability, csi_dataset, data_grid_offline, ap_n, ant_d, freq)
            fp_database = zeros(size(data_grid_offline, 1), size(data_grid_offline, 2), 2, 181, ap_n);
            for y = 1:size(fp_database, 1)
                for x = 1:size(fp_database, 2)
                    if grid_availability(y, x) == 1
                        for ap_i = 1:ap_n
                            sector_idx = data.get_sector_idx(data_grid_offline(:, :, ap_i), x, y);
                            csi_sector = csi_dataset(sector_idx, :, :, ap_i);
                            fingerprint = mobloc_v1.build_fp(csi_sector, ant_d, freq);
                            fp_database(y, x, :, :, ap_i) = fingerprint;
                        end
                    end
                end
            end
            
            % Normalize inverse weights across all sectors
            fp_database(:, :, 2, :, :) = rescale(fp_database(:, :, 2, :, :), 0, 1);
        end

        function [fpdb_online_fps, fpdb_online_labels, fpdb_online_labels_grid] = build_test_fingerprints(csi, labels, data_grid_test, grid_availability, ap_n, ant_d, freq, fp_packet_count)
            % Build online fingerprints using a pre-set number of packets.
            fpdb_online_fps = zeros(0, 2, 181, ap_n);
            fpdb_online_labels = zeros(0, 2);
            fpdb_online_labels_grid = zeros(0, 2);
            
            for y = 1:size(data_grid_test, 1)
                for x = 1:size(data_grid_test, 2)
                    if grid_availability(y, x) == 1
                        % Get indexes of test samples within this sector
                        sector_idx = data.get_sector_idx(data_grid_test, x, y);
            
                        % Calculate how many fingerprints can we build
                        sector_fp_count = floor(length(sector_idx) / fp_packet_count);
            
                        % Manually cap amount of tile test cases to 5 (to
                        % reduce bias from tiles with more data)
                        if sector_fp_count > 10
                            sector_fp_count = 10;
                        end

                        % Go ahead and build these FPs
                        for fp_idx = 1:sector_fp_count
                            start_idx = (fp_idx - 1) * fp_packet_count + 1;
                            end_idx = start_idx + fp_packet_count - 1;
            
                            % Indexes of samples in the dataset for this fingerprint
                            sector_sample_idx = sector_idx(start_idx : end_idx);
            
                            % Retrieve labels for samples included in this fingerprint
                            % Then, take a mean of all labels (to find the
                            % center for all coordinates)
                            sector_labels = labels(sector_sample_idx, :);
                            sector_labels_mean = mean(sector_labels, 1);
            
                            % Build FPs separately for each AP involved
                            sector_fp = zeros(2, 181, ap_n);
                            for ap_idx = 1:ap_n
                                sector_fp_csi = csi(sector_sample_idx, :, :, ap_idx);
                                fingerprint = mobloc_v1.build_fp(sector_fp_csi, ant_d, freq);
                                sector_fp(:, :, ap_idx) = fingerprint;
                            end
            
                            % Save FP into the list
                            fpdb_online_fps(end + 1, :, :, :) = sector_fp;
                            fpdb_online_labels(end + 1, :) = sector_labels_mean;
                            fpdb_online_labels_grid(end + 1, :) = [y, x];
                        end
                    end
                end
            end
            
            % Normalize inverse weights across all sectors
            fpdb_online_fps(:, 2, :, :) = rescale(fpdb_online_fps(:, 2, :, :), 0, 1);
        end

        function [fp_merged] = merge_fps(fps)
            fp_len = size(fps, 2);
            ap_n = size(fps, 3);
        
            fp_merged = zeros(2, fp_len * ap_n);
            fp_idx_start = 1;
            for ap_idx = 1:ap_n
                fp_iter_start = fp_idx_start;
                fp_iter_end = fp_iter_start + fp_len;
        
                fp_merged(:, fp_iter_start:fp_iter_end-1) = fps(:, :, ap_idx);
        
                fp_idx_start = fp_idx_start + fp_len;
            end
        end
    end
end