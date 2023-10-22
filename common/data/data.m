classdef data
    methods (Static)

        % Extracts all valuables from a given file in DLoc dataset
        function [csi, labels, ap, freq, lambda, ant_d, ap_n] = load_dataset(filepath)
            S = load(filepath);

            opt = getfield(S, 'opt');
            ap = getfield(S, 'ap');
            csi = getfield(S, 'channels');
            labels = getfield(S, 'labels');
        
            freq = getfield(opt, 'freq');
            lambda = getfield(opt, 'lambda');
            ant_d = getfield(opt, 'ant_sep');

            ap_n = size(csi, 4);
        end

        function [fov_test_idx, fov_train_idx, non_fov_test_idx, non_fov_train_idx] = load_split_idx(filepath)
            S = load(filepath);

            fov_test_idx = getfield(S, 'fov_test_idx');
            fov_train_idx = getfield(S, 'fov_train_idx');
            non_fov_test_idx = getfield(S, 'non_fov_test_idx');
            non_fov_train_idx = getfield(S, 'non_fov_train_idx');
        end

        % Removes redundant datapoints from the dataset (outside of DoA boundaries)
        function [csi_filtered, labels_filtered] = filter_dataset(csi, labels, ap_coords)
            idx_filtered = [];
            for packet_idx = 1:size(labels, 1)
                doa_ap_1 = data.ground_aoa(cell2mat(ap_coords(1, 1)), labels(packet_idx, :), 1);
                doa_ap_2 = data.ground_aoa(cell2mat(ap_coords(1, 2)), labels(packet_idx, :), 2);
                doa_ap_3 = data.ground_aoa(cell2mat(ap_coords(1, 3)), labels(packet_idx, :), 3);
                doa_ap_4 = data.ground_aoa(cell2mat(ap_coords(1, 4)), labels(packet_idx, :), 4);

                if data.is_doa_valid(doa_ap_1) && data.is_doa_valid(doa_ap_2) && data.is_doa_valid(doa_ap_3) && data.is_doa_valid(doa_ap_4)
                    idx_filtered = [idx_filtered; packet_idx];
                end
            end

            csi_filtered = csi(idx_filtered, :, :, :);
            labels_filtered = labels(idx_filtered, :);
        end

        % Plots given labels & AP locations on the map
        function [] = plot_dataset(ap_coords, labels, plot_figure, color)
            if plot_figure
                figure;
            end

            hold on
            
            for i = 1:size(labels, 1)
                plot(labels(i, 1), labels(i, 2), color);
            end

            colors = {'redo', 'greeno', 'magentao', 'blueo'};

            for i = 1:size(ap_coords, 2)
                AP = cell2mat(ap_coords(1, i));
                plot(AP(1, 1), AP(1, 2), string(colors(i)));
            end
            
            xlim([-1, 25]);
            ylim([-1, 16]);
            hold off
        end

        % Find IDs of data points for each of the sectors in the area
        function [data_grid] = find_sectors(labels, area_height_m, area_width_m, grid_height, grid_width)            
            sector_height = area_height_m / grid_height;
            sector_width = area_width_m / grid_width;
            data_grid = cell(grid_height, grid_width);
            
            for x_i = 1:grid_width
                for y_i = 1:grid_height
                    sector_idx = data.find_sector_data_idx(labels, ...
                        sector_width * (x_i - 1), sector_width * x_i, ...
                        sector_height * (y_i - 1), sector_height * y_i);
                    data_grid(y_i, x_i) = {sector_idx};
                end
            end
        end

        % Retrieves indexes of data points for a certain grid sector 
        % (x - columns, y - rows)
        function [sector_idx] = get_sector_idx(data_grid, sector_x, sector_y)
            try
                sector_idx = cell2mat(data_grid(sector_y, sector_x));
            catch
                sector_idx = [];
            end
        end

        function [] = plot_sector_counts(grid)
            counts = zeros(size(grid, 1), size(grid, 2));
            for i = 1:size(counts, 1)
                for j = 1:size(counts, 2)
                    items = data.get_sector_idx(grid, j, i);
                    counts(i, j) = length(items);
                end
            end
        
            figure
            h = heatmap(counts);
            h.FontSize = 10;
            h.NodeChildren(3).YDir='normal';      
        end

        % Checks if calculated DoA is within permitted boundaries
        function [valid] = is_doa_valid(doa)
            valid = 25 <= doa && doa <= (180 - 25);
        end

        % Calculates center frequency from array of frequencies for each
        % subcarrier
        function [freq_center] = calculate_center_freq(freq)
            freq_center = (freq(1, 117) + freq(1, 118))/2;
        end

        % Calculates ground truth angle
        function [aoa] = ground_aoa(ap, p, ap_idx)
            ap_x = ap(1, 1);
            ap_y = ap(1, 2);
            p_x = p(1, 1);
            p_y = p(1, 2);

            % Process cases when point is in the second quarter of the circle
            switch ap_idx
                case 1
                    aoa = rad2deg(atan(abs(ap_x - p_x) / abs(ap_y - p_y)));
                    if ap_y > p_y
                        aoa = 180 - aoa;
                    end
                case 2
                    aoa = rad2deg(atan(abs(ap_x - p_x) / abs(ap_y - p_y)));
                    if ap_y < p_y
                        aoa = 180 - aoa;
                    end
                case 3
                    aoa = rad2deg(atan(abs(ap_y - p_y) / abs(ap_x - p_x)));
                    if ap_x > p_x
                        aoa = 180 - aoa;
                    end
                case 4
                    aoa = rad2deg(atan(abs(ap_y - p_y) / abs(ap_x - p_x)));
                    if ap_x < p_x
                        aoa = 180 - aoa;
                    end
            end
        end
    
        % Finds all data points that are located within a certain sector
        function [sector_idx] = find_sector_data_idx(labels, x_min, x_max, y_min, y_max)
            sector_idx = [];
            for label_i = 1:size(labels, 1)
                point = labels(label_i, :);
        
                if (x_min <= point(1) && point(1) <= x_max && y_min <= point(2) && point(2) <= y_max)
                    sector_idx = [sector_idx; label_i];
                end
            end
        end

        function [availability_map] = build_availability_map(grid)
            availability_map = zeros(size(grid));
            for y = 1:size(grid, 1)
                for x = 1:size(grid, 2)
                    sector_idx = data.get_sector_idx(grid, x, y);
                    
                    if ~isempty(sector_idx)
                        availability_map(y, x) = 1;
                    end
                end
            end
        end

        function [availability_map] = build_availability_map_thresh(grid, threshold)
            availability_map = zeros(size(grid));
            for y = 1:size(grid, 1)
                for x = 1:size(grid, 2)
                    sector_idx = data.get_sector_idx(grid, x, y);
                    
                    if length(sector_idx) >= threshold
                        availability_map(y, x) = 1;
                    end
                end
            end
        end

        % Perform test/train split of the dataset
        function [X_train, X_test, y_train, y_test] = train_test_split(X, y, test_size, random_state)
            % If split info isn't provided -- take 1/3 of data for testing
            if nargin<3
                test_size = 1/3;
            end
            
            if nargin>3
                rng('default');
                rng(random_state);
            end
            
            n = size(X,1);
            n_test = floor(test_size*n);
            n_train = n - n_test;
            
            indp = randperm(n);
            ind_train = indp(1:n_train);
            ind_test = indp(n_train+1:n);
            
            X_train = X(ind_train,:,:,:);
            X_test = X(ind_test,:,:,:);
            y_train = y(ind_train,:);
            y_test = y(ind_test,:);
        end

        function [csi_merged, labels_merged] = merge_datasets(csi_1, csi_2, labels_1, labels_2)
            csi_merged = zeros(size(csi_1, 1) + size(csi_2, 1), size(csi_1, 2), size(csi_1, 3), size(csi_1, 4));
            csi_merged(1:size(csi_1, 1), :, :, :) = csi_1;
            csi_merged(size(csi_1, 1) + 1:size(csi_merged, 1), :, :, :) = csi_2;
        
            labels_merged = zeros(size(labels_1, 1) + size(labels_2, 1), size(labels_1, 2));
            labels_merged(1:size(labels_1, 1), :) = labels_1;
            labels_merged(size(labels_1, 1) + 1:size(labels_merged, 1), :) = labels_2;
        end

        % Convert grid with labels from 1-sensor to N-sensor format
        function [data_grid_ready] = adapt_grid(data_grid, ap_n)
            data_grid_ready = cell(size(data_grid, 1), size(data_grid, 2), ap_n);
            for ap_i = 1:ap_n
                data_grid_ready(:, :, ap_i) = data_grid;
            end
        end

        % Find mutually available sectors from two availability grids
        % (values must be 0 or 1)
        function [grid_common] = find_common_sectors(grid1, grid2)
            grid_common = zeros(size(grid1));
            for y = 1:size(grid_common, 1)
                for x = 1:size(grid_common, 2)
                    grid_common(y, x) = grid1(y, x) && grid2(y, x);
                end
            end
        end
    end
end