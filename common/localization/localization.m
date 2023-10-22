classdef localization
    methods (Static)
        function [x, y] = grid_idx_to_loc(grid_x, grid_y, tile_height, tile_width)
            y = (grid_y - 1) * tile_height + 0.5 * tile_height;
            x = (grid_x - 1) * tile_width + 0.5 * tile_width;
        end

        function [test_grid] = test_entire_grid_accuracy(grid_availability, fpdb_offline, fpdb_online, plot_result, plot_title)
            test_grid = zeros(size(grid_availability));
            for y = 1:size(test_grid, 1)
                for x = 1:size(test_grid, 2)
                    % Make sure there's data for this sector
                    if grid_availability(y, x) == 1
                        fp_online = squeeze(fpdb_online(y, x, :, :, :));
       
                        % Run FP matching across the whole offline FPDB
                        sector_result = localization.test_localization_merged(grid_availability, fpdb_offline, fp_online, false, '', false);
            
                        % Determine the most likely sector based on fingerprinting
                        [y_res, x_res] = localization.get_x_min(sector_result, 1); 
            
                        % If the sector matches ground truth -- set 1
                        if y_res == y && x_res == x
                            test_grid(y, x) = 1;
                        else % Otherwise, calculate Euclid distance to ground truth
                            test_grid(y, x) = sqrt((y-y_res)^2 + (x-x_res)^2);
                        end
                    else % Otherwise, set sector to NaN (so that color scheme wouldn't interfere with real data)
                        test_grid(y, x) = NaN;
                    end
                end
            end
            
            % Plot data if necessary
            if plot_result
                f = figure;
                f.Position = [0 0 1000 700];
                h = heatmap(test_grid);
                h.FontSize = 10;
                h.NodeChildren(3).YDir='normal';  
                title(plot_title);
            end
        end

        function [grid] = grid_trim_above_x(grid, x)
            grid(grid >= x) = NaN;

            figure
            h = heatmap(grid, 'ColorLimits',[0 10]);
            colormap(autumn)
            h.FontSize = 10;
            h.NodeChildren(3).YDir='normal';
        end

        % Calculate kernel function using FP matching distances & kernel coefficient
        function [k_grid] = calculate_kernel_function(grid_loc_result, grid_availability, k_coef, plot_grid)
            k_grid = zeros(size(grid_availability));
            
            for i = 1:size(k_grid, 1)
                for j = 1:size(k_grid, 2)
                    if grid_availability(i, j) == 1
                        d = grid_loc_result(i, j);
                        k_grid(i, j) = exp(-1 * k_coef * d);
                    else
                        k_grid(i, j) = NaN;
                    end
                end
            end
            
            if plot_grid
                figure
                h = heatmap(k_grid);
                h.FontSize = 10;
                h.NodeChildren(3).YDir='normal'; 
            end
        end

        % Calculate weighted centroid using kernel function grid
        function [centroid_coord] = calculate_centroid(k_grid, tile_height, tile_width)
            centroid_coord = zeros(2, 1);
            
            for i = 1:size(k_grid, 1)
                for j = 1:size(k_grid, 2)
                    if ~isnan(k_grid(i, j))
                        % Convert tile idx to coordinates
                        [x_tile, y_tile] = localization.grid_idx_to_loc(j, i, tile_height, tile_width);
            
                        % Get kernel coefficient
                        k = k_grid(i, j);
            
                        % Apply iteration to coordinates
                        centroid_coord(1, 1) = centroid_coord(1, 1) + x_tile * k;
                        centroid_coord(2, 1) = centroid_coord(2, 1) + y_tile * k;
                    end
                end
            end
            
            centroid_coord = centroid_coord / sum(k_grid, 'all', 'omitnan');
        end

        % Retrieves indexes of X smallest elements in the 2D array
        function [row, col] = get_x_min(arr, x)
            [~, idx] = mink(arr(:), x);
            [row, col] = ind2sub(size(arr), idx);
        end
    end
end