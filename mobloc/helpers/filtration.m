classdef filtration
    methods (Static)
        
        % Remove outliers from data grid based on IAO analysis
        function [data_grid_filtered] = apply_iao_filtration(csi, data_grid, ap_n)
            data_grid_filtered = cell(size(data_grid, 1), size(data_grid, 2), ap_n);
            
            % Apply IAO to sectors to remove outliers
            for y = 1:size(data_grid, 1)
                for x = 1:size(data_grid, 2)
                    sector_idx = data.get_sector_idx(data_grid, x, y);
            
                    for ap_idx = 1:ap_n
                        if isempty(sector_idx)
                            data_grid_filtered(y, x, ap_idx) = {[]};
                        else
                            sector_csi = csi(sector_idx, :, :, :);
                            sector_offsets = sigproc.estimate_iao(sector_csi, false);
                            sector_offsets_ap = sector_offsets(:, :, ap_idx);
                            sector_good_idx = filtration.evaluate_iao_outliers(sector_offsets_ap, ap_idx, 10, 1, false);
                            sector_idx_filtered = sector_idx(sector_good_idx);
                            data_grid_filtered(y, x, ap_idx) = {sector_idx_filtered};
                        end
                    end
                end
            end
        end

        function [good_idx] = evaluate_iao_outliers(offsets, ap_idx, window_size, window_step, plot_outliers)
            window_best_mean = Inf;
            window_best_eval = [];
            
            % Use a shifting window to go through all possible base idx
            for window_start = 1:window_step:(size(offsets, 1)-window_size)
                window_idx = window_start:window_start+window_size-1;
                eval = filtration.analyze_outliers(offsets, window_idx);
                window_mean = mean(eval);
            
                if window_best_mean > window_mean && window_mean ~= 0
                    window_best_mean = window_mean;
                    window_best_eval = eval;
                end
            end
            
            if plot_outliers
                figure
                hold on
                plot(1:length(window_best_eval), window_best_eval);
                plot([1, length(window_best_eval)], [window_best_mean, window_best_mean])
                hold off
                title(strcat("AP", num2str(ap_idx),": IAO L2PCA Outlier Analysis"));
            end
        
            % Determine which samples to drop
            good_idx = window_best_eval < window_best_mean;
        end

        function [eval] = analyze_outliers(samples, base_sample_idx)
            X = samples(base_sample_idx, :);
            r = 1;
            [P, ~, ~] = filtration.mSVD(X, r);
        
            eval = zeros(size(samples, 1), 1);
            for i = 1:length(eval)
                Xn = samples(i, :);
                eval(i) = norm((Xn - P'*P*Xn), 2);
            end
        end
        
        function [U, S, V] = mSVD(X, r)
            [UK, SK, VK] = svd(X,'econ');
            S = SK(1:r, 1:r);
            U = UK(:, 1:r);
            V = VK(:, 1:r);
        end

        function [csi_filtered] = apply_csi_filter_tap(csi, threshold)
            csi_filtered = zeros(size(csi));
            for packet_idx = 1:size(csi_filtered, 1)
                for ant_idx = 1:size(csi_filtered, 3)
                    for ap_idx = 1:size(csi_filtered, 4)
                        csi_packet_raw = permute(squeeze(csi(packet_idx, :, ant_idx, ap_idx)), [2, 1]);
                        csi_filtered(packet_idx, :, ant_idx, ap_idx) = sigproc.apply_tap_filter(csi_packet_raw, threshold, false);
                    end
                end
            end
        end

    end
end