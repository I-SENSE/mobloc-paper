classdef dloc_split
    methods (Static)

        function [csi_train, csi_test, labels_train, labels_test, ap_coords, freq_center, freq, ant_d, ap_n, cdf_title] = prepare_data_test_figure10b(DATASETS, SPLIT_IDX)
            % Test case: Figure 10(b)
            % Train: dataset_jacobs_July28 -- looks like this dataset is used entirely for training
            % Train: dataset_non_fov_train_jacobs_July28_2
            % Train: dataset_fov_train_jacobs_July28_2
            % Test: dataset_fov_test_jacobs_July28_2
            % Test: dataset_non_fov_test_jacobs_July28_2
        
            disp('..Load datasets & indexes.');
            [csi_1, labels_1, ap_coords, freq, ~, ant_d, ap_n] = data.load_dataset(DATASETS(2));
            [csi_2, labels_2, ~, ~, ~, ~, ~] = data.load_dataset(DATASETS(3));
        
            [fov_test_idx_2, fov_train_idx_2, non_fov_test_idx_2, non_fov_train_idx_2] = data.load_split_idx(SPLIT_IDX(3));
        
            % For some reason dimensions are flipped. Needs to be fixed.
            non_fov_test_idx_2 = permute(non_fov_test_idx_2, [2, 1]);
            non_fov_train_idx_2 = permute(non_fov_train_idx_2, [2, 1]);
        
            % Calculate center frequency
            freq_center = data.calculate_center_freq(freq);
        
            disp('..Perform test/train splitting.');
            
            % Extract CSI data for training
            train_idx_2 = [fov_train_idx_2; non_fov_train_idx_2];
            [csi_train, labels_train] = data.merge_datasets(csi_1, csi_2(train_idx_2, :, :, :), labels_1, labels_2(train_idx_2, :));
            
            % Extract CSI data for testing
            test_idx_2 = [fov_test_idx_2; non_fov_test_idx_2];
            csi_test = csi_2(test_idx_2, :, :, :);
            labels_test = labels_2(test_idx_2, :);
        
            cdf_title = 'DLoc Evaluation. Figure 10(b): complex environment.';
        end
        
        function [csi_train, csi_test, labels_train, labels_test, ap_coords, freq_center, freq, ant_d, ap_n, cdf_title] = prepare_data_test_figure10a(DATASETS, SPLIT_IDX)
            % Test case: Figure 10(a)
            % Train: dataset_non_fov_train_July18
            % Train: dataset_fov_train_July18
            % Test: dataset_non_fov_test_July18
            % Test: dataset_fov_test_July18
        
            disp('..Load datasets & indexes.');
            [csi, labels, ap_coords, freq, ~, ant_d, ap_n] = data.load_dataset(DATASETS(1));
            [fov_test_idx, fov_train_idx, non_fov_test_idx, non_fov_train_idx] = data.load_split_idx(SPLIT_IDX(1));
            
            % For some reason dimensions are flipped. Needs to be fixed.
            non_fov_test_idx = permute(non_fov_test_idx, [2, 1]);
            non_fov_train_idx = permute(non_fov_train_idx, [2, 1]);
            
            % Calculate center frequency
            freq_center = data.calculate_center_freq(freq);
            
            disp('..Perform test/train splitting.');
            % Merge fov & non-fov splits as they are used together for this case.
            train_idx = [fov_train_idx; non_fov_train_idx];
            test_idx = [fov_test_idx; non_fov_test_idx];
            
            csi_train = csi(train_idx, :, :, :);
            csi_test = csi(test_idx, :, :, :);
            labels_train = labels(train_idx, :);
            labels_test = labels(test_idx, :);
        
            cdf_title = 'DLoc Evaluation. Figure 10(a): simple environment.';
        end

    end
end