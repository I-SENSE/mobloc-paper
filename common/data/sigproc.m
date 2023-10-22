classdef sigproc
    methods (Static)
        function [offsets_full] = estimate_iao(csi, plot_iao)
            offsets_full = zeros(size(csi, 1), 3, size(csi, 4));
        
            % Loop through all APs (for each there will be a plot
            for ap_i = 1:size(csi, 4)
                % Prepare arrays to store offsets
                offsets = zeros(size(csi, 1), 3);
            
                % Loop through all packets and calculate offsets
                for packet_i = 1:size(offsets, 1)
                    packet_csi = squeeze(csi(packet_i, :, :, ap_i));
                    [offset_a12, offset_a13, offset_a14] = sigproc.calculate_antenna_offsets(packet_csi);
                    offsets(packet_i, 1) = offset_a12;
                    offsets(packet_i, 2) = offset_a13;
                    offsets(packet_i, 3) = offset_a14;
                end
            
                offsets(:, 1) = unwrap(offsets(:, 1));
                offsets(:, 2) = unwrap(offsets(:, 2));
                offsets(:, 3) = unwrap(offsets(:, 3));

                % Save offsets for returning from function
                offsets_full(:, :, ap_i) = offsets;
            
                % Plot offsets
                if plot_iao
                    figure
                    hold on;
                    plot(1:size(offsets, 1), offsets(:, 1), 'red.');
                    plot(1:size(offsets, 1), offsets(:, 2), 'blue.');
                    plot(1:size(offsets, 1), offsets(:, 3), 'green.');
                    hold off;
                
                    title(strcat('Inter-antenna Phase Offsets (unwrapped). AP:', num2str(ap_i)));
                end
            end
        end

        function [offset_a12, offset_a13, offset_a14] = calculate_antenna_offsets(csi)
            offset_a12 = angle(csi(:, 1)' * csi(:, 2));
            offset_a13 = angle(csi(:, 1)' * csi(:, 3));
            offset_a14 = angle(csi(:, 1)' * csi(:, 4));
        end

        function [csi_packet_trimmed] = apply_tap_filter(csi_packet, threshold, plot_steps)
            % Plot CSI amplitude for this packet
            if plot_steps
                figure
                plot(1:size(csi_packet, 1), abs(csi_packet));
                title('CSI Amplitude: before tap filtering');
            end
            
            % 4. Apply IFFT to the packet
            cir_packet = ifft(csi_packet);
            
            if plot_steps
                figure
                plot(1:size(cir_packet, 1), abs(cir_packet));
                title('CIR: before tap filtering');
            end
            
            % 5. Calculate average channel power
            % How: bring each real CIR component to power
            U_packet = real(cir_packet).^2;
            
            C_packet = zeros(size(U_packet));
            for i = 1:size(C_packet, 1)
                C_packet(i) = sum(U_packet(1:i));
            end
            C_packet = C_packet / sum(U_packet);
            
            if plot_steps
                [idx, ~] = find(C_packet > threshold/100);

                figure
                hold on;
                plot(1:size(C_packet, 1), abs(C_packet));
                plot([idx(2), idx(2)], [0, 1], '--red');
                hold off;

                xlabel('IFFT (Tap) Index');
                ylabel('Cumulative Contribution Rate, %');
                title('Packet Cumulative Contribution Rate');

                xlim([0, 40]);
                ylim([0, 1]);

                legend('', 'Threshold: 99%');
            end
            
            % 6. Use threshold to truncate CIR, using cumulative contribution rate
            tap_count = sum(C_packet <= threshold);
            cir_packet_trimmed = cir_packet(1:tap_count);
            
            if plot_steps
                figure
                plot(1:size(cir_packet_trimmed, 1), abs(cir_packet_trimmed));
                title('CIR: after tap filtering');
            end
            
            % 7. Apply FFT to bring trimmed CIR back to CFR form
            csi_packet_trimmed = fft(cir_packet_trimmed, 234);
            
            if plot_steps
                figure
                plot(1:size(csi_packet_trimmed, 1), abs(csi_packet_trimmed));
                title('CSI Amplitude: after tap filtering');
            end
        end

        function [csi_adjusted] = remove_sto(csi)
            % Input matrix: LxMxN, where L - # of packets, M - # of antennas, N - # of subcarriers
            csi_adjusted = zeros(size(csi));
        
            nfft = size(csi, 3);
        
            % Perform STO estimation separately for each antenna
            for ant_idx = 1:size(csi, 2)
                % 1. Obtain PDP from CIR for each of L packets & calculate Nsto
                Nsto = zeros(size(csi, 1), 1);
                for packet_idx = 1:size(csi, 1)
                    packet = squeeze(csi(packet_idx, ant_idx, :));
                    pdp = abs(ifft(packet, nfft));
                    [~, Nsto_packet] = max(pdp.^2);
                    Nsto(packet_idx, 1) = Nsto_packet;
                end
            
                % 2. Apply STO offset: -2 * pi * k * mean(Nsto) / K, where K is the
                % number of taps (or nfft parameter from FFT)
                sto_offset = permute(exp(-2 * pi * (0:233) * mean(Nsto) / nfft), [2, 1]);
        
                % 3. Create adjusted CSI measurements
                for packet_idx = 1:size(csi, 1)
                    csi_adjusted(packet_idx, ant_idx, :) = permute(csi(packet_idx, ant_idx, :), [3, 2, 1]) .* sto_offset;
                end
            end
        end

        function [csi_adjusted] = remove_sfo(csi)
            % Input matrix: MxN, where M - # of antennas, N - # of subcarriers
        
            % Unwrap phase from CSI matrix (with tolerance of pi, and against
            % subcarrier dimension)
            phase_matrix = unwrap(angle(csi), pi, 2);
        
            % Formulate X values for fitting
            ant_n = size(csi, 1);
            subcarr_n = size(csi, 2);
        
            for ant_idx = 1:ant_n
                for subcarr_idx = 1:subcarr_n
                    start_idx = 1 + subcarr_n * (ant_idx - 1);
                    fit_X(start_idx:start_idx + subcarr_n - 1) = 1:1:subcarr_n;
                end
            end
        
            % Formulate Y values for fitting
            fit_Y = zeros(subcarr_n * ant_n, 1);
            fit_Y_idx = 1;
            for i = 1:size(phase_matrix, 1)
                for j = 1:size(phase_matrix, 2)
                    fit_Y(fit_Y_idx) = phase_matrix(i, j);
                    fit_Y_idx = fit_Y_idx + 1;
                end
            end
        
            % Perform linear fitting
            % Note: linear fit is common across all antennas
            result = polyfit(fit_X, fit_Y, 1);
            tau = result(1);
        
            % Calculate corrected phase using calculated tau
            phase_corrected = phase_matrix - repmat((0:(subcarr_n - 1)) * tau, ant_n, 1);
            
            % Reconstruct the CSI matrix with the adjusted phase
            csi_adjusted = abs(csi) .* exp(1i * phase_corrected);
        end
    end
end