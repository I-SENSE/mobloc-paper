classdef music
    methods (Static)
        
        %% 1D MUSIC implementation
        function [theta, Pmusic] = run_music(X, d, freq, ant_n)
            % Input (X) must have a shape: (M, N), where M -- # of antennas, N -- # of samples
        
            K = 1; % The number of signal sources
            lambda = physconst('LightSpeed')/freq; % Wavelength (meters)
            
            % Run MUSIC
            Rx = cov(X'); %Data covarivance matrix 
            [eigenVec, ~] = eig(Rx); %Find the eigenvalues and eigenvectors of Rx 
        
            Vn = eigenVec(:,1:ant_n-K);
            
            theta = -90:1:90;     

            % Perform search through all available values of DoA (theta)
            for i=1:length(theta)
                % Compose steering vector [4x1]
                SS = zeros(ant_n,1); 
                SS = exp(-1j*2*pi*d*(0:ant_n-1)'*sind(theta(i))/lambda);

                % Obtain pseudospectrum value for this theta value
                Pmusic(i) = 1/(SS'*(Vn*Vn')*SS); 
            end
            Pmusic = real(10*log10(Pmusic)); % Spatial Spectrum function (converts to dB)
        end
    end
end