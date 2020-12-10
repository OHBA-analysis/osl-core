function [ hpi ] = osl_headpos( infile, to_plot )
%function hpi = osl_headpos( infile, to_plot )
%
% Function to extract cHPI information from the output of a maxfilter run.
%
% infile - .pos file maxfilter output, normally found next to the sss file.
% to_plot - optionally set to 1 to produce diagnostic image, default is 0
%
% Outputs are normalised relative to the start position for that scan.
%
% hpi.hpi_fit values should be very high ie > .99
%

if nargin < 2
    to_plot = 0;
end

fid = fopen(infile);

header = fgetl(fid);
C = textscan(fid,'%d%f%f%f%f%f%f%f%f%f','MultipleDelimsAsOne',true);

% transform translation from meters to mm
C{5} = C{5} * 1000;
C{6} = C{7} * 1000;
C{7} = C{8} * 1000;

% find the mean
x_mean = mean(C{5});
y_mean = mean(C{6});
z_mean = mean(C{7});
% normalise to start of recording
C{5} = C{5} - C{5}(1);
C{6} = C{6} - C{6}(1);
C{7} = C{7} - C{7}(1);

translation = [C{5}, C{6}, C{7}];

% decompose rotation matrix
C_0 = sqrt( 1 - C{2}.^2 - C{3}.^2 - C{4}.^2 );
r_11 = C_0.^2 + C{2}.^2 - C{3}.^2 - C{4}.^2;
r_21 = 2 * ( C{2}.*C{3} + C_0.*C{4});
r_31 = 2 * ( C{2}.*C{4} - C_0.*C{3});
r_32 = 2 * ( C{3}.*C{4} + C_0.*C{2});
r_33 = C_0.^2 + C{4}.^2 - C{2}.^2 - C{3}.^2;
theta_x = atan2( r_32 , r_33);
theta_y = atan2(-r_31, sqrt( r_32.^2 + r_33.^2 ) );
theta_z = atan2( r_21, r_11);

%normalise to start of recording
theta_x = theta_x - theta_x(1);
theta_y = theta_y - theta_y(1);
theta_z = theta_z - theta_z(1);

rotation = [theta_x, theta_y, theta_z];

% delta rotation

delta_rot = sqrt( C{2}.^2 + ...
                 C{3}.^2 + ...
                 C{4}.^2 );

% rotation drift

drift_rot = sqrt( cumsum(C{2}).^2 + ...
                 cumsum(C{3}).^2 + ...
                 cumsum(C{4}).^2 );

% delta translation

delta_trans = sqrt( ( C{5}(1:end-1)-C{5}(2:end) ).^2 + ...
                   ( C{6}(1:end-1)-C{6}(2:end) ).^2 + ...
                   ( C{7}(1:end-1)-C{7}(2:end) ).^2 );

% translation drift

drift_trans = sqrt( ( cumsum(C{5}(1:end-1)-C{5}(2:end)) ).^2 + ...
                   ( cumsum(C{6}(1:end-1)-C{6}(2:end)) ).^2 + ...
                   ( cumsum(C{7}(1:end-1)-C{7}(2:end)) ).^2 );

% Build output structure
hpi = [];
hpi.translation = translation;
hpi.rotation = rotation;
hpi.delta_trans = delta_trans;
hpi.drift_trans = drift_trans;
hpi.delta_rot = delta_rot;
hpi.drift_rot = drift_rot;
hpi.hpi_fit = C{8};
hpi.hpi_fit_err = C{9};
hpi.time_vect = C{1};

if to_plot==1;
    h = figure('name','Head position','tag','headpos');

    subplot(311); grid on; hold on
    plot(hpi.time_vect, hpi.hpi_fit,'linewidth',2)
    plot(hpi.time_vect, hpi.hpi_fit + hpi.hpi_fit_err,'linewidth',.5);
    plot(hpi.time_vect, hpi.hpi_fit - hpi.hpi_fit_err,'linewidth',.5);
    %plot(hpi.time_vect, hpi.hpi_fit)
    title(['HPI Goodness of Fit File'])
    ylim([.98 1])

    subplot(312); grid on
    plot(hpi.time_vect, hpi.translation)
    title('Translation over time')
    mx = max(max(abs(hpi.translation)));
    ylim([-mx mx])

    subplot(313); grid on
    plot(hpi.time_vect, hpi.rotation)
    title('Rotation over time')
    mx = max(max(abs(hpi.rotation)));
    ylim([-mx mx]) 
end
