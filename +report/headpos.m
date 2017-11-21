function h = headpos(pos_file)
    hpi = osl_headpos(pos_file);

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

