function [PKS,LOCS] = find_max(data,X,min_peak_width,max_peak_width,visualize)
    if nargin<=4
        visualize = false;
    end
    %FIND_MIN Find the minimum peaks of the data
    [PKS,LOCS] = findpeaks(data,X,'MinPeakWidth',min_peak_width,"MaxPeakWidth",max_peak_width); % Find the minimum
    if visualize
        plot(X,data);
        hold on
        scatter(LOCS,PKS,'r^','filled');
        hold off
    end
end

