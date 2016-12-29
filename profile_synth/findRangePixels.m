function [R_idx,R_theta] = findRangePixels(zz_dist,R)
% Given an array of distances from a specified point, finds pixels that are
% a distance R from the point. The script also outputs the angle from the 
% point. 
%
% Inputs: 
% - zz_dist: Array of distances from a point. It is assumed that the point
%            is located within the array as a distance of 0. 
% - R: Range vector. Currently scripted to work only with 1 value.
%
% Outputs: R_idx: Logical array that identifies pixels that are a distance
%                 R from the point. 
%
% TJ Young
% 07 December 2016

% Pre-allocate matrices
R_idx = false(size(zz_dist));

% Run through each column to find first pixel nearest to specified depth
[R_rmag,R_ridx] = min(abs(zz_dist-R)); 

% Checks: shallow ranges
R_ridx(zz_dist(1,:) > R) = NaN; 

% Checks: deep ranges
R_dr = zz_dist(end,:) - zz_dist(end-1,:);
R_ridx(R-zz_dist(end,:) > R_dr) = NaN; 

for jj = 1:size(zz_dist,2) 
    % Assign indices to master index
    if ~isnan(R_ridx(jj))
        R_idx(R_ridx(jj),jj) = 1;
    end
end

% % Run through each column to find first pixel nearest to specified depth
% for jj = 1:size(zz_dist,2) 
%     
%     % Magnitude and index of cell nearest to R
%     [R_rmag(jj),R_ridx(jj)] = min(abs(zz_dist(:,jj)-R)); 
%     
%     % Checks: shallow ranges
%     if zz_dist(1,jj) > R
%         R_ridx(jj) = NaN; % All pixels in column too far from R
%     end
%     % Checks: deep ranges
%     if zz_dist(end,jj) < R
%         R_dr = zz_dist(end,jj) - zz_dist(end-1,jj);
%         if R - zz_dist(end,jj) > R_dr
%             R_ridx(jj) = NaN; % R is deeper than all pixels in column
%         end
%     end
%     
%     % Assign indices to master index
%     if ~isnan(R_ridx(jj))
%         R_idx(R_ridx(jj),jj) = 1;
%     end
%     
% end