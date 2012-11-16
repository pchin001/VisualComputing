function [ avgerr ] = CalculatePixelError( predicted, actual )
% Uses euclidian distance between two pixels as the error.
% Each point is aligned to each column in 'predicted' and 'actual'.

[~,numPoints] = size(predicted);
diff = abs(actual - predicted);
dist = sqrt(diff(1,:).^2 + diff(2,:).^2);
avgerr = sum(dist) ./ numPoints;

end

