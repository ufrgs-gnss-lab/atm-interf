function N3 = logavg (N1, N2)
    N3 = exp( (log(N1) + log(N2)) ./ 2 );  % logarithmic mean
    %N3 = (N1 + N2) ./ 2;  % arithmetic mean
    % for quantities decaying logarithmically with altitude, 
    % the layer vertical average \int(N)dh/\int(1)dh corresponds 
    % to the value at the layer center of mass (centroid), 
    % which is below the midpoint.
end
