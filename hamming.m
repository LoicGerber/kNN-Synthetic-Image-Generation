function hammingDist = hamming(x,y)

index = ~(isnan(x) | isnan(y));
validX = x(index);
validY = y(index);

binaryX = (validX>0);
binaryY = (validY>0);
hammingDist = mean(binaryX(:) ~= binaryY(:));

end
