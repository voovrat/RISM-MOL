function [INDEX_MATRIX, DISTANCES ]=build_distance_index2(XYZ,tape,level)

XYZ=XYZ(:,1:3);

D =mydist(XYZ');

D = roundDistanceMatrix2(D,tape,level);
[INDEX_MATRIX,DISTANCES] = getUniqueIndeces2(D);
