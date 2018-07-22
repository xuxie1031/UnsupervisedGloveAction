%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main script of clustering pipeline %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dir='hand_data/';									% the path to hand data directory
data_name='sample_hand_data';                           % hand data file name
dim=26*3;												% dimension of collected hand data, given (Fx, Fy, Fz) on each finger phalanx

% dump hand feature
% feature_dump_vec(data_dir, [data_name,'.csv'], dim);

% hierarchical clustering
minC = 9;
maxC = 11;												% specify the minimum and maximum cluster numbers
% hClustering(minC, maxC);

% prepare ACA binaries input
file_names = {data_name};								% data file names need to generate ACA binaries
ACAUtils(file_names, minC, maxC);

