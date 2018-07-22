close all
clear

data_dir='hand_data_ground_truth/';                    % path to hand data ground truth label directory
data_name='sample_hand_data';                          % ground truth file name prefix
dim=26*3;                                              % dimension of collected hand data, given (Fx, Fy, Fz) on each finger phalanx
nlabel=6;

% fit gaussian to compute annealing likelihood
gaussian_fit(data_dir, data_name, dim, nlabel);