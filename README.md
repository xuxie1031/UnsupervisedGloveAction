Unsupervised Learning of Hierarchical Models for Hand-Object Interactions (ICRA 2018)
======
Source Code of **Clustering** pipeline / **ACA** optimization / **Inference** labeling

## Usage

$ git clone https://github.com/xiaozhuchacha/AtomicAction.git


## Clustering

* use 'main.m' as the entry in 'Clustering' folder to start clustering pipeline
* sample hand data is placed in 'Clustering/hand_data'


## ACA

* use 'ACA.cpp' as the entry in 'ACA' folder to start clustering optimization
* to compile: g++ ACA.cpp -o ACA
* execute: ./ACA [NC] [data_name], eg. ./ACA 9 sample_hand_data


## Inference

* use 'AnnealGibbs.py' as the entry in 'Inference' folder to start Gibbs annealing
* parser is placed in 'Inference/induced_grammar' to calculate prior
* gaussians are placed in 'Inference/Gaussians' to calculate likelihood
* use 'main.m' as the entry in 'Inference/Gaussians' to compute gauusian parameters w.r.t ground truth labeled data
* labeled motion sequence as input can be found either in 'Clustering/hcluster' as hierarchical clustering result or 'ACA/ACAbin' as ACA optimization result
* execute: python AnnealGibbs.py [--nlabel] [--data-name], eg. python AnnealGibbs.py --nlabel 6 --data-name sample_hand_data


## Contact

* xiexu@ucla.edu


## Paper Link

* see gloveaction2018icra.pdf