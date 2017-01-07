# Efficient Wasserstein barycenter in MATLAB

This repository includes the propotype codes for computing the Wasserstein barycenter.
The codes are released to reproduce parts of the experiments reported in the following journal paper:

- Jianbo Ye, Panruo Wu, James Z. Wang and Jia Li, "Fast Discrete Distribution Clustering Using Wasserstein Barycenter with Sparse Support", *IEEE Transactions on Signal Processing* ([arXiv:1510.00012](http://arxiv.org/abs/1510.00012) [stat.CO], September 2015), accepted with mandatory minor revisions.

If you need more efficient and scalable implementation in MPI and C (patent pending), please contact the author.

The codes also implements algorithms in the following papers:

- Jean-David Benamou, Guillaume Carlier, Marco Cuturi, Luca Nenna, and Gabriel Peyré, "Iterative Bregman projections for regularized transportation problems." *SIAM Journal on Scientific Computing* 37.2 (2015): A1111-A1138.

- Jianbo Ye and Jia Li, "Scaling Up Discrete Distribution Clustering Using ADMM", *IEEE International Conference on Image Processing*, Paris, France, October 2014.

- Marco Cuturi and Arnaud Doucet. "Fast Computation of Wasserstein Barycenters.", *International Conference on Machine Learning*, Atlanta, USA, June 2013.

- Jia Li and James Z. Wang. "Real-time computerized annotation of pictures." *IEEE Transactions on Pattern Analysis and Machine Intelligence* 30.6 (2008): 985-1002.




## Algorithm Prototypes

Directory `matlab` contains prototype version in Matlab R2015b.
  
- `profile_centroid.m` --- profiling the convergence of centroid updates
- `Wasserstein_Barycenter.m` --- main function
- `profile_kantorovich.m` --- profiling LP solution of transportation problem
- `centroid_sph*.m` --- computing centroid of a single phase


----
Copyright 2017 Jianbo Ye

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
