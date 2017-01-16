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
Copyright (C) 2017 The Pennsylvania State University, USA - All Rights Reserved 

Installation and use of this software for academic, non-profit, or government-sponsored research purposes is hereby granted. Use of the software under this license is restricted to non-commercial purposes. COMMERCIAL USE OF THE SOFTWARE REQUIRES A SEPARATELY EXECUTED WRITTEN LICENSE AGREEMENT with the copyright holder.

Any use of the software, and any modifications, improvements, or derivatives to the software the user(s) may create (collectively, “Improvements”) must be solely for internal, non-commercial purposes. A user of the software shall not distribute or transfer the software or Improvements to any person or third parties without prior written permission from copyright holder.

Any publication of results obtained with the software must acknowledge its use by an appropriate citation to the article titled "Fast Discrete Distribution Clustering Using Wasserstein Barycenter with Sparse Support" published in the IEEE Transactions on Signal Processing in 2017.

ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. THE COPYRIGHT HOLDER MAKES NO REPRESENTATION OR WARRANTY THAT THE SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Author: Jianbo Ye, College of Information Sciences and Technology, yelpoo@gmail.com
