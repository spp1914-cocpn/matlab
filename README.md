# MATLAB

All not directly OMNeT++ co-simulation related MATLAB code (libraries, generic NCS and stuff) resides here.

For MATLAB-OMNeT++ co-simulation API, see ../libncs_matlab

Reference implementations of the following control and estimation algorithms can be found here:
* F. Rosenthal, B. Noack, and U.D. Hanebeck, [*State Estimation in Networked Control Systems with Delayed and Lossy Acknowledgments,*](https://doi.org/10.1007/978-3-319-90509-9_2) in: Lee S., Ko H., Oh S. (eds) Multisensor Fusion and Integration in the Wake of Big Data, Deep Learning and Cyber Physical System, Lecture Notes in Electrical Engineering, Volume 501, Springer, Cham, 2018.
* F. Rosenthal, B. Noack, and U.D. Hanebeck, [*Sequence-Based Receding Horizon Control Over Networks with Delays and Data Losses,*](https://doi.org/10.23919/ACC.2019.8815149) Proceedings of the 2019 American Control Conference (ACC), Philadelphia, PA, USA, 2019. 
* F. Rosenthal and U. D. Hanebeck, [*Sequence-Based Stochastic Receding Horizon Control Using IMM Filtering and Value Function Approximation,*](https://doi.org/10.1109/CDC40024.2019.9029717) Proceedings of the 58th IEEE Conference on Decision and Control (CDC), Nice, France, 2019.
* F. Rosenthal and U. D. Hanebeck, [*Stability Analysis of Polytopic Markov Jump Linear Systems with Applications to Sequence-Based Control over Networks,*](https://doi.org/10.1016/j.ifacol.2020.12.1030) in: IFAC-PapersOnLine, Volume 53, Issue 2, 2020.  


## Externals
The following external libraries/functions are used and included in the folder 'external'.
* [Nonlinear Estimation Toolbox (GPLv3)](https://nonlinearestimation.bitbucket.io/) by Jannik Steinbring, only the required subset is included
* [YALMIP](https://yalmip.github.io/) by Johan Löfberg
* [DiscreteSample (FreeBSD)](https://de.mathworks.com/matlabcentral/fileexchange/21912-sampling-from-a-discrete-distribution) by Dahua Lin
* [mtimesx (FreeBSD)](https://de.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support) by James Tursa
* [Armadillo (Apache)](http://arma.sourceforge.net/) by Conrad Sanderson and Ryan Curtin
* [SDPT3 (GPLv2)](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/) by Kim-Chuan Toh, Michael J. Todd, and Reha H. Tutuncu, only the required subset is included

<br>
<br>

For additional information regarding YALMIP, Armadillo, and SDPT3, please refer to the corresponding papers:
* Johan Löfberg, [YALMIP: a toolbox for modeling and optimization in MATLAB](https://doi.org/10.1109/CACSD.2004.1393890), Proceedings of the 2004 IEEE International Symposium on Computer Aided Control Systems Design, Taipei, Taiwan, 2004.
* Conrad Sanderson and Ryan Curtin, [Armadillo: a template-based C++ library for linear algebra](http://arma.sourceforge.net/armadillo_joss_2016.pdf), Journal of Open Source Software, Vol. 1, pp. 26, 2016.
* Kim-Chuan Toh, Michael J. Todd, and Reha H. Tutuncu, [SDPT3 — A Matlab software package for semidefinite programming](https://doi.org/10.1080/10556789908805762), Optimization Methods and Software, 11 (1999), pp. 545–581.


