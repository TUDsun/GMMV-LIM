GMMV-LIM
========

Introduction
------------

GMMV-LIM is a MATLAB-based package of the **Linear Inversion Method based on the Generalized Multiple Measurement Vectors model**. This package inverts the [Fresnel data](http://www.fresnel.fr/3Ddatabase/) with multiple frequency components. In addition, it also contains the linear sampling method (LSM). Comparison can be easily made to verify the advantage of GMMV-LIM. 

- Basic formulation

	<div align=center><img src="http://latex.codecogs.com/gif.latex?%5Cmin%5C%20%5Ckappa%28J%29%20%5Cquad%20%5Ctext%7Bs.%20t.%7D%5C%20%5Cleft%5C%7C%5CPhi%20%5Ccdot%20J%20-%20Y%5Cright%5C%7C_F%20%5Cleq%20%5Ctilde%7B%5Csigma%7D"/></div>

	where,

	<div align=center><img src="http://latex.codecogs.com/gif.latex?%5Ckappa%28J%29%3D%5C%7CJ%5C%7C_%7B1%2C2%7D%3A%3D%5Csum_%7Bn%3D1%7D%5EN%5Cleft%5C%7CJ_%7Bn%2C%3A%7D%5ET%5Cright%5C%7C_2"/></div>

- See `INSTALL.md` for installation instruction.


Feedback
--------
We would be delighted to hear from you if you find GMMV-LIM useful, or if you have any suggestions, contributions, or bug reports. Please send these to

Shilong Sun (shilongsun@icloud.com)


Citation
--------

- [*A Linear Model for Microwave Imaging of Highly Conductive Scatterers*](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8123524), S. Sun, B. J. Kooij, and A. G. Yarovoy, IEEE Transactions on Microwave Theory and Techniques, 66 (3), 1149 - 1164, 2018

- [*A linear method for shape reconstruction based on the generalized multiple measurement vectors model*](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8292840), S. Sun, B. J. Kooij, A. G. Yarovoy, and T. Jin, IEEE Transactions on Antennas and Propagation, 66 (4), 2016-2025, 2018
