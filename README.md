# Simple direct solution to Perspective-3-Point problem

&copy; 2020 NEC Corporation

This repository is an official MATLAB implementation of the paper "A Simple Direct Solution to the Perspective-Three-Point Problem", BMVC2019 [(pdf)](https://bmvc2019.org/wp-content/uploads/papers/0533-paper.pdf) [\[1\]](#reference).
Other p3p solvers [\[2-6\]](#reference) evaluated in the paper, ported from C++ to MATLAB, are also included.
For a fair comparison, a quartic equation solver by L. Kneip is used for all solvers except for LambdaTwist [\[6\]](#reference).  
For Japanese readers: 日本語によるP3P問題の解説は私の博士論文 [\[7\]](#reference) をご一読ください．

## License

This software is released under the NEC Corporation License.
See [LICENSE](https://github.com/g9nkn/p3p_problem/LICENSE) before using the code. If you use this code, please cite the paper.

```bibtex
@inproceedings{nakano2019simple,
  title={A Simple Direct Solution to the Perspective-Three-Point Problem},
  author={Nakano, Gaku},
  numpages={12},
  booktitle={Proceedings of the British Machine Vision Conference (BMVC)},
  publisher={BMVA Press},
  year={2019}
}
```

For commercial use, please contact Gaku Nakano \<g-nakano@nec.com\>.

## Usage

### Install

Copy `p3p_problem` folder and set a path to the folder by `addpath(genpath('p3p_problem'))`. If you use p3p solvers in `p3p_problem/other_p3p3_solvers`, please follow their licenses.

### Function API

```
[R, t] = p3p_nakano_bmvc2019(m, X, polishing)

INPUTS:
  m - 3x3 matrix of 2D points represented by homogeneous coordinates.
      Each column m(:,i) corresponds to the 3D point X(:,i),
      [u1, u2, u3
       v1, v2, v3
       w1, w2, w3]
      where each column is normalized by sqrt(u^2+v^2+w^2)=1.
  X - 3x3 matrix of 3D points.
      Each column X(:,i) corresponds to the 2D point m(:,i),
      [x1, x2, x3
       y1, y2, y3
       z1, z2, z3].
  polishing - (optional) an integer to set the number of iterations of
              root polishing. If <= 0, the root polishing is not performed.
              (default: 1)
OUTPUS:
  R - 3x3xM rotation matrix (1<= M <= 4).
      R(:,:,i) corresponds to t(:,i). 
  t - 3xM translation vector.
```

Run `demo_p3p_nakano_bmvc2019.m` to understand how it works.

## Reference

1. Gaku Nakano, "A Simple Direct Solution to the Perspective-Three-Point Problem," BMVC2019.  
<https://bmvc2019.org/wp-content/uploads/papers/0533-paper.pdf>

2. X. Gao et al., "Completesolution classification for the erspective-three-point problem," IEEE PAMI, 25(8):930–943, 2003.  
Code: <https://github.com/opencv/opencv/blob/4.1.0/modules/calib3d/src/p3p.cpp#L132>

3. L. Kneip et al., "A Novel Parametrization of the P3P-Problem for a Direct Computation of Absolute Camera Position and Orientation," CVPR2011.  
Code: ~~<https://dl.dropboxusercontent.com/u/23966023/home_page_files/p3p_code_final.zip>~~ (dead link)

4. T. Ke and S. I. Roumeliotis, "An Efficient Algebraic Solution to the Perspective-Three-Point Problem," CVPR2017.  
Code: <https://github.com/opencv/opencv/blob/4.1.0/modules/calib3d/src/ap3p.cpp#L157>

5. A. Banno, "A P3P problem solver representing all parameters as a linear combination," Image Vision and Computing, 70 (2018) 55-62.  
Code: <https://github.com/atsuhikobanno/p3p>

6. M. Persson and K. Nordberg, "Lambda Twist: An Accurate Fast Robust Perspective Three Point (P3P) Solver," ECCV2018.  
Code: <https://github.com/midjji/pnp/blob/master/lambdatwist/lambdatwist.p3p.h>

7. 中野学，Perspective-n-Point問題とその派生問題に対する安定かつ高速な解法に関する研究，博士論文，筑波大学，2021年3月．
<https://jpn.nec.com/rd/people/docs/doctoral_thesis_nakano.pdf>

## Contributors

- Gaku Nakano, Central Research Laboratories, NEC Corporation.  
<g-nakano@nec.com>  
(ENG) <https://www.nec.com/en/global/rd/people/gaku_nakano.html>  
(JPN) <https://jpn.nec.com/rd/people/gaku_nakano.html>