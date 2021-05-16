Opt-IGFEM-2D ========

This project optimizes parallel and serpentine channels in 2D.

Features
--------

- Self-contained: does not need Abaqus to generate input file 
- Efficient (typically an hour to several hours for optimizing one design)
- Independent optimizations (limited only by the resources available)

Installation
------------

The relevant mex functions for the project IGFEM-Curves-2D first need to be
install (see the README file in the IGFEM-Curves-2D directory). 
Then the mex function mx_assemble_pseudo_adjoint_forces.mex* in the directory
mx_sensitivity need to be installed. The installation process is the same as
that of mx_assemble_sparse.mex* (see the README file in the IGFEM-Curves-2D directory)

Directories/subdirectories
--------------------------
Opt-IGFEM-2D
- M_opt_postprocessing
- M_optimization
- M_tests
- mx_sensitivity
IGFEM-Curves-2D
- Armadillo_packages
- ChannelFiles
- M_channels
- M_error_analysis
- M_FEM
- M_geom_toolbox
- M_postprocessing
- M_preFEM
- M_tests
- mx_FEM
NURBS
- nurbs-toolbox
SISL
- SISL-masters

Support
-------
If you are having issues or would like to develop these codes further, please
contact the author via email.

Walkthrough
-----------
A few walkthrough examples are provided in the file User_manual.pdf in the
directory Manuals

Copyright
---------
2013-2017 University of Illinois, Urbana-Champaign. All rights reserved.

Author
------
Marcus Hwai Yik Tan (marcus8tan@gmail.com)
PhD, Theoretical and Applied Mechanics