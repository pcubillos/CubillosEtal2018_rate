topdir=`pwd`

# Make sure you have installed sympy:
# (pip install sympy)


# Clone rate:
cd $topdir/
git clone https://github.com/pcubillos/rate
cd $topdir/rate
git checkout 6d364b3

# Clone TEA
#git clone https://github.com/dzesmin/TEA
cd $topdir/
git clone https://github.com/pcubillos/TEA
cd $topdir/TEA
git checkout e585f3e


# Run VULCAN code (see results in $topdir/plots/ folder):
cd $topdir/code
python $topdir/code/vulcan.py


# 'Do' algebra to get the polynomial coefficients:
cd $topdir/code
python $topdir/code/algebra.py


# Run TEA grid for C>=O parameter space:
cd $topdir/run02_tea
python $topdir/code/tea_grid.py

# Fit turn-over pressure with quartic polynomial:
cd $topdir/run02_tea
python $topdir/code/turnover.py


