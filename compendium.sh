topdir=`pwd`

# Clone rate:
cd $topdir/
git clone https://github.com/pcubillos/rate
#git checkout e585f3e

# Clone TEA
#git clone https://github.com/dzesmin/TEA
cd $topdir/
git clone https://github.com/pcubillos/TEA
cd $topdir/TEA
git checkout e585f3e


# Run VULCAN code:
cd $topdir/code
python $topdir/code/vulcan.py

# 'Do' algebra to get the polynomial coefficients:
python $topdir/code/algebra.py

# Run TEA grid for C>=O parameter space:
cd $topdir/run02_tea
python $topdir/code/tea_grid.py

# Fit turn-over pressure with quartic polynomial:
cd $topdir/run02_tea
python $topdir/code/turnover.py


