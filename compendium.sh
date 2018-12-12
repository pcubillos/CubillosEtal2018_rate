topdir=`pwd`

# Make sure you have installed sympy:
# (e.g., pip install sympy==0.7.6)


# Clone rate:
cd $topdir/
git clone https://github.com/pcubillos/rate
cd $topdir/rate
git checkout ee7b5aa
# Clone TEA
#git clone https://github.com/dzesmin/TEA
cd $topdir/
git clone https://github.com/pcubillos/TEA
cd $topdir/TEA
git checkout 4ad1537

# Run VULCAN code (see results in $topdir/plots/ folder):
cd $topdir
python $topdir/code/vulcan.py

# 'Do' algebra to get the polynomial coefficients:
python $topdir/code/algebra.py


# Run TEA grid for C>=O parameter space:
# Note this will take between ~1 to a few hours and create 500M of data.
cd $topdir/run
python $topdir/code/tea_grid.py

# Fit turn-over pressure with quartic polynomial:
cd $topdir/run
python $topdir/code/turnover.py

# Explore parameter space to find most adequate approach for each case:
cd $topdir/run
python $topdir/code/explore.py

# Benchmark the rate code, explore grid of [T, p, [M/H], C/O] for each species:
cd $topdir/run
python $topdir/code/benchmark.py


# Article plots:
cd $topdir
python $topdir/code/fig_domain.py
python $topdir/code/fig_regimes.py
cd $topdir/run
python $topdir/code/fig_benchmark.py
python $topdir/code/fig_hot.py
