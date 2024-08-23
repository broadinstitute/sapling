# Building locally

Install Conan 1.59.0 (we do not yet support using Conan 2.0):
```
pip3 install 'conan==1.59.0'   # See https://docs.conan.io/en/latest/installation.html for details
```

Make some config adjustment to your Conan profile (see https://docs.conan.io/en/latest/getting_started.html)
```
conan profile new default --detect
conan profile update settings.compiler.libcxx=libstdc++11 default
```

Make sure git submodules are checked out in the `third-party` directory:
```
git submodule update --init
```

Set up build directory and install dependencies
```
mkdir -p build/debug && cd build/debug
conan install ../..
```

Then compile as usual:
```
cmake ../.. -DCMAKE_BUILD_TYPE=Debug  # Or Release
make -j 6
```

This results in the executable `sapling`.


# Demo

After succesful compilation, you should be able to run the command
```
sapling --version
```
and obtain a version string like the following:
``
Sapling Version 0.0.1 (build 1, commit 467948f)
```

If that works, you can run a short simulation to verify that everything is working.
```
sapling --seed 1 -n 10 --exp-pop-n0 6.0 --exp-pop-g 10 --t0 2024-07-31 --min-tip-t 2024-01-01 --max-tip-t 2024-07-31 --hky-kappa 5.0 --hky-pi-A 0.3 --hky-pi-C 0.18 --hky-pi-G 0.20 --hky-pi-T 0.32 --mu 0.1 -L 10 --out-info demo_info.json --out-newick demo_tree.nwk --out-fasta demo.fasta --out-maple demo.maple
```
That simulates an exponentially growing population from `min_tip_t = 2024-01-01` to `max_tip_t = 2024-07-31`, with effective population size times generation time in years given by `N(t)*rho = 6.0 exp[10.0 * (t - max_tip_t)]`.  It samples `n = 10` members of that population uniformly (i.e., their times are distributed proportional to `N(t)`) and runs a coalescent simulation to sample an ancestry for them.  It then simulates a random sequence of length `L = 10` at the root of the tree, with state probabilities of `pi_A = 0.3`, `pi_C = 0.18`, `pi_G = 0.20` and `pi_T = 0.32`.  It then evolves that sequence randomly over the tree using an HKY evolution model with transition-transversion ratio `kappa = 5.0` and mutation rate `mu = 0.1` mutations per site per year.  That produces a detailed tree with concrete sequences at the tips.  A summary of the properties of that tree is written into `demo_info.json`, while the tree itself is written in Newick format to `demo_tree.nwk`.  The tip sequences are written in FASTA format to `demo.fasta` and in MAPLE format to `demo.maple`.  The tip names are in the form `TIP_N|YYYY-MM-DD`, where N is a tip ID starting from 1, and YYYY-MM-DD is the tip date.

That should produce the following output in the console:
```
Seed: 1
Wrote Info to demo_info.json
Wrote Newick to demo_tree.nwk
Wrote FASTA to demo.fasta
Wrote Maple to demo.maple
Total number of mutations: 1
Total branch length: 2.03464 years
```
The contents of the demo files produced should match those in the `demo-data` folder.
