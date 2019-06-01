# willard-chandler
Python-based tool to calculate instantaneous interfaces and concentration/orientation profiles from molecular simulation trajectories in slab geometry.

The method is described by A.P. Willard and D. Chandler in [_Instantaneous Liquid Interfaces_](https://pubs.acs.org/doi/10.1021/jp909219k). The algorithm is an adaptation of parts of the [`pytim` code](https://github.com/Marcello-Sega/pytim) by M. Sega.

# usage

python ./willardchandler.py -s <top> -t <traj> -o <output> -b <begin> -e <every> -m <molecules> -l <layers>

for example

python ./willardchandler.py -s npt.gro -t npt.xtc -o out.p -b 1000 -e 10 -m 'H1 O' 'NA' 'S1 C2 N3' -l -0.1 0.1 0.3 0.5
