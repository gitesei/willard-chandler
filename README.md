# willard-chandler
Python-based tool to calculate instantaneous interfaces and concentration/orientation profiles from molecular simulation trajectories in slab geometry.

The method is described by A.P. Willard and D. Chandler in [_Instantaneous Liquid Interfaces_](https://pubs.acs.org/doi/10.1021/jp909219k). The algorithm is an adaptation of parts of the [`pytim` code](https://github.com/Marcello-Sega/pytim) by M. Sega.

# usage

To calculate concentration profiles, orientational distributions of linear ions and free OH groups w.r.t. the rough instantaneous interface:

python ./willardchandler.py -s TOPOLOGY -t TRAJECTORY -o OUTPUT PICKLE FILE -b FIRST FRAME -e FRAME INTERVAL -m MOLECULES -l LAYERS

for example (NaSCN aqueous solution)

python ./willardchandler.py -s npt.gro -t npt.xtc -o out.p -b 1000 -e 10 -m 'H1 O' 'NA' 'S1 C2 N3' -l -0.1 0.1 0.3 0.5

To calculate concentration profiles, orientational distributions of linear ions and free OH groups w.r.t. the planar Gibbs dividing surface:

python ./gds.py -s TOPOLOGY -t TRAJECTORY -o OUTPUT PICKLE FILE -b FIRST FRAME -e FRAME INTERVAL -m MOLECULES -l LAYERS

for example (NaSCN aqueous solution)

python ./willardchandler.py -s npt.gro -t npt.xtc -o out.p -b 1000 -e 10 -m 'H1 O' 'NA' 'S1 C2 N3' -l -0.1 0.1 0.3 0.5

