#!/bin/bash
SEL='(simHit_pabs>50)'

python3 simple_plot.py --infile NTuples_fromGENSIM_v0.root --sel $SEL --var simHit_R --imgname simHit_R
python3 simple_plot.py --infile NTuples_fromGENSIM_v0.root --sel $SEL --var simHit_x --imgname simHit_x
python3 simple_plot.py --infile NTuples_fromGENSIM_v0.root --sel $SEL --var simHit_y --imgname simHit_y
python3 simple_plot.py --infile NTuples_fromGENSIM_v0.root --sel $SEL --var simHit_z --imgname simHit_z
python3 simple_plot.py --infile NTuples_fromGENSIM_v0.root --sel $SEL --var simHit_tof --imgname simHit_tof
python3 simple_plot.py --infile NTuples_fromGENSIM_v0.root --sel $SEL --var simHit_x simHit_y --imgname 2D_simHit_xy
