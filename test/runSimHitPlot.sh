#!/bin/bash
SEL='(simHit_isHardProcess*(simHit_genLxy>=100))'

python3 simple_plot.py --infile $EOS/NTuples_fromGENSIM_v0.root --sel $SEL --var simHit_R --imgname simHit_R_genLxy_ge_100
python3 simple_plot.py --infile $EOS/NTuples_fromGENSIM_v0.root --sel $SEL --var simHit_x --imgname simHit_x_genLxy_ge_100
python3 simple_plot.py --infile $EOS/NTuples_fromGENSIM_v0.root --sel $SEL --var simHit_y --imgname simHit_y_genLxy_ge_100
python3 simple_plot.py --infile $EOS/NTuples_fromGENSIM_v0.root --sel $SEL --var simHit_z --imgname simHit_z_genLxy_ge_100
python3 simple_plot.py --infile $EOS/NTuples_fromGENSIM_v0.root --sel $SEL --var simHit_tof --imgname simHit_tof_genLxy_ge_100
python3 simple_plot.py --infile $EOS/NTuples_fromGENSIM_v0.root --sel $SEL --var simHit_x simHit_y --imgname 2D_simHit_xy_genLxy_ge_100
