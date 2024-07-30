#!/bin/bash
reference='/eos/user/r/rlopezru/MuonReco/TKAlIssue/NTuples_FinalTest_RelValZMM_14_nominal_v1.root'
target='/eos/user/r/rlopezru/MuonReco/TKAlIssue/NTuples_FinalTest_RelValZMM_14_target_v5_mixed.root'

#mkdir -p /eos/user/r/rlopezru/www/MuonReco/TkAlIssue/finalTest/
python3 simple_plot.py --infile $reference --var mu_eta mu_iso03_hadEt/mu_pt --imgname finalTest/ZMM_eta_iso03hadEt_nominal
python3 simple_plot.py --infile $reference --var mu_eta mu_iso03_emEt/mu_pt --imgname finalTest/ZMM_eta_iso03emEt_nominal
python3 simple_plot.py --infile $reference --var mu_eta mu_iso03_sumPt/mu_pt --imgname finalTest/ZMM_eta_iso03sumPt_nominal
python3 simple_plot.py --infile $reference --var mu_phi mu_iso03_emEt/mu_pt --imgname finalTest/ZMM_phi_iso03emEt_nominal
python3 simple_plot.py --infile $reference --var mu_phi mu_iso03_hadEt/mu_pt --imgname finalTest/ZMM_phi_iso03hadEt_nominal
python3 simple_plot.py --infile $reference --var mu_phi mu_iso03_sumPt/mu_pt --imgname finalTest/ZMM_phi_iso03sumPt_nominal
python3 simple_plot.py --infile $reference --var mu_phi mu_eta mu_iso03_emEt/mu_pt --imgname finalTest/ZMM_phi_eta_iso03emEt_nominal
python3 simple_plot.py --infile $reference --var mu_phi mu_eta mu_iso03_hadEt/mu_pt --imgname finalTest/ZMM_phi_eta_iso03hadEt_nominal
python3 simple_plot.py --infile $reference --var mu_phi mu_eta mu_iso03_sumPt/mu_pt --imgname finalTest/ZMM_phi_eta_iso03sumPt_nominal
python3 simple_plot.py --infile $target --var mu_phi mu_iso03_hadEt/mu_pt --imgname finalTest/ZMM_phi_iso03hadEt_target_mixed_v5
python3 simple_plot.py --infile $target --var mu_phi mu_iso03_emEt/mu_pt --imgname finalTest/ZMM_phi_iso03emEt_target_mixed_v5
python3 simple_plot.py --infile $target --var mu_phi mu_iso03_sumPt/mu_pt --imgname finalTest/ZMM_phi_iso03sumPt_target_mixed_v5
python3 simple_plot.py --infile $target --var mu_eta mu_iso03_emEt/mu_pt --imgname finalTest/ZMM_eta_iso03emEt_target_mixed_v5
python3 simple_plot.py --infile $target --var mu_eta mu_iso03_hadEt/mu_pt --imgname finalTest/ZMM_eta_iso03hadEt_target_mixed_v5
python3 simple_plot.py --infile $target --var mu_eta mu_iso03_sumPt/mu_pt --imgname finalTest/ZMM_eta_iso03sumPt_target_mixed_v5
python3 simple_plot.py --infile $target --var mu_phi mu_eta mu_iso03_emEt/mu_pt --imgname finalTest/ZMM_phi_eta_iso03emEt_target_mixed_v5
python3 simple_plot.py --infile $target --var mu_phi mu_eta mu_iso03_hadEt/mu_pt --imgname finalTest/ZMM_phi_eta_iso03hadEt_target_mixed_v5
python3 simple_plot.py --infile $target --var mu_phi mu_eta mu_iso03_sumPt/mu_pt --imgname finalTest/ZMM_phi_eta_iso03sumPt_target_mixed_v5
