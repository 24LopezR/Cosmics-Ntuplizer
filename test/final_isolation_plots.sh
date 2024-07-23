#!/bin/bash
reference='/eos/user/r/rlopezru/MuonReco/TKAlIssue/NTuples_FinalTest_RelValZMM_14_nominal_v1.root'
target='/eos/user/r/rlopezru/MuonReco/TKAlIssue/NTuples_FinalTest_RelValZMM_14_target_v2.root'

mkdir -p /eos/user/r/rlopezru/www/MuonReco/TkAlIssue/finalTest_HS/
python3 simple_plot.py --infile $reference --var mu_eta mu_iso03_hadEt/mu_pt --imgname finalTest_HS/ZMM_eta_iso03hadEt_nominal
python3 simple_plot.py --infile $reference --var mu_eta mu_iso03_emEt/mu_pt --imgname finalTest_HS/ZMM_eta_iso03emEt_nominal
python3 simple_plot.py --infile $reference --var mu_phi mu_iso03_emEt/mu_pt --imgname finalTest_HS/ZMM_phi_iso03emEt_nominal
python3 simple_plot.py --infile $reference --var mu_phi mu_iso03_hadEt/mu_pt --imgname finalTest_HS/ZMM_phi_iso03hadEt_nominal
python3 simple_plot.py --infile $reference --var mu_phi mu_eta mu_iso03_emEt/mu_pt --imgname finalTest_HS/ZMM_phi_eta_iso03emEt_nominal
python3 simple_plot.py --infile $reference --var mu_phi mu_eta mu_iso03_hadEt/mu_pt --imgname finalTest_HS/ZMM_phi_eta_iso03hadEt_nominal
python3 simple_plot.py --infile $target --var mu_phi mu_iso03_hadEt/mu_pt --imgname finalTest_HS/ZMM_phi_iso03hadEt_target
python3 simple_plot.py --infile $target --var mu_phi mu_iso03_emEt/mu_pt --imgname finalTest_HS/ZMM_phi_iso03emEt_target
python3 simple_plot.py --infile $target --var mu_eta mu_iso03_emEt/mu_pt --imgname finalTest_HS/ZMM_eta_iso03emEt_target
python3 simple_plot.py --infile $target --var mu_eta mu_iso03_hadEt/mu_pt --imgname finalTest_HS/ZMM_eta_iso03hadEt_target
python3 simple_plot.py --infile $target --var mu_phi mu_eta mu_iso03_emEt/mu_pt --imgname finalTest_HS/ZMM_phi_eta_iso03emEt_target
python3 simple_plot.py --infile $target --var mu_phi mu_eta mu_iso03_hadEt/mu_pt --imgname finalTest_HS/ZMM_phi_eta_iso03hadEt_target
