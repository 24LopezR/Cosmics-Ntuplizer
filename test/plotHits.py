import numpy as np
import ROOT as R
from style import setStyle
import math
from tqdm import tqdm
import os
import pickle 

Xs = np.array([-2.74609,-6.53641,-10.6466,-15.7231,-26.8481,-27.0807,-35.3487,-35.5812,-40.3124,-50.6929,-51.6188,-58.1599,-58.5919,-70.6757,-71.1064,-75.5727,-88.1971,-93.9702,-108.578])
Ys = np.array([0.307657,0.738272,1.20632,1.78576,3.06145,3.08825,4.04121,4.06825,4.61531,5.82229,5.92942,6.69269,6.74339,8.16311,8.21325,8.74054,10.2355,10.9233,12.6715])

gen_particle_hits = {}
DR_TRESH = 0.01

def deltaR(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = abs(phi1 - phi2)
    if dphi > math.pi: dphi -= 2*math.pi
    dR = np.sqrt(deta*deta + dphi*dphi)
    return dR

if __name__=='__main__':
    #R.gROOT.ProcessLine('.L ./tdrstyle2D.C')
    R.gROOT.SetBatch(1)
    #R.setTDRStyle()
    setStyle(width=800, height=800, font=42, fontsize=0.04)


    if not os.path.exists('./gentracks.pkl'):
        # Build hits
        f = R.TFile.Open("NTuples_fromGENSIM_v0.root", "READ")
        tree = f.Get("Events")

        for ev in tqdm(tree, total=tree.GetEntries()):
            for i in range(len(ev.gen_pdgId)):
                if ev.gen_status[i] != 1: continue
                gen_eta = ev.gen_eta[i]
                gen_phi = ev.gen_phi[i]
                for j in range(len(ev.simTrack_trackId)):
                    simTrack_eta = ev.simTrack_eta[j]
                    simTrack_phi = ev.simTrack_phi[j]
                    if deltaR(gen_eta, gen_phi, simTrack_eta, simTrack_phi) > DR_TRESH: continue
                    key = str(ev.simTrack_trackId[j])
                    if not key in gen_particle_hits: 
                        gen_particle_hits[key] = {"id": ev.simTrack_trackId[j], "x": [], "y": [], "z": []}
                    for k in range(len(ev.simHit_trackId)):
                        if ev.simHit_trackId[k] != ev.simTrack_trackId[j]: continue
                        else:
                            gen_particle_hits[key]["x"].append(ev.simHit_globalPosition_x[k])
                            gen_particle_hits[key]["y"].append(ev.simHit_globalPosition_y[k])
                            gen_particle_hits[key]["z"].append(ev.simHit_globalPosition_z[k])
                    
        with open('gentracks.pkl', 'wb') as out:
            pickle.dump(gen_particle_hits, out)

    with open('gentracks.pkl', 'rb') as iin:
        gen_particle_hits = pickle.load(iin)
    print(f'Number of muon SimTracks with hits: {len(gen_particle_hits.keys())}')

    graphs_xy = []
    graphs_zx = []
    graphs_zy = []
    for key in gen_particle_hits:
        gen_track = gen_particle_hits[key]
        if len(gen_track["x"]) == 0: continue
        gxy_temp = R.TGraph(len(gen_track["x"]), np.array(gen_track["x"]), np.array(gen_track["y"]))
        gxy_temp.SetTitle(f"Sim Track {key};x (cm);y (cm)")
        graphs_xy.append(gxy_temp)
        gzx_temp = R.TGraph(len(gen_track["x"]), np.array(gen_track["z"]), np.array(gen_track["x"]))
        gzx_temp.SetTitle(f"Sim Track {key};z (cm);x (cm)")
        graphs_zx.append(gzx_temp)
        gzy_temp = R.TGraph(len(gen_track["x"]), np.array(gen_track["z"]), np.array(gen_track["y"]))
        gzy_temp.SetTitle(f"Sim Track {key};z (cm);y (cm)")
        graphs_zy.append(gzy_temp)

    c = R.TCanvas("cxy","cxy")
    c.cd()
    mg = R.TMultiGraph()
    mg.SetTitle(";x (cm);y (cm)")
    for g in graphs_xy: 
        mg.Add(g)
    mg.Draw("AL")
    c.SaveAs("/eos/user/r/rlopezru/www/SignalModelStudies/Stop_testplots/tracks_xy.png")
    c = R.TCanvas("czx","czx")
    c.cd()
    mg = R.TMultiGraph()
    mg.SetTitle(";z (cm);x (cm)")
    for g in graphs_zx: 
        mg.Add(g)
    mg.Draw("AL")
    c.SaveAs("/eos/user/r/rlopezru/www/SignalModelStudies/Stop_testplots/tracks_zx.png")
    c = R.TCanvas("czy","czy")
    c.cd()
    mg = R.TMultiGraph()
    mg.SetTitle(";z (cm);y (cm)")
    for g in graphs_zy: 
        mg.Add(g)
    mg.Draw("AL")
    c.SaveAs("/eos/user/r/rlopezru/www/SignalModelStudies/Stop_testplots/tracks_zy.png")
