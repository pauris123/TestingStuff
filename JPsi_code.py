# Code for J/Psi with uproot 2 and 3
import time
bigBang=time.time()
import uproot as up
import uproot3 as up3
from ROOT import TLorentzVector
import numpy as np
import pandas as pn
import sys 
import awkward as ak

def wall_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "{0:.0f}h:{1:.0f}min:{2:.0f}s".format(h,m,s)
def l_wall_time(seconds):
    m, s = divmod(seconds, 60)
    return "{0:.0f}min:{1:.0f}s".format(m,s)

File_Names = pn.read_csv("FileNames_2017_v2_JPsi_MuonEG.txt",sep='\t',names=["FileNames"])

print("Got all the libraries and filenames")

From = 0#int(sys.argv[1])
To = len(File_Names["FileNames"])#int(sys.argv[2])

cloumn_names = ["index","Electron_pt","nElectron","Electron_eta","Electron_charge","Electron_phi","Electron_cutBased","Electron_r9"]
cloumn_names_extra = ["Electron_pt","Electron_eta","Electron_charge","Electron_phi","Electron_cutBased","Electron_r9","event_i","inv_mass"]

DF_FULL = pn.DataFrame(columns = cloumn_names)

DF_RESULTS = pn.DataFrame(columns = cloumn_names_extra)



for i in range(From,To):
    
    start=time.time()
    
    print("Working on file number {0}/{1}, {2}".format(i+1,To,wall_time(time.time()-bigBang)))
    
    rootfile="root://cmsxrootd.fnal.gov///"+File_Names["FileNames"][i]
    tree = up.open(rootfile)["Events"]
    #up3_tree = up3.open(rootfile)["Events"]
    print("Tree has been read in! {0}".format(l_wall_time(time.time()-start)))
    branches = ["Electron_pt","nElectron","Electron_eta","Electron_charge","Electron_phi","Electron_cutBased","Electron_r9"]
    
    DataTable = tree.pandas.df(branches, entrystop=-1)
    print("I have the branches! {0}".format(l_wall_time(time.time()-start)))
    df = DataTable.query("nElectron==2")
    df = df.reset_index()
    print("I have {0} diE events! {1}".format(len(df["nElectron"]),l_wall_time(time.time()-start)))
    DF_FULL = DF_FULL.append(df)
    
    el_InvMass = []
    el_event_i = []
    
    el_Pt = []
    el_Eta = []
    el_Charge = []
    el_Phi = []
    el_CutBased = []
    el_R9 = []
    
    
    for y in range(len(df["nElectron"])):
        if (df["Electron_charge"][y][0]+df["Electron_charge"][y][1] == 0):
            if ((df["Electron_pt"][y][0] > 7 and abs(df["Electron_eta"][y][0]) < 1.479 and df["Electron_cutBased"][y][0] in [2,3,4]) or \
            (df["Electron_pt"][y][1] > 7 and abs(df["Electron_eta"][y][1]) < 1.479 and df["Electron_cutBased"][y][1] in [2,3,4])):
            
            
                el1 = TLorentzVector()
                el2 = TLorentzVector()

                el1.SetPtEtaPhiM(df["Electron_pt"][y][0],
                                 df["Electron_eta"][y][0],
                                 df["Electron_phi"][y][0],0.000511)
                el2.SetPtEtaPhiM(df["Electron_pt"][y][1],
                                 df["Electron_eta"][y][1],
                                 df["Electron_phi"][y][1],0.000511)
                if ((el1+el2).M() < 20):
                    el_InvMass.append((el1+el2).M())
                    el_event_i.append(y)
    
                    el_Pt.append(df["Electron_pt"][y])
                    el_Eta.append(df["Electron_eta"][y])
                    el_Charge.append(df["Electron_charge"][y])
                    el_Phi.append(df["Electron_phi"][y])
                    el_CutBased.append(df["Electron_cutBased"][y])
                    el_R9.append(df["Electron_r9"][y])
    
    results = pn.DataFrame(list(zip(el_Pt,el_Eta,el_Charge,el_Phi,el_CutBased,el_R9,el_event_i,el_InvMass)),columns = cloumn_names_extra)
    DF_RESULTS = DF_RESULTS.append(results)
    
    print("I have {0} raw data events! {1}".format(len(el_InvMass),l_wall_time(time.time()-start)))
    print(" Each event took {0:.2f} ms".format(1000*(time.time()-start)/len(el_InvMass)))

print("ROOT file writing has started! {0}".format(wall_time(time.time()-bigBang)))     
# Arrays for FULL data    
ak0_array_pt = ak.to_awkward0(ak.Array(DF_FULL["Electron_pt"]))
ak0_array_eta = ak.to_awkward0(ak.Array(DF_FULL["Electron_eta"]))
ak0_array_charge = ak.to_awkward0(ak.Array(DF_FULL["Electron_charge"]))
ak0_array_phi = ak.to_awkward0(ak.Array(DF_FULL["Electron_phi"]))
ak0_array_cut = ak.to_awkward0(ak.Array(DF_FULL["Electron_cutBased"]))
ak0_array_r9 = ak.to_awkward0(ak.Array(DF_FULL["Electron_r9"]))

#Arrays for RESULT data
ak0_array_pt_R = ak.to_awkward0(ak.Array(DF_RESULTS["Electron_pt"]))
ak0_array_eta_R = ak.to_awkward0(ak.Array(DF_RESULTS["Electron_eta"]))
ak0_array_charge_R = ak.to_awkward0(ak.Array(DF_RESULTS["Electron_charge"]))
ak0_array_phi_R = ak.to_awkward0(ak.Array(DF_RESULTS["Electron_phi"]))
ak0_array_cut_R = ak.to_awkward0(ak.Array(DF_RESULTS["Electron_cutBased"]))
ak0_array_r9_R = ak.to_awkward0(ak.Array(DF_RESULTS["Electron_r9"]))  

array_i_R = np.array(DF_RESULTS["event_i"])
array_InvM_R = np.array(DF_RESULTS["inv_mass"],dtype = np.float32)

#ROOT file for FULL data
ROOT_File_FULL = up3.recreate("2017_JPsi_MuonEG_test_FULL.root")
ROOT_File_FULL["tree1"] = up3.newtree({"Electron_Pt":up3.newbranch(np.dtype(np.float32), size="n"),
                                        "Electron_Eta":up3.newbranch(np.dtype(np.float32), size="n"),
                                        "Electron_Charge":up3.newbranch(np.dtype(np.int32), size="n"),
                                        "Electron_Phi":up3.newbranch(np.dtype(np.float32), size="n"),
                                        "Electron_CutBased":up3.newbranch(np.dtype(np.int32), size="n"),
                                        "Electron_R9":up3.newbranch(np.dtype(np.float32), size="n")})
        
ROOT_File_FULL["tree1"].extend({"Electron_Pt": ak0_array_pt, "n": ak0_array_pt.counts,
                                "Electron_Eta": ak0_array_eta, "n": ak0_array_eta.counts,
                                "Electron_Charge": ak0_array_charge, "n": ak0_array_charge.counts,
                                "Electron_Phi": ak0_array_phi, "n": ak0_array_phi.counts,
                                "Electron_CutBased": ak0_array_cut, "n": ak0_array_cut.counts,
                                "Electron_R9": ak0_array_r9, "n": ak0_array_r9.counts})
ROOT_File_FULL.close 

# ROOT file for RESULT data
ROOT_File_RESULTS = up3.recreate("2017_JPsi_MuonEG_test_RESULTS.root")
ROOT_File_RESULTS["tree1"] = up3.newtree({"Electron_Pt":up3.newbranch(np.dtype(np.float32), size="n"),
                                        "Electron_Eta":up3.newbranch(np.dtype(np.float32), size="n"),
                                        "Electron_Charge":up3.newbranch(np.dtype(np.int32), size="n"),
                                        "Electron_Phi":up3.newbranch(np.dtype(np.float32), size="n"),
                                        "Electron_CutBased":up3.newbranch(np.dtype(np.int32), size="n"),
                                        "Electron_R9":up3.newbranch(np.dtype(np.float32), size="n")})

ROOT_File_RESULTS["tree2"] = up3.newtree({"Event_i": np.int32,"Inv_Mass": np.float32})
        
ROOT_File_RESULTS["tree1"].extend({"Electron_Pt": ak0_array_pt_R, "n": ak0_array_pt_R.counts,
                                "Electron_Eta": ak0_array_eta_R, "n": ak0_array_eta_R.counts,
                                "Electron_Charge": ak0_array_charge_R, "n": ak0_array_charge_R.counts,
                                "Electron_Phi": ak0_array_phi_R, "n": ak0_array_phi_R.counts,
                                "Electron_CutBased": ak0_array_cut_R, "n": ak0_array_cut_R.counts,
                                "Electron_R9": ak0_array_r9_R, "n": ak0_array_r9_R.counts})

ROOT_File_RESULTS["tree2"].extend({"Event_i": array_i_R,"Inv_Mass": array_InvM_R})


ROOT_File_RESULTS.close 

print("All done! {0}".format(wall_time(time.time()-bigBang)))  
