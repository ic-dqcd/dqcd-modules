import os

from analysis_tools.utils import import_root, randomize
from Base.Modules.baseModules import JetLepMetSyst

ROOT = import_root()


class DQCDULBDTNoJetProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        scenario = kwargs.pop("scenario", "A")
        default_name = "bdt_scenario%s" % scenario
        use_sigmoid = "false"
        if scenario == "A":
            default_model_path = os.path.expandvars(
                # "$CMSSW_BASE/src/DQCD/Modules/data/model_ul_saved.model")
                #"$CMSSW_BASE/src/DQCD/Modules/data/model_ul_saved_scenarioA_no_jet.model")
                "$CMSSW_BASE/src/DQCD/Modules/data/XGB_NOJETS_45_scenarioA_final.model")
                #"$CMSSW_BASE/src/DQCD/Modules/data/model_ul_saved_scenarioA_new.model")
        elif scenario == "B1":
            default_model_path = os.path.expandvars(
                # "$CMSSW_BASE/src/DQCD/Modules/data/model_ul_saved_scenarioB1.model")
                "$CMSSW_BASE/src/DQCD/Modules/data/XGB_NOJETS_46_scenarioB1_final.model")
        elif scenario == "B2":
            raise ValueError("B2 no jet not implemented")
            default_model_path = os.path.expandvars(
                "$CMSSW_BASE/src/DQCD/Modules/data/model_ul_saved_scenarioB2_new.model")
        elif scenario == "C":
            raise ValueError("C no jet not implemented")
            default_model_path = os.path.expandvars(
                "$CMSSW_BASE/src/DQCD/Modules/data/model_ul_saved_scenarioC_new.model")
        elif scenario == "vector":
            default_model_path = os.path.expandvars(
                "$CMSSW_BASE/src/DQCD/Modules/data/XGB_NOJETS_102_vector_final.model")
            default_name = "bdt_vector"
            # use_sigmoid = "true"
        else:
            raise ValueError("Only BDTs for scenarios A, B1, B2, and C are already implemented")

        self.model_path = kwargs.pop("model_path", default_model_path)
        self.model = self.model_path.replace("/", "_").replace(".", "_")
        self.bdt_name = kwargs.pop("bdt_name", default_name)
        # self.model_m = kwargs.pop("model_m", 2.0)
        # self.model_ctau = kwargs.pop("model_ctau", 10.0)
        # self.model_xi0 = kwargs.pop("model_xi0", 1.0)
        # self.model_xiL = kwargs.pop("model_xiL", 1.0)

        super(DQCDULBDTNoJetProducer, self).__init__(*args, **kwargs)

        base = "{}/{}/src/DQCD/Modules".format(
            os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))

        if not os.getenv("_DQCDBDT"):
            os.environ["_DQCDBDT"] = "_DQCDBDT"

            ROOT.gSystem.Load("libDQCDModules.so")
            ROOT.gROOT.ProcessLine(".L {}/interface/BDTinference.h".format(base))

        if not os.getenv("_DQCDBDT_%s" % self.model):
            os.environ["_DQCDBDT_%s" % self.model] = "_DQCDBDT_%s" % self.model

            ROOT.gInterpreter.Declare("""
                auto bdt%s = BDTinference("%s", %s);
            """ % (self.model, self.model_path, use_sigmoid))

            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                std::vector<float> get_bdt_outputs_%s (
                    int nMuon, int nSV, int nmuonSV,
                    Vfloat Muon_eta, Vfloat Muon_phi, Vfloat Muon_pt, Vfloat Muon_ptErr,
                    Vfloat Muon_dxy, Vfloat Muon_dxyErr, Vfloat Muon_dz, Vfloat Muon_dzErr,
                    Vfloat Muon_ip3d, Vfloat Muon_sip3d, Vfloat Muon_charge, Vfloat Muon_tightId,
                    Vfloat Muon_softMva, Vfloat Muon_pfRelIso03_all,
                    Vfloat Muon_miniPFRelIso_all, Vfloat Muon_jetIdx,
                    Vfloat muonSV_chi2, Vfloat muonSV_pAngle, Vfloat muonSV_dlen, Vfloat muonSV_dlenSig,
                    Vfloat muonSV_dxy, Vfloat muonSV_dxySig,
                    Vfloat muonSV_mu1pt, Vfloat muonSV_mu1eta, Vfloat muonSV_mu1phi,
                    Vfloat muonSV_mu2pt, Vfloat muonSV_mu2eta, Vfloat muonSV_mu2phi,
                    Vfloat muonSV_x, Vfloat muonSV_y, Vfloat muonSV_z,
                    Vfloat SV_pt, Vfloat SV_eta, Vfloat SV_phi, 
                    Vfloat SV_x, Vfloat SV_y, Vfloat SV_z,
                    Vfloat SV_dxy, Vfloat SV_dxySig, Vfloat SV_dlen, Vfloat SV_dlenSig,
                    Vfloat SV_pAngle, Vfloat SV_chi2, Vfloat SV_ndof
                )
                {
                    std::vector<muon_legacy_t> muons(std::max(8, nMuon), muon_legacy_t());
                    for (int i = 0; i < nMuon; i++) {
                        muons[i] = muon_legacy_t({
                            Muon_eta[i], Muon_phi[i], Muon_pt[i],
                            Muon_ptErr[i], Muon_dxy[i], Muon_dxyErr[i], Muon_dz[i], Muon_dzErr[i],
                            Muon_ip3d[i], Muon_sip3d[i], Muon_charge[i], Muon_tightId[i],
                            Muon_softMva[i], Muon_pfRelIso03_all[i], Muon_miniPFRelIso_all[i],
                            Muon_jetIdx[i]
                        });
                    }
                    std::stable_sort(muons.begin(), muons.end(), muon_legacy_sort);

                    std::vector<muonsv_legacy_t> muonsvs(std::max(8, nmuonSV), muonsv_legacy_t());
                    for (int i = 0; i < nmuonSV; i++) {
                        auto muonsv_deltar = reco::deltaR(muonSV_mu1eta[i], muonSV_mu1phi[i],
                            muonSV_mu2eta[i], muonSV_mu2phi[i]);
                        muonsvs[i] = muonsv_legacy_t({
                            muonSV_chi2[i], muonSV_pAngle[i], muonSV_dlen[i], muonSV_dlenSig[i],
                            muonSV_dxy[i], muonSV_dxySig[i],
                            muonSV_mu1pt[i], muonSV_mu1eta[i], muonSV_mu1phi[i],
                            muonSV_mu2pt[i], muonSV_mu2eta[i], muonSV_mu2phi[i],
                            muonSV_x[i], muonSV_y[i], muonSV_z[i]
                        });
                    }
                    std::stable_sort(muonsvs.begin(), muonsvs.end(), muonsv_ul_sort);

                    std::vector<sv_legacy_t> svs(std::max(8, nSV), sv_legacy_t());
                    for (int i = 0; i < nSV; i++) {
                        svs[i] = sv_legacy_t({
                            SV_pt[i], SV_eta[i], SV_phi[i],
                            SV_x[i], SV_y[i], SV_z[i], SV_dxy[i], SV_dxySig[i],
                            SV_dlen[i], SV_dlenSig[i], SV_pAngle[i], SV_chi2[i], SV_ndof[i]
                        });
                    }
                    std::stable_sort(svs.begin(), svs.end(), sv_legacy_sort);

                    return bdt%s.get_bdt_outputs({
                        (float) nMuon, (float) nSV, (float) nmuonSV,
                        muons[0].eta, muons[1].eta, muons[2].eta, muons[3].eta, 
                        muons[4].eta, muons[5].eta, muons[6].eta, muons[7].eta, 
                        muons[0].phi, muons[1].phi, muons[2].phi, muons[3].phi, 
                        muons[4].phi, muons[5].phi, muons[6].phi, muons[7].phi, 
                        muons[0].pt, muons[1].pt, muons[2].pt, muons[3].pt, 
                        muons[4].pt, muons[5].pt, muons[6].pt, muons[7].pt, 
                        muons[0].ptErr, muons[1].ptErr, muons[2].ptErr, muons[3].ptErr, 
                        muons[4].ptErr, muons[5].ptErr, muons[6].ptErr, muons[7].ptErr, 
                        muons[0].dxy, muons[1].dxy, muons[2].dxy, muons[3].dxy, 
                        muons[4].dxy, muons[5].dxy, muons[6].dxy, muons[7].dxy, 
                        muons[0].dxyErr, muons[1].dxyErr, muons[2].dxyErr, muons[3].dxyErr, 
                        muons[4].dxyErr, muons[5].dxyErr, muons[6].dxyErr, muons[7].dxyErr, 
                        muons[0].dz, muons[1].dz, muons[2].dz, muons[3].dz, 
                        muons[4].dz, muons[5].dz, muons[6].dz, muons[7].dz, 
                        muons[0].dzErr, muons[1].dzErr, muons[2].dzErr, muons[3].dzErr, 
                        muons[4].dzErr, muons[5].dzErr, muons[6].dzErr, muons[7].dzErr, 
                        muons[0].ip3d, muons[1].ip3d, muons[2].ip3d, muons[3].ip3d, 
                        muons[4].ip3d, muons[5].ip3d, muons[6].ip3d, muons[7].ip3d, 
                        muons[0].sip3d, muons[1].sip3d, muons[2].sip3d, muons[3].sip3d, 
                        muons[4].sip3d, muons[5].sip3d, muons[6].sip3d, muons[7].sip3d, 
                        muons[0].charge, muons[1].charge, muons[2].charge, muons[3].charge, 
                        muons[4].charge, muons[5].charge, muons[6].charge, muons[7].charge, 
                        muons[0].tightId, muons[1].tightId, muons[2].tightId, muons[3].tightId, 
                        muons[4].tightId, muons[5].tightId, muons[6].tightId, muons[7].tightId, 
                        muons[0].softMva, muons[1].softMva, muons[2].softMva, muons[3].softMva, 
                        muons[4].softMva, muons[5].softMva, muons[6].softMva, muons[7].softMva, 
                        muons[0].pfRelIso03_all, muons[1].pfRelIso03_all,
                        muons[2].pfRelIso03_all, muons[3].pfRelIso03_all, 
                        muons[4].pfRelIso03_all, muons[5].pfRelIso03_all,
                        muons[6].pfRelIso03_all, muons[7].pfRelIso03_all, 
                        muons[0].miniPFRelIso_all, muons[1].miniPFRelIso_all,
                        muons[2].miniPFRelIso_all, muons[3].miniPFRelIso_all, 
                        muons[4].miniPFRelIso_all, muons[5].miniPFRelIso_all,
                        muons[6].miniPFRelIso_all, muons[7].miniPFRelIso_all, 
                        muons[0].jetIdx, muons[1].jetIdx, muons[2].jetIdx, muons[3].jetIdx, 
                        muons[4].jetIdx, muons[5].jetIdx, muons[6].jetIdx, muons[7].jetIdx, 
                        muonsvs[0].chi2, muonsvs[1].chi2, muonsvs[2].chi2, muonsvs[3].chi2, 
                        muonsvs[4].chi2, muonsvs[5].chi2, muonsvs[6].chi2, muonsvs[7].chi2, 
                        muonsvs[0].pAngle, muonsvs[1].pAngle, muonsvs[2].pAngle, muonsvs[3].pAngle, 
                        muonsvs[4].pAngle, muonsvs[5].pAngle, muonsvs[6].pAngle, muonsvs[7].pAngle, 
                        muonsvs[0].dlen, muonsvs[1].dlen, muonsvs[2].dlen, muonsvs[3].dlen, 
                        muonsvs[4].dlen, muonsvs[5].dlen, muonsvs[6].dlen, muonsvs[7].dlen, 
                        muonsvs[0].dlenSig, muonsvs[1].dlenSig, muonsvs[2].dlenSig, muonsvs[3].dlenSig, 
                        muonsvs[4].dlenSig, muonsvs[5].dlenSig, muonsvs[6].dlenSig, muonsvs[7].dlenSig, 
                        muonsvs[0].dxy, muonsvs[1].dxy, muonsvs[2].dxy, muonsvs[3].dxy, 
                        muonsvs[4].dxy, muonsvs[5].dxy, muonsvs[6].dxy, muonsvs[7].dxy, 
                        muonsvs[0].dxySig, muonsvs[1].dxySig, muonsvs[2].dxySig, muonsvs[3].dxySig, 
                        muonsvs[4].dxySig, muonsvs[5].dxySig, muonsvs[6].dxySig, muonsvs[7].dxySig, 
                        muonsvs[0].mu1pt, muonsvs[1].mu1pt, muonsvs[2].mu1pt, muonsvs[3].mu1pt, 
                        muonsvs[4].mu1pt, muonsvs[5].mu1pt, muonsvs[6].mu1pt, muonsvs[7].mu1pt, 
                        muonsvs[0].mu1eta, muonsvs[1].mu1eta, muonsvs[2].mu1eta, muonsvs[3].mu1eta, 
                        muonsvs[4].mu1eta, muonsvs[5].mu1eta, muonsvs[6].mu1eta, muonsvs[7].mu1eta, 
                        muonsvs[0].mu1phi, muonsvs[1].mu1phi, muonsvs[2].mu1phi, muonsvs[3].mu1phi, 
                        muonsvs[4].mu1phi, muonsvs[5].mu1phi, muonsvs[6].mu1phi, muonsvs[7].mu1phi, 
                        muonsvs[0].mu2pt, muonsvs[1].mu2pt, muonsvs[2].mu2pt, muonsvs[3].mu2pt, 
                        muonsvs[4].mu2pt, muonsvs[5].mu2pt, muonsvs[6].mu2pt, muonsvs[7].mu2pt, 
                        muonsvs[0].mu2eta, muonsvs[1].mu2eta, muonsvs[2].mu2eta, muonsvs[3].mu2eta, 
                        muonsvs[4].mu2eta, muonsvs[5].mu2eta, muonsvs[6].mu2eta, muonsvs[7].mu2eta, 
                        muonsvs[0].mu2phi, muonsvs[1].mu2phi, muonsvs[2].mu2phi, muonsvs[3].mu2phi, 
                        muonsvs[4].mu2phi, muonsvs[5].mu2phi, muonsvs[6].mu2phi, muonsvs[7].mu2phi, 
                        muonsvs[0].x, muonsvs[1].x, muonsvs[2].x, muonsvs[3].x, 
                        muonsvs[4].x, muonsvs[5].x, muonsvs[6].x, muonsvs[7].x, 
                        muonsvs[0].y, muonsvs[1].y, muonsvs[2].y, muonsvs[3].y, 
                        muonsvs[4].y, muonsvs[5].y, muonsvs[6].y, muonsvs[7].y, 
                        muonsvs[0].z, muonsvs[1].z, muonsvs[2].z, muonsvs[3].z, 
                        muonsvs[4].z, muonsvs[5].z, muonsvs[6].z, muonsvs[7].z, 
                        svs[0].pt, svs[1].pt, svs[2].pt, svs[3].pt, 
                        svs[4].pt, svs[5].pt, svs[6].pt, svs[7].pt, 
                        svs[0].eta, svs[1].eta, svs[2].eta, svs[3].eta, 
                        svs[4].eta, svs[5].eta, svs[6].eta, svs[7].eta, 
                        svs[0].phi, svs[1].phi, svs[2].phi, svs[3].phi, 
                        svs[4].phi, svs[5].phi, svs[6].phi, svs[7].phi, 
                        svs[0].x, svs[1].x, svs[2].x, svs[3].x, 
                        svs[4].x, svs[5].x, svs[6].x, svs[7].x, 
                        svs[0].y, svs[1].y, svs[2].y, svs[3].y, 
                        svs[4].y, svs[5].y, svs[6].y, svs[7].y, 
                        svs[0].z, svs[1].z, svs[2].z, svs[3].z, 
                        svs[4].z, svs[5].z, svs[6].z, svs[7].z, 
                        svs[0].dxy, svs[1].dxy, svs[2].dxy, svs[3].dxy, 
                        svs[4].dxy, svs[5].dxy, svs[6].dxy, svs[7].dxy, 
                        svs[0].dxySig, svs[1].dxySig, svs[2].dxySig, svs[3].dxySig, 
                        svs[4].dxySig, svs[5].dxySig, svs[6].dxySig, svs[7].dxySig, 
                        svs[0].dlen, svs[1].dlen, svs[2].dlen, svs[3].dlen, 
                        svs[4].dlen, svs[5].dlen, svs[6].dlen, svs[7].dlen, 
                        svs[0].dlenSig, svs[1].dlenSig, svs[2].dlenSig, svs[3].dlenSig, 
                        svs[4].dlenSig, svs[5].dlenSig, svs[6].dlenSig, svs[7].dlenSig, 
                        svs[0].pAngle, svs[1].pAngle, svs[2].pAngle, svs[3].pAngle, 
                        svs[4].pAngle, svs[5].pAngle, svs[6].pAngle, svs[7].pAngle, 
                        svs[0].chi2, svs[1].chi2, svs[2].chi2, svs[3].chi2, 
                        svs[4].chi2, svs[5].chi2, svs[6].chi2, svs[7].chi2, 
                        svs[0].ndof, svs[1].ndof, svs[2].ndof, svs[3].ndof, 
                        svs[4].ndof, svs[5].ndof, svs[6].ndof, svs[7].ndof
                    });
                }
            """ % (self.model, self.model))


    def run(self, df):
        s = randomize("bdt")

        # Some fixes needed
        #   - muonSV_mu*pt systematics (probably need to use mu*index and extract the pt from it)
        #   - SV systematics?

        df = df.Define(s, f"""!((muonSV_mu1pt[muonSV_mu1pt > 5].size() + muonSV_mu2pt[muonSV_mu2pt > 5].size() > 0))
           ? std::vector<float>(1, -1.)
           : get_bdt_outputs_{self.model}(
                nMuon, nSV, nmuonSV,
                Muon_eta, Muon_phi, Muon_pt{self.muon_syst}, Muon_ptErr,
                Muon_dxy, Muon_dxyErr, Muon_dz, Muon_dzErr,
                Muon_ip3d, Muon_sip3d, Muon_charge, Muon_tightId,
                Muon_softMva, Muon_pfRelIso03_all,
                Muon_miniPFRelIso_all, Muon_jetIdx,
                muonSV_chi2, muonSV_pAngle, muonSV_dlen, muonSV_dlenSig,
                muonSV_dxy, muonSV_dxySig,
                muonSV_mu1pt, muonSV_mu1eta, muonSV_mu1phi,
                muonSV_mu2pt, muonSV_mu2eta, muonSV_mu2phi,
                muonSV_x, muonSV_y, muonSV_z,
                SV_pt, SV_eta, SV_phi, 
                SV_x, SV_y, SV_z,
                SV_dxy, SV_dxySig, SV_dlen, SV_dlenSig,
                SV_pAngle, SV_chi2, SV_ndof
        )""")

        # p = [self.bdt_name, self.model_m, self.model_ctau, self.model_xi0, self.model_xiL]
        p = [self.bdt_name]
        p = [str(param).replace(".", "p") for param in p]

        # b_name = (f"{p[0]}_m_{p[1]}_ctau_{p[2]}_xi0_{p[3]}_xiL_{p[4]}{self.systs}")
        b_name = (f"{p[0]}{self.systs}")
        df = df.Define(b_name, " %s.at(0)" % s)

        return df, [b_name]

def DQCDULBDTNoJet(*args, **kwargs):
    """
    Returns the DQCD UL BDT output.

    Lepton and jet systematics (used for pt and mass variables) can be modified using the parameters
    from :ref:`BaseModules_JetLepMetSyst`.

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: DQCDULBDT
            path: DQCD.Modules.BDTULinference
            parameters:
                isMC: self.dataset.process.isMC
                scenario: "A"
                # model_m: 0
                # model_ctau: 0
                # model_xi0: 0
                # model_xiL: 0

    """
    return lambda: DQCDULBDTNoJetProducer(*args, **kwargs)
