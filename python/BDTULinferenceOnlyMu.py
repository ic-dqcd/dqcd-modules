import os

from analysis_tools.utils import import_root, randomize
from Base.Modules.baseModules import JetLepMetSyst

ROOT = import_root()


class DQCDULBDTOnlyMuProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        scenario = kwargs.pop("scenario", "A")
        default_name = "bdt_scenario%s" % scenario
        use_sigmoid = "false"
        if scenario == "A":
            default_model_path = os.path.expandvars(
                "$CMSSW_BASE/src/DQCD/Modules/data/XGB-ONLYMU_46_scenarioA.model")
        elif scenario == "B1":
            default_model_path = os.path.expandvars(
                # "$CMSSW_BASE/src/DQCD/Modules/data/model_ul_saved_scenarioB1.model")
                "$CMSSW_BASE/src/DQCD/Modules/data/XGB-ONLYMU_41_scenarioB1.model")
        elif scenario == "B2":
            raise ValueError("B2 onlymu not implemented")
            default_model_path = os.path.expandvars(
                "$CMSSW_BASE/src/DQCD/Modules/data/model_ul_saved_scenarioB2_new.model")
        elif scenario == "C":
            raise ValueError("C onlymu not implemented")
            default_model_path = os.path.expandvars(
                "$CMSSW_BASE/src/DQCD/Modules/data/model_ul_saved_scenarioC_new.model")
        elif scenario == "vector":
            default_model_path = os.path.expandvars(
                "$CMSSW_BASE/src/DQCD/Modules/data/XGB-ONLYMU_108_vector.model")
            default_name = "bdt_vector"
            # use_sigmoid = "true"
        else:
            raise ValueError("Only BDT for scenario A, B1, and vector portal are already implemented")

        self.model_path = kwargs.pop("model_path", default_model_path)
        self.model = self.model_path.replace("/", "_").replace(".", "_")
        self.bdt_name = kwargs.pop("bdt_name", default_name)
        # self.model_m = kwargs.pop("model_m", 2.0)
        # self.model_ctau = kwargs.pop("model_ctau", 10.0)
        # self.model_xi0 = kwargs.pop("model_xi0", 1.0)
        # self.model_xiL = kwargs.pop("model_xiL", 1.0)

        super(DQCDULBDTOnlyMuProducer, self).__init__(*args, **kwargs)

        base = "{}/{}/src/DQCD/Modules".format(
            os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))

        if not os.getenv("_DQCDBDT"):
            os.environ["_DQCDBDT"] = "_DQCDBDT"

            ROOT.gSystem.Load("libDQCDModules.so")
            ROOT.gROOT.ProcessLine(".L {}/interface/BDTinference.h".format(base))

        if not os.getenv("_DQCDBDTOnlyMu_%s" % self.model):
            os.environ["_DQCDBDTOnlyMu_%s" % self.model] = "_DQCDBDT_%s" % self.model

            ROOT.gInterpreter.Declare("""
                auto bdt%s = BDTinference("%s", %s);
            """ % (self.model, self.model_path, use_sigmoid))

            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                std::vector<float> get_bdt_outputs_%s (
                    int nMuon, int nmuonSV,
                    Vfloat Muon_eta, Vfloat Muon_phi, Vfloat Muon_pt, Vfloat Muon_ptErr,
                    Vfloat Muon_dxy, Vfloat Muon_dxyErr, Vfloat Muon_dz, Vfloat Muon_dzErr,
                    Vfloat Muon_ip3d, Vfloat Muon_sip3d, Vfloat Muon_charge, Vfloat Muon_tightId,
                    Vfloat Muon_softMva, Vfloat Muon_pfRelIso03_all,
                    Vfloat Muon_miniPFRelIso_all, Vfloat Muon_jetIdx,
                    Vfloat muonSV_chi2, Vfloat muonSV_pAngle,
                    Vfloat muonSV_dlen, Vfloat muonSV_dlenSig,
                    Vfloat muonSV_dxy, Vfloat muonSV_dxySig,
                    Vfloat muonSV_mu1pt, Vfloat muonSV_mu1eta, Vfloat muonSV_mu1phi,
                    Vfloat muonSV_mu2pt, Vfloat muonSV_mu2eta, Vfloat muonSV_mu2phi,
                    Vfloat muonSV_x, Vfloat muonSV_y, Vfloat muonSV_z
                )
                {
                    std::vector<muon_legacy_t> muons(std::max(4, nMuon), muon_legacy_t());
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

                    std::vector<muonsv_legacy_t> muonsvs(std::max(2, nmuonSV), muonsv_legacy_t());
                    for (int i = 0; i < nmuonSV; i++) {
                        // auto muonsv_deltar = reco::deltaR(muonSV_mu1eta[i], muonSV_mu1phi[i],
                        //    muonSV_mu2eta[i], muonSV_mu2phi[i]);
                        muonsvs[i] = muonsv_legacy_t({
                            muonSV_chi2[i], muonSV_pAngle[i], muonSV_dlen[i], muonSV_dlenSig[i],
                            muonSV_dxy[i], muonSV_dxySig[i],
                            muonSV_mu1pt[i], muonSV_mu1eta[i], muonSV_mu1phi[i],
                            muonSV_mu2pt[i], muonSV_mu2eta[i], muonSV_mu2phi[i],
                            muonSV_x[i], muonSV_y[i], muonSV_z[i], 1.
                        });
                    }
                    std::stable_sort(muonsvs.begin(), muonsvs.end(), muonsv_ul_sort);

                    return bdt%s.get_bdt_outputs({
                        (float) nMuon,
                        muons[0].eta, muons[1].eta, muons[2].eta, muons[3].eta,
                        muons[0].phi, muons[1].phi, muons[2].phi, muons[3].phi,
                        muons[0].pt, muons[1].pt, muons[2].pt, muons[3].pt,
                        muons[0].ptErr, muons[1].ptErr, muons[2].ptErr, muons[3].ptErr,
                        muons[0].dxy, muons[1].dxy, muons[2].dxy, muons[3].dxy,
                        muons[0].dxyErr, muons[1].dxyErr, muons[2].dxyErr, muons[3].dxyErr,
                        muons[0].dz, muons[1].dz, muons[2].dz, muons[3].dz,
                        muons[0].dzErr, muons[1].dzErr, muons[2].dzErr, muons[3].dzErr,
                        muons[0].ip3d, muons[1].ip3d, muons[2].ip3d, muons[3].ip3d,
                        muons[0].sip3d, muons[1].sip3d, muons[2].sip3d, muons[3].sip3d,
                        muons[0].charge, muons[1].charge, muons[2].charge, muons[3].charge,
                        muons[0].tightId, muons[1].tightId, muons[2].tightId, muons[3].tightId,
                        muons[0].softMva, muons[1].softMva, muons[2].softMva, muons[3].softMva,
                        muons[0].pfRelIso03_all, muons[1].pfRelIso03_all,
                        muons[2].pfRelIso03_all, muons[3].pfRelIso03_all,
                        muons[0].miniPFRelIso_all, muons[1].miniPFRelIso_all,
                        muons[2].miniPFRelIso_all, muons[3].miniPFRelIso_all,
                        muons[0].jetIdx, muons[1].jetIdx, muons[2].jetIdx, muons[3].jetIdx,
                        muonsvs[0].chi2, muonsvs[1].chi2,
                        muonsvs[0].pAngle, muonsvs[1].pAngle,
                        muonsvs[0].dlen, muonsvs[1].dlen,
                        muonsvs[0].dlenSig, muonsvs[1].dlenSig,
                        muonsvs[0].dxy, muonsvs[1].dxy,
                        muonsvs[0].dxySig, muonsvs[1].dxySig,
                        muonsvs[0].mu1pt, muonsvs[1].mu1pt,
                        muonsvs[0].mu1eta, muonsvs[1].mu1eta,
                        muonsvs[0].mu1phi, muonsvs[1].mu1phi,
                        muonsvs[0].mu2pt, muonsvs[1].mu2pt,
                        muonsvs[0].mu2eta, muonsvs[1].mu2eta,
                        muonsvs[0].mu2phi, muonsvs[1].mu2phi,
                        muonsvs[0].x, muonsvs[1].x,
                        muonsvs[0].y, muonsvs[1].y,
                        muonsvs[0].z, muonsvs[1].z
                    });
                }
            """ % (self.model, self.model))


    def run(self, df):
        s = randomize("bdt")

        # Some fixes needed
        #   - muonSV_mu*pt systematics (probably need to use mu*index and extract the pt from it)
        #   - SV systematics?

        df = df.Define(s, f"""!(((HLT_Mu9_IP6_part0 == 1) ||
                (HLT_Mu9_IP6_part1 == 1) || (HLT_Mu9_IP6_part2 == 1) || (HLT_Mu9_IP6_part3 == 1) ||
                (HLT_Mu9_IP6_part4 == 1)) && 
            (muonSV_mu1pt[muonSV_mu1pt > 5].size() + muonSV_mu2pt[muonSV_mu2pt > 5].size() > 0))
            ? std::vector<float>(1, -1.)
            : get_bdt_outputs_{self.model}(
                nMuon, nmuonSV,
                Muon_eta, Muon_phi, Muon_pt{self.muon_syst}, Muon_ptErr,
                Muon_dxy, Muon_dxyErr, Muon_dz, Muon_dzErr,
                Muon_ip3d, Muon_sip3d, Muon_charge, Muon_tightId,
                Muon_softMva, Muon_pfRelIso03_all,
                Muon_miniPFRelIso_all, Muon_jetIdx,
                muonSV_chi2, muonSV_pAngle, muonSV_dlen, muonSV_dlenSig,
                muonSV_dxy, muonSV_dxySig,
                muonSV_mu1pt, muonSV_mu1eta, muonSV_mu1phi,
                muonSV_mu2pt, muonSV_mu2eta, muonSV_mu2phi,
                muonSV_x, muonSV_y, muonSV_z
        )""")

        # p = [self.bdt_name, self.model_m, self.model_ctau, self.model_xi0, self.model_xiL]
        p = [self.bdt_name]
        p = [str(param).replace(".", "p") for param in p]

        # b_name = (f"{p[0]}_m_{p[1]}_ctau_{p[2]}_xi0_{p[3]}_xiL_{p[4]}{self.systs}")
        b_name = (f"{p[0]}{self.systs}")
        df = df.Define(b_name, " %s.at(0)" % s)

        return df, [b_name]


def DQCDULBDTOnlyMu(*args, **kwargs):
    """
    Returns the DQCD UL BDT output.

    Lepton and jet systematics (used for pt and mass variables) can be modified using the parameters
    from :ref:`BaseModules_JetLepMetSyst`.

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: DQCDULBDTOnlyMu
            path: DQCD.Modules.BDTULinference
            parameters:
                isMC: self.dataset.process.isMC
                scenario: "A"
                # model_m: 0
                # model_ctau: 0
                # model_xi0: 0
                # model_xiL: 0

    """
    return lambda: DQCDULBDTOnlyMuProducer(*args, **kwargs)
