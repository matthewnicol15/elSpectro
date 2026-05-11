#include "Interface.h"
#include "LorentzVector.h"
#include "DistTF1.h"
#include "FunctionsForElectronScattering.h"
#include "LundWriter.h"
#include "TwoBody_stu.h"
#include "DecayModelst_NoPS.h"

#include <TBenchmark.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

// ---------------------------------------------------------------------------
// Diagnostic histograms
// ---------------------------------------------------------------------------

double minMass = 0.2;
double maxMass = 3;
TH1F hQ2("Q2", "Q2", 1000, 0, 5);
TH1D heE("eE", "eE", 1000, 0, 20);
TH1D heTh("eTh", "eTh", 1000, 0, 180);
TH1D hKp1Th("hKp1Th", "hKp1Th", 100, 0, 180);
TH1D hKp2Th("hKp2Th", "hKp2Th", 100, 0, 180);
TH2D hKpTh("hKpTh", "hKpTh", 100, 0, 180, 100, 0, 180);
TH1D hYTh("YTh", "YTh", 1000, 0, 180);
TH1F hW("W", "W", 1000, 2, 5);
TH1F hV("Vertex", "Vertex", 100, 0, 20);
TH1F ht("t", "t", 1000, 0, 10);
TH1F hgE("gE", "gE", 1000, 0, 20);
TH1F hgTh("gTh", "gTh", 1000, 0, 180);
TH1F hLM("HyperonM", "; Y* to #Xi K+ Mass (GeV)", 1000, 1.0, 4);
TH2F hHypW("Hyperon_v_W", "; Y* to #Xi K+ Mass (GeV)", 1000, 1.0, 4, 1000, 1.0, 4);
TH1F hPKp("PKp", "K+", 100, 0, 10);
TH1F hPXi("PXi", "#Xi-", 100, 0, 10);
TH1F hPKp2("PK", "K+", 100, 0, 10);
TH1F hMcap("Mcap", ";W - m_{K} (GeV)", 1000, 1, 5);

void ep_to_KKXi_broad(double ebeamE, int nEvents, int fileno)
{
  using namespace elSpectro;
  elSpectro::Manager::Instance();
  gRandom->SetSeed(0);

  // define e- beam, pdg =11 momentum = _beamP
  auto elBeam = initial(11, ebeamE);
  auto elbeam = elBeam->GetInteracting4Vector();
  // proton target at rest
  auto prTarget = initial(2212, 0);
  auto prbeam = prTarget->GetInteracting4Vector();

  // ep -> e' K+ Y* (~3 GeV with 2 GeV width) -> e' K+ K+ Xi-

  // Hyperon Y* -> K+ Xi-
  double thr = 1.32171 + 0.493677;
  double wmax = (*elbeam + *prbeam).M();
  std::cout << "wmax " << wmax << std::endl;
  auto dist1 = new DistTF1{TF1(Form("Ybw_%d", (int)ebeamE), "TMath::BreitWigner(x,3.0,2.0)", thr, wmax)};
  mass_distribution(9995, new DistTF1{TF1("hh", "TMath::BreitWigner(x,3.0,2.0)", thr, (*elbeam + *prbeam).M())});
  auto hyperon = static_cast<DecayingParticle *>(particle(9995, model(new GenericModelst{dist1, {}, {3312, 321}})));

  // decay of pGamma* to K+ Y*
  auto pGammaStarDecay = static_cast<DecayModelst *>(model(new DecayModelst_NoPS{{hyperon}, {321}}));
  //
  // create mesonex electroproduction of X + neutron
  // TwoBody_stu{0.1, 0.9, 3 ,0,0} //0.1 strength  s distribution (flat angular dist.),  0.9 strength t distribution with slope b = 3
  mesonex(elBeam, prTarget, new DecayModelQ2W{0, pGammaStarDecay, new TwoBody_stu{0.1, 0.9, 3, 0, 0}});

  // give limits to the detected electron
  auto production = dynamic_cast<ElectronScattering *>(generator().Reaction());
  // production->SetLimitTarRest_eThmin(3.5*TMath::DegToRad());
  // production->SetLimitTarRest_eThmax(5.5*TMath::DegToRad());
  // production->SetLimitTarRest_ePmin(0.4);
  // production->SetLimitTarRest_ePmax(6);

  // get pointers to produced particles fror diagnostic histos
  auto Xi = hyperon->Model()->Product(1);
  auto Kp2 = hyperon->Model()->Product(0);

  auto Kp = pGammaStarDecay->GetMeson();
  auto electron = dynamic_cast<DecayModelQ2W *>(production->Model())->GetScatteredElectron();

  // ---------------------------------------------------------------------------
  // Initialize LUND
  // ---------------------------------------------------------------------------

  writer(new LundWriter{Form("/home/nics/work/York/elSpectro/scripts/outputs/ep_to_KKXi_2000MeV_Mcap_%d_%d.dat", (int)ebeamE, fileno)});

  // initilase the generator, may take some time for making distribution tables
  initGenerator();

  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------

  gBenchmark->Start("e");

  for (int i = 0; i < nEvents; i++)
  {
    // for(int i=0;i<0;i++){
    nextEvent();

    // fill diagnostic histograms
    auto photon = *elbeam - electron->P4();
    double Q2 = -photon.M2();
    double W = (photon + *prbeam).M();
    double t = -1 * (Kp2->P4() - *prbeam).M2();

    double Mcap = W - 0.493677;
    if (Mcap < 4.0)
    {
      i--;
      continue;
    }

    if (i % 1000 == 0)
      std::cout << "event number " << i << std::endl;

    hQ2.Fill(Q2);
    hW.Fill(W);
    hMcap.Fill(Mcap);
    ht.Fill(t);
    hgE.Fill(photon.E());

    auto elec = electron->P4();
    heTh.Fill(elec.Theta() * TMath::RadToDeg());
    heE.Fill(elec.E());

    hKp1Th.Fill(Kp->P4().Theta() * TMath::RadToDeg());
    hKp2Th.Fill(Kp2->P4().Theta() * TMath::RadToDeg());
    hKpTh.Fill(Kp->P4().Theta() * TMath::RadToDeg(), Kp2->P4().Theta() * TMath::RadToDeg());

    auto Lrec = Kp2->P4() + Xi->P4();
    hLM.Fill(Lrec.M());
    hHypW.Fill(Lrec.M(), W);

    hPKp.Fill(Kp->P4().P());
    hPXi.Fill(Xi->P4().P());
    hPKp2.Fill(Kp2->P4().P());

    hV.Fill(Kp2->VertexPosition()->P() / 10);
  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");

  // internally stored histograms for total ep cross section
  TH1D *hWdist = (TH1D *)gDirectory->FindObject("Wdist");
  TH1D *hGenWdist = (TH1D *)gDirectory->FindObject("genWdist");

  TFile *fout = TFile::Open(Form("/home/nics/work/York/elSpectro/scripts/outputs/ep_to_KKXi_2000MeV_Mcap_%d_%d.root", (int)ebeamE, fileno), "recreate");
  // total ep cross section inputs
  if (hWdist)
    hWdist->Write();
  if (hGenWdist)
    hGenWdist->Write();

  // generated event distributions
  hQ2.Write();
  hW.Write();
  ht.Write();
  hgE.Write();
  heTh.Write();
  heE.Write();
  hLM.Write();
  hPKp.Write();
  hPXi.Write();
  hPKp2.Write();
  hHypW.Write();
  hMcap.Write();
  hKp1Th.Write();
  hKp2Th.Write();
  hKpTh.Write();
  fout->Close();

  generator().Summary();
}
