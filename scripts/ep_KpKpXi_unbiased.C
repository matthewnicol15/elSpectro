#include "Interface.h"
#include "LorentzVector.h"
#include "DistTF1.h"
#include "FunctionsForElectronScattering.h"
#include "LundWriter.h"
#include "TwoBody_stu.h"
#include "DecayModelst_NoPS.h"
#include "PhaseSpaceDecay.h"

#include <TBenchmark.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

// ---------------------------------------------------------------------------
// Diagnostic histograms
// ---------------------------------------------------------------------------

double minMass = 0.2;
double maxMass = 3;
TH1F hQ2("Q2", "Q2", 1000, -2, 12);
TH2F hQ2_gE("Q2_gE", "", 1000, -2, 12, 1000, -2, 12);
TH1D heE("eE", "eE", 1000, 0, 20);
TH1D heTh("eTh", "eTh", 1000, 0, 180);
TH1D hKp1Th("hKp1Th", "hKp1Th", 100, 0, 180);
TH1D hKp2Th("hKp2Th", "hKp2Th", 100, 0, 180);
TH1D hKpAcceptance("KpAcceptance", "K^{+} Acceptance", 4, 0, 3);
TH1D hYTh("YTh", "YTh", 1000, 0, 180);
TH1F hW("W", "W", 1000, 2, 5);
TH1F hV("Vertex", "Vertex", 100, 0, 20);
TH1F ht("t", "t", 1000, -2, 12);
TH1F hgE("gE", "gE", 1000, -2, 12);
TH1F hgTh("gTh", "gTh", 1000, 0, 180);
TH1F hHyperon1Rec("Hyperon1Rec", "; Y*_{1} to Y*_{2} #K^{+}; Y*_{1} Mass [GeV]", 1000, 1.0, 4);
TH1F hHyperon2Rec("Hyperon2Rec", "; Y*_{2} to #Xi^{-} #pi^{0}; Y*_{2} Mass [GeV]", 1000, 1.0, 4);
TH1F hPKp("PKp", ";P_{K^{+}_{1}} [GeV]", 100, 0, 10);
TH1F hPKp2("PKp2", ";P_{K^{+}_{2}} [GeV]", 100, 0, 10);
TH1F hPPi0("PPi0", ";P_{#pi^{0}} [GeV]", 100, 0, 10);
TH1F hPXi("PXi", ";p_{#Xi^{-}} [GeV]", 100, 0, 10);
TH1F hPHyperon1("PHyperon1", ";P_{Y_{1}} [GeV]", 100, 0, 10);
TH1F hPHyperon2("PHyperon2", ";P_{Y_{2}} [GeV]", 100, 0, 10);

void ep_KpKpXi_unbiased(double ebeamE, int nEvents, int fileno)
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

  const double Xi_rest = 1.32171;
  const double Pi_rest = 0.134976;
  const double Kp_rest = 0.493677;

  // ep -> e' K+ Y(9995, broad) -> e' K+ K+ Y2(9996, broad)
  //    -> e' K+ K+ Xi- pi0
  double wmax = (*elbeam + *prbeam).M();
  std::cout << "wmax " << wmax << std::endl;

  // Xi(9996) -> K+ (9995)[Xi pi0]
  const double thr2 = Xi_rest + Pi_rest;
  double max2 = wmax - 2.0 * Kp_rest;
  if (max2 <= thr2)
    max2 = thr2 + 1E-3;
  auto dist2 = new DistTF1{TF1(Form("XiStarBW_%d", (int)ebeamE), "TMath::BreitWigner(x,3.0,3.0)", thr2, max2)};
  mass_distribution(9996, new DistTF1{TF1("XiStarMass", "TMath::BreitWigner(x,3.0,3.0)", thr2, max2)});
  auto hyperon2 = static_cast<DecayingParticle *>(particle(9996, model(new GenericModelst{dist2, {}, {3312, 111}})));

  // Y*(9995) -> K+ Xi*(9996)
  const double thr1 = Kp_rest + thr2;
  double max1 = wmax - 1.0 * Kp_rest;
  std::cout << "thr1 " << thr1 << "  max  " << max1 << std::endl;
  if (max1 <= thr1)
    max1 = thr1 + 1E-3;
  auto dist1 = new DistTF1{TF1(Form("YStarBW_%d", (int)ebeamE), "TMath::BreitWigner(x,3.0,3.0)", thr1, max1)};
  mass_distribution(9995, new DistTF1{TF1("YStarMass", "TMath::BreitWigner(x,3.0,3.0)", thr1, max1)});
  auto hyperon1 = static_cast<DecayingParticle *>(particle(9995, model(new GenericModelst{dist1, {hyperon2}, {321}})));

  // decay of pGamma* to K+ Y*
  auto pGammaStarDecay = static_cast<DecayModelst *>(model(new DecayModelst_NoPS{{hyperon1}, {321}}));
  //
  // TwoBody_stu{0.1, 0.9, 3 ,0,0} //0.1 strength  s distribution (flat angular dist.),  0.9 strength t distribution with slope b = 3
  double s_strength = 0.1;
  double t_strength = 0.9;
  double t_slope = 0.1;
  mesonex(elBeam, prTarget, new DecayModelQ2W{0, pGammaStarDecay, new TwoBody_stu{s_strength, t_strength, t_slope, 0, 0}});

  // give limits to the detected electron
  auto production = dynamic_cast<ElectronScattering *>(generator().Reaction());
  auto q2wModel = dynamic_cast<DecayModelQ2W *>(production->Model());

  production->SetLimitTarRest_eThmin(4.5 * TMath::DegToRad());
  production->SetLimitTarRest_eThmax(40 * TMath::DegToRad());

  // get pointers to produced particles fror diagnostic histos
  auto Hyper1 = pGammaStarDecay->Product(0);
  auto Kp = pGammaStarDecay->Product(1);

  auto Hyper2 = hyperon1->Model()->Product(0);
  auto Kp2 = hyperon1->Model()->Product(1);

  auto Xi = hyperon2->Model()->Product(0);
  auto Pi0 = hyperon2->Model()->Product(1);

  auto electron = dynamic_cast<DecayModelQ2W *>(production->Model())->GetScatteredElectron();

  // ---------------------------------------------------------------------------
  // Initialize LUND
  // ---------------------------------------------------------------------------

  writer(new LundWriter{Form("/home/nics/work/York/elSpectro/scripts/outputs/ep_to_KpKpXi_unbiased_%d_%d.dat", (int)ebeamE, fileno)});

  // initilase the generator, may take some time for making distribution tables
  initGenerator();

  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------

  gBenchmark->Start("e");

  int Acceptance = 0;
  int percentage = (nEvents >= 100) ? (nEvents / 100) : 1;
  int total = 0;

  for (int i = 0; i < nEvents; i++)
  {
    nextEvent();

    total++;
    // fill diagnostic histograms
    auto photon = *elbeam - electron->P4();
    double Q2 = -photon.M2();

    double W = (photon + *prbeam).M();
    double t = -1.0 * (photon - Kp->P4()).M2();

    if (percentage > 0 && i % percentage == 0)
      std::cout << "event number " << i << std::endl;

    Acceptance = 0;

    double Kp1Theta = Kp->P4().Theta() * TMath::RadToDeg();
    double Kp2Theta = Kp2->P4().Theta() * TMath::RadToDeg();
    hKp1Th.Fill(Kp1Theta);
    hKp2Th.Fill(Kp2Theta);

    if (Kp1Theta >= 5 && Kp1Theta <= 45)
      Acceptance++;
    if (Kp2Theta >= 5 && Kp2Theta <= 45)
      Acceptance++;

    hKpAcceptance.Fill(Acceptance);

    hQ2.Fill(Q2);
    hW.Fill(W);
    ht.Fill(t);
    hgE.Fill(photon.E());
    hQ2_gE.Fill(Q2, photon.E());

    auto elec = electron->P4();
    heTh.Fill(elec.Theta() * TMath::RadToDeg());
    heE.Fill(elec.E());

    auto hyperon1Rec = Hyper2->P4() + Kp2->P4();
    hHyperon1Rec.Fill(hyperon1Rec.M());

    auto hyperon2Rec = Xi->P4() + Pi0->P4();
    hHyperon2Rec.Fill(hyperon2Rec.M());

    hPKp.Fill(Kp->P4().P());
    hPKp2.Fill(Kp2->P4().P());
    hPHyperon1.Fill(hyperon1->P4().P());
    hPHyperon2.Fill(hyperon2->P4().P());
    hPXi.Fill(Xi->P4().P());
  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");

  // internally stored histograms for total ep cross section
  TH1D *hWdist = (TH1D *)gDirectory->FindObject("Wdist");
  TH1D *hGenWdist = (TH1D *)gDirectory->FindObject("genWdist");

  TFile *fout = TFile::Open(Form("/home/nics/work/York/elSpectro/scripts/outputs/ep_to_KpKpXi_unbiased_%d_%d.root", (int)ebeamE, fileno), "recreate");
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
  hHyperon1Rec.Write();
  hHyperon2Rec.Write();
  hPKp.Write();
  hPPi0.Write();
  hKp1Th.Write();
  hKpAcceptance.Write();
  hQ2_gE.Write();
  hPHyperon1.Write();
  hPHyperon2.Write();
  hPXi.Write();
  hPKp2.Write();
  fout->Close();

  generator().Summary();
}
