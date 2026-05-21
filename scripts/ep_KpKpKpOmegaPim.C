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
TH1F hQ2("Q2", "Q2", 1000, 0, 5);
TH1D heE("eE", "eE", 1000, 0, 20);
TH1D heTh("eTh", "eTh", 1000, 0, 180);
TH1D hKp1Th("hKp1Th", "hKp1Th", 100, 0, 180);
TH1D hKp2Th("hKp2Th", "hKp2Th", 100, 0, 180);
TH1D hKp2Th("hKp3Th", "hKp3Th", 100, 0, 180);
TH1D hKpAcceptance("hKpAcceptance", "hKpAcceptance", 4, 0, 3);
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

void ep_KpKpKpOmegaPim(double ebeamE, int nEvents, int fileno)
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

  auto Kp_rest = LorentzVector(0, 0, 0, 0.493677); // K+ at rest

  // ep -> e' K+ Y(9995, broad) -> e' K+ K+ Y2(9996, broad)
  //    -> e' K+ K+ K+ Y3(9997, Omega pi-) -> e' K+ K+ K+ Omega pi-
  double wmax = (*elbeam + *prbeam).M();
  std::cout << "wmax " << wmax << std::endl;

  // Intermediate state for [Omega pi-] to keep all decays two-body
  // 9997 -> Omega- pi-
  const double thr3 = 1.67245 + 0.13957;
  double max3 = wmax - 3.0 * 0.493677;
  if (max3 <= thr3)
    max3 = thr3 + 1E-3;
  auto dist3 = new DistTF1{TF1(Form("OmegaPiBW_%d", (int)ebeamE), "TMath::BreitWigner(x,2.3,1.5)", thr3, max3)};
  mass_distribution(9997, new DistTF1{TF1("OmegaPiMass", "TMath::BreitWigner(x,2.3,1.5)", thr3, max3)});
  auto omegaPi = static_cast<DecayingParticle *>(particle(9997, model(new GenericModelst{dist3, {}, {3334, -211}})));

  // Xi*(9996) -> K+ (9995)[Omega pi-]
  const double thr2 = 0.493677 + thr3;
  double max2 = wmax - 2.0 * 0.493677;
  if (max2 <= thr2)
    max2 = thr2 + 1E-3;
  auto dist2 = new DistTF1{TF1(Form("XiStarBW_%d", (int)ebeamE), "TMath::BreitWigner(x,3.0,2.0)", thr2, max2)};
  mass_distribution(9996, new DistTF1{TF1("XiStarMass", "TMath::BreitWigner(x,3.0,2.0)", thr2, max2)});
  auto hyperon2 = static_cast<DecayingParticle *>(particle(9996, model(new GenericModelst{dist2, {omegaPi}, {321}})));

  // Y*(9995) -> K+ Xi*(9996)
  const double thr1 = 0.493677 + thr2;
  double max1 = wmax - 1.0 * 0.493677;
  if (max1 <= thr1)
    max1 = thr1 + 1E-3;
  auto dist1 = new DistTF1{TF1(Form("YStarBW_%d", (int)ebeamE), "TMath::BreitWigner(x,4.0,3.0)", thr1, max1)};
  mass_distribution(9995, new DistTF1{TF1("YStarMass", "TMath::BreitWigner(x,4.0,3.0)", thr1, max1)});
  auto hyperon1 = static_cast<DecayingParticle *>(particle(9995, model(new GenericModelst{dist1, {hyperon2}, {321}})));

  // decay of pGamma* to K+ Y*
  auto pGammaStarDecay = static_cast<DecayModelst *>(model(new DecayModelst_NoPS{{hyperon1}, {321}}));
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
  auto XiStar = hyperon1->Model()->Product(0);
  auto Kp2 = hyperon1->Model()->Product(1);

  auto omegaPiReco = dynamic_cast<const DecayingParticle *>(hyperon2->Model()->Product(0));
  auto Kp3 = hyperon2->Model()->Product(1);

  auto Omega = omegaPiReco->Model()->Product(0);
  auto pim = omegaPiReco->Model()->Product(1);

  auto Kp = pGammaStarDecay->GetMeson();
  auto electron = dynamic_cast<DecayModelQ2W *>(production->Model())->GetScatteredElectron();

  // ---------------------------------------------------------------------------
  // Initialize LUND
  // ---------------------------------------------------------------------------

  writer(new LundWriter{Form("/home/nics/work/York/elSpectro/scripts/outputs/ep_to_KpKpKpOmegaPim_%d_%d.dat", (int)ebeamE, fileno)});

  // initilase the generator, may take some time for making distribution tables
  initGenerator();

  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------

  gBenchmark->Start("e");

  int Acceptance = 0;

  for (int i = 0; i < nEvents; i++)
  {

    nextEvent();

    Acceptance = 0;
    // fill diagnostic histograms
    auto photon = *elbeam - electron->P4();
    double Q2 = -photon.M2();
    double W = (photon + *prbeam).M();
    double t = -1.0 * (hyperon1->P4() - *prbeam).M2();

    double Mcap = W - 0.493677;
    if (Mcap < 3.0)
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

    double Kp1Theta = Kp->P4().Theta() * TMath::RadToDeg();
    double Kp2Theta = Kp2->P4().Theta() * TMath::RadToDeg();
    double Kp3Theta = Kp3->P4().Theta() * TMath::RadToDeg();
    hKp1Th.Fill(Kp1Theta);
    hKp2Th.Fill(Kp2Theta);
    hKp3Th.Fill(Kp3Theta);

    if (Kp1Theta >= 5 && Kp1Theta <= 45)
      Acceptance++;
    if (Kp2Theta >= 5 && Kp2Theta <= 45)
      Acceptance++;
    if (Kp3Theta >= 5 && Kp3Theta <= 45)
      Acceptance++;

    hKpAcceptance.Fill(Acceptance);

    auto Lrec = Kp2->P4() + XiStar->P4();
    hLM.Fill(Lrec.M());
    hHypW.Fill(Lrec.M(), W);

    hPKp.Fill(Kp->P4().P());
    hPXi.Fill(XiStar->P4().P());
    hPKp2.Fill(Kp2->P4().P());

    hV.Fill(Kp2->VertexPosition()->P() / 10);
  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");

  // internally stored histograms for total ep cross section
  TH1D *hWdist = (TH1D *)gDirectory->FindObject("Wdist");
  TH1D *hGenWdist = (TH1D *)gDirectory->FindObject("genWdist");

  TFile *fout = TFile::Open(Form("/home/nics/work/York/elSpectro/scripts/outputs/ep_to_KpKpKpOmegaPim_%d_%d.root", (int)ebeamE, fileno), "recreate");
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
