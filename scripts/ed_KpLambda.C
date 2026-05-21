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
TH1D hKpAcceptance("hKpAcceptance", "hKpAcceptance", 4, 0, 3);
TH1D hYTh("YTh", "YTh", 1000, 0, 180);
TH1F hW("W", "W", 1000, 2, 5);
TH1F hV("Vertex", "Vertex", 100, 0, 20);
TH1F ht("t", "t", 1000, 0, 10);
TH1F hgE("gE", "gE", 1000, 0, 20);
TH1F hgTh("gTh", "gTh", 1000, 0, 180);
TH1F hHyperon1Rec("hHyperon1Rec", "; Y* to p #pi^{-} Mass (GeV)", 1000, 1.0, 4);
TH2F hKp1Th_v_LambdaRec("hKp1Th_v_LambdaRec", "", 1000, 1.0, 4, 1000, 0, 180);
TH2F hHypW("Hyperon_v_W", "; Y* to #Xi K+ Mass (GeV)", 1000, 1.0, 4, 1000, 1.0, 4);
TH1F hPKp("PKp", "K+", 100, 0, 10);
TH1F hPLambda("PLambda", "#Lambda", 100, 0, 10);
TH1F hPProton("PProton", "P", 100, 0, 10);
TH1F hPPim("PPim", "#pi^{-}", 100, 0, 10);
TH1F hMcap("Mcap", ";W - m_{K} (GeV)", 1000, 1, 5);

void ed_KpLambda(double ebeamE, int nEvents, int fileno)
{
  using namespace elSpectro;
  elSpectro::Manager::Instance();
  gRandom->SetSeed(0);

  // define e- beam, pdg =11 momentum = _beamP
  auto elBeam = initial(11, ebeamE);
  auto elbeam = elBeam->GetInteracting4Vector();

  // deuteron target at rest
  // Use CDBonn potential to give nucleon momentum
  auto qfTarget = initial(2212, 0, 1000010020, model(new NuclearBreakup(2212, 2112)), new QuasiFreeNucleon(QuasiFree::CDBonnMomentum())); // CDBonn potential

  // get 4-momentum of target nucleon
  auto nucleon = qfTarget->GetInteracting4Vector();

  auto Kp_rest = LorentzVector(0, 0, 0, 0.493677); // K+ at rest

  // ed -> e' n K+ Y(9995, broad)
  double wmax = (*elbeam + *nucleon).M();
  std::cout << "wmax " << wmax << std::endl;

  // Y*(9995) -> p pi-
  const double thr1 = 0.938272 + 0.13957;
  double max1 = wmax - 1.0 * 0.493677;
  if (max1 <= thr1)
    max1 = thr1 + 1E-3;

  auto dist_Y = new DistTF1{TF1(Form("YStarBW_%d", (int)ebeamE), "TMath::BreitWigner(x,4.0,3.0)", thr1, max1)};
  mass_distribution(9995, new DistTF1{TF1("YStarMass", "TMath::BreitWigner(x,4.0,3.0)", thr1, max1)});
  auto hyperon1 = static_cast<DecayingParticle *>(particle(9995, model(new GenericModelst{dist_Y, {}, {2212, -211}})));

  // decay of pGamma* to K+ + Y*
  auto pGammaStarDecay = static_cast<DecayModelst *>(model(new DecayModelst_NoPS{{hyperon1}, {321}}));
  //
  // create mesonex electroproduction of X + neutron
  // TwoBody_stu{0.1, 0.9, 3 ,0,0} //0.1 strength  s distribution (flat angular dist.),  0.9 strength t distribution with slope b = 3
  mesonex(elBeam, qfTarget, new DecayModelQ2W{0, pGammaStarDecay, new TwoBody_stu{0.1, 0.9, 3, 0, 0}});

  // give limits to the detected electron
  auto production = dynamic_cast<ElectronScattering *>(generator().Reaction());
  // production->SetLimitTarRest_eThmin(3.5*TMath::DegToRad());
  // production->SetLimitTarRest_eThmax(5.5*TMath::DegToRad());
  // production->SetLimitTarRest_ePmin(0.4);
  // production->SetLimitTarRest_ePmax(6);

  // get pointers to produced particles fror diagnostic histos
  auto Proton = hyperon1->Model()->Product(0);
  auto Pim = hyperon1->Model()->Product(1);

  auto Kp = pGammaStarDecay->GetMeson();

  auto spec_nucleon = pGammaStarDecay->GetBaryon();
  auto electron = dynamic_cast<DecayModelQ2W *>(production->Model())->GetScatteredElectron();

  // ---------------------------------------------------------------------------
  // Initialize LUND
  // ---------------------------------------------------------------------------

  writer(new LundWriter{Form("/home/nics/work/York/elSpectro/scripts/outputs/ed_to_KpLambda_%d_%d.dat", (int)ebeamE, fileno)});

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
    double W = (photon + *nucleon).M();
    double t = -1.0 * (hyperon1->P4() - *nucleon).M2();

    double Mcap = W - 0.493677;
    if (Mcap < 3.2)
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
    hKp1Th.Fill(Kp1Theta);

    if (Kp1Theta >= 5 && Kp1Theta <= 45)
      Acceptance++;

    hKpAcceptance.Fill(Acceptance);

    auto LambdaRec = Proton->P4() + Pim->P4();
    hHyperon1Rec.Fill(LambdaRec.M());
    hHypW.Fill(LambdaRec.M(), W);

    hKp1Th_v_LambdaRec.Fill(LambdaRec.M(), Kp1Theta);

    hPKp.Fill(Kp->P4().P());
    hPLambda.Fill(hyperon1->P4().P());
    hPProton.Fill(Proton->P4().P());
    hPPim.Fill(Pim->P4().P());
  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");

  // internally stored histograms for total ep cross section
  TH1D *hWdist = (TH1D *)gDirectory->FindObject("Wdist");
  TH1D *hGenWdist = (TH1D *)gDirectory->FindObject("genWdist");

  TFile *fout = TFile::Open(Form("/home/nics/work/York/elSpectro/scripts/outputs/ed_to_KpLambda_%d_%d.root", (int)ebeamE, fileno), "recreate");
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
  hKp1Th_v_LambdaRec.Write();
  hPKp.Write();
  hPLambda.Write();
  hPProton.Write();
  hPPim.Write();
  hHypW.Write();
  hMcap.Write();
  hKp1Th.Write();
  hKpAcceptance.Write();
  fout->Close();

  generator().Summary();
}
