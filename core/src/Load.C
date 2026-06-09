namespace elSpectro{};
using namespace elSpectro;

namespace elSpectro{namespace escat{}};
using namespace escat;

namespace jpacPhoto{};
using namespace jpacPhoto;



void Load(){

 
 
  TString JPAC = gSystem->Getenv("JPACPHOTO");
  
  TString ELSPECTRO = gSystem->Getenv("ELSPECTRO");
  if(ELSPECTRO.Length()==0){
    Fatal("elSpectro::Load","environment variable ELSPECTRO not set");
  }
  //if not defined jpacPhoto look for elSpectro submodule
  if(JPAC.Length()==0) JPAC = ELSPECTRO+"/jpacPhoto";

  // Make ROOT autoload aware of local build/install library locations.
  gSystem->AddDynamicPath(ELSPECTRO+"/build/core");
  gSystem->AddDynamicPath(ELSPECTRO+"/lib");
  gSystem->AddDynamicPath(ELSPECTRO+"/build/jpacPhoto/src/core");
  gSystem->AddDynamicPath(JPAC+"/lib");
  
  gInterpreter->AddIncludePath(JPAC+"/include/");
  // Prefer in-tree build output when developing, then fall back to installed libs.
  auto jlib=gSystem->Load(ELSPECTRO+"/build/jpacPhoto/src/core/libjpacPhoto."+gSystem->GetSoExt());
  if(jlib!=0) jlib=gSystem->Load(JPAC+"/lib/libjpacPhoto."+gSystem->GetSoExt());
  // If not, check LD_LIBRARY_PATH
  //if(jlib!=0) jlib=gSystem->Load(TString("libjpacPhoto.")+gSystem->GetSoExt());
  if(jlib!=0) Warning("elSpectro::Load","libjpacPhoto not found");
  
  gInterpreter->AddIncludePath(ELSPECTRO+"/core");
  // Prefer in-tree build output when developing, then fall back to installed libs.
  auto ellib=gSystem->Load(ELSPECTRO+"/build/core/libelSpectro."+gSystem->GetSoExt());
  if(ellib!=0) ellib=gSystem->Load(ELSPECTRO+"/lib/libelSpectro."+gSystem->GetSoExt());
  // If not, check LD_LIBRARY_PATH
  if(ellib!=0)ellib=gSystem->Load(TString("libelSpectro.")+gSystem->GetSoExt());
  if(ellib!=0) Fatal("elSpectro::Load","libelSpectro not found");

  //libs loaded can continue
  gROOT->ProcessLine("elSpectro::Manager::Instance();");

}
