///////////////////////////////////////////////////////////////////////////////
// DeltaFinder debugging tool
//
// takes hit printout from a .txt file and plots them
// assume we're passing data from one station, one combohit per face
//
// x = new TDisplayHits("a.txt")
// x->Draw()
// x->evaluateSeeds2()     // prints all intersections
///////////////////////////////////////////////////////////////////////////////
#include "TCanvas.h"
#include "TObjArray.h"
#include "TH2.h"
#include "TLine.h"

// #include "Offline/CalPatRec/inc/DeltaFinder_types.hh"
// #include "Offline/CalPatRec/inc/DeltaSeed.hh"

#include "Stntuple/geom/tracker_geom.hh"
#include "Stntuple/gui/TEvdComboHit.hh"

using stntuple::TEvdComboHit;
using mu2e::ComboHit;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
class TDisplayHits : public TObject {
public:
  TObjArray* fListOfComboHits;

  TObjArray* fListOfEvdHits;  // combo hits
  TObjArray* fListOfWires;    // wires

  TH2F*      fHist;

  TDisplayHits(const char* Fn = "");
  ~TDisplayHits();

  int nHits() { return fListOfComboHits->GetEntriesFast(); }

  ComboHit*     comboHit(int I) { return (ComboHit*)     fListOfComboHits->At(I); }
  TEvdComboHit* evdHit  (int I) { return (TEvdComboHit*) fListOfEvdHits->At(I); }

  void          calculateCogAndChi2(float RCore=3.0, float SigmaR2=3.0);
  int           evaluateSeeds2     ();
  int           evaluateSeed       (int I0, int I1, int I2, int I3);
  int           readData           (const char* Fn);
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  void Paint(Option_t* Opt);
};


//-----------------------------------------------------------------------------
TDisplayHits::TDisplayHits(const char* Fn) {
  fListOfComboHits  = new TObjArray();
  fListOfComboHits->SetOwner(kFALSE);
  
  fListOfEvdHits  = new TObjArray();
  fListOfEvdHits->SetOwner(kTRUE);

  fListOfWires  = new TObjArray();
  fListOfWires->SetOwner(kTRUE);

  readData(Fn);

  // fCanvas = new TCanvas("c","c",900,900);

  fHist = new TH2F("h2","h2",160,-800,800,160,-800,800);
  fHist->SetStats(0);
  fHist->Draw();

  gPad->SetGridx(1);
  gPad->SetGridy(1);

  Draw();

  gPad->Modified();
  gPad->Update();
}

//-----------------------------------------------------------------------------
TDisplayHits::~TDisplayHits() {
  delete fListOfEvdHits;
  
  int nh = fListOfComboHits->GetEntriesFast();

  for (int i=0; i<nh; i++) {
    ComboHit* ch = comboHit(i);
    delete ch;
  }

  delete fListOfComboHits;
  delete fListOfWires;
  delete fHist;
}

/*            hex
  145  2062 0x00000042  2   1:0:0: 2 0 0 14    1384.65  1378.40  -1.000   0.002758    375.410  160.611         11          852    4.121  252.539 -507.984 -1341.614    33     1     0
  211  2725 0x00000042  1   1:1:2: 2 5 1 37    1465.51  1435.82  -1.000   0.001654    -72.034   39.960         11          852    4.121  197.857 -460.093 -1321.786    34     1     0
  212  2726 0x00000042  2   1:1:2: 2 5 0 38    1390.55  1375.64  -1.000   0.003157    -37.913   28.493         11          854    0.123  166.112 -473.452 -1319.080    33     1     0
  257  3482 0x00000002  1   1:2:1: 3 3 0 26    1397.91  1382.63  -1.000   0.002606   -268.633   44.590         11          852    4.121  136.201 -516.105 -1293.626    -1     1     0
*/
//-----------------------------------------------------------------------------
// read hits as they are printed, w/o reformatting, use part of that to intialize
// TEvdComboHits
//-----------------------------------------------------------------------------
int TDisplayHits::readData(const char* Fn) {
  int rc(0);

  // ChannelID cx, co;

  FILE* f = fopen(Fn,"r");
  if (f == nullptr) {
    printf("ERROR cant open %s. BAIL OUT\n",Fn);
    return -1;
  }

  char   c[1000];

  char   c_flag[100], c_code[100];
  int    ih, sid, flag, nsh, station, zface, iplane, plane, panel, layer(0), straw; 
  float  t, tcorr, dt, edep, wdist, wres;
  int    pdg_id, sim_id, mom_pdg, mom_id;
  float  pstart, p, pz, x , y, z;
  int    delta_id, radOK, edepOK;

  int color[4] = { kRed, kGreen, kBlue, kMagenta };  // rgbm

  while ((c[0]=getc(f)) != EOF) {

					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);
					// parse line
      fscanf(f,"%i %i %s %i"      ,&ih,&sid,c_flag,&nsh );
      sscanf(c_flag,"0x%x"        ,&flag);
      
      fscanf(f,"%s",c_code);
      sscanf(c_code,"%i:%i:%2i"    ,&station,&zface,&plane);

      fscanf(f,"%i %i"            ,&panel,&straw  );
      fscanf(f,"%f %f %f %f %f %f",&t,&tcorr,&dt,&edep,&wdist,&wres);
      fscanf(f,"%i %i %i %i"      ,&pdg_id,&sim_id,&mom_pdg,&mom_id);
      fscanf(f,"%f %f %f %f %f %f",&pstart,&p,&pz,&x,&y,&z);
      fscanf(f,"%i"               ,&delta_id);
//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------
      printf("c_flag: %s\n",c_flag);
      printf("c_code: %s ",c_code);
      printf("station:zface:plane:  panel,straw: %02i:%i:%02i  %2i %2i\n",station,zface,plane,panel,straw);

      printf("%10.3f %10.3f %10.3f %10.5f %10.3f %10.3f\n",t,tcorr,dt,edep,wdist,wres);
      printf("pdg_id, sim_id: %5i %5i\n",pdg_id,sim_id);
      printf("%10.3f %10.3f %10.3f\n",x,y,z);
//-----------------------------------------------------------------------------
// create a new combo hit
//-----------------------------------------------------------------------------
      ComboHit* ch = new ComboHit();
      ch->_pos.SetXYZ(x,y,z);
      
      int iplane = plane   % 2;
      int is     = station % 2;
      // int ipanel = panel   / 2;

      stntuple::LayerGeom_t* lg = &stntuple::lgeom[is][iplane][panel][0];

      float wx = lg->wx;
      float wy = lg->wy;

      ch->_wdir.SetXYZ(wx,wy,0);

      ch->_sid    = mu2e::StrawId(sid);
      ch->_wres   = wres;
      ch->_wdist  = wdist;
      ch->_time   = t;
      ch->_nsh    = nsh;
      ch->_ncombo = nsh;
      ch->_edep   = edep;

      printf("-- iplane:panel,wx,wy,wres,wdist,t,tcorr,edep : %i:%i %8.5f %8.5f %8.3f %8.3f %8.3f %8.3f %8.5f\n",
             iplane,panel,wx,wy,wres,wdist,t,tcorr,edep);

      fListOfComboHits->Add((TObject*) ch);
//-----------------------------------------------------------------------------
// now create a TEvdComboHit
//-----------------------------------------------------------------------------
      TEvdComboHit* ech = new TEvdComboHit(ch,nullptr,nullptr,-1, p,p);
      fListOfEvdHits->Add(ech);
      ech->SetColor(color[zface]);
      ech->SetLineWidth(3);

      TLine* wire  = new TLine(x-wx*1000,y-wy*1000,x+wx*1000,y+wy*1000);
      fListOfWires->Add(wire);
    }
					// skip (the end of ?) the line
    fgets(c,1000,f);
  }

  fclose(f);

  printf("n(EVD/Combo hits) = %i\n",fListOfEvdHits->GetEntries());
  
  return rc;
}

//-----------------------------------------------------------------------------
// loops over all pairs of hits
// for each pair, calculates an intersection of the wires and prints its parameters
//-----------------------------------------------------------------------------
int TDisplayHits::evaluateSeeds2() {

  int nh = nHits();

  printf("i1 i2    x0      y0       x1      y1      nx1      ny1     x2      y2       nx2      ny2     wd1     wd2    wres1   wres2  chi21  chi22\n");
  printf("---------------------------------------------------------------------------------------------------------------------------------------\n");

  for (int ih1=0; ih1<nh-1; ih1++) {
    ComboHit* ch1 = comboHit(ih1);

    double x1    = ch1->pos().x();
    double y1    = ch1->pos().y();
    double nx1   = ch1->wdir().x();
    double ny1   = ch1->wdir().y();

    for (int ih2=ih1+1; ih2<nh; ih2++) {
      ComboHit* ch2 = comboHit(ih2);
//-----------------------------------------------------------------------------
// do not intersect parallel wires
//-----------------------------------------------------------------------------
      if ((ch1->strawId().panel() == ch2->strawId().panel()) and 
	  (ch1->strawId().plane() == ch2->strawId().plane())     )    continue;

      double x2    = ch2->pos().x();
      double y2    = ch2->pos().y();
      double nx2   = ch2->wdir().x();
      double ny2   = ch2->wdir().y();

      double n1n2  = nx1*nx2+ny1*ny2;
      double q12   = 1-n1n2*n1n2;

      double r12n1 = (x1-x2)*nx1+(y1-y2)*ny1;
      double r12n2 = (x1-x2)*nx2+(y1-y2)*ny2;

      double wd1   = -(r12n2*n1n2-r12n1)/q12;

      double x0 = x1-nx1*wd1;
      double y0 = y1-ny1*wd1;

      double wd2   = -(r12n2-n1n2*r12n1)/q12;
//-----------------------------------------------------------------------------
// require both hits to be close enough to the intersection point
//-----------------------------------------------------------------------------
      float wres1 = ch1->wireRes();
      float wres2 = ch2->wireRes();

      float hd1_chi2 = wd1*wd1/ch1->wireErr2();
      float hd2_chi2 = wd2*wd2/ch2->wireErr2();

      printf("%2i %2i",ih1, ih2);
      printf(" %7.2f %7.2f %7.2f %7.2f %8.5f %8.5f %7.2f %7.2f %8.5f %8.5f",x0, y0, x1, y1, nx1, ny1,x2, y2,nx2,ny2);
      printf(" %7.2f %7.2f %7.2f %7.2f %6.1f %6.1f",wd1, wd2, wres1, wres2, hd1_chi2,hd2_chi2);
      printf("\n");
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
// calculate chi2 of the seed constructed of more than two hits
//-----------------------------------------------------------------------------
int TDisplayHits::evaluateSeed(int I0, int I1, int I2, int I3) {
//-----------------------------------------------------------------------------
// intersect the two straws, we need coordinates of the intersection point and
// two distances from hits to the intersection point, 4 numbers in total
//-----------------------------------------------------------------------------
  int nh = nHits();

  int ind[4];

  ind[0] = I0;
  ind[1] = I1;
  ind[2] = I2;
  ind[3] = I3;

  double sn(0), snx2(0), snxy(0), sny2(0), snxr(0), snyr(0);

  for (int i=0; i<4; i++) {
    if (ind[i] < 0)                                 continue;

    ComboHit* ch = comboHit(ind[i]);

    double x    = ch->pos().x();
    double y    = ch->pos().y();
    double nx   = ch->wdir().x();
    double ny   = ch->wdir().y();
    
    double nr = x*ny-y*nx;

    sn    += 1;
    snx2  += nx*nx;
    snxy  += nx*ny;
    sny2  += ny*ny;
    snxr  += nx*nr;
    snyr  += ny*nr;
  }

  double nxym, nx2m, ny2m, nxrm, nyrm;

  nxym = snxy/sn;
  nx2m = snx2/sn;
  ny2m = sny2/sn;
  nxrm = snxr/sn;
  nyrm = snyr/sn;
  
  double d  = nx2m*ny2m-nxym*nxym;
  
  double xc = (nyrm*nx2m-nxrm*nxym)/d;
  double yc = (nyrm*nxym-nxrm*ny2m)/d;
//-----------------------------------------------------------------------------
// calculate seed chi2 - can this be optimized ?
//-----------------------------------------------------------------------------
  double sum_chi2_all  = 0;
  double sum_chi2_perp = 0;

  double  sigmaR2(30*30.);         // 30mm^2

  printf(" i    xc      yc       x       y      dx      dy    dr_par dr_perp dr2_perp    chi2_par chi2_perp\n");
  printf("-------------------------------------------------------------------------------------------------\n");
  for (int i=0; i<4; i++) {
    if (ind[i] < 0) continue;
    ComboHit* ch = comboHit(ind[i]);

    double dx = ch->pos().x()-xc;
    double dy = ch->pos().y()-yc;
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
//-----------------------------------------------------------------------------
    const XYZVectorF& wdir = ch->wdir();

    double dxy_dot_wdir    = dx*wdir.x()+dy*wdir.y();
      
    double dx_par          = dxy_dot_wdir*wdir.x();
    double dy_par          = dxy_dot_wdir*wdir.y();
      
    double dx_perp         = dx-dx_par;
    double dy_perp         = dy-dy_par;

    double dxy_perp        = dx*wdir.y()-dy*wdir.x();

    float  dxy2_perp       = dx_perp*dx_perp+dy_perp*dy_perp;
    
    float  chi2_par        = dxy_dot_wdir*dxy_dot_wdir/ch->wireErr2();
    float  chi2_perp       = dxy2_perp/sigmaR2;
    float  chi2            = chi2_par + chi2_perp;
    sum_chi2_all          += chi2;
    sum_chi2_perp         += chi2_perp;

    printf("%2i %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f",ind[i],xc,yc,ch->pos().x(),ch->pos().y(),dx,dy);
    printf(" %7.2f %7.2f %8.2f %7.2f    %7.2f\n",dxy_dot_wdir, dxy_perp, dxy2_perp, chi2_par, chi2_perp);
  }

  printf("sum_chi2_all, sum_chi2_perp : %7.2f  %7.2f\n",sum_chi2_all,sum_chi2_perp);

  return 0;
}

//-----------------------------------------------------------------------------
void TDisplayHits::Paint(Option_t* Opt) {
					// 2D histogram sets the scale

  //   fHist->Paint();


  fListOfWires->Paint();

  int nh = fListOfEvdHits->GetEntriesFast();

  for (int ih=0; ih<nh; ih++) {
    stntuple::TEvdComboHit* hit = evdHit(ih);
    hit->PaintXY();
  }
}


void TDisplayHits::calculateCogAndChi2(float RCore, float SigmaR) {
//-----------------------------------------------------------------------------
// update seed time and X and Y coordinates, accurate knowledge of Z is not very relevant
// if the seed has only two hits from initial intersection, there is no coordinates
// to redefine, chi2perp = 0 and chi2w is the sum of the two ...
//-----------------------------------------------------------------------------
  double fSnx2(0), fSny2(0), fSnxy(0), fSnxr(0), fSnyr(0), qn(0);

  int nh = nHits();
  printf(" i qn    x0      y0      nx    ny      nr    fSnx2     fSnxy    fSny2    fSnxr    fSnyr\n");
  printf("---------------------------------------------------------------------------------------\n");
  for (int i=0; i<nh; i++) {
    ComboHit* ch = comboHit(i);
      
    double x0 = ch->pos().x();
    double y0 = ch->pos().y();
    double nx = ch->wdir().x();
    double ny = ch->wdir().y();

    double nr = x0*ny-y0*nx;

    fSnx2    += nx*nx;
    fSnxy    += nx*ny;
    fSny2    += ny*ny;
    fSnxr    += nx*nr;
    fSnyr    += ny*nr;

    qn++;

    printf("%2i %.0f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
           i,qn,x0,y0,nx,ny,nr,fSnx2,fSnxy,fSny2,fSnxr,fSnyr);
  }

  double d  = fSnx2*fSny2-fSnxy*fSnxy;

  double xc = (fSnyr*fSnx2-fSnxr*fSnxy)/d;
  double yc = (fSnyr*fSnxy-fSnxr*fSny2)/d;

  printf("snx2,snxy,sny2,snxr,snyr,d,xc,yc : %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
         fSnx2,fSnxy,fSny2,fSnxr,fSnyr,d,xc,yc);
  // CofM.SetX(xc);
  // CofM.SetY(yc);
//-----------------------------------------------------------------------------
// calculate seed chi2 - can this be optimized ?
//-----------------------------------------------------------------------------
  float fChi2Par  = 0;
  float fChi2Perp = 0;
  printf("  i      x        y     wx      wy   dxy_dot_w dxy_dot_n   wres     drr  chi2_par  chi2_perp\n");
  printf("--------------------------------------------------------------------------------------------\n");
  for (int i=0; i<nh; i++) {
    ComboHit* ch = comboHit(i);

    double dx = ch->pos().x()-xc;
    double dy = ch->pos().y()-yc;
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
//-----------------------------------------------------------------------------
    const XYZVectorF& wdir = ch->wdir();

    float dxy_dot_w = dx*wdir.x()+dy*wdir.y();
    float dxy_dot_n = dx*wdir.y()-dy*wdir.x();

    float wres2     = ch->wireErr2();

    float sigR2     = SigmaR*SigmaR;
    float chi2_par  = (dxy_dot_w*dxy_dot_w)/(sigR2+wres2);
    float drr       = fmax(fabs(dxy_dot_n)-RCore,0);
    float chi2_perp = (drr*drr)/sigR2;
    fChi2Par       += chi2_par;
    fChi2Perp      += chi2_perp;

    printf(" %2i %8.3f %8.3f %6.3f %6.3f  %8.3f %8.3f  %8.3f %7.3f %8.3f %8.3f\n",
           i,ch->pos().x(),ch->pos().y(),wdir.x(),wdir.y(), 
           dxy_dot_w, dxy_dot_n, sqrt(wres2), drr, chi2_par,chi2_perp);
           
  }
  printf("qn: %3.0f chi2_par: %10.3f chi2_perp: %10.3f\n",
         qn,fChi2Par,fChi2Perp);
}
