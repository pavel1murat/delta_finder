///////////////////////////////////////////////////////////////////////////////
// stations
#include "TCanvas.h"
#include "TMarker.h"
#include "TGaxis.h"
//-----------------------------------------------------------------------------
void draw_panels(int Station = 0) {
  double rmax(80.), r(70.), rmin(38.);

  TCanvas* c = new TCanvas(Form("c_%i",Station),Form("station %i",Station),1000,1000);

  c->Divide(2,2);

  TEllipse* e1[4];

  double len  = sqrt(r*r-rmin*rmin);
    
  struct data_t {
    double phi;
    double phi_n;  // normal pointing outwards;
  };

  double z[2][4] = {
    -1518., -1490., -1467., -1539.,
    -1342., -1319., -1291., -1268.
  };
//-----------------------------------------------------------------------------
// ordered in Z
//-----------------------------------------------------------------------------
  data_t data[2][4][3] = {
    data_t{-75  , -75.- 90}, data_t{  45.,  45-90.},  data_t{ 165.,  165-90.},   // station:0 plane:0 face:1 z:-1518
    data_t{ 15. ,  15.+ 90}, data_t{ 135., 135+90.},  data_t{ 255.,  255+90.},   // station:0 plane:0 face:0 z:-1490
    data_t{-15. , -15 - 90}, data_t{-135.,-135-90.},  data_t{ 105.,  105-90.},   // station:0 plane:1 face:0 z:-1467
    data_t{ 75. ,  75.+ 90}, data_t{ -45., -45+90.},  data_t{-165., -165+90.},   // station:0 plane:1 face:1 z:-1539

    data_t{-15. , -15.- 90}, data_t{-135.,-135-90.},  data_t{ 105.,  105-90.},   // station:1 plane:2 face:0 z:-1342
    data_t{ 75. ,  75.+ 90}, data_t{ -45., -45+90.},  data_t{-165., -165+90.},   // station:1 plane:2 face:1 z:-1319
    data_t{-75. , -75.- 90}, data_t{ +45., +45-90.},  data_t{ 165.,  165-90.},   // station:1 plane:3 face:0 z:-1290
    data_t{ 15. ,  15.+ 90}, data_t{ 135., 135+90.},  data_t{-105., -105+90.},   // station:1 plane:3 face:1 z:-1268
  } ;

  TH2F* h2[4];
  
  for (int face=0; face<4; face++) {
    
    c->cd(face+1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    h2[face] = new TH2F(Form("h_face_%i",face),Form("zface:%i z:%8.2f",face,z[Station][face]),160,-80,80,160,-80,80);
    h2[face]->SetStats(0);
    h2[face]->Draw();

    e1[face]= new TEllipse(0,0,rmax-10,rmax-10,0.,2*M_PI*360,0);
    e1[face]->SetLineColor(kRed);
    e1[face]->Draw();

    for (int i=0; i<3; i++) { 
      double phi   = data[Station][face][i].phi*M_PI/180.;
      double phi_n = data[Station][face][i].phi_n*M_PI/180.;
      
      double x1    = rmin*cos(phi_n)+len*cos(phi);
      double y1    = rmin*sin(phi_n)+len*sin(phi);
      double x2    = rmin*cos(phi_n)-len*cos(phi);
      double y2    = rmin*sin(phi_n)-len*sin(phi);

      TArrow* a = new TArrow(x2,y2,x1,y1,0.015);
      a->Draw();

      double r1  = rmin+10;
      double xt1 = r1*cos(phi_n);
      double yt1 = r1*sin(phi_n);

      TText* t = new TText(xt1,yt1,Form("panel %i",i));
      t->SetTextAngle(data[Station][face][i].phi);
      t->SetTextFont(42);
      t->SetTextSize(0.035);
      t->Draw();
    }
  }
    
                              
  gPad->Modified();
  gPad->Update();
}
