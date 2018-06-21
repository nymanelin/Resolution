
//#define NH 3
#define NC 6 // Number of Centrality bin

//https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/jonderwa/2014-Apr-29-analysis_note-analysis_note_event_plane_calibration.pdf
//https://www.hepdata.net/record/78365
//https://indico.cern.ch/event/703569/contributions/2886293/attachments/1597688/2531598/2018.02.08-Slupecki-ToyFlow.pdf
//https://aliceinfo.cern.ch/ArtSubmission/sites/aliceinfo.cern.ch.ArtSubmission/files/draft/cholm/2017-Jun-15-paper_draft-cds-rb-adraft-20170109-1.pdf

//#define ETADST_N 34
#define ETADST_N 36
static double etadst[ETADST_N] = {
  -3.8, //underflow
  -3.375,-3.125,-2.875,-2.625,-2.375,-2.125,-1.875,-1.625,
  -1.375,-1.125,-0.875,-0.625,-0.375,-0.125, 0.125, 0.375,
  0.625, 0.875, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375,
  2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.125, 4.375,
  4.625, 4.875,
  5.2, //overflow
};

static double etanch[NC+2][ETADST_N] = {
  {1640,1643,1670,1718,1787,1835,1912,1968,2001,2021,2017,1995,1970,1943,1929,1929,1943,1970,1995,2017,2021,2001,1968,1912,1835,1787,1718,1670,1643,1563,1474,1370,1324,1281,1244,1240},
  {1360,1364,1391,1424,1474,1507,1569,1644,1679,1682,1672,1646,1621,1597,1583,1583,1597,1621,1646,1672,1682,1679,1644,1569,1507,1474,1424,1391,1364,1292,1218,1132,1093,1062,1032,1030},
  {1000,1038,1061,1080,1114,1136,1178,1229,1253,1256,1247,1229,1210,1191,1181,1181,1191,1210,1229,1247,1256,1253,1229,1178,1136,1114,1080,1061,1038,977,921.3,857.7,829.6,807.4,787,780},
  {700,714,726,738,759,772,797,827,842,844,838,826,811.9,799.2,792.4,792.4,799.2,811.9,826,838,844,842,827,797,772,759,738,726,714,665,625.4,582.6,565.5,551.4,538,530},
  {460,475,482.7,489.7,502.6,510.6,522,539.9,549,549.3,545.5,537.5,527.6,519.3,514.7,514.7,519.3,527.6,537.5,545.5,549.3,549,539.9,522,510.6,502.6,489.7,482.7,475,440,413.6,386.7,375.6,368,359.9,350},
  {300,302,306.3,310.1,317.9,322.3,327.6,335.1,340,340.2,337.7,332.5,326.3,320.7,317.5,317.5,320.7,326.3,332.5,337.7,340.2,340,335.1,327.6,322.3,317.9,310.1,306.3,302,277.5,261.3,244.7,238.4,233.8,229.4,225},
  {177,178,179.9,181.7,186,188.2,189.8,193.5,196.4,196.5,194.8,191.4,187.5,184.3,182.5,182.5,184.3,187.5,191.4,194.8,196.5,196.4,193.5,189.8,188.2,186,181.7,179.9,178,163.2,153.4,143.8,140.3,138.7,136,135},
  {93,94.9,96.1,96.8,98.3,98.8,99.1,101.2,102.7,103.1,102,100.3,98,96.1,95.2,95.2,96.1,98,100.3,102,103.1,102.7,101.2,99.1,98.8,98.3,96.8,96.1,94.9,86.8,81.9,77.3,75.8,75.1,73.8,73}
};

static double CentBins[NC+3] = {0,5,10,20,30,40,50, 60, 70};

enum DETECTOR{
  D_TPC, //TPC full coverage
  D_TPC_ETAA, //TPC with eta gap
  D_TPC_ETAC,
  D_V0A,
  D_V0C,
  D_V0P, //V0+
  D_COUNT
};
static double cov[D_COUNT][2] = {
  {-1.5,1.5},
  {-1.5,-0.4},
  {0.4,1.5},
  {2.8,5.1},
  {-3.7,-1.7},
  {2.19,5.08}
};

enum RESOLUTION{
  R_V0A,
  R_V0C,
  R_V0P,
  R_COUNT
};
static const char *presn[] = {"V0A","V0C","V0P"};

//correlation histograms
TH1D *pah[R_COUNT][NC], *pbh[R_COUNT][NC], *pch[R_COUNT][NC];
TH1D *evph[R_COUNT][NC];
TH2D *contami2d[R_COUNT][NC];
TH2D *highcontami2d[R_COUNT][NC];
TH2D *tracks2d[R_COUNT][NC];
TH2D *evpcorr2d[R_COUNT][NC];
TH2D *evpcorrvsdet2d[NC];
TH1D *evpdifference[R_COUNT][NC];
TH2D *jetcorr2d[R_COUNT][NC];
TH1D *h_jetdirection[R_COUNT][NC][2];

TProfile *h_v2[R_COUNT];
TProfile *h_resolution[R_COUNT];
TH1D *p_v2[R_COUNT];
TH1D *p_reso[R_COUNT];
//TH2D *samecorr2dtrue[R_COUNT][NC][6];
TH2D *samecorr2d[R_COUNT][NC][6][6];

double pi = TMath::Pi();

int checkplane(double evp, double phi) {
  double diff = phi- evp;
  if (diff < 0) diff += 2*pi;

  //  if (diff < 1./6*pi && diff >= 0) return 0;
  //  else if (diff < 2./6*pi && diff >= 1./6*pi ) return 1;
  //  else if (diff < 4./6*pi && diff >= 2./6*pi ) return 2;
  //  else if (diff < 5./6*pi && diff >= 4./6*pi ) return 3;
  //  else if (diff < 7./6*pi && diff >= 5./6*pi ) return 4;
  //  else if (diff < 8./6*pi && diff >= 7./6*pi ) return 5;
  //  else if (diff < 10./6*pi && diff >= 8./6*pi ) return 6;
  //  else if (diff < 11./6*pi && diff >= 10./6*pi ) return 7;
  //  else if (diff < 12./6*pi && diff >= 11./6*pi ) return 0;

  //  if (diff < 1./8*pi && diff >= 0) return 0;
  //  else if (diff < 3./8*pi && diff >= 1./8*pi ) return 1;
  //  else if (diff < 5./8*pi && diff >= 3./8*pi ) return 2;
  //  else if (diff < 7./8*pi && diff >= 5./8*pi ) return 3;
  //  else if (diff < 9./8*pi && diff >= 7./8*pi ) return 4;
  //  else if (diff < 11./8*pi && diff >= 9./8*pi ) return 5;
  //  else if (diff < 13./8*pi && diff >= 11./8*pi ) return 6;
  //  else if (diff < 15./8*pi && diff >= 13./8*pi ) return 7;
  //  else if (diff < 16./8*pi && diff >= 15./8*pi ) return 0;

  //  if (diff < 1./12*pi && diff >= 0) return 0;
  //  else if (diff < 4./12*pi && diff >= 2./12*pi ) return 1;
  //  else if (diff < 7./12*pi && diff >= 5./12*pi ) return 2;
  //  else if (diff < 10./12*pi && diff >= 8./12*pi ) return 3;
  //  else if (diff < 13./12*pi && diff >= 11./12*pi ) return 4;
  //  else if (diff < 16./12*pi && diff >= 14./12*pi ) return 5;
  //  else if (diff < 19./12*pi && diff >= 17./12*pi ) return 6;
  //  else if (diff < 22./12*pi && diff >= 20./12*pi ) return 7;
  //  else if (diff < 24./12*pi && diff >= 23./12*pi ) return 0;

  if (diff < 1./12*pi && diff >= 0) return 0;
  else if (diff < 2./12*pi && diff >= 1./12*pi) return 1;
  else if (diff < 3./12*pi && diff >= 2./12*pi) return 2;
  else if (diff < 4./12*pi && diff >= 3./12*pi) return 3;
  else if (diff < 5./12*pi && diff >= 4./12*pi) return 4;
  else if (diff < 7./12*pi && diff >= 5./12*pi) return 5;
  else if (diff < 8./12*pi && diff >= 7./12*pi) return 4;
  else if (diff < 9./12*pi && diff >= 8./12*pi) return 3;
  else if (diff < 10./12*pi && diff >= 9./12*pi) return 2;
  else if (diff < 11./12*pi && diff >= 10./12*pi) return 1;
  else if (diff < 13./12*pi && diff >= 11./12*pi) return 0;
  else if (diff < 14./12*pi && diff >= 13./12*pi) return 1;
  else if (diff < 15./12*pi && diff >= 14./12*pi) return 2;
  else if (diff < 16./12*pi && diff >= 15./12*pi) return 3;
  else if (diff < 17./12*pi && diff >= 16./12*pi) return 4;
  else if (diff < 19./12*pi && diff >= 17./12*pi) return 5;
  else if (diff < 20./12*pi && diff >= 19./12*pi) return 4;
  else if (diff < 21./12*pi && diff >= 20./12*pi) return 3;
  else if (diff < 22./12*pi && diff >= 21./12*pi) return 2;
  else if (diff < 23./12*pi && diff >= 22./12*pi) return 1;
  else if (diff < 24./12*pi && diff >= 23./12*pi) return 0;

  return -9;
}

double CheckDetectorPhi(double phi) {

  double angle[8] = {-3*pi/4, -1*pi/2, -1*pi/4, 0, pi/4, pi/2, 3*pi/4, pi};
  double medianangle[8] = {-7*pi/8, -5*pi/8, -3*pi/8, -1*pi/8, pi/8, 3*pi/8, 5*pi/8, 7*pi/8};
  int i = 0;
  for ( ; i < 8; i++) {
    if (phi < angle[i]) break;
  }
  return medianangle[i];
}

