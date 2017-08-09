#define nd 10
extern "C" {
  void gmuini_();
  void zduini_();

  void gmucha_();

  void gmubeg_();

  void zduevt_(int* iwant);
  int luchge_(int&);
  extern struct {
    double inpe, inpp;
    int intge, intgp, gpdf, spdf, pmod, emod, ipair, nquark;
  } beam_;
  extern struct {
    int n, k[5][4000];
    float p[5][4000], v[5][4000];
  } lujets_;
  extern struct {
    double t1min,t1max;
    double t2min,t2max;
    double d3;
  } photons_;
  extern struct {
    double s1,s2,t1,t2;
  } extra_;
  extern struct {
    float s1,s2,s3,s4;
  } vgres_;
  extern struct {
    int ndim,ncvg,itmx,nprn,igraph,npoin,nprin,ntreat,ibeg,iend,ngen;
  } vegpar_;
  extern struct {
    double w,valtreat, x[nd], z[nd];
  } treatb_;
  extern struct {
    double u1,u2,v1,v2;
    double t11,t12,t21,t22;
  } peric_;
  extern struct {
    double tmx;
  } mykin_;
  extern struct {
    double cotth1,cotth2,ecut,ptcutmin,ptcutmax,mxmin2,mxmax2;
    float thmax,thmin;
    double qp2min,qp2max;
    int modcut;
    float mxmn,mxmx,q2mn,q2mx;
  } cuts_;
}
