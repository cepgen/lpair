#include <iostream>

using namespace std;

extern "C" {
  void zduini_();
  void zduevt_(int* iwant);
  extern struct {
    double inpe, inpp;
    int intge, intgp, gpdf, spdf, pmod, emod, ipair, nquark;
  } beam_;
  void printlhe_(int* mode);
  extern struct {
    int n, npad, k[5][4000];
    double p[5][4000], v[5][4000];
  } pyjets_;
}

main() {
  int one = 1;
  int two = 2;
  int three = 3;

  zduini_();
  printlhe_(&one);
  for (int i=0; i<1e2; i++) {
    zduevt_(&one);
    for (int p=0; p<pyjets_.n; p++) {
      cout << "--> " << p << " = " << pyjets_.p[0][p] << "\t" << pyjets_.p[1][p] << "\t" << pyjets_.p[2][p] << "\t" << pyjets_.p[3][p] << endl;
    }
    printlhe_(&two);
  }
  printlhe_(&three);
  cout << "-> " << beam_.pmod << endl;
  cout << "-> " << beam_.ipair << endl;
  cout << "-> " << beam_.inpe << endl;
  cout << "-> " << beam_.inpp << endl;
}
