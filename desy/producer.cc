#include <iostream>

extern "C" {
  void get_seeds_(int&);
  void put_seeds_(int&);
  void zduini_();
  void zduevt_(int&);
  void zduend_();
}

int main() {
  int i;

  /*int IEND = 3;
  int NCVG = 1e5;
  int ITVG = 9;
  double INPP = 3500.;
  int PMOD = 11;
  double INPT = 3500.;
  int EMOD = 2;
  int PAIR = 13;
  int MCUT = 2;
  double PTCT = 15;*/

  i = 1;
  get_seeds_(i);
  
  zduini_();
  for (int n=0; n<10; n++) {
    i = 0;
    zduevt_(i);
  }
  zduend_();

  return 0;
}
