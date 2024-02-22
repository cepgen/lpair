#ifndef cdf_lpair_h
#define cdf_lpair_h

#ifdef __cplusplus
extern "C" {
#endif

void fileini_();
void integrate_();
void generate_(int&);
void fragmentation_();
int luchge_(int&);

extern struct {
  std::array<int, 20> ipar;
  std::array<double, 20> lpar;
} datapar_;

extern struct {
  int accepted, ndim;
  std::array<double, 10> x;
} event_;

extern struct {
  int n;
  std::array<int, 5> k[4000];
  std::array<float, 5> p[4000];
} lujets_;

extern struct {
  double s1, s2, s3, s4;
} result_;


#ifdef __cplusplus
}
#endif

#endif
