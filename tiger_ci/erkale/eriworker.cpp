#include "eriworker.h"
#include "solidharmonics.h"
#include "mathf.h"

static void transform_i0(size_t Nj, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(1*Nj*Nk*Nl,0.0);
  for(size_t jj=0;jj<Nj;jj++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] +=  2.8209479177387814e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
      }
}

static void transform_i1(size_t Nj, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(3*Nj*Nk*Nl,0.0);
  for(size_t jj=0;jj<Nj;jj++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] +=  4.8860251190291992e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] +=  4.8860251190291992e-01 * (*input)[(( 1*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 1*Nj+jj)*Nk+kk)*Nl+ll] +=  4.8860251190291992e-01 * (*input)[(( 2*Nj+jj)*Nk+kk)*Nl+ll];
      }
}

static void transform_i2(size_t Nj, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(5*Nj*Nk*Nl,0.0);
  for(size_t jj=0;jj<Nj;jj++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] += -3.1539156525252005e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  5.4627421529603959e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] +=  1.0925484305920792e+00 * (*input)[(( 1*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] +=  1.0925484305920792e+00 * (*input)[(( 2*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] += -3.1539156525252005e-01 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] += -5.4627421529603959e-01 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 1*Nj+jj)*Nk+kk)*Nl+ll] +=  1.0925484305920792e+00 * (*input)[(( 4*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] +=  6.3078313050503998e-01 * (*input)[(( 5*Nj+jj)*Nk+kk)*Nl+ll];
      }
}

static void transform_i3(size_t Nj, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(7*Nj*Nk*Nl,0.0);
  for(size_t jj=0;jj<Nj;jj++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] += -4.5704579946446572e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] +=  5.9004358992664352e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] +=  1.7701307697799304e+00 * (*input)[(( 1*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] += -4.5704579946446572e-01 * (*input)[(( 1*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] += -1.1195289977703462e+00 * (*input)[(( 2*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] +=  1.4453057213202769e+00 * (*input)[(( 2*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] += -4.5704579946446572e-01 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] += -1.7701307697799304e+00 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 1*Nj+jj)*Nk+kk)*Nl+ll] +=  2.8906114426405538e+00 * (*input)[(( 4*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  1.8281831978578631e+00 * (*input)[(( 5*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] += -5.9004358992664352e-01 * (*input)[(( 6*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] += -4.5704579946446572e-01 * (*input)[(( 6*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] += -1.1195289977703462e+00 * (*input)[(( 7*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] += -1.4453057213202769e+00 * (*input)[(( 7*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] +=  1.8281831978578631e+00 * (*input)[(( 8*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] +=  7.4635266518023080e-01 * (*input)[(( 9*Nj+jj)*Nk+kk)*Nl+ll];
      }
}

static void transform_i4(size_t Nj, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(9*Nj*Nk*Nl,0.0);
  for(size_t jj=0;jj<Nj;jj++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  3.1735664074561293e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] += -4.7308734787877998e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] +=  6.2583573544917614e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] +=  2.5033429417967046e+00 * (*input)[(( 1*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] += -9.4617469575755997e-01 * (*input)[(( 1*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] += -2.0071396306718676e+00 * (*input)[(( 2*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 7*Nj+jj)*Nk+kk)*Nl+ll] +=  1.7701307697799307e+00 * (*input)[(( 2*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  6.3471328149122586e-01 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] += -3.7550144126950569e+00 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 1*Nj+jj)*Nk+kk)*Nl+ll] +=  5.3103923093397922e+00 * (*input)[(( 4*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] += -2.0071396306718676e+00 * (*input)[(( 4*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] += -2.5388531259649034e+00 * (*input)[(( 5*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] +=  2.8385240872726798e+00 * (*input)[(( 5*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] += -2.5033429417967046e+00 * (*input)[(( 6*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] += -9.4617469575755997e-01 * (*input)[(( 6*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] += -2.0071396306718676e+00 * (*input)[(( 7*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 7*Nj+jj)*Nk+kk)*Nl+ll] += -5.3103923093397922e+00 * (*input)[(( 7*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] +=  5.6770481745453596e+00 * (*input)[(( 8*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] +=  2.6761861742291573e+00 * (*input)[(( 9*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  3.1735664074561293e-01 * (*input)[((10*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] +=  4.7308734787877998e-01 * (*input)[((10*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] +=  6.2583573544917614e-01 * (*input)[((10*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 1*Nj+jj)*Nk+kk)*Nl+ll] += -1.7701307697799307e+00 * (*input)[((11*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] += -2.0071396306718676e+00 * (*input)[((11*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] += -2.5388531259649034e+00 * (*input)[((12*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] += -2.8385240872726798e+00 * (*input)[((12*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] +=  2.6761861742291573e+00 * (*input)[((13*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  8.4628437532163425e-01 * (*input)[((14*Nj+jj)*Nk+kk)*Nl+ll];
      }
}

static void transform_i5(size_t Nj, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(11*Nj*Nk*Nl,0.0);
  for(size_t jj=0;jj<Nj;jj++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] +=  4.5294665119569688e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] += -4.8923829943525038e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((10*Nj+jj)*Nk+kk)*Nl+ll] +=  6.5638205684017015e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] +=  3.2819102842008507e+00 * (*input)[(( 1*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] += -1.4677148983057511e+00 * (*input)[(( 1*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  4.5294665119569688e-01 * (*input)[(( 1*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] +=  1.7542548368013540e+00 * (*input)[(( 2*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 7*Nj+jj)*Nk+kk)*Nl+ll] += -2.3967683924866621e+00 * (*input)[(( 2*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 9*Nj+jj)*Nk+kk)*Nl+ll] +=  2.0756623148810416e+00 * (*input)[(( 2*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] +=  9.0589330239139376e-01 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] +=  9.7847659887050065e-01 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((10*Nj+jj)*Nk+kk)*Nl+ll] += -6.5638205684017015e+00 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 1*Nj+jj)*Nk+kk)*Nl+ll] +=  8.3026492595241663e+00 * (*input)[(( 4*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] += -4.7935367849733241e+00 * (*input)[(( 4*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] += -5.4353598143483630e+00 * (*input)[(( 5*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] +=  3.9139063954820035e+00 * (*input)[(( 5*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] += -6.5638205684017015e+00 * (*input)[(( 6*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] += -9.7847659887050065e-01 * (*input)[(( 6*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  9.0589330239139376e-01 * (*input)[(( 6*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] +=  3.5085096736027079e+00 * (*input)[(( 7*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 9*Nj+jj)*Nk+kk)*Nl+ll] += -1.2453973889286249e+01 * (*input)[(( 7*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] +=  1.1741719186446010e+01 * (*input)[(( 8*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] += -5.4353598143483630e+00 * (*input)[(( 8*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] += -4.6780128981369433e+00 * (*input)[(( 9*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 7*Nj+jj)*Nk+kk)*Nl+ll] +=  4.7935367849733250e+00 * (*input)[(( 9*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] +=  4.5294665119569688e-01 * (*input)[((10*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] +=  1.4677148983057511e+00 * (*input)[((10*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((10*Nj+jj)*Nk+kk)*Nl+ll] +=  3.2819102842008507e+00 * (*input)[((10*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 1*Nj+jj)*Nk+kk)*Nl+ll] += -8.3026492595241663e+00 * (*input)[((11*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] += -4.7935367849733241e+00 * (*input)[((11*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] += -5.4353598143483630e+00 * (*input)[((12*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] += -1.1741719186446010e+01 * (*input)[((12*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] +=  9.5870735699466501e+00 * (*input)[((13*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] +=  3.6235732095655746e+00 * (*input)[((14*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] +=  6.5638205684017015e-01 * (*input)[((15*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] +=  4.8923829943525038e-01 * (*input)[((15*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  4.5294665119569688e-01 * (*input)[((15*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] +=  1.7542548368013540e+00 * (*input)[((16*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 7*Nj+jj)*Nk+kk)*Nl+ll] +=  2.3967683924866621e+00 * (*input)[((16*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 9*Nj+jj)*Nk+kk)*Nl+ll] +=  2.0756623148810416e+00 * (*input)[((16*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] += -3.9139063954820035e+00 * (*input)[((17*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] += -5.4353598143483630e+00 * (*input)[((17*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] += -4.6780128981369433e+00 * (*input)[((18*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 7*Nj+jj)*Nk+kk)*Nl+ll] += -4.7935367849733250e+00 * (*input)[((18*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  3.6235732095655746e+00 * (*input)[((19*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] +=  9.3560257962738991e-01 * (*input)[((20*Nj+jj)*Nk+kk)*Nl+ll];
      }
}

static void transform_i6(size_t Nj, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(13*Nj*Nk*Nl,0.0);
  for(size_t jj=0;jj<Nj;jj++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] += -3.1784601133814211e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] +=  4.6060262975746175e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((10*Nj+jj)*Nk+kk)*Nl+ll] += -5.0456490072872417e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((12*Nj+jj)*Nk+kk)*Nl+ll] +=  6.8318410519191419e-01 * (*input)[(( 0*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] +=  4.0991046311514854e+00 * (*input)[(( 1*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] += -2.0182596029148967e+00 * (*input)[(( 1*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  9.2120525951492349e-01 * (*input)[(( 1*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 7*Nj+jj)*Nk+kk)*Nl+ll] +=  2.9131068125936568e+00 * (*input)[(( 2*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 9*Nj+jj)*Nk+kk)*Nl+ll] += -2.7636157785447706e+00 * (*input)[(( 2*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((11*Nj+jj)*Nk+kk)*Nl+ll] +=  2.3666191622317521e+00 * (*input)[(( 2*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] += -9.5353803401442638e-01 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] +=  4.6060262975746175e-01 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((10*Nj+jj)*Nk+kk)*Nl+ll] +=  2.5228245036436205e+00 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((12*Nj+jj)*Nk+kk)*Nl+ll] += -1.0247761577878713e+01 * (*input)[(( 3*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 1*Nj+jj)*Nk+kk)*Nl+ll] +=  1.1833095811158760e+01 * (*input)[(( 4*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] += -8.2908473356343109e+00 * (*input)[(( 4*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] +=  2.9131068125936568e+00 * (*input)[(( 4*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] +=  5.7212282040865574e+00 * (*input)[(( 5*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] += -7.3696420761193870e+00 * (*input)[(( 5*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((10*Nj+jj)*Nk+kk)*Nl+ll] +=  5.0456490072872420e+00 * (*input)[(( 5*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] += -1.3663682103838283e+01 * (*input)[(( 6*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  1.8424105190298470e+00 * (*input)[(( 6*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 7*Nj+jj)*Nk+kk)*Nl+ll] +=  5.8262136251873136e+00 * (*input)[(( 7*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 9*Nj+jj)*Nk+kk)*Nl+ll] +=  5.5272315570895403e+00 * (*input)[(( 7*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((11*Nj+jj)*Nk+kk)*Nl+ll] += -2.3666191622317520e+01 * (*input)[(( 7*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] +=  2.0182596029148968e+01 * (*input)[(( 8*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] += -1.4739284152238774e+01 * (*input)[(( 8*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 7*Nj+jj)*Nk+kk)*Nl+ll] += -1.1652427250374625e+01 * (*input)[(( 9*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 9*Nj+jj)*Nk+kk)*Nl+ll] +=  7.3696420761193879e+00 * (*input)[(( 9*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] += -9.5353803401442638e-01 * (*input)[((10*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] += -4.6060262975746175e-01 * (*input)[((10*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((10*Nj+jj)*Nk+kk)*Nl+ll] +=  2.5228245036436205e+00 * (*input)[((10*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((12*Nj+jj)*Nk+kk)*Nl+ll] +=  1.0247761577878713e+01 * (*input)[((10*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 1*Nj+jj)*Nk+kk)*Nl+ll] += -2.3666191622317520e+01 * (*input)[((11*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] += -5.5272315570895403e+00 * (*input)[((11*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] +=  5.8262136251873136e+00 * (*input)[((11*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] +=  1.1442456408173115e+01 * (*input)[((12*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((10*Nj+jj)*Nk+kk)*Nl+ll] += -3.0273894043723452e+01 * (*input)[((12*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] +=  2.2108926228358165e+01 * (*input)[((13*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] += -1.1652427250374625e+01 * (*input)[((13*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] += -7.6283042721154128e+00 * (*input)[((14*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] +=  7.3696420761193888e+00 * (*input)[((14*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 0*Nj+jj)*Nk+kk)*Nl+ll] +=  4.0991046311514854e+00 * (*input)[((15*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] +=  2.0182596029148967e+00 * (*input)[((15*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  9.2120525951492349e-01 * (*input)[((15*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 7*Nj+jj)*Nk+kk)*Nl+ll] +=  2.9131068125936568e+00 * (*input)[((16*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 9*Nj+jj)*Nk+kk)*Nl+ll] +=  8.2908473356343109e+00 * (*input)[((16*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((11*Nj+jj)*Nk+kk)*Nl+ll] +=  1.1833095811158760e+01 * (*input)[((16*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 2*Nj+jj)*Nk+kk)*Nl+ll] += -2.0182596029148968e+01 * (*input)[((17*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] += -1.4739284152238774e+01 * (*input)[((17*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 7*Nj+jj)*Nk+kk)*Nl+ll] += -1.1652427250374625e+01 * (*input)[((18*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 9*Nj+jj)*Nk+kk)*Nl+ll] += -2.2108926228358165e+01 * (*input)[((18*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 4*Nj+jj)*Nk+kk)*Nl+ll] +=  1.4739284152238778e+01 * (*input)[((19*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 7*Nj+jj)*Nk+kk)*Nl+ll] +=  4.6609709001498523e+00 * (*input)[((20*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] += -3.1784601133814211e-01 * (*input)[((21*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] += -4.6060262975746175e-01 * (*input)[((21*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((10*Nj+jj)*Nk+kk)*Nl+ll] += -5.0456490072872417e-01 * (*input)[((21*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((12*Nj+jj)*Nk+kk)*Nl+ll] += -6.8318410519191419e-01 * (*input)[((21*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 1*Nj+jj)*Nk+kk)*Nl+ll] +=  2.3666191622317521e+00 * (*input)[((22*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] +=  2.7636157785447706e+00 * (*input)[((22*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] +=  2.9131068125936568e+00 * (*input)[((22*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] +=  5.7212282040865574e+00 * (*input)[((23*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] +=  7.3696420761193870e+00 * (*input)[((23*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[((10*Nj+jj)*Nk+kk)*Nl+ll] +=  5.0456490072872420e+00 * (*input)[((23*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 3*Nj+jj)*Nk+kk)*Nl+ll] += -7.3696420761193879e+00 * (*input)[((24*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] += -1.1652427250374625e+01 * (*input)[((24*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] += -7.6283042721154128e+00 * (*input)[((25*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 8*Nj+jj)*Nk+kk)*Nl+ll] += -7.3696420761193888e+00 * (*input)[((25*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 5*Nj+jj)*Nk+kk)*Nl+ll] +=  4.6609709001498523e+00 * (*input)[((26*Nj+jj)*Nk+kk)*Nl+ll];
        (*output)[(( 6*Nj+jj)*Nk+kk)*Nl+ll] +=  1.0171072362820530e+00 * (*input)[((27*Nj+jj)*Nk+kk)*Nl+ll];
      }
}

void IntegralWorker::transform_i(int am, size_t Nj, size_t Nk, size_t Nl) {
  static void (*f[7])(size_t, size_t, size_t, const std::vector<double> *, std::vector<double> *)={
    transform_i0,
    transform_i1,
    transform_i2,
    transform_i3,
    transform_i4,
    transform_i5,
    transform_i6
  };
    f[am](Nj,Nk,Nl,input,output);
  std::swap(input,output);
}
static void transform_j0(size_t Ni, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(1*Ni*Nk*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii* 1+ 0)*Nk+kk)*Nl+ll] +=  2.8209479177387814e-01 * (*input)[((ii* 1+ 0)*Nk+kk)*Nl+ll];
      }
}

static void transform_j1(size_t Ni, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(3*Ni*Nk*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii* 3+ 2)*Nk+kk)*Nl+ll] +=  4.8860251190291992e-01 * (*input)[((ii* 3+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii* 3+ 0)*Nk+kk)*Nl+ll] +=  4.8860251190291992e-01 * (*input)[((ii* 3+ 1)*Nk+kk)*Nl+ll];
        (*output)[((ii* 3+ 1)*Nk+kk)*Nl+ll] +=  4.8860251190291992e-01 * (*input)[((ii* 3+ 2)*Nk+kk)*Nl+ll];
      }
}

static void transform_j2(size_t Ni, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(5*Ni*Nk*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii* 5+ 2)*Nk+kk)*Nl+ll] += -3.1539156525252005e-01 * (*input)[((ii* 6+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii* 5+ 4)*Nk+kk)*Nl+ll] +=  5.4627421529603959e-01 * (*input)[((ii* 6+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii* 5+ 0)*Nk+kk)*Nl+ll] +=  1.0925484305920792e+00 * (*input)[((ii* 6+ 1)*Nk+kk)*Nl+ll];
        (*output)[((ii* 5+ 3)*Nk+kk)*Nl+ll] +=  1.0925484305920792e+00 * (*input)[((ii* 6+ 2)*Nk+kk)*Nl+ll];
        (*output)[((ii* 5+ 2)*Nk+kk)*Nl+ll] += -3.1539156525252005e-01 * (*input)[((ii* 6+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii* 5+ 4)*Nk+kk)*Nl+ll] += -5.4627421529603959e-01 * (*input)[((ii* 6+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii* 5+ 1)*Nk+kk)*Nl+ll] +=  1.0925484305920792e+00 * (*input)[((ii* 6+ 4)*Nk+kk)*Nl+ll];
        (*output)[((ii* 5+ 2)*Nk+kk)*Nl+ll] +=  6.3078313050503998e-01 * (*input)[((ii* 6+ 5)*Nk+kk)*Nl+ll];
      }
}

static void transform_j3(size_t Ni, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(7*Ni*Nk*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii* 7+ 4)*Nk+kk)*Nl+ll] += -4.5704579946446572e-01 * (*input)[((ii*10+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 6)*Nk+kk)*Nl+ll] +=  5.9004358992664352e-01 * (*input)[((ii*10+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 0)*Nk+kk)*Nl+ll] +=  1.7701307697799304e+00 * (*input)[((ii*10+ 1)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 2)*Nk+kk)*Nl+ll] += -4.5704579946446572e-01 * (*input)[((ii*10+ 1)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 3)*Nk+kk)*Nl+ll] += -1.1195289977703462e+00 * (*input)[((ii*10+ 2)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 5)*Nk+kk)*Nl+ll] +=  1.4453057213202769e+00 * (*input)[((ii*10+ 2)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 4)*Nk+kk)*Nl+ll] += -4.5704579946446572e-01 * (*input)[((ii*10+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 6)*Nk+kk)*Nl+ll] += -1.7701307697799304e+00 * (*input)[((ii*10+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 1)*Nk+kk)*Nl+ll] +=  2.8906114426405538e+00 * (*input)[((ii*10+ 4)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 4)*Nk+kk)*Nl+ll] +=  1.8281831978578631e+00 * (*input)[((ii*10+ 5)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 0)*Nk+kk)*Nl+ll] += -5.9004358992664352e-01 * (*input)[((ii*10+ 6)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 2)*Nk+kk)*Nl+ll] += -4.5704579946446572e-01 * (*input)[((ii*10+ 6)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 3)*Nk+kk)*Nl+ll] += -1.1195289977703462e+00 * (*input)[((ii*10+ 7)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 5)*Nk+kk)*Nl+ll] += -1.4453057213202769e+00 * (*input)[((ii*10+ 7)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 2)*Nk+kk)*Nl+ll] +=  1.8281831978578631e+00 * (*input)[((ii*10+ 8)*Nk+kk)*Nl+ll];
        (*output)[((ii* 7+ 3)*Nk+kk)*Nl+ll] +=  7.4635266518023080e-01 * (*input)[((ii*10+ 9)*Nk+kk)*Nl+ll];
      }
}

static void transform_j4(size_t Ni, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(9*Ni*Nk*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii* 9+ 4)*Nk+kk)*Nl+ll] +=  3.1735664074561293e-01 * (*input)[((ii*15+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 6)*Nk+kk)*Nl+ll] += -4.7308734787877998e-01 * (*input)[((ii*15+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 8)*Nk+kk)*Nl+ll] +=  6.2583573544917614e-01 * (*input)[((ii*15+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 0)*Nk+kk)*Nl+ll] +=  2.5033429417967046e+00 * (*input)[((ii*15+ 1)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 2)*Nk+kk)*Nl+ll] += -9.4617469575755997e-01 * (*input)[((ii*15+ 1)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 5)*Nk+kk)*Nl+ll] += -2.0071396306718676e+00 * (*input)[((ii*15+ 2)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 7)*Nk+kk)*Nl+ll] +=  1.7701307697799307e+00 * (*input)[((ii*15+ 2)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 4)*Nk+kk)*Nl+ll] +=  6.3471328149122586e-01 * (*input)[((ii*15+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 8)*Nk+kk)*Nl+ll] += -3.7550144126950569e+00 * (*input)[((ii*15+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 1)*Nk+kk)*Nl+ll] +=  5.3103923093397922e+00 * (*input)[((ii*15+ 4)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 3)*Nk+kk)*Nl+ll] += -2.0071396306718676e+00 * (*input)[((ii*15+ 4)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 4)*Nk+kk)*Nl+ll] += -2.5388531259649034e+00 * (*input)[((ii*15+ 5)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 6)*Nk+kk)*Nl+ll] +=  2.8385240872726798e+00 * (*input)[((ii*15+ 5)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 0)*Nk+kk)*Nl+ll] += -2.5033429417967046e+00 * (*input)[((ii*15+ 6)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 2)*Nk+kk)*Nl+ll] += -9.4617469575755997e-01 * (*input)[((ii*15+ 6)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 5)*Nk+kk)*Nl+ll] += -2.0071396306718676e+00 * (*input)[((ii*15+ 7)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 7)*Nk+kk)*Nl+ll] += -5.3103923093397922e+00 * (*input)[((ii*15+ 7)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 2)*Nk+kk)*Nl+ll] +=  5.6770481745453596e+00 * (*input)[((ii*15+ 8)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 5)*Nk+kk)*Nl+ll] +=  2.6761861742291573e+00 * (*input)[((ii*15+ 9)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 4)*Nk+kk)*Nl+ll] +=  3.1735664074561293e-01 * (*input)[((ii*15+10)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 6)*Nk+kk)*Nl+ll] +=  4.7308734787877998e-01 * (*input)[((ii*15+10)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 8)*Nk+kk)*Nl+ll] +=  6.2583573544917614e-01 * (*input)[((ii*15+10)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 1)*Nk+kk)*Nl+ll] += -1.7701307697799307e+00 * (*input)[((ii*15+11)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 3)*Nk+kk)*Nl+ll] += -2.0071396306718676e+00 * (*input)[((ii*15+11)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 4)*Nk+kk)*Nl+ll] += -2.5388531259649034e+00 * (*input)[((ii*15+12)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 6)*Nk+kk)*Nl+ll] += -2.8385240872726798e+00 * (*input)[((ii*15+12)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 3)*Nk+kk)*Nl+ll] +=  2.6761861742291573e+00 * (*input)[((ii*15+13)*Nk+kk)*Nl+ll];
        (*output)[((ii* 9+ 4)*Nk+kk)*Nl+ll] +=  8.4628437532163425e-01 * (*input)[((ii*15+14)*Nk+kk)*Nl+ll];
      }
}

static void transform_j5(size_t Ni, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(11*Ni*Nk*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii*11+ 6)*Nk+kk)*Nl+ll] +=  4.5294665119569688e-01 * (*input)[((ii*21+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 8)*Nk+kk)*Nl+ll] += -4.8923829943525038e-01 * (*input)[((ii*21+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+10)*Nk+kk)*Nl+ll] +=  6.5638205684017015e-01 * (*input)[((ii*21+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 0)*Nk+kk)*Nl+ll] +=  3.2819102842008507e+00 * (*input)[((ii*21+ 1)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 2)*Nk+kk)*Nl+ll] += -1.4677148983057511e+00 * (*input)[((ii*21+ 1)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 4)*Nk+kk)*Nl+ll] +=  4.5294665119569688e-01 * (*input)[((ii*21+ 1)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 5)*Nk+kk)*Nl+ll] +=  1.7542548368013540e+00 * (*input)[((ii*21+ 2)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 7)*Nk+kk)*Nl+ll] += -2.3967683924866621e+00 * (*input)[((ii*21+ 2)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 9)*Nk+kk)*Nl+ll] +=  2.0756623148810416e+00 * (*input)[((ii*21+ 2)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 6)*Nk+kk)*Nl+ll] +=  9.0589330239139376e-01 * (*input)[((ii*21+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 8)*Nk+kk)*Nl+ll] +=  9.7847659887050065e-01 * (*input)[((ii*21+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+10)*Nk+kk)*Nl+ll] += -6.5638205684017015e+00 * (*input)[((ii*21+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 1)*Nk+kk)*Nl+ll] +=  8.3026492595241663e+00 * (*input)[((ii*21+ 4)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 3)*Nk+kk)*Nl+ll] += -4.7935367849733241e+00 * (*input)[((ii*21+ 4)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 6)*Nk+kk)*Nl+ll] += -5.4353598143483630e+00 * (*input)[((ii*21+ 5)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 8)*Nk+kk)*Nl+ll] +=  3.9139063954820035e+00 * (*input)[((ii*21+ 5)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 0)*Nk+kk)*Nl+ll] += -6.5638205684017015e+00 * (*input)[((ii*21+ 6)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 2)*Nk+kk)*Nl+ll] += -9.7847659887050065e-01 * (*input)[((ii*21+ 6)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 4)*Nk+kk)*Nl+ll] +=  9.0589330239139376e-01 * (*input)[((ii*21+ 6)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 5)*Nk+kk)*Nl+ll] +=  3.5085096736027079e+00 * (*input)[((ii*21+ 7)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 9)*Nk+kk)*Nl+ll] += -1.2453973889286249e+01 * (*input)[((ii*21+ 7)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 2)*Nk+kk)*Nl+ll] +=  1.1741719186446010e+01 * (*input)[((ii*21+ 8)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 4)*Nk+kk)*Nl+ll] += -5.4353598143483630e+00 * (*input)[((ii*21+ 8)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 5)*Nk+kk)*Nl+ll] += -4.6780128981369433e+00 * (*input)[((ii*21+ 9)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 7)*Nk+kk)*Nl+ll] +=  4.7935367849733250e+00 * (*input)[((ii*21+ 9)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 6)*Nk+kk)*Nl+ll] +=  4.5294665119569688e-01 * (*input)[((ii*21+10)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 8)*Nk+kk)*Nl+ll] +=  1.4677148983057511e+00 * (*input)[((ii*21+10)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+10)*Nk+kk)*Nl+ll] +=  3.2819102842008507e+00 * (*input)[((ii*21+10)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 1)*Nk+kk)*Nl+ll] += -8.3026492595241663e+00 * (*input)[((ii*21+11)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 3)*Nk+kk)*Nl+ll] += -4.7935367849733241e+00 * (*input)[((ii*21+11)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 6)*Nk+kk)*Nl+ll] += -5.4353598143483630e+00 * (*input)[((ii*21+12)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 8)*Nk+kk)*Nl+ll] += -1.1741719186446010e+01 * (*input)[((ii*21+12)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 3)*Nk+kk)*Nl+ll] +=  9.5870735699466501e+00 * (*input)[((ii*21+13)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 6)*Nk+kk)*Nl+ll] +=  3.6235732095655746e+00 * (*input)[((ii*21+14)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 0)*Nk+kk)*Nl+ll] +=  6.5638205684017015e-01 * (*input)[((ii*21+15)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 2)*Nk+kk)*Nl+ll] +=  4.8923829943525038e-01 * (*input)[((ii*21+15)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 4)*Nk+kk)*Nl+ll] +=  4.5294665119569688e-01 * (*input)[((ii*21+15)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 5)*Nk+kk)*Nl+ll] +=  1.7542548368013540e+00 * (*input)[((ii*21+16)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 7)*Nk+kk)*Nl+ll] +=  2.3967683924866621e+00 * (*input)[((ii*21+16)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 9)*Nk+kk)*Nl+ll] +=  2.0756623148810416e+00 * (*input)[((ii*21+16)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 2)*Nk+kk)*Nl+ll] += -3.9139063954820035e+00 * (*input)[((ii*21+17)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 4)*Nk+kk)*Nl+ll] += -5.4353598143483630e+00 * (*input)[((ii*21+17)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 5)*Nk+kk)*Nl+ll] += -4.6780128981369433e+00 * (*input)[((ii*21+18)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 7)*Nk+kk)*Nl+ll] += -4.7935367849733250e+00 * (*input)[((ii*21+18)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 4)*Nk+kk)*Nl+ll] +=  3.6235732095655746e+00 * (*input)[((ii*21+19)*Nk+kk)*Nl+ll];
        (*output)[((ii*11+ 5)*Nk+kk)*Nl+ll] +=  9.3560257962738991e-01 * (*input)[((ii*21+20)*Nk+kk)*Nl+ll];
      }
}

static void transform_j6(size_t Ni, size_t Nk, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(13*Ni*Nk*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t kk=0;kk<Nk;kk++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii*13+ 6)*Nk+kk)*Nl+ll] += -3.1784601133814211e-01 * (*input)[((ii*28+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 8)*Nk+kk)*Nl+ll] +=  4.6060262975746175e-01 * (*input)[((ii*28+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+10)*Nk+kk)*Nl+ll] += -5.0456490072872417e-01 * (*input)[((ii*28+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+12)*Nk+kk)*Nl+ll] +=  6.8318410519191419e-01 * (*input)[((ii*28+ 0)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 0)*Nk+kk)*Nl+ll] +=  4.0991046311514854e+00 * (*input)[((ii*28+ 1)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 2)*Nk+kk)*Nl+ll] += -2.0182596029148967e+00 * (*input)[((ii*28+ 1)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 4)*Nk+kk)*Nl+ll] +=  9.2120525951492349e-01 * (*input)[((ii*28+ 1)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 7)*Nk+kk)*Nl+ll] +=  2.9131068125936568e+00 * (*input)[((ii*28+ 2)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 9)*Nk+kk)*Nl+ll] += -2.7636157785447706e+00 * (*input)[((ii*28+ 2)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+11)*Nk+kk)*Nl+ll] +=  2.3666191622317521e+00 * (*input)[((ii*28+ 2)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 6)*Nk+kk)*Nl+ll] += -9.5353803401442638e-01 * (*input)[((ii*28+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 8)*Nk+kk)*Nl+ll] +=  4.6060262975746175e-01 * (*input)[((ii*28+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+10)*Nk+kk)*Nl+ll] +=  2.5228245036436205e+00 * (*input)[((ii*28+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+12)*Nk+kk)*Nl+ll] += -1.0247761577878713e+01 * (*input)[((ii*28+ 3)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 1)*Nk+kk)*Nl+ll] +=  1.1833095811158760e+01 * (*input)[((ii*28+ 4)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 3)*Nk+kk)*Nl+ll] += -8.2908473356343109e+00 * (*input)[((ii*28+ 4)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 5)*Nk+kk)*Nl+ll] +=  2.9131068125936568e+00 * (*input)[((ii*28+ 4)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 6)*Nk+kk)*Nl+ll] +=  5.7212282040865574e+00 * (*input)[((ii*28+ 5)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 8)*Nk+kk)*Nl+ll] += -7.3696420761193870e+00 * (*input)[((ii*28+ 5)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+10)*Nk+kk)*Nl+ll] +=  5.0456490072872420e+00 * (*input)[((ii*28+ 5)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 0)*Nk+kk)*Nl+ll] += -1.3663682103838283e+01 * (*input)[((ii*28+ 6)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 4)*Nk+kk)*Nl+ll] +=  1.8424105190298470e+00 * (*input)[((ii*28+ 6)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 7)*Nk+kk)*Nl+ll] +=  5.8262136251873136e+00 * (*input)[((ii*28+ 7)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 9)*Nk+kk)*Nl+ll] +=  5.5272315570895403e+00 * (*input)[((ii*28+ 7)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+11)*Nk+kk)*Nl+ll] += -2.3666191622317520e+01 * (*input)[((ii*28+ 7)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 2)*Nk+kk)*Nl+ll] +=  2.0182596029148968e+01 * (*input)[((ii*28+ 8)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 4)*Nk+kk)*Nl+ll] += -1.4739284152238774e+01 * (*input)[((ii*28+ 8)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 7)*Nk+kk)*Nl+ll] += -1.1652427250374625e+01 * (*input)[((ii*28+ 9)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 9)*Nk+kk)*Nl+ll] +=  7.3696420761193879e+00 * (*input)[((ii*28+ 9)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 6)*Nk+kk)*Nl+ll] += -9.5353803401442638e-01 * (*input)[((ii*28+10)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 8)*Nk+kk)*Nl+ll] += -4.6060262975746175e-01 * (*input)[((ii*28+10)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+10)*Nk+kk)*Nl+ll] +=  2.5228245036436205e+00 * (*input)[((ii*28+10)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+12)*Nk+kk)*Nl+ll] +=  1.0247761577878713e+01 * (*input)[((ii*28+10)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 1)*Nk+kk)*Nl+ll] += -2.3666191622317520e+01 * (*input)[((ii*28+11)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 3)*Nk+kk)*Nl+ll] += -5.5272315570895403e+00 * (*input)[((ii*28+11)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 5)*Nk+kk)*Nl+ll] +=  5.8262136251873136e+00 * (*input)[((ii*28+11)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 6)*Nk+kk)*Nl+ll] +=  1.1442456408173115e+01 * (*input)[((ii*28+12)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+10)*Nk+kk)*Nl+ll] += -3.0273894043723452e+01 * (*input)[((ii*28+12)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 3)*Nk+kk)*Nl+ll] +=  2.2108926228358165e+01 * (*input)[((ii*28+13)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 5)*Nk+kk)*Nl+ll] += -1.1652427250374625e+01 * (*input)[((ii*28+13)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 6)*Nk+kk)*Nl+ll] += -7.6283042721154128e+00 * (*input)[((ii*28+14)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 8)*Nk+kk)*Nl+ll] +=  7.3696420761193888e+00 * (*input)[((ii*28+14)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 0)*Nk+kk)*Nl+ll] +=  4.0991046311514854e+00 * (*input)[((ii*28+15)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 2)*Nk+kk)*Nl+ll] +=  2.0182596029148967e+00 * (*input)[((ii*28+15)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 4)*Nk+kk)*Nl+ll] +=  9.2120525951492349e-01 * (*input)[((ii*28+15)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 7)*Nk+kk)*Nl+ll] +=  2.9131068125936568e+00 * (*input)[((ii*28+16)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 9)*Nk+kk)*Nl+ll] +=  8.2908473356343109e+00 * (*input)[((ii*28+16)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+11)*Nk+kk)*Nl+ll] +=  1.1833095811158760e+01 * (*input)[((ii*28+16)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 2)*Nk+kk)*Nl+ll] += -2.0182596029148968e+01 * (*input)[((ii*28+17)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 4)*Nk+kk)*Nl+ll] += -1.4739284152238774e+01 * (*input)[((ii*28+17)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 7)*Nk+kk)*Nl+ll] += -1.1652427250374625e+01 * (*input)[((ii*28+18)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 9)*Nk+kk)*Nl+ll] += -2.2108926228358165e+01 * (*input)[((ii*28+18)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 4)*Nk+kk)*Nl+ll] +=  1.4739284152238778e+01 * (*input)[((ii*28+19)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 7)*Nk+kk)*Nl+ll] +=  4.6609709001498523e+00 * (*input)[((ii*28+20)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 6)*Nk+kk)*Nl+ll] += -3.1784601133814211e-01 * (*input)[((ii*28+21)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 8)*Nk+kk)*Nl+ll] += -4.6060262975746175e-01 * (*input)[((ii*28+21)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+10)*Nk+kk)*Nl+ll] += -5.0456490072872417e-01 * (*input)[((ii*28+21)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+12)*Nk+kk)*Nl+ll] += -6.8318410519191419e-01 * (*input)[((ii*28+21)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 1)*Nk+kk)*Nl+ll] +=  2.3666191622317521e+00 * (*input)[((ii*28+22)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 3)*Nk+kk)*Nl+ll] +=  2.7636157785447706e+00 * (*input)[((ii*28+22)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 5)*Nk+kk)*Nl+ll] +=  2.9131068125936568e+00 * (*input)[((ii*28+22)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 6)*Nk+kk)*Nl+ll] +=  5.7212282040865574e+00 * (*input)[((ii*28+23)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 8)*Nk+kk)*Nl+ll] +=  7.3696420761193870e+00 * (*input)[((ii*28+23)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+10)*Nk+kk)*Nl+ll] +=  5.0456490072872420e+00 * (*input)[((ii*28+23)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 3)*Nk+kk)*Nl+ll] += -7.3696420761193879e+00 * (*input)[((ii*28+24)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 5)*Nk+kk)*Nl+ll] += -1.1652427250374625e+01 * (*input)[((ii*28+24)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 6)*Nk+kk)*Nl+ll] += -7.6283042721154128e+00 * (*input)[((ii*28+25)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 8)*Nk+kk)*Nl+ll] += -7.3696420761193888e+00 * (*input)[((ii*28+25)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 5)*Nk+kk)*Nl+ll] +=  4.6609709001498523e+00 * (*input)[((ii*28+26)*Nk+kk)*Nl+ll];
        (*output)[((ii*13+ 6)*Nk+kk)*Nl+ll] +=  1.0171072362820530e+00 * (*input)[((ii*28+27)*Nk+kk)*Nl+ll];
      }
}

void IntegralWorker::transform_j(int am, size_t Ni, size_t Nk, size_t Nl) {
  static void (*f[7])(size_t, size_t, size_t, const std::vector<double> *, std::vector<double> *)={
    transform_j0,
    transform_j1,
    transform_j2,
    transform_j3,
    transform_j4,
    transform_j5,
    transform_j6
  };
    f[am](Ni,Nk,Nl,input,output);
  std::swap(input,output);
}
static void transform_k0(size_t Ni, size_t Nj, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(1*Ni*Nj*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii*Nj+jj)* 1+ 0)*Nl+ll] +=  2.8209479177387814e-01 * (*input)[((ii*Nj+jj)* 1+ 0)*Nl+ll];
      }
}

static void transform_k1(size_t Ni, size_t Nj, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(3*Ni*Nj*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii*Nj+jj)* 3+ 2)*Nl+ll] +=  4.8860251190291992e-01 * (*input)[((ii*Nj+jj)* 3+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)* 3+ 0)*Nl+ll] +=  4.8860251190291992e-01 * (*input)[((ii*Nj+jj)* 3+ 1)*Nl+ll];
        (*output)[((ii*Nj+jj)* 3+ 1)*Nl+ll] +=  4.8860251190291992e-01 * (*input)[((ii*Nj+jj)* 3+ 2)*Nl+ll];
      }
}

static void transform_k2(size_t Ni, size_t Nj, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(5*Ni*Nj*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii*Nj+jj)* 5+ 2)*Nl+ll] += -3.1539156525252005e-01 * (*input)[((ii*Nj+jj)* 6+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)* 5+ 4)*Nl+ll] +=  5.4627421529603959e-01 * (*input)[((ii*Nj+jj)* 6+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)* 5+ 0)*Nl+ll] +=  1.0925484305920792e+00 * (*input)[((ii*Nj+jj)* 6+ 1)*Nl+ll];
        (*output)[((ii*Nj+jj)* 5+ 3)*Nl+ll] +=  1.0925484305920792e+00 * (*input)[((ii*Nj+jj)* 6+ 2)*Nl+ll];
        (*output)[((ii*Nj+jj)* 5+ 2)*Nl+ll] += -3.1539156525252005e-01 * (*input)[((ii*Nj+jj)* 6+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)* 5+ 4)*Nl+ll] += -5.4627421529603959e-01 * (*input)[((ii*Nj+jj)* 6+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)* 5+ 1)*Nl+ll] +=  1.0925484305920792e+00 * (*input)[((ii*Nj+jj)* 6+ 4)*Nl+ll];
        (*output)[((ii*Nj+jj)* 5+ 2)*Nl+ll] +=  6.3078313050503998e-01 * (*input)[((ii*Nj+jj)* 6+ 5)*Nl+ll];
      }
}

static void transform_k3(size_t Ni, size_t Nj, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(7*Ni*Nj*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii*Nj+jj)* 7+ 4)*Nl+ll] += -4.5704579946446572e-01 * (*input)[((ii*Nj+jj)*10+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 6)*Nl+ll] +=  5.9004358992664352e-01 * (*input)[((ii*Nj+jj)*10+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 0)*Nl+ll] +=  1.7701307697799304e+00 * (*input)[((ii*Nj+jj)*10+ 1)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 2)*Nl+ll] += -4.5704579946446572e-01 * (*input)[((ii*Nj+jj)*10+ 1)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 3)*Nl+ll] += -1.1195289977703462e+00 * (*input)[((ii*Nj+jj)*10+ 2)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 5)*Nl+ll] +=  1.4453057213202769e+00 * (*input)[((ii*Nj+jj)*10+ 2)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 4)*Nl+ll] += -4.5704579946446572e-01 * (*input)[((ii*Nj+jj)*10+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 6)*Nl+ll] += -1.7701307697799304e+00 * (*input)[((ii*Nj+jj)*10+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 1)*Nl+ll] +=  2.8906114426405538e+00 * (*input)[((ii*Nj+jj)*10+ 4)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 4)*Nl+ll] +=  1.8281831978578631e+00 * (*input)[((ii*Nj+jj)*10+ 5)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 0)*Nl+ll] += -5.9004358992664352e-01 * (*input)[((ii*Nj+jj)*10+ 6)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 2)*Nl+ll] += -4.5704579946446572e-01 * (*input)[((ii*Nj+jj)*10+ 6)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 3)*Nl+ll] += -1.1195289977703462e+00 * (*input)[((ii*Nj+jj)*10+ 7)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 5)*Nl+ll] += -1.4453057213202769e+00 * (*input)[((ii*Nj+jj)*10+ 7)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 2)*Nl+ll] +=  1.8281831978578631e+00 * (*input)[((ii*Nj+jj)*10+ 8)*Nl+ll];
        (*output)[((ii*Nj+jj)* 7+ 3)*Nl+ll] +=  7.4635266518023080e-01 * (*input)[((ii*Nj+jj)*10+ 9)*Nl+ll];
      }
}

static void transform_k4(size_t Ni, size_t Nj, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(9*Ni*Nj*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii*Nj+jj)* 9+ 4)*Nl+ll] +=  3.1735664074561293e-01 * (*input)[((ii*Nj+jj)*15+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 6)*Nl+ll] += -4.7308734787877998e-01 * (*input)[((ii*Nj+jj)*15+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 8)*Nl+ll] +=  6.2583573544917614e-01 * (*input)[((ii*Nj+jj)*15+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 0)*Nl+ll] +=  2.5033429417967046e+00 * (*input)[((ii*Nj+jj)*15+ 1)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 2)*Nl+ll] += -9.4617469575755997e-01 * (*input)[((ii*Nj+jj)*15+ 1)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 5)*Nl+ll] += -2.0071396306718676e+00 * (*input)[((ii*Nj+jj)*15+ 2)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 7)*Nl+ll] +=  1.7701307697799307e+00 * (*input)[((ii*Nj+jj)*15+ 2)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 4)*Nl+ll] +=  6.3471328149122586e-01 * (*input)[((ii*Nj+jj)*15+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 8)*Nl+ll] += -3.7550144126950569e+00 * (*input)[((ii*Nj+jj)*15+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 1)*Nl+ll] +=  5.3103923093397922e+00 * (*input)[((ii*Nj+jj)*15+ 4)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 3)*Nl+ll] += -2.0071396306718676e+00 * (*input)[((ii*Nj+jj)*15+ 4)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 4)*Nl+ll] += -2.5388531259649034e+00 * (*input)[((ii*Nj+jj)*15+ 5)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 6)*Nl+ll] +=  2.8385240872726798e+00 * (*input)[((ii*Nj+jj)*15+ 5)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 0)*Nl+ll] += -2.5033429417967046e+00 * (*input)[((ii*Nj+jj)*15+ 6)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 2)*Nl+ll] += -9.4617469575755997e-01 * (*input)[((ii*Nj+jj)*15+ 6)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 5)*Nl+ll] += -2.0071396306718676e+00 * (*input)[((ii*Nj+jj)*15+ 7)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 7)*Nl+ll] += -5.3103923093397922e+00 * (*input)[((ii*Nj+jj)*15+ 7)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 2)*Nl+ll] +=  5.6770481745453596e+00 * (*input)[((ii*Nj+jj)*15+ 8)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 5)*Nl+ll] +=  2.6761861742291573e+00 * (*input)[((ii*Nj+jj)*15+ 9)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 4)*Nl+ll] +=  3.1735664074561293e-01 * (*input)[((ii*Nj+jj)*15+10)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 6)*Nl+ll] +=  4.7308734787877998e-01 * (*input)[((ii*Nj+jj)*15+10)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 8)*Nl+ll] +=  6.2583573544917614e-01 * (*input)[((ii*Nj+jj)*15+10)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 1)*Nl+ll] += -1.7701307697799307e+00 * (*input)[((ii*Nj+jj)*15+11)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 3)*Nl+ll] += -2.0071396306718676e+00 * (*input)[((ii*Nj+jj)*15+11)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 4)*Nl+ll] += -2.5388531259649034e+00 * (*input)[((ii*Nj+jj)*15+12)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 6)*Nl+ll] += -2.8385240872726798e+00 * (*input)[((ii*Nj+jj)*15+12)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 3)*Nl+ll] +=  2.6761861742291573e+00 * (*input)[((ii*Nj+jj)*15+13)*Nl+ll];
        (*output)[((ii*Nj+jj)* 9+ 4)*Nl+ll] +=  8.4628437532163425e-01 * (*input)[((ii*Nj+jj)*15+14)*Nl+ll];
      }
}

static void transform_k5(size_t Ni, size_t Nj, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(11*Ni*Nj*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii*Nj+jj)*11+ 6)*Nl+ll] +=  4.5294665119569688e-01 * (*input)[((ii*Nj+jj)*21+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 8)*Nl+ll] += -4.8923829943525038e-01 * (*input)[((ii*Nj+jj)*21+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+10)*Nl+ll] +=  6.5638205684017015e-01 * (*input)[((ii*Nj+jj)*21+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 0)*Nl+ll] +=  3.2819102842008507e+00 * (*input)[((ii*Nj+jj)*21+ 1)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 2)*Nl+ll] += -1.4677148983057511e+00 * (*input)[((ii*Nj+jj)*21+ 1)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 4)*Nl+ll] +=  4.5294665119569688e-01 * (*input)[((ii*Nj+jj)*21+ 1)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 5)*Nl+ll] +=  1.7542548368013540e+00 * (*input)[((ii*Nj+jj)*21+ 2)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 7)*Nl+ll] += -2.3967683924866621e+00 * (*input)[((ii*Nj+jj)*21+ 2)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 9)*Nl+ll] +=  2.0756623148810416e+00 * (*input)[((ii*Nj+jj)*21+ 2)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 6)*Nl+ll] +=  9.0589330239139376e-01 * (*input)[((ii*Nj+jj)*21+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 8)*Nl+ll] +=  9.7847659887050065e-01 * (*input)[((ii*Nj+jj)*21+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+10)*Nl+ll] += -6.5638205684017015e+00 * (*input)[((ii*Nj+jj)*21+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 1)*Nl+ll] +=  8.3026492595241663e+00 * (*input)[((ii*Nj+jj)*21+ 4)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 3)*Nl+ll] += -4.7935367849733241e+00 * (*input)[((ii*Nj+jj)*21+ 4)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 6)*Nl+ll] += -5.4353598143483630e+00 * (*input)[((ii*Nj+jj)*21+ 5)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 8)*Nl+ll] +=  3.9139063954820035e+00 * (*input)[((ii*Nj+jj)*21+ 5)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 0)*Nl+ll] += -6.5638205684017015e+00 * (*input)[((ii*Nj+jj)*21+ 6)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 2)*Nl+ll] += -9.7847659887050065e-01 * (*input)[((ii*Nj+jj)*21+ 6)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 4)*Nl+ll] +=  9.0589330239139376e-01 * (*input)[((ii*Nj+jj)*21+ 6)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 5)*Nl+ll] +=  3.5085096736027079e+00 * (*input)[((ii*Nj+jj)*21+ 7)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 9)*Nl+ll] += -1.2453973889286249e+01 * (*input)[((ii*Nj+jj)*21+ 7)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 2)*Nl+ll] +=  1.1741719186446010e+01 * (*input)[((ii*Nj+jj)*21+ 8)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 4)*Nl+ll] += -5.4353598143483630e+00 * (*input)[((ii*Nj+jj)*21+ 8)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 5)*Nl+ll] += -4.6780128981369433e+00 * (*input)[((ii*Nj+jj)*21+ 9)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 7)*Nl+ll] +=  4.7935367849733250e+00 * (*input)[((ii*Nj+jj)*21+ 9)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 6)*Nl+ll] +=  4.5294665119569688e-01 * (*input)[((ii*Nj+jj)*21+10)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 8)*Nl+ll] +=  1.4677148983057511e+00 * (*input)[((ii*Nj+jj)*21+10)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+10)*Nl+ll] +=  3.2819102842008507e+00 * (*input)[((ii*Nj+jj)*21+10)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 1)*Nl+ll] += -8.3026492595241663e+00 * (*input)[((ii*Nj+jj)*21+11)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 3)*Nl+ll] += -4.7935367849733241e+00 * (*input)[((ii*Nj+jj)*21+11)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 6)*Nl+ll] += -5.4353598143483630e+00 * (*input)[((ii*Nj+jj)*21+12)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 8)*Nl+ll] += -1.1741719186446010e+01 * (*input)[((ii*Nj+jj)*21+12)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 3)*Nl+ll] +=  9.5870735699466501e+00 * (*input)[((ii*Nj+jj)*21+13)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 6)*Nl+ll] +=  3.6235732095655746e+00 * (*input)[((ii*Nj+jj)*21+14)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 0)*Nl+ll] +=  6.5638205684017015e-01 * (*input)[((ii*Nj+jj)*21+15)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 2)*Nl+ll] +=  4.8923829943525038e-01 * (*input)[((ii*Nj+jj)*21+15)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 4)*Nl+ll] +=  4.5294665119569688e-01 * (*input)[((ii*Nj+jj)*21+15)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 5)*Nl+ll] +=  1.7542548368013540e+00 * (*input)[((ii*Nj+jj)*21+16)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 7)*Nl+ll] +=  2.3967683924866621e+00 * (*input)[((ii*Nj+jj)*21+16)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 9)*Nl+ll] +=  2.0756623148810416e+00 * (*input)[((ii*Nj+jj)*21+16)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 2)*Nl+ll] += -3.9139063954820035e+00 * (*input)[((ii*Nj+jj)*21+17)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 4)*Nl+ll] += -5.4353598143483630e+00 * (*input)[((ii*Nj+jj)*21+17)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 5)*Nl+ll] += -4.6780128981369433e+00 * (*input)[((ii*Nj+jj)*21+18)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 7)*Nl+ll] += -4.7935367849733250e+00 * (*input)[((ii*Nj+jj)*21+18)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 4)*Nl+ll] +=  3.6235732095655746e+00 * (*input)[((ii*Nj+jj)*21+19)*Nl+ll];
        (*output)[((ii*Nj+jj)*11+ 5)*Nl+ll] +=  9.3560257962738991e-01 * (*input)[((ii*Nj+jj)*21+20)*Nl+ll];
      }
}

static void transform_k6(size_t Ni, size_t Nj, size_t Nl, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(13*Ni*Nj*Nl,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t ll=0;ll<Nl;ll++) {
        (*output)[((ii*Nj+jj)*13+ 6)*Nl+ll] += -3.1784601133814211e-01 * (*input)[((ii*Nj+jj)*28+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 8)*Nl+ll] +=  4.6060262975746175e-01 * (*input)[((ii*Nj+jj)*28+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+10)*Nl+ll] += -5.0456490072872417e-01 * (*input)[((ii*Nj+jj)*28+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+12)*Nl+ll] +=  6.8318410519191419e-01 * (*input)[((ii*Nj+jj)*28+ 0)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 0)*Nl+ll] +=  4.0991046311514854e+00 * (*input)[((ii*Nj+jj)*28+ 1)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 2)*Nl+ll] += -2.0182596029148967e+00 * (*input)[((ii*Nj+jj)*28+ 1)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 4)*Nl+ll] +=  9.2120525951492349e-01 * (*input)[((ii*Nj+jj)*28+ 1)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 7)*Nl+ll] +=  2.9131068125936568e+00 * (*input)[((ii*Nj+jj)*28+ 2)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 9)*Nl+ll] += -2.7636157785447706e+00 * (*input)[((ii*Nj+jj)*28+ 2)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+11)*Nl+ll] +=  2.3666191622317521e+00 * (*input)[((ii*Nj+jj)*28+ 2)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 6)*Nl+ll] += -9.5353803401442638e-01 * (*input)[((ii*Nj+jj)*28+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 8)*Nl+ll] +=  4.6060262975746175e-01 * (*input)[((ii*Nj+jj)*28+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+10)*Nl+ll] +=  2.5228245036436205e+00 * (*input)[((ii*Nj+jj)*28+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+12)*Nl+ll] += -1.0247761577878713e+01 * (*input)[((ii*Nj+jj)*28+ 3)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 1)*Nl+ll] +=  1.1833095811158760e+01 * (*input)[((ii*Nj+jj)*28+ 4)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 3)*Nl+ll] += -8.2908473356343109e+00 * (*input)[((ii*Nj+jj)*28+ 4)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 5)*Nl+ll] +=  2.9131068125936568e+00 * (*input)[((ii*Nj+jj)*28+ 4)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 6)*Nl+ll] +=  5.7212282040865574e+00 * (*input)[((ii*Nj+jj)*28+ 5)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 8)*Nl+ll] += -7.3696420761193870e+00 * (*input)[((ii*Nj+jj)*28+ 5)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+10)*Nl+ll] +=  5.0456490072872420e+00 * (*input)[((ii*Nj+jj)*28+ 5)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 0)*Nl+ll] += -1.3663682103838283e+01 * (*input)[((ii*Nj+jj)*28+ 6)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 4)*Nl+ll] +=  1.8424105190298470e+00 * (*input)[((ii*Nj+jj)*28+ 6)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 7)*Nl+ll] +=  5.8262136251873136e+00 * (*input)[((ii*Nj+jj)*28+ 7)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 9)*Nl+ll] +=  5.5272315570895403e+00 * (*input)[((ii*Nj+jj)*28+ 7)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+11)*Nl+ll] += -2.3666191622317520e+01 * (*input)[((ii*Nj+jj)*28+ 7)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 2)*Nl+ll] +=  2.0182596029148968e+01 * (*input)[((ii*Nj+jj)*28+ 8)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 4)*Nl+ll] += -1.4739284152238774e+01 * (*input)[((ii*Nj+jj)*28+ 8)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 7)*Nl+ll] += -1.1652427250374625e+01 * (*input)[((ii*Nj+jj)*28+ 9)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 9)*Nl+ll] +=  7.3696420761193879e+00 * (*input)[((ii*Nj+jj)*28+ 9)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 6)*Nl+ll] += -9.5353803401442638e-01 * (*input)[((ii*Nj+jj)*28+10)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 8)*Nl+ll] += -4.6060262975746175e-01 * (*input)[((ii*Nj+jj)*28+10)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+10)*Nl+ll] +=  2.5228245036436205e+00 * (*input)[((ii*Nj+jj)*28+10)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+12)*Nl+ll] +=  1.0247761577878713e+01 * (*input)[((ii*Nj+jj)*28+10)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 1)*Nl+ll] += -2.3666191622317520e+01 * (*input)[((ii*Nj+jj)*28+11)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 3)*Nl+ll] += -5.5272315570895403e+00 * (*input)[((ii*Nj+jj)*28+11)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 5)*Nl+ll] +=  5.8262136251873136e+00 * (*input)[((ii*Nj+jj)*28+11)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 6)*Nl+ll] +=  1.1442456408173115e+01 * (*input)[((ii*Nj+jj)*28+12)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+10)*Nl+ll] += -3.0273894043723452e+01 * (*input)[((ii*Nj+jj)*28+12)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 3)*Nl+ll] +=  2.2108926228358165e+01 * (*input)[((ii*Nj+jj)*28+13)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 5)*Nl+ll] += -1.1652427250374625e+01 * (*input)[((ii*Nj+jj)*28+13)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 6)*Nl+ll] += -7.6283042721154128e+00 * (*input)[((ii*Nj+jj)*28+14)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 8)*Nl+ll] +=  7.3696420761193888e+00 * (*input)[((ii*Nj+jj)*28+14)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 0)*Nl+ll] +=  4.0991046311514854e+00 * (*input)[((ii*Nj+jj)*28+15)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 2)*Nl+ll] +=  2.0182596029148967e+00 * (*input)[((ii*Nj+jj)*28+15)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 4)*Nl+ll] +=  9.2120525951492349e-01 * (*input)[((ii*Nj+jj)*28+15)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 7)*Nl+ll] +=  2.9131068125936568e+00 * (*input)[((ii*Nj+jj)*28+16)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 9)*Nl+ll] +=  8.2908473356343109e+00 * (*input)[((ii*Nj+jj)*28+16)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+11)*Nl+ll] +=  1.1833095811158760e+01 * (*input)[((ii*Nj+jj)*28+16)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 2)*Nl+ll] += -2.0182596029148968e+01 * (*input)[((ii*Nj+jj)*28+17)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 4)*Nl+ll] += -1.4739284152238774e+01 * (*input)[((ii*Nj+jj)*28+17)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 7)*Nl+ll] += -1.1652427250374625e+01 * (*input)[((ii*Nj+jj)*28+18)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 9)*Nl+ll] += -2.2108926228358165e+01 * (*input)[((ii*Nj+jj)*28+18)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 4)*Nl+ll] +=  1.4739284152238778e+01 * (*input)[((ii*Nj+jj)*28+19)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 7)*Nl+ll] +=  4.6609709001498523e+00 * (*input)[((ii*Nj+jj)*28+20)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 6)*Nl+ll] += -3.1784601133814211e-01 * (*input)[((ii*Nj+jj)*28+21)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 8)*Nl+ll] += -4.6060262975746175e-01 * (*input)[((ii*Nj+jj)*28+21)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+10)*Nl+ll] += -5.0456490072872417e-01 * (*input)[((ii*Nj+jj)*28+21)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+12)*Nl+ll] += -6.8318410519191419e-01 * (*input)[((ii*Nj+jj)*28+21)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 1)*Nl+ll] +=  2.3666191622317521e+00 * (*input)[((ii*Nj+jj)*28+22)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 3)*Nl+ll] +=  2.7636157785447706e+00 * (*input)[((ii*Nj+jj)*28+22)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 5)*Nl+ll] +=  2.9131068125936568e+00 * (*input)[((ii*Nj+jj)*28+22)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 6)*Nl+ll] +=  5.7212282040865574e+00 * (*input)[((ii*Nj+jj)*28+23)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 8)*Nl+ll] +=  7.3696420761193870e+00 * (*input)[((ii*Nj+jj)*28+23)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+10)*Nl+ll] +=  5.0456490072872420e+00 * (*input)[((ii*Nj+jj)*28+23)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 3)*Nl+ll] += -7.3696420761193879e+00 * (*input)[((ii*Nj+jj)*28+24)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 5)*Nl+ll] += -1.1652427250374625e+01 * (*input)[((ii*Nj+jj)*28+24)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 6)*Nl+ll] += -7.6283042721154128e+00 * (*input)[((ii*Nj+jj)*28+25)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 8)*Nl+ll] += -7.3696420761193888e+00 * (*input)[((ii*Nj+jj)*28+25)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 5)*Nl+ll] +=  4.6609709001498523e+00 * (*input)[((ii*Nj+jj)*28+26)*Nl+ll];
        (*output)[((ii*Nj+jj)*13+ 6)*Nl+ll] +=  1.0171072362820530e+00 * (*input)[((ii*Nj+jj)*28+27)*Nl+ll];
      }
}

void IntegralWorker::transform_k(int am, size_t Ni, size_t Nj, size_t Nl) {
  static void (*f[7])(size_t, size_t, size_t, const std::vector<double> *, std::vector<double> *)={
    transform_k0,
    transform_k1,
    transform_k2,
    transform_k3,
    transform_k4,
    transform_k5,
    transform_k6
  };
    f[am](Ni,Nj,Nl,input,output);
  std::swap(input,output);
}
static void transform_l0(size_t Ni, size_t Nj, size_t Nk, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(1*Ni*Nj*Nk,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t kk=0;kk<Nk;kk++) {
        (*output)[((ii*Nj+jj)*Nk+kk)* 1+ 0] +=  2.8209479177387814e-01 * (*input)[((ii*Nj+jj)*Nk+kk)* 1+ 0];
      }
}

static void transform_l1(size_t Ni, size_t Nj, size_t Nk, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(3*Ni*Nj*Nk,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t kk=0;kk<Nk;kk++) {
        (*output)[((ii*Nj+jj)*Nk+kk)* 3+ 2] +=  4.8860251190291992e-01 * (*input)[((ii*Nj+jj)*Nk+kk)* 3+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)* 3+ 0] +=  4.8860251190291992e-01 * (*input)[((ii*Nj+jj)*Nk+kk)* 3+ 1];
        (*output)[((ii*Nj+jj)*Nk+kk)* 3+ 1] +=  4.8860251190291992e-01 * (*input)[((ii*Nj+jj)*Nk+kk)* 3+ 2];
      }
}

static void transform_l2(size_t Ni, size_t Nj, size_t Nk, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(5*Ni*Nj*Nk,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t kk=0;kk<Nk;kk++) {
        (*output)[((ii*Nj+jj)*Nk+kk)* 5+ 2] += -3.1539156525252005e-01 * (*input)[((ii*Nj+jj)*Nk+kk)* 6+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)* 5+ 4] +=  5.4627421529603959e-01 * (*input)[((ii*Nj+jj)*Nk+kk)* 6+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)* 5+ 0] +=  1.0925484305920792e+00 * (*input)[((ii*Nj+jj)*Nk+kk)* 6+ 1];
        (*output)[((ii*Nj+jj)*Nk+kk)* 5+ 3] +=  1.0925484305920792e+00 * (*input)[((ii*Nj+jj)*Nk+kk)* 6+ 2];
        (*output)[((ii*Nj+jj)*Nk+kk)* 5+ 2] += -3.1539156525252005e-01 * (*input)[((ii*Nj+jj)*Nk+kk)* 6+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)* 5+ 4] += -5.4627421529603959e-01 * (*input)[((ii*Nj+jj)*Nk+kk)* 6+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)* 5+ 1] +=  1.0925484305920792e+00 * (*input)[((ii*Nj+jj)*Nk+kk)* 6+ 4];
        (*output)[((ii*Nj+jj)*Nk+kk)* 5+ 2] +=  6.3078313050503998e-01 * (*input)[((ii*Nj+jj)*Nk+kk)* 6+ 5];
      }
}

static void transform_l3(size_t Ni, size_t Nj, size_t Nk, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(7*Ni*Nj*Nk,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t kk=0;kk<Nk;kk++) {
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 4] += -4.5704579946446572e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 6] +=  5.9004358992664352e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 0] +=  1.7701307697799304e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 1];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 2] += -4.5704579946446572e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 1];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 3] += -1.1195289977703462e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 2];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 5] +=  1.4453057213202769e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 2];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 4] += -4.5704579946446572e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 6] += -1.7701307697799304e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 1] +=  2.8906114426405538e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 4];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 4] +=  1.8281831978578631e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 5];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 0] += -5.9004358992664352e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 6];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 2] += -4.5704579946446572e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 6];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 3] += -1.1195289977703462e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 7];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 5] += -1.4453057213202769e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 7];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 2] +=  1.8281831978578631e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 8];
        (*output)[((ii*Nj+jj)*Nk+kk)* 7+ 3] +=  7.4635266518023080e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*10+ 9];
      }
}

static void transform_l4(size_t Ni, size_t Nj, size_t Nk, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(9*Ni*Nj*Nk,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t kk=0;kk<Nk;kk++) {
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 4] +=  3.1735664074561293e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 6] += -4.7308734787877998e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 8] +=  6.2583573544917614e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 0] +=  2.5033429417967046e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 1];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 2] += -9.4617469575755997e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 1];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 5] += -2.0071396306718676e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 2];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 7] +=  1.7701307697799307e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 2];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 4] +=  6.3471328149122586e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 8] += -3.7550144126950569e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 1] +=  5.3103923093397922e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 4];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 3] += -2.0071396306718676e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 4];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 4] += -2.5388531259649034e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 5];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 6] +=  2.8385240872726798e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 5];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 0] += -2.5033429417967046e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 6];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 2] += -9.4617469575755997e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 6];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 5] += -2.0071396306718676e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 7];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 7] += -5.3103923093397922e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 7];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 2] +=  5.6770481745453596e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 8];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 5] +=  2.6761861742291573e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+ 9];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 4] +=  3.1735664074561293e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*15+10];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 6] +=  4.7308734787877998e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*15+10];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 8] +=  6.2583573544917614e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*15+10];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 1] += -1.7701307697799307e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+11];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 3] += -2.0071396306718676e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+11];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 4] += -2.5388531259649034e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+12];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 6] += -2.8385240872726798e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+12];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 3] +=  2.6761861742291573e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*15+13];
        (*output)[((ii*Nj+jj)*Nk+kk)* 9+ 4] +=  8.4628437532163425e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*15+14];
      }
}

static void transform_l5(size_t Ni, size_t Nj, size_t Nk, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(11*Ni*Nj*Nk,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t kk=0;kk<Nk;kk++) {
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 6] +=  4.5294665119569688e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 8] += -4.8923829943525038e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+10] +=  6.5638205684017015e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 0] +=  3.2819102842008507e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 1];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 2] += -1.4677148983057511e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 1];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 4] +=  4.5294665119569688e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 1];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 5] +=  1.7542548368013540e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 2];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 7] += -2.3967683924866621e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 2];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 9] +=  2.0756623148810416e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 2];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 6] +=  9.0589330239139376e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 8] +=  9.7847659887050065e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+10] += -6.5638205684017015e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 1] +=  8.3026492595241663e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 4];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 3] += -4.7935367849733241e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 4];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 6] += -5.4353598143483630e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 5];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 8] +=  3.9139063954820035e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 5];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 0] += -6.5638205684017015e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 6];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 2] += -9.7847659887050065e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 6];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 4] +=  9.0589330239139376e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 6];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 5] +=  3.5085096736027079e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 7];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 9] += -1.2453973889286249e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 7];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 2] +=  1.1741719186446010e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 8];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 4] += -5.4353598143483630e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 8];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 5] += -4.6780128981369433e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 9];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 7] +=  4.7935367849733250e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+ 9];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 6] +=  4.5294665119569688e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+10];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 8] +=  1.4677148983057511e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+10];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+10] +=  3.2819102842008507e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+10];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 1] += -8.3026492595241663e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+11];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 3] += -4.7935367849733241e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+11];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 6] += -5.4353598143483630e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+12];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 8] += -1.1741719186446010e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+12];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 3] +=  9.5870735699466501e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+13];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 6] +=  3.6235732095655746e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+14];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 0] +=  6.5638205684017015e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+15];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 2] +=  4.8923829943525038e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+15];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 4] +=  4.5294665119569688e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+15];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 5] +=  1.7542548368013540e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+16];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 7] +=  2.3967683924866621e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+16];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 9] +=  2.0756623148810416e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+16];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 2] += -3.9139063954820035e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+17];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 4] += -5.4353598143483630e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+17];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 5] += -4.6780128981369433e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+18];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 7] += -4.7935367849733250e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+18];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 4] +=  3.6235732095655746e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*21+19];
        (*output)[((ii*Nj+jj)*Nk+kk)*11+ 5] +=  9.3560257962738991e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*21+20];
      }
}

static void transform_l6(size_t Ni, size_t Nj, size_t Nk, const std::vector<double> *input, std::vector<double> *output) {
  (*output).clear();
  (*output).resize(13*Ni*Nj*Nk,0.0);
  for(size_t ii=0;ii<Ni;ii++)
    for(size_t jj=0;jj<Nj;jj++)
      for(size_t kk=0;kk<Nk;kk++) {
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 6] += -3.1784601133814211e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 8] +=  4.6060262975746175e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+10] += -5.0456490072872417e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+12] +=  6.8318410519191419e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 0];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 0] +=  4.0991046311514854e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 1];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 2] += -2.0182596029148967e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 1];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 4] +=  9.2120525951492349e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 1];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 7] +=  2.9131068125936568e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 2];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 9] += -2.7636157785447706e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 2];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+11] +=  2.3666191622317521e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 2];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 6] += -9.5353803401442638e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 8] +=  4.6060262975746175e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+10] +=  2.5228245036436205e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+12] += -1.0247761577878713e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 3];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 1] +=  1.1833095811158760e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 4];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 3] += -8.2908473356343109e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 4];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 5] +=  2.9131068125936568e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 4];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 6] +=  5.7212282040865574e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 5];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 8] += -7.3696420761193870e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 5];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+10] +=  5.0456490072872420e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 5];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 0] += -1.3663682103838283e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 6];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 4] +=  1.8424105190298470e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 6];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 7] +=  5.8262136251873136e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 7];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 9] +=  5.5272315570895403e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 7];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+11] += -2.3666191622317520e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 7];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 2] +=  2.0182596029148968e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 8];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 4] += -1.4739284152238774e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 8];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 7] += -1.1652427250374625e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 9];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 9] +=  7.3696420761193879e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+ 9];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 6] += -9.5353803401442638e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+10];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 8] += -4.6060262975746175e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+10];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+10] +=  2.5228245036436205e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+10];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+12] +=  1.0247761577878713e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+10];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 1] += -2.3666191622317520e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+11];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 3] += -5.5272315570895403e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+11];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 5] +=  5.8262136251873136e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+11];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 6] +=  1.1442456408173115e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+12];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+10] += -3.0273894043723452e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+12];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 3] +=  2.2108926228358165e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+13];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 5] += -1.1652427250374625e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+13];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 6] += -7.6283042721154128e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+14];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 8] +=  7.3696420761193888e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+14];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 0] +=  4.0991046311514854e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+15];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 2] +=  2.0182596029148967e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+15];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 4] +=  9.2120525951492349e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+15];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 7] +=  2.9131068125936568e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+16];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 9] +=  8.2908473356343109e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+16];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+11] +=  1.1833095811158760e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+16];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 2] += -2.0182596029148968e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+17];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 4] += -1.4739284152238774e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+17];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 7] += -1.1652427250374625e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+18];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 9] += -2.2108926228358165e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+18];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 4] +=  1.4739284152238778e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+19];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 7] +=  4.6609709001498523e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+20];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 6] += -3.1784601133814211e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+21];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 8] += -4.6060262975746175e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+21];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+10] += -5.0456490072872417e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+21];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+12] += -6.8318410519191419e-01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+21];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 1] +=  2.3666191622317521e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+22];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 3] +=  2.7636157785447706e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+22];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 5] +=  2.9131068125936568e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+22];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 6] +=  5.7212282040865574e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+23];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 8] +=  7.3696420761193870e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+23];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+10] +=  5.0456490072872420e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+23];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 3] += -7.3696420761193879e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+24];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 5] += -1.1652427250374625e+01 * (*input)[((ii*Nj+jj)*Nk+kk)*28+24];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 6] += -7.6283042721154128e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+25];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 8] += -7.3696420761193888e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+25];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 5] +=  4.6609709001498523e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+26];
        (*output)[((ii*Nj+jj)*Nk+kk)*13+ 6] +=  1.0171072362820530e+00 * (*input)[((ii*Nj+jj)*Nk+kk)*28+27];
      }
}

void IntegralWorker::transform_l(int am, size_t Ni, size_t Nj, size_t Nk) {
  static void (*f[7])(size_t, size_t, size_t, const std::vector<double> *, std::vector<double> *)={
    transform_l0,
    transform_l1,
    transform_l2,
    transform_l3,
    transform_l4,
    transform_l5,
    transform_l6
  };
    f[am](Ni,Nj,Nk,input,output);
  std::swap(input,output);
}

/*
 *                This source code is part of
 *
 *                     E  R  K  A  L  E
 *                             -
 *                       DFT from Hel
 *
 * Written by Susi Lehtola, 2010-2013
 * Copyright (c) 2010-2013, Susi Lehtola
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 */

#include "eriworker.h"
#include "linalg.h"
#include "mathf.h"
#include "boys.h"
#include "integrals.h"
#include <cfloat>

// No Boys function interpolation?
#define BOYSNOINTERP

IntegralWorker::IntegralWorker() {
  input=&arrone;
  output=&arrtwo;

#ifndef BOYSNOINTERP
  // Maximum value of m for Boys function is
#ifdef _OPENMP
#pragma omp critical
#endif
  BoysTable::fill(4*LIBINT_MAX_AM);
#endif
}

IntegralWorker::~IntegralWorker() {
}

void IntegralWorker::spherical_transform(const GaussianShell *is, const GaussianShell *js, const GaussianShell *ks, const GaussianShell *ls) {
  const bool is_lm=is->lm_in_use();
  const bool js_lm=js->lm_in_use();
  const bool ks_lm=ks->lm_in_use();
  const bool ls_lm=ls->lm_in_use();

  const int am_i=is->get_am();
  const int am_j=js->get_am();
  const int am_k=ks->get_am();
  const int am_l=ls->get_am();

  const size_t Nj_tgt=js->get_Nbf();
  const size_t Nk_tgt=ks->get_Nbf();
  const size_t Nl_tgt=ls->get_Nbf();

  const size_t Ni_cart=is->get_Ncart();
  const size_t Nj_cart=js->get_Ncart();
  const size_t Nk_cart=ks->get_Ncart();

  // Do transforms
  if(ls_lm)
    transform_l(am_l,Ni_cart,Nj_cart,Nk_cart);
  if(ks_lm)
    transform_k(am_k,Ni_cart,Nj_cart,Nl_tgt);
  if(js_lm)
    transform_j(am_j,Ni_cart,Nk_tgt,Nl_tgt);
  if(is_lm)
    transform_i(am_i,Nj_tgt,Nk_tgt,Nl_tgt);
}

ERIWorker::ERIWorker(int maxam, int maxcontr) {
  if(maxam>=LIBINT_MAX_AM) {
    ERROR_INFO();
    throw std::domain_error("You need a version of LIBINT that supports larger angular momentum.\n");
  }

  // Initialize evaluator
  init_libint(&libint,maxam,pow(maxcontr,4));
}

ERIWorker::~ERIWorker() {
  free_libint(&libint);
}

dERIWorker::dERIWorker(int maxam, int maxcontr) {
  if(maxam>=LIBDERIV_MAX_AM1) {
    ERROR_INFO();
    throw std::domain_error("You need a version of LIBDERIV that supports larger angular momentum.\n");
  }

  // Maximum amount of basis functions is
  int nbf=(maxam+1)*(maxam+2)/2;

  // Initialize evaluator
  init_libderiv1(&libderiv,maxam,pow(maxcontr,4),pow(nbf,4));
}

dERIWorker::~dERIWorker() {
  free_libderiv(&libderiv);
}

void ERIWorker::compute_cartesian(const GaussianShell *is, const GaussianShell *js, const GaussianShell *ks, const GaussianShell *ls) {
  eri_precursor_t ip=compute_precursor(is,js);
  eri_precursor_t jp=compute_precursor(ks,ls);

  // Compute shell of cartesian ERIs using libint

  // Libint computes (ab|cd) for
  // l(a)>=l(b), l(c)>=l(d) and l(a)+l(b)>=l(c)+l(d)
  // where l(a) is the angular momentum type of the a shell.
  // THIS ROUTINE ASSUMES THAT THE SHELLS ARE ALREADY IN CORRECT ORDER!

  // The sum of angular momenta is
  int mmax=is->get_am()+js->get_am()+ks->get_am()+ls->get_am();

  // Figure out the number of contractions
  size_t Ncomb=is->get_Ncontr()*js->get_Ncontr()*ks->get_Ncontr()*ls->get_Ncontr();

  // Check that all is in order
  if(is->get_am()<js->get_am()) {
    ERROR_INFO();
    throw std::runtime_error("lambda_i < lambda_j\n");
  }

  if(ks->get_am()<ls->get_am()) {
    ERROR_INFO();
    throw std::runtime_error("lambda_k < lambda_l\n");
  }

  if( (is->get_am()+js->get_am()) > (ks->get_am()+ls->get_am())) {
    ERROR_INFO();
    throw std::runtime_error("lambda_k + lambda_l < lambda_i + lambda_j\n");
  }

  // Compute data for LIBINT
  compute_libint_data(ip,jp,mmax);

  // Pointer to integrals table
  double *ints;

  // Special handling of (ss|ss) integrals:
  if(mmax==0) {
    double tmp=0.0;
    for(size_t i=0;i<Ncomb;i++)
      tmp+=libint.PrimQuartet[i].F[0];

    // Plug in normalizations
    tmp*=is->get_cart()[0].relnorm;
    tmp*=js->get_cart()[0].relnorm;
    tmp*=ks->get_cart()[0].relnorm;
    tmp*=ls->get_cart()[0].relnorm;

    (*input).resize(1);
    (*input)[0]=tmp;
  } else {
    //printf("Computing shell %i %i %i %i\n",is->get_am(),js->get_am(),ks->get_am(),ls->get_am());
    //    printf("which consists of basis functions (%i-%i)x(%i-%i)x(%i-%i)x(%i-%i).\n",(int) shells[is].get_first_ind(),(int) shells[is].get_last_ind(),(int) shells[js].get_first_ind(),(int) shells[js].get_last_ind(),(int) shells[ks].get_first_ind(),(int) shells[ks].get_last_ind(),(int) shells[ls].get_first_ind(),(int) shells[ls].get_last_ind());

    // Now we can compute the integrals using libint:
    ints=build_eri[is->get_am()][js->get_am()][ks->get_am()][ls->get_am()](&libint,Ncomb);

    // and collect the results, plugging in the normalization factors
    size_t ind_ij, ind_ijk, ind;
    double norm_i, norm_ij, norm_ijk, norm;

    // Numbers of functions on each shell
    std::vector<shellf_t> ci=is->get_cart();
    std::vector<shellf_t> cj=js->get_cart();
    std::vector<shellf_t> ck=ks->get_cart();
    std::vector<shellf_t> cl=ls->get_cart();
    (*input).resize(ci.size()*cj.size()*ck.size()*cl.size());

    for(size_t ii=0;ii<ci.size();ii++) {
      norm_i=ci[ii].relnorm;
      for(size_t ji=0;ji<cj.size();ji++) {
	ind_ij=ii*cj.size()+ji;
	norm_ij=cj[ji].relnorm*norm_i;
	for(size_t ki=0;ki<ck.size();ki++) {
	  ind_ijk=ind_ij*ck.size()+ki;
	  norm_ijk=ck[ki].relnorm*norm_ij;
	  for(size_t li=0;li<cl.size();li++) {
	    // Index in computed integrals table
	    ind=ind_ijk*cl.size()+li;
	    // Total norm factor
	    norm=cl[li].relnorm*norm_ijk;
	    // Store scaled integral
	    //printf("Scaling integral (%i %i %i %i) = %i by % e % e % e % e\n",(int) ii,(int) ji, (int) ki, (int) li, (int) ind, ci[ii].relnorm,cj[ji].relnorm,ck[ki].relnorm,cl[li].relnorm);
	    (*input)[ind]=norm*ints[ind];
	  }
	}
      }
    }
  }
}

void dERIWorker::compute_cartesian() {
  eri_precursor_t ip=compute_precursor(is,js);
  eri_precursor_t jp=compute_precursor(ks,ls);

  // Compute shell of cartesian ERI derivatives using libderiv

  // Libint computes (ab|cd) for
  // l(a)>=l(b), l(c)>=l(d) and l(a)+l(b)>=l(c)+l(d)
  // where l(a) is the angular momentum type of the a shell.
  // THIS ROUTINE ASSUMES THAT THE SHELLS ARE ALREADY IN CORRECT ORDER!

  // The sum of angular momenta is
  int mmax=is->get_am()+js->get_am()+ks->get_am()+ls->get_am();

  // Figure out the number of contractions
  size_t Ncomb=is->get_Ncontr()*js->get_Ncontr()*ks->get_Ncontr()*ls->get_Ncontr();

  // Check that all is in order
  if(is->get_am()<js->get_am()) {
    ERROR_INFO();
    throw std::runtime_error("lambda_i < lambda_j\n");
  }

  if(ks->get_am()<ls->get_am()) {
    ERROR_INFO();
    throw std::runtime_error("lambda_k < lambda_l\n");
  }

  if( (is->get_am()+js->get_am()) > (ks->get_am()+ls->get_am())) {
    ERROR_INFO();
    throw std::runtime_error("lambda_k + lambda_l < lambda_i + lambda_j\n");
  }

  // Compute data for LIBINT
  compute_libderiv_data(ip,jp,mmax);

  // Now we can compute the integrals using libint:
  build_deriv1_eri[is->get_am()][js->get_am()][ks->get_am()][ls->get_am()](&libderiv,Ncomb);

  // and plug in the normalization factors.
  size_t ind_i, ind_ij, ind_ijk, ind;
  double norm_i, norm_ij, norm_ijk, norm;

  // Numbers of functions on each shell
  std::vector<shellf_t> ci=is->get_cart();
  std::vector<shellf_t> cj=js->get_cart();
  std::vector<shellf_t> ck=ks->get_cart();
  std::vector<shellf_t> cl=ls->get_cart();

  // Integrals computed by libderiv
  const int idx[]={0, 1, 2, 6, 7, 8, 9, 10, 11};

  for(size_t ii=0;ii<ci.size();ii++) {
    ind_i=ii*cj.size();
    norm_i=ci[ii].relnorm;
    for(size_t ji=0;ji<cj.size();ji++) {
      ind_ij=(ind_i+ji)*ck.size();
      norm_ij=cj[ji].relnorm*norm_i;
      for(size_t ki=0;ki<ck.size();ki++) {
	  ind_ijk=(ind_ij+ki)*cl.size();
	  norm_ijk=ck[ki].relnorm*norm_ij;
	  for(size_t li=0;li<cl.size();li++) {
	    // Index in computed integrals table
	    ind=ind_ijk+li;
	    // Total norm factor
	    norm=cl[li].relnorm*norm_ijk;
	    // Normalize integrals
	    for(size_t iidx=0;iidx<sizeof(idx)/sizeof(idx[0]);iidx++)
	      libderiv.ABCD[idx[iidx]][ind]*=norm;
	  }
      }
    }
  }
}

void dERIWorker::get_idx(int idx) {
  // Amount of integrals is
  size_t N=(is->get_Ncart())*(js->get_Ncart())*(ks->get_Ncart())*(ls->get_Ncart());
  (*input).resize(N);

  // Determine what is the real center requested (juggling of integrals)
  if(idx>=0 && idx<3) {
    if(swap_ij && swap_ijkl)
      // i-> j -> l
      idx+=9;
    else if(swap_ij)
      // i -> j
      idx+=3;
    else if(swap_ijkl)
      // i -> k
      idx+=6;

  } else if(idx>=3 && idx<6) {
    if(swap_ij && swap_ijkl)
      // j -> i -> k
      idx+=3;
    else if(swap_ij)
      // j -> i
      idx-=3;
    else if(swap_ijkl)
      // j -> l
      idx+=6;

  } else if(idx>=6 && idx<9) {
    if(swap_kl && swap_ijkl)
      // k -> l -> j
      idx-=3;
    else if(swap_kl)
      // k -> l
      idx+=3;
    else if(swap_ijkl)
      // k -> i
      idx-=6;

  } else if(idx>=9 && idx<12) {

    if(swap_kl && swap_ijkl)
      // l -> k -> i
      idx-=9;
    else if(swap_kl)
      // l -> k
      idx-=3;
    else if(swap_ijkl)
      // l -> j
      idx-=6;
  }

  if((idx>=0 && idx<3) || (idx>5 && idx<12)) {
    // Integrals have been computed explicitely.
    for(size_t i=0;i<N;i++)
      (*input)[i]=libderiv.ABCD[idx][i];
  } else if(idx>=3 && idx<=5) {
    // Derivate wrt B_i: d/dB_i = - d/dA_i - d/dC_i - d/dD_i
    for(size_t i=0;i<N;i++)
      (*input)[i]= - libderiv.ABCD[idx-3][i] - libderiv.ABCD[idx+3][i] - libderiv.ABCD[idx+6][i];
  } else {
    ERROR_INFO();
    throw std::runtime_error("Invalid derivative index requested!\n");
  }

  // Reorder integrals
  reorder(is_orig,js_orig,ks_orig,ls_orig,swap_ij,swap_kl,swap_ijkl);
  // and transform them into the spherical basis
  spherical_transform(is_orig,js_orig,ks_orig,ls_orig);
}

std::vector<double> dERIWorker::get(int idx) {
  // Get the integrals
  get_idx(idx);

#ifdef DEBUGDERIV
  // Evaluate other integrals
  std::vector<double> eris=get_debug(idx);

  // Check the integrals
  size_t Nfail=0;
  for(size_t i=0;i<N;i++)
    if(fabs(ints[i]-eris[i])>1e-6) {
      Nfail++;
      printf("(%c %c %c %c) integral, idx = %i, first functions (%i %i %i %i)\n",shell_types[is_orig->get_am()],shell_types[js_orig->get_am()],shell_types[ks_orig->get_am()],shell_types[ls_orig->get_am()],idx,(int) is_orig->get_first_ind(),(int) js_orig->get_first_ind(),(int) ks_orig->get_first_ind(),(int) ls_orig->get_first_ind());
      printf("%i % e % e % e\n",(int) i,ints[i],eris[i],ints[i]-eris[i]);
    }

  if(Nfail) {
    /*
      is_orig->print();
      js_orig->print();
    ks_orig->print();
    ls_orig->print();
    */

    /*
    for(size_t i=0;i<N;i++)
      printf("%i % e % e % e\n",(int) i,ints[i],eris[i],ints[i]-eris[i]);
    */
    printf("%i integrals failed\n",(int) Nfail);
  }

#endif

  return *input;
}

const std::vector<double> * dERIWorker::getp(int idx) {
  get_idx(idx);
  return input;
}


std::vector<double> dERIWorker::get_debug(int idx) {
  // Amount of integrals
  size_t N=(is->get_Ncart())*(js->get_Ncart())*(ks->get_Ncart())*(ls->get_Ncart());
  (*input).resize(N);

  // Form cartesian shells
  GaussianShell isc(is_orig->get_am(),false,is_orig->get_contr());
  isc.set_center(is_orig->get_center(),is_orig->get_center_ind());
  isc.set_first_ind(is_orig->get_first_ind());
  isc.normalize(); // Need to normalize, otherwise cartesian relative factors aren't initalized

  GaussianShell jsc(js_orig->get_am(),false,js_orig->get_contr());
  jsc.set_center(js_orig->get_center(),js_orig->get_center_ind());
  jsc.set_first_ind(js_orig->get_first_ind());
  jsc.normalize(); // Need to normalize, otherwise cartesian relative factors aren't initalized

  GaussianShell ksc(ks_orig->get_am(),false,ks_orig->get_contr());
  ksc.set_center(ks_orig->get_center(),ks_orig->get_center_ind());
  ksc.set_first_ind(ks_orig->get_first_ind());
  ksc.normalize(); // Need to normalize, otherwise cartesian relative factors aren't initalized

  GaussianShell lsc(ls_orig->get_am(),false,ls_orig->get_contr());
  lsc.set_center(ls_orig->get_center(),ls_orig->get_center_ind());
  lsc.set_first_ind(ls_orig->get_first_ind());
  lsc.normalize(); // Need to normalize, otherwise cartesian relative factors aren't initalized

  // Sizes
  size_t Ni=isc.get_Ncart();
  size_t Nj=jsc.get_Ncart();
  size_t Nk=ksc.get_Ncart();
  size_t Nl=lsc.get_Ncart();

  // ERI evaluator
  ERIWorker eri(std::max(std::max(isc.get_am(),jsc.get_am()),std::max(ksc.get_am(),lsc.get_am()))+1,std::max(std::max(isc.get_Ncontr(),jsc.get_Ncontr()),std::max(ksc.get_Ncontr(),lsc.get_Ncontr())));

  // Zero out (*input) array
  (*input).assign(N,0.0);

  // Which derivative are we taking?
  if(idx>=0 && idx<3) {
    // Get normalized contraction.
    std::vector<contr_t> ic=isc.get_contr_normalized();

    // Get cartesians
    std::vector<shellf_t> icart=isc.get_cart();

    // Loop over contraction
    for(size_t iic=0;iic<ic.size();iic++) {
      // Dummy contraction
      std::vector<contr_t> dumcontr(1);
      dumcontr[0].c=1.0;
      dumcontr[0].z=ic[iic].z;

      // Form helpers
      GaussianShell icp(is_orig->get_am()+1,false,dumcontr);
      icp.set_center(is_orig->get_center(),is_orig->get_center_ind());
      icp.set_first_ind(is_orig->get_first_ind());
      icp.normalize();

      // Evaluate ERI
      const std::vector<double> * erip;
      eri.compute(&icp,&jsc,&ksc,&lsc);
      erip=eri.getp();

      // Collect terms
      for(size_t ica=0;ica<icart.size();ica++) {
	int l=icart[ica].l;
	int m=icart[ica].m;
	int n=icart[ica].n;

	int il=l;
	int im=m;
	int in=n;

	double fac=ic[iic].c*sqrt(ic[iic].z);
	if(idx==0) {
	  fac*=sqrt(2*l+1);
	  il++;
	} else if(idx==1) {
	  fac*=sqrt(2*m+1);
	  im++;
	} else if(idx==2) {
	  fac*=sqrt(2*n+1);
	  in++;
	}

	// i index of target integral is
	size_t idt=getind(l,m,n);
	// i index of source integral is
	size_t idp=getind(il,im,in);

	// Loop over target integrals
	for(size_t jj=0;jj<Nj;jj++)
	  for(size_t kk=0;kk<Nk;kk++)
	    for(size_t ll=0;ll<Nl;ll++)
	      (*input)[((idt*Nj+jj)*Nk+kk)*Nl+ll]+=fac*(*erip)[((idp*Nj+jj)*Nk+kk)*Nl+ll];
      }

      if(is_orig->get_am()>0) {
	GaussianShell icm=GaussianShell(is_orig->get_am()-1,false,dumcontr);
	icm.set_center(is_orig->get_center(),is_orig->get_center_ind());
	icm.set_first_ind(is_orig->get_first_ind());
	icm.normalize();

	// Evaluate ERI
	eri.compute(&icm,&jsc,&ksc,&lsc);
	erip=eri.getp();

	// Collect terms
	for(size_t ica=0;ica<icart.size();ica++) {
	  int l=icart[ica].l;
	  int m=icart[ica].m;
	  int n=icart[ica].n;

	  // Skip nonexistent integrals
	  if(idx==0 && l==0)
	    continue;
	  if(idx==1 && m==0)
	    continue;
	  if(idx==2 && n==0)
	    continue;

	  int il=l;
	  int im=m;
	  int in=n;

	  double fac=-2*ic[iic].c*sqrt(ic[iic].z);
	  if(idx==0) {
	    fac*=l/sqrt(2*l-1);
	    il--;
	  } else if(idx==1) {
	    fac*=m/sqrt(2*m-1);
	    im--;
	  } else if(idx==2) {
	    fac*=n/sqrt(2*n-1);
	    in--;
	  }

	  // i index of target integral is
	  size_t idt=getind(l,m,n);
	  // i index of source integral is
	  size_t idm=getind(il,im,in);

	  // Loop over target integrals
	  for(size_t jj=0;jj<Nj;jj++)
	    for(size_t kk=0;kk<Nk;kk++)
	      for(size_t ll=0;ll<Nl;ll++)
		(*input)[((idt*Nj+jj)*Nk+kk)*Nl+ll]+=fac*(*erip)[((idm*Nj+jj)*Nk+kk)*Nl+ll];
	}
      }
    }

  } else if(idx>=3 && idx<6) {
    // Get normalized contraction.
    std::vector<contr_t> jc=jsc.get_contr_normalized();

    // Get cartesians
    std::vector<shellf_t> jcart=jsc.get_cart();

    // Loop over contraction
    for(size_t jic=0;jic<jc.size();jic++) {
      // Dummy contraction
      std::vector<contr_t> dumcontr(1);
      dumcontr[0].c=1.0;
      dumcontr[0].z=jc[jic].z;

      // Form helpers
      GaussianShell jcp(js_orig->get_am()+1,false,dumcontr);
      jcp.set_center(js_orig->get_center(),js_orig->get_center_ind());
      jcp.set_first_ind(js_orig->get_first_ind());
      jcp.normalize();
      size_t Njp=jcp.get_Ncart();

      // Evaluate ERI
      const std::vector<double> * erip;
      eri.compute(&isc,&jcp,&ksc,&lsc);
      erip=eri.getp();

      // Collect terms
      for(size_t jca=0;jca<jcart.size();jca++) {
	int l=jcart[jca].l;
	int m=jcart[jca].m;
	int n=jcart[jca].n;

	int il=l;
	int im=m;
	int in=n;

	double fac=jc[jic].c*sqrt(jc[jic].z);
	if(idx==3) {
	  fac*=sqrt(2*l+1);
	  il++;
	} else if(idx==4) {
	  fac*=sqrt(2*m+1);
	  im++;
	} else if(idx==5) {
	  fac*=sqrt(2*n+1);
	  in++;
	}

	// j index of target integral is
	size_t jdt=getind(l,m,n);
	// j index of source integral is
	size_t jdp=getind(il,im,in);

	// Loop over target integrals
	for(size_t ii=0;ii<Ni;ii++)
	  for(size_t kk=0;kk<Nk;kk++)
	    for(size_t ll=0;ll<Nl;ll++)
	      (*input)[((ii*Nj+jdt)*Nk+kk)*Nl+ll]+=fac*(*erip)[((ii*Njp+jdp)*Nk+kk)*Nl+ll];
      }

      if(js_orig->get_am()>0) {
	GaussianShell jcm=GaussianShell(js_orig->get_am()-1,false,dumcontr);
	jcm.set_center(js_orig->get_center(),js_orig->get_center_ind());
	jcm.set_first_ind(js_orig->get_first_ind());
	jcm.normalize();
	size_t Njm=jcm.get_Ncart();

	// Evaluate ERI
	eri.compute(&isc,&jcm,&ksc,&lsc);
	erip=eri.getp();

	// Collect terms
	for(size_t jca=0;jca<jcart.size();jca++) {
	  int l=jcart[jca].l;
	  int m=jcart[jca].m;
	  int n=jcart[jca].n;

	  // Skip nonexistent integrals
	  if(idx==3 && l==0)
	    continue;
	  if(idx==4 && m==0)
	    continue;
	  if(idx==5 && n==0)
	    continue;

	  int il=l;
	  int im=m;
	  int in=n;

	  double fac=-jc[jic].c*sqrt(jc[jic].z);
	  if(idx==3) {
	    fac*=2*l/sqrt(2*l-1);
	    il--;
	  } else if(idx==4) {
	    fac*=2*m/sqrt(2*m-1);
	    im--;
	  } else if(idx==5) {
	    fac*=2*n/sqrt(2*n-1);
	    in--;
	  }

	  // j index of target integral is
	  size_t jdt=getind(l,m,n);
	  // j index of source integral is
	  size_t jdm=getind(il,im,in);

	  // Loop over target integrals
	  for(size_t ii=0;ii<Ni;ii++)
	    for(size_t kk=0;kk<Nk;kk++)
	      for(size_t ll=0;ll<Nl;ll++)
		(*input)[((ii*Nj+jdt)*Nk+kk)*Nl+ll]+=fac*(*erip)[((ii*Njm+jdm)*Nk+kk)*Nl+ll];
	}
      }
    }

  } else if(idx>=6 && idx<9) {
    // Get normalized contraction.
    std::vector<contr_t> kc=ksc.get_contr_normalized();

    // Get cartesians
    std::vector<shellf_t> kcart=ksc.get_cart();

    // Loop over contraction
    for(size_t kic=0;kic<kc.size();kic++) {
      // Dummy contraction
      std::vector<contr_t> dumcontr(1);
      dumcontr[0].c=1.0;
      dumcontr[0].z=kc[kic].z;

      // Form helpers
      GaussianShell kcp(ks_orig->get_am()+1,false,dumcontr);
      kcp.set_center(ks_orig->get_center(),ks_orig->get_center_ind());
      kcp.set_first_ind(ks_orig->get_first_ind());
      kcp.normalize();
      size_t Nkp=kcp.get_Ncart();

      // Evaluate ERI
      const std::vector<double> * erip;
      eri.compute(&isc,&jsc,&kcp,&lsc);
      erip=eri.getp();

      // Collect terms
      for(size_t kca=0;kca<kcart.size();kca++) {
	int l=kcart[kca].l;
	int m=kcart[kca].m;
	int n=kcart[kca].n;

	int il=l;
	int im=m;
	int in=n;

	double fac=kc[kic].c*sqrt(kc[kic].z);
	if(idx==6) {
	  fac*=sqrt(2*l+1);
	  il++;
	} else if(idx==7) {
	  fac*=sqrt(2*m+1);
	  im++;
	} else if(idx==8) {
	  fac*=sqrt(2*n+1);
	  in++;
	}

	// k index of target integral is
	size_t kdt=getind(l,m,n);
	// k index of source integral is
	size_t kdp=getind(il,im,in);

	// Loop over target integrals
	for(size_t ii=0;ii<Ni;ii++)
	  for(size_t jj=0;jj<Nj;jj++)
	    for(size_t ll=0;ll<Nl;ll++)
	      (*input)[((ii*Nj+jj)*Nk+kdt)*Nl+ll]+=fac*(*erip)[((ii*Nj+jj)*Nkp+kdp)*Nl+ll];
      }

      if(ks_orig->get_am()>0) {
	GaussianShell kcm=GaussianShell(ks_orig->get_am()-1,false,dumcontr);
	kcm.set_center(ks_orig->get_center(),ks_orig->get_center_ind());
	kcm.set_first_ind(ks_orig->get_first_ind());
	kcm.normalize();
	size_t Nkm=kcm.get_Ncart();

	// Evaluate ERI
	eri.compute(&isc,&jsc,&kcm,&lsc);
	erip=eri.getp();

	// Collect terms
	for(size_t kca=0;kca<kcart.size();kca++) {
	  int l=kcart[kca].l;
	  int m=kcart[kca].m;
	  int n=kcart[kca].n;

	  // Skip nonexistent integrals
	  if(idx==6 && l==0)
	    continue;
	  if(idx==7 && m==0)
	    continue;
	  if(idx==8 && n==0)
	    continue;

	  int il=l;
	  int im=m;
	  int in=n;

	  double fac=-kc[kic].c*sqrt(kc[kic].z);
	  if(idx==6) {
	    fac*=2*l/sqrt(2*l-1);
	    il--;
	  } else if(idx==7) {
	    fac*=2*m/sqrt(2*m-1);
	    im--;
	  } else if(idx==8) {
	    fac*=2*n/sqrt(2*n-1);
	    in--;
	  }

	  // k index of target integral is
	  size_t kdt=getind(l,m,n);
	  // k index of source integral is
	  size_t kdm=getind(il,im,in);

	  // Loop over target integrals
	  for(size_t ii=0;ii<Ni;ii++)
	    for(size_t jj=0;jj<Nj;jj++)
	      for(size_t ll=0;ll<Nl;ll++)
		(*input)[((ii*Nj+jj)*Nk+kdt)*Nl+ll]+=fac*(*erip)[((ii*Nj+jj)*Nkm+kdm)*Nl+ll];
	}
      }
    }
  } else if(idx>=9 && idx<12) {
    // Get normalized contraction.
    std::vector<contr_t> lc=lsc.get_contr_normalized();

    // Get cartesians
    std::vector<shellf_t> lcart=lsc.get_cart();

    // Loop over contraction
    for(size_t lic=0;lic<lc.size();lic++) {
      // Dummy contraction
      std::vector<contr_t> dumcontr(1);
      dumcontr[0].c=1.0;
      dumcontr[0].z=lc[lic].z;

      // Form helpers
      GaussianShell lcp(ls_orig->get_am()+1,false,dumcontr);
      lcp.set_center(ls_orig->get_center(),ls_orig->get_center_ind());
      lcp.set_first_ind(ls_orig->get_first_ind());
      lcp.normalize();
      size_t Nlp=lcp.get_Ncart();

      // Evaluate ERI
      const std::vector<double> * erip;
      eri.compute(&isc,&jsc,&ksc,&lcp);
      erip=eri.getp();

      // Collect terms
      for(size_t lca=0;lca<lcart.size();lca++) {
	int l=lcart[lca].l;
	int m=lcart[lca].m;
	int n=lcart[lca].n;

	int il=l;
	int im=m;
	int in=n;

	double fac=lc[lic].c*sqrt(lc[lic].z);
	if(idx==9) {
	  fac*=sqrt(2*l+1);
	  il++;
	} else if(idx==10) {
	  fac*=sqrt(2*m+1);
	  im++;
	} else if(idx==11) {
	  fac*=sqrt(2*n+1);
	  in++;
	}

	// l index of target integral is
	size_t ldt=getind(l,m,n);
	// l index of source integral is
	size_t ldp=getind(il,im,in);

	// Loop over target integrals
	for(size_t ii=0;ii<Ni;ii++)
	  for(size_t jj=0;jj<Nj;jj++)
	    for(size_t kk=0;kk<Nk;kk++)
	      (*input)[((ii*Nj+jj)*Nk+kk)*Nl+ldt]+=fac*(*erip)[((ii*Nj+jj)*Nk+kk)*Nlp+ldp];
      }

      if(ls_orig->get_am()>0) {
	GaussianShell lcm=GaussianShell(ls_orig->get_am()-1,false,dumcontr);
	lcm.set_center(ls_orig->get_center(),ls_orig->get_center_ind());
	lcm.set_first_ind(ls_orig->get_first_ind());
	lcm.normalize();
	size_t Nlm=lcm.get_Ncart();

	// Evaluate ERI
	eri.compute(&isc,&jsc,&ksc,&lcm);
	erip=eri.getp();

	// Collect terms
	for(size_t lca=0;lca<lcart.size();lca++) {
	  int l=lcart[lca].l;
	  int m=lcart[lca].m;
	  int n=lcart[lca].n;

	  // Skip nonexistent integrals
	  if(idx==9 && l==0)
	    continue;
	  if(idx==10 && m==0)
	    continue;
	  if(idx==11 && n==0)
	    continue;

	  int il=l;
	  int im=m;
	  int in=n;

	  double fac=-lc[lic].c*sqrt(lc[lic].z);
	  if(idx==9) {
	    fac*=2*l/sqrt(2*l-1);
	    il--;
	  } else if(idx==10) {
	    fac*=2*m/sqrt(2*m-1);
	    im--;
	  } else if(idx==11) {
	    fac*=2*n/sqrt(2*n-1);
	    in--;
	  }

	  // l index of target integral is
	  size_t ldt=getind(l,m,n);
	  // l index of source integral is
	  size_t ldm=getind(il,im,in);

	  // Loop over target integrals
	  for(size_t ii=0;ii<Ni;ii++)
	    for(size_t jj=0;jj<Nj;jj++)
	      for(size_t kk=0;kk<Nk;kk++)
		(*input)[((ii*Nj+jj)*Nk+kk)*Nl+ldt]+=fac*(*erip)[((ii*Nj+jj)*Nk+kk)*Nlm+ldm];
	}
      }
    }
  }

  // Convert the integrals to the shperical harmonics basis
  spherical_transform(is_orig,js_orig,ks_orig,ls_orig);

  // Return the integrals
  return *input;
}

eri_precursor_t IntegralWorker::compute_precursor(const GaussianShell *is, const GaussianShell *js) {
  // Returned precursor
  eri_precursor_t r;

  // Initialize arrays
  r.AB.zeros(3);

  r.zeta.zeros(is->get_Ncontr(),js->get_Ncontr());
  r.P.zeros(is->get_Ncontr(),js->get_Ncontr(),3);
  r.PA.zeros(is->get_Ncontr(),js->get_Ncontr(),3);
  r.PB.zeros(is->get_Ncontr(),js->get_Ncontr(),3);
  r.S.zeros(is->get_Ncontr(),js->get_Ncontr());

  // Get data
  r.ic=is->get_contr();
  r.jc=js->get_contr();

  coords_t Ac=is->get_center();
  arma::vec A(3);
  A(0)=Ac.x;
  A(1)=Ac.y;
  A(2)=Ac.z;

  coords_t Bc=js->get_center();
  arma::vec B(3);
  B(0)=Bc.x;
  B(1)=Bc.y;
  B(2)=Bc.z;

  // Compute AB
  r.AB=A-B;
  double rabsq=arma::dot(r.AB,r.AB);

  // Compute zeta
  for(size_t i=0;i<r.ic.size();i++)
    for(size_t j=0;j<r.jc.size();j++)
      r.zeta(i,j)=r.ic[i].z+r.jc[j].z;

  // Form P
  for(size_t i=0;i<r.ic.size();i++)
    for(size_t j=0;j<r.jc.size();j++)
      for(int k=0;k<3;k++)
	r.P(i,j,k)=(r.ic[i].z*A(k) + r.jc[j].z*B(k))/r.zeta(i,j);

  // Compute PA and PB
  for(size_t i=0;i<r.ic.size();i++)
    for(size_t j=0;j<r.jc.size();j++)
      for(int k=0;k<3;k++) {
	r.PA(i,j,k)=r.P(i,j,k)-A(k);
	r.PB(i,j,k)=r.P(i,j,k)-B(k);
      }

  // Compute S
  for(size_t i=0;i<r.ic.size();i++)
    for(size_t j=0;j<r.jc.size();j++)
      r.S(i,j)=r.ic[i].c*r.jc[j].c*(M_PI/r.zeta(i,j))*sqrt(M_PI/r.zeta(i,j))*exp(-r.ic[i].z*r.jc[j].z/r.zeta(i,j)*rabsq);

  return r;
}

void IntegralWorker::compute_G(double rho, double T, int nmax) {
  // Evaluate Boys' function
  (void) rho;
#ifdef BOYSNOINTERP
  boysF_arr(nmax,T,Gn);
#else
  BoysTable::eval(nmax,T,Gn);
#endif
}

void ERIWorker::compute_libint_data(const eri_precursor_t & ip, const eri_precursor_t & jp, int mmax) {
  // Store AB and CD
  for(int i=0;i<3;i++) {
    libint.AB[i]=ip.AB(i);
    libint.CD[i]=jp.AB(i);
  }

  size_t ind=0;

  // Two-product centers
  double P[3], Q[3];
  // Four-product center
  double W[3];
  // Distances
  double PQ[3], WP[3], WQ[3];

  // Compute primitive data
  for(size_t p=0;p<ip.ic.size();p++)
    for(size_t q=0;q<ip.jc.size();q++) {
      double zeta=ip.zeta(p,q);

      // Compute overlaps for auxiliary integrals
      double S12=ip.S(p,q);

      for(size_t r=0;r<jp.ic.size();r++)
	for(size_t s=0;s<jp.jc.size();s++) {
	  double eta=jp.zeta(r,s);

	  // Overlap for auxiliary integral
	  double S34=jp.S(r,s);

	  // Reduced exponent
          double rho=zeta*eta/(zeta+eta);

	  P[0]=ip.P(p,q,0);
	  P[1]=ip.P(p,q,1);
	  P[2]=ip.P(p,q,2);

	  Q[0]=jp.P(r,s,0);
	  Q[1]=jp.P(r,s,1);
	  Q[2]=jp.P(r,s,2);

	  W[0]=(zeta*P[0]+eta*Q[0])/(zeta+eta);
	  W[1]=(zeta*P[1]+eta*Q[1])/(zeta+eta);
	  W[2]=(zeta*P[2]+eta*Q[2])/(zeta+eta);

	  PQ[0]=P[0]-Q[0];
	  PQ[1]=P[1]-Q[1];
	  PQ[2]=P[2]-Q[2];

	  WP[0]=W[0]-P[0];
	  WP[1]=W[1]-P[1];
	  WP[2]=W[2]-P[2];

	  WQ[0]=W[0]-Q[0];
	  WQ[1]=W[1]-Q[1];
	  WQ[2]=W[2]-Q[2];

	  double rpqsq=pow(PQ[0],2) + pow(PQ[1],2) + pow(PQ[2],2);

	  // Helper variable
	  prim_data data;

          // Store PA, QC, WP and WQ
          for(int i=0;i<3;i++) {
            data.U[0][i]=ip.PA(p,q,i); // PA
            data.U[2][i]=jp.PA(r,s,i); // QC
            data.U[4][i]=WP[i]; // WP
            data.U[5][i]=WQ[i]; // WQ
          }

	  // Store exponents
	  data.oo2z=0.5/zeta;
	  data.oo2n=0.5/eta;
	  data.oo2zn=0.5/(eta+zeta);
	  data.oo2p=0.5/rho;
	  data.poz=rho/zeta;
	  data.pon=rho/eta;

	  // Kernel prefactor is
	  double prefac=2.0*sqrt(rho/M_PI)*S12*S34;

	  // Compute the kernel
	  compute_G(rho,rho*rpqsq,mmax);

	  // Store the integrals
	  for(int i=0;i<=mmax;i++)
	    data.F[i]=prefac*Gn(i);
	  
	  // We have all necessary data; store quartet.
	  libint.PrimQuartet[ind++]=data;
	}
    }
}

void dERIWorker::compute_libderiv_data(const eri_precursor_t & ip, const eri_precursor_t & jp, int mmax) {
  // Store AB and CD
  for(int i=0;i<3;i++) {
    libderiv.AB[i]=ip.AB(i);
    libderiv.CD[i]=jp.AB(i);
  }

  size_t ind=0;

  // Two-product centers
  double P[3], Q[3];
  // Four-product center
  double W[3];
  // Distances
  double PQ[3], WP[3], WQ[3];


  // Compute primitive data
  for(size_t p=0;p<ip.ic.size();p++)
    for(size_t q=0;q<ip.jc.size();q++) {
      double zeta=ip.zeta(p,q);

      // Compute overlaps for auxiliary integrals
      double S12=ip.S(p,q);

      for(size_t r=0;r<jp.ic.size();r++)
	for(size_t s=0;s<jp.jc.size();s++) {
	  double eta=jp.zeta(r,s);

	  // Overlap for auxiliary integral
	  double S34=jp.S(r,s);

	  // Reduced exponent
          double rho=zeta*eta/(zeta+eta);

	  P[0]=ip.P(p,q,0);
	  P[1]=ip.P(p,q,1);
	  P[2]=ip.P(p,q,2);

	  Q[0]=jp.P(r,s,0);
	  Q[1]=jp.P(r,s,1);
	  Q[2]=jp.P(r,s,2);

	  W[0]=(zeta*P[0]+eta*Q[0])/(zeta+eta);
	  W[1]=(zeta*P[1]+eta*Q[1])/(zeta+eta);
	  W[2]=(zeta*P[2]+eta*Q[2])/(zeta+eta);

	  PQ[0]=P[0]-Q[0];
	  PQ[1]=P[1]-Q[1];
	  PQ[2]=P[2]-Q[2];

	  WP[0]=W[0]-P[0];
	  WP[1]=W[1]-P[1];
	  WP[2]=W[2]-P[2];

	  WQ[0]=W[0]-Q[0];
	  WQ[1]=W[1]-Q[1];
	  WQ[2]=W[2]-Q[2];

	  double rpqsq=pow(PQ[0],2) + pow(PQ[1],2) + pow(PQ[2],2);

	  // Helper variable
	  prim_data data;

          // Store PA, PB, QC, QD, WP and WQ
          for(int i=0;i<3;i++) {
            data.U[0][i]=ip.PA(p,q,i); // PA
	    data.U[1][i]=ip.PB(p,q,i); // PB
            data.U[2][i]=jp.PA(r,s,i); // QC
	    data.U[3][i]=jp.PB(r,s,i); // QD
            data.U[4][i]=WP[i]; // WP
            data.U[5][i]=WQ[i]; // WQ
          }

	  // Store exponents
	  data.twozeta_a=2.0*ip.ic[p].z;
	  data.twozeta_b=2.0*ip.jc[q].z;
	  data.twozeta_c=2.0*jp.ic[r].z;
	  data.twozeta_d=2.0*jp.jc[s].z;

	  data.oo2z=0.5/zeta;
	  data.oo2n=0.5/eta;
	  data.oo2zn=0.5/(eta+zeta);
	  data.oo2p=0.5/rho;
	  data.poz=rho/zeta;
	  data.pon=rho/eta;

	  // Kernel prefactor is
	  double prefac=2.0*sqrt(rho/M_PI)*S12*S34;

	  // Compute the kernel
	  compute_G(rho,rho*rpqsq,mmax+1);

	  // Store auxiliary integrals
	  for(int i=0;i<=mmax+1;i++)
	    data.F[i]=prefac*Gn[i];

	  // We have all necessary data; store quartet.
	  libderiv.PrimQuartet[ind++]=data;
	}
    }
}

void ERIWorker::compute(const GaussianShell *is_orig, const GaussianShell *js_orig, const GaussianShell *ks_orig, const GaussianShell *ls_orig) {
  // Calculate ERIs and transform them to spherical harmonics basis, if necessary.

  // Helpers
  const GaussianShell *is=is_orig;
  const GaussianShell *js=js_orig;
  const GaussianShell *ks=ks_orig;
  const GaussianShell *ls=ls_orig;

  // Did we need to swap the indices?
  bool swap_ij=(is->get_am()<js->get_am());
  bool swap_kl=(ks->get_am()<ls->get_am());
  bool swap_ijkl=(is->get_am()+js->get_am() > ks->get_am() + ls->get_am());

  // Swap shells if necessary
  if(swap_ij)
    std::swap(is,js);
  if(swap_kl)
    std::swap(ks,ls);
  if(swap_ijkl) {
    std::swap(is,ks);
    std::swap(js,ls);
  }

  // Get the cartesian ERIs
  compute_cartesian(is,js,ks,ls);
  // Restore the original order
  reorder(is_orig,js_orig,ks_orig,ls_orig,swap_ij,swap_kl,swap_ijkl);
  // and transform them into the spherical basis
  spherical_transform(is_orig,js_orig,ks_orig,ls_orig);
}

void ERIWorker::compute_debug(const GaussianShell *is, const GaussianShell *js, const GaussianShell *ks, const GaussianShell *ls) {
  // Get the cartesian ERIs
  compute_cartesian_debug(is,js,ks,ls);
  // and transform them into the spherical basis
  spherical_transform(is,js,ks,ls);
}

void ERIWorker::compute_cartesian_debug(const GaussianShell *is, const GaussianShell *js, const GaussianShell *ks, const GaussianShell *ls) {
  std::vector<shellf_t> carti(is->get_cart());
  std::vector<shellf_t> cartj(js->get_cart());
  std::vector<shellf_t> cartk(ks->get_cart());
  std::vector<shellf_t> cartl(ls->get_cart());

  std::vector<contr_t> contri(is->get_contr());
  std::vector<contr_t> contrj(js->get_contr());
  std::vector<contr_t> contrk(ks->get_contr());
  std::vector<contr_t> contrl(ls->get_contr());
  
  coords_t Ri(is->get_center());
  coords_t Rj(js->get_center());
  coords_t Rk(ks->get_center());
  coords_t Rl(ls->get_center());

  input->assign(carti.size()*cartj.size()*cartk.size()*cartl.size(),0.0);

  for(size_t ic=0;ic<carti.size();ic++)
    for(size_t jc=0;jc<cartj.size();jc++)
      for(size_t kc=0;kc<cartk.size();kc++)
	for(size_t lc=0;lc<cartl.size();lc++) {
	  int li(carti[ic].l);
	  int mi(carti[ic].m);
	  int ni(carti[ic].n);
	  double reli(carti[ic].relnorm);
	  
	  int lj(cartj[jc].l);
	  int mj(cartj[jc].m);
	  int nj(cartj[jc].n);
	  double relj(cartj[jc].relnorm);
	  
	  int lk(cartk[kc].l);
	  int mk(cartk[kc].m);
	  int nk(cartk[kc].n);
	  double relk(cartk[kc].relnorm);
	  
	  int ll(cartl[lc].l);
	  int ml(cartl[lc].m);
	  int nl(cartl[lc].n);
	  double rell(cartl[lc].relnorm);

	  double el=0.0;
	  
	  for(size_t xi=0;xi<contri.size();xi++)
	    for(size_t xj=0;xj<contrj.size();xj++)
	      for(size_t xk=0;xk<contrk.size();xk++)
		for(size_t xl=0;xl<contrl.size();xl++) {
		  double zi(contri[xi].z);
		  double zj(contrj[xj].z);
		  double zk(contrk[xk].z);
		  double zl(contrl[xl].z);

		  double ci(contri[xi].c);
		  double cj(contrj[xj].c);
		  double ck(contrk[xk].c);
		  double cl(contrl[xl].c);

		  el+=ci*cj*ck*cl*ERI_int(li,mi,ni,Ri.x,Ri.y,Ri.z,zi,	\
					  lj,mj,nj,Rj.x,Rj.y,Rj.z,zj,	\
					  lk,mk,nk,Rk.x,Rk.y,Rk.z,zk,	\
					  ll,ml,nl,Rl.x,Rl.y,Rl.z,zl);
		}
	  
	  (*input)[((ic*cartj.size()+jc)*cartk.size()+kc)*cartl.size()+lc]=reli*relj*relk*rell*el;
	}
}

std::vector<double> ERIWorker::get() const {
  return *input;
}

const std::vector<double> * ERIWorker::getp() const {
  return input;
}

void dERIWorker::compute(const GaussianShell *is_origv, const GaussianShell *js_origv, const GaussianShell *ks_origv, const GaussianShell *ls_origv) {
  // Calculate derivative ERIs

  is_orig=is_origv;
  js_orig=js_origv;
  ks_orig=ks_origv;
  ls_orig=ls_origv;

  is=is_orig;
  js=js_orig;
  ks=ks_orig;
  ls=ls_orig;

  // Did we need to swap the indices?
  swap_ij=false;
  swap_kl=false;
  swap_ijkl=false;

  // Check order and swap shells if necessary
  if(is->get_am()<js->get_am()) {
    swap_ij=true;
    std::swap(is,js);
  }

  if(ks->get_am()<ls->get_am()) {
    swap_kl=true;
    std::swap(ks,ls);
  }

  if(is->get_am()+js->get_am() > ks->get_am() + ls->get_am()) {
    swap_ijkl=true;
    std::swap(is,ks);
    std::swap(js,ls);
  }

  // Get the cartesian ERIs
  compute_cartesian();
}

void IntegralWorker::reorder(const GaussianShell *is, const GaussianShell *js, const GaussianShell *ks, const GaussianShell *ls, bool swap_ij, bool swap_kl, bool swap_ijkl) {
  // Integrals are now in the (*input) array.
  // Get them in the original order.

  if(!swap_ijkl && !swap_kl && !swap_ij) // 000
    // Already in correct order.
    return;

  // Numbers of functions on each shell
  const size_t Ni=is->get_Ncart();
  const size_t Nj=js->get_Ncart();
  const size_t Nk=ks->get_Ncart();
  const size_t Nl=ls->get_Ncart();

  // Check that we have enough memory
  (*output).resize((*input).size());

  if(!swap_ijkl && !swap_kl && swap_ij) { // 001
    // Need two switch i and j.
    for(size_t ii=0;ii<Ni;ii++)
      for(size_t jj=0;jj<Nj;jj++)
	for(size_t kk=0;kk<Nk;kk++)
	  for(size_t ll=0;ll<Nl;ll++) {
	    size_t iout(((ii*Nj+jj)*Nk+kk)*Nl+ll);
	    size_t iin(((jj*Ni+ii)*Nk+kk)*Nl+ll);
	    (*output)[iout]=(*input)[iin];
	  }

  } else if(!swap_ijkl && swap_kl && !swap_ij) { // 010
    // Need to switch k and l
    
    for(size_t ii=0;ii<Ni;ii++)
      for(size_t jj=0;jj<Nj;jj++)
	for(size_t kk=0;kk<Nk;kk++)
	  for(size_t ll=0;ll<Nl;ll++) {
	    size_t iout(((ii*Nj+jj)*Nk+kk)*Nl+ll);
	    size_t iin(((ii*Nj+jj)*Nl+ll)*Nk+kk);
	    (*output)[iout]=(*input)[iin];
	  }

  } else if(!swap_ijkl && swap_kl && swap_ij) { // 011
    // Switch i <-> j, and k <-> l

    for(size_t ii=0;ii<Ni;ii++)
      for(size_t jj=0;jj<Nj;jj++)
	for(size_t kk=0;kk<Nk;kk++)
	  for(size_t ll=0;ll<Nl;ll++) {
	    size_t iout(((ii*Nj+jj)*Nk+kk)*Nl+ll);
	    size_t iin(((jj*Ni+ii)*Nl+ll)*Nk+kk);
	    (*output)[iout]=(*input)[iin];
	  }

  } else if(swap_ijkl && !swap_kl && !swap_ij) { // 100
    // Switch i <-> k, and j <-> l

    for(size_t ii=0;ii<Ni;ii++)
      for(size_t jj=0;jj<Nj;jj++)
	for(size_t kk=0;kk<Nk;kk++)
	  for(size_t ll=0;ll<Nl;ll++) {
	    size_t iout(((ii*Nj+jj)*Nk+kk)*Nl+ll);
	    size_t iin(((kk*Nl+ll)*Ni+ii)*Nj+jj);
	    (*output)[iout]=(*input)[iin];
	  }

  } else if(swap_ijkl && !swap_kl && swap_ij) { // 101
    // i -> k, j -> l, k -> j, l -> i

    for(size_t ii=0;ii<Ni;ii++)
      for(size_t jj=0;jj<Nj;jj++)
	for(size_t kk=0;kk<Nk;kk++)
	  for(size_t ll=0;ll<Nl;ll++) {
	    size_t iout(((ii*Nj+jj)*Nk+kk)*Nl+ll);
	    size_t iin(((kk*Nl+ll)*Nj+jj)*Ni+ii);
	    (*output)[iout]=(*input)[iin];
	  }

  } else if(swap_ijkl && swap_kl && !swap_ij) { // 110
    // i -> l, j -> k, k -> i, l -> j
    
    for(size_t ii=0;ii<Ni;ii++)
      for(size_t jj=0;jj<Nj;jj++)
	for(size_t kk=0;kk<Nk;kk++)
	  for(size_t ll=0;ll<Nl;ll++) {
	    size_t iout(((ii*Nj+jj)*Nk+kk)*Nl+ll);
	    size_t iin(((ll*Nk+kk)*Ni+ii)*Nj+jj);
	    (*output)[iout]=(*input)[iin];
	  }
    

  } else if(swap_ijkl && swap_kl && swap_ij) { // 111
    // i <-> l, j <-> k
    
    for(size_t ii=0;ii<Ni;ii++)
      for(size_t jj=0;jj<Nj;jj++)
	for(size_t kk=0;kk<Nk;kk++)
	  for(size_t ll=0;ll<Nl;ll++) {
	    size_t iout(((ii*Nj+jj)*Nk+kk)*Nl+ll);
	    size_t iin(((ll*Nk+kk)*Nj+jj)*Ni+ii);
	    (*output)[iout]=(*input)[iin];
	  }
    
  } else
    throw std::logic_error("Should not be here!\n");

  // Swap arrays
  std::swap(input,output);
}

ERIWorker_srlr::ERIWorker_srlr(int maxam, int maxcontr, double w, double a, double b) : ERIWorker(maxam,maxcontr) {
  omega=w;
  alpha=a;
  beta=b;
}

ERIWorker_srlr::~ERIWorker_srlr() {
}

dERIWorker_srlr::dERIWorker_srlr(int maxam, int maxcontr, double w, double a, double b) : dERIWorker(maxam,maxcontr) {
  omega=w;
  alpha=a;
  beta=b;
}

dERIWorker_srlr::~dERIWorker_srlr() {
}

void ERIWorker_srlr::compute_G(double rho, double T, int nmax) {
  // Compute helpers
  double omegasq=omega*omega;
  double rhomegasq=omegasq / (omegasq + rho);

  // Evaluate Boys' functions
#ifdef BOYSNOINTERP
  boysF_arr(nmax,T,bf_long);
  boysF_arr(nmax,T*rhomegasq,bf_short);
#else
  BoysTable::eval(nmax,T,bf_long);
  BoysTable::eval(nmax,T*rhomegasq,bf_short);
#endif
  
  // Store values
  Gn.zeros(nmax+1);

  // [ w^2 / (w^2 + r) ]^(n+0.5)
  double w2_wrN=sqrt(rhomegasq);
  for(int i=0;i<=nmax;i++) {
    Gn(i)=(alpha+beta)*bf_long(i) - beta*w2_wrN*bf_short(i);
    w2_wrN*=rhomegasq;
  }
}

void dERIWorker_srlr::compute_G(double rho, double T, int nmax) {
  // Compute helpers
  double omegasq=omega*omega;
  double rhomegasq=omegasq / (omegasq + rho);

  // Evaluate Boys' functions
#ifdef BOYSNOINTERP
  boysF_arr(nmax,T,bf_long);
  boysF_arr(nmax,T*rhomegasq,bf_short);
#else
  BoysTable::eval(nmax,T,bf_long);
  BoysTable::eval(nmax,T*rhomegasq,bf_short);
#endif

  // Store values
  Gn.zeros(nmax+1);

  // [ w^2 / (w^2 + r) ]^(n+0.5)
  double w2_wrN=sqrt(rhomegasq);
  for(int i=0;i<=nmax;i++) {
    Gn(i)=(alpha+beta)*bf_long(i) - beta*w2_wrN*bf_short(i);
    w2_wrN*=rhomegasq;
  }
}
