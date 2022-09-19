#include "NtuAnalysis/Common/interface/NtuData.h"

template <>
NtuData<float>::NtuData() {
  numberType = "F";
}


template <>
NtuData<double>::NtuData() {
  numberType = "D";
}


