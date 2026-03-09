#ifndef GR_H
#define GR_H

#include <NTL/ZZ_pE.h>

using namespace NTL;

// Random element in R*
ZZ_pE randomInvertible();

ZZ_pE indexedElementInExceptionalSet(long index);

Vec<ZZ_pE> getExceptionalSubset(long size);

// Random element in A
ZZ_pE randomInExceptionalSet();

// Random element in A*
ZZ_pE randomNonZeroInExceptionalSet();

ZZ_pE getInverse(ZZ_pE element);

ZZ_pX primitiveIrredPoly(long degree);

// ============================================================
// 新增：判断当前 GR(2^k, d) 是否为 large ring
// large ring: 2^d >= 2^lambda, 即 d >= lambda
// ============================================================
inline bool isLargeRing(long d, long lambda = 128) {
    return d >= lambda;
}

inline long smallRingPackingFactor(long d, long lambda = 128) {
    if (d >= lambda) return 1;
    return (lambda + d - 1) / d;
}

#endif