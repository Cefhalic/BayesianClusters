#pragma once

/* ===== C++ ===== */
#include <vector>

enum sort_direction { up, down };

template < sort_direction aDir , typename T >
void BitonicSort( const typename std::vector<T>::iterator& aDataStart, const typename std::vector<T>::iterator& aDataEnd );

template < sort_direction aDir , typename T >
void BitonicMerge( const typename std::vector<T>::iterator& aDataStart, const typename std::vector<T>::iterator& aDataEnd);


template < sort_direction aDir , typename T >
void BitonicSort( const typename std::vector<T>::iterator& aDataStart, const typename std::vector<T>::iterator& aDataEnd) {
  uint32_t lSize(aDataEnd - aDataStart);
  if (lSize > 1) {
    typename std::vector<T>::iterator lMidpoint(aDataStart + (lSize >> 1));
    if (aDir == down) {
      BitonicSort<up,T>(aDataStart, lMidpoint);
      BitonicSort<down,T>( lMidpoint, aDataEnd);
    } else {
      BitonicSort<down,T>(aDataStart, lMidpoint);
      BitonicSort<up,T>( lMidpoint, aDataEnd);
    }
    BitonicMerge<aDir,T>(aDataStart, aDataEnd);
  }
}

template < sort_direction aDir , typename T >
void BitonicMerge( const typename std::vector<T>::iterator& aDataStart, const typename std::vector<T>::iterator& aDataEnd) {
  uint32_t lSize(aDataEnd - aDataStart);
  if (lSize > 1) {
    uint32_t lPower2(1);
    while (lPower2 < lSize) lPower2 <<= 1;

    typename std::vector<T>::iterator lMidpoint(aDataStart + (lPower2 >> 1)) , lFirst(aDataStart) , lSecond(lMidpoint);

    for (; lSecond != aDataEnd; ++lFirst, ++lSecond) {
      if (((*lSecond) < (*lFirst)) == (aDir == up)) {
        std::swap(*lFirst, *lSecond);
      }
    }

    BitonicMerge<aDir,T>( aDataStart, lMidpoint);
    BitonicMerge<aDir,T>( lMidpoint, aDataEnd);
  }
}
