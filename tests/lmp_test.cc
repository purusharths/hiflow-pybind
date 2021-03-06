// Copyright (C) 2011-2021 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for
// more details.
//
// You should have received a copy of the European Union Public Licence (EUPL)
// v1.2 along with HiFlow3.  If not, see
// <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

#define BOOST_TEST_MODULE lmp

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "hiflow.h"
#include <iostream>
#include <stdexcept>

using namespace hiflow::la;
using namespace hiflow;

void SimpleVector_CopyConstructor() {
  CPUsimple_lVector< double > v1(3, "v1");
  v1.add_value(0, 0.0);
  v1.add_value(1, 0.0);
  v1.add_value(2, 0.0);

  // No death test available using CppUnit++. Otherwise following check would be
  // possible:
  // CPUsimple_lVector<double> v2 = v1;
  // v1.~CPUsimple_lVector<double>();
  // CHECK_EXIT(v2.~CPUsimple_lVector<double>(), -1);
};

void CompareVectorMultSimpleOpenMP() {
  CPUsimple_lVector< double > svi(3, "svi");
  svi.add_value(0, 1.0);
  svi.add_value(1, -2.0);
  svi.add_value(2, 0.0);

  CPUsimple_CSR_lMatrix< double > sm(7, 3, 3, "sm");
  std::vector< int > row, col;
  std::vector< double > val;

  // coo data
  row.push_back(0), col.push_back(0), val.push_back(2.0);
  row.push_back(0), col.push_back(1), val.push_back(-1.0);
  row.push_back(1), col.push_back(0), val.push_back(-1.0);
  row.push_back(1), col.push_back(1), val.push_back(2.0);
  row.push_back(1), col.push_back(2), val.push_back(-1.0);
  row.push_back(2), col.push_back(1), val.push_back(-1.0);
  row.push_back(2), col.push_back(2), val.push_back(2.0);

  sm.TransformFromCOO(vec2ptr(row), vec2ptr(col), vec2ptr(val), 3, 3,
                      val.size());

  CPUsimple_lVector< double > svo(3, "svo");
  sm.VectorMult(svi, &svo);

#ifdef WITH_OPENMP
  CPUopenmp_lVector< double > ovi(3, "ovi");
  ovi.CopyFrom(svi);

  CPUopenmp_CSR_lMatrix< double > om;
  om.CopyStructureFrom(sm);
  om.CopyFrom(sm);
  om.set_num_threads(1);

  CPUopenmp_lVector< double > ovo(3, "ovo");
  om.VectorMult(ovi, &ovo);
#endif

  // expected, actual, tolerance
  BOOST_TEST(4.0 - svo.buffer[0] == 0.0);
  BOOST_TEST(-5.0 - svo.buffer[1] == 0.0);
  BOOST_TEST(2.0 - svo.buffer[2] == 0.0);

#ifdef WITH_OPENMP
  CHECK_TEST(4.0 - ovo.buffer[0] == 0.0);
  CHECK_TEST(-5.0- ovo.buffer[1] == 0.0);
  CHECK_TEST(2.0 - ovo.buffer[2] == 0.0);
#endif
};

BOOST_AUTO_TEST_CASE(lmp, *utf::tolerance(1.0e-12)) {

  SimpleVector_CopyConstructor();

  CompareVectorMultSimpleOpenMP();

}
