/* *************************************************************

Copyright 2019 Terry Lyons.

Distributed under the terms of the GNU General Public License,
Version 3. (See accompanying file License.txt)

************************************************************* */
// a debugging tool - SHOW(X) outputs variable name X and its content to a stream (e.g. cout) 
#include "SHOW.h"

// libalgebra functionality
#include "alg_framework.h"

// std:: dependencies for current tests
#include <iostream>
#include <vector>
#include <random>
#include <string>

// the unit test framework
#include <UnitTest++/UnitTest++.h>

// Local frameworks
#include "brown_path_increments.h"
#include "memfile.h"
#include "time_and_details.h"

typedef brown_path_increments<5, 5, 60> pathsetup5560;

TEST_FIXTURE(pathsetup5560, logsignature_versus_cbh)
{
	TEST_DETAILS();
	// collect the signature
	TENSOR sig = signature(increments.begin(), increments.end());

	// collect input for cbh (vector of pointers to Lie increments)
	std::vector<const LIE*> vec_of_ptr_to_lie;

	for (std::vector<LIE>::const_iterator i = increments.cbegin(); i != increments.cend(); i++)
		vec_of_ptr_to_lie.push_back(&(*i));

	// make logsignatures
	LIE logsig1 = maps.t2l(log(sig));
	LIE logsig2 = cbh.full(vec_of_ptr_to_lie);

	// compare logsignatures
	LIE err = logsig1 - logsig2;
	for (auto k : err) {
		CHECK_CLOSE(k.second, 0., 7.0e-16);
	}

	// check dimension of log signature
	CHECK_EQUAL(logsig1.size(), 829);
}

TEST_FIXTURE(pathsetup5560, simple_multiplication)
{
	TEST_DETAILS();
	auto begin = increments.cbegin();
	auto end = increments.cend();
	auto sig1 = signature(begin, begin + (end - begin) / 2);
	auto sig2 = signature(begin + (end - begin) / 2, end);
	auto sig = signature(begin, end);

	TENSOR err = sig - sig1 * sig2;
	for (auto k : err) {
		CHECK_CLOSE(k.second, 0., 2.0e-15);
	}
}

TEST_FIXTURE(pathsetup5560, long_multiplication)
{
	TEST_DETAILS();
	auto begin = increments.cbegin();
	auto end = increments.cend();
	TENSOR sig = signature(begin, end);
	for (auto i = begin; i != end; i++) {
		TENSOR err = sig - signature(begin, i) * signature(i, end);
		for (auto k : err) {
			CHECK_CLOSE(k.second, 0., 2.0e-15);
		}
	}
}

TEST_FIXTURE(pathsetup5560, fine_changes_to_arithmetic_using_memory_mapped_file)
{
	TEST_DETAILS();
	UNITTEST_TIME_CONSTRAINT(5000);
	// the data from the framework
	auto cbegin_framework = increments.cbegin();
	auto cend_framework = increments.cend();
		
	TENSOR sig = signature(cbegin_framework, cend_framework);
	CHECK_compare_with_file(sig, "signature.raw");

	LIE logsig = logsignature(cbegin_framework, cend_framework);
	CHECK_compare_with_file(logsig, "logsignature.raw");
}