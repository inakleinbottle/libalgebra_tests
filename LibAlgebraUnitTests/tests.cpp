// libalgebra functionality
#include "libalgebra/alg_types.h"

// boost dependencies for current tests
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp> // exists
#include <boost/filesystem/fstream.hpp>

// std:: dependencies for current tests
#include <iostream>
#include <vector>
#include <random>

// the unit test framework
#include <UnitTest++/UnitTest++.h>

// a debugging tool - SHOW(X) outputs variable name X and its content to a stream (e.g. cout) 
#include "SHOW.h"

// the determining template variables
const unsigned DEPTH = 5;
const unsigned ALPHABET_SIZE = 5;
const unsigned STEPS = 60;

/// brown_path - a class exposing a fixed Brownian path of dimension "ALPHABET_SIZE" 
/// the path is evenly sampled at times 0, 1/Steps, 2/Steps, ..., 1.0 within the 
/// unit interval. The path is presented as a one dimensional stream of doubles [begin(),end())
struct brown_path
{
// state
	const size_t steps;
	const size_t width;	
	std::vector<double> path;

// accessors
	// dimensions
	size_t Steps() const { return steps; }
	size_t Width() const { return width; }

	// the path data
	const double* begin() const { return &(*path.begin()); }
	const double* end() const { return begin() + path.size(); }

// constructor
	brown_path() : 
		steps(STEPS),
		width(ALPHABET_SIZE),
		path((Steps() + 1) * Width(), 0.)
	{
	// set up random number generation
		// seeded generator	
		unsigned int seed = (const unsigned int&)0x6d35f0e5b8f6c603;//std::random_device seed; unsigned int seed = seed();
		std::mt19937 generator;
		generator.seed(seed);

		// rng
		std::normal_distribution<double> distribution(0., 1. / sqrt(Steps()));//distribution(mean, std deviation)

		// create random path with dimension width and steps increments, and so a (steps+1) x width C matrix. 
		for (int i = 0; i < Steps(); ++i)
			for (int j = 0; j < Width(); ++j) {
				double increment = distribution(generator);
				path[(std::size_t(i) + 1) * Width() + j] = path[std::size_t(i) * width + j] + increment;
			}
	}

// destructor
	~brown_path() {}
};

/// brown_path_increments - This class is intended to be derived from in unit tests.
/// It exposes a fixed Brownian path of dimension "ALPHABET_SIZE", and the path is 
/// presented via the LIE representation of the (truncated to level "DEPTH") of the 
/// path increments over the intervals 
///         [0, 1/Steps), [1/Steps, 2/Steps), ..., [(Steps-1)/Steps, 1.0)
/// within the unit interval. The path is presented as a stream of size "STEPS" of LIE elements 
/// in a range [increments.begin(),increments.end())
///
/// In addition to the LIE representation of the path increments it provides the basic tensor objects and the maps 
/// as well as helper functions that can calculate signatures and log signatures of the path over subintervals [b,e)
///
struct brown_path_increments : private brown_path
{
// types
	typedef alg_types<DEPTH, ALPHABET_SIZE, DPReal>  ALG;
	typedef typename ALG::TENSOR TENSOR;
	typedef typename ALG::LIE LIE;
	typedef typename ALG::MAPS MAPS;
	typedef typename ALG::CBH CBH;
	typedef typename ALG::S S;

// state
	mutable MAPS maps; // has cached state
	mutable CBH cbh; // has cached state
	using brown_path::Steps;
	using brown_path::Width;
	std::vector<LIE> increments;

// helper functions
	/// takes flat vector to increment of LIE type
	template<class const_scalar_iterator>
	LIE LieDifference(const_scalar_iterator i) const
	{
		LIE increment;
		for (int j = 0; j < Width(); j++)
			// Hall basis elements start at index 1 with zero as a reserved parent index		
			increment += LIE(j + 1, *(i + j + Width()) - *(i + j));
		return increment;
	}

	/// brown_path_increments inherits its constructor
	template<typename ...Ts>
	brown_path_increments(Ts...) 
		: brown_path(Ts...) {
		// convert sampled path to LIE increments
		for (auto i = path.cbegin(); (i + Width()) != path.cend(); i += Width()) {
			increments.emplace_back(LieDifference(i));
		}
		// remove path data
		brown_path::path.clear();
	}

	/// computes a signature from an iterable sequence of lie elements
	template<class ITERATOR_T>
	TENSOR signature(ITERATOR_T begin, ITERATOR_T end) const
	{
		TENSOR signature(S(1));
		for (ITERATOR_T i = begin; i != end; i++)
			signature *= exp(maps.l2t(*i));
		return signature;
	}

	/// computes the logsignature from the signature
	template<class ITERATOR_T>
	LIE logsignature(ITERATOR_T begin, ITERATOR_T end) const
	{
		TENSOR signature = signature(begin, end);
		return maps.t2l(log(signature));
	}
};

TEST_FIXTURE(brown_path_increments, logsignature_versus_cbh)
{
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

TEST_FIXTURE(brown_path_increments, simple_multiplication)
{
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

TEST_FIXTURE(brown_path_increments, long_multiplication)
{
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

TEST_FIXTURE(brown_path_increments, fine_changes_to_artithmetic_using_memory_mapped_file)
{
	auto begin = increments.cbegin();
	auto end = increments.cend();
	std::string filename("signature.raw");
	int numberOfElements = 3906;
	typedef std::pair<TENSOR::KEY, TENSOR::SCALAR> value_type;
	int numberOfBytes = numberOfElements * sizeof(value_type);
	int count = 10;
	boost::filesystem::path p(filename.c_str());
	while (!boost::filesystem::exists(p) && (count-- > 0)) {
		{
			boost::filesystem::ofstream file;
			file.open(p);
		}// file closed
		boost::filesystem::resize_file(p, 62496);
		// create file
		boost::iostreams::mapped_file_sink file;
		file.open(p);
		if (file.is_open()) {
			value_type* data_begin = (value_type*)(file.data());
			value_type* data_end = data_begin + numberOfElements;
			TENSOR sig = signature(begin, end);
			std::copy(sig.begin(), sig.end(), data_begin);

			double c = data_begin[0].second;
			std::sort(data_begin, data_end, [](const value_type lhs, const value_type rhs) {return lhs.first < rhs.first; });
			file.close();
		}
	}
	if (count == 0) throw;

	// compare the calculated with the stored data
	{
		TENSOR sig;
		boost::iostreams::mapped_file_source file;
		file.open("signature.raw", numberOfBytes);
		if (file.is_open()) {
			const value_type* data_begin = (value_type*)file.data();
			const value_type* data_end = data_begin + numberOfElements;
			for (auto a = data_begin; a < data_end; a++)
				sig[a->first] = a->second;
			file.close();
		}
		TENSOR sig1 = signature(begin, end);
		TENSOR err = sig1 - sig;
		//SHOW(err);
		CHECK_EQUAL(TENSOR(), err);
	}
}
