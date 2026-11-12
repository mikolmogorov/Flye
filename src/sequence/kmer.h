//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <cstddef>
#include <cstdint>
#if defined(_MSC_VER)
#include <cstdlib>
#endif
#include <functional>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>
#include <deque>

#include "../common/config.h"
#include "../sequence/sequence.h"

static_assert(sizeof(size_t) == 8, "32-bit architectures are not supported");

class Kmer
{
public:
	typedef size_t KmerRepr;

	explicit Kmer(KmerRepr repr=0): _representation(repr) {}

	Kmer(const DnaSequence& dnaString, 
		   size_t start, size_t length):
		_representation(0)
	{
		if (length != Parameters::get().kmerSize)
		{
			throw std::runtime_error("Kmer length inconsistency");
		}

		for (size_t i = start; i < start + length; ++i)	
		{
			_representation <<= 2;
			_representation += dnaString.atRaw(i);
		}
	}

	static KmerRepr bswap64(KmerRepr v)
	{
#if defined(__GNUC__) || defined(__clang__)
		return __builtin_bswap64(v);
#elif defined(_MSC_VER)
		return _byteswap_uint64(v);
#else
		v = ((v >> 8) & 0x00FF00FF00FF00FFULL) | ((v & 0x00FF00FF00FF00FFULL) << 8);
		v = ((v >> 16) & 0x0000FFFF0000FFFFULL) | ((v & 0x0000FFFF0000FFFFULL) << 16);
		return (v >> 32) | (v << 32);
#endif
	}

	Kmer reverseComplement() const
	{
		KmerRepr v = ~_representation;
		v = ((v >> 2) & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
		v = ((v >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
		v = bswap64(v);
		Kmer newKmer;
		newKmer._representation = v >> (64 - (Parameters::get().kmerSize << 1));
		return newKmer;
	}

	bool standardForm()
	{
		Kmer complKmer = this->reverseComplement();
		if (complKmer._representation < _representation)
		{
			_representation = complKmer._representation;
			return true;
		}
		return false;
	}

	void appendRight(DnaSequence::NuclType dnaSymbol)
	{
		_representation <<= 2;
		_representation += dnaSymbol;

		KmerRepr kmerSize = Parameters::get().kmerSize;
		KmerRepr kmerMask = (kmerSize >= 32) ? ~KmerRepr(0) : (((KmerRepr)1 << (kmerSize * 2)) - 1);
		_representation &= kmerMask;
	}

	void appendLeft(DnaSequence::NuclType dnaSymbol)
	{
		_representation >>= 2;

		KmerRepr kmerSize = Parameters::get().kmerSize;
		KmerRepr shift = kmerSize * 2 - 2;
		_representation += dnaSymbol << shift;
	}


	bool operator == (const Kmer& other) const
		{return this->_representation == other._representation;}

	bool operator != (const Kmer& other) const
		{return !(*this == other);}

	size_t hash() const
	{
		size_t x = _representation;
		size_t z = (x += 0x9E3779B97F4A7C15ULL);
		z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
		z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
		return z ^ (z >> 31);
	}

	bool operator< (const Kmer& other)
	{
		return _representation < other._representation;
	}

	size_t numRepr() {return _representation;}

private:
	KmerRepr _representation;
};

namespace std
{
	template <>
	struct hash<Kmer>
	{
		std::size_t operator()(const Kmer& kmer) const
		{
			return kmer.hash();
		}
	};
}

struct KmerPosition
{
	KmerPosition(Kmer kmer, int32_t position):
		kmer(kmer), position(position) {}
	Kmer kmer;
	int32_t position;
};

class KmerIterator
{
public:
    typedef std::forward_iterator_tag iterator_category;

	KmerIterator(const DnaSequence* readSeq, size_t position):
		_readSeq(readSeq),
		_position(position)
	{
		if (position != readSeq->length() - Parameters::get().kmerSize)
		{
			//_kmer = Kmer(readSeq->substr(0, Parameters::get().kmerSize));
			_kmer = Kmer(*readSeq, 0, Parameters::get().kmerSize);
		}
	}

	bool operator==(const KmerIterator& other) const
	{
		return _readSeq == other._readSeq && _position == other._position;
	}

	bool operator!=(const KmerIterator& other) const
	{
		return !(*this == other);
	}

	KmerPosition operator*() const
	{
		return KmerPosition(_kmer, _position);
	}

	KmerIterator& operator++()
	{
		size_t appendPos = _position + Parameters::get().kmerSize;
		_kmer.appendRight(_readSeq->atRaw(appendPos));
		++_position;
		return *this;
	}

protected:
	const DnaSequence* _readSeq;
	size_t 	_position;
	Kmer 	_kmer;
};


class IterKmers
{
public:
	IterKmers(const DnaSequence& sequence, size_t start = 0,
			  size_t length = std::string::npos):
		_sequence(sequence), _start(start), _length(length)
	{}

	KmerIterator begin()
	{
		if (_sequence.length() < Parameters::get().kmerSize + _start)
			return this->end();

		return KmerIterator(&_sequence, _start);
	}

	KmerIterator end()
	{
		size_t end = _length == std::string::npos ?
						_sequence.length() : _length + _start;
		return KmerIterator(&_sequence, end - Parameters::get().kmerSize);
	}

private:
	const DnaSequence& _sequence;
	const size_t _start;
	const size_t _length;
};

inline std::vector<KmerPosition> yieldMinimizers(const DnaSequence& sequence, int window)
{
	if (window < 1) throw std::runtime_error("wrong minimizer length");

	struct KmerAndHash
	{
		KmerPosition kp;
		size_t hash;
	};
	thread_local std::deque<KmerAndHash> miniQueue;
	miniQueue.clear();

	std::vector<KmerPosition> minimizers;
	const size_t expectedSize = sequence.length() / window * 2;
	minimizers.reserve(1.5 * expectedSize);

	if (window == 1)
	{
		for (auto kmerPos : IterKmers(sequence))
		{
			minimizers.push_back(kmerPos);
		}
		return minimizers;
	}

	for (auto kmerPos : IterKmers(sequence))
	{
		auto stdKmer = kmerPos.kmer;
		stdKmer.standardForm();
		size_t curHash = stdKmer.hash();
		
		while (!miniQueue.empty() && miniQueue.back().hash > curHash)
		{
			miniQueue.pop_back();
		}
		miniQueue.push_back({kmerPos, curHash});
		if (miniQueue.front().kp.position <= kmerPos.position - window)
		{
			while (miniQueue.front().kp.position <= kmerPos.position - window)
			{
				miniQueue.pop_front();
			}
			while (miniQueue.size() >= 2 && miniQueue[0].hash == miniQueue[1].hash)
			{
				miniQueue.pop_front();
			}
		}
		if (minimizers.empty() || minimizers.back().position != 
								  miniQueue.front().kp.position)
		{
			minimizers.push_back(miniQueue.front().kp);
		}
	}

	//Logger::get().debug() << _seqContainer.seqLen(seqId) << " " << minimizers.size();
	return minimizers;
}
