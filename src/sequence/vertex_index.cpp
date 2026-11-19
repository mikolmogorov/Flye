//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <cmath>
#include <cstdlib>
#include <thread>
#include <new>
#include <sys/mman.h>

#include "vertex_index.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "../common/config.h"
#include "../common/memory_info.h"

//FLYE_FREQ_DIGEST=1: sum a hash of every (k-mer, frequency) pair returned by
//KmerCounter::getFreq(). The sum ignores call order, so two runs that return
//the same frequencies for the same calls print the same digest at clear().
namespace
{
	bool freqDigestOn = false;
	std::atomic<uint64_t> freqDigest(0);
	std::atomic<uint64_t> freqCalls(0);

	uint64_t mixFreq(uint64_t kmer, uint64_t freq)
	{
		uint64_t z = kmer * 0x9E3779B97F4A7C15ULL + freq + 0x632BE59BD9B4E019ULL;
		z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
		z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
		return z ^ (z >> 31);
	}

	size_t recordFreq(Kmer kmer, size_t freq)
	{
		if (freqDigestOn)
		{
			freqDigest.fetch_add(mixFreq(kmer.numRepr(), freq), std::memory_order_relaxed);
			freqCalls.fetch_add(1, std::memory_order_relaxed);
		}
		return freq;
	}
}


void VertexIndex::countKmers()
{
	_kmerCounter.count(/*use flat counter*/ true);
}


void VertexIndex::buildIndexUnevenCoverage(int globalMinFreq, float selectRate,
										   int tandemFreq)
{
	//this->setRepeatCutoff(globalMinFreq);

	//_solidMultiplier = 1;

	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}

	//first, count the number of k-mers that will be actually stored in the index
	_kmerIndex.reserve(_kmerCounter.getKmerNum() / 10);
	if (_outputProgress) Logger::get().info() << "Filling index table (1/2)";
	std::function<void(const FastaRecord::Id&)> initializeIndex = 
	[this, globalMinFreq, selectRate, tandemFreq] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		auto topKmers = this->yieldFrequentKmers(readId, selectRate, tandemFreq);
		for (auto kmerFreq : topKmers)
		{
			kmerFreq.kmer.standardForm();

			if (kmerFreq.freq < (size_t)globalMinFreq) continue;

			ReadVector defVec((uint32_t)1, (uint32_t)0);
			_kmerIndex.upsert(kmerFreq.kmer, 
							  [](ReadVector& rv){++rv.capacity;}, defVec);
		}
	};
	processInParallel(allReads, initializeIndex, 
					  Parameters::get().numThreads, _outputProgress);
	
	this->filterFrequentKmers(globalMinFreq, (float)Config::get("repeat_kmer_rate"));
	this->allocateIndexMemory();

	if (_outputProgress) Logger::get().info() << "Filling index table (2/2)";
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, globalMinFreq, selectRate, tandemFreq] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		auto topKmers = this->yieldFrequentKmers(readId, selectRate, tandemFreq);
		for (auto kmerFreq : topKmers)
		{
			if (kmerFreq.freq < (size_t)globalMinFreq ||
				kmerFreq.freq > _repetitiveFrequency) continue;

			KmerPosition kmerPos(kmerFreq.kmer, kmerFreq.position);
			FastaRecord::Id targetRead = readId;
			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
				targetRead = targetRead.rc();
			}

			//will not trigger update for k-mer not in the index
			_kmerIndex.update_fn(kmerPos.kmer, 
				[targetRead, &kmerPos, this](ReadVector& rv)
				{
					if (rv.size == rv.capacity) 
					{
						Logger::get().warning() << "Index size mismatch " << rv.capacity;
						return;
					}
					size_t globPos = _seqContainer
							.globalPosition(targetRead, kmerPos.position);
					rv.data[rv.size].set(globPos);
					++rv.size;
				});
		}
	};
	processInParallel(allReads, indexUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	_kmerCounter.clear();

	Logger::get().debug() << "Sorting k-mer index";
	for (const auto& kmerVec : _kmerIndex.lock_table())
	{
		std::sort(kmerVec.second.data, kmerVec.second.data + kmerVec.second.size,
				  [](const IndexChunk& p1, const IndexChunk& p2)
				  	{return p1.get() < p2.get();});
	}
	
	size_t totalEntries = 0;
	for (const auto& kmerRec : _kmerIndex.lock_table())
	{
		totalEntries += kmerRec.second.size;
	}
	Logger::get().debug() << "Selected k-mers: " << _kmerIndex.size();
	Logger::get().debug() << "Index size: " << totalEntries;
	Logger::get().debug() << "Mean k-mer index frequency: " 
		<< (float)totalEntries / _kmerIndex.size();
}

namespace
{
	template <class T>
	size_t getFreq(T& histIter)
		{return histIter->first;};

	template <class T>
	size_t getCount(T& histIter)
		{return histIter->second;};

}

void VertexIndex::filterFrequentKmers(int minCoverage, float rate)
{
	size_t totalKmers = 0;
	size_t uniqueKmers = 0;
	for (const auto& kmer : _kmerIndex.lock_table())
	{
		if (kmer.second.capacity >= (size_t)minCoverage)
		{
			totalKmers += kmer.second.capacity;
			uniqueKmers += 1;
		}
	}
	float meanFrequency = (float)totalKmers / (uniqueKmers + 1);
	_repetitiveFrequency = rate * meanFrequency;
	
	size_t repetitiveKmers = 0;
	for (const auto& kmer : _kmerIndex.lock_table())
	{
		if (kmer.second.capacity > _repetitiveFrequency)
		{
			//++repetitiveKmers;
			repetitiveKmers += kmer.second.capacity;
			_repetitiveKmers.insert(kmer.first, true);
		}
	}

	for (const auto& kmer : _repetitiveKmers.lock_table())
	{
		_kmerIndex.erase(kmer.first);
	}

	float filteredRate = (float)repetitiveKmers / totalKmers;
	Logger::get().debug() << "Mean k-mer frequency: " 
						  << meanFrequency;
	Logger::get().debug() << "Repetitive k-mer frequency: " 
						  << _repetitiveFrequency;
	Logger::get().debug() << "Filtered " << repetitiveKmers 
						  << " repetitive k-mers (" <<
						  filteredRate << ")";
}

std::vector<VertexIndex::KmerFreq>
	VertexIndex::yieldFrequentKmers(const FastaRecord::Id& seqId,
									float selectRate, int tandemFreq)
{
	struct FlatSlot
	{
		Kmer key;
		uint32_t count;
		uint32_t epoch;
	};
	thread_local std::vector<FlatSlot> flatTable;
	thread_local uint32_t curEpoch = 1;

	const DnaSequence& seq = _seqContainer.getSeq(seqId);
	const size_t seqLen = seq.length();
	if (seqLen < Parameters::get().kmerSize) return {};
	const size_t numKmers = seqLen - Parameters::get().kmerSize + 1;

	size_t tableMask = 0;
	if (tandemFreq > 0)
	{
		size_t cap = 1024;
		while (cap < numKmers * 2) cap <<= 1;
		if (flatTable.size() < cap)
		{
			flatTable.assign(cap, FlatSlot{Kmer(0), 0, 0});
		}
		tableMask = flatTable.size() - 1;
		if (++curEpoch == 0)
		{
			for (auto& s : flatTable) s.epoch = 0;
			curEpoch = 1;
		}
	}

	std::vector<KmerFreq> topKmers;
	topKmers.reserve(numKmers);

	//Collect the read's k-mers first, so the random flat-counter lookups can
	//be prefetched ahead and several cache misses stay in flight at once.
	//getFreq() and the per-read table still see the k-mers in the original order.
	thread_local std::vector<KmerPosition> readKmers;
	thread_local std::vector<Kmer> stdKmers;
	readKmers.clear();
	stdKmers.clear();
	for (const auto& kmerPos : IterKmers(seq))
	{
		auto stdKmer = kmerPos.kmer;
		stdKmer.standardForm();
		readKmers.push_back(kmerPos);
		stdKmers.push_back(stdKmer);
	}

	const size_t PREFETCH_AHEAD = 16;
	const size_t numRead = stdKmers.size();
	for (size_t i = 0; i < std::min(PREFETCH_AHEAD, numRead); ++i)
	{
		_kmerCounter.prefetch(stdKmers[i]);
	}
	for (size_t i = 0; i < numRead; ++i)
	{
		if (i + PREFETCH_AHEAD < numRead) _kmerCounter.prefetch(stdKmers[i + PREFETCH_AHEAD]);
		const Kmer& stdKmer = stdKmers[i];
		size_t freq = _kmerCounter.getFreq(stdKmer);

		if (tandemFreq > 0)
		{
			size_t idx = stdKmer.hash() & tableMask;
			while (flatTable[idx].epoch == curEpoch && flatTable[idx].key != stdKmer)
			{
				idx = (idx + 1) & tableMask;
			}
			if (flatTable[idx].epoch != curEpoch)
			{
				flatTable[idx] = {stdKmer, 1, curEpoch};
			}
			else
			{
				++flatTable[idx].count;
			}
		}
		topKmers.push_back({readKmers[i].kmer, readKmers[i].position, freq});
	}

	if (topKmers.empty()) return {};
	const size_t maxKmers = std::min<size_t>(selectRate * topKmers.size(), topKmers.size() - 1);
	std::nth_element(topKmers.begin(), topKmers.begin() + maxKmers, topKmers.end(),
					 [](const KmerFreq& k1, const KmerFreq& k2)
					 {return k1.freq > k2.freq;});
	const size_t minFreq = topKmers[maxKmers].freq;

	topKmers.erase(std::remove_if(topKmers.begin(), topKmers.end(),
						[minFreq, tandemFreq, tableMask](KmerFreq kf)
						{
							if (kf.freq < minFreq) return true;
							if (tandemFreq > 0)
							{
								kf.kmer.standardForm();
								size_t idx = kf.kmer.hash() & tableMask;
								while (flatTable[idx].key != kf.kmer)
								{
									idx = (idx + 1) & tableMask;
								}
								return flatTable[idx].count > (uint32_t)tandemFreq;
							}
							return false;
						}),
				   topKmers.end());

	return topKmers;
}



void VertexIndex::allocateIndexMemory()
{
	_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
	size_t chunkOffset = 0;
	//Important: since packed structures are apparently not thread-safe,
	//make sure that adjacent k-mer index arrays (that are accessed in parallel)
	//do not overlap within 8-byte window
	const size_t PADDING = 1;
	for (auto& kmer : _kmerIndex.lock_table())
	{
		if (MEM_CHUNK < kmer.second.capacity + PADDING) 
		{
			throw std::runtime_error("k-mer is too frequent");
		}
		if (MEM_CHUNK - chunkOffset < kmer.second.capacity + PADDING)
		{
			_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
			chunkOffset = 0;
		}
		kmer.second.data = _memoryChunks.back() + chunkOffset;
		chunkOffset += kmer.second.capacity + PADDING;
	}

	//Logger::get().debug() << "Total chunks " << _memoryChunks.size()
	//	<< " wasted space: " << wasted;
}

void VertexIndex::buildIndexMinimizers(int minCoverage, int wndLen)
{
	if (_outputProgress) Logger::get().info() << "Building minimizer index";

	std::vector<FastaRecord::Id> allReads;
	size_t totalLen = 0;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
		if (seq.id.strand()) totalLen += seq.sequence.length();
	}

	_kmerIndex.reserve(1000000);
	if (_outputProgress) Logger::get().info() << "Pre-calculating index storage";
	std::function<void(const FastaRecord::Id&)> initializeIndex = 
	[this, minCoverage, wndLen] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		auto minimizers = yieldMinimizers(_seqContainer.getSeq(readId), wndLen);
		for (auto kmerPos : minimizers)
		{
			auto stdKmer = kmerPos.kmer;
			stdKmer.standardForm();
			ReadVector defVec((uint32_t)1, (uint32_t)0);
			_kmerIndex.upsert(stdKmer, 
							  [](ReadVector& rv){++rv.capacity;}, defVec);
		}
	};
	processInParallel(allReads, initializeIndex, 
					  Parameters::get().numThreads, _outputProgress);

	this->filterFrequentKmers(minCoverage, (float)Config::get("repeat_kmer_rate"));
	this->allocateIndexMemory();
	
	if (_outputProgress) Logger::get().info() << "Filling index";
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, minCoverage, wndLen] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;
		auto minimizers = yieldMinimizers(_seqContainer.getSeq(readId), wndLen);
		for (auto kmerPos : minimizers)
		{
			FastaRecord::Id targetRead = readId;
			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
				targetRead = targetRead.rc();
			}

			if (_repetitiveKmers.contains(kmerPos.kmer)) continue;

			_kmerIndex.update_fn(kmerPos.kmer, 
				[targetRead, &kmerPos, this](ReadVector& rv)
				{
					if (rv.size == rv.capacity) 
					{
						Logger::get().warning() << "Index size mismatch " << rv.capacity;
						return;
					}
					size_t globPos = _seqContainer
							.globalPosition(targetRead, kmerPos.position);
					rv.data[rv.size].set(globPos);
					++rv.size;
				});
		}
	};
	processInParallel(allReads, indexUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	Logger::get().debug() << "Sorting k-mer index";
	for (const auto& kmerVec : _kmerIndex.lock_table())
	{
		std::sort(kmerVec.second.data, kmerVec.second.data + kmerVec.second.size,
				  [](const IndexChunk& p1, const IndexChunk& p2)
				  	{return p1.get() < p2.get();});
	}
	
	size_t totalEntries = 0;
	for (const auto& kmerRec : _kmerIndex.lock_table())
	{
		totalEntries += kmerRec.second.size;
	}
	Logger::get().debug() << "Selected k-mers: " << _kmerIndex.size();
	Logger::get().debug() << "K-mer index size: " << totalEntries;
	Logger::get().debug() << "Mean k-mer frequency: " 
		<< (float)totalEntries / _kmerIndex.size();

	float minimizerRate = (float)totalLen / totalEntries;
	Logger::get().debug() << "Minimizer rate: " << minimizerRate;
	_sampleRate = minimizerRate;
}


void VertexIndex::clear()
{
	for (auto& chunk : _memoryChunks) delete[] chunk;
	_memoryChunks.clear();

	_kmerIndex.clear();
	_kmerIndex.reserve(0);

	_kmerCounter.clear();
	//_kmerCounts.reserve(0);
}


void KmerCounter::count(bool useFlatCounter)
{
	//Logger::get().debug() << "Before counter: " 
	//	<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	if (useFlatCounter && Parameters::get().kmerSize > 17)
	{
		throw std::runtime_error("Can't use flat counter for k-mer size > 17");
	}
	_useFlatCounter = useFlatCounter;

	//flat array for all possible k-mers, 4 bits for each
	//in case of k=17, takes 8Gb
	static const size_t COUNTER_LEN = std::pow(4, Parameters::get().kmerSize) / 2;
	if (useFlatCounter)
	{
		//2 MiB-aligned and marked for transparent huge pages: lookups into this
		//8 GB array are random, and on 4 KiB pages nearly every one also needs a
		//page-table walk. Freed with free() in clear() and ~KmerCounter().
		const size_t HUGE_PAGE = 2 * 1024 * 1024;
		const size_t allocLen = (COUNTER_LEN + HUGE_PAGE - 1) / HUGE_PAGE * HUGE_PAGE;
		void* flatMem = nullptr;
		if (posix_memalign(&flatMem, HUGE_PAGE, allocLen) != 0) throw std::bad_alloc();
		madvise(flatMem, allocLen, MADV_HUGEPAGE);
		std::memset(flatMem, 0, allocLen);
		_flatCounter = static_cast<std::atomic<uint8_t>*>(flatMem);
	}
 
	if (_outputProgress) Logger::get().info() << "Counting k-mers:";
	std::function<void(const FastaRecord::Id&)> readUpdate = 
	[this] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;
		
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			kmerPos.kmer.standardForm();
			bool addOne = true;
			if (_useFlatCounter)
			{
				size_t arrayPos = kmerPos.kmer.numRepr() / 2;
				bool highBits = kmerPos.kmer.numRepr() % 2;

				while (true)
				{
					uint8_t expected = _flatCounter[arrayPos]; 
					uint8_t count = highBits ? (expected >> 4) : (expected & 15);
					if (count == 15)
					{
						break;
					}

					uint8_t updated = highBits ? (expected + 16) : (expected + 1);
					if (_flatCounter[arrayPos].compare_exchange_weak(expected,  updated))
					{
						if (count == 0) ++_numKmers;
						addOne = false; //not saturated yet, don't update hash counter
						break;
					}
				}
			}

			if (addOne)
			{
				_hashCounter.upsert(kmerPos.kmer, [](size_t& num){++num;}, 1);
			}
		}
	};
	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}
	processInParallel(allReads, readUpdate, Parameters::get().numThreads, _outputProgress);

	Logger::get().debug() << "Updating k-mer histogram";
	if (_useFlatCounter)
	{
		const size_t numChunks = 256;
		std::vector<size_t> chunks(numChunks);
		for (size_t i = 0; i < numChunks; ++i) chunks[i] = i;

		std::mutex histMutex;
		size_t globalSmallHist[16] = {0};
		const uint8_t* rawFlat = reinterpret_cast<const uint8_t*>(_flatCounter);

		std::function<void(const size_t&)> scanChunk =
		[rawFlat, numChunks, &histMutex, &globalSmallHist](const size_t& chunkId)
		{
			const size_t startByte = (COUNTER_LEN * chunkId) / numChunks;
			const size_t endByte = (COUNTER_LEN * (chunkId + 1)) / numChunks;
			size_t localHist[16] = {0};

			size_t i = startByte;
			for (; i + 8 <= endByte; i += 8)
			{
				uint64_t word;
				std::memcpy(&word, rawFlat + i, sizeof(uint64_t));
				if (word == 0) continue;
				for (int b = 0; b < 16; ++b)
				{
					++localHist[(word >> (b * 4)) & 0xF];
				}
			}
			for (; i < endByte; ++i)
			{
				uint8_t byte = rawFlat[i];
				if (byte == 0) continue;
				++localHist[byte & 0xF];
				++localHist[byte >> 4];
			}

			histMutex.lock();
			for (int f = 1; f < 16; ++f) globalSmallHist[f] += localHist[f];
			histMutex.unlock();
		};
		processInParallel(chunks, scanChunk, Parameters::get().numThreads, false);

		for (const auto& kv : _hashCounter.lock_table())
		{
			--globalSmallHist[15];
			_kmerDistribution[15 + kv.second] += 1;
		}
		for (size_t f = 1; f < 16; ++f)
		{
			if (globalSmallHist[f] > 0) _kmerDistribution[f] += globalSmallHist[f];
		}

		//FLYE_HIST_VERIFY=1: rebuild the histogram with the original serial loop
		//and abort on any difference.
		const char* verifyEnv = std::getenv("FLYE_HIST_VERIFY");
		if (verifyEnv && verifyEnv[0] == '1')
		{
			KmerDistribution reference;
			for (size_t kmerId = 0; kmerId < COUNTER_LEN * 2; ++kmerId)
			{
				Kmer kmer(kmerId);
				size_t freq = this->getFreq(kmer);
				if (freq > 0) reference[freq] += 1;
			}
			if (reference != _kmerDistribution)
			{
				Logger::get().error() << "FLYE_HIST_VERIFY: k-mer histogram mismatch";
				std::abort();
			}
			Logger::get().debug() << "FLYE_HIST_VERIFY: histogram identical ("
				<< reference.size() << " distinct frequencies)";
		}
	}
	else
	{
		for (const auto& kmer : _hashCounter.lock_table())
		{
			_kmerDistribution[kmer.second] += 1;
		}
	}

	//Logger::get().debug() << "After counter: " 
	//	<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	Logger::get().debug() << "Hash size: " << _hashCounter.size();
	Logger::get().debug() << "Total k-mers " << _numKmers;

	//Counting is finished; nothing writes _hashCounter until clear().
	//FLYE_FREQ_NOFREEZE=1 keeps the original locking lookups (for A/B checks).
	const char* noFreezeEnv = std::getenv("FLYE_FREQ_NOFREEZE");
	if (!(noFreezeEnv && noFreezeEnv[0] == '1'))
	{
		_frozenHash.reset(new cuckoohash_map<Kmer, size_t>::locked_table(
							_hashCounter.lock_table()));
	}
	const char* digestEnv = std::getenv("FLYE_FREQ_DIGEST");
	freqDigestOn = digestEnv && digestEnv[0] == '1';
}

size_t KmerCounter::getFreq(Kmer kmer) const
{
	//kmer.standardForm();

	size_t addCount = 0;
	if (_useFlatCounter)
	{
		size_t arrayPos = kmer.numRepr() / 2;
		bool highBits = kmer.numRepr() % 2;
		//Counting threads have been joined, so a relaxed load returns the same
		//byte. The default seq_cst load is a load-acquire (ldarb) on Arm, and
		//no ordering is needed here (no measurable speed change on its own).
		uint8_t bits = _flatCounter[arrayPos].load(std::memory_order_relaxed);
		uint8_t count = highBits ? (bits >> 4) : (bits & 15);
		if (count < 15)
		{
			return recordFreq(kmer, count);
		}
		else
		{
			addCount = count;
		}
	}

	size_t freq = 0;
	if (_frozenHash)
	{
		auto it = _frozenHash->find(kmer);	//no locks: all are held by _frozenHash
		if (it != _frozenHash->end()) freq = it->second;
	}
	else
	{
		_hashCounter.find(kmer, freq);
	}
	return recordFreq(kmer, freq + addCount);
}

void KmerCounter::clear()
{
	if (freqDigestOn)
	{
		Logger::get().debug() << "FLYE_FREQ_DIGEST: calls " << freqCalls.load()
			<< " digest " << std::hex << freqDigest.load() << std::dec
			<< (_frozenHash ? " (frozen)" : " (locking)");
	}
	_frozenHash.reset();	//release the locks before clearing
	_hashCounter.clear();
	_hashCounter.reserve(0);
	if (_flatCounter)
	{
		free(_flatCounter);
		_flatCounter = nullptr;
	}
}

size_t KmerCounter::getKmerNum() const
{
	if (!_useFlatCounter) return _hashCounter.size();
	return _numKmers;
}
