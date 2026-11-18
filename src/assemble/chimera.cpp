//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <vector>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <iomanip>
#include <cmath>

#include "../common/config.h"
#include "../common/logger.h"
#include "chimera.h"



/*bool ChimeraDetector::isChimeric(FastaRecord::Id readId)
{
	if (!_chimeras.contains(readId))
	{
		const auto& ovlps = _ovlpContainer.lazySeqOverlaps(readId);
		bool result = this->testReadByCoverage(readId, ovlps);
					  //_ovlpContainer.hasSelfOverlaps(readId);
		_chimeras.insert(readId, result);
		_chimeras.insert(readId.rc(), result);
	}
	return _chimeras.find(readId);
}*/

bool ChimeraDetector::isChimeric(FastaRecord::Id readId,
								 const std::vector<OverlapRange>& readOvlps)
{
	// One entry per physical read: separate strand entries can be populated
	// by different racing callers. Coverage windows below use the same
	// forward coordinate system regardless of which strand arrives first.
	auto cacheId = readId.strand() ? readId : readId.rc();
	if (!_chimeras.contains(cacheId))
	{
		bool result = this->testReadByCoverage(readId, readOvlps);
		/*for (const auto& ovlp : IterNoOverhang(readOvlps))
		{
			if (ovlp.curId == ovlp.extId.rc()) 
			{
				int32_t projEnd = ovlp.extLen - ovlp.extEnd - 1;
				if (abs(ovlp.curEnd - projEnd) < JUMP)
				{
					result = true;
				}
			}
		}*/
		_chimeras.insert(cacheId, result);
	}
	return _chimeras.find(cacheId);
}

void ChimeraDetector::estimateGlobalCoverage()
{
	Logger::get().debug() << "Estimating overlap coverage";

	int numSamples = std::min(1000, (int)_seqContainer.iterSeqs().size());
	int sampleRate = (int)_seqContainer.iterSeqs().size() / numSamples;
	//int minCoverage = _inputCoverage / 
	//				(int)Config::get("max_coverage_drop_rate") + 1;
	//int maxCoverage = _inputCoverage * 
	//				(int)Config::get("max_coverage_drop_rate");
	int flankSize = 0;

	std::unordered_map<int32_t, int32_t> readHist;
	std::vector<int32_t> covList;
	
	//std::ofstream fout("../cov_hist.txt");

	int64_t sum = 0;
	int64_t num = 0;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		if (rand() % sampleRate) continue;
		auto coverage = this->getReadCoverage(seq.id, _ovlpContainer.lazySeqOverlaps(seq.id));
		bool nonZero = false;
		for (auto c : coverage) nonZero |= (c != 0);
		if (!nonZero) continue;

		for (size_t i = flankSize; i < coverage.size() - flankSize; ++i)
		{
			{
				++readHist[coverage[i]];
				sum += coverage[i];
				++num;
				covList.push_back(coverage[i]);
			}
		}
	}

	if (readHist.empty())
	{
		Logger::get().warning() << "No overlaps found!";
		_overlapCoverage = 0;
	}
	else
	{
		_overlapCoverage = median(covList);
	}

	Logger::get().info() << "Overlap-based coverage: " << _overlapCoverage;
}

std::vector<int32_t> 
	ChimeraDetector::getReadCoverage(FastaRecord::Id readId,
									 const std::vector<OverlapRange>& readOverlaps)
{
	static const int WINDOW = Config::get("chimera_window");
	const int FLANK = 1;

	std::vector<int> coverage;
	int numWindows = std::ceil((float)_seqContainer.seqLen(readId) / WINDOW) + 1;
	if (numWindows - 2 * FLANK <= 0) return {0};

	coverage.assign(numWindows - 2 * FLANK, 0);
	for (const auto& ovlp : IterNoOverhang(readOverlaps))
	{
		if (ovlp.curId == ovlp.extId.rc() ||
			ovlp.curId == ovlp.extId) continue;

		//skip 2 first/last windows of overlap to be more robust to
		//possible coorinate shifts
		// Window rounding is not symmetric for arbitrary read lengths. Bin
		// both strands in forward coordinates so the cached verdict does not
		// depend on the strand of the first caller.
		int begin = readId.strand() ? ovlp.curBegin : ovlp.curLen - ovlp.curEnd - 1;
		int end = readId.strand() ? ovlp.curEnd : ovlp.curLen - ovlp.curBegin - 1;
		for (int pos = begin / WINDOW + FLANK;
			 pos <= end / WINDOW - FLANK; ++pos)
		{
			//assert(pos - FLANK >= 0 && pos - FLANK < (int)coverage.size());
			++coverage.at(pos - FLANK);
		}
	}

	return coverage;
}


float ChimeraDetector::maxCoverageDrop(FastaRecord::Id readId,
									   const std::vector<OverlapRange>& readOvlps)
{
	auto coverage = this->getReadCoverage(readId, readOvlps);
	if (coverage.empty()) return 0;

	const int CHIMERA_OVERHANG = (int)Config::get("chimera_overhang");
	const int MAX_FLANK = CHIMERA_OVERHANG / Config::get("chimera_window");
	int32_t goodStart = MAX_FLANK;
	int32_t goodEnd = coverage.size() - MAX_FLANK - 1;

	if (goodEnd <= goodStart) return 0;
	float maxDrop = 0;
	for (int32_t i = goodStart; i <= goodEnd - 1; ++i)
	{
		float rate = (float)(std::max(coverage[i], coverage[i + 1]) + 1) / (std::min(coverage[i], coverage[i + 1]) + 1);
		if (maxDrop < rate) maxDrop = rate;
	}

	return maxDrop;
}

bool ChimeraDetector::testReadByCoverage(FastaRecord::Id readId,
										 const std::vector<OverlapRange>& readOvlps)
{
	const float MAX_DROP_RATE = Config::get("max_coverage_drop_rate");

	auto coverage = this->getReadCoverage(readId, readOvlps);
	if (coverage.empty()) return false;

	const int CHIMERA_OVERHANG = (int)Config::get("chimera_overhang");
	const int MAX_FLANK = CHIMERA_OVERHANG / Config::get("chimera_window");
	int32_t goodStart = MAX_FLANK;
	int32_t goodEnd = coverage.size() - MAX_FLANK - 1;

	int32_t maxCov = 0;
	int64_t sumCov = 0;
	for (int32_t i = goodStart; i <= goodEnd; ++i)
	{
		maxCov = std::max(maxCov, coverage[i]);
		sumCov += coverage[i];
	}
	int32_t medianCoverage = median(coverage);
	if (sumCov == 0) return true;

	int threshold = 0;	
	if (!Parameters::get().unevenCoverage)
	{
		threshold = std::max(1L, std::lround((float)_overlapCoverage / 
											 MAX_DROP_RATE));
	}
	else
	{
		threshold = std::min(10L, std::max(1L, std::lround(medianCoverage / MAX_DROP_RATE)));
	}

	bool lowCoverage = false;
	if (goodEnd <= goodStart) lowCoverage = true;
	for (int32_t i = goodStart; i <= goodEnd; ++i)
	{
		if (coverage[i] < threshold)
		{
			lowCoverage = true;
			break;
		}
	}

	/*static std::ofstream fout("chimera_dump.txt");
	static std::mutex logLock;
	logLock.lock();
	std::string covStr;
	for (int32_t i = 0; i < (int)coverage.size() - 1; ++i)
	{
		if (i == goodStart) covStr += " { ";
		covStr += std::to_string(coverage[i]) + " ";
		if (i == goodEnd) covStr += " } ";
	}
	fout << _seqContainer.seqName(readId) << " " 
		<< _seqContainer.seqLen(readId) << std::endl;
	fout << covStr << std::endl;
	fout << "max: " << maxCov << " median: " << medianCoverage
		<< " threshold: " << threshold << " chim:" << lowCoverage << std::endl << std::endl;
	logLock.unlock();*/

	return lowCoverage;
}

bool ChimeraDetector::isRepetitiveRegion(FastaRecord::Id readId, int32_t start, 
										 int32_t end, bool debug)
{
	const float HANG_END_RATE = 0.75f;
	const float REPEAT_WINDOW_RATE = 0.75f;
	const int WINDOW = Config::get("chimera_window");
	
	/*int numWindows = std::ceil((float)_seqContainer.seqLen(readId) / WINDOW) + 1;
	int vecSize = numWindows - 2 * FLANK;
	if (vecSize <= 0) return false;

	std::vector<int> coverage(vecSize, 0);
	std::vector<int> junctions(vecSize, 0);
	for (const auto& ovlp : this->getLocalOverlaps(readId))
	//for (const auto& ovlp : _containNoOverhangs.lazySeqOverlaps(readId))
	{
		if (ovlp.curId == ovlp.extId.rc() ||
			ovlp.curId == ovlp.extId) continue;

		for (int pos = ovlp.curBegin / WINDOW + FLANK; 
			 pos <= ovlp.curEnd / WINDOW - FLANK; ++pos)
		{
			if (ovlp.lrOverhang() > MAX_OVERHANG)
			{
				++junctions.at(pos - FLANK);
			}
			else
			{
				++coverage.at(pos - FLANK);
			}
		}
	}*/

	auto cachedCoverage = this->getCachedCoverage(readId);
	const std::vector<int32_t>& coverage = *cachedCoverage.coverageFullAln;
	const std::vector<int32_t>& junctions = *cachedCoverage.coverageIncomleteAln;

	int numSuspicious = 0;
	int rangeLen = 0;
	for (int32_t pos = std::max(0, start / WINDOW); 
		 pos < std::min((int32_t)coverage.size(), end / WINDOW); ++pos)
	{
		if (HANG_END_RATE * coverage[pos] <= junctions[pos])	//if both are 0, it's suspicious
		{
			++numSuspicious;
		}
		++rangeLen;
	}

	if (debug)
	{
		//static std::mutex logLock;
		//logLock.lock();
		Logger::get().debug() << "Checking repeat";
		std::string covStr;
		std::string juncStr;
		for (int32_t i = 0; i < (int)coverage.size() - 1; ++i)
		{
			covStr += std::to_string(coverage[i]) + " ";
			juncStr += std::to_string(junctions[i]) + " ";
		}
		//Logger::get().debug() << _seqContainer.seqName(readId);
		Logger::get().debug() << _seqContainer.seqLen(readId) << " " << start << " " << end;
		Logger::get().debug() << covStr;
		Logger::get().debug() << juncStr;
		
		if ((float)numSuspicious / rangeLen > REPEAT_WINDOW_RATE) Logger::get().debug() << "Flagged";
	}

	if ((float)numSuspicious / rangeLen > REPEAT_WINDOW_RATE)
	{
		return true;
	}
	return false;
}

ChimeraDetector::CachedCoverage
	ChimeraDetector::getCachedCoverage(FastaRecord::Id readId)
{
	//upsert creates default value if it does not exist
	CachedCoverage cached;
	_localOvlpsStorage.upsert(readId, 	
		[&cached](CachedCoverage& val)
			{cached = val;}, CachedCoverage());
	if (cached.cached)
	{
		return cached;
	}

	//not cached - need to copmute
	const int WINDOW = Config::get("chimera_window");
	const int MAX_OVERHANG = Config::get("maximum_overhang");
	const int FLANK = 1;

	int numWindows = std::ceil((float)_seqContainer.seqLen(readId) / WINDOW) + 1;
	int vecSize = numWindows - 2 * FLANK;
	if (vecSize <= 0) throw std::runtime_error("Zero-sized coverage vector");

	std::vector<int> coverage(vecSize, 0);
	std::vector<int> junctions(vecSize, 0);
	auto overlaps = _ovlpContainer.quickSeqOverlaps(readId, /*max ovlps*/ 0, /*force local*/ true);
	for (const auto& ovlp : overlaps)
	{
		if (ovlp.curId == ovlp.extId.rc() ||
			ovlp.curId == ovlp.extId) continue;

		for (int pos = ovlp.curBegin / WINDOW + FLANK; 
			 pos <= ovlp.curEnd / WINDOW - FLANK; ++pos)
		{
			if (ovlp.lrOverhang() > MAX_OVERHANG)
			{
				++junctions.at(pos - FLANK);
			}
			else
			{
				++coverage.at(pos - FLANK);
			}
		}
	}
	//

	//updating cache
	_localOvlpsStorage.update_fn(readId,
		[&cached, &coverage, &junctions, this]
		(CachedCoverage& val)
		{
			if (!val.cached)
			{
				val.coverageFullAln = new std::vector<int32_t>;
				val.coverageFullAln->swap(coverage);
				val.coverageIncomleteAln = new std::vector<int32_t>;
				val.coverageIncomleteAln->swap(junctions);
				val.cached = true;
			}
			cached = val;
		});

	return cached;

}
