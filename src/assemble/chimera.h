//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "../sequence/overlap.h"
#include "../sequence/sequence_container.h"
#include <unordered_map>

class ChimeraDetector
{
public:
	ChimeraDetector(const SequenceContainer& readContainer,
					OverlapContainer& ovlpContainer):
		_seqContainer(readContainer),
		_ovlpContainer(ovlpContainer), 
		_overlapCoverage(0)
	{}

	void estimateGlobalCoverage();
	// Requires the complete overlap list, not a capped quickSeqOverlaps
	// prefix: this verdict is cached and reused by subsequent callers.
	bool isChimeric(FastaRecord::Id readId, 
					const std::vector<OverlapRange>& readOvlps);

	float maxCoverageDrop(FastaRecord::Id readId, 
						  const std::vector<OverlapRange>& readOvlps);

	int  getOverlapCoverage() const {return _overlapCoverage;}
	int  getRightTrim(FastaRecord::Id readId);
	bool isRepetitiveRegion(FastaRecord::Id readId, int32_t start, int32_t end, bool debug=false);

private:
	std::vector<int32_t> getReadCoverage(FastaRecord::Id readId,
										 const std::vector<OverlapRange>& readOvlps);

	bool testReadByCoverage(FastaRecord::Id readId,
							const std::vector<OverlapRange>& readOvlps);

	struct CachedCoverage
	{
		CachedCoverage():
			cached(false) {}

		std::vector<int32_t>* coverageFullAln;
		std::vector<int32_t>* coverageIncomleteAln;
		bool cached;
	};
	CachedCoverage getCachedCoverage(FastaRecord::Id readId);

	const SequenceContainer& _seqContainer;
	OverlapContainer& 		 _ovlpContainer;
	cuckoohash_map<FastaRecord::Id, bool> _chimeras;
	cuckoohash_map<FastaRecord::Id, CachedCoverage> _localOvlpsStorage;
	int _overlapCoverage;
};
