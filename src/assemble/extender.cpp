//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <algorithm>
#include <iomanip>
#include <stack>
#include <cmath>

#include "../common/config.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "extender.h"

namespace
{
	OverlapRange getOverlapBetween(OverlapContainer& ovlpCnt, FastaRecord::Id readOne, FastaRecord::Id readTwo)
	{
		bool found = false;
		OverlapRange readsOvlp;

		for (const auto& ovlp : ovlpCnt.lazySeqOverlaps(readOne))
		{
			if (ovlp.extId == readTwo) 
			{
				readsOvlp = ovlp;
				found = true;
				break;
			}
		}
		for (const auto& ovlp : ovlpCnt.lazySeqOverlaps(readTwo))
		{
			if (ovlp.extId == readOne) 
			{
				if (!found || readsOvlp.minRange() < ovlp.minRange())
				{
					readsOvlp = ovlp.reverse();
					found = true;
					break;
				}
			}
		}
		if (!found) throw std::runtime_error("Ovlp not found!");

		return readsOvlp;
	}
}

Extender::ExtensionInfo Extender::extendDisjointig(FastaRecord::Id startRead)
{

	std::unordered_set<FastaRecord::Id> currentReads;
	currentReads.insert(startRead);
	currentReads.insert(startRead.rc());

	bool rightExtension = true;
	FastaRecord::Id currentRead = startRead;
	std::vector<int> numExtensions;
	std::vector<int> overlapSizes;
	ExtensionInfo exInfo;
	exInfo.reads.push_back(startRead);
	exInfo.assembledLength = _readsContainer.seqLen(startRead);

	auto startOverlaps = _ovlpContainer.lazySeqOverlaps(startRead);
	auto leftExtendsStart = [startRead, this, &startOverlaps](const FastaRecord::Id readId)
	{
		for (const auto& ovlp : IterNoOverhang(startOverlaps))
		{
			if (ovlp.extId == readId && this->extendsLeft(ovlp)) return true;
		}
		return false;
	};

	while(true)
	{
		const std::vector<OverlapRange>& curOverlaps = _ovlpContainer.lazySeqOverlaps(currentRead);
		std::vector<OverlapRange> extensions;
		for (const auto& ovlp : IterNoOverhang(curOverlaps))
		{
			if (this->extendsRight(ovlp)) extensions.push_back(ovlp);
		}
		numExtensions.push_back(extensions.size());

		//sort from longes to shortest overlap
		std::sort(extensions.begin(), extensions.end(), 
				  [](const OverlapRange& a, const OverlapRange& b)
					 {return a.curRange() > b.curRange();});

		//bool foundExtension = false;
		const float COV_DROP = Config::get("max_extensions_drop_rate");
		int minExtensions = std::roundf((float)median(numExtensions) / COV_DROP);
		minExtensions = std::min(10, std::max(1, minExtensions));

		/*Logger::get().debug() << _readsContainer.seqName(currentRead) 
			<< "\t" << overlaps.size();
		for (auto& ovlp : overlaps)
		{
			if (this->extendsRight(ovlp))
			{
				Logger::get().debug() << "\t" << 
					_readsContainer.seqName(ovlp.extId) << "\t" 
					<< "o:" << ovlp.curRange() << "\tchim:"
					<< _chimDetector.isChimeric(ovlp.extId)
					<< "\trep:" << this->isRightRepeat(ovlp.extId)
					<< "\text:" << this->countRightExtensions(ovlp.extId);
			}
		}*/

		const OverlapRange* bestPreferered = nullptr;
		const OverlapRange* bestSuspicious = nullptr;
		const OverlapRange* bestDeadEnd = nullptr;
		for (const auto& ovlp : extensions)
		{
			//need to check this condition, otherwise there will
			//be complications when we initiate extension of startRead to the left
			if (leftExtendsStart(ovlp.extId)) continue;

			//only using reads longer than safeOverlap
			if (ovlp.extLen < _safeOverlap) continue;

			//if overlap is shorter than minOverlap parameter (which happens rarely)
			//do an additional repeat check. 
			if (ovlp.minRange() < _safeOverlap)
			{
				bool curRepeat = _chimDetector.isRepetitiveRegion(ovlp.curId, ovlp.curBegin, ovlp.curEnd);
				bool extRepeat = _chimDetector.isRepetitiveRegion(ovlp.extId, ovlp.extBegin, ovlp.extEnd);
				if (curRepeat && extRepeat) continue;
			}

			const std::vector<OverlapRange>& extOverlaps = _ovlpContainer.lazySeqOverlaps(ovlp.extId);

			const float MAX_COVERAGE_DROP = 5.0f;
			if (_chimDetector.isChimeric(ovlp.extId, extOverlaps) &&
				_chimDetector.maxCoverageDrop(ovlp.extId, extOverlaps) > MAX_COVERAGE_DROP) continue;

			//optimistically, pick the first available highly reliable extenion (which will
			//also be with the longest overlap as extensions are sorted based on that)
			if (!_chimDetector.isChimeric(ovlp.extId, extOverlaps) &&
				this->countRightExtensions(extOverlaps) >= minExtensions &&
				ovlp.minRange() > _safeOverlap)
			{
				bestPreferered = &ovlp;
				break;
			}

			//alternatives, in the decreasing order of preference
			else if (this->countRightExtensions(extOverlaps) > 0)
			{
				if (!bestSuspicious) bestSuspicious = &ovlp;
				if (ovlp.minRange() < _safeOverlap) break;
			}
			else
			{
				if (!bestDeadEnd || bestDeadEnd->rightShift() < ovlp.rightShift())
				{
					bestDeadEnd = &ovlp;
				}
			}
		}

		const OverlapRange* selectedExtension = nullptr;
		if (bestPreferered)
		{
			selectedExtension = bestPreferered;
		}
		else if (bestSuspicious)
		{
			selectedExtension = bestSuspicious;
			++exInfo.numSuspicious;
		}
		else if (bestDeadEnd)
		{
			selectedExtension = bestDeadEnd;
			++exInfo.numSuspicious;
		}

		if (selectedExtension)
		{
			exInfo.assembledLength += selectedExtension->rightShift();
			currentRead = selectedExtension->extId;
			if (selectedExtension->minRange() < _safeOverlap) ++exInfo.shortExtensions;
			exInfo.reads.push_back(currentRead);
			overlapSizes.push_back(selectedExtension->curRange());

			//_chimDetector.isRepetitiveRegion(selectedExtension->curId, selectedExtension->curBegin, 
			//								 selectedExtension->curEnd, true);
			//_chimDetector.isRepetitiveRegion(selectedExtension->extId, selectedExtension->extBegin, 
			//								 selectedExtension->extEnd, true);
		}
		else
		{
			rightExtension ? exInfo.leftTip = true : exInfo.rightTip = true;
		}

		if (!selectedExtension || _innerReads.contains(currentRead) ||
			currentReads.count(currentRead))
		{
			//Logger::get().debug() << "Not found: " << !foundExtension << 
			//	" overlaps visited: " << overlapsVisited;

			//right extension done, try to extend left from start read
			if (rightExtension && !exInfo.reads.empty())
			{
				exInfo.stepsToTurn = exInfo.reads.size();
				rightExtension = false;
				currentRead = exInfo.reads.front().rc();
				std::reverse(exInfo.reads.begin(), exInfo.reads.end());
				for (size_t i = 0; i < exInfo.reads.size(); ++i) 
				{
					exInfo.reads[i] = exInfo.reads[i].rc();
				}
			}

			//both left and right are done
			else
			{
				break;
			}
		}

		currentReads.insert(currentRead);
		currentReads.insert(currentRead.rc());
	}

	if (!numExtensions.empty())
	{
		exInfo.meanOverlaps = median(numExtensions);
	}
	if (!overlapSizes.empty())
	{
		exInfo.avgOverlapSize = median(overlapSizes);
		exInfo.minOverlapSize = *std::min_element(overlapSizes.begin(), 
												  overlapSizes.end());
	}

	/*if (_innerReads.contains(exInfo.reads.front())) 
	{
		exInfo.leftAsmOverlap = _readsContainer.seqLen(exInfo.reads.front());
	}
	if (_innerReads.contains(exInfo.reads.back())) 
	{
		exInfo.rightAsmOverlap = _readsContainer.seqLen(exInfo.reads.back());
	}*/

	return exInfo;
}


void Extender::assembleDisjointigs()
{
	Logger::get().info() << "Extending reads";
	_chimDetector.estimateGlobalCoverage();
	_ovlpContainer.overlapDivergenceStats();
	_innerReads.clear();
	cuckoohash_map<FastaRecord::Id, size_t> coveredReads;
	
	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _readsContainer.iterSeqs())
	{
		if (seq.sequence.length() > (size_t)Parameters::get().minimumOverlap &&
			seq.id.strand())
		{
			allReads.push_back(seq.id);
		}
	}
	int totalReads = allReads.size() * 2;	//counting both strands
	
	std::mutex indexMutex;
	ProgressPercent progress(totalReads);
	progress.setValue(0);
	auto processRead = [this, &indexMutex, &coveredReads, totalReads, &progress] 
		(FastaRecord::Id startRead)
	{
		//most of the reads will fall into the inner categoty -
		//so no further processing will be needed
		if (_innerReads.contains(startRead)) return;

		coveredReads.insert(startRead);
		coveredReads.insert(startRead.rc());

		//getting overlaps without caching first - so we don't
		//store overlap information for many trashy reads
		//that won't result into disjointig extension
		auto startOvlps = _ovlpContainer.quickSeqOverlaps(startRead, 
														  /*max overlaps*/ 0);
		int numInnerOvlp = 0;
		int totalOverlaps = 0;
		for (const auto& ovlp : IterNoOverhang(startOvlps))
		{
			if (_innerReads.contains(ovlp.extId)) ++numInnerOvlp;
			++totalOverlaps;
		}

		//int maxStartExt = _chimDetector.getOverlapCoverage() * 10;
		//int minStartExt = 1;
		//int extLeft = this->countLeftExtensions(startOvlps);
		//int extRight = this->countRightExtensions(startOvlps);

		if (_chimDetector.isChimeric(startRead, startOvlps) ||
			_readsContainer.seqLen(startRead) < _safeOverlap) return;

		const bool aggressiveDupFilt = (int)Config::get("aggressive_dup_filter");
		if (aggressiveDupFilt && numInnerOvlp > totalOverlaps / 2) return;
		
		//Good to go!
		ExtensionInfo exInfo = this->extendDisjointig(startRead);

		//Exclusive part - updating the overall assembly
		std::lock_guard<std::mutex> guard(indexMutex);

		/*if (exInfo.reads.size() - exInfo.numSuspicious < 
			(size_t)Config::get("min_reads_in_disjointig"))
		{
			//Logger::get().debug() << "Thrown away: " << exInfo.reads.size() << " " << exInfo.numSuspicious
			//	<< " " << exInfo.leftTip << " " << exInfo.rightTip;
			return;
		}*/
		
		int innerCount = 0;
		//do not count first and last reads - they are inner by defalut
		for (size_t i = 1; i < exInfo.reads.size() - 1; ++i)
		{
			if (_innerReads.contains(exInfo.reads[i])) ++innerCount;
		}
		int innerThreshold = std::min((int)Config::get("max_inner_reads"),
									  int((float)Config::get("max_inner_fraction") * 
										  exInfo.reads.size()));
		if (innerCount > innerThreshold)
		{
			Logger::get().debug() << "Discarded disjointig with "
				<< exInfo.reads.size() << " reads and "
				<< innerCount << " inner overlaps";
			return;
		}

		Logger::get().debug() << "Assembled disjointig " 
			<< std::to_string(_readLists.size() + 1)
			<< "\n\tWith " << exInfo.reads.size() << " reads"
			<< "\n\tStart read: " << _readsContainer.seqName(startRead)
			<< "\n\tAt position: " << exInfo.stepsToTurn
			<< "\n\tleftTip: " << exInfo.leftTip 
			<< " rightTip: " << exInfo.rightTip
			<< "\n\tSuspicious: " << exInfo.numSuspicious
			<< "\n\tShort ext: " << exInfo.shortExtensions
			<< "\n\tMean extensions: " << exInfo.meanOverlaps
			<< "\n\tAvg overlap len: " << exInfo.avgOverlapSize
			<< "\n\tMin overlap len: " << exInfo.minOverlapSize
			<< "\n\tInner reads: " << innerCount
			<< "\n\tLength: " << exInfo.assembledLength;

		//Logger::get().debug() << "Ovlp index size: " << _ovlpContainer.indexSize();
		
		//update inner read index
		std::vector<OverlapRange> allOverlaps;
		for (const auto& readId : exInfo.reads)
		{
			coveredReads.insert(readId, true);
			coveredReads.insert(readId.rc(), true);
			_innerReads.insert(readId, true);
			_innerReads.insert(readId.rc(), true);

			for (const auto& ovlp : IterNoOverhang(_ovlpContainer.lazySeqOverlaps(readId)))
			{
				allOverlaps.push_back(ovlp);
				if (ovlp.minRange() > _safeOverlap)
				{
					coveredReads.insert(ovlp.extId, true);
					coveredReads.insert(ovlp.extId.rc(), true);
				}
			}
		}

		auto innerReads = this->getInnerReads(allOverlaps);
		for (const auto& read : innerReads)
		{
			_innerReads.insert(read, true);
			_innerReads.insert(read.rc(), true);
		}

		Logger::get().debug() << "Inner: " << 
			_innerReads.size() << " covered: " << coveredReads.size()
			<< " total: "<< totalReads;
		progress.setValue(coveredReads.size());
		
		_readLists.push_back(std::move(exInfo));
	};

	std::function<void(const FastaRecord::Id&)> threadWorker = 
		[processRead] (const FastaRecord::Id& readId)
	{
		processRead(readId);
	};

	//deterministic shuffling
	std::sort(allReads.begin(), allReads.end(), 
			  [](const FastaRecord::Id& id1, const FastaRecord::Id& id2)
			  {return id1.hash() < id2.hash();});
	processInParallel(allReads, threadWorker,
					  Parameters::get().numThreads, /*progress*/ false);
	progress.setDone();

	/*bool addSingletons = (bool)Config::get("add_unassembled_reads");
	if (addSingletons)
	{
		std::vector<FastaRecord::Id> sortedByLength;
		for (const auto& seq : _readsContainer.iterSeqs())
		{
			if (seq.id.strand() && !_innerReads.contains(seq.id) &&
				_readsContainer.seqLen(seq.id) > _safeOverlap)
			{
				sortedByLength.push_back(seq.id);
			}
		}
		std::sort(sortedByLength.begin(), sortedByLength.end(),
				  [this](FastaRecord::Id idOne, FastaRecord::Id idTwo)
				  	{return _readsContainer.seqLen(idOne) > 
							_readsContainer.seqLen(idTwo);});

		int singletonsAdded = 0;
		std::unordered_set<FastaRecord::Id> coveredLocal;
		for (const auto& readId : sortedByLength)
		{
			if (!coveredLocal.count(readId))
			{
				for (const auto& ovlp : IterNoOverhang(_ovlpContainer.lazySeqOverlaps(readId)))
				{
					if (ovlp.leftShift() >= 0 && ovlp.rightShift() <= 0)
					{
						coveredLocal.insert(ovlp.extId);
						coveredLocal.insert(ovlp.extId.rc());
					}
				}
				ExtensionInfo path;
				path.singleton = true;
				path.reads.push_back(readId);
				_readLists.push_back(path);
				++singletonsAdded;
			}
		}
		Logger::get().info() << "Added " << singletonsAdded << " singleton reads";
	}*/

	this->convertToDisjointigs();
	Logger::get().info() << "Assembled " << _disjointigPaths.size() 
		<< " disjointigs";
}

std::vector<FastaRecord::Id> 
	Extender::getInnerReads(const std::vector<OverlapRange>& ovlps)
{
	static const int WINDOW = Config::get("chimera_window");
	static const int OVERHANG = Config::get("maximum_overhang");

	std::unordered_map<FastaRecord::Id, 
					   std::vector<int32_t>> readsCoverage;
	for (const auto& ovlp: ovlps)
	{
		auto& coverage = readsCoverage[ovlp.extId];
		if (coverage.empty())
		{
			int numWindows = _readsContainer.seqLen(ovlp.extId) / WINDOW;
			if (numWindows < 1)
			{
				Logger::get().warning() << "Wrong read length: " << numWindows;
				numWindows = 1;
			}
			coverage.assign(numWindows, 0);
		}

		for (int pos = ovlp.extBegin / WINDOW + 1; 
			 pos < ovlp.extEnd / WINDOW; ++pos)
		{
			++coverage[pos];
		}
	}

	std::vector<FastaRecord::Id> innerReads;
	for (auto rc : readsCoverage)
	{
		int leftZeros = 0;
		for (size_t i = 0; i < rc.second.size(); ++i)
		{
			if (rc.second[i] != 0) break;
			++leftZeros;
		}
		int rightZeros = 0;
		for (size_t i = 0; i < rc.second.size(); ++i)
		{
			if (rc.second[rc.second.size() - i - 1] != 0) break;
			++rightZeros;
		}
		bool middleZero = false;
		for (size_t i = leftZeros + 1; i < rc.second.size() - rightZeros; ++i)
		{
			if (rc.second[i] == 0) middleZero = true;
		}

		if (!middleZero && leftZeros < OVERHANG / WINDOW && 
			rightZeros < OVERHANG / WINDOW) innerReads.push_back(rc.first);

		/*if (middleZero)
		{
			std::string covStr;
			for (int cov : rc.second)
			{
				covStr += std::to_string(cov) + " ";
			}
			Logger::get().debug() << covStr;
		}*/
	}

	return innerReads;
}

void Extender::convertToDisjointigs()
{
	for (const auto& exInfo : _readLists)
	{
		ContigPath path;
		if (!exInfo.singleton)
		{
			path.name = "disjointig_" + std::to_string(_disjointigPaths.size() + 1);
		}
		else
		{
			path.name = "read_" + std::to_string(_disjointigPaths.size() + 1);
		}
		//path.trimLeft = std::max(0, exInfo.leftAsmOverlap - 
		//							2 * Parameters::get().minimumOverlap);
		//path.trimRight = std::max(0, exInfo.rightAsmOverlap - 
		//							 2 * Parameters::get().minimumOverlap);

		for (size_t i = 0; i < exInfo.reads.size() - 1; ++i)
		{
			auto readsOvlp = getOverlapBetween(_ovlpContainer, exInfo.reads[i], 
											   exInfo.reads[i + 1]);
			path.sequences.push_back(_readsContainer.getSeq(exInfo.reads[i]));
			path.overlaps.push_back(readsOvlp);
		}
		path.sequences.push_back(_readsContainer.getSeq(exInfo.reads.back()));
		_disjointigPaths.push_back(std::move(path));
	}
}

int Extender::countRightExtensions(const std::vector<OverlapRange>& ovlps) const
{
	int count = 0;
	for (const auto& ovlp : IterNoOverhang(ovlps))
	{
		if (this->extendsRight(ovlp)) ++count;
	}
	return count;
}

/*int Extender::countRightExtensions(FastaRecord::Id readId) const
{
	return this->countRightExtensions(_ovlpContainer.lazySeqOverlaps(readId));
}*/

bool Extender::extendsRight(const OverlapRange& ovlp) const
{
	static const int MAX_JUMP = Config::get("maximum_jump");
	return ovlp.rightShift() > MAX_JUMP;
}

int Extender::countLeftExtensions(const std::vector<OverlapRange>& ovlps) const
{
	int count = 0;
	for (const auto& ovlp : IterNoOverhang(ovlps))
	{
		if (this->extendsLeft(ovlp)) ++count;
	}
	return count;
}

bool Extender::extendsLeft(const OverlapRange& ovlp) const
{
	static const int MAX_JUMP = Config::get("maximum_jump");
	return ovlp.leftShift() < -MAX_JUMP;
}
