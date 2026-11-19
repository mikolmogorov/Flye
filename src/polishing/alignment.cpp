//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "alignment.h"
#include <algorithm>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>
#include "../polishing/subs_matrix.h"

Alignment::Alignment(size_t size, const SubstitutionMatrix& sm):
	_forwardScores(size),
	_reverseScores(size),
	_subsMatrix(sm)
{ 
}


AlnScoreType Alignment::globalAlignment(const std::string& consensus,
							 			const std::vector<std::string>& reads)
{
	AlnScoreType finalScore = 0;
	std::string revConsensus(consensus.rbegin(), consensus.rend());
	std::string revRead;

	for (size_t readId = 0; readId < _forwardScores.size(); ++readId)
	{
		unsigned int x = consensus.size() + 1;
		unsigned int y = reads[readId].size() + 1;
		_forwardScores[readId].resize(x, y);
		AlnScoreType score = this->getScoringMatrix(consensus, reads[readId], 
													_forwardScores[readId]);

		revRead.assign(reads[readId].rbegin(), reads[readId].rend());
		_reverseScores[readId].resize(x, y);
		this->getScoringMatrix(revConsensus, revRead, _reverseScores[readId]);

		finalScore += score;
	}
	return finalScore;
}

AlnScoreType Alignment::addDeletion(unsigned int letterIndex) const
{
	AlnScoreType finalScore = 0;
	for (size_t readId = 0; readId < _forwardScores.size(); ++readId)
	{
		const ScoreMatrix& forwardScore = _forwardScores[readId];
		const ScoreMatrix& reverseScore = _reverseScores[readId];

		size_t frontRow = letterIndex - 1;
		size_t revRow = reverseScore.nrows() - 1 - letterIndex;
		size_t ncols = forwardScore.ncols();
		const AlnScoreType* fwdPtr = forwardScore.rowPtr(frontRow);
		const AlnScoreType* revPtr = reverseScore.rowPtr(revRow);
		
		AlnScoreType maxVal = std::numeric_limits<AlnScoreType>::lowest();
		for (size_t col = 0; col < ncols; ++col)
		{
			AlnScoreType sum = fwdPtr[col] + revPtr[ncols - col - 1];
			maxVal = std::max(maxVal, sum);
		}
		finalScore += maxVal;
	}
	return finalScore;
}


AlnScoreType Alignment::addSubstitution(unsigned int letterIndex, char base, 
										const std::vector<std::string>& reads) const
{
	AlnScoreType finalScore = 0;
	const AlnScoreType baseGapScore = _subsMatrix.getScore(base, '-');
	for (size_t readId = 0; readId < reads.size(); ++readId)
	{
		const ScoreMatrix& forwardScore = _forwardScores[readId];
		const ScoreMatrix& reverseScore = _reverseScores[readId];
		const std::string& read = reads[readId];
		const size_t readLen = read.size();

		size_t frontRow = letterIndex - 1;
		size_t revRow = reverseScore.nrows() - 1 - letterIndex;
		const AlnScoreType* fwdPtr = forwardScore.rowPtr(frontRow);
		const AlnScoreType* revPtr = reverseScore.rowPtr(revRow);

		AlnScoreType maxVal = (fwdPtr[0] + baseGapScore) + revPtr[readLen];
		for (size_t i = 0; i < readLen; ++i)
		{
			AlnScoreType match = fwdPtr[i] + _subsMatrix.getScore(base, read[i]);
			AlnScoreType ins = fwdPtr[i + 1] + baseGapScore;
			AlnScoreType subVal = std::max(match, ins);
			maxVal = std::max(maxVal, subVal + revPtr[readLen - 1 - i]);
		}
		finalScore += maxVal;
	}
	return finalScore;
}


AlnScoreType Alignment::addInsertion(unsigned int pos, char base, 
									 const std::vector<std::string>& reads) const
{
	AlnScoreType finalScore = 0;
	const AlnScoreType baseGapScore = _subsMatrix.getScore(base, '-');
	for (size_t readId = 0; readId < reads.size(); ++readId)
	{
		const ScoreMatrix& forwardScore = _forwardScores[readId];
		const ScoreMatrix& reverseScore = _reverseScores[readId];
		const std::string& read = reads[readId];
		const size_t readLen = read.size();

		size_t frontRow = pos - 1;
		size_t revRow = reverseScore.nrows() - pos;
		const AlnScoreType* fwdPtr = forwardScore.rowPtr(frontRow);
		const AlnScoreType* revPtr = reverseScore.rowPtr(revRow);

		AlnScoreType maxVal = (fwdPtr[0] + baseGapScore) + revPtr[readLen];
		for (size_t i = 0; i < readLen; ++i)
		{
			AlnScoreType match = fwdPtr[i] + _subsMatrix.getScore(base, read[i]);
			AlnScoreType ins = fwdPtr[i + 1] + baseGapScore;
			AlnScoreType subVal = std::max(match, ins);
			maxVal = std::max(maxVal, subVal + revPtr[readLen - 1 - i]);
		}
		finalScore += maxVal;
	}
	return finalScore;
}


AlnScoreType Alignment::getScoringMatrix(const std::string& v, 
										 const std::string& w,
								  		 ScoreMatrix& scoreMat) 
{
	AlnScoreType score = 0;
	scoreMat.at(0, 0) = 0;
	
	for (size_t i = 0; i < v.size(); i++) 
	{
		AlnScoreType score = _subsMatrix.getScore(v[i], '-');
		scoreMat.at(i + 1, 0) = scoreMat.at(i, 0) + score;
	}


	for (size_t i = 0; i < w.size(); i++) {
		AlnScoreType score = _subsMatrix.getScore('-', w[i]);
		scoreMat.at(0, i + 1) = scoreMat.at(0, i) + score;
	}


	for (size_t i = 1; i < v.size() + 1; i++)
	{
		char key1 = v[i - 1];
		AlnScoreType key1Gap = _subsMatrix.getScore(key1, '-');
		const AlnScoreType* prevRow = scoreMat.rowPtr(i - 1);
		AlnScoreType* curRow = scoreMat.rowPtr(i);
		for (size_t j = 1; j < w.size() + 1; j++) 
		{
			char key2 = w[j - 1];

			AlnScoreType left = curRow[j - 1] + _subsMatrix.getScore('-', key2);
			AlnScoreType up = prevRow[j] + key1Gap;
			score = std::max(left, up);

			AlnScoreType cross = prevRow[j - 1] + _subsMatrix.getScore(key1, key2);
			score = std::max(score, cross);
			curRow[j] = score;
		}
	}

	return score;
}
