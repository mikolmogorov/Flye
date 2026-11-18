#include <iostream>
#include <random>
#include "assemble/extender.h"
#include "common/config.h"

std::string dna(std::mt19937& rng, size_t size)
{
    std::string result(size, 'A');
    for (auto& base : result) base = "ACGT"[rng() % 4];
    return result;
}

int main(int argc, char** argv)
{
    if (argc != 3) return 2;
    Config::load(argv[1]);
    bool fullFirst = std::string(argv[2]) == "full-first";
    Parameters::get().unevenCoverage = true;
    Parameters::get().minimumOverlap = 1000;
    Parameters::get().kmerSize = 17;
    Parameters::get().numThreads = 1;
    std::mt19937 rng(42);
    std::string query = dna(rng, 10000);
    std::vector<std::pair<std::string, std::string>> partial, complete;
    for (int i = 0; i < 100; ++i)
        partial.emplace_back(dna(rng, 5500) + query.substr(0, 4500), "left_" + std::to_string(i));
    for (int i = 0; i < 31; ++i)
        complete.emplace_back(query, "full_" + std::to_string(i));
    SequenceContainer reads;
    auto add = [&](const std::vector<std::pair<std::string, std::string>>& records) {
        for (const auto& record : records)
            reads.addSequence(DnaSequence(record.first), record.second);
    };
    // No early full-length query read: every complete seed must encounter
    // the same 100 partial overlaps first in the partial-first permutation.
    if (fullFirst) { add(complete); add(partial); }
    else { add(partial); add(complete); }
    reads.buildPositionIndex();
    VertexIndex index(reads);
    index.buildIndexMinimizers(1, 5);
    OverlapDetector detector(reads, index, 1500, 1000, 1500,
                             false, true, 1.0, false, false, false);
    OverlapContainer overlaps(detector, reads);
    Extender extender(reads, overlaps, 1000);
    extender.assembleDisjointigs();
    ConsensusGenerator consensus;
    auto sequences = consensus.generateConsensuses(extender.getDisjointigPaths());
    bool recovered = false;
    std::string reverse = DnaSequence(query).complement().str();
    for (const auto& sequence : sequences)
    {
        auto bases = sequence.sequence.str();
        recovered |= bases.find(query) != std::string::npos ||
                     bases.find(reverse) != std::string::npos;
    }
    std::cout << argv[2] << ": disjointigs=" << sequences.size()
              << " query recovered=" << recovered << '\n';
    return recovered ? 0 : 1;
}
