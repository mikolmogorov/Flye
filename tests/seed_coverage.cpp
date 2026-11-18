#include <iostream>
#include <random>
#include "assemble/chimera.h"
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
    for (int i = 0; i < 30; ++i)
        complete.emplace_back(query, "full_" + std::to_string(i));
    SequenceContainer reads;
    auto id = reads.addSequence(DnaSequence(query), "query").id;
    auto add = [&](const std::vector<std::pair<std::string, std::string>>& records) {
        for (const auto& record : records)
            reads.addSequence(DnaSequence(record.first), record.second);
    };
    if (fullFirst) { add(complete); add(partial); }
    else { add(partial); add(complete); }
    reads.buildPositionIndex();
    VertexIndex index(reads);
    index.buildIndexMinimizers(1, 5);
    OverlapDetector detector(reads, index, 1500, 1000, 1500,
                             false, true, 1.0, false, false, false);
    OverlapContainer overlaps(detector, reads);
    auto capped = overlaps.quickSeqOverlaps(id, 100);
    // Same complete, uncached evidence required by Extender's seed screening.
    auto full = overlaps.quickSeqOverlaps(id, 0);
    if (overlaps.indexSize() != 0) return 3;
    ChimeraDetector seedFirst(reads, overlaps), extensionFirst(reads, overlaps);
    ChimeraDetector incomplete(reads, overlaps);
    bool seed = seedFirst.isChimeric(id, full);
    bool cached = seedFirst.isChimeric(id, overlaps.lazySeqOverlaps(id));
    bool extension = extensionFirst.isChimeric(id, overlaps.lazySeqOverlaps(id));
    bool seedAfterExtension = extensionFirst.isChimeric(id, full);
    bool prefix = incomplete.isChimeric(id, capped);
    std::cout << argv[2] << ": full=" << full.size() << " capped=" << capped.size()
              << " complete verdicts=" << seed << cached << extension << seedAfterExtension
              << " prefix verdict=" << prefix << '\n';
    // Also demonstrate that this fixture is sensitive to the original bug:
    // the ordered prefix is falsely chimeric only in the partial-first order.
    return full.size() != 130 || capped.size() != 100 || seed || cached ||
           extension || seedAfterExtension || prefix == fullFirst;
}
