#include <atomic>
#include <iostream>
#include <thread>
#include "assemble/chimera.h"
#include "common/config.h"

bool check(int length, int end)
{
    SequenceContainer reads;
    auto id = reads.addSequence(DnaSequence(std::string(length, 'A')), "query").id;
    auto ext = reads.addSequence(DnaSequence(std::string(length, 'A')), "other").id;
    VertexIndex index(reads);
    OverlapDetector detector(reads, index, 1500, 1000, 1500,
                             false, true, 1.0, false, false, false);
    OverlapContainer overlaps(detector, reads);
    OverlapRange o(id, ext, 0, length - end - 1, length, length);
    o.curEnd = end;
    o.extEnd = length - 1;
    std::vector<OverlapRange> forward{o}, reverse{o.complement()};
    ChimeraDetector fwdFirst(reads, overlaps), revFirst(reads, overlaps);
    bool expected = fwdFirst.isChimeric(id, forward);
    bool ok = fwdFirst.isChimeric(id.rc(), reverse) == expected &&
              revFirst.isChimeric(id.rc(), reverse) == expected &&
              revFirst.isChimeric(id, forward) == expected;
    if (end == length - 1) ok &= !expected;
    if (length == 10050 && end == 9000) ok &= expected;
    // Race forward and reverse callers on a fresh cache, repeatedly.
    for (int trial = 0; trial < 10; ++trial)
    {
        ChimeraDetector concurrent(reads, overlaps);
        std::atomic<bool> start(false);
        std::vector<std::thread> workers;
        int results[8] = {};
        for (int i = 0; i < 8; ++i)
            workers.emplace_back([&, i]() {
                while (!start.load()) std::this_thread::yield();
                results[i] = i % 2 ? concurrent.isChimeric(id, forward) :
                                    concurrent.isChimeric(id.rc(), reverse);
            });
        start.store(true);
        for (auto& worker : workers) worker.join();
        for (int result : results) ok &= result == expected;
    }
    if (!ok) std::cerr << "FAIL length=" << length << " overlapEnd=" << end << '\n';
    return ok;
}

int main(int argc, char** argv)
{
    if (argc != 2) return 2;
    Config::load(argv[1]);
    Parameters::get().unevenCoverage = true;
    Parameters::get().minimumOverlap = 1000;
    Parameters::get().kmerSize = 17;
    int failures = 0;
    for (int length : {10000, 10050, 10100})
        for (int end : {8999, 9000, 9049, 9100, length - 1})
            failures += !check(length, end);
    std::cout << "Chimera-strand cases: 15; failures: " << failures << '\n';
    return failures ? 1 : 0;
}
