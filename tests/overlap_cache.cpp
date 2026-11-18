#include <atomic>
#include <iostream>
#include <thread>
#include "sequence/overlap.h"
#include "common/config.h"

int main(int argc, char** argv)
{
    if (argc != 2) return 2;
    Config::load(argv[1]);
    Parameters::get().minimumOverlap = 1000;
    Parameters::get().kmerSize = 17;
    SequenceContainer reads;
    auto a = reads.addSequence(DnaSequence(std::string(10000, 'A')), "a").id;
    auto b = reads.addSequence(DnaSequence(std::string(10000, 'A')), "b").id;
    VertexIndex index(reads);
    OverlapDetector detector(reads, index, 1500, 1000, 1500,
                             false, true, 1.0, false, false, false);
    OverlapContainer cache(detector, reads);
    std::vector<OverlapRange> copied;
    if (cache.copyCachedSeqOverlaps(a, copied) || cache.indexSize()) return 1;
    OverlapRange overlap(a, b, 200, 0, 10000, 10000);
    overlap.curEnd = 9999;
    overlap.extEnd = 9799;
    overlap.score = 123;
    overlap.kmerMatches = new std::vector<std::pair<int32_t, int32_t>>{
        {200, 0}, {9999, 9799}};
    std::vector<OverlapRange> complete{overlap};
    std::atomic<bool> start(false);
    std::vector<std::thread> workers;
    for (int i = 0; i < 8; ++i)
        workers.emplace_back([&]() {
            while (!start.load()) std::this_thread::yield();
            cache.cacheForwardOverlaps(a, complete);
        });
    start.store(true);
    for (auto& worker : workers) worker.join();
    if (cache.indexSize() != 1 || !cache.copyCachedSeqOverlaps(a, copied) ||
        copied.size() != 1 || copied.front().score != 123) return 2;
    copied.front().kmerMatches->front().first = 999;
    // The index is deliberately not built: lazy lookups must reuse the
    // published list, not invoke overlap detection again.
    const auto& forward = cache.lazySeqOverlaps(a);
    const auto& reverse = cache.lazySeqOverlaps(a.rc());
    auto expectedRC = overlap.complement();
    if (forward.size() != 1 || reverse.size() != 1 ||
        reverse.front().curId != expectedRC.curId ||
        reverse.front().extId != expectedRC.extId ||
        reverse.front().curBegin != expectedRC.curBegin ||
        reverse.front().curEnd != expectedRC.curEnd ||
        *reverse.front().kmerMatches != *expectedRC.kmerMatches ||
        *forward.front().kmerMatches != *overlap.kmerMatches) return 3;
    const auto* address = forward.data();
    complete.front().score = 456;
    cache.cacheForwardOverlaps(a, complete);
    if (cache.indexSize() != 1 || forward.data() != address || forward.front().score != 123)
        return 4;
    if (!cache.copyCachedSeqOverlaps(a.rc(), copied) || copied.front().curId != a.rc())
        return 5;
    cache.cacheForwardOverlaps(b, {});
    if (!cache.copyCachedSeqOverlaps(b, copied) || !copied.empty() || cache.indexSize() != 1)
        return 6;
    bool rejectedReverse = false;
    try { cache.cacheForwardOverlaps(a.rc(), {}); }
    catch (const std::runtime_error&) { rejectedReverse = true; }
    if (!rejectedReverse) return 7;
    std::cout << "Overlap-cache reuse, immutable publication, empty results and concurrent writers: PASS\n";
}
