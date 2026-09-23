#include <iostream>
#include <random>
#include "sequence/sequence_container.h"
#include "common/config.h"
void removeContainedDisjointigs(std::vector<FastaRecord>&, float);

std::string dna(std::mt19937& rng, size_t size)
{
    std::string result(size, 'A');
    for (auto& c : result) c = "ACGT"[rng() % 4];
    return result;
}

bool check(const std::vector<DnaSequence>& sequences, size_t expected)
{
    SequenceContainer records;
    std::vector<FastaRecord> disjointigs;
    for (size_t i = 0; i < sequences.size(); ++i)
        disjointigs.push_back(records.addSequence(sequences[i], "seq" + std::to_string(i)));
    removeContainedDisjointigs(disjointigs, 0.1f);
    bool ok = disjointigs.size() == expected;
    if (!ok) std::cerr << "FAIL inputs=" << sequences.size()
                       << " expected=" << expected << " retained=" << disjointigs.size() << '\n';
    return ok;
}

int main(int argc, char** argv)
{
    if (argc != 2) return 2;
    Config::load(argv[1]);
    Parameters::get().numThreads = 1;
    Parameters::get().minimumOverlap = 1000;
    Parameters::get().kmerSize = 17;
    std::mt19937 rng(42);
    std::string seq = dna(rng, 10000);
    DnaSequence a(seq), b(dna(rng, 10000)), shortSeq(seq.substr(2000, 6000));
    int failures = 0;
    failures += !check({a, a}, 1);
    failures += !check({a, a, a}, 1);
    failures += !check({a, a.complement()}, 1);
    failures += !check({a.complement(), a}, 1);
    failures += !check({a, b}, 2);
    failures += !check({a, shortSeq}, 1);
    failures += !check({shortSeq, a}, 1);
    failures += !check({a}, 1);
    std::cout << "Contained-disjointig cases: 8; failures: " << failures << '\n';
    return failures ? 1 : 0;
}
