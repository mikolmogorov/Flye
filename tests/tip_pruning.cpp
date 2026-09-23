#include <iostream>
#include "repeat_graph/multiplicity_inferer.h"
#include "common/config.h"

// Exercise the real pruning routine on paired-strand graphs. Swapping the
// numerical IDs of an edge and its complement cannot change its survival.
bool check(bool flip, bool interior, int length, int cov, bool onlyTips)
{
    SequenceContainer assembly, sequences, reads;
    RepeatGraph graph(assembly, &sequences);
    GraphNode* n[12];
    for (auto& node : n) node = graph.addNode();
    auto edge = [&](int left, int right, int id, int coverage, int len) {
        GraphEdge e(n[left], n[right], FastaRecord::Id(id));
        e.meanCoverage = coverage;
        e.seqSegments.emplace_back(FastaRecord::ID_NONE, len);
        graph.addEdge(std::move(e));
    };
    edge(0, 1, 2, 20, 2000000); edge(5, 4, 3, 20, 2000000);
    edge(1, 2, 4, 20, 2000000); edge(6, 5, 5, 20, 2000000);
    edge(1, 3, flip ? 1 : 0, cov, length);
    edge(7, 5, flip ? 0 : 1, cov, length);
    if (interior)
    {
        // Add a second high-coverage backbone at the other end of the
        // candidate, turning it from a tip into an internal branch.
        edge(8, 3, 6, 20, 2000000); edge(7, 10, 7, 20, 2000000);
        edge(3, 9, 8, 20, 2000000); edge(11, 7, 9, 20, 2000000);
    }
    ReadAligner aligner(graph, reads);
    MultiplicityInferer inferer(graph, aligner, assembly);
    int removed = inferer.removeUnsupportedEdges(onlyTips);
    bool expectedRemoval = cov < 3 && (!onlyTips || !interior);
    bool retained = graph.getEdge(FastaRecord::Id(0)) != nullptr;
    bool rcRetained = graph.getEdge(FastaRecord::Id(1)) != nullptr;
    bool ok = retained == !expectedRemoval && rcRetained == retained &&
              removed == int(expectedRemoval);
    if (!ok)
        std::cerr << "FAIL flip=" << flip << " interior=" << interior
                  << " length=" << length << " cov=" << cov
                  << " onlyTips=" << onlyTips << " removed=" << removed
                  << " retained=" << retained << '/' << rcRetained << '\n';
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
    for (bool flip : {false, true})
        for (bool interior : {false, true})
            for (int length : {10000, 1200000})
                for (int cov : {1, 3, 20})
                    for (bool onlyTips : {false, true})
                        failures += !check(flip, interior, length, cov, onlyTips);
    std::cout << "Tip-pruning cases: 48; failures: " << failures << '\n';
    return failures ? 1 : 0;
}
