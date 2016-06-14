
#pragma once

#include <tuple>
#include <vector>

#include <pacbio/consensus/Mutation.h>

namespace PacBio {
namespace Consensus {

// forward declaration
class AbstractIntegrator;

struct PolishConfig
{
    size_t MaximumIterations;
    size_t MutationSeparation;
    size_t MutationNeighborhood;

    PolishConfig(size_t iterations = 40, size_t separation = 10, size_t neighborhood = 20);
};

std::tuple<bool, size_t, size_t> Polish(AbstractIntegrator* ai, const PolishConfig& cfg);

struct QualityValues
{
    std::vector<int> Qualities;
    std::vector<int> DeletionQVs;
    std::vector<int> InsertionQVs;
    std::vector<int> SubstitutionQVs;
};

std::vector<int> ConsensusQualities(AbstractIntegrator& ai);

QualityValues ConsensusQVs(AbstractIntegrator& ai);

std::vector<Mutation> Mutations(const AbstractIntegrator& ai);

}  // namespace Consensus
}  // namespace PacBio
