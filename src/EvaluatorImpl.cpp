
#include <cmath>

#include <boost/optional.hpp>

#include <pacbio/consensus/Evaluator.h>

#include "EvaluatorImpl.h"

namespace PacBio {
namespace Consensus {
namespace {  // anonymous

constexpr size_t EXTEND_BUFFER_COLUMNS = 8;

#if 0
std::ostream& operator<<(std::ostream& out, const std::pair<size_t, size_t>& x)
{
    return out << '(' << x.first << ", " << x.second << ')';
}

void WriteMatrix(const ScaledMatrix& mat)
{
    std::cerr << std::pair<size_t, size_t>(mat.Rows(), mat.Columns()) << std::endl;

    for (size_t j = 0; j < mat.Columns(); ++j)
        std::cerr << " " << mat.UsedRowRange(j);
    std::cerr << std::endl;

    std::cerr << "lg: ";
    for (size_t j = 0; j < mat.Columns(); ++j)
        std::cerr << "\t" << std::fixed << std::setprecision(3) << mat.GetLogScale(j);
    std::cerr << std::endl;

    std::cerr << "lgS: ";
    double lgS = 0.0;
    for (size_t j = 0; j < mat.Columns(); ++j)
        std::cerr << "\t" << std::fixed << std::setprecision(3) << (lgS += mat.GetLogScale(j));
    std::cerr << std::endl;

    for (size_t i = 0; i < mat.Rows(); ++i)
    {
        for (size_t j = 0; j < mat.Columns(); ++j)
        {
            std::cerr << "\t" << std::fixed << std::setprecision(3) << std::log(mat.Get(i, j)) + mat.GetLogScale(j);
        }
        std::cerr << std::endl;
    }
}
#endif

}  // namespace anonymous

EvaluatorImpl::EvaluatorImpl(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                             const double scoreDiff)
    : recursor_(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff)
    , alpha_(mr.Length() + 1, recursor_.tpl_->Length() + 1)
    , beta_(mr.Length() + 1, recursor_.tpl_->Length() + 1)
    , extendBuffer_(mr.Length() + 1, EXTEND_BUFFER_COLUMNS)
{
    ValidateAssumptions();
    recursor_.FillAlphaBeta(alpha_, beta_);
    if (!std::isfinite(LL())) throw AlphaBetaMismatch();
}
    
void EvaluatorImpl::ValidateAssumptions() {
    if (TemplateSpan() == 0) {
        throw std::runtime_error("Cannot score a read that aligns to an empty template sequence");
    }
}

double EvaluatorImpl::LL(const Mutation& mut_)
{
    boost::optional<Mutation> mut(recursor_.tpl_->Mutate(mut_));

    // apply the virtual mutation
    if (!mut) return LL();

    size_t betaLinkCol = 1 + mut->End();
    size_t absoluteLinkColumn = 1 + mut->End() + mut->LengthDiff();

    double score;

    bool atBegin = mut->Start() < 3;
    bool atEnd = (mut->End() + 3) > beta_.Columns();

    if (!atBegin && !atEnd) {
        int extendStartCol, extendLength = 2;

        if (mut->Type == MutationType::DELETION) {
            // Future thought: If we revise the semantic of Extra,
            // we can remove the extend and just link alpha and
            // beta directly.
            extendStartCol = mut->Start() - 1;
        } else {
            extendStartCol = mut->Start();
            assert(extendLength <= EXTEND_BUFFER_COLUMNS);
        }

        recursor_.ExtendAlpha(alpha_, extendStartCol, extendBuffer_, extendLength);
        score = recursor_.LinkAlphaBeta(extendBuffer_, extendLength, beta_, betaLinkCol,
                                        absoluteLinkColumn) +
                alpha_.GetLogProdScales(0, extendStartCol);
    } else if (!atBegin && atEnd) {
        //
        // Extend alpha to end
        //
        size_t extendStartCol = mut->Start() - 1;
        assert(recursor_.tpl_->Length() + 1 > extendStartCol);
        size_t extendLength = recursor_.tpl_->Length() - extendStartCol + 1;

        recursor_.ExtendAlpha(alpha_, extendStartCol, extendBuffer_, extendLength);
        score = std::log(extendBuffer_(recursor_.read_.Length(), extendLength - 1)) +
                alpha_.GetLogProdScales(0, extendStartCol) +
                extendBuffer_.GetLogProdScales(0, extendLength);
    } else if (atBegin && !atEnd) {
        // If the mutation occurs at positions 0 - 2
        size_t extendLastCol = mut->End();
        // We duplicate this math inside the function
        size_t extendLength = 1 + mut->End() + mut->LengthDiff();

        recursor_.ExtendBeta(beta_, extendLastCol, extendBuffer_, mut->LengthDiff());
        score = std::log(extendBuffer_(0, 0)) +
                beta_.GetLogProdScales(extendLastCol + 1, beta_.Columns()) +
                extendBuffer_.GetLogProdScales(0, extendLength);
    } else {
        assert(atBegin && atEnd);
        /* This should basically never happen...
           and is a total disaster if it does.  The basic idea is that
           FillAlpha and FillBeta use the real "template" while we test
           mutations using "virtual" template positions and the Extend/Link
           methods.  Trying to call FillAlpha to calculate the likelihood of a
        virtual
           mutation is therefore going to fail, as it calculates using the
           "real" template.
        throw TooSmallTemplateException();
         */

        /*
         Just do the whole fill, but guard against the case where we have deleted 
         away the final base in the template that this aligned to (leads to 
         nasty SIGSEGV */
        if (TemplateSpan() == 0) {
            recursor_.tpl_->Reset();
            return 0;
        } else {
            ScaledMatrix alphaP(recursor_.read_.Length() + 1, recursor_.tpl_->Length() + 1);
            recursor_.FillAlpha(ScaledMatrix::Null(), alphaP);
            score = std::log(alphaP(recursor_.read_.Length(), recursor_.tpl_->Length())) +
                alphaP.GetLogProdScales();
        }
    }

    // reset the virtual mutation
    recursor_.tpl_->Reset();

    return score + recursor_.tpl_->UndoCounterWeights(recursor_.read_.Length());
}

double EvaluatorImpl::LL() const
{
    return std::log(beta_(0, 0)) + beta_.GetLogProdScales() +
           recursor_.tpl_->UndoCounterWeights(recursor_.read_.Length());
}

size_t EvaluatorImpl::TemplateSpan() const {
    return recursor_.tpl_->Length();
}

std::pair<double, double> EvaluatorImpl::NormalParameters() const
{
    return recursor_.tpl_->NormalParameters();
}

double EvaluatorImpl::ZScore() const
{
    double mean, var;
    std::tie(mean, var) = NormalParameters();
    return (LL() - mean) / std::sqrt(var);
}

void EvaluatorImpl::Recalculate()
{
    ValidateAssumptions();
    size_t I = recursor_.read_.Length() + 1;
    size_t J = recursor_.tpl_->Length() + 1;
    alpha_.Reset(I, J);
    beta_.Reset(I, J);
    extendBuffer_.Reset(I, EXTEND_BUFFER_COLUMNS);
    recursor_.FillAlphaBeta(alpha_, beta_);
}

void EvaluatorImpl::ApplyMutation(const Mutation& mut)
{
    recursor_.tpl_->ApplyMutation(mut);
}

void EvaluatorImpl::ApplyMutations(std::vector<Mutation>* muts)
{
    recursor_.tpl_->ApplyMutations(muts);
}

/*
Matrix* AlphaMatrix();
Matrix* BetaMatrix();
*/

}  // namespace Consensus
}  // namespace PacBio
