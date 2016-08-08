// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/consensus/Read.h>

#include "../ModelFactory.h"
#include "../Recursor.h"

namespace PacBio {
namespace Consensus {
namespace {
    
constexpr size_t OUTCOME_NUMBER = 4;
constexpr size_t CONTEXT_NUMBER = 8;

constexpr double kCounterWeight = 2.0;

class P6C4NoCovModel : public ModelConfig
{
    REGISTER_MODEL(P6C4NoCovModel);

public:
    static std::set<std::string> Names() { return {"P6-C4::Snr"}; }
    P6C4NoCovModel(const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                     const MappedRead& mr, double scoreDiff) const;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double ExpectedLLForEmission(MoveType move, uint8_t prev, uint8_t curr,
                                 MomentType moment) const;


private:
    SNR snr_;
    double cachedEmissionExpectations_[CONTEXT_NUMBER][3][2];

};

REGISTER_MODEL_IMPL(P6C4NoCovModel);

// TODO(lhepler) comments regarding the CRTP
class P6C4NoCovRecursor : public Recursor<P6C4NoCovRecursor>
{
public:
    P6C4NoCovRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                      double scoreDiff);
    static inline std::vector<uint8_t> EncodeRead(const MappedRead& read);
    static inline double EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr);
    virtual double UndoCounterWeights(size_t nEmissions) const;
private:
    double ctxTrans_[CONTEXT_NUMBER][4];
   
};

double P6C4NoCovParams[4][2][3][4] = {
    { // A
     {// NA
         { -3.16995059942652, 0.124783365029361, -0.00690196141925598, 0.000177995862249767  },
         { -4.01783083321503, 0.318063447634267, -0.0300729852540775, 0.0007831389331307  },
         { 0.109068191364145, -0.57698988101572, 0.0301752645136841, -0.000543103717822943  } },
     {// AA
         { -3.48395040742402, 0.0126593051973385, 0.00325840287267268, -0.000107077000106521  },
         { -4.56345047094783, 0.499490687180124, -0.0439352725030146, 0.00105119070507155  },
         { -0.859831905262017, -0.243301363980332, 0.00733992534660296, 4.26157389357674e-05  } }},
    { // C
     {// NC
         { -4.48381290560576, 0.594014860084145, -0.0734335219550843, 0.00284890591211796  },
         { -0.534933599351947, -0.650210598148757, 0.049903690079979, -0.00120375227890548  },
         { 2.47066086087946, -1.6802237943306, 0.156887347129726, -0.00496467932771181  } },
     {// CC
         { -10.3097185116682, 1.79562668245579, -0.192673135145183, 0.00706787328507119  },
         { -2.56400604276783, 0.104651248402936, -0.0442677508098647, 0.00254383834736708  },
         { 1.63597793376433, -1.0906295165667, 0.0964704509493554, -0.0029393576406705  } }},
    { // G
     {// NG
         { -3.52685181915902, -0.116760910611688, 0.0105665414951042, 1.77090807527055e-05  },
         { -3.69177532118875, 0.40918026835223, -0.0751139249916928, 0.00339016991998504  },
         { 1.2466774545486, -1.16043019260752, 0.0947547076789284, -0.00279983324959412  }},
     {// GG
         { -2.20496179118271, -0.617206594653536, 0.068601552864775, -0.00237925284887096  },
         { -3.21187089890749, 0.134544816905904, -0.0423057091608044, 0.00218649088584648  },
         { 0.503909617791076, -0.761260981497719, 0.0622724688170413, -0.00181947652105094  } }},
    { // T
     {// NT
         { -5.03617456615196, 0.596997287662421, -0.0729281692099805, 0.00292661830349639  },
         { -2.76376571774453, 0.224388055841335, -0.0661786464522392, 0.00346615080359996  },
         { 3.18648305314085, -1.87512264762453, 0.193090866939792, -0.00741463951401265  } },
     {// TT
         { -2.94466862205134, -0.0254074517613081, 0.00627679565298507, -0.000111876999992817  },
         { -2.56927108557671, 0.151325602973565, -0.0597420172307808, 0.00315569090189143  },
         { 0.854896636773977, -0.830591794425688, 0.059047379681373, -0.00147406054452713  } }}};
    
    
    constexpr double emissionPmf[3][CONTEXT_NUMBER][OUTCOME_NUMBER] = {
    // AA CC GG TT NA NC NG NT
    {   // Match PMF
        
        {    0.993857288,   0.00494559296,  0.000385300035,  0.000575618025},
        {   0.0465109567,     0.952012096,  0.000600314495,  0.000666174201},
        {   0.0411554441,  0.000544930781,     0.956760221,   0.00133578784},
        {   0.0388533948,  0.000575700341,   0.00166993904,     0.958664439},
        {    0.992865151,   0.00656866703,  0.000164777704,  0.000322755116},
        {  0.00433138659,     0.994100935,  0.000396422343,   0.00110404528},
        {  0.00119537112,  0.000440240277,     0.997277117,   0.00102056588},
        { 0.000972007313,   0.00072722415,   0.00173750926,     0.996481506}
    },
        
    {   // Branch PMF
        
        {    0.993899765,  0.000381264688,  0.000381264688,  0.000381264688},
        {  0.00159160579,     0.974555317,   0.00159020515,   0.00159020515},
        { 0.000790343196,   0.00078946654,     0.987367659,   0.00078946654},
        { 0.000302541831,  0.000302233939,  0.000302233939,     0.995163949},
        {    0.999034758,   6.0327653e-05,   6.0327653e-05,   6.0327653e-05},
        { 0.000107011829,     0.998287811,  0.000107011829,  0.000107011829},
        { 0.000203215502,  0.000203215502,     0.996748552,  0.000203215502},
        { 0.000194014393,  0.000194014393,  0.000194014393,      0.99689577}
        },

    { // Stick PMF
        
        { 0.000429870922,     0.306480035,     0.332319633,     0.355182139},
        {    0.263877239,  0.000359740517,     0.401227098,     0.329967236},
        {    0.236686656,     0.398616081,  0.000573262408,     0.356968582},
        {    0.308995996,     0.317470903,     0.366711811,  0.000509608509},
        { 0.000148241682,     0.358717055,     0.318864058,     0.320343504},
        {    0.262894962,  9.96324496e-05,     0.382794886,     0.352915298},
        {    0.314950856,     0.301650343,   0.00014629655,     0.381350649},
        {    0.327107321,     0.332843922,     0.337905548,  0.000153086352}
    }};
    
inline double CalculateExpectedLogLikelihoodOfOutcomeRow(const int index, const uint8_t row, const bool secondMoment)  {
    double expectedLL = 0;
    for(size_t i = 0; i < OUTCOME_NUMBER; i++) {
        double curProb = emissionPmf[index][row][i];
        double lgCurProb = std::log(curProb);
        if(!secondMoment) {
            expectedLL +=  curProb * lgCurProb;
        } else {
            expectedLL += curProb * pow(lgCurProb, 2.0);
        }
    }
    return expectedLL;
}

// For P6-C4 we cap SNR at 20.0 (19.0 for C); as the training set only went that
// high; extrapolation beyond this cap goes haywire because of the higher-order
// terms in the regression model.  See bug 31423.
P6C4NoCovModel::P6C4NoCovModel(const SNR& snr)
    : snr_(ClampSNR(snr, SNR(2.87, 2.0, 2.1, 2.28), SNR(24, 15, 16.3, 15.5)))
{
    for(int ctx = 0; ctx < CONTEXT_NUMBER; ctx++) {
        for (int index = 0; index < 3; index++) {
            cachedEmissionExpectations_[ctx][index][0] = CalculateExpectedLogLikelihoodOfOutcomeRow(index, ctx, false);
            cachedEmissionExpectations_[ctx][index][1] = CalculateExpectedLogLikelihoodOfOutcomeRow(index, ctx, true);
        }
    }
}
    



std::vector<TemplatePosition> P6C4NoCovModel::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    uint8_t prev = detail::TranslationTable[static_cast<uint8_t>(tpl[0])];
    if (prev > 3) throw std::invalid_argument("invalid character in sequence!");

    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t curr = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];
        if (curr > 3) throw std::invalid_argument("invalid character in sequence!");
        const bool hp = tpl[i - 1] == tpl[i];  // NA -> 0, AA -> 1
        const auto params = P6C4NoCovParams[curr][hp];
        const double snr = snr_[curr], snr2 = snr * snr, snr3 = snr2 * snr;
        double tprobs[3];
        double sum = 1.0;

        for (size_t j = 0; j < 3; ++j) {
            double xb =
                params[j][0] + snr * params[j][1] + snr2 * params[j][2] + snr3 * params[j][3];
            xb = std::exp(xb);
            tprobs[j] = xb;
            sum += xb;
        }

        for (size_t j = 0; j < 3; ++j)
            tprobs[j] /= sum;

        result.emplace_back(TemplatePosition{
            tpl[i - 1], prev,
            tprobs[1],  // match
            1.0 / sum,  // branch
            tprobs[2],  // stick
            tprobs[0]   // deletion
        });

        prev = curr;
    }

    result.emplace_back(TemplatePosition{tpl.back(), prev, 1.0, 0.0, 0.0, 0.0});

    return result;
}

std::unique_ptr<AbstractRecursor> P6C4NoCovModel::CreateRecursor(
    std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff) const
{
    return std::unique_ptr<AbstractRecursor>(
        new P6C4NoCovRecursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff));
}

inline int GetRow(uint8_t prev, uint8_t curr) {
        auto toAdd = prev == curr ? 0 : 4;
        const auto row = curr + toAdd;
        return row;
}

double P6C4NoCovModel::ExpectedLLForEmission(const MoveType move, const uint8_t prev,
                                                 const uint8_t curr, const MomentType moment) const

{
        auto row = GetRow(prev, curr);
        return cachedEmissionExpectations_[row][static_cast<uint8_t>(move)][static_cast<uint8_t>(moment)];
}
    

P6C4NoCovRecursor::P6C4NoCovRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                                     double scoreDiff)
    : Recursor<P6C4NoCovRecursor>(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr,
                                  scoreDiff)
{
}

std::vector<uint8_t> P6C4NoCovRecursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;

    for (const char bp : read.Seq) {
        const uint8_t em = detail::TranslationTable[static_cast<uint8_t>(bp)];
        if (em > 3) throw std::invalid_argument("invalid character in read!");
        result.emplace_back(em);
    }

    return result;
}

    
double P6C4NoCovRecursor::EmissionPr(MoveType move, const uint8_t emission, const uint8_t prev,
                                     const uint8_t curr)
{
        assert(move != MoveType::DELETION);
        const auto row = GetRow(prev, curr);
        return emissionPmf[static_cast<uint8_t>(move)][row][emission] * kCounterWeight;
}

double P6C4NoCovRecursor::UndoCounterWeights(const size_t nEmissions) const
{
    return -std::log(kCounterWeight) * nEmissions;
}

}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
