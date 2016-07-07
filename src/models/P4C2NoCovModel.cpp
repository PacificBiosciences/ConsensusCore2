
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

constexpr double kCounterWeight = 1.894736842105264607;
const size_t OUTCOME_NUMBER = 4;
const size_t CONTEXT_NUMBER = 8;

class P4C2NoCovModel : public ModelConfig
{
    REGISTER_MODEL(P4C2NoCovModel);

public:
    static std::set<std::string> Names() { return {"P4-C2"}; }
    P4C2NoCovModel(const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                     const MappedRead& mr, double scoreDiff) const;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double ExpectedLogLikelihoodForMatchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const;
    double ExpectedLogLikelihoodForStickEmission(uint8_t prev, uint8_t curr, bool secondMoment) const;
    double ExpectedLogLikelihoodForBranchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const;
private:
    SNR snr_;
    double cachedEmissionExpectations_[CONTEXT_NUMBER][3][2];
    double ExpectedLogLikelihoodOfOutcomeRow(const int index, const uint8_t prev, const uint8_t curr, const bool secondMoment) const;
};

REGISTER_MODEL_IMPL(P4C2NoCovModel);

// TODO(lhepler) comments regarding the CRTP
class P4C2NoCovRecursor : public Recursor<P4C2NoCovRecursor>
{
public:
    P4C2NoCovRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                      double scoreDiff);
    static inline std::vector<uint8_t> EncodeRead(const MappedRead& read);
    static inline double EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr);
    virtual double UndoCounterWeights(size_t nEmissions) const;
};
constexpr double emissionPmf[3][CONTEXT_NUMBER][OUTCOME_NUMBER] = {
    {// matchPmf
        {     0.99704287,   0.00142953706,    0.0011356778,   0.00028849568},
        {   0.0650431997,     0.930724503,   0.00101193446,   0.00296335123},
        {   0.0452593033,   0.00161481595,     0.951263803,   0.00160737687},
        {   0.0336642993,   0.00171904484,   0.00341177001,     0.961098766},
        {     0.99638081,   0.00169677024,   0.00157124048,  0.000295862976},
        {    0.019536231,      0.97861974,  0.000431533394,   0.00133781704},
        {  0.00200866256,   0.00115526791,     0.994426014,   0.00233404165},
        {  0.00144478772,   0.00311183879,   0.00458133666,     0.990803541}},

    {// branchPmf
        {    0.996618788,  0.000211325763,  0.000211325763,  0.000211325763},
        {  0.00098634885,     0.984298187,  0.000981030928,  0.000981030928},
        { 0.000886336125,  0.000877586924,      0.98594986,  0.000877586924},
        { 0.000134741142,  0.000131163082,  0.000131163082,     0.997897813},
        {    0.999030846,  6.05721377e-05,  6.05721377e-05,  6.05721377e-05},
        { 8.29135034e-05,     0.998673384,  8.29135034e-05,  8.29135034e-05},
        { 0.000191546593,  0.000191546593,     0.996935255,  0.000191546593},
        { 0.000162127506,  0.000162127506,  0.000162127506,      0.99740596}  },

    {// stickPmf
        { 0.000138828593,     0.493085447,     0.312770192,     0.192200761},
        {    0.259663553,  0.000535303321,     0.564413554,     0.168649035},
        {    0.172939367,     0.617816768,  0.000589887862,     0.201334046},
        {    0.209255859,     0.452518724,     0.336452865,  0.000129835235},
        { 8.77384914e-05,     0.542515781,     0.296979203,     0.159276677},
        {    0.289400267,   0.00015477583,     0.472487967,     0.235944905},
        {    0.202007653,     0.561388144,  0.000138893928,     0.234659688},
        {    0.270003409,     0.432586014,     0.296507605,  6.44980403e-05}}};

double P6C4NoCovParams[4][2][3][4] = {
    { // A
     {// NA
    { -3.13042706599906, 0.0928738507007798, -0.0068259313168565, 0.000226837967915887  },
    { -2.77149311404184, -0.0492894388427295, 0.00318374495525272, -0.000104028211974967  },
    { 0.803656427026809, -1.08745490878257, 0.0904491441724511, -0.00249024971528355  } },
     {// AA
    { -2.8862977387348, -0.0604183411870013, 0.00285464541363591, -4.38116647355821e-05  },
    { -2.44610461607789, -0.0388917733637863, -0.00222395923003413, 0.000137857833075156  },
    { 0.0706812431146588, -0.601308221089493, 0.0430115439090431, -0.00100391276309512  } }},
    { // C
     {// NC
    { -3.03668861956175, -0.0110203910326558, 0.00941716538820111, -0.000405787810964155  },
    { -1.53034380578289, -0.486569780370429, 0.0421533848926535, -0.00118412226386447  },
    { 0.471377907209154, -0.905944629261466, 0.0712407677290752, -0.00189737449589616  } },
     {// CC
    { -6.88951775784362, 0.674256382650639, -0.0512910847933769, 0.0013541298330032  },
    { -2.18796852210547, -0.292641361281024, 0.0219538047721538, -0.000474020751617844  },
    { 0.992018689451749, -0.748909802260926, 0.054507350350456, -0.00133801987078912  }}},
    { // G
     {// NG
    { -3.90567763701861, 0.0406355883281611, 0.00391483493872408, -0.000283192473088248  },
    { -2.91097982225128, -0.077806479416406, 0.0072256793887255, -0.000218398644451469  },
    { 1.340853920229, -1.18434278912557, 0.0963594819199681, -0.00267347984061916  } },
     {// GG
    { -2.4408707618026, -0.380970014269318, 0.0317109174018829, -0.000873934344669762  },
    { -3.01731812712857, -0.0344980496324789, -0.00279895538245104, 0.00012815289915934  },
    { 0.542282202259098, -0.632302123281132, 0.0339754460352515, -0.000485581827886687  } }},
    { // T
     {// NT
    { -4.21082767596584, 0.198022950775399, -0.0169634974224188, 0.000393834751568083  },
    { -1.74550760023557, -0.290457036946731, 0.0235405082981186, -0.000546093207126016  },
    { 2.7143463332484, -1.83113271380589, 0.179619418446414, -0.0062680998438969  } },
     {// TT
    { -2.32064310123281, -0.19005586609383, 0.0213592374921381, -0.000624000871181434  },
    { -1.8591597497171, -0.277708318688102, 0.0217401482689649, -0.000513096774869628  },
    { 0.913919454768021, -1.07856774020178, 0.0909136744823634, -0.00261508580897891  } }}};

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

P4C2NoCovModel::P4C2NoCovModel(const SNR& snr)
    : snr_(ClampSNR(snr, SNR(2.18,2.0,2.0,2.0), SNR(22.4, 22.9, 29.9, 26.5)))
{
    for(int ctx = 0; ctx < CONTEXT_NUMBER; ctx++) {
        for (int index = 0; index < 3; index++) {
            cachedEmissionExpectations_[ctx][index][0] = CalculateExpectedLogLikelihoodOfOutcomeRow(index, ctx, false);
            cachedEmissionExpectations_[ctx][index][1] = CalculateExpectedLogLikelihoodOfOutcomeRow(index, ctx, true);
        }
    }

}

std::vector<TemplatePosition> P4C2NoCovModel::Populate(const std::string& tpl) const
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
            1.0 / sum,  // match
            tprobs[0],  // branch
            tprobs[1],  // stick
            tprobs[2]   // deletion
        });

        prev = curr;
    }

    result.emplace_back(TemplatePosition{tpl.back(), prev, 1.0, 0.0, 0.0, 0.0});

    return result;
}

std::unique_ptr<AbstractRecursor> P4C2NoCovModel::CreateRecursor(
    std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff) const
{
    return std::unique_ptr<AbstractRecursor>(
        new P4C2NoCovRecursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff));
}



inline int GetRow(uint8_t prev, uint8_t curr) {
    auto toAdd = prev == curr ? 0 : 4;
    const auto row = curr + toAdd;
    return row;
}

double P4C2NoCovModel::ExpectedLogLikelihoodOfOutcomeRow(const int index, const uint8_t prev, const uint8_t curr, const bool secondMoment) const {
    auto row = GetRow(prev, curr);
    const auto moment = secondMoment ? 1 : 0;
    return cachedEmissionExpectations_[row][index][moment];
}

double P4C2NoCovModel::ExpectedLogLikelihoodForMatchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const {
    return ExpectedLogLikelihoodOfOutcomeRow(static_cast<uint8_t>(MoveType::MATCH), prev, curr, secondMoment);
}
double P4C2NoCovModel::ExpectedLogLikelihoodForStickEmission(uint8_t prev, uint8_t curr, bool secondMoment) const {
    return ExpectedLogLikelihoodOfOutcomeRow(static_cast<uint8_t>(MoveType::STICK), prev, curr, secondMoment);
}
double P4C2NoCovModel::ExpectedLogLikelihoodForBranchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const {
    return ExpectedLogLikelihoodOfOutcomeRow(static_cast<uint8_t>(MoveType::BRANCH), prev, curr, secondMoment);
}    
P4C2NoCovRecursor::P4C2NoCovRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                                     double scoreDiff)
    : Recursor<P4C2NoCovRecursor>(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr,
                                  scoreDiff)
{
}


std::vector<uint8_t> P4C2NoCovRecursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;

    for (const char bp : read.Seq) {
        const uint8_t em = detail::TranslationTable[static_cast<uint8_t>(bp)];
        if (em > 3) throw std::invalid_argument("invalid character in read!");
        result.emplace_back(em);
    }

    return result;
}

double P4C2NoCovRecursor::EmissionPr(MoveType move, const uint8_t emission, const uint8_t prev,
                                     const uint8_t curr)
{
    assert(move != MoveType::DELETION);    assert(move != MoveType::DELETION);
    const auto row = GetRow(prev, curr);
    return emissionPmf[static_cast<uint8_t>(move)][row][emission] * kCounterWeight;
}

double P4C2NoCovRecursor::UndoCounterWeights(const size_t nEmissions) const
{
    return -std::log(kCounterWeight) * nEmissions;
}

}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
