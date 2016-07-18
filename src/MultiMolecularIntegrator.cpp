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
#include <limits>
#include <utility>

#include <pacbio/consensus/MultiMolecularIntegrator.h>
#include <pacbio/consensus/Sequence.h>

#include "ModelFactory.h"

#ifdef DEBUG_BAM_OUTPUT
#include <fstream>
#include <cram/md5.h>
#include <pbbam/BamWriter.h>
#include <pbbam/BamFile.h>
#include <pbbam/MD5.h>
#include <pacbio/consensus/align/LinearAlignment.h>
#endif

namespace PacBio {
namespace Consensus {

MultiMolecularIntegrator::MultiMolecularIntegrator(const std::string& tpl,
                                                   const IntegratorConfig& cfg)
    : AbstractIntegrator(cfg), fwdTpl_{tpl}, revTpl_{::PacBio::Consensus::ReverseComplement(tpl)}
{
}

State MultiMolecularIntegrator::AddRead(const MappedRead& read)
{
    try {
        return AbstractIntegrator::AddRead(GetTemplate(read, read.SignalToNoise), read);
    } catch (const TemplateTooSmall& e) {
        return State::TEMPLATE_TOO_SMALL;
    }
}

size_t MultiMolecularIntegrator::TemplateLength() const { return fwdTpl_.length(); }
char MultiMolecularIntegrator::operator[](const size_t i) const { return fwdTpl_[i]; }
MultiMolecularIntegrator::operator std::string() const { return fwdTpl_; }
void MultiMolecularIntegrator::ApplyMutation(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));

    std::vector<Mutation> fwdMuts = {fwdMut};
    std::vector<Mutation> revMuts = {revMut};

    fwdTpl_ = ::PacBio::Consensus::ApplyMutations(fwdTpl_, &fwdMuts);
    revTpl_ = ::PacBio::Consensus::ApplyMutations(revTpl_, &revMuts);

    for (auto& eval : evals_) {
        if (eval.Strand() == StrandType::FORWARD)
            eval.ApplyMutation(fwdMut);
        else if (eval.Strand() == StrandType::REVERSE)
            eval.ApplyMutation(revMut);
    }

    assert(fwdTpl_.length() == revTpl_.length());
    assert(fwdTpl_ == ::PacBio::Consensus::ReverseComplement(revTpl_));
}

void MultiMolecularIntegrator::ApplyMutations(std::vector<Mutation>* fwdMuts)
{
    std::vector<Mutation> revMuts;

    for (auto it = fwdMuts->crbegin(); it != fwdMuts->crend(); ++it)
        revMuts.emplace_back(ReverseComplement(*it));

    fwdTpl_ = ::PacBio::Consensus::ApplyMutations(fwdTpl_, fwdMuts);
    revTpl_ = ::PacBio::Consensus::ApplyMutations(revTpl_, &revMuts);

    for (auto& eval : evals_) {
        if (eval.Strand() == StrandType::FORWARD)
            eval.ApplyMutations(fwdMuts);
        else if (eval.Strand() == StrandType::REVERSE)
            eval.ApplyMutations(&revMuts);
    }

    assert(fwdTpl_.length() == revTpl_.length());
    assert(fwdTpl_ == ::PacBio::Consensus::ReverseComplement(revTpl_));
}

std::unique_ptr<AbstractTemplate> MultiMolecularIntegrator::GetTemplate(const MappedRead& read,
                                                                        const SNR& snr)
{
    const size_t len = read.TemplateEnd - read.TemplateStart;

    if (read.Strand == StrandType::FORWARD) {
        const size_t start = read.TemplateStart;
        const size_t end = read.TemplateEnd;

        return std::unique_ptr<AbstractTemplate>(
            new Template(fwdTpl_.substr(start, len), ModelFactory::Create(read.Model, snr), start,
                         end, read.PinStart, read.PinEnd));
    } else if (read.Strand == StrandType::REVERSE) {
        const size_t start = revTpl_.size() - read.TemplateEnd;
        const size_t end = revTpl_.size() - read.TemplateStart;

        return std::unique_ptr<AbstractTemplate>(
            new Template(revTpl_.substr(start, len), ModelFactory::Create(read.Model, snr), start,
                         end, read.PinEnd, read.PinStart));
    }

    throw std::invalid_argument("read is unmapped!");
}

#ifdef DEBUG_BAM_OUTPUT
std::string AlignmentCigar(const std::string& tpl, const std::string& query)
{
    PairwiseAlignment* p = AlignLinear(tpl, query);
    return p->Cigar();
}

size_t TemplateSpan(const MappedRead& read)
{
    return read.TemplateEnd - read.TemplateStart;
}

int32_t ReadToFlag(const MappedRead& read)
{
    uint32_t flag = 0;
    if (read.Strand == StrandType::REVERSE)
        flag |= static_cast<uint32_t>(16);
    return flag;
}

std::string ReadToSequence(const MappedRead& read)
{
    if (read.Strand == StrandType::REVERSE) {
        return ::PacBio::Consensus::ReverseComplement(read.Seq);
    } else {
        return read.Seq;
    }
}

std::string MultiMolecularIntegrator::ReadToCigar(const MappedRead& read) const
{
    const size_t tlen = read.TemplateEnd - read.TemplateStart;
    const std::string tpl = std::string(fwdTpl_).substr(read.TemplateStart, tlen);
    if (read.Strand == StrandType::REVERSE) {
        const std::string rcSeq = ::PacBio::Consensus::ReverseComplement(read.Seq);
        return AlignmentCigar(tpl, rcSeq);
    } else {
        return AlignmentCigar(tpl, read.Seq);
    }
}

std::vector<PacBio::BAM::BamRecordImpl> MultiMolecularIntegrator::ToBamRecords(const int32_t rid) const
{
    using namespace PacBio::BAM;
    std::vector<BamRecordImpl> records;
    for (auto& eval : evals_) {
        // Create a BamRecord for every (A) valid eval with (B) a valid mapped read
        if (eval) {
            if (const auto& read = eval.Read()) {
                records.emplace_back();
                records.back().Name(read->Name);                                   // QNAME
                records.back().Flag(ReadToFlag(*read));                            // FLAG
                records.back().ReferenceId(rid);                                   // RNAME
                records.back().Position(read->TemplateStart);                      // POS
                records.back().CigarData(ReadToCigar(*read));                      // CIGAR
                records.back().InsertSize(TemplateSpan(*read));                    // TLEN
                records.back().SetSequenceAndQualities(ReadToSequence(*read));     // SEQ
                records.back().AddTag("ip", Tag{read->IPD});                       // TAG "ip" for IPD
                records.back().AddTag("mo", Tag{read->Model});                     // TAG "mo" for sequencing model
                records.back().AddTag("pw", Tag{read->PulseWidth});                // TAG "pw" for PulseWidth
                records.back().AddTag("sn", Tag{read->SignalToNoise.ToVector()});  // TAG "sn" for SNR
            }
        }
    }

    return records;
}

void MultiMolecularIntegrator::WriteBamFile(const std::string& filepath,
                                            const std::string& name,
                                            const bool sorted,
                                            const bool indexed) const
{
    using namespace PacBio::BAM;

    BamHeader header;
    // Add SQ entry with the template information
    SequenceInfo si(name, std::to_string(fwdTpl_.size()));
    si.Checksum(MD5Hash(fwdTpl_));
    header.AddSequence(si);
    // Add CO entry with the raw sequence
    header.AddComment(name + "\t" + std::string(fwdTpl_));

    {   // Blocked to ensure flushing and completion of file I/O
        BamWriter writer(filepath, header);
        const int32_t sid = header.SequenceId(name);
        std::vector<BamRecordImpl> records = ToBamRecords(sid);
    
        if (sorted)
            std::stable_sort(records.begin(), records.end(),
                [](const BamRecordImpl& lhs, const BamRecordImpl& rhs) { 
                    return lhs.Position() < rhs.Position(); 
                });

        for (const auto& record : records)
            writer.Write(record);
    }

    if (indexed) {
        BamFile f(filepath);
        f.CreateStandardIndex();
    }
}

void MultiMolecularIntegrator::WriteReferenceFasta(const std::string& filepath, 
                                                   const std::string& name) const
{
    const size_t COLUMN_SIZE = 60;
    std::string tpl(fwdTpl_);
    std::ofstream fastaFile;
    fastaFile.open(filepath);
    fastaFile << ">" << name << std::endl;
    for (size_t i = 0; i < tpl.size();) {
        fastaFile << tpl.substr(i, COLUMN_SIZE) << std::endl;
        i += COLUMN_SIZE;
    }
    fastaFile.close();
}

void MultiMolecularIntegrator::WriteAlignments(const std::string& filePrefix) const
{
    WriteReferenceFasta(filePrefix + ".fasta");
    WriteBamFile(filePrefix + ".bam");
}

#endif

}  // namespace Consensus
}  // namespace PacBio
