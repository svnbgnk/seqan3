// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include <seqan3/search/all.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#if __has_include(<seqan/align.h>)
    #define SEQAN3_HAS_SEQAN2 1
#endif

#ifdef SEQAN3_HAS_SEQAN2
    #include <seqan/align.h>
    #include <seqan/basic.h>
    #include <seqan/sequence.h>
#endif

using namespace seqan3;



template <typename alphabet_t>
auto generate_sequence_seqan3(size_t const len = 500,
                              size_t const variance = 0,
                              size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, alphabet_size_v<alphabet_t> - 1);
    std::uniform_int_distribution<size_t> dis_length(len - variance, len + variance);

    std::vector<alphabet_t> sequence;

    size_t length = dis_length(gen);
    for (size_t l = 0; l < length; ++l)
        sequence.push_back(alphabet_t{}.assign_rank(dis_alpha(gen)));

    return sequence;
}

template<typename alphabet_t>
void mutateInsertion(std::vector<alphabet_t> & seq, size_t const  errors){
    auto const alphabetSize = alphabet_t::value_size;
    uint32_t rpos = rand() % (seq.size() - errors);
    int rValue = rand() % alphabetSize;
    alphabet_t cbase;
    cbase.assign_rank(rValue);
    seq.insert(seq.begin() + rpos, cbase);
}

template<typename alphabet_t>
void mutateDeletion(std::vector<alphabet_t> & seq, size_t const  errors){
    uint32_t rpos = rand() % (seq.size() - errors);
    seq.erase(seq.begin() + rpos);
}

template<typename alphabet_t>
void mutateSubstitution(std::vector<alphabet_t> & seq, size_t const  errors){
    auto const alphabetSize = alphabet_t::value_size;
    uint32_t rpos = rand() % (seq.size() - errors);
    int rValue = rand() % (alphabetSize - 1);
    alphabet_t & cbase = seq[rpos];
    int cValue = to_rank(cbase);
    if(rValue >=  cValue)
        ++rValue;
    cbase.assign_rank(rValue);
}

template<typename alphabet_t>
void generate_reads(std::vector<std::vector<alphabet_t> > & reads,
                    std::vector<alphabet_t> & ref,
                    size_t const number_of_Reads,
                    size_t const read_length,
                    size_t const simulated_errors,
                    float const probInsertion,
                    float const probDeletion
                   )
{
    for(size_t i = 0; i < number_of_Reads; ++i){
        size_t rpos = rand() % (ref.size() - read_length- simulated_errors);
        std::vector<alphabet_t> query_tmp{ref.begin() + rpos, ref.begin() + rpos + read_length + simulated_errors};
        for(size_t j = 0; j < simulated_errors; ++j)
        {
            //Substitution
            float prob = (float) rand()/RAND_MAX;
            if(probInsertion + probDeletion < prob)
            {
                mutateSubstitution(query_tmp, simulated_errors);
            }
            //Insertion
            else if(probInsertion < prob)
            {
                mutateInsertion(query_tmp, simulated_errors);
            }
            //Deletion
            else
            {
                mutateDeletion(query_tmp, simulated_errors);
            }
        }
        std::vector<alphabet_t>  read{query_tmp.begin(), query_tmp.begin() + read_length};
        reads.push_back(read);
    }
}
/*
#ifdef SEQAN3_HAS_SEQAN2
template <typename alphabet_t>
auto generate_sequence_seqan2(size_t const len = 500,
                              size_t const variance = 0,
                              size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, seqan::ValueSize<alphabet_t>::VALUE - 1);
    std::uniform_int_distribution<size_t> dis_length(len - variance, len + variance);

    seqan::String<alphabet_t> sequence;
    size_t length = dis_length(gen);
    for (size_t l = 0; l < length; ++l)
        appendValue(sequence, alphabet_t{dis_alpha(gen)});

    return sequence;
}
#endif // generate seqan2 data.

// ============================================================================
//  affine; score; dna4; single
// ============================================================================

void seqan3_affine_dna4(benchmark::State & state)
{
    auto cfg = align_cfg::mode{global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}} |
               align_cfg::result{with_score};

    auto seq1 = generate_sequence_seqan3<seqan3::dna4>(500, 0, 0);
    auto seq2 = generate_sequence_seqan3<seqan3::dna4>(500, 0, 1);

    for (auto _ : state)
    {
        auto rng = align_pairwise(std::tie(seq1, seq2), cfg);
        *seqan3::begin(rng);
    }
}

BENCHMARK(seqan3_affine_dna4);

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4(benchmark::State & state)
{
    auto seq1 = generate_sequence_seqan2<seqan::Dna>(500, 0, 0);
    auto seq2 = generate_sequence_seqan2<seqan::Dna>(500, 0, 1);

    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan::globalAlignmentScore(seq1, seq2, seqan::Score<int>{4, -5, -1, -11});
    }
}

BENCHMARK(seqan2_affine_dna4);
#endif // SEQAN3_HAS_SEQAN2*/

// ============================================================================
//  affine; trace; dna4; single
// ============================================================================
/*
void prepareSearch(size_t const searchedErrors, size_t const simulated_errors)
{
    uint32_t index_length = 100000;
    int number_of_Reads = 1;
    uint32_t read_length = 100;
    float probInsertion = 0.18;
    float probDeletion = 0.18;
    
    auto ref = generate_sequence_seqan3<seqan3::dna4>(index_length, 0, 0);
    fm_index<std::vector<seqan3::dna4>> index{ref};
    std::vector<std::vector<seqan3::dna4> > reads;
//     std::vector<std::vector<std::vector<seqan3::dna4> > > collection_reads;
//     collection_reads.resize(simulated_errors + 1);
//     for(int sE = 0; sE <= simulated_errors; ++sE){
//         generate_reads(collection_reads[sE], ref, number_of_Reads, read_length, sE, probInsertion, probDeletion);
//     }
    generate_reads(reads, ref, number_of_Reads, read_length, simulated_errors, probInsertion, probDeletion);
//      std::cout << "Generated reads: " << collection_reads[0].size() << "\n";
//     seqan3::debug_stream << "Generated reads: " << collection_reads[0][0] << "\n";
    configuration const cfg = search_cfg::max_error{search_cfg::total{searchedErrors},
                                                search_cfg::substitution{searchedErrors},
                                                search_cfg::insertion{searchedErrors},
                                                search_cfg::deletion{searchedErrors}};
    
}*/

void seqan3_unidirectional_search(benchmark::State & state/*, size_t const simulated_errors, size_t const searchedErrors*/)
{
//     typedef seqan3::dna4 dna_type;
    uint8_t const simulated_errors = 3;
    uint8_t const searchedErrors = 3;
    uint8_t const searchedErrors_tmp = searchedErrors;
    uint32_t index_length = 100000;
    int number_of_Reads = 1;
    uint32_t read_length = 100;
    float probInsertion = 0.18;
    float probDeletion = 0.18;
    
    auto ref = generate_sequence_seqan3<seqan3::dna4>(index_length, 0, 0);
    fm_index<std::vector<seqan3::dna4>> index{ref};
    std::vector<std::vector<seqan3::dna4> > reads;
    generate_reads(reads, ref, number_of_Reads, read_length, simulated_errors, probInsertion, probDeletion);
    configuration cfg = search_cfg::max_error{search_cfg::total{searchedErrors_tmp},
                                                search_cfg::substitution{searchedErrors_tmp},
                                                search_cfg::insertion{searchedErrors_tmp},
                                                search_cfg::deletion{searchedErrors_tmp}};
//     std::cout << "Searching with " << searchedErrors << "\t" << simulated_errors << "\n";
    for (auto _ : state)
    {
        auto results = search(index, reads, cfg);
    }
    

    
}

BENCHMARK(seqan3_unidirectional_search);

/*
#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4_trace(benchmark::State & state)
{
    auto seq1 = generate_sequence_seqan2<seqan::Dna>(500, 0, 0);
    auto seq2 = generate_sequence_seqan2<seqan::Dna>(500, 0, 1);

    seqan::Gaps<decltype(seq1)> gap1{seq1};
    seqan::Gaps<decltype(seq2)> gap2{seq2};
    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan::globalAlignment(gap1, gap2, seqan::Score<int>{4, -5, -1, -11});
    }
}

BENCHMARK(seqan2_affine_dna4_trace);
#endif // SEQAN3_HAS_SEQAN2*/

// ============================================================================
//  affine; score; dna4; collection
// ============================================================================
/*
void seqan3_affine_dna4_collection(benchmark::State & state)
{
    auto cfg = align_cfg::mode{global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}} |
               align_cfg::result{with_score};

    using sequence_t = decltype(generate_sequence_seqan3<seqan3::dna4>());

    std::vector<std::pair<sequence_t, sequence_t>> vec;
    for (unsigned i = 0; i < 100; ++i)
    {
        sequence_t seq1 = generate_sequence_seqan3<seqan3::dna4>(100, 0, i);
        sequence_t seq2 = generate_sequence_seqan3<seqan3::dna4>(100, 0, i + 100);
        vec.push_back(std::pair{seq1, seq2});
    }

    for (auto _ : state)
    {
        for (auto && rng : align_pairwise(vec, cfg))
            rng.score();
    }
}

BENCHMARK(seqan3_affine_dna4_collection);*/

/*
#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4_collection(benchmark::State & state)
{
    using sequence_t = decltype(generate_sequence_seqan2<seqan::Dna>());

    seqan::StringSet<sequence_t> vec1;
    seqan::StringSet<sequence_t> vec2;
    for (unsigned i = 0; i < 100; ++i)
    {
        sequence_t seq1 = generate_sequence_seqan2<seqan::Dna>(100, 0, i);
        sequence_t seq2 = generate_sequence_seqan2<seqan::Dna>(100, 0, i + 100);
        appendValue(vec1, seq1);
        appendValue(vec2, seq2);
    }

    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan::globalAlignmentScore(vec1, vec2, seqan::Score<int>{4, -5, -1, -11});
    }
}

BENCHMARK(seqan2_affine_dna4_collection);
#endif // SEQAN3_HAS_SEQAN2*/

// ============================================================================
//  affine; trace; dna4; collection
// ============================================================================
/*
void seqan3_affine_dna4_trace_collection(benchmark::State & state)
{
    auto cfg = align_cfg::mode{global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}} |
               align_cfg::result{with_alignment};

    using sequence_t = decltype(generate_sequence_seqan3<seqan3::dna4>());

    std::vector<std::pair<sequence_t, sequence_t>> vec;
    for (unsigned i = 0; i < 100; ++i)
    {
        sequence_t seq1 = generate_sequence_seqan3<seqan3::dna4>(100, 0, i);
        sequence_t seq2 = generate_sequence_seqan3<seqan3::dna4>(100, 0, i + 100);
        vec.push_back(std::pair{seq1, seq2});
    }

    for (auto _ : state)
    {
        for (auto && rng : align_pairwise(vec, cfg))
            rng.score();
    }
}

BENCHMARK(seqan3_affine_dna4_trace_collection);
*/

/*

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4_trace_collection(benchmark::State & state)
{
    using sequence_t = decltype(generate_sequence_seqan2<seqan::Dna>());

    seqan::StringSet<sequence_t> vec1;
    seqan::StringSet<sequence_t> vec2;
    for (unsigned i = 0; i < 100; ++i)
    {
        sequence_t seq1 = generate_sequence_seqan2<seqan::Dna>(100, 0, i);
        sequence_t seq2 = generate_sequence_seqan2<seqan::Dna>(100, 0, i + 100);
        appendValue(vec1, seq1);
        appendValue(vec2, seq2);
    }


    seqan::StringSet<seqan::Gaps<sequence_t>> gap1;
    seqan::StringSet<seqan::Gaps<sequence_t>> gap2;

    for (unsigned i = 0; i < 100; ++i)
    {
        appendValue(gap1, seqan::Gaps<sequence_t>{vec1[i]});
        appendValue(gap2, seqan::Gaps<sequence_t>{vec2[i]});
    }

    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan::globalAlignment(gap1, gap2, seqan::Score<int>{4, -5, -1, -11});
    }
}

BENCHMARK(seqan2_affine_dna4_trace_collection);
#endif // SEQAN3_HAS_SEQAN2

*/

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
