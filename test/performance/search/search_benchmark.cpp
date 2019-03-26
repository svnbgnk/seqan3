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
#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

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

//============================================================================
//  undirectional; trivial_search, single, dna4, simulated error = 1
//============================================================================

void unidirectional_search_sim_e1(benchmark::State & state)
{
    uint8_t const searchedErrors_tmp = state.range(0);
    uint8_t const simulated_errors = 1;
    uint32_t index_length = 100000;
    int number_of_Reads = 1;
    uint32_t read_length = 100;
    float probInsertion = 0.18;
    float probDeletion = 0.18;

    auto ref = generate_sequence_seqan3<seqan3::dna4>(index_length, 0, 0);
    fm_index<std::vector<seqan3::dna4>> index{ref};
    std::vector<std::vector<seqan3::dna4> > reads;
    generate_reads(reads, ref, number_of_Reads, read_length, simulated_errors, probInsertion, probDeletion);
    configuration cfg =  search_cfg::max_error{search_cfg::total{searchedErrors_tmp}};

    for (auto _ : state)
    {
        auto results = search(index, reads, cfg);
    }
}

BENCHMARK(unidirectional_search_sim_e1)
    ->Arg(1)
    ->Arg(2)
    ->Arg(3);

//============================================================================
//  undirectional; trivial_search, single, dna4, searched error = 1
//============================================================================

void unidirectional_search_e3(benchmark::State & state)
{
    uint8_t const searchedErrors_tmp = 3;
    uint8_t const simulated_errors = state.range(0);
    uint32_t index_length = 100000;
    int number_of_Reads = 1;
    uint32_t read_length = 100;
    float probInsertion = 0.18;
    float probDeletion = 0.18;

    auto ref = generate_sequence_seqan3<seqan3::dna4>(index_length, 0, 0);
    fm_index<std::vector<seqan3::dna4>> index{ref};
    std::vector<std::vector<seqan3::dna4> > reads;
    generate_reads(reads, ref, number_of_Reads, read_length, simulated_errors, probInsertion, probDeletion);
    configuration cfg =  search_cfg::max_error{search_cfg::total{searchedErrors_tmp}};
    
    for (auto _ : state)
    {
        auto results = search(index, reads, cfg);
    }
}

BENCHMARK(unidirectional_search_e3)
    ->Arg(1)
    ->Arg(2)
    ->Arg(3);
    
// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
