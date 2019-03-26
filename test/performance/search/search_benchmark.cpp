// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/algorithm/all.hpp>
//
using namespace seqan3;

template<typename alphabet_t>
void generate_text(std::vector<alphabet_t> & text,  unsigned const length)
{
    size_t const alphabet_size = alphabet_t::value_size; 
    text.resize(length);
    for (unsigned i = 0; i < length; ++i)
    {
        alphabet_t r;
        r.assign_rank(std::rand() % alphabet_size);
        text[i] = r;
    }
}

template<typename alphabet_t>
void mutate_insertion(std::vector<alphabet_t> & seq, size_t const  errors){
    size_t const alphabet_size = alphabet_t::value_size;
    size_t rpos = rand() % (seq.size() - errors);
    size_t rrank = rand() % alphabet_size;
    alphabet_t cbase;
    cbase.assign_rank(rrank);
    seq.insert(seq.begin() + rpos, cbase);
}

template<typename alphabet_t>
void mutate_deletion(std::vector<alphabet_t> & seq, size_t const  errors){
    uint32_t rpos = rand() % (seq.size() - errors);
    seq.erase(seq.begin() + rpos);
}

template<typename alphabet_t>
void mutate_substitution(std::vector<alphabet_t> & seq, size_t const  errors){
    size_t const alphabet_size = alphabet_t::value_size;
    size_t rpos = rand() % (seq.size() - errors);
    size_t rrank = rand() % (alphabet_size - 1);
    alphabet_t & cbase = seq[rpos];
    size_t crank = to_rank(cbase);
    if(rrank >=  crank)
        ++rrank;
    cbase.assign_rank(rrank);
}

template<typename alphabet_t>
void generate_reads(std::vector<std::vector<alphabet_t> > & reads,
                    std::vector<alphabet_t> & ref,
                    size_t const number_of_Reads,
                    size_t const read_length,
                    size_t const simulated_errors,
                    float const prob_insertion,
                    float const prob_deletion
                   )
{
    for(size_t i = 0; i < number_of_Reads; ++i){
        size_t rpos = rand() % (ref.size() - read_length - simulated_errors);
        std::vector<alphabet_t> read_tmp{ref.begin() + rpos,
            ref.begin() + rpos + read_length + simulated_errors};
        for(size_t j = 0; j < simulated_errors; ++j)
        {
            //Substitution
            float prob = (float) rand()/RAND_MAX;
            if(prob_insertion + prob_deletion < prob)
            {
                mutate_substitution(read_tmp, simulated_errors);
            }
            //Insertion
            else if(prob_insertion < prob)
            {
                mutate_insertion(read_tmp, simulated_errors);
            }
            //Deletion
            else
            {
                mutate_deletion(read_tmp, simulated_errors);
            }
        }
        std::vector<alphabet_t>  read{read_tmp.begin(), read_tmp.begin() + read_length};
        reads.push_back(read);
    }
}

//============================================================================
//  undirectional; trivial_search, single, dna4, searched_error = 3
//============================================================================

void unidirectional_search(benchmark::State & state)
{
    uint8_t const simulated_errors = state.range(0);
    uint8_t const searched_errors = state.range(1);
    uint32_t index_length = 100000;
    int number_of_reads = 500;
    uint32_t read_length = 100;
    float prob_insertion = 0.18;
    float prob_deletion = 0.18;

    std::vector<seqan3::dna4> ref;
    generate_text(ref, index_length);
    fm_index<std::vector<seqan3::dna4> > index{ref};
    std::vector<std::vector<seqan3::dna4> > reads;
    generate_reads(reads, ref, number_of_reads, read_length, simulated_errors, prob_insertion, prob_deletion);
    configuration cfg =  search_cfg::max_error{search_cfg::total{searched_errors}};
    
    for (auto _ : state)
    {
        auto results = search(index, reads, cfg);
    }
}

BENCHMARK(unidirectional_search)
    ->Args({0, 1})
    ->Args({1, 1})
    ->Args({0, 3})
    ->Args({1, 3})
    ->Args({2, 3})
    ->Args({3, 3});
    
// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
