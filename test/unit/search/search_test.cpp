// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <type_traits>

#include "helper.hpp"

#include <seqan3/search/algorithm/all.hpp>
#include <gtest/gtest.h>

using namespace seqan3;
using namespace seqan3::search_cfg;
using namespace std::string_literals;

template <typename T>
class search_test : public ::testing::Test
{
public:
    std::vector<dna4> text{"ACGTACGTACGT"_dna4};
    T index{text};
};

template <typename T>
class search_string_test : public ::testing::Test
{
public:
    std::string text{"Garfield the fat cat."};
    T index{text};
};

template <Alphabet alphabet_t>
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

template<Alphabet alphabet_t>
void mutate_insertion(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed = 0){
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, alphabet_size_v<alphabet_t> - 1);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(seq) - overlap);
    alphabet_t cbase;
    seq.insert(seq.begin() + random_pos(gen), alphabet_t{}.assign_rank(dis_alpha(gen)));
}

template<Alphabet alphabet_t>
void mutate_deletion(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed = 0){
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(seq) - overlap);
    seq.erase(seq.begin() + random_pos(gen));
}

template<Alphabet alphabet_t>
void mutate_substitution(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed = 0){
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha_short(0, alphabet_size_v<alphabet_t> - 2);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(seq) - overlap);
    alphabet_t & cbase = seq[random_pos(gen)];
    uint8_t crank = to_rank(cbase);
    uint8_t rrank = dis_alpha_short(gen);
    if (rrank >=  crank)
        ++rrank;
    cbase.assign_rank(rrank);
}

template<Alphabet alphabet_t>
auto generate_reads(std::vector<alphabet_t> & ref,
                    size_t const number_of_reads,
                    size_t const read_length,
                    auto const max_error,
                    float const prob_insertion,
                    float const prob_deletion,
                    size_t const seed = 0)
{
    std::vector<std::vector<alphabet_t> > reads;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> seeds (0, SIZE_MAX);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(ref) - read_length - simulated_errors);
    for (size_t i = 0; i < number_of_reads; ++i){
        size_t rpos = random_pos(gen);
        std::vector<alphabet_t> read_tmp{ref.begin() + rpos,
            ref.begin() + rpos + read_length + max_error.deletion};
        while (max_error.total)
        {
            //Substitution
            float prob = (float) rand()/RAND_MAX;
            if (prob_insertion + prob_deletion < prob)
            {
                mutate_substitution(read_tmp, simulated_errors, seeds(gen));
            }
            //Insertion
            else if (prob_insertion < prob)
            {
                mutate_insertion(read_tmp, simulated_errors, seeds(gen));
            }
            //Deletion
            else
            {
                mutate_deletion(read_tmp, simulated_errors, seeds(gen));
            }
        }
        read_tmp.erase(read_tmp.begin() + read_length, read_tmp.end());
        reads.push_back(read_tmp);
    }
    return reads;
}


using fm_index_types        = ::testing::Types<fm_index<std::vector<dna4>>, bi_fm_index<std::vector<dna4>>>;
using fm_index_string_types = ::testing::Types<fm_index<std::string>, bi_fm_index<std::string>>;

TYPED_TEST_CASE(search_test, fm_index_types);
TYPED_TEST_CASE(search_string_test, fm_index_string_types);

TYPED_TEST(search_test, error_levenshtein)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;
    using namespace search_cfg;
    configuration cfg = max_error{total{3}, deletion{3}, insertion{3}, substitution{3}};
    
    {
        configuration const cfg = max_error{total{1}};
        EXPECT_EQ(uniquify(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }

    {
        configuration const cfg = max_error{total{2}};
        EXPECT_EQ(uniquify(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    }
}

TYPED_TEST(search_test, error_free)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search with empty cfg
        configuration const cfg;
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using empty max_total_error
        configuration const cfg = max_error{};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error
        configuration const cfg = max_error{total{0}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using max_total_error
        configuration const cfg = max_error{total{0}, substitution{0}, insertion{0}, deletion{0}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using empty max_total_error_rate
        configuration const cfg = max_error_rate{};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error_rate
        configuration const cfg = max_error_rate{total{.0}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using max_total_error_rate
        configuration const cfg = max_error_rate{total{.0}, substitution{.0}, insertion{.0}, deletion{.0}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }
}

TYPED_TEST(search_test, multiple_queries)
{
    using hits_result_t = std::vector<std::vector<typename TypeParam::size_type>>;
    std::vector<std::vector<dna4>> const queries{{"GG"_dna4, "ACGTACGTACGT"_dna4, "ACGTA"_dna4}};

    configuration const cfg = max_error_rate{total{.0}, substitution{.0}, insertion{.0}, deletion{.0}};
    EXPECT_EQ(uniquify(search(this->index, queries, cfg)), (hits_result_t{{}, {0}, {0, 4}})); // 0, 1 and 2 hits
}

TYPED_TEST(search_test, invalid_error_configuration)
{
    configuration const cfg = max_error{total{0}, substitution{1}};
    EXPECT_THROW(search(this->index, "A"_dna4, cfg), std::invalid_argument);
}

TYPED_TEST(search_test, error_substitution)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error_rate{total{.25}, substitution{.25}};

        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(uniquify(search(this->index, "CGG"_dna4     , cfg)), (hits_result_t{}));        // not enough mismatches
        EXPECT_EQ(uniquify(search(this->index, "CGTC"_dna4    , cfg)), (hits_result_t{1, 5}));    // 1 mismatch
        EXPECT_EQ(uniquify(search(this->index, "ACGGACG"_dna4 , cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(uniquify(search(this->index, "ACGGACGG"_dna4, cfg)), (hits_result_t{0, 4}));    // 2 mismatches
    }

    {
        configuration const cfg = max_error_rate{total{.25}, substitution{.25}, insertion{.0}, deletion{.0}};

        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(uniquify(search(this->index, "CGG"_dna4     , cfg)), (hits_result_t{}));        // not enough mismatches
        EXPECT_EQ(uniquify(search(this->index, "CGTC"_dna4    , cfg)), (hits_result_t{1, 5}));    // 1 mismatch
        EXPECT_EQ(uniquify(search(this->index, "ACGGACG"_dna4 , cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(uniquify(search(this->index, "ACGGACGG"_dna4, cfg)), (hits_result_t{0, 4}));    // 2 mismatches
    }

    {
        configuration const cfg = max_error{total{1}, substitution{1}};

        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(uniquify(search(this->index, "CGTTT"_dna4   , cfg)), (hits_result_t{}));        // not enough mismatches
        EXPECT_EQ(uniquify(search(this->index, "CGG"_dna4     , cfg)), (hits_result_t{1, 5, 9})); // 1 mismatch
        EXPECT_EQ(uniquify(search(this->index, "ACGGACG"_dna4 , cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(uniquify(search(this->index, "CGTCCGTA"_dna4, cfg)), (hits_result_t{1}));       // 1 mismatch
    }

    {
        configuration const cfg = max_error{total{1}, substitution{1}, insertion{0}, deletion{0}};

        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(uniquify(search(this->index, "CGTTT"_dna4   , cfg)), (hits_result_t{}));        // not enough mismatches
        EXPECT_EQ(uniquify(search(this->index, "CGG"_dna4     , cfg)), (hits_result_t{1, 5, 9})); // 1 mismatch
        EXPECT_EQ(uniquify(search(this->index, "ACGGACG"_dna4 , cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(uniquify(search(this->index, "CGTCCGTA"_dna4, cfg)), (hits_result_t{1}));       // 1 mismatch
    }
}

TYPED_TEST(search_test, error_configuration_types)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        uint8_t s = 1, t = 1;
        configuration const cfg = max_error{total{t}, substitution{s}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
    }

    {
        configuration const cfg = max_error{substitution{1}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
    }
}

TYPED_TEST(search_test, error_insertion)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error_rate{total{.25}, insertion{.25}};

        // exact match and insertion at the beginning of the query
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
        // 1 insertion
        EXPECT_EQ(uniquify(search(this->index, "CCGT"_dna4    , cfg)), (hits_result_t{1, 5, 9}));
        // 2 insertions
        EXPECT_EQ(uniquify(search(this->index, "ACCGGTAC"_dna4, cfg)), (hits_result_t{0, 4}));
        // 2 insertions necessary, only 1 allowed
        EXPECT_EQ(uniquify(search(this->index, "ACCGG"_dna4   , cfg)), (hits_result_t{}));
        // deletion necessary, not allowed
        EXPECT_EQ(uniquify(search(this->index, "ACTACGT"_dna4 , cfg)), (hits_result_t{}));
    }

    {
        configuration const cfg = max_error{total{1}, insertion{1}};

        // exact match and insertion at the beginning of the query
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
        // 1 insertion
        EXPECT_EQ(uniquify(search(this->index, "CCGT"_dna4    , cfg)), (hits_result_t{1, 5, 9}));
        // 2 insertions necessary, only 1 allowed
        EXPECT_EQ(uniquify(search(this->index, "ACCGGTAC"_dna4, cfg)), (hits_result_t{}));
        // deletion necessary, not allowed
        EXPECT_EQ(uniquify(search(this->index, "ACTACGT"_dna4 , cfg)), (hits_result_t{}));
    }
}

TYPED_TEST(search_test, error_deletion)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error_rate{total{.25}, deletion{.25}};

        // exact match, no deletion
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 4, 8}));
        // not enough max errors
        EXPECT_EQ(uniquify(search(this->index, "AGT"_dna4     , cfg)), (hits_result_t{}));
        // one deletion (C)
        EXPECT_EQ(uniquify(search(this->index, "AGTA"_dna4    , cfg)), (hits_result_t{0, 4}));
        // two deletion (C)
        EXPECT_EQ(uniquify(search(this->index, "AGTAGTAC"_dna4, cfg)), (hits_result_t{0}));
        // no deletion at beginning. 0 and 4 cannot be reported
        EXPECT_EQ(uniquify(search(this->index, "CGTACGT"_dna4 , cfg)), (hits_result_t{1, 5}));
    }

    {
        configuration const cfg = max_error{total{1}, deletion{1}};

        // exact match, no deletion
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 4, 8}));
        // one deletion (C)
        EXPECT_EQ(uniquify(search(this->index, "AGTA"_dna4    , cfg)), (hits_result_t{0, 4}));
        // 2 deletions necessary, only 1 allowed
        EXPECT_EQ(uniquify(search(this->index, "AGTAGTAC"_dna4, cfg)), (hits_result_t{}));
        // no deletion at beginning. 0 and 4 cannot be reported
        EXPECT_EQ(uniquify(search(this->index, "CGTACGT"_dna4 , cfg)), (hits_result_t{1, 5}));
    }
}

TYPED_TEST(search_test, error_levenshtein)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error{total{1}};
        EXPECT_EQ(uniquify(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }

    {
        configuration const cfg = max_error{total{2}};
        EXPECT_EQ(uniquify(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    }
}

TYPED_TEST(search_test, error_indel_no_substitution)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error{total{2}, deletion{2}, insertion{2}};
        EXPECT_EQ(uniquify(search(this->index, "GTACCTAC"_dna4, cfg)), (hits_result_t{2}));
    }

    {
        configuration const cfg = max_error{total{3}, deletion{3}, insertion{3}};
        EXPECT_EQ(uniquify(search(this->index, "GTATCCTAC"_dna4, cfg)), (hits_result_t{2}));
    }
}

TYPED_TEST(search_test, search_strategy_all)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error{total{1}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }

    {
        configuration const cfg = max_error{total{1}} | mode{all};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }
}

TYPED_TEST(search_test, search_strategy_best)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error{total{1}} | mode{best};

        hits_result_t possible_hits{0, 4, 8}; // any of 0, 4, 8 ... 1, 5, 9 are not best hits
        hits_result_t result = search(this->index, "ACGT"_dna4, cfg);
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());

        EXPECT_EQ(search(this->index, "AAAA"_dna4, cfg), (hits_result_t{})); // no hit
    }
}

TYPED_TEST(search_test, search_strategy_all_best)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error{total{1}} | mode{all_best};

        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8})); // 1, 5, 9 are not best hits

        EXPECT_EQ(search(this->index, "AAAA"_dna4, cfg), (hits_result_t{})); // no hit
    }
}

TYPED_TEST(search_test, search_strategy_strata)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error{total{1}} | mode{strata{0}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
    }

    {
        configuration const cfg = max_error{total{1}} | mode{strata{1}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }

    {
        configuration const cfg = max_error{total{1}} | mode{strata{1}};
        EXPECT_EQ(search(this->index, "AAAA"_dna4, cfg), (hits_result_t{})); // no hit
    }

    // {
    //     // best hit ACGT with 1 error, i.e. 1+1
    //     configuration const cfg = max_total_error{1} | strategy_strata{1};
    //     EXPECT_EQ(uniquify(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    // }

    // {
    //     // best hit ACGT with 1 error, i.e. 1+1
    //     configuration const cfg = max_total_error{1} | strategy_strata{1};
    //     EXPECT_EQ(uniquify(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    // }
}

TYPED_TEST(search_string_test, error_free_string)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search(this->index, "at"s)), (hits_result_t{14, 18}));
        EXPECT_EQ(uniquify(search(this->index, "Jon"s)), (hits_result_t{}));
    }
}

TYPED_TEST(search_string_test, error_free_raw)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search(this->index, "at")), (hits_result_t{14, 18}));
        EXPECT_EQ(uniquify(search(this->index, "Jon")), (hits_result_t{}));
    }
}

TYPED_TEST(search_string_test, multiple_queries_string)
{
    using hits_result_t = std::vector<std::vector<typename TypeParam::size_type>>;
    std::vector<std::string> const queries{"at", "Jon"};

    EXPECT_EQ(uniquify(search(this->index, queries)), (hits_result_t{{14, 18}, {}})); // 2 and 0 hits
}

TYPED_TEST(search_string_test, multiple_queries_raw)
{
    using hits_result_t = std::vector<std::vector<typename TypeParam::size_type>>;

    EXPECT_EQ(uniquify(search(this->index, {"at", "Jon"})), (hits_result_t{{14, 18}, {}})); // 2 and 0 hits
}

// TYPED_TEST(search_test, return_iterator_index)
// {
// }
//
// TYPED_TEST(search_test, on_hit)
// {
// }
