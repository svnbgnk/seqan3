#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/repeat.hpp>
#include <seqan3/std/ranges>

int main()
{
    auto v = seqan3::view::repeat('A');

    seqan3::debug_stream << *std::ranges::begin(v) << '\n'; // prints 'A'
    seqan3::debug_stream << v[12355] << '\n';               // also prints 'A'. It always prints 'A'

    v[1345] = 'C';

    // Now it always prints 'C'
    seqan3::debug_stream << *std::ranges::begin(v) << '\n'; // prints 'C'
    seqan3::debug_stream << v[12355] << '\n';               // prints 'C'
}
