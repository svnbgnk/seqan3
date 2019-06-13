// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::deferred_crtp_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief An invocable wrapper that defers the instantiation of a crtp_base class.
 * \ingroup metafunction
 * \tparam crtp_base The crtp base class to be deferred. Must be a template template parameter.
 * \tparam args_t    A variadic template pack used to augment the `crtp_base` class.
 *
 * \details
 *
 * This metafunction wrapper allows to defer the template instantiation of crtp-base classes. This can be useful
 * if the crtp_base class should be augmented with traits or other templates, especially when using variadic
 * crtp_bases. The help function seqan3::detail::invoke_deferred_crtp_base can be used to instantiate the
 * deferred crtp base with the respective derived type.
 *
 * ### Example
 *
 * The following snippet demonstrates the use of the deferred crtp base class instantiation.
 *
 * \include test/snippet/core/metafunction/deferred_crtp_base.cpp
 *
 * \see seqan3::detail::invoke_deferred_crtp_base
 */
template <template <typename ...> typename crtp_base, typename ...args_t>
struct deferred_crtp_base
{
    /*!\brief Invokes the deferred crtp_base with the corresponding derived type.
     * \tparam derived_t The derived type to instantiate the crtp_base with.
     */
    template <typename derived_t>
    using invoke = crtp_base<derived_t, args_t...>;
};

/*!\brief Template alias to instantiate the deferred crtp base with the derived class.
 * \ingroup metafunction
 * \tparam deferred_crtp_base_t The deferred crtp base class.
 * \tparam derived_t            The derived type to instantiate the crtp base class with.
 *
 * \details
 *
 * Effectively declares the type resulting from `deferred_crtp_base_t::template invoke<derived_t>`.
 *
 * \see seqan3::detail::deferred_crtp_base
 */
template <typename deferred_crtp_base_t, typename derived_t>
//!\cond
    requires requires { typename deferred_crtp_base_t::template invoke<derived_t>; }
//!\endcond
using invoke_deferred_crtp_base = typename deferred_crtp_base_t::template invoke<derived_t>;

} // namespace seqan3::detail