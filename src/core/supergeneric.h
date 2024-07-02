/*!
  \file   supergeneric.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Thu Jun 13 14:44:02 2002
  
  \brief  Declarations for supergeneric methods.

  The class Any can be used to mark supergenerio IO.
*/

#ifndef supergeneric_h
#define supergeneric_h

#include <format_tags.h>

#include <ostream>

//! A placeholder for any type.
/*!  Used to mark supergeneric methods in file
  methods.cc.
*/
class Any {
  // Nothing to do here.
  friend std::ostream& operator<<(std::ostream& os, const Any&) { return os; }
};

template <>
struct std::formatter<Any> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x._make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Any&, FmtContext& ctx) const {
    return ctx.out();
  }
};

#endif  // supergeneric_h
