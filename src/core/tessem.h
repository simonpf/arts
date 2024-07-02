/*!
  \file   tessem.h

  \brief  This file contains functions that are adapted from TESSEM
  code which is used to calculate surface emissivity.
*/

#ifndef tessem_h
#define tessem_h

#include <fstream>

#include "matpack_data.h"

struct TessemNN {
  Index nb_inputs;
  Index nb_outputs;
  Index nb_cache;
  Vector b1;
  Vector b2;
  Matrix w1;
  Matrix w2;
  Vector x_min;
  Vector x_max;
  Vector y_min;
  Vector y_max;

  friend std::ostream& operator<<(std::ostream& os, const TessemNN&) {
    return os;
  }
};

void tessem_read_ascii(std::ifstream& is, TessemNN& net);

void tessem_prop_nn(VectorView ny, const TessemNN& net, ConstVectorView nx);

template <>
struct std::formatter<TessemNN> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const TessemNN& v, FmtContext& ctx) const {
    const std::string_view sep = tags.comma ? ",\n"sv : "\n"sv;

    std::format_to(ctx,
                   "{}{}{}{}{}{}\n",
                   v.nb_inputs,
                   sep,
                   v.nb_outputs,
                   sep,
                   v.nb_cache,
                   sep);
    std::formatter<Vector> vec;
    std::formatter<Matrix> mat;
    make_compat(vec, mat);

    vec.format(v.b1, ctx);
    std::format_to(ctx, "{}", sep);
    vec.format(v.b2, ctx);
    std::format_to(ctx, "{}", sep);
    mat.format(v.w1, ctx);
    std::format_to(ctx, "{}", sep);
    mat.format(v.w2, ctx);
    std::format_to(ctx, "{}", sep);
    vec.format(v.x_min, ctx);
    std::format_to(ctx, "{}", sep);
    vec.format(v.x_max, ctx);
    std::format_to(ctx, "{}", sep);
    vec.format(v.y_min, ctx);
    std::format_to(ctx, "{}", sep);
    vec.format(v.y_max, ctx);
    std::format_to(ctx, "{}", sep);

    return ctx.out();
  }
};

#endif /* tessem_h */
