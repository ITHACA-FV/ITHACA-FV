#include <ostream>

/// \file
/// Simple header and source file of the Color::Modifier class to change color to the output stream.

#ifndef colormod_H
#define colormod_H

namespace Color
{
enum Code
{
    FG_RED      = 31,
    FG_GREEN    = 32,
    FG_BLUE     = 34,
    FG_DEFAULT  = 39,
    BG_RED      = 41,
    BG_GREEN    = 42,
    BG_BLUE     = 44,
    BG_DEFAULT  = 49
};
/// Class to change color to the output stream
class Modifier
{
        Code code;
    public:
        Modifier(Code pCode) : code(pCode) {}
        friend std::ostream&
        operator<<(std::ostream& os, const Modifier& mod)
        {
            return os << "\033[" << mod.code << "m";
        }
};
}

#endif
