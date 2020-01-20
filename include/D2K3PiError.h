/*
 * Analysis-specific exceptions
 */

#ifndef D2K3PIERROR_H
#define D2K3PIERROR_H

#include <exception>

struct D2K3PiException : public std::exception {
    const char* what() const throw() { return "Base D2K3Pi Exception"; }
};

#endif // D2K3PIERROR_h
