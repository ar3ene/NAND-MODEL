#pragma once

#include <cstdarg>
#include <cstdio>

inline void LOG_DEBUG(int level, const char* fmt, ...) {
    (void)level;
    va_list ap;
    va_start(ap, fmt);
    std::vfprintf(stderr, fmt, ap);
    std::fprintf(stderr, "\n");
    va_end(ap);
}
