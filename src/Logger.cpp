#include "Logger.h"

#include <cstdarg>
#include <cstdio>

static bool s_logger_initialized = false;
static int g_debug_level = 10;
static std::ostream* g_out = &std::cerr;

void LogInit(int level, std::ostream& out) {
    g_debug_level = level;
    g_out = &out;
    s_logger_initialized = true;
}

void LogSetDbgLevel(int level) {
    g_debug_level = level;
}

static void vlogf(const char* tag, int level, const char* fmt, va_list ap) {
    if (level > g_debug_level) {
        return;
    }
    std::fprintf(stderr, "[%s] ", tag);
    std::vfprintf(stderr, fmt, ap);
    std::fprintf(stderr, "\n");
    if (g_out) {
        (*g_out) << ""; // keep stream in sync
    }
}

void LogDebug(int level, const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vlogf("DEBUG", level, fmt, ap);
    va_end(ap);
}

void LogInfo(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vlogf("INFO", 1, fmt, ap);
    va_end(ap);
}

void LogError(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vlogf("ERROR", 0, fmt, ap);
    va_end(ap);
}

