#pragma once

#include <iostream>

void LogInit(int level, std::ostream& out);
void LogSetDbgLevel(int level);

void LogDebug(int level, const char* fmt, ...);
void LogInfo(const char* fmt, ...);
void LogError(const char* fmt, ...);

#define LOG_DEBUG(level, fmt, ...) LogDebug(level, fmt, ##__VA_ARGS__)
#define LOG_INFO(fmt, ...) LogInfo(fmt, ##__VA_ARGS__)
#define LOG_ERROR(fmt, ...) LogError(fmt, ##__VA_ARGS__)
