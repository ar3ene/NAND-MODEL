#include "CConfiguration.h"

#include <fstream>
#include <sstream>
#include <cctype>

static std::string Trim(const std::string& s) {
    size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
        ++start;
    }
    size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
        --end;
    }
    return s.substr(start, end - start);
}

CConfiguration::CConfiguration(const std::string& path) {
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        // Remove comments
        size_t hash = line.find('#');
        if (hash != std::string::npos) {
            line = line.substr(0, hash);
        }
        line = Trim(line);
        if (line.empty()) {
            continue;
        }
        // Expect key = values
        size_t eq = line.find('=');
        if (eq == std::string::npos) {
            continue;
        }
        std::string key = Trim(line.substr(0, eq));
        std::string rhs = Trim(line.substr(eq + 1));
        if (key.empty() || rhs.empty()) {
            continue;
        }
        std::istringstream iss(rhs);
        std::vector<double> vals;
        double v;
        while (iss >> v) {
            vals.push_back(v);
            // Skip commas if present
            if (iss.peek() == ',') {
                iss.get();
            }
        }
        if (!vals.empty()) {
            m_values[key] = vals;
        }
    }
}

bool CConfiguration::GetParamValueInt(const std::string& key, int& out) const {
    auto it = m_values.find(key);
    if (it == m_values.end() || it->second.empty()) {
        return false;
    }
    out = static_cast<int>(it->second[0]);
    return true;
}

bool CConfiguration::GetParamValueDouble(const std::string& key, double& out) const {
    auto it = m_values.find(key);
    if (it == m_values.end() || it->second.empty()) {
        return false;
    }
    out = it->second[0];
    return true;
}

bool CConfiguration::GetParamValueIntVec(const std::string& key, std::vector<int>& out, int expected) const {
    auto it = m_values.find(key);
    if (it == m_values.end()) {
        return false;
    }
    out.clear();
    for (double v : it->second) {
        out.push_back(static_cast<int>(v));
    }
    if (expected > 0 && static_cast<int>(out.size()) < expected) {
        return false;
    }
    return true;
}

bool CConfiguration::GetParamValueDoubleVec(const std::string& key, std::vector<double>& out, int expected) const {
    auto it = m_values.find(key);
    if (it == m_values.end()) {
        return false;
    }
    out = it->second;
    if (expected > 0 && static_cast<int>(out.size()) < expected) {
        return false;
    }
    return true;
}

bool CConfiguration::GetParamValueDoubleMatrix(const std::string& key, std::vector<std::vector<double>>& out, int n) const {
    auto it = m_values.find(key);
    if (it == m_values.end()) {
        return false;
    }
    const std::vector<double>& vals = it->second;
    if (n <= 0) {
        return false;
    }
    if (static_cast<int>(vals.size()) < n * n) {
        return false;
    }
    out.assign(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            out[i][j] = vals[i * n + j];
        }
    }
    return true;
}

bool CConfiguration::GetParamVectorMask(const std::string& key, std::vector<int>& out, int expected) const {
    auto it = m_values.find(key);
    if (it == m_values.end()) {
        return false;
    }
    out.clear();
    for (double v : it->second) {
        out.push_back(static_cast<int>(v));
    }
    if (expected > 0 && static_cast<int>(out.size()) < expected) {
        return false;
    }
    return true;
}

bool CConfiguration::GetParamMatrixMask(const std::string& key, std::vector<std::vector<int>>& out, int n) const {
    auto it = m_values.find(key);
    if (it == m_values.end()) {
        return false;
    }
    const std::vector<double>& vals = it->second;
    if (n <= 0) {
        return false;
    }
    if (static_cast<int>(vals.size()) < n * n) {
        return false;
    }
    out.assign(n, std::vector<int>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            out[i][j] = static_cast<int>(vals[i * n + j]);
        }
    }
    return true;
}
