#pragma once

#include <string>
#include <unordered_map>
#include <vector>

class CConfiguration {
public:
    explicit CConfiguration(const std::string& path);

    bool GetParamValueInt(const std::string& key, int& out) const;
    bool GetParamValueDouble(const std::string& key, double& out) const;

    bool GetParamValueIntVec(const std::string& key, std::vector<int>& out, int expected) const;
    bool GetParamValueDoubleVec(const std::string& key, std::vector<double>& out, int expected) const;

    bool GetParamValueDoubleMatrix(const std::string& key, std::vector<std::vector<double>>& out, int n) const;
    bool GetParamVectorMask(const std::string& key, std::vector<int>& out, int expected) const;
    bool GetParamMatrixMask(const std::string& key, std::vector<std::vector<int>>& out, int n) const;

private:
    std::unordered_map<std::string, std::vector<double>> m_values;
};
