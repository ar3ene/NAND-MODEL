#ifndef _NELDER_MEAD_ALGORITHM_H
#define _NELDER_MEAD_ALGORITHM_H

#include <iostream>
#include <vector>
#include <functional>

#include "SysTypes.h"

using Point = std::vector<double>;
using Points = std::vector<Point>;
using FUNCTION = std::function<double(std::vector<double>&)>;

class CNelderMeadAlgorithm {
public:
    CNelderMeadAlgorithm(int num_params, int num_vertices, int trace_path);
    ~CNelderMeadAlgorithm();

    void Init(const std::vector<std::vector<double>>& initial_simplex);

    void Curvefitting0(FUNCTION function, std::vector<double>& optimal,
                       int num_iterations, double tolerance, int num_params,
                       double reflection_coeff, double expansion_coeff, double contraction_coeff);

    void Curvefitting1(FUNCTION function, std::vector<double>& optimal,
                       int num_iterations, double tolerance, int num_params,
                       double reflection_coeff, double expansion_coeff, double contraction_coeff);

    void Curvefitting2(FUNCTION function, std::vector<double>& optimal,
                       int num_iterations, double tolerance, int num_params,
                       double reflection_coeff, double expansion_coeff, double contraction_coeff);

    void Curvefitting3(FUNCTION function, std::vector<double>& optimal,
                       int num_iterations, double tolerance, int num_params,
                       double reflection_coeff, double expansion_coeff, double contraction_coeff);

private:
    void Reflect(const Point& centroid, const Point& p, Point& p_reflected, double reflection_coeff);
    void Expand(const Point& centroid, Point& p_expanded, const Point& p_reflected, double expansion_coeff);

    void ContractOneDim(const Points& vertexes_ensemble, const Point& centroid,
                        Point& p_contracted, int max_idx, double contraction_coeff);
    void ContractMultiDim(Points& vertexes_ensemble, int min_idx);
    void ExpandMultiDim(int min_idx, double expansion_coeff);

    void SimulatedAnnealing();

    void ComputeVertexValues(const Points& vertexes_ensemble, std::vector<double>& f_values, FUNCTION function);
    void ComputeCentroid(const Points& vertexes_ensemble, Point& centroid, int num_params, int max_idx);

    void MiddlePoint(const Point& a, const Point& b, Point& mid);
    void RelativeAdd(const Point& o, const Point& a, const Point& b, Point& s);
    void RelativeSub(const Point& o, const Point& a, const Point& b, Point& s);

    int FindMinIndex(const std::vector<double>& values);
    int FindMaxIndex(const std::vector<double>& values);

    bool IsConverged(double tolerance);
    void TraceNewMin(const Point& p, int num_params, double fmin, int idx, const char* info);

    Points m_simplex;

    int m_num_params = 0;
    int m_num_vertexes = 0;
    int m_trace_path = 0;
};

#endif

