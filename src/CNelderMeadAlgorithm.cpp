#include "CNelderMeadAlgorithm.h"

#include <cmath>
#include <cstdio>

#include "Logger.h"

CNelderMeadAlgorithm::CNelderMeadAlgorithm(int num_params, int num_vertices, int trace_path)
    : m_num_params(num_params), m_num_vertexes(num_vertices), m_trace_path(trace_path) {
}

CNelderMeadAlgorithm::~CNelderMeadAlgorithm() = default;

void CNelderMeadAlgorithm::Init(const std::vector<std::vector<double>>& initial_simplex) {
    m_simplex = initial_simplex;
}

void CNelderMeadAlgorithm::Curvefitting0(FUNCTION function, std::vector<double>& optimal,
                                         int num_iterations, double tolerance, int num_params,
                                         double reflection_coeff, double expansion_coeff, double contraction_coeff) {
    Curvefitting3(function, optimal, num_iterations, tolerance, num_params, reflection_coeff, expansion_coeff, contraction_coeff);
}

void CNelderMeadAlgorithm::Curvefitting1(FUNCTION function, std::vector<double>& optimal,
                                         int num_iterations, double tolerance, int num_params,
                                         double reflection_coeff, double expansion_coeff, double contraction_coeff) {
    Curvefitting3(function, optimal, num_iterations, tolerance, num_params, reflection_coeff, expansion_coeff, contraction_coeff);
}

void CNelderMeadAlgorithm::Curvefitting2(FUNCTION function, std::vector<double>& optimal,
                                         int num_iterations, double tolerance, int num_params,
                                         double reflection_coeff, double expansion_coeff, double contraction_coeff) {
    Curvefitting3(function, optimal, num_iterations, tolerance, num_params, reflection_coeff, expansion_coeff, contraction_coeff);
}

void CNelderMeadAlgorithm::Curvefitting3(FUNCTION function, std::vector<double>& optimal,
                                         int max_iterations, double tolerance, int num_params,
                                         double reflection_coeff, double expansion_coeff, double contraction_coeff) {
    int min_idx = 0;
    int max_idx = 0;

    std::vector<double> f_values(m_num_vertexes, 0);
    Point centroid(num_params, 0);
    Point p_reflected(num_params, 0);
    Point p_expanded(num_params, 0);
    Point p_contracted(num_params, 0);

    ComputeVertexValues(m_simplex, f_values, function);
    min_idx = FindMinIndex(f_values);
    TraceNewMin(m_simplex[min_idx], m_num_params, f_values[min_idx], min_idx, "init");

    for (int iters = 0; iters < max_iterations; ++iters) {
        min_idx = FindMinIndex(f_values);
        max_idx = FindMaxIndex(f_values);

        double fmin = f_values[min_idx];
        double fmax = f_values[max_idx];

        LogInfo("iteration=%d, fmin=%.9g, fmax=%.9g", iters, fmin, fmax);

        if (IsConverged(tolerance)) {
            break;
        }

        ComputeCentroid(m_simplex, centroid, num_params, max_idx);

        Reflect(centroid, m_simplex[max_idx], p_reflected, reflection_coeff);
        double reflected_value = function(p_reflected);

        if (reflected_value < fmin) {
            Expand(centroid, p_expanded, p_reflected, expansion_coeff);
            double expanded_value = function(p_expanded);
            if (expanded_value < reflected_value) {
                m_simplex[max_idx] = p_expanded;
                f_values[max_idx] = expanded_value;
                TraceNewMin(p_expanded, m_num_params, expanded_value, max_idx, "expd");
            } else {
                m_simplex[max_idx] = p_reflected;
                f_values[max_idx] = reflected_value;
                TraceNewMin(p_reflected, m_num_params, reflected_value, max_idx, "refl");
            }
        } else if (reflected_value < f_values[max_idx]) {
            m_simplex[max_idx] = p_reflected;
            f_values[max_idx] = reflected_value;
        } else {
            ContractOneDim(m_simplex, centroid, p_contracted, max_idx, contraction_coeff);
            double contracted_value = function(p_contracted);
            if (contracted_value < f_values[max_idx]) {
                m_simplex[max_idx] = p_contracted;
                f_values[max_idx] = contracted_value;
                TraceNewMin(p_contracted, m_num_params, contracted_value, max_idx, "contr");
            } else {
                // Shrink towards best
                ContractMultiDim(m_simplex, min_idx);
                ComputeVertexValues(m_simplex, f_values, function);
            }
        }
    }

    optimal = m_simplex[min_idx];
}

void CNelderMeadAlgorithm::Reflect(const Point& centroid, const Point& p, Point& p_reflected, double reflection_coeff) {
    p_reflected.resize(m_num_params);
    for (int i = 0; i < m_num_params; ++i) {
        double element = centroid[i] + reflection_coeff * (centroid[i] - p[i]);
        p_reflected[i] = element;
    }
}

void CNelderMeadAlgorithm::Expand(const Point& centroid, Point& p_expanded, const Point& p_reflected, double expansion_coeff) {
    p_expanded.resize(m_num_params);
    for (int i = 0; i < m_num_params; ++i) {
        double element = centroid[i] + expansion_coeff * (p_reflected[i] - centroid[i]);
        p_expanded[i] = element;
    }
}

void CNelderMeadAlgorithm::ContractOneDim(const Points& vertexes_ensemble, const Point& centroid,
                                         Point& p_contracted, int max_idx, double contraction_coeff) {
    p_contracted.resize(m_num_params);
    for (int i = 0; i < m_num_params; ++i) {
        double element = centroid[i] + contraction_coeff * (vertexes_ensemble[max_idx][i] - centroid[i]);
        p_contracted[i] = element;
    }
}

void CNelderMeadAlgorithm::ContractMultiDim(Points& vertexes_ensemble, int min_idx) {
    for (int i = 0; i < m_num_vertexes; ++i) {
        for (int j = 0; j < m_num_params; ++j) {
            vertexes_ensemble[i][j] = 0.5 * (vertexes_ensemble[i][j] + vertexes_ensemble[min_idx][j]);
        }
    }
}

void CNelderMeadAlgorithm::ExpandMultiDim(int min_idx, double expansion_coeff) {
    Point pmin = m_simplex[min_idx];
    for (int i = 0; i < m_num_vertexes; ++i) {
        for (int j = 0; j < m_num_params; ++j) {
            m_simplex[i][j] = m_simplex[i][j] + expansion_coeff * (m_simplex[i][j] - pmin[j]);
        }
    }
}

void CNelderMeadAlgorithm::SimulatedAnnealing() {
    // Placeholder for optional extension.
}

void CNelderMeadAlgorithm::ComputeVertexValues(const Points& vertexes_ensemble, std::vector<double>& f_values, FUNCTION function) {
    for (int i = 0; i < m_num_vertexes; ++i) {
        double f_value = function(vertexes_ensemble[i]);
        f_values[i] = f_value;
    }
}

void CNelderMeadAlgorithm::ComputeCentroid(const Points& vertexes_ensemble, Point& centroid, int num_params, int max_idx) {
    for (int i = 0; i < num_params; ++i) {
        double tmp = 0.0;
        for (int j = 0; j < m_num_vertexes; ++j) {
            tmp += vertexes_ensemble[j][i];
        }
        double element = (tmp - vertexes_ensemble[max_idx][i]) / (m_num_vertexes - 1);
        centroid[i] = element;
    }
}

void CNelderMeadAlgorithm::MiddlePoint(const Point& a, const Point& b, Point& mid) {
    int size = static_cast<int>(a.size());
    mid.resize(size);
    for (int i = 0; i < size; i++) {
        mid[i] = (a[i] + b[i]) * 0.5;
    }
}

void CNelderMeadAlgorithm::RelativeAdd(const Point& o, const Point& a, const Point& b, Point& s) {
    int size = static_cast<int>(o.size());
    s.resize(size);
    for (int i = 0; i < size; i++) {
        s[i] = a[i] + b[i] - o[i];
    }
}

void CNelderMeadAlgorithm::RelativeSub(const Point& o, const Point& a, const Point& b, Point& s) {
    int size = static_cast<int>(o.size());
    s.resize(size);
    for (int i = 0; i < size; i++) {
        s[i] = a[i] - b[i] + o[i];
    }
}

int CNelderMeadAlgorithm::FindMinIndex(const std::vector<double>& values) {
    double min = values[0];
    int min_idx = 0;
    for (unsigned int i = 1; i < values.size(); ++i) {
        if (values[i] < min) {
            min = values[i];
            min_idx = static_cast<int>(i);
        }
    }
    return min_idx;
}

int CNelderMeadAlgorithm::FindMaxIndex(const std::vector<double>& values) {
    double max = values[0];
    int max_idx = 0;
    for (unsigned int i = 1; i < values.size(); ++i) {
        if (values[i] > max) {
            max = values[i];
            max_idx = static_cast<int>(i);
        }
    }
    return max_idx;
}

bool CNelderMeadAlgorithm::IsConverged(double tolerance) {
    bool converged = true;
    for (int i = 0; i < m_num_params; i++) {
        double pmin = m_simplex[0][i];
        double pmax = m_simplex[0][i];
        for (int j = 1; j < m_num_vertexes; j++) {
            double p = m_simplex[j][i];
            if (p < pmin) {
                pmin = p;
            }
            if (p > pmax) {
                pmax = p;
            }
        }
        if ((pmax - pmin) > fabs(pmin) * tolerance) {
            converged = false;
        }
    }
    return converged;
}

void CNelderMeadAlgorithm::TraceNewMin(const Point& p, int num_params, double fmin, int idx, const char* action) {
    if (!m_trace_path) {
        return;
    }

    printf("NEWMIN(%s %d): ", action, idx);
    for (int i = 0; i < num_params; i++) {
        printf("%9g ", p[i]);
    }
    printf("\n");
    (void)fmin;
}

