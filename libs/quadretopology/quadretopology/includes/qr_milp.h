/***************************************************************************/
/* Copyright(C) 2021


The authors of

Reliable Feature-Line Driven Quad-Remeshing
Siggraph 2021


 All rights reserved.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****************************************************************************/

/* Small MILP modelling layer on top of HiGHS (MIT), replacing the Gurobi C++ API.
 * HiGHS does not solve mixed-integer quadratic problems, so:
 * - |x| terms become an epigraph variable with two linear inequalities (exact, since they are
 *   only minimized with positive weights);
 * - squared terms (a*x + b)^2 of an integer variable x are replaced by their piecewise linear
 *   interpolation at integer breakpoints (exact at every integer within MILP_SQUARE_EXACT_RANGE of
 *   the minimum, coarser breakpoints further away). */

#ifndef QR_MILP_H
#define QR_MILP_H

#include <Highs.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#define MILP_INFINITY std::numeric_limits<double>::infinity()
#define MILP_SQUARE_EXACT_RANGE 16

namespace QuadRetopology {
namespace internal {

class MILPException : public std::runtime_error
{
public:
    explicit MILPException(const std::string& message) : std::runtime_error(message) {}
};

struct MILPVar {
    int id = -1;
};

class MILPLinExpr
{
public:
    MILPLinExpr(double constant = 0.0) : constant(constant) {}
    MILPLinExpr(int constant) : constant(constant) {}
    MILPLinExpr(const MILPVar& var) : terms{{var.id, 1.0}}, constant(0.0) {}

    MILPLinExpr& operator+=(const MILPLinExpr& other) {
        terms.insert(terms.end(), other.terms.begin(), other.terms.end());
        constant += other.constant;
        return *this;
    }
    MILPLinExpr& operator-=(const MILPLinExpr& other) {
        for (const std::pair<int, double>& t : other.terms)
            terms.emplace_back(t.first, -t.second);
        constant -= other.constant;
        return *this;
    }
    MILPLinExpr& operator*=(const double factor) {
        for (std::pair<int, double>& t : terms)
            t.second *= factor;
        constant *= factor;
        return *this;
    }
    MILPLinExpr& operator/=(const double factor) {
        return *this *= (1.0 / factor);
    }

    //Merge duplicated variables and drop zero coefficients
    void compact() {
        std::sort(terms.begin(), terms.end());
        std::vector<std::pair<int, double>> merged;
        for (const std::pair<int, double>& t : terms) {
            if (!merged.empty() && merged.back().first == t.first)
                merged.back().second += t.second;
            else
                merged.push_back(t);
        }
        merged.erase(std::remove_if(merged.begin(), merged.end(),
                                    [](const std::pair<int, double>& t) { return t.second == 0.0; }),
                     merged.end());
        terms = std::move(merged);
    }

    std::vector<std::pair<int, double>> terms;
    double constant;
};

inline MILPLinExpr operator+(MILPLinExpr a, const MILPLinExpr& b) { return a += b; }
inline MILPLinExpr operator-(MILPLinExpr a, const MILPLinExpr& b) { return a -= b; }
inline MILPLinExpr operator-(MILPLinExpr a) { return a *= -1.0; }
inline MILPLinExpr operator*(MILPLinExpr a, const double b) { return a *= b; }
inline MILPLinExpr operator*(const double a, MILPLinExpr b) { return b *= a; }
inline MILPLinExpr operator/(MILPLinExpr a, const double b) { return a /= b; }

//Constraint expr in [lower, upper]
struct MILPConstr {
    MILPLinExpr expr;
    double lower;
    double upper;
};

inline MILPConstr operator<=(const MILPLinExpr& a, const MILPLinExpr& b) { return {a - b, -MILP_INFINITY, 0.0}; }
inline MILPConstr operator>=(const MILPLinExpr& a, const MILPLinExpr& b) { return {a - b, 0.0, MILP_INFINITY}; }
inline MILPConstr operator==(const MILPLinExpr& a, const MILPLinExpr& b) { return {a - b, 0.0, 0.0}; }

enum class MILPResult { OPTIMAL, FEASIBLE, NOSOLUTION };

class MILPModel
{
public:
    MILPModel() : timeLimit(MILP_INFINITY), gapLimit(0.0), verbose(true), mipGap(MILP_INFINITY) {}

    MILPVar addVar(const double lower, const double upper, const bool integer) {
        MILPVar var;
        var.id = static_cast<int>(colLower.size());
        colLower.push_back(lower);
        colUpper.push_back(upper);
        colInteger.push_back(integer);
        startValue.push_back(kHighsUndefined);
        return var;
    }

    void addConstr(MILPConstr constr) {
        constr.expr.compact();
        //The constant is moved to the bounds
        const double lower = constr.lower - constr.expr.constant;
        const double upper = constr.upper - constr.expr.constant;
        if (constr.expr.terms.empty()) {
            if (lower > 1e-9 || upper < -1e-9)
                throw MILPException("Infeasible constant constraint");
            return;
        }
        rowLower.push_back(lower);
        rowUpper.push_back(upper);
        rowStart.push_back(static_cast<HighsInt>(rowIndex.size()));
        for (const std::pair<int, double>& t : constr.expr.terms) {
            rowIndex.push_back(t.first);
            rowValue.push_back(t.second);
        }
    }

    //Returns a variable equal to |value| at the optimum (value only minimized with positive weight).
    //As with Gurobi's addGenConstrAbs on an integer variable diff == value, an integer abs also
    //constrains value to be integer.
    MILPLinExpr addAbs(const MILPLinExpr& value, const bool integer) {
        MILPLinExpr diff = value;
        if (integer) {
            MILPVar diffVar = addVar(-MILP_INFINITY, MILP_INFINITY, true);
            addConstr(MILPLinExpr(diffVar) == value);
            diff = diffVar;
        }
        MILPVar abs = addVar(0, MILP_INFINITY, integer);
        addConstr(MILPLinExpr(abs) >= diff);
        addConstr(MILPLinExpr(abs) >= -diff);
        return abs;
    }

    //Returns a variable equal to value^2 (at integer values of value) at the optimum.
    //value must be integer valued, or a*x + b with x an integer variable.
    MILPLinExpr addSquare(MILPLinExpr value) {
        value.compact();
        if (value.terms.empty())
            return value.constant * value.constant;

        MILPVar x;
        double a = 1.0;
        double b = 0.0;
        if (value.terms.size() == 1 && colInteger[value.terms[0].first]) {
            x.id = value.terms[0].first;
            a = value.terms[0].second;
            b = value.constant;
        }
        else {
            x = addVar(-MILP_INFINITY, MILP_INFINITY, true);
            addConstr(MILPLinExpr(x) == value);
        }

        //f(x) = (a*x + b)^2 is convex and is replaced by its piecewise linear interpolation at
        //integer breakpoints b_0 < ... < b_n (incremental formulation): x = b_0 + sum_i d_i with
        //0 <= d_i <= b_(i+1) - b_i and cost f(b_0) + sum_i slope_i * d_i. The slopes increase, so
        //at the optimum the segments are filled in order and the cost equals the interpolation.
        const double center = -b / a;
        const double lower = colLower[x.id];
        const double upper = colUpper[x.id];
        const double range = MILP_SQUARE_EXACT_RANGE + std::fabs(center);

        std::vector<double> breakpoints;
        for (double k = std::floor(center) - MILP_SQUARE_EXACT_RANGE; k <= std::ceil(center) + MILP_SQUARE_EXACT_RANGE; k++)
            breakpoints.push_back(k);
        for (double step = 2; step <= 4 * range; step *= 2) {
            breakpoints.push_back(std::ceil(center) + MILP_SQUARE_EXACT_RANGE + step);
            breakpoints.push_back(std::floor(center) - MILP_SQUARE_EXACT_RANGE - step);
        }
        breakpoints.push_back(lower);
        breakpoints.push_back(upper);
        breakpoints.erase(std::remove_if(breakpoints.begin(), breakpoints.end(),
                                         [&](const double k) { return !std::isfinite(k) || k < lower || k > upper; }),
                          breakpoints.end());
        std::sort(breakpoints.begin(), breakpoints.end());
        breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end()), breakpoints.end());

        auto f = [&](const double k) { return (a * k + b) * (a * k + b); };

        const double first = breakpoints.front();
        const double last = breakpoints.back();

        MILPLinExpr square = f(first);
        MILPLinExpr position = first;
        for (size_t i = 0; i + 1 < breakpoints.size(); i++) {
            const double k0 = breakpoints[i];
            const double k1 = breakpoints[i + 1];
            MILPVar d = addVar(0, k1 - k0, false);
            square += MILPLinExpr(d) * ((f(k1) - f(k0)) / (k1 - k0));
            position += d;
        }
        //Unbounded variable: the outermost segments are extended with unit slopes
        if (!std::isfinite(lower)) {
            MILPVar d = addVar(0, MILP_INFINITY, false);
            square += MILPLinExpr(d) * (f(first - 1) - f(first));
            position -= d;
        }
        if (!std::isfinite(upper)) {
            MILPVar d = addVar(0, MILP_INFINITY, false);
            square += MILPLinExpr(d) * (f(last + 1) - f(last));
            position += d;
        }
        addConstr(MILPLinExpr(x) == position);

        return square;
    }

    void setObjective(MILPLinExpr obj) {
        obj.compact();
        objective = std::move(obj);
    }

    //MIP start: variables without a start value are completed by HiGHS
    void setStart(const MILPVar& var, const double value) {
        startValue[var.id] = value;
    }

    //Early stop rule: stop when the gap is below gaps[i] after times[i] seconds
    void setEarlyStop(const std::vector<float>& times, const std::vector<float>& gaps) {
        stopTimes = times;
        stopGaps = gaps;
    }

    MILPResult optimize() {
        const HighsInt numCol = static_cast<HighsInt>(colLower.size());
        const HighsInt numRow = static_cast<HighsInt>(rowLower.size());

        HighsLp lp;
        lp.num_col_ = numCol;
        lp.num_row_ = numRow;
        lp.sense_ = ObjSense::kMinimize;
        lp.offset_ = objective.constant;
        lp.col_cost_.assign(numCol, 0.0);
        for (const std::pair<int, double>& t : objective.terms)
            lp.col_cost_[t.first] = t.second;
        lp.col_lower_ = colLower;
        lp.col_upper_ = colUpper;
        lp.row_lower_ = rowLower;
        lp.row_upper_ = rowUpper;
        lp.a_matrix_.format_ = MatrixFormat::kRowwise;
        lp.a_matrix_.num_col_ = numCol;
        lp.a_matrix_.num_row_ = numRow;
        lp.a_matrix_.start_ = rowStart;
        lp.a_matrix_.start_.push_back(static_cast<HighsInt>(rowIndex.size()));
        lp.a_matrix_.index_ = rowIndex;
        lp.a_matrix_.value_ = rowValue;
        lp.integrality_.resize(numCol);
        for (HighsInt i = 0; i < numCol; i++)
            lp.integrality_[i] = colInteger[i] ? HighsVarType::kInteger : HighsVarType::kContinuous;

        highs.setOptionValue("output_flag", verbose);
        highs.setOptionValue("time_limit", std::isfinite(timeLimit) ? timeLimit : kHighsInf);
        highs.setOptionValue("mip_rel_gap", gapLimit);

        if (highs.passModel(lp) == HighsStatus::kError)
            throw MILPException("HiGHS: invalid model");

        std::vector<HighsInt> startIndex;
        std::vector<double> startValues;
        for (HighsInt i = 0; i < numCol; i++) {
            if (startValue[i] != kHighsUndefined) {
                startIndex.push_back(i);
                startValues.push_back(startValue[i]);
            }
        }
        if (!startIndex.empty())
            highs.setSolution(static_cast<HighsInt>(startIndex.size()), startIndex.data(), startValues.data());

        size_t stopIndex = 0;
        bool stopped = false;
        if (!stopTimes.empty() && stopTimes.size() == stopGaps.size()) {
            highs.setCallback([this, &stopIndex, &stopped](int callbackType, const std::string&,
                                                 const HighsCallbackOutput* dataOut,
                                                 HighsCallbackInput* dataIn, void*) {
                if (callbackType != kCallbackMipInterrupt || stopped)
                    return;

                while (stopIndex < stopTimes.size() && stopTimes[stopIndex] <= dataOut->running_time)
                    stopIndex++;

                if (stopIndex > 0) {
                    const double objbst = dataOut->mip_primal_bound;
                    const double objbnd = dataOut->mip_dual_bound;
                    if (!std::isfinite(objbst))
                        return;

                    const double currentGap = std::fabs(objbnd - objbst) / std::fabs(objbst);
                    if (currentGap < stopGaps[stopIndex - 1]) {
                        std::cout << "Stop early - Gap: " << currentGap << " < " << stopGaps[stopIndex - 1]
                                  << " gap achieved after " << stopTimes[stopIndex - 1] << " seconds" << std::endl;
                        dataIn->user_interrupt = true;
                        stopped = true;
                    }
                }
            }, nullptr);
            highs.startCallback(kCallbackMipInterrupt);
        }

        const HighsStatus runStatus = highs.run();
        if (runStatus == HighsStatus::kError)
            throw MILPException("HiGHS: error during the optimization");

        const HighsInfo& info = highs.getInfo();
        mipGap = info.mip_gap;

        if (info.primal_solution_status != kSolutionStatusFeasible) {
            throw MILPException("HiGHS: no feasible solution found (" +
                                highs.modelStatusToString(highs.getModelStatus()) + ")");
        }

        solution = highs.getSolution().col_value;
        return highs.getModelStatus() == HighsModelStatus::kOptimal ? MILPResult::OPTIMAL : MILPResult::FEASIBLE;
    }

    double getValue(const MILPVar& var) const {
        return solution.at(var.id);
    }

    double getValue(const MILPLinExpr& expr) const {
        double value = expr.constant;
        for (const std::pair<int, double>& t : expr.terms)
            value += t.second * solution.at(t.first);
        return value;
    }

    double getMIPGap() const {
        return mipGap;
    }

    double timeLimit;
    double gapLimit;
    bool verbose;

private:
    Highs highs;

    std::vector<double> colLower;
    std::vector<double> colUpper;
    std::vector<bool> colInteger;
    std::vector<double> startValue;

    std::vector<double> rowLower;
    std::vector<double> rowUpper;
    std::vector<HighsInt> rowStart;
    std::vector<HighsInt> rowIndex;
    std::vector<double> rowValue;

    MILPLinExpr objective;

    std::vector<float> stopTimes;
    std::vector<float> stopGaps;

    std::vector<double> solution;
    double mipGap;
};

}
}

#endif // QR_MILP_H
