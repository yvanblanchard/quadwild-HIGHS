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

#include "qr_ilp.h"

#include "qr_milp.h"

#define MIN_SUBDIVISION_VALUE 1
#define FEASIBILITY_FIX_COST 1000000.0

namespace QuadRetopology {
namespace internal {

void getChartSubsideSum(
        const ChartData& chartData,
        const size_t& cId,
        const std::vector<MILPVar>& vars,
        const std::vector<bool>& isFixed,
        const std::vector<int>& ilpResults,
        const bool hardParityConstraint,
        std::vector<MILPLinExpr>& chartSubsideSum);

void getChartSubsideSumResults(
        const ChartData& chartData,
        const size_t& cId,
        const std::vector<int>& results,
        std::vector<int>& chartSubsideSum);

MILPLinExpr getMILPCostTermInteger(MILPModel& model, const ILPMethod& method, const MILPLinExpr& value);
MILPLinExpr getMILPCostTermContinuous(MILPModel& model, const ILPMethod& method, const MILPLinExpr& value);

MILPLinExpr getMILPAbsInteger(MILPModel& model, const MILPLinExpr& value);
MILPLinExpr getMILPAbsContinuous(MILPModel& model, const MILPLinExpr& value);

inline void solveILP(
        const ChartData& chartData,
        const std::vector<double>& chartEdgeLength,
        const ILPMethod& method,
        const double alpha,
        const bool isometry,
        const bool regularityQuadrilaterals,
        const bool regularityNonQuadrilaterals,
        const double regularityNonQuadrilateralsWeight,
        const bool alignSingularities,
        const double alignSingularitiesWeight,
        const int repeatLosingConstraintsIterations,
        const bool repeatLosingConstraintsQuads,
        const bool repeatLosingConstraintsNonQuads,
        const bool repeatLosingConstraintsAlign,
        const bool feasibilityFix,
        const bool hardParityConstraint,
        const double timeLimit,
        const double gapLimit,
        const std::vector<float>& callbackTimeLimit,
        const std::vector<float>& callbackGapLimit,
        double& gap,
        ILPStatus& status,
        std::vector<int>& ilpResults)
{
    using namespace std;

    vector<int> result;

    std::vector<bool> constraintRespected;

    bool allRespectConstraints = false;
    int it = 0;
    do {
        std::cout << std::endl << std::endl << " --- OPTIMIZATION: iteration " << it << " ----" << std::endl;
        try {
            std::vector<bool> isFixed(chartData.subsides.size(), false);
            std::vector<bool> isComputable(chartData.subsides.size(), true);
            for (size_t subsideId = 0; subsideId < chartData.subsides.size(); ++subsideId) {
                if (ilpResults[subsideId] == ILP_IGNORE) {
                    isComputable[subsideId] = false;
                }
                else if (ilpResults[subsideId] >= 0) {
                    isFixed[subsideId] = true;
                }
            }

            result = ilpResults;

            MILPModel model;

            if (!callbackTimeLimit.empty() && !callbackGapLimit.empty() && callbackTimeLimit.size() == callbackGapLimit.size()) {
                model.setEarlyStop(callbackTimeLimit, callbackGapLimit);
            }

#if defined(MILP_NON_VERBOSE) || defined(GUROBI_NON_VERBOSE)
            model.verbose = false;
#endif

            model.timeLimit = timeLimit;
            model.gapLimit = gapLimit;

            // Create variables
            MILPLinExpr obj = 0;
            MILPLinExpr supportObj = 0;

            vector<MILPVar> vars(chartData.subsides.size());

            size_t constraintRespectedId = 0;

            for (size_t subsideId = 0; subsideId < chartData.subsides.size(); subsideId++) {
                if (!isComputable[subsideId]) {
                    continue;
                }

                const ChartSubside& subside = chartData.subsides[subsideId];

                //If it is not a border (free)
                if (!isFixed[subsideId]) {
                    vars[subsideId] = model.addVar(MIN_SUBDIVISION_VALUE, MILP_INFINITY, true);

                    //MIP start: (even) subdivision closest to the isometric one
                    const int startChartId = subside.incidentCharts[0] >= 0 ? subside.incidentCharts[0] : subside.incidentCharts[1];
                    if (startChartId >= 0) {
                        const double sideSubdivision = subside.length / chartEdgeLength[startChartId];
                        if (!hardParityConstraint) {
                            model.setStart(vars[subsideId], std::max(static_cast<double>(MIN_SUBDIVISION_VALUE), std::round(sideSubdivision / 2.0)));
                        }
                        else {
                            model.setStart(vars[subsideId], std::max(2.0, 2.0 * std::round(sideSubdivision / 2.0)));
                        }
                    }
                }
                else if (feasibilityFix && subside.isOnBorder) {
                    const int fixedSize = hardParityConstraint ? ilpResults[subsideId] : std::round(ilpResults[subsideId] / 2.0);

                    vars[subsideId] = model.addVar(std::max(fixedSize - 1, MIN_SUBDIVISION_VALUE), fixedSize + 1, true);
                    model.setStart(vars[subsideId], fixedSize);
                    MILPLinExpr value = vars[subsideId] - fixedSize;

                    supportObj += getMILPAbsInteger(model, value) * FEASIBILITY_FIX_COST * (1.0 / fixedSize);

                    isFixed[subsideId] = false;
                }
            }

            std::cout << chartData.subsides.size() << " subsides!" << std::endl;

            for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
                const Chart& chart = chartData.charts[cId];

                bool computable = true;

                for (size_t i = 0; i < chart.chartSubsides.size(); i++) {
                    const size_t subsideId = chart.chartSubsides[i];
                    if (!isComputable[subsideId])
                        computable = false;
                }

                if (computable && chart.faces.size() > 0) {
                    size_t nSides = chart.chartSides.size();


                    size_t numRegularityTerms = 0;
                    size_t numIsometryTerms = 0;

                    MILPLinExpr regExpr = 0;
                    MILPLinExpr isoExpr = 0;

                    /* ------------------------ ISOMETRY ------------------------ */

                    if (isometry) {
                        for (size_t i = 0; i < chart.chartSubsides.size(); i++) {
                            const size_t subsideId = chart.chartSubsides[i];
                            const ChartSubside& subside = chartData.subsides[subsideId];

                            //If it is not fixed (free)
                            if (!isFixed[subsideId]) {
                                double edgeLength = chartEdgeLength[cId];

                                double sideSubdivision = subside.length / edgeLength;
                                if (!hardParityConstraint) {
                                    sideSubdivision /= 2.0;
                                }

                                sideSubdivision = std::max(static_cast<double>(MIN_SUBDIVISION_VALUE), sideSubdivision);

                                MILPLinExpr value = vars[subsideId] - sideSubdivision;

                                isoExpr += getMILPCostTermContinuous(model, method, value);
                                numIsometryTerms++;
                            }
                        }
                    }


                    /* ------------------------ REGULARITY ------------------------ */

                    std::vector<MILPLinExpr> chartSubsideSum;
                    getChartSubsideSum(chartData, cId, vars, isFixed, ilpResults, hardParityConstraint, chartSubsideSum);

                    for (size_t j = 0; j < nSides; j++) {
                        MILPLinExpr value = 0.0;
                        bool valueComputed = false;

                        //Regularity for quad case
                        if (nSides == 4 && regularityQuadrilaterals) {
                            const MILPLinExpr& subside0Sum = chartSubsideSum[j];
                            const MILPLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];

                            value = subside0Sum - subside2Sum;
                            valueComputed = true;
                        }
                        //Regularity for triangular case
                        else if (nSides == 3 && regularityNonQuadrilaterals) {
                            const MILPLinExpr& subside0Sum = chartSubsideSum[j];
                            const MILPLinExpr& subside1Sum = chartSubsideSum[(j+1)%nSides];
                            const MILPLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];

                            MILPVar c = model.addVar(0, MILP_INFINITY, true);
                            model.addConstr(subside0Sum + 1 <= subside1Sum + subside2Sum + c);

                            value = c;
                            valueComputed = true;
                        }
                        //Regularity for pentagonal case
                        else if (nSides == 5 && regularityNonQuadrilaterals) {
                            const MILPLinExpr& subside0Sum = chartSubsideSum[j];
                            const MILPLinExpr& subside1Sum = chartSubsideSum[(j+1)%nSides];
                            const MILPLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];
                            const MILPLinExpr& subside3Sum = chartSubsideSum[(j+3)%nSides];
                            const MILPLinExpr& subside4Sum = chartSubsideSum[(j+4)%nSides];

                            MILPVar c = model.addVar(0, MILP_INFINITY, true);
                            model.addConstr(subside0Sum + subside1Sum + 1 <= subside2Sum + subside3Sum + subside4Sum + c);

                            value = c;
                            valueComputed = true;
                        }
                        //Regularity for hexagonal case
                        else if (nSides == 6 && regularityNonQuadrilaterals) {
                            const MILPLinExpr& subside0Sum = chartSubsideSum[j];
                            const MILPLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];
                            const MILPLinExpr& subside4Sum = chartSubsideSum[(j+4)%nSides];

                            MILPVar c = model.addVar(0, MILP_INFINITY, true);
                            model.addConstr(subside0Sum + 1 <= subside2Sum + subside4Sum + c);

                            MILPLinExpr parityEquation = subside0Sum + subside2Sum + subside4Sum;
                            MILPVar hexParity = model.addVar(3, MILP_INFINITY, true);
                            MILPVar hexFree = model.addVar(0, 1, true);
                            model.addConstr(hexParity * 2 == parityEquation + hexFree);

                            value = (c + hexFree) / 2.0;
                            valueComputed = true;
                        }

                        if (valueComputed) {
                            if (constraintRespected.empty() || constraintRespected[constraintRespectedId]) {
                                if (nSides == 4) {
                                    regExpr += getMILPCostTermInteger(model, ILPMethod::ABS, value) / nSides;
                                }
                                else {
                                    regExpr += getMILPCostTermInteger(model, ILPMethod::ABS, value) / nSides * regularityNonQuadrilateralsWeight;
                                }

                                numRegularityTerms++;
                            }

                            constraintRespectedId++;
                        }

                    }



                    /* ------------------------ ALIGN SINGULARITIES ------------------------ */

                    if (alignSingularities && (nSides == 3 || nSides == 5 || nSides == 6)) {
                        for (size_t j = 0; j < nSides; j++) {
                            MILPLinExpr value1 = 0.0;
                            MILPLinExpr value2 = 0.0;

                            int currentChartId = cId;
                            int currentChartSideId = j;
                            size_t currentChartNSides = chartData.charts[currentChartId].chartSides.size();

                            bool valueComputed = false;
                            do {
                                const Chart& currentChart = chartData.charts[currentChartId];
                                const ChartSide& currentChartSide = currentChart.chartSides[currentChartSideId];

                                if (currentChartSide.subsides.size() != 1) {
                                    currentChartId = -1;
                                    currentChartSideId = -1;
                                    currentChartNSides = 0;
                                }
                                else {
                                    int adjOppositeSideId = currentChartSideId;
                                    if (currentChartId != static_cast<int>(cId)) {
                                        adjOppositeSideId = (currentChartSideId + 2) % chartData.charts[currentChartId].chartSides.size();
                                    }
                                    const ChartSide& adjOppositeSide = currentChart.chartSides[adjOppositeSideId];

                                    const std::array<int, 2>& incidentCharts = chartData.subsides[adjOppositeSide.subsides[0]].incidentCharts;
                                    const std::array<int, 2>& incidentChartSides = chartData.subsides[adjOppositeSide.subsides[0]].incidentChartSideId;

                                    if (incidentCharts[0] == static_cast<int>(currentChartId)) {
                                        currentChartId = incidentCharts[1];
                                        currentChartSideId = incidentChartSides[1];
                                    }
                                    else if (incidentCharts[1] == static_cast<int>(currentChartId)) {
                                        currentChartId = incidentCharts[0];
                                        currentChartSideId = incidentChartSides[0];
                                    }

                                    if (currentChartId > -1) {
                                        currentChartNSides = chartData.charts[currentChartId].chartSides.size();
                                    }
                                }

                            } while (currentChartId != static_cast<int>(cId) && currentChartId > -1 && currentChartNSides == 4);

                            if (currentChartId != static_cast<int>(cId) && currentChartId > -1 && (currentChartNSides == 3 || currentChartNSides == 5 || currentChartNSides == 6)) {
                                bool currentComputable = true;
                                for (size_t i = 0; i < chartData.charts[currentChartId].chartSubsides.size(); i++) {
                                    const size_t subsideId = chartData.charts[currentChartId].chartSubsides[i];
                                    if (!isComputable[subsideId])
                                        currentComputable = false;
                                }

                                if (currentComputable) {
                                    MILPLinExpr valueDown1 = 0.0;
                                    MILPLinExpr valueUp1 = 0.0;
                                    MILPLinExpr valueDown2 = 0.0;
                                    MILPLinExpr valueUp2 = 0.0;

                                    //Singularity alignment for triangular case
                                    if (nSides == 3) {
                                        const MILPLinExpr& subside0Sum = chartSubsideSum[j];
                                        const MILPLinExpr& subside1Sum = chartSubsideSum[(j+1)%nSides];
                                        const MILPLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];

                                        MILPLinExpr down = subside0Sum + subside2Sum - subside1Sum;
                                        MILPLinExpr up = subside1Sum + subside0Sum - subside2Sum;

                                        valueDown1 = down;
                                        valueUp1 = up;
                                    }
                                    //Singularity alignment for pentagonal case
                                    else if (nSides == 5) {
                                        const MILPLinExpr& subside0Sum = chartSubsideSum[j];
                                        const MILPLinExpr& subside1Sum = chartSubsideSum[(j+1)%nSides];
                                        const MILPLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];
                                        const MILPLinExpr& subside3Sum = chartSubsideSum[(j+3)%nSides];
                                        const MILPLinExpr& subside4Sum = chartSubsideSum[(j+4)%nSides];

                                        MILPLinExpr down = (subside0Sum + subside1Sum + subside2Sum) - (subside3Sum + subside4Sum);
                                        MILPLinExpr up = (subside3Sum + subside4Sum + subside0Sum) - (subside1Sum + subside2Sum);

                                        valueDown1 = down;
                                        valueUp1 = up;
                                    }
                                    //Singularity alignment for hexagonal case
                                    else if (nSides == 6) {
                                        const MILPLinExpr& subside0Sum = chartSubsideSum[j];
                                        const MILPLinExpr& subside2Sum = chartSubsideSum[(j+2)%nSides];
                                        const MILPLinExpr& subside4Sum = chartSubsideSum[(j+4)%nSides];

                                        MILPLinExpr down = subside0Sum + subside2Sum - subside4Sum;
                                        MILPLinExpr up = subside4Sum + subside0Sum - subside2Sum;

                                        valueDown1 = down;
                                        valueUp1 = up;
                                    }




                                    std::vector<MILPLinExpr> adjChartSubsideSum;
                                    getChartSubsideSum(chartData, currentChartId, vars, isFixed, ilpResults, hardParityConstraint, adjChartSubsideSum);

                                    //Singularity alignment for triangular case
                                    if (currentChartNSides == 3) {
                                        const MILPLinExpr& subside0Sum = adjChartSubsideSum[currentChartSideId];
                                        const MILPLinExpr& subside1Sum = adjChartSubsideSum[(currentChartSideId+1)%currentChartNSides];
                                        const MILPLinExpr& subside2Sum = adjChartSubsideSum[(currentChartSideId+2)%currentChartNSides];

                                        MILPLinExpr down = subside0Sum + subside2Sum - subside1Sum;
                                        MILPLinExpr up = subside1Sum + subside0Sum - subside2Sum;

                                        valueDown2 = down;
                                        valueUp2 = up;
                                    }
                                    //Singularity alignment for pentagonal case
                                    else if (currentChartNSides == 5) {
                                        const MILPLinExpr& subside0Sum = adjChartSubsideSum[currentChartSideId];
                                        const MILPLinExpr& subside1Sum = adjChartSubsideSum[(currentChartSideId+1)%currentChartNSides];
                                        const MILPLinExpr& subside2Sum = adjChartSubsideSum[(currentChartSideId+2)%currentChartNSides];
                                        const MILPLinExpr& subside3Sum = adjChartSubsideSum[(currentChartSideId+3)%currentChartNSides];
                                        const MILPLinExpr& subside4Sum = adjChartSubsideSum[(currentChartSideId+4)%currentChartNSides];

                                        MILPLinExpr down = (subside0Sum + subside1Sum + subside2Sum) - (subside3Sum + subside4Sum);
                                        MILPLinExpr up = (subside3Sum + subside4Sum + subside0Sum) - (subside1Sum + subside2Sum);

                                        valueDown2 = down;
                                        valueUp2 = up;
                                    }
                                    //Singularity alignment for hexagonal case
                                    else if (currentChartNSides == 6) {
                                        const MILPLinExpr& subside0Sum = adjChartSubsideSum[currentChartSideId];
                                        const MILPLinExpr& subside2Sum = adjChartSubsideSum[(currentChartSideId+2)%currentChartNSides];
                                        const MILPLinExpr& subside4Sum = adjChartSubsideSum[(currentChartSideId+4)%currentChartNSides];

                                        MILPLinExpr down = subside0Sum + subside2Sum - subside4Sum;
                                        MILPLinExpr up = subside4Sum + subside0Sum - subside2Sum;

                                        valueDown2 = down;
                                        valueUp2 = up;
                                    }

                                    value1 = valueUp1 - valueDown2;
                                    value2 = valueDown1 - valueUp2;
                                    valueComputed = true;
                                }
                            }

                            if (valueComputed) {
                                if (constraintRespected.empty() || constraintRespected[constraintRespectedId]) {
                                    regExpr += (
                                        0.5 * getMILPCostTermInteger(model, ILPMethod::ABS, value1) +
                                        0.5 * getMILPCostTermInteger(model, ILPMethod::ABS, value2)
                                        )  / nSides * alignSingularitiesWeight;
                                }

                                constraintRespectedId++;
                            }
                        }
                    }

                    const double isoWeight = alpha;
                    const double regWeight = (1 - alpha);

                    if (numIsometryTerms > 0)
                        obj += isoWeight * isoExpr / numIsometryTerms;
                    if (numRegularityTerms > 0)
                        obj += regWeight * regExpr / numRegularityTerms;


                    //Even side size sum constraint in a chart
                    if (hardParityConstraint) {
                        MILPLinExpr sumExp = 0;
                        for (const size_t& subsideId : chart.chartSubsides) {
                            const ChartSubside& subside = chartData.subsides[subsideId];

                            if (!isFixed[subsideId]) {
                                assert(isComputable[subsideId] && !isFixed[subsideId]);
                                sumExp += vars[subsideId];
                            }
                            else {                                
                                const int fixedSize = hardParityConstraint ? ilpResults[subsideId] : std::round(ilpResults[subsideId] / 2.0);
                                sumExp += fixedSize;
                            }
                        }

                        MILPVar free = model.addVar(2, MILP_INFINITY, true);
                        model.addConstr(free * 2 == sumExp);
                    }

                    if (chart.chartSides.size() < 3 || chart.chartSides.size() > 6) {
                        std::cout << "Chart " << cId << " has " << chart.chartSides.size() << " sides." << std::endl;
#ifdef ASSERT_FOR_NUMBER_SIDES
                        assert(chart.chartSides.size() >= 3 && chart.chartSides.size() <= 6);
#endif
                    }

                }
            }

            //Set objective function
            model.setObjective(obj + supportObj);

            //model.write("out.lp");

            //Optimize model
            model.optimize();

            for (size_t subsideId = 0; subsideId < chartData.subsides.size(); subsideId++) {
                const ChartSubside& subside = chartData.subsides[subsideId];
                if (isComputable[subsideId]) {
                    if (!isFixed[subsideId]) {
                        assert(isComputable[subsideId] && !isFixed[subsideId]);
                        result[subsideId] = static_cast<int>(std::round(model.getValue(vars[subsideId])));
                    }
                    else {
                        const int fixedSize = hardParityConstraint ? ilpResults[subsideId] : std::round(ilpResults[subsideId] / 2.0);
                        result[subsideId] = fixedSize;
                    }
    
                    if (!hardParityConstraint) {
                        result[subsideId] *= 2;
                    }
                }
            }


            if (it < repeatLosingConstraintsIterations) {
                int numLostConstraintsQuad = 0;
                int numLostConstraintsNonQuad = 0;
                int numLostConstraintsAlign = 0;

                if (constraintRespected.empty()) {
                    constraintRespected.resize(constraintRespectedId, true);
                }

                allRespectConstraints = true;

                constraintRespectedId = 0;

                for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
                    const Chart& chart = chartData.charts[cId];

                    bool computable = true;

                    for (size_t i = 0; i < chart.chartSubsides.size(); i++) {
                        const size_t subsideId = chart.chartSubsides[i];
                        if (!isComputable[subsideId])
                            computable = false;
                    }

                    if (computable && chart.faces.size() > 0) {
                        size_t nSides = chart.chartSides.size();

                        std::vector<int> chartSubsideSumResults;
                        getChartSubsideSumResults(chartData, cId, result, chartSubsideSumResults);

                        for (size_t j = 0; j < nSides; j++) {
                            int value = 0.0;
                            bool valueComputed = false;
                            bool respected = false;

                            //Regularity for quad case
                            if (nSides == 4 && regularityQuadrilaterals) {
                                const int& subside0Sum = chartSubsideSumResults[j];
                                const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];

                                value = subside0Sum - subside2Sum;
                                respected = value == 0;
                                valueComputed = true;
                            }
                            //Regularity for triangular case
                            else if (nSides == 3 && regularityNonQuadrilaterals) {
                                const int& subside0Sum = chartSubsideSumResults[j];
                                const int& subside1Sum = chartSubsideSumResults[(j+1)%nSides];
                                const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];

                                respected = subside0Sum < subside1Sum + subside2Sum;
                                valueComputed = true;
                            }
                            //Regularity for pentagonal case
                            else if (nSides == 5 && regularityNonQuadrilaterals) {
                                const int& subside0Sum = chartSubsideSumResults[j];
                                const int& subside1Sum = chartSubsideSumResults[(j+1)%nSides];
                                const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];
                                const int& subside3Sum = chartSubsideSumResults[(j+3)%nSides];
                                const int& subside4Sum = chartSubsideSumResults[(j+4)%nSides];

                                respected = subside0Sum + subside1Sum < subside2Sum + subside3Sum + subside4Sum;
                                valueComputed = true;
                            }
                            //Regularity for hexagonal case
                            else if (nSides == 6 && regularityNonQuadrilaterals) {
                                const int& subside0Sum = chartSubsideSumResults[j];
                                const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];
                                const int& subside4Sum = chartSubsideSumResults[(j+4)%nSides];

                                respected = subside0Sum < subside2Sum + subside4Sum && ((subside0Sum + subside2Sum + subside4Sum) % 2 == 0);
                                valueComputed = true;
                            }

                            if (valueComputed) {
                                if (!respected) {
                                    if (nSides == 4 && repeatLosingConstraintsQuads) {
                                        numLostConstraintsQuad++;
                                        allRespectConstraints = false;
                                        constraintRespected[constraintRespectedId] = false;
                                    }
                                    else if (((nSides == 3 || nSides == 5 || nSides == 6) && repeatLosingConstraintsNonQuads)) {
                                        numLostConstraintsNonQuad++;
                                        allRespectConstraints = false;
                                        constraintRespected[constraintRespectedId] = false;
                                    }
                                }

                                constraintRespectedId++;
                            }
                        }

                        if (alignSingularities && (nSides == 3 || nSides == 5 || nSides == 6)) {
                            for (size_t j = 0; j < nSides; j++) {
                                int value1 = 0.0;
                                int value2 = 0.0;

                                int currentChartId = cId;
                                int currentChartSideId = j;
                                size_t currentChartNSides = chartData.charts[currentChartId].chartSides.size();

                                bool valueComputed = false;
                                bool respected = false;
                                do {
                                    const Chart& currentChart = chartData.charts[currentChartId];
                                    const ChartSide& currentChartSide = currentChart.chartSides[currentChartSideId];

                                    if (currentChartSide.subsides.size() != 1) {
                                        currentChartId = -1;
                                        currentChartSideId = -1;
                                        currentChartNSides = 0;
                                    }
                                    else {
                                        int adjOppositeSideId = currentChartSideId;
                                        if (currentChartId != static_cast<int>(cId)) {
                                            adjOppositeSideId = (currentChartSideId + 2) % chartData.charts[currentChartId].chartSides.size();
                                        }
                                        const ChartSide& adjOppositeSide = currentChart.chartSides[adjOppositeSideId];

                                        const std::array<int, 2>& incidentCharts = chartData.subsides[adjOppositeSide.subsides[0]].incidentCharts;
                                        const std::array<int, 2>& incidentChartSides = chartData.subsides[adjOppositeSide.subsides[0]].incidentChartSideId;

                                        if (incidentCharts[0] == static_cast<int>(currentChartId)) {
                                            currentChartId = incidentCharts[1];
                                            currentChartSideId = incidentChartSides[1];
                                        }
                                        else if (incidentCharts[1] == static_cast<int>(currentChartId)) {
                                            currentChartId = incidentCharts[0];
                                            currentChartSideId = incidentChartSides[0];
                                        }

                                        if (currentChartId > -1) {
                                            currentChartNSides = chartData.charts[currentChartId].chartSides.size();
                                        }
                                    }

                                } while (currentChartId != static_cast<int>(cId) && currentChartId > -1 && currentChartNSides == 4);

                                if (currentChartId != static_cast<int>(cId) && currentChartId > -1 && (currentChartNSides == 3 || currentChartNSides == 5 || currentChartNSides == 6)) {
                                    bool currentComputable = true;
                                    for (size_t i = 0; i < chartData.charts[currentChartId].chartSubsides.size(); i++) {
                                        const size_t subsideId = chartData.charts[currentChartId].chartSubsides[i];
                                        if (!isComputable[subsideId])
                                            currentComputable = false;
                                    }

                                    if (currentComputable) {

                                        int valueDown1 = 0.0;
                                        int valueUp1 = 0.0;
                                        int valueDown2 = 0.0;
                                        int valueUp2 = 0.0;

                                        //Singularity alignment for triangular case
                                        if (nSides == 3) {
                                            const int& subside0Sum = chartSubsideSumResults[j];
                                            const int& subside1Sum = chartSubsideSumResults[(j+1)%nSides];
                                            const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];

                                            int down = subside0Sum + subside2Sum - subside1Sum;
                                            int up = subside1Sum + subside0Sum - subside2Sum;

                                            valueDown1 = down;
                                            valueUp1 = up;
                                        }
                                        //Singularity alignment for pentagonal case
                                        else if (nSides == 5) {
                                            const int& subside0Sum = chartSubsideSumResults[j];
                                            const int& subside1Sum = chartSubsideSumResults[(j+1)%nSides];
                                            const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];
                                            const int& subside3Sum = chartSubsideSumResults[(j+3)%nSides];
                                            const int& subside4Sum = chartSubsideSumResults[(j+4)%nSides];

                                            int down = (subside0Sum + subside1Sum + subside2Sum) - (subside3Sum + subside4Sum);
                                            int up = (subside3Sum + subside4Sum + subside0Sum) - (subside1Sum + subside2Sum);

                                            valueDown1 = down;
                                            valueUp1 = up;
                                        }
                                        //Singularity alignment for hexagonal case
                                        else if (nSides == 6) {
                                            const int& subside0Sum = chartSubsideSumResults[j];
                                            const int& subside2Sum = chartSubsideSumResults[(j+2)%nSides];
                                            const int& subside4Sum = chartSubsideSumResults[(j+4)%nSides];

                                            int down = subside0Sum + subside2Sum - subside4Sum;
                                            int up = subside4Sum + subside0Sum - subside2Sum;

                                            valueDown1 = down;
                                            valueUp1 = up;
                                        }




                                        std::vector<int> adjChartSubsideSumResults;
                                        getChartSubsideSumResults(chartData, currentChartId, result, adjChartSubsideSumResults);

                                        //Singularity alignment for triangular case
                                        if (currentChartNSides == 3) {
                                            const int& subside0Sum = adjChartSubsideSumResults[currentChartSideId];
                                            const int& subside1Sum = adjChartSubsideSumResults[(currentChartSideId+1)%currentChartNSides];
                                            const int& subside2Sum = adjChartSubsideSumResults[(currentChartSideId+2)%currentChartNSides];

                                            int down = subside0Sum + subside2Sum - subside1Sum;
                                            int up = subside1Sum + subside0Sum - subside2Sum;

                                            valueDown2 = down;
                                            valueUp2 = up;
                                        }
                                        //Singularity alignment for pentagonal case
                                        else if (currentChartNSides == 5) {
                                            const int& subside0Sum = adjChartSubsideSumResults[currentChartSideId];
                                            const int& subside1Sum = adjChartSubsideSumResults[(currentChartSideId+1)%currentChartNSides];
                                            const int& subside2Sum = adjChartSubsideSumResults[(currentChartSideId+2)%currentChartNSides];
                                            const int& subside3Sum = adjChartSubsideSumResults[(currentChartSideId+3)%currentChartNSides];
                                            const int& subside4Sum = adjChartSubsideSumResults[(currentChartSideId+4)%currentChartNSides];

                                            int down = (subside0Sum + subside1Sum + subside2Sum) - (subside3Sum + subside4Sum);
                                            int up = (subside3Sum + subside4Sum + subside0Sum) - (subside1Sum + subside2Sum);

                                            valueDown2 = down;
                                            valueUp2 = up;
                                        }
                                        //Singularity alignment for hexagonal case
                                        else if (currentChartNSides == 6) {
                                            const int& subside0Sum = adjChartSubsideSumResults[currentChartSideId];
                                            const int& subside2Sum = adjChartSubsideSumResults[(currentChartSideId+2)%currentChartNSides];
                                            const int& subside4Sum = adjChartSubsideSumResults[(currentChartSideId+4)%currentChartNSides];

                                            int down = subside0Sum + subside2Sum - subside4Sum;
                                            int up = subside4Sum + subside0Sum - subside2Sum;

                                            valueDown2 = down;
                                            valueUp2 = up;
                                        }

                                        value1 = valueUp1 - valueDown2;
                                        value2 = valueDown1 - valueUp2;

                                        respected = value1 == 0 && value2 == 0;

                                        valueComputed = true;
                                    }
                                }

                                if (valueComputed) {
                                    if (!respected) {
                                        if (((nSides == 3 || nSides == 5 || nSides == 6) && repeatLosingConstraintsAlign)) {
                                            allRespectConstraints = false;
                                            constraintRespected[constraintRespectedId] = false;
                                            numLostConstraintsAlign++;
                                        }
                                    }

                                    constraintRespectedId++;
                                }
                            }
                        }
                    }
                }

                std::cout << "Lost contraints quad: " << numLostConstraintsQuad << std::endl;
                std::cout << "Lost contraints non quad: " << numLostConstraintsNonQuad << std::endl;
                std::cout << "Lost contraints alignment: " << numLostConstraintsAlign << std::endl << std::endl << std::endl;
            }

#if !defined(MILP_NON_VERBOSE) && !defined(GUROBI_NON_VERBOSE)
            cout << "Support obj: " << model.getValue(supportObj) << std::endl;
            cout << "Obj: " << model.getValue(obj) << std::endl;
#endif
            gap = model.getMIPGap();

            status = ILPStatus::SOLUTIONFOUND;



            for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
                const Chart& chart = chartData.charts[cId];

                bool computable = true;

                for (size_t i = 0; i < chart.chartSubsides.size(); i++) {
                    const size_t subsideId = chart.chartSubsides[i];
                    if (!isComputable[subsideId])
                        computable = false;
                }

                if (computable && chart.faces.size() > 0) {
                    int sizeSum = 0;
                    for (const size_t& subsideId : chart.chartSubsides) {
                        sizeSum += result[subsideId];
                    }

                    if (sizeSum % 2 == 1) {
                        std::cout << "Error not even, chart: " << cId << " -> ";
                        for (const size_t& subsideId : chart.chartSubsides) {
                            std::cout << result[subsideId] << " ";
                        }
                        std::cout << " = " << sizeSum << std::endl;

                        status = ILPStatus::SOLUTIONWRONG;
                    }
                }
            }

        } catch(const MILPException& e) {
            cout << e.what() << std::endl;

            status = ILPStatus::INFEASIBLE;
        }

        it++;
    }
    while (status == ILPStatus::SOLUTIONFOUND && it < repeatLosingConstraintsIterations + 1 && !allRespectConstraints);

    ilpResults = result;
}

inline void getChartSubsideSum(
        const ChartData& chartData,
        const size_t& cId,
        const std::vector<MILPVar>& vars,
        const std::vector<bool>& isFixed,
        const std::vector<int>& ilpResults,
        const bool hardParityConstraint,
        std::vector<MILPLinExpr>& chartSubsideSum)
{
    const Chart& chart = chartData.charts[cId];

    size_t nSides = chart.chartSides.size();

    chartSubsideSum.resize(nSides, 0.0);

    for (size_t i = 0; i < nSides; i++) {
        const ChartSide& side = chart.chartSides[i];
        MILPLinExpr subsideSum = 0;
        for (const size_t& subsideId : side.subsides) {
            const ChartSubside& subside = chartData.subsides[subsideId];

            if (!isFixed[subsideId]) {
                assert(ilpResults[subsideId] != ILP_IGNORE && !isFixed[subsideId]);
                subsideSum += vars[subsideId];
            }
            else {                
                const int fixedSize = hardParityConstraint ? ilpResults[subsideId] : std::round(ilpResults[subsideId] / 2.0);
                subsideSum += fixedSize;
            }
        }

        chartSubsideSum[i] = subsideSum;
    }
}

inline void getChartSubsideSumResults(
        const ChartData& chartData,
        const size_t& cId,
        const std::vector<int>& results,
        std::vector<int>& chartSubsideSum)
{
    const Chart& chart = chartData.charts[cId];

    size_t nSides = chart.chartSides.size();

    chartSubsideSum.resize(nSides, 0.0);

    for (size_t i = 0; i < nSides; i++) {        
        const ChartSide& side = chart.chartSides[i];
        int subsideSum = 0;
        for (const size_t& subsideId : side.subsides) {
            subsideSum += results[subsideId];
        }

        chartSubsideSum[i] = subsideSum;
    }
}

inline MILPLinExpr getMILPCostTermInteger(MILPModel& model, const ILPMethod& method, const MILPLinExpr& value)
{
    MILPLinExpr expr;

    if (method == LEASTSQUARES) {
        expr = model.addSquare(value);
    }
    else {
        expr = getMILPAbsInteger(model, value);
    }

    return expr;
}

inline MILPLinExpr getMILPCostTermContinuous(MILPModel& model, const ILPMethod& method, const MILPLinExpr& value)
{
    MILPLinExpr expr;

    if (method == LEASTSQUARES) {
        //value is (subdivision variable - constant): the square is exact at the integer values
        expr = model.addSquare(value);
    }
    else {
        expr = getMILPAbsContinuous(model, value);
    }

    return expr;
}

inline MILPLinExpr getMILPAbsInteger(MILPModel& model, const MILPLinExpr& value) {
    return model.addAbs(value, true);
}

inline MILPLinExpr getMILPAbsContinuous(MILPModel& model, const MILPLinExpr& value) {
    return model.addAbs(value, false);
}

}
}


