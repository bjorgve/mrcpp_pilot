#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MRCPP/Plotter"

#include "MRCPP/MWFunctions"
#include "MRCPP/Gaussians"
#include "MRCPP/MWOperators"
#include "MRCPP/core/LegendreBasis.h"
#include "MRCPP/core/QuadratureCache.h"
#include "MRCPP/core/InterpolatingBasis.h"
#include "MRCPP/core/ScalingBasis.h"
#include "../include/config.h"

#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// This is your personal sandbox to test out ideas.
// See the examples/ directory for inspiration.
// Please do not commit your pilot code to git.

using namespace mrcpp;
const auto min_scale = 0;
const auto max_depth = 25;

const auto order = 7;
const auto prec = 1.0e-5;
const auto D = 3; // Dimensions
int main(int argc, char **argv) {
    Timer timer;

    // Plots
    // Initialize printing
    int printlevel = 0;
    Printer::init(printlevel);
    Printer::printEnvironment();
    Printer::printHeader(0, "MRCPP pilot code");



    std::string file;
    std::string path = MW_FILTER_DIR;
    file = path + "/L_ph_deriv_1.txt";

    std::ifstream ifs(file.c_str());
    std::string line;


    std::ofstream file_w("test.txt");
    file_w.precision(16);
    file_w.setf(std::ios::scientific);

// auto give_basis(int order) {
//     if (order == 0) return 1.0;
//     return LegendreBasis(order);
// }

    for (auto kp1 = 2; kp1 < 5; kp1++) {

        auto basis = LegendreBasis(kp1-1);

        getQuadratureCache(qCache);
        const VectorXd &roots = qCache.getRoots(kp1);
        const VectorXd &weights = qCache.getWeights(kp1);

        getline(ifs, line);
        if (kp1 != std::stoi(line)) MSG_FATAL("Something wrong");
        MatrixXd data = MatrixXd::Zero(3*kp1, kp1);
        for (int i = 0; i < 3*kp1; i++) {
            getline(ifs, line);
            std::istringstream iss(line);
            for (int j = 0; j < kp1; j++) {
                iss >> data(i, j);
            }
        }

        MatrixXd S_p1 = data.block(0*kp1, 0, kp1, kp1);
        MatrixXd S_0 = data.block(1*kp1, 0, kp1, kp1);
        MatrixXd S_m1 = data.block(2*kp1, 0, kp1, kp1);

        MatrixXd S = MatrixXd::Zero(kp1, kp1);

        for (auto row = 0; row < kp1; row++) {
            for (auto column = 0; column < kp1; column++) {
                S(row, column) = basis.getFunc(row).evalf(roots[column])*weights[column];
            }
        }

        // Writing to new file
        file_w << kp1 << '\n';
        file_w << S.inverse()*S_p1*S << '\n';
        file_w << S.inverse()*S_0*S<< '\n';
        file_w << S.inverse()*S_m1*S << '\n';


        // println(0, data);

    }
    return 0;
}
